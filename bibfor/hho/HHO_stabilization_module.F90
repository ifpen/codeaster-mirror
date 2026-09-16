! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

module HHO_stabilization_module
!
    use HHO_basis_module
    use HHO_massmat_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_tracemat_module
    use HHO_type
    use HHO_utils_module
    use HHO_matrix_module
    use HHO_algebra_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterf_debug.h"
#include "asterfort/assert.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/utmess.h"
#include "blas/dgemm.h"
#include "blas/dposv.h"
#include "blas/dpotrf.h"
#include "blas/dpotrs.h"
#include "blas/dsyrk.h"
!
!---------------------------------------------------------------------------------------------------
!  HHO - stabilization
!
!  This module contains all the routines to compute the stabilzation operator for HHO methods
!
!---------------------------------------------------------------------------------------------------
!
    public :: hhoStabScal, hhoStabVec, hdgStabScal, hdgStabVec, hhoStabSymVec, MatScal2Vec
!    private ::
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoStabScal(hhoCell, hhoData, gradrec, stab)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(in) :: gradrec
        type(HHO_matrix), intent(out) :: stab
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the hho stabilization of a scalar function
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   In gradrec      : matrix of the gradient reconstruction
!   Out stab        : matrix of stabilization (lhs member for laplacian problem)
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_Face) :: hhoFace
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_massmat_cell) :: massMat
        type(HHO_massmat_face) :: faceMass
        type(HHO_matrix) :: proj1
        real(kind=8) :: invH
        real(kind=8), dimension(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL) :: M1, M2
        real(kind=8), dimension(MSIZE_FACE_SCAL, MSIZE_FACE_SCAL) :: piKF
        real(kind=8), dimension(MSIZE_FACE_SCAL, MSIZE_CELL_SCAL) :: MR1, traceMat
        real(kind=8), dimension(MSIZE_FACE_SCAL, MSIZE_TDOFS_SCAL) :: proj2, proj3, TMP
        integer(kind=8) :: dimMassMat, ifromM2, itoM2, colsM2, i, j
        integer(kind=8) :: cbs, fbs, total_dofs, iface, fromFace, toFace, cell_offset
        blas_int :: b_n, b_nrhs, b_lda, b_ldb, info
        blas_int :: b_k, b_ldc, b_m
! --------------------------------------------------------------------------------------------------
!
        ASSERT(hhoCell%l_face_init)
! -- init cell basis
        call hhoBasisCell%initialize(hhoCell)
!
! -- number of dofs
        call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        cell_offset = total_dofs-cbs+1
!
! -- compute cell mass matrix
        call massMat%compute(hhoCell, 0, hhoData%face_degree()+1)
        dimMassMat = massMat%nrows
!
! -- Range
        call hhoBasisCell%BSRange(1, hhoData%face_degree()+1, ifromM2, itoM2)
        colsM2 = hhoBasisCell%BSSize(1, hhoData%face_degree()+1)
!
! -- extract M2:
        M2 = 0.d0
        M2(1:cbs, 1:colsM2) = massMat%m(1:cbs, ifromM2:itoM2)
!
! -- Verif size
        ASSERT(MSIZE_CELL_SCAL >= colsM2 .and. MSIZE_TDOFS_SCAL >= cbs)
!
        call stab%initialize(total_dofs, total_dofs, 0.0)
        call proj1%initialize(cbs, total_dofs, 0.0)
!
! -- Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)
! -- compute proj1: Step 1: compute \pi_T^k p_T^k v (third term).
!
! -- Compute proj1 = -M2 * gradrec
        b_ldc = to_blas_int(proj1%max_nrows)
        b_ldb = to_blas_int(gradrec%max_nrows)
        b_lda = to_blas_int(MSIZE_CELL_SCAL)
        b_m = to_blas_int(cbs)
        b_n = to_blas_int(total_dofs)
        b_k = to_blas_int(colsM2)
        call dgemm('N', 'N', b_m, b_n, b_k, &
                   -1.d0, M2, b_lda, gradrec%m, b_ldb, &
                   0.d0, proj1%m, b_ldc)
!
        if (.not. massMat%isIdentity) then
!
! -- extract M1:
            M1 = 0.d0
            M1(1:cbs, 1:cbs) = massMat%m(:cbs, 1:cbs)
! -- Solve proj1 = M1^-1 * proj1
! -- Verif strange bug if info neq 0 in entry
            info = 0
            b_n = to_blas_int(cbs)
            b_nrhs = to_blas_int(total_dofs)
            b_lda = to_blas_int(MSIZE_CELL_SCAL)
            b_ldb = to_blas_int(proj1%max_nrows)
            call dposv('U', b_n, b_nrhs, M1, b_lda, &
                       proj1%m, b_ldb, info)
!
! - Sucess ?
            if (info .ne. 0) then
                call utmess('F', 'HHO1_4')
            end if
!
        end if
!
! --  Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
! -- Compute proj1 = proj1 + I_Cell
        do i = 1, cbs
            proj1%m(i, cell_offset-1+i) = proj1%m(i, cell_offset-1+i)+1.d0
        end do
!
! Step 3: project on faces (eqn. 21)
!
        toFace = 0
! -- Loop on the faces
        do iface = 1, hhoCell%nbfaces
            hhoFace = hhoCell%faces(iface)
            invH = 1.d0/hhoFace%diameter
            fromFace = toFace+1
            toFace = fromFace+fbs-1
!
! ----- Compute face mass matrix
            call faceMass%compute(hhoFace, 0, hhoData%face_degree())
!
! ----- Compute trace mass matrix
            call hhoTraceMatScal(hhoCell, 0, hhoData%face_degree()+1, hhoFace, 0, &
                                 hhoData%face_degree(), traceMat)
!
            if (.not. faceMass%isIdentity) then
!
! ---- Factorize face Mass
                piKF = 0.d0
                piKF(1:fbs, 1:fbs) = faceMass%m(1:fbs, 1:fbs)
! ---- Verif strange bug if info neq 0 in entry
                info = 0
                b_n = to_blas_int(fbs)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                call dpotrf('U', b_n, piKF, b_lda, info)
!
! --- Sucess ?
                if (info .ne. 0) then
                    call utmess('F', 'HHO1_4')
                end if
            end if
!
! ----  Step 3a: \pi_F^k( v_F - p_T^k v )
            MR1 = 0.d0
            MR1(1:fbs, 1:colsM2) = traceMat(1:fbs, ifromM2:itoM2)
!
! ----  compute proj2 = MR1 * gradrec
            proj2 = 0.d0
            b_ldc = to_blas_int(MSIZE_FACE_SCAL)
            b_ldb = to_blas_int(gradrec%max_nrows)
            b_lda = to_blas_int(MSIZE_FACE_SCAL)
            b_m = to_blas_int(fbs)
            b_n = to_blas_int(total_dofs)
            b_k = to_blas_int(colsM2)
            call dgemm('N', 'N', b_m, b_n, b_k, &
                       1.d0, MR1, b_lda, gradrec%m, b_ldb, &
                       0.d0, proj2, b_ldc)
!
            if (.not. faceMass%isIdentity) then
!
! ---- Solve proj2 = pikF^-1 * proj2
! ---- Verif strange bug if info neq 0 in entry
                info = 0
                b_n = to_blas_int(fbs)
                b_nrhs = to_blas_int(total_dofs)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                b_ldb = to_blas_int(MSIZE_FACE_SCAL)
                call dpotrs('U', b_n, b_nrhs, piKF, b_lda, &
                            proj2, b_ldb, info)
!
! --- Sucess ?
                if (info .ne. 0) then
                    call utmess('F', 'HHO1_4')
                end if
            end if
!
! ---- Compute proj2 -= I_F
            i = 1
            do j = fromFace, toFace
                proj2(i, j) = proj2(i, j)-1.d0
                i = i+1
            end do
!
! ---- Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
!
!
! ---- Compute proj3 = traceMat * proj1
            proj3 = 0.d0
            b_ldc = to_blas_int(MSIZE_FACE_SCAL)
            b_ldb = to_blas_int(proj1%max_nrows)
            b_lda = to_blas_int(MSIZE_FACE_SCAL)
            b_m = to_blas_int(fbs)
            b_n = to_blas_int(total_dofs)
            b_k = to_blas_int(cbs)
            call dgemm('N', 'N', b_m, b_n, b_k, &
                       1.d0, traceMat, b_lda, proj1%m, b_ldb, &
                       0.d0, proj3, b_ldc)
!
            if (.not. faceMass%isIdentity) then
!
! ---- Solve proj3 = pikF^-1 * proj3
                info = 0
                b_n = to_blas_int(fbs)
                b_nrhs = to_blas_int(total_dofs)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                b_ldb = to_blas_int(MSIZE_FACE_SCAL)
                call dpotrs('U', b_n, b_nrhs, piKF, b_lda, &
                            proj3, b_ldb, info)
!
! --- -Success ?
! ---- Verif strange bug if info neq 0 in entry
                if (info .ne. 0) then
                    call utmess('F', 'HHO1_4')
                end if
            end if
!
! ---- proj3 = proj3 + proj2
            proj3(1:fbs, 1:total_dofs) = proj3(1:fbs, 1:total_dofs)+proj2(1:fbs, 1:total_dofs)
!
            if (.not. faceMass%isIdentity) then
!
! ---- Compute TMP = faceMass * proj3
                TMP = 0.d0
                b_ldc = to_blas_int(MSIZE_FACE_SCAL)
                b_ldb = to_blas_int(MSIZE_FACE_SCAL)
                b_lda = to_blas_int(faceMass%max_nrows)
                b_m = to_blas_int(fbs)
                b_n = to_blas_int(total_dofs)
                b_k = to_blas_int(fbs)
                call dgemm('N', 'N', b_m, b_n, b_k, &
                           1.d0, faceMass%m, b_lda, proj3, b_ldb, &
                           0.d0, TMP, b_ldc)
!
! ---- Compute stab += invH * proj3**T * TMP
                b_ldc = to_blas_int(stab%max_nrows)
                b_ldb = to_blas_int(MSIZE_FACE_SCAL)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                b_m = to_blas_int(total_dofs)
                b_n = to_blas_int(total_dofs)
                b_k = to_blas_int(fbs)
                call dgemm('T', 'N', b_m, b_n, b_k, &
                           invH, proj3, b_lda, TMP, b_ldb, &
                           1.d0, stab%m, b_ldc)
            else
!
! ---- Compute stab += invH * proj3**T * proj3 since faceMass = Id
                b_ldc = to_blas_int(stab%max_nrows)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                b_n = to_blas_int(total_dofs)
                b_k = to_blas_int(fbs)
                call dsyrk('U', 'T', b_n, b_k, &
                           invH, proj3, b_lda, &
                           1.d0, stab%m, b_ldc)
            end if
!
        end do
!
        if (faceMass%isIdentity) then
            call stab%copySymU()
        end if
!
        call proj1%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoStabVec(hhoCell, hhoData, gradrec_scal, stab)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(in) :: gradrec_scal
        type(HHO_matrix), intent(out) :: stab
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the hho stabilization of a vectorial function
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   In gradrec_scal : matrix of the gradient reconstruction of a scalar function
!   Out stab        : matrix of stabilization (lhs member for vectorial laplacian problem)
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_matrix) :: stab_scal
! --------------------------------------------------------------------------------------------------
!
! -- compute scalar stabilization
        call hhoStabScal(hhoCell, hhoData, gradrec_scal, stab_scal)
!
! -- copy the scalar stabilization in the vectorial stabilization
        call MatScal2Vec(hhoCell, hhoData, stab_scal, stab)
        call stab_scal%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoStabSymVec(hhoCell, hhoData, gradrec, stab)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(in) :: gradrec
        type(HHO_matrix), intent(out) :: stab
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the hho stabilization of a vectorial function
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   In gradrec      : matrix of the symmetric gradient reconstruction
!   Out stab        : matrix of stabilization (lhs member for laplacian problem)
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_Face) :: hhoFace
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_massmat_cell) :: massMat
        type(HHO_massmat_face) :: faceMass
        real(kind=8) :: invH
        real(kind=8), dimension(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL) :: M1, M2
        real(kind=8), dimension(MSIZE_FACE_SCAL, MSIZE_FACE_SCAL) :: piKF
        type(HHO_matrix) :: proj1, proj2, proj3, TMP
        real(kind=8), dimension(MSIZE_FACE_SCAL, MSIZE_CELL_SCAL) :: MR1, MR2, traceMat
        integer(kind=8) :: dimMassMat, ifromM1, itoM1, ifromM2, itoM2, colsM2, i, j, idir
        integer(kind=8) :: cbs, fbs, total_dofs, iface, fromFace, toFace, faces_dofs_comp
        integer(kind=8) :: ifromGrad, itoGrad, ifromProj, itoProj, fbs_comp, faces_dofs
        integer(kind=8) :: dimM1
        blas_int :: b_n, b_nrhs, b_lda, b_ldb, info
        blas_int :: b_k, b_ldc, b_m
! --------------------------------------------------------------------------------------------------
!
        ASSERT(hhoCell%l_face_init)
! -- init cell basis
        call hhoBasisCell%initialize(hhoCell)
!
! -- number of dofs
        call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
        fbs_comp = fbs/hhoCell%ndim
        faces_dofs = total_dofs-cbs
        faces_dofs_comp = faces_dofs/hhoCell%ndim
!
! -- compute cell mass matrix
        call massMat%compute(hhoCell, 0, hhoData%face_degree()+1)
        dimMassMat = massMat%nrows
!
! -- Range
        call hhoBasisCell%BSRange(0, hhoData%cell_degree(), ifromM1, itoM1)
        dimM1 = hhoBasisCell%BSSize(0, hhoData%cell_degree())
        ASSERT(dimM1*hhoCell%ndim == cbs)
        call hhoBasisCell%BSRange(1, hhoData%face_degree()+1, ifromM2, itoM2)
        colsM2 = hhoBasisCell%BSSize(1, hhoData%face_degree()+1)
!
! -- factorize M1
        if (.not. massMat%isIdentity) then
!
! -- extract M1:
            M1 = 0.d0
            M1(1:dimM1, 1:dimM1) = massMat%m(ifromM1:itoM1, ifromM1:itoM1)
            info = 0
            b_n = to_blas_int(dimM1)
            b_lda = to_blas_int(MSIZE_CELL_SCAL)
            call dpotrf('U', b_n, M1, b_lda, info)
!
! -- Sucess ?
            if (info .ne. 0) then
                call utmess('F', 'HHO1_4')
            end if
        end if
!
! -- extract M2:
        M2 = 0.d0
        M2(1:dimM1, 1:colsM2) = massMat%m(ifromM1:itoM1, ifromM2:itoM2)
!
! -- Verif size
        ASSERT(MSIZE_CELL_SCAL >= colsM2 .and. MSIZE_TDOFS_SCAL >= cbs)
!
        call stab%initialize(total_dofs, total_dofs, 0.0)
        call proj1%initialize(cbs+3, total_dofs, 0.d0)
        call proj2%initialize(fbs_comp, total_dofs, 0.d0)
        call proj3%initialize(fbs_comp, total_dofs, 0.d0)
        call TMP%initialize(fbs_comp, total_dofs, 0.d0)
!
! -- Build \pi_F^k (v_F - P_T^K v) equations (21) and (22)
! -- compute proj1: Step 1: compute \pi_T^k p_T^k v (third term).
!
! -- Compute proj1 = -M2 * gradrec
!
        do idir = 1, hhoCell%ndim
!
            ifromGrad = (idir-1)*colsM2+1
            itoGrad = ifromGrad+colsM2-1
            ifromProj = (idir-1)*dimM1+1
            itoProj = ifromProj+dimM1-1
!
            b_ldc = to_blas_int(dimM1)
            b_ldb = to_blas_int(colsM2)
            b_lda = to_blas_int(MSIZE_CELL_SCAL)
            b_m = to_blas_int(dimM1)
            b_n = to_blas_int(total_dofs)
            b_k = to_blas_int(colsM2)
            call dgemm('N', 'N', b_m, b_n, b_k, &
                       -1.d0, M2, b_lda, gradrec%m(ifromGrad:itoGrad, 1:total_dofs), b_ldb, &
                       0.d0, proj1%m(ifromProj:itoProj, 1:total_dofs), b_ldc)
!
            if (.not. massMat%isIdentity) then
! -- Solve proj1 = M1^-1 * proj1
! -- Verif strange bug if info neq 0 in entry
                info = 0
                b_n = to_blas_int(dimM1)
                b_nrhs = to_blas_int(total_dofs)
                b_lda = to_blas_int(MSIZE_CELL_SCAL)
                b_ldb = to_blas_int(dimM1)
                call dpotrs('U', b_n, b_nrhs, M1, b_lda, &
                            proj1%m(ifromProj:itoProj, 1:total_dofs), b_ldb, info)
!
! -- Sucess ?
                if (info .ne. 0) then
                    call utmess('F', 'HHO1_4')
                end if
            end if
!
        end do
!
! --  Step 2: v_T - \pi_T^k p_T^k v (first term minus third term)
! -- Compute proj1 = proj1 + I_Cell
        do i = 1, cbs
            proj1%m(i, faces_dofs+i) = proj1%m(i, faces_dofs+i)+1.d0
        end do
!
! Step 3: project on faces (eqn. 21)
!
! -- Loop on the faces
        do iface = 1, hhoCell%nbfaces
            hhoFace = hhoCell%faces(iface)
            invH = 1.d0/hhoFace%diameter
!
! ----- Compute face mass matrix
            call faceMass%compute(hhoFace, 0, hhoData%face_degree())
!
! ----- Compute trace mass matrix
            call hhoTraceMatScal(hhoCell, 0, hhoData%face_degree()+1, hhoFace, 0, &
                                 hhoData%face_degree(), traceMat)
!
            if (.not. faceMass%isIdentity) then
!
! ---- Factorize face Mass
                piKF = 0.d0
                piKF(1:fbs_comp, 1:fbs_comp) = faceMass%m(1:fbs_comp, 1:fbs_comp)
! ---- Verif strange bug if info neq 0 in entry
                info = 0
                b_n = to_blas_int(fbs_comp)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                call dpotrf('U', b_n, piKF, b_lda, info)
!
! --- Sucess ?
                if (info .ne. 0) then
                    call utmess('F', 'HHO1_4')
                end if
            end if
!
! ----  Step 3a: \pi_F^k( v_F - p_T^k v )
            MR1 = 0.d0
            MR1(1:fbs_comp, 1:colsM2) = traceMat(1:fbs_comp, ifromM2:itoM2)
!
            MR2 = 0.d0
            MR2(1:fbs_comp, 1:dimM1) = traceMat(1:fbs_comp, ifromM1:itoM1)
!
            do idir = 1, hhoCell%ndim
!
                call proj2%setValue(0.0)
                call proj3%setValue(0.0)
!
                ifromGrad = (idir-1)*colsM2+1
                itoGrad = ifromGrad+colsM2-1
!
! ----  compute proj2 = MR1 * gradrec
                b_ldc = to_blas_int(proj2%max_nrows)
                b_ldb = to_blas_int(colsM2)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                b_m = to_blas_int(fbs_comp)
                b_n = to_blas_int(total_dofs)
                b_k = to_blas_int(colsM2)
                call dgemm('N', 'N', b_m, b_n, b_k, &
                           1.d0, MR1, b_lda, gradrec%m(ifromGrad:itoGrad, 1:total_dofs), b_ldb, &
                           0.d0, proj2%m, b_ldc)
!
                if (.not. faceMass%isIdentity) then
!
! ---- Solve proj2 = pikF^-1 * proj2
! ---- Verif strange bug if info neq 0 in entry
                    info = 0
                    b_n = to_blas_int(fbs_comp)
                    b_nrhs = to_blas_int(total_dofs)
                    b_lda = to_blas_int(MSIZE_FACE_SCAL)
                    b_ldb = to_blas_int(proj2%max_nrows)
                    call dpotrs('U', b_n, b_nrhs, piKF, b_lda, &
                                proj2%m, b_ldb, info)
!
! --- Sucess ?
                    if (info .ne. 0) then
                        call utmess('F', 'HHO1_4')
                    end if
                end if
!
! ---- Compute proj2 -= I_F
                fromFace = (idir-1)*fbs_comp+(iface-1)*fbs+1
                toFace = fromFace+fbs_comp-1
                i = 1
                do j = fromFace, toFace
                    proj2%m(i, j) = proj2%m(i, j)-1.d0
                    i = i+1
                end do
!
! ---- Step 3b: \pi_F^k( v_T - \pi_T^k p_T^k v )
! ---- Compute proj3 = MR2 * proj1
!
                ifromProj = (idir-1)*dimM1+1
                itoProj = ifromProj+dimM1-1
!
                b_ldc = to_blas_int(proj3%max_nrows)
                b_ldb = to_blas_int(dimM1)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                b_m = to_blas_int(fbs_comp)
                b_n = to_blas_int(total_dofs)
                b_k = to_blas_int(dimM1)
                call dgemm('N', 'N', b_m, b_n, b_k, &
                           1.d0, MR2, b_lda, proj1%m(ifromProj:itoProj, 1:total_dofs), b_ldb, &
                           0.d0, proj3%m, b_ldc)
!
                if (.not. faceMass%isIdentity) then
!
! ---- Solve proj3 = pikF^-1 * proj3
                    info = 0
                    b_n = to_blas_int(fbs_comp)
                    b_nrhs = to_blas_int(total_dofs)
                    b_lda = to_blas_int(MSIZE_FACE_SCAL)
                    b_ldb = to_blas_int(proj3%max_nrows)
                    call dpotrs('U', b_n, b_nrhs, piKF, b_lda, &
                                proj3%m, b_ldb, info)
!
! --- -Success ?
! ---- Verif strange bug if info neq 0 in entry
                    if (info .ne. 0) then
                        call utmess('F', 'HHO1_4')
                    end if
                end if
!
! ---- proj3 = proj3 + proj2
                call proj3%add(proj2)
!
                if (.not. faceMass%isIdentity) then
!
! ---- Compute TMP = faceMass * proj3
                    call TMP%setValue(0.d0)
                    b_ldc = to_blas_int(TMP%max_nrows)
                    b_ldb = to_blas_int(proj3%max_nrows)
                    b_lda = to_blas_int(MSIZE_FACE_SCAL)
                    b_m = to_blas_int(fbs_comp)
                    b_n = to_blas_int(total_dofs)
                    b_k = to_blas_int(fbs_comp)
                    call dgemm('N', 'N', b_m, b_n, b_k, &
                               1.d0, faceMass%m, b_lda, proj3%m, b_ldb, &
                               0.d0, TMP%m, b_ldc)
!
! ---- Compute stab += invH * proj3**T * TMP
                    call hho_dgemm_TN(invH, proj3, TMP, 1.d0, stab)
                else
!
! ---- Compute stab += invH * proj3**T * proj3
                    call hho_dgemm_TN(invH, proj3, proj3, 1.d0, stab)
                end if

!
            end do
        end do
!
        call proj1%free()
        call proj2%free()
        call proj3%free()
        call TMP%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hdgStabScal(hhoCell, hhoData, stab)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(out) :: stab
!
! --------------------------------------------------------------------------------------------------
!   HHO - HDG type stabilisation 1/h_F(v_F - pi^k_F(vT))_F
!
!   Compute the hdg stabilization of a scalar function
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   Out stab        : matrix of stabilization (lhs member for laplacian problem)
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_Face) :: hhoFace
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_massmat_face) :: faceMass
        real(kind=8) :: invH
        real(kind=8), dimension(MSIZE_FACE_SCAL, MSIZE_FACE_SCAL) :: invM
        real(kind=8), dimension(MSIZE_FACE_SCAL, MSIZE_CELL_SCAL) :: traceMat, piKF
        real(kind=8), dimension(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL) :: S_TT
        integer(kind=8) :: cbs, fbs, total_dofs, iface, fromFace, toFace
        integer(kind=8) :: cell_offset
        blas_int :: b_n, b_nrhs, b_lda, b_ldb, info
        blas_int :: b_k, b_ldc, b_m
! --------------------------------------------------------------------------------------------------
!
        ASSERT(hhoCell%l_face_init)
! -- init cell basis
        call hhoBasisCell%initialize(hhoCell)
!
! -- number of dofs
        call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        cell_offset = total_dofs-cbs+1
!
        call stab%initialize(total_dofs, total_dofs, 0.0)
        piKF = 0.d0
        S_TT = 0.d0
!
! Step 3: project on faces (eqn. 21)
        toFace = 0
!
! -- Loop on the faces
        do iface = 1, hhoCell%nbfaces
            hhoFace = hhoCell%faces(iface)
            invH = 1.d0/hhoFace%diameter
            fromFace = toFace+1
            toFace = fromFace+fbs-1
!
! ----- Compute face mass matrix
            call faceMass%compute(hhoFace, 0, hhoData%face_degree())
!
! ----- Compute trace mass matrix
            call hhoTraceMatScal(hhoCell, 0, hhoData%cell_degree(), hhoFace, 0, &
                                 hhoData%face_degree(), traceMat)
!
            piKF(1:fbs, 1:cbs) = traceMat(1:fbs, 1:cbs)
!
            if (.not. faceMass%isIdentity) then
!
! ---- Factorize face Mass
                invM(1:fbs, 1:fbs) = faceMass%m(1:fbs, 1:fbs)
!
                b_n = to_blas_int(fbs)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                call dpotrf('U', b_n, invM, b_lda, info)
!
! --- Sucess ?
                if (info .ne. 0) then
                    call utmess('F', 'HHO1_4')
                end if
!
! ---- Solve piKF = invM^-1 * piKF
                info = 0
                b_n = to_blas_int(fbs)
                b_nrhs = to_blas_int(cbs)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                b_ldb = to_blas_int(MSIZE_FACE_SCAL)
                call dpotrs('U', b_n, b_nrhs, invM, b_lda, &
                            piKF, b_ldb, info)
!
! --- -Success ?
! ---- Verif strange bug if info neq 0 in entry
                if (info .ne. 0) then
                    call utmess('F', 'HHO1_4')
                end if
            end if
!
            stab%m(fromFace:toFace, cell_offset:total_dofs) = (-invH)*traceMat(1:fbs, 1:cbs)
            stab%m(fromFace:toFace, fromFace:toFace) = invH*faceMass%m(1:fbs, 1:fbs)
!
            if (.not. faceMass%isIdentity) then
!
! ---- Compute stab += invH * traceMat**T * piKF
                b_ldc = to_blas_int(MSIZE_CELL_SCAL)
                b_ldb = to_blas_int(MSIZE_FACE_SCAL)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                b_m = to_blas_int(cbs)
                b_n = to_blas_int(cbs)
                b_k = to_blas_int(fbs)
                call dgemm('T', 'N', b_m, b_n, b_k, &
                           invH, traceMat, b_lda, piKF, b_ldb, &
                           1.d0, S_TT, b_ldc)
!
            else
!
! ---- Compute stab += invH * traceMat**T * traceMat since faceMass = Id
                b_ldc = to_blas_int(MSIZE_CELL_SCAL)
                b_lda = to_blas_int(MSIZE_FACE_SCAL)
                b_n = to_blas_int(cbs)
                b_k = to_blas_int(fbs)
                call dsyrk('U', 'T', b_n, b_k, &
                           invH, traceMat, b_lda, &
                           1.d0, S_TT, b_ldc)
            end if
        end do
!
        stab%m(cell_offset:total_dofs, cell_offset:total_dofs) = S_TT(1:cbs, 1:cbs)
!
        call stab%copySymU()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hdgStabVec(hhoCell, hhoData, stab)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(out) :: stab
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the hdg stabilization of a vectorial function
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   Out stab        : matrix of stabilization (lhs member for vectorial laplacian problem)
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_matrix) :: stab_scal
! --------------------------------------------------------------------------------------------------
!
! -- compute scalar stabilization
        call hdgStabScal(hhoCell, hhoData, stab_scal)
!
! -- copy the scalar stabilization in the vectorial stabilization
        call MatScal2Vec(hhoCell, hhoData, stab_scal, stab)
        call stab_scal%free()
!
    end subroutine
!
end module
