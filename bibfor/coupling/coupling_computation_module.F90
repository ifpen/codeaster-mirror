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
!
module coupling_computation_module
!
    use coupling_quadrature_module, only: cplGetQuadOrderFEM
    use coupling_type
    use FE_Basis_module
    use FE_quadrature_module
    use FE_topo_module
    use HHO_basis_module
    use HHO_geometry_module
    use HHO_init_module
    use HHO_massmat_module
    use HHO_matrix_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
    use HHO_utils_module, only: hhoCopySymPartMat
!
    implicit none
!
    private
!
#include "asterf_types.h"
#include "asterfort/apnorm.h"
#include "asterfort/assert.h"
#include "asterfort/coupling_type.h"
#include "asterfort/HHO_basis_module.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/matrHooke3d.h"
#include "asterfort/matrHookePlaneStrain.h"
#include "asterfort/reereg.h"
#include "asterfort/separ_RI_elas_3D.h"
#include "asterfort/jevech.h"
#include "asterfort/utmess.h"
#include "asterfort/readMatrix.h"
#include "blas/dger.h"
#include "blas/dposv.h"
#include "blas/dsyr.h"
#include "FE_basis_module.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! Coupling methods - computation
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: cplCmpPenaFEMHHOMatScal, cplCmpPenaFEMHHOMatVec
    public :: cplCmpPenaFEMFEMMatScal, cplCmpPenaFEMFEMMatVec
    public :: cplCmpLagrFEMFEMMatScal, cplCmpLagrFEMFEMMatVec
    public :: cplCmpStressFEMHHOMat, cplCmpGradLiftFEMHHOMatVec
    private :: cplEvalElasParaFace, cplEvalElasParaCell
    private :: cplEvalMatrHB, cplEvalMatrHBn
    private :: cplTensorizeFEMHHOMat
    private :: cplCmpProjFEM2HHOFaceScalRead
    public :: cplCmpProjFEM2HHOFaceResScal
    public :: cplCmpProjFEM2HHOFaceResVec
!
! -- Section 6 projected cG unknowns.
!    .true.  => use Pi_F^k(u_CG) - u_F in the penalty term.
!    .false. => recover the original local penalty term u_CG - u_F.
    aster_logical, parameter :: USE_L2_HHO_PROJ = ASTER_TRUE
!
!
contains
!
    subroutine cplTensorizeFEMHHOMat(matScal, ndim, nFEScal, nHHOScal, matVec)
!
        implicit none
!
        type(HHO_matrix), intent(in) :: matScal
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: nFEScal
        integer(kind=8), intent(in) :: nHHOScal
        type(HHO_matrix), intent(out) :: matVec
!
!===================================================================================================
! Tensorize a scalar FEM/HHO matrix.
!
! Scalar ordering:
! [CG node 1, ..., CG node nFE,
! HHO basis 1, ..., HHO basis nHHO]
!
! Vector ordering:
! CG:
! [u1_x, u1_y, (u1_z), u2_x, u2_y, (u2_z), ...]
!
! HHO:
! [phi1_x, ..., phiN_x,
! phi1_y, ..., phiN_y,
! (phi1_z, ..., phiN_z)]
!
! The resulting operator is block diagonal with respect to the physical
! components:
!
! A_vec = I_ndim tensor A_scal
!
! taking into account the different CG and HHO storage conventions.
!===================================================================================================
!
        integer(kind=8) :: nScal, nVec
        integer(kind=8) :: idim, iScal, jScal
        integer(kind=8) :: iVec, jVec
!
        nScal = nFEScal+nHHOScal
        nVec = ndim*nScal
!
        ASSERT(ndim >= 1)
        ASSERT(matScal%nrows == nScal)
        ASSERT(matScal%ncols == nScal)
!
        call matVec%initialize(nVec, nVec, 0.d0)
!
        do idim = 1, ndim
        do jScal = 1, nScal
!
! -------- Map the scalar column to the vector column.
            if (jScal <= nFEScal) then
! ----------- CG ordering: node first, then component.
                jVec = ndim*(jScal-1)+idim
            else
! ----------- HHO ordering: component first, then polynomial basis.
                jVec = ndim*nFEScal &
                       +(idim-1)*nHHOScal &
                       +(jScal-nFEScal)
            end if
!
            do iScal = 1, nScal
!
! ------------ Map the scalar row to the vector row.
                if (iScal <= nFEScal) then
                    iVec = ndim*(iScal-1)+idim
                else
                    iVec = ndim*nFEScal &
                           +(idim-1)*nHHOScal &
                           +(iScal-nFEScal)
                end if
!
                matVec%m(iVec, jVec) = matScal%m(iScal, jScal)
            end do
        end do
        end do
!
    end subroutine cplTensorizeFEMHHOMat
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplCmpPenaFEMHHOMatScal(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                       coef_pena, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuadMa
        real(kind=8), intent(in) :: coef_pena
        type(HHO_matrix), intent(out) :: lhs
!
!===================================================================================================
!    FEM/HHO - Compute penalty matrix.
!
!    Standard term:
!        gamma/hF * (u_CG - u_F, v_CG - v_F)_S
!
!    Projected term, section 6 interpretation for the aster pair-by-pair assembly:
!        - Pi_F^k is defined on the whole HHO face F.
!        - The elementary contribution is integrated only on the current intersection S.
!
!    Hence, for each current FEM basis phi_i, we first compute the coefficients
!
!        c_i = Pi_F^k(phi_i)
!
!    using a full HHO-face quadrature on F. Then we assemble on the coupled quadrature S:
!
!        gamma/hF * (Pi_F^k(phi_i) - psi_a, Pi_F^k(phi_j) - psi_b)_S.
!===================================================================================================
!
        type(HHO_basis_face) :: hhoBasisFace
        type(FE_Basis) :: FEBasis
        real(kind=8), dimension(MAX_BS_CG) :: BSEval, BSEvalProj
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: basisScalEval
        real(kind=8), dimension(MSIZE_FACE_SCAL, MAX_BS_CG) :: projFEM2HHO
        real(kind=8), dimension(MAX_BS_CG+MSIZE_FACE_SCAL) :: BS
        real(kind=8) :: coeff
        integer(kind=8) :: fbs, nbDoFsFE, nbDoFs, ipg, iBasis
        blas_int, parameter :: b_one = to_blas_int(1)
        blas_int :: b_lda, b_n

!
! -- init face basis
        call hhoBasisFace%initialize(hhoFaceMa)
        call FEBasis%initFace(FEFaceSl)
!
! -- number of dofs
        call hhoTherFaceDofs(hhoFaceMa, hhoDataMa, fbs)
        nbDoFsFE = FEBasis%size
        nbDoFs = nbDoFsFE+fbs
!
        call lhs%initialize(nbDoFs, nbDoFs, 0.d0)
        b_n = to_blas_int(lhs%nrows)
        b_lda = to_blas_int(lhs%max_nrows)
!
        if (USE_L2_HHO_PROJ) then
!
! ----- Compute Pi_F^k(phi_i^CG) on the intersection and not whole HHO face F.
! Projected term on the current intersection S:
! - the HHO polynomial space is restricted from the complete face F to S;
! - the L2 projection is computed using the mass matrix restricted to S;
! - the resulting operator depends only on the CG dofs of the current intersection.
!
            call cplCmpProjFEM2HHOFaceScalRead(FEFaceSl, hhoFaceMa, hhoDataMa, &
                                               projFEM2HHO)
        end if
!
! ----- Assemble the elementary contribution only on the current intersection S.
        do ipg = 1, FEQuadSl%nbQuadPoints
!
! --------- HHO face basis on the current intersection point.
            call hhoBasisFace%BSEval(hhoQuadMa%points(1:3, ipg), 0, &
                                     hhoDataMa%face_degree(), basisScalEval)
!
            if (USE_L2_HHO_PROJ) then
!
! --------- Evaluate projected FEM basis functions at the same HHO point:
!           Pi_F^k(phi_i^CG)(x_q) = sum_a c_i(a) psi_a(x_q).
                do iBasis = 1, nbDoFsFE
                    BSEvalProj(iBasis) = dot_product(projFEM2HHO(1:fbs, iBasis), &
                                                     basisScalEval(1:fbs))
                end do
!
! --------- Assemble [Pi_F^k(u_CG), -u_F] on S.
                BS(1:nbDoFsFE) = BSEvalProj(1:nbDoFsFE)
            else
!
! --------- Compute u_CG.
                BSEval = FEBasis%func(FEQuadSl%points_param(1:3, ipg))
! --------- Assemble [u_CG, -u_F].
                BS(1:nbDoFsFE) = BSEval(1:nbDoFsFE)
            end if
! ---------- [-uF] part
            BS(nbDoFsFE+1:nbDoFs) = -basisScalEval(1:fbs)
!
            coeff = FEQuadSl%weights(ipg)*coef_pena/hhoFaceMa%diameter
            call dsyr('U', b_n, coeff, BS, b_one, lhs%m, b_lda)
        end do
!
        call lhs%copySymU()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine cplCmpProjFEM2HHOFaceScalRead(FEFaceSl, hhoFaceMa, hhoDataMa, &
                                             projFEM2HHO)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        real(kind=8), intent(out) :: projFEM2HHO(MSIZE_FACE_SCAL, MAX_BS_CG)
!
!===================================================================================================
!    Read the L2 projection coefficients Pi_F^k(phi_i^CG) on the intersection and
!    not the whole HHO face F. It could be computed directly since it is local.
!===================================================================================================
!
        integer(kind=8) :: fbs_cmp
!
        call hhoTherFaceDofs(hhoFaceMa, hhoDataMa, fbs_cmp)
!
        projFEM2HHO = 0.d0
        call readMatrix("PMPHHOR", fbs_cmp, FEFaceSl%nbnodes, ASTER_FALSE, projFEM2HHO)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine cplCmpProjFEM2HHOFaceResScal(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                            projFEM2HHO)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuadMa
        real(kind=8), intent(out) :: projFEM2HHO(MSIZE_FACE_SCAL, MAX_BS_CG)
!
!===================================================================================================
!    Compute the local L2 projection coefficients Pi_F^k(phi_i^CG) on the intersection
!    and not on the whole HHO face F. Since, it is not local.
!===================================================================================================
!
        type(FE_Basis) :: FEBasis
        type(HHO_basis_face) :: hhoBasisFace
        real(kind=8) :: BSFEEval(MAX_BS_CG), BSHHOEval(MSIZE_FACE_SCAL)
        real(kind=8) :: massMat(MSIZE_FACE_SCAL, MSIZE_FACE_SCAL)
        integer(kind=8) :: face_degree, ipg, feOrder, fbs
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n, b_nrhs, b_ldb, info
!
! -- init face basis
        call FEBasis%initFace(FEFaceSl)
        call hhoBasisFace%initialize(hhoFaceMa)
!
! -- number of dofs and projection degree
        call hhoTherFaceDofs(hhoFaceMa, hhoDataMa, fbs)
        face_degree = hhoDataMa%face_degree()
!
! -- HHO face quadrature for the L2 projection on the intersection.
!    The RHS of the projection contains FEM basis functions times HHO basis functions.
!    Therefore the quadrature order must be high enough for both the HHO face
!    mass matrix check and the FEM/HHO product in the RHS.
        feOrder = cplGetQuadOrderFEM(FEFaceSl%typemas)
        ASSERT((face_degree+feOrder+1) <= hhoQuadMa%order)
!
        b_lda = to_blas_int(MSIZE_FACE_SCAL)
        b_m = to_blas_int(fbs)
        b_n = to_blas_int(FEBasis%size)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
!
! -- Compute local RHS of the projection (uF, u_cG) on the intersection part
!
        projFEM2HHO = 0.d0
        massMat = 0.d0
        ASSERT(hhoQuadMa%nbQuadPoints == FEQuadSl%nbQuadPoints)
        do ipg = 1, hhoQuadMa%nbQuadPoints
            BSFEEval = FEBasis%func(FEQuadSl%points_param(1:3, ipg))
            call hhoBasisFace%BSEval(hhoQuadMa%points(1:3, ipg), 0, face_degree, BSHHOEval)
!
            call dger(b_m, b_n, hhoQuadMa%weights(ipg), BSHHOEval, b_incx, BSFEEval, b_incy, &
                      projFEM2HHO, b_lda)

            call dsyr("U", b_m, hhoQuadMa%weights(ipg), BSHHOEval, b_incx, massMat, b_lda)
        end do
!
! -- Solve the system - Local mass matrix too.
!
        b_n = to_blas_int(fbs)
        b_nrhs = to_blas_int(FEBasis%size)
        b_lda = to_blas_int(MSIZE_FACE_SCAL)
        b_ldb = to_blas_int(MSIZE_FACE_SCAL)
        info = 0
        call dposv('U', b_n, b_nrhs, massMat, b_lda, &
                   projFEM2HHO, b_ldb, info)
!
! -- Sucess ?
!
        if (info .ne. 0) then
            call utmess('F', 'HHO1_4')
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine cplCmpProjFEM2HHOFaceResVec(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                           projFEM2HHO)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuadMa
        type(HHO_matrix), intent(out) :: projFEM2HHO
!
!
!===================================================================================================
!    Compute the local L2 projection coefficients Pi_F^k(phi_i^CG) on the intersection
!    and not on the whole HHO face F. Since, it is not local.
!    Vectorial case by tensorization.
!===================================================================================================
!
        real(kind=8) :: rhs_scal(MSIZE_FACE_SCAL, MAX_BS_CG)
        integer(kind=8) :: ndim, fbs_scal, fbs, nFE_scal, nFE
        integer(kind=8) :: i, j, idim, irow, jcol

        call cplCmpProjFEM2HHOFaceResScal(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                          rhs_scal)
!
        ndim = hhoFaceMa%ndim+1
        call hhoTherFaceDofs(hhoFaceMa, hhoDataMa, fbs_scal)
        fbs = fbs_scal*ndim
        nFE_scal = FEFaceSl%nbnodes
        nFE = nFE_scal*ndim
!
        call projFEM2HHO%initialize(fbs, nFE, 0.d0)
!
        irow = 0
        do idim = 1, ndim
            do i = 1, fbs_scal
                irow = irow+1
                jcol = idim
                do j = 1, nFE_scal
                    projFEM2HHO%m(irow, jcol) = rhs_scal(i, j)
                    jcol = jcol+ndim
                end do
            end do
        end do
!
        ASSERT(irow == fbs)
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine cplCmpGradLiftFEMHHOMatVec(FEFaceSl, FECellSl, hhoFaceMa, hhoDataMa, &
                                          FEQuadSl, hhoQuadMa, cplData, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(FE_Cell), intent(in) :: FECellSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuadMa
        type(CouplingData), intent(in) :: cplData
        type(HHO_matrix), intent(out) :: lhs
!
!===================================================================================================
!    FEM/HHO - Elastic energy of the local CG-side gradient lifting.
!
!    On a CG cell T adjacent to the coupling interface, define
!
!        R_T = G_T^{k,CG} - grad(u_CG)
!
!    through
!
!        (R_T, tau)_T = (u_F^HHO - Pi_F^k(u_CG), tau*n_TF)_F.
!
!    This routine computes for small-strain elasticity,
!
!        (C:sym(R_T(u)), sym(R_T(v)))_T.
!
!    The output ordering is the same as for cplCmpPenaFEMHHOMatVec:
!      - CG face dofs are sorted by node, then by component;
!      - HHO face dofs are sorted by component, then by polynomial basis.
!
!    The volume integrals are evaluated with the native HHO cell quadrature,
!    exactly as for the HHO gradient reconstruction. The current implementation
!    therefore assumes a straight-sided triangular or quadrilateral CG cell.
!===================================================================================================
!
        type(HHO_Cell) :: hhoCellSl
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_basis_face) :: hhoBasisFace
        type(HHO_Quadrature) :: hhoQuadCellSl
        type(FE_Basis) :: FEBasisFaceSl, FEBasisCellSl
        type(HHO_massmat_cell) :: hhoMassMatCellSl
        type(HHO_matrix) :: rhs
        real(kind=8), dimension(MSIZE_FACE_SCAL, MAX_BS_CG) :: projFEM2HHO
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: BSCEvalGradSl
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: BSEvalFaceMa
        real(kind=8), dimension(MAX_BS_CG) :: BSEvalFaceSlProj
        real(kind=8), dimension(MAX_BS_CG) :: BSEvalCellSl
        real(kind=8), dimension(6, MAX_BV_CG+MSIZE_FACE_VEC) :: strain, cstrain
        real(kind=8) :: normalSl(3), coeff, weight
        real(kind=8) :: E, nu, lambda, mu, trace
        real(kind=8), parameter :: un_rac2 = 1.d0/sqrt(2.d0)
        integer(kind=8), parameter ::  XX = 1, YY = 2, XY = 3, ZZ = 4, XZ = 5, YZ = 6
        integer(kind=8) :: ndim, faceDegree, gradDegree, gbs
        integer(kind=8) :: fbsCmp, nbDoFsFE
        integer(kind=8) :: nbDoFs, ipg, i_gb, iBasis, iDoF, jDoF
        integer(kind=8) :: colX, colY, rhsXX, rhsYY, rhsXY, nbSigm
        integer(kind=8) :: colZ, rhsZZ, rhsXZ, rhsYZ
        character(len=8) :: typema
        blas_int :: b_n, b_nrhs, b_lda, b_ldb, info
!
        ndim = FECellSl%ndim
        ASSERT(hhoFaceMa%ndim+1 == ndim)
        if (ndim == 2) then
            nbSigm = 3
        else
            nbSigm = 6
        end if
!
        faceDegree = hhoDataMa%face_degree()
        gradDegree = hhoDataMa%grad_degree()
        ASSERT(faceDegree == gradDegree)
        ASSERT(cplGetQuadOrderFEM(FEFaceSl%typemas) <= faceDegree+1)
!
        select case (FECellSl%typemas)
        case ("TR3", "TR6", "TR7", "TR1")
            typema = "TRIA7"
        case ("QU4", "QU8", "QU9", "Q12")
            typema = "QUAD9"
        case ("HE8", "H20", "H27")
            typema = "HEXA27"
        case ("TE4", "T10", "T20")
            typema = "TETRA15"
        case default
            ASSERT(ASTER_FALSE)
        end select
!
        call hhoCellInit(typema, FECellSl%coorno, hhoCellSl)
        call hhoQuadCellSl%getQuadCell(hhoCellSl, 2*gradDegree, param=ASTER_TRUE)
!
        call FEBasisFaceSl%initFace(FEFaceSl)
        call FEBasisCellSl%initCell(FECellSl)
!
! ----- Use inertial basis since orthonormal does not exist for FEM side
        call hhoBasisCell%initialize(hhoCellSl, BASIS_INERTIAL)
        gbs = hhoBasisCell%BSSize(0, gradDegree)
        ASSERT(gbs <= MSIZE_CELL_SCAL)
!
! ----- Compute mass matrix on FEM side
        call hhoMassMatCellSl%compute(hhoCellSl, 0, gradDegree, &
                                      hhoQuadCellSl, BASIS_INERTIAL)
!
        call hhoBasisFace%initialize(hhoFaceMa)
        call hhoTherFaceDofs(hhoFaceMa, hhoDataMa, fbsCmp)
        call cplCmpProjFEM2HHOFaceScalRead(FEFaceSl, hhoFaceMa, hhoDataMa, &
                                           projFEM2HHO)
        nbDoFsFE = FEBasisFaceSl%size
        nbDoFs = ndim*(nbDoFsFE+fbsCmp)
!
        call rhs%initialize(gbs, nbSigm*nbDoFs, 0.d0)
!
        do ipg = 1, FEQuadSl%nbQuadPoints
            call hhoBasisCell%BSEval(FEQuadSl%points(1:3, ipg), 0, gradDegree, &
                                     BSCEvalGradSl)
            call hhoBasisFace%BSEval(hhoQuadMa%points(1:3, ipg), 0, faceDegree, &
                                     BSEvalFaceMa)
!
            BSEvalFaceSlProj = 0.d0
            do iBasis = 1, nbDoFsFE
                BSEvalFaceSlProj(iBasis) = dot_product(projFEM2HHO(1:fbsCmp, iBasis), &
                                                       BSEvalFaceMa(1:fbsCmp))
            end do
!
            normalSl = FEFaceSl%normal(FEQuadSl%points_param(1:3, ipg))
            weight = FEQuadSl%weights(ipg)
!
! -------- CG contribution: -Pi_F^k(u_CG).
!
! CG dof ordering:
!
! 2D: ux_1, uy_1, ux_2, uy_2, ...
! 3D: ux_1, uy_1, uz_1, ux_2, uy_2, uz_2, ...
!
            do iBasis = 1, nbDoFsFE

                colX = ndim*(iBasis-1)+1
                colY = ndim*(iBasis-1)+2

                if (ndim == 3) then
                    colZ = ndim*(iBasis-1)+3
                end if

                do i_gb = 1, gbs

                    coeff = -weight &
                            *BSCEvalGradSl(i_gb) &
                            *BSEvalFaceSlProj(iBasis)

! ---------------------------------------------------
! Contribution of the x-component CG degree of freedom
! ---------------------------------------------------

                    rhsXX = (XX-1)*nbDoFs+colX

                    rhs%m(i_gb, rhsXX) = rhs%m(i_gb, rhsXX) &
                                         +coeff*normalSl(1)

                    rhsXY = (XY-1)*nbDoFs+colX

                    rhs%m(i_gb, rhsXY) = rhs%m(i_gb, rhsXY) &
                                         +coeff*normalSl(2)*un_rac2

                    if (ndim == 3) then
                        rhsXZ = (XZ-1)*nbDoFs+colX

                        rhs%m(i_gb, rhsXZ) = rhs%m(i_gb, rhsXZ) &
                                             +coeff*normalSl(3)*un_rac2
                    end if

! ---------------------------------------------------
! Contribution of the y-component CG degree of freedom
! ---------------------------------------------------

                    rhsYY = (YY-1)*nbDoFs+colY

                    rhs%m(i_gb, rhsYY) = rhs%m(i_gb, rhsYY) &
                                         +coeff*normalSl(2)

                    rhsXY = (XY-1)*nbDoFs+colY

                    rhs%m(i_gb, rhsXY) = rhs%m(i_gb, rhsXY) &
                                         +coeff*normalSl(1)*un_rac2

                    if (ndim == 3) then
                        rhsYZ = (YZ-1)*nbDoFs+colY

                        rhs%m(i_gb, rhsYZ) = rhs%m(i_gb, rhsYZ) &
                                             +coeff*normalSl(3)*un_rac2
                    end if

! ---------------------------------------------------
! Contribution of the z-component CG degree of freedom
! ---------------------------------------------------

                    if (ndim == 3) then
                        rhsZZ = (ZZ-1)*nbDoFs+colZ

                        rhs%m(i_gb, rhsZZ) = rhs%m(i_gb, rhsZZ) &
                                             +coeff*normalSl(3)

                        rhsXZ = (XZ-1)*nbDoFs+colZ

                        rhs%m(i_gb, rhsXZ) = rhs%m(i_gb, rhsXZ) &
                                             +coeff*normalSl(1)*un_rac2

                        rhsYZ = (YZ-1)*nbDoFs+colZ

                        rhs%m(i_gb, rhsYZ) = rhs%m(i_gb, rhsYZ) &
                                             +coeff*normalSl(2)*un_rac2
                    end if

                end do
            end do

! -------- HHO contribution: +u_F^HHO.
!
! HHO face dof ordering:
!
! 2D:
! ux_phi1, ..., ux_phiN,
! uy_phi1, ..., uy_phiN
!
! 3D:
! ux_phi1, ..., ux_phiN,
! uy_phi1, ..., uy_phiN,
! uz_phi1, ..., uz_phiN
!
            do iBasis = 1, fbsCmp

! First HHO dof comes after all CG dofs.
                colX = ndim*nbDoFsFE+iBasis
                colY = ndim*nbDoFsFE+fbsCmp+iBasis

                if (ndim == 3) then
                    colZ = ndim*nbDoFsFE+2*fbsCmp+iBasis
                end if

                do i_gb = 1, gbs

                    coeff = weight &
                            *BSCEvalGradSl(i_gb) &
                            *BSEvalFaceMa(iBasis)

! ---------------------------------------------------
! Contribution of the x-component HHO degree of freedom
! ---------------------------------------------------

                    rhsXX = (XX-1)*nbDoFs+colX

                    rhs%m(i_gb, rhsXX) = rhs%m(i_gb, rhsXX) &
                                         +coeff*normalSl(1)

                    rhsXY = (XY-1)*nbDoFs+colX

                    rhs%m(i_gb, rhsXY) = rhs%m(i_gb, rhsXY) &
                                         +coeff*normalSl(2)*un_rac2

                    if (ndim == 3) then
                        rhsXZ = (XZ-1)*nbDoFs+colX

                        rhs%m(i_gb, rhsXZ) = rhs%m(i_gb, rhsXZ) &
                                             +coeff*normalSl(3)*un_rac2
                    end if

! ---------------------------------------------------
! Contribution of the y-component HHO degree of freedom
! ---------------------------------------------------

                    rhsYY = (YY-1)*nbDoFs+colY

                    rhs%m(i_gb, rhsYY) = rhs%m(i_gb, rhsYY) &
                                         +coeff*normalSl(2)

                    rhsXY = (XY-1)*nbDoFs+colY

                    rhs%m(i_gb, rhsXY) = rhs%m(i_gb, rhsXY) &
                                         +coeff*normalSl(1)*un_rac2

                    if (ndim == 3) then
                        rhsYZ = (YZ-1)*nbDoFs+colY

                        rhs%m(i_gb, rhsYZ) = rhs%m(i_gb, rhsYZ) &
                                             +coeff*normalSl(3)*un_rac2
                    end if

! ---------------------------------------------------
! Contribution of the z-component HHO degree of freedom
! ---------------------------------------------------

                    if (ndim == 3) then
                        rhsZZ = (ZZ-1)*nbDoFs+colZ

                        rhs%m(i_gb, rhsZZ) = rhs%m(i_gb, rhsZZ) &
                                             +coeff*normalSl(3)

                        rhsXZ = (XZ-1)*nbDoFs+colZ

                        rhs%m(i_gb, rhsXZ) = rhs%m(i_gb, rhsXZ) &
                                             +coeff*normalSl(1)*un_rac2

                        rhsYZ = (YZ-1)*nbDoFs+colZ

                        rhs%m(i_gb, rhsYZ) = rhs%m(i_gb, rhsYZ) &
                                             +coeff*normalSl(2)*un_rac2
                    end if

                end do
            end do
        end do
!
! --- Compute gradient reconstruction
        if (.not. hhoMassMatCellSl%isIdentity) then
            info = 0
            b_n = to_blas_int(gbs)
            b_nrhs = to_blas_int(nbSigm*nbDoFs)
            b_lda = to_blas_int(hhoMassMatCellSl%max_nrows)
            b_ldb = to_blas_int(rhs%max_nrows)
            call dposv('U', b_n, b_nrhs, hhoMassMatCellSl%m, b_lda, rhs%m, b_ldb, info)
            if (info /= 0) then
                call utmess('F', 'HHO1_4')
            end if
        end if
!
!
        strain = 0.d0
        cstrain = 0.d0
        call lhs%initialize(nbDoFs, nbDoFs, 0.d0)
!
        do ipg = 1, hhoQuadCellSl%nbQuadPoints
            call hhoBasisCell%BSEval(hhoQuadCellSl%points(1:3, ipg), 0, gradDegree, &
                                     BSCEvalGradSl)
            strain = 0.d0
            do iDoF = 1, nbDoFs
                do i_gb = 1, gbs
                    strain(XX, iDoF) = strain(XX, iDoF)+ &
                                       BSCEvalGradSl(i_gb)*rhs%m(i_gb, (XX-1)*nbDoFs+iDoF)
                    strain(YY, iDoF) = strain(YY, iDoF)+ &
                                       BSCEvalGradSl(i_gb)*rhs%m(i_gb, (YY-1)*nbDoFs+iDoF)
                    strain(XY, iDoF) = strain(XY, iDoF)+ &
                                       BSCEvalGradSl(i_gb)*rhs%m(i_gb, (XY-1)*nbDoFs+iDoF)
                    if (ndim == 3) then
                        strain(ZZ, iDoF) = strain(ZZ, iDoF)+ &
                                           BSCEvalGradSl(i_gb)*rhs%m(i_gb, (ZZ-1)*nbDoFs+iDoF)
                        strain(XZ, iDoF) = strain(XZ, iDoF)+ &
                                           BSCEvalGradSl(i_gb)*rhs%m(i_gb, (XZ-1)*nbDoFs+iDoF)
                        strain(YZ, iDoF) = strain(YZ, iDoF)+ &
                                           BSCEvalGradSl(i_gb)*rhs%m(i_gb, (YZ-1)*nbDoFs+iDoF)
                    end if
                end do
            end do
!
! ------- Compute elastic parameters
!
            BSEvalCellSl = FEBasisCellSl%func(hhoQuadCellSl%points_param(1:3, ipg))
            call cplEvalElasParaCell(cplData, FEBasisCellSl%size, BSEvalCellSl, E, nu, lambda, mu)
!
! ------- Compute elastic rigidity
!
            do iDoF = 1, nbDoFs
                trace = strain(XX, iDoF)+strain(YY, iDoF)+strain(ZZ, iDoF)
                cstrain(1:nbSigm, iDoF) = 2.d0*mu*strain(1:nbSigm, iDoF)
                cstrain(XX, iDOF) = cstrain(XX, iDOF)+lambda*trace
                cstrain(YY, iDOF) = cstrain(YY, iDOF)+lambda*trace
!
                if (ndim == 3) then
                    cstrain(ZZ, iDOF) = cstrain(ZZ, iDOF)+lambda*trace
                end if
            end do
!
            weight = hhoQuadCellSl%weights(ipg)
            do jDoF = 1, nbDoFs
                do iDoF = 1, jDoF
                    lhs%m(iDoF, jDoF) = lhs%m(iDoF, jDoF)+weight* &
                                        dot_product(strain(1:nbSigm, iDoF), cstrain(1:nbSigm, jDoF))
                end do
            end do
        end do
!
        call lhs%copySymU()
        call rhs%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine cplCmpPenaFEMHHOMatVec(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                      coef_pena, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuadMa
        real(kind=8), intent(in) :: coef_pena
        type(HHO_matrix), intent(out) :: lhs
!
!===================================================================================================
!    FEM/HHO - Compute matrix: gamma/hF*(uh-uF, vh-vF)
!              unknowns (uh, uF)
!===================================================================================================
!
        type(HHO_matrix) :: lhs_scal
        integer(kind=8) :: ndim, fbs_scal

        call cplCmpPenaFEMHHOMatScal(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                     coef_pena, lhs_scal)
!
        ndim = hhoFaceMa%ndim+1
        call hhoTherFaceDofs(hhoFaceMa, hhoDataMa, fbs_scal)
!
        call cplTensorizeFEMHHOMat(lhs_scal, ndim, FEFaceSl%nbnodes, fbs_scal, lhs)
!
        call lhs_scal%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplCmpPenaFEMFEMMatScal(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, &
                                       coef_pena, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        real(kind=8), intent(in) :: coef_pena
        real(kind=8), intent(out) :: lhs(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
!
!===================================================================================================
!    FEM/FEM - Compute matrix: gamma/hF*(uh1-uh2, vh1-vh2)
!              unknowns (uh1, uh2)
!===================================================================================================
!
        type(FE_Basis) :: FEBasisSl, FEBasisMa
        real(kind=8), dimension(MAX_BS_CG) :: BSEvalSl, BSEvalMa
        real(kind=8), dimension(MAX_BS_CG+MAX_BS_CG) :: BS
        real(kind=8) :: coeff, hF
        integer(kind=8) :: nbDoFsFEMa, nbDoFsFESl, nbDoFs, ipg
        blas_int, parameter :: b_one = to_blas_int(1)
        blas_int :: b_lda, b_n
!
! -- init face basis
        call FEBasisSl%initFace(FEFaceSl)
        call FEBasisMa%initFace(FEFaceMa)
!
! -- number of dofs
        nbDoFsFESl = FEBasisSl%size
        nbDoFsFEMa = FEBasisMa%size
        nbDoFs = nbDoFsFESl+nbDoFsFEMa
!
        b_n = to_blas_int(nbDoFs)
        b_lda = to_blas_int(MSIZE_CPL_PENA_FEM)
!
        hF = FEFaceSl%diameter()
!
        lhs = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuadSl%nbQuadPoints
! ----- Eval face basis function at the quadrature point
! -------- Compute uh1
            BSEvalSl = FEBasisSl%func(FEQuadSl%points_param(1:3, ipg))
            BSEvalMa = FEBasisMa%func(FEQuadMa%points_param(1:3, ipg))
!
! -------- Assemble [uh1, -uh2]
            BS(1:nbDoFsFESl) = BSEvalSl(1:nbDoFsFESl)
            BS(nbDoFsFESl+1:nbDoFs) = -BSEvalMa(1:nbDoFsFEMa)
!
! --------  Eval

            coeff = FEQuadSl%weights(ipg)*coef_pena/hF
            call dsyr('U', b_n, coeff, BS, b_one, lhs, b_lda)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine cplCmpPenaFEMFEMMatVec(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, coef_pena, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        real(kind=8), intent(in) :: coef_pena
        real(kind=8), intent(out) :: lhs(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
!
!===================================================================================================
!    FEM/FEM - Compute matrix: gamma/hF*(uh1-uh2, vh1-vh2)
!              unknowns (uh1, uh2)
!===================================================================================================
!
        real(kind=8) :: lhs_scal(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
        integer(kind=8) :: ndim, ndofsSl_scal, ndofsMa_scal, ndofs_scal
        integer(kind=8) :: idim, i, j, deca_row, deca_col, ndofs

        call cplCmpPenaFEMFEMMatScal(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, &
                                     coef_pena, lhs_scal)
!
        ndim = FEFaceSl%ndim+1
        ndofsSl_scal = FEFaceSl%nbnodes
        ndofsMa_scal = FEFaceMa%nbnodes
        ndofs_scal = ndofsSl_scal+ndofsMa_scal
        ndofs = ndim*ndofs_scal
!
        lhs = 0.d0
!
        do idim = 1, ndim
            deca_col = idim
            do j = 1, ndofs_scal
                deca_row = idim
                do i = 1, j
                    lhs(deca_row, deca_col) = lhs_scal(i, j)
                    deca_row = deca_row+ndim
                end do
                deca_col = deca_col+ndim
            end do
        end do
!
! ----- Copy the lower part
!
        call hhoCopySymPartMat('U', lhs(1:ndofs, 1:ndofs))
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplCmpLagrFEMFEMMatScal(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa, FEFaceLagSl
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        real(kind=8), intent(out) :: lhs(MSIZE_LAGR_FEM, MSIZE_CPL_PENA_FEM)
!
!===================================================================================================
!    FEM/FEM - Compute matrix: (mu^s, u^s-u^m)
!              unknowns ((u^s, lag^s), u^m)
!===================================================================================================
!
        type(FE_Basis) :: FEBasisSl, FEBasisMa, FEBasisLagSl
        real(kind=8), dimension(MAX_BS_CG) :: BSEvalSl, BSEvalMa, BSEvalLagSl
        real(kind=8), dimension(MAX_BS_CG+MAX_BS_CG) :: BS
        real(kind=8) :: coeff
        integer(kind=8) :: nbDoFsFEMa, nbDoFsFESl, nbDoFs, ipg, nbDoFsFELagSl
        blas_int, parameter :: b_one = to_blas_int(1)
        blas_int :: b_lda, b_n, b_m
!
! -- init face basis
        call FEBasisSl%initFace(FEFaceSl)
        call FEBasisMa%initFace(FEFaceMa)
        call FEBasisLagSl%initFace(FEFaceLagSl)
!
! -- number of dofs
        nbDoFsFELagSl = FEBasisLagSl%size
        nbDoFsFESl = FEBasisSl%size
        nbDoFsFEMa = FEBasisMa%size
        nbDoFs = nbDoFsFESl+nbDoFsFEMa
!
        b_n = to_blas_int(nbDoFs)
        b_m = to_blas_int(nbDoFsFELagSl)
        b_lda = to_blas_int(MSIZE_LAGR_FEM)
!
        lhs = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuadSl%nbQuadPoints
! ----- Eval face basis function at the quadrature point
! -------- Compute uh1
            BSEvalSl = FEBasisSl%func(FEQuadSl%points_param(1:3, ipg))
            BSEvalMa = FEBasisMa%func(FEQuadMa%points_param(1:3, ipg))
            BSEvalLagSl = FEBasisLagSl%func(FEQuadSl%points_param(1:3, ipg))
!
! -------- Assemble [u^s, -u^m]
            BS(1:nbDoFsFESl) = BSEvalSl(1:nbDoFsFESl)
            BS(nbDoFsFESl+1:nbDoFs) = -BSEvalMa(1:nbDoFsFEMa)
!
! --------  Eval

            coeff = FEQuadSl%weights(ipg)
            call dger(b_m, b_n, coeff, BSEvalLagSl, b_one, BS, b_one, lhs, b_lda)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine cplCmpLagrFEMFEMMatVec(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa, FEFaceLagSl
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        real(kind=8), intent(out) :: lhs(MSIZE_LAGR_FEM, MSIZE_CPL_PENA_FEM)
!
!===================================================================================================
!    FEM/FEM - Compute matrix: (mu^s, u^s-u^m)
!              unknowns ((u^s, lag^s), u^m)
!===================================================================================================
!
        real(kind=8) :: lhs_scal(MSIZE_LAGR_FEM, MSIZE_CPL_PENA_FEM)
        integer(kind=8) :: ndim, ndofsSl_scal, ndofsMa_scal, ndofs_scal, ndofsLagSl_scal
        integer(kind=8) :: idim, i, j, deca_row, deca_col, ndofs

        call cplCmpLagrFEMFEMMatScal(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, &
                                     lhs_scal)
!
        ndim = FEFaceSl%ndim+1
        ndofsLagSl_scal = FEFaceLagSl%nbnodes
        ndofsSl_scal = FEFaceSl%nbnodes
        ndofsMa_scal = FEFaceMa%nbnodes
        ndofs_scal = ndofsSl_scal+ndofsMa_scal
        ndofs = ndim*(ndofs_scal+ndofsLagSl_scal)
!
        lhs = 0.d0
!
        do idim = 1, ndim
            deca_col = idim
            do j = 1, ndofs_scal
                deca_row = idim
                do i = 1, ndofsLagSl_scal
                    lhs(deca_row, deca_col) = lhs_scal(i, j)
                    deca_row = deca_row+ndim
                end do
                deca_col = deca_col+ndim
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplCmpStressFEMHHOMat(FEFaceSl, FECellSl, hhoFaceMa, hhoDataMa, &
                                     FEQuadSl, hhoQuadMa, cplMap, cplData, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(FE_Cell), intent(in) :: FECellSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuadMa
        type(CouplingMap), intent(in) :: cplMap
        type(CouplingData), intent(in) :: cplData
        type(HHO_matrix), intent(out) :: lhs
!
!===================================================================================================
!    FEM/HHO - Compute matrix: (sigma(u^s).n^s, v^s-v^m)
!              unknowns (u^s, u^m)
!===================================================================================================
!
        type(HHO_basis_face) :: hhoBasisFace
        type(FE_Basis) :: FEBasisFaceSl, FEBasisCellSl
        real(kind=8), dimension(MAX_BS_CG) :: BSEvalFaceSl, BSEvalFaceSlProj
        real(kind=8), dimension(3, MAX_BS_CG) :: BSGEvalCellSl
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: BSEvalFaceMa
        real(kind=8), dimension(MSIZE_FACE_SCAL, MAX_BS_CG) :: projFEM2HHO
        real(kind=8) :: coor_qp_vo(3), normalSl(3), E, nu, lambda, mu
        real(kind=8) :: stress_n(3, 3), weight, matHB(6, 3)
        integer(kind=8) :: nbDoFsFace, nbDoFsFEFaceSl, nbDoFsFECellSl
        integer(kind=8) :: ndim, ipg, fbs, fbs_cmp, idim, iBasis, iRow, iret, jCol
        integer(kind=8) :: jdim, jBasis
!
! -- init face basis
!
        call FEBasisFaceSl%initFace(FEFaceSl)
        call FEBasisCellSl%initCell(FECellSl)
        call hhoBasisFace%initialize(hhoFaceMa)
!
! -- number of dofs
!
        ndim = FECellSl%ndim
        call hhoTherFaceDofs(hhoFaceMa, hhoDataMa, fbs_cmp)
        fbs = ndim*fbs_cmp
        nbDoFsFECellSl = ndim*FEBasisCellSl%size
        nbDoFsFEFaceSl = ndim*FEBasisFaceSl%size
        nbDoFsFace = nbDoFsFEFaceSl+fbs
!
        call lhs%initialize(nbDoFsFECellSl, nbDoFsFace, 0.d0)
!
! -- Section 6 projected consistency term.
!    The projection operator Pi_F^k is defined on the complete HHO face F,
!    while the final elementary integration remains local to the current intersection S.
!    Only the CG trace part in (v_CG - v_F) is replaced by Pi_F^k(v_CG).
        if (USE_L2_HHO_PROJ) then
            call cplCmpProjFEM2HHOFaceScalRead(FEFaceSl, hhoFaceMa, hhoDataMa, &
                                               projFEM2HHO)
        end if
!
        BSEvalFaceSlProj = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuadSl%nbQuadPoints
            weight = FEQuadSl%weights(ipg)
!
! ----- Projection of node on volumic slave cell (volumic parametric space)
!
            coor_qp_vo = 0.d0
            call reereg('S', FECellSl%typemas, FECellSl%nbnodes, FECellSl%coorno, &
                        FEQuadSl%points(1:3, ipg), ndim, coor_qp_vo, iret, &
                        ndim_coor_=3)
            ASSERT(iret == 0)
!
! ----- Eval basis function at the quadrature point
! -------- Compute uh and grad(uh)
!
            BSGEvalCellSl = FEBasisCellSl%grad(coor_qp_vo)
            BSEvalFaceSl = FEBasisFaceSl%func(FEQuadSl%points_param(1:3, ipg))
!
! -------- Compute uF
!
            call hhoBasisFace%BSEval(hhoQuadMa%points(1:3, ipg), 0, hhoDataMa%face_degree(), &
                                     BSEvalFaceMa)
!
! -------- For the projected Nitsche term, evaluate Pi_F^k(phi_j^CG) at the
!          current intersection quadrature point. The HHO part remains -v_F.
            if (USE_L2_HHO_PROJ) then
                do jBasis = 1, FEBasisFaceSl%size
                    BSEvalFaceSlProj(jBasis) = dot_product(projFEM2HHO(1:fbs_cmp, jBasis), &
                                                           BSEvalFaceMa(1:fbs_cmp))
                end do
            else
                BSEvalFaceSlProj(1:FEBasisFaceSl%size) = BSEvalFaceSl(1:FEBasisFaceSl%size)
            end if
!
! ------- Compute outward normal
!
            normalSl = FEFaceSl%normal(FEQuadSl%points_param(1:3, ipg))
!
! ------- Compute elastic parameters
!
            call cplEvalElasParaFace(cplData, cplMap, FEBasisFaceSl%size, BSEvalFaceSl, &
                                     E, nu, lambda, mu)
!
! --------  Product (sigma(u^s).n^s, v^s-v^m)
!
            iRow = 0
            do iBasis = 1, FEBasisCellSl%size
! ------------- Compute [H] : {eps}
                call cplEvalMatrHB(BSGEvalCellSl(1:3, iBasis), lambda, mu, matHB)
! ------------- Compute normal component
                call cplEvalMatrHBn(matHB, normalSl, stress_n)
!
                do idim = 1, ndim
                    iRow = iRow+1
!
! ----------------- (sigma(u^s).n^s, v^s) - sorted by node
                    jCol = 0
                    do jBasis = 1, FEBasisFaceSl%size
                        do jdim = 1, ndim
                            jCol = jCol+1
                            lhs%m(iRow, jCol) = lhs%m(iRow, jCol)+ &
                                                weight*stress_n(jdim, idim)*BSEvalFaceSlProj(jBasis)
                        end do
                    end do
!
! ----------------- (sigma(u^s).n^s, -v^m) - sorted by dimension
                    do jdim = 1, ndim
                        do jBasis = 1, fbs_cmp
                            jCol = jCol+1
                            lhs%m(iRow, jCol) = lhs%m(iRow, jCol)- &
                                                weight*stress_n(jdim, idim)*BSEvalFaceMa(jBasis)
                        end do
                    end do
                end do
            end do
        end do
!
        ASSERT(iRow == nbDoFsFECellSl)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplEvalElasParaCell(cplData, nbNodes, BSEvalCellSl, E, nu, lambda, mu)
!
        implicit none
!
        type(CouplingData), intent(in) :: cplData
        integer(kind=8), intent(in) :: nbNodes
        real(kind=8), dimension(MAX_BS_CG), intent(in) :: BSEvalCellSl
        real(kind=8), intent(out) :: E, nu, lambda, mu
!
!===================================================================================================
!    Evaluate elastic parameter by interpolation
!===================================================================================================
!
        integer(kind=8) :: iNode
!
        E = 0.d0
        nu = 0.d0
!
        do iNode = 1, nbNodes
            E = E+cplData%E(iNode)*BSEvalCellSl(iNode)
            nu = nu+cplData%nu(iNode)*BSEvalCellSl(iNode)
        end do
!
        lambda = E*nu/(1.d0+nu)/(1.d0-2.d0*nu)
        mu = 0.5d0*E/(1.d0+nu)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplEvalElasParaFace(cplData, cplMap, nbNodes, BSEvalFaceSl, E, nu, lambda, mu)
!
        implicit none
!
        type(CouplingData), intent(in) :: cplData
        type(CouplingMap), intent(in) :: cplMap
        integer(kind=8), intent(in) :: nbNodes
        real(kind=8), dimension(MAX_BS_CG), intent(in) :: BSEvalFaceSl
        real(kind=8), intent(out) :: E, nu, lambda, mu
!
!===================================================================================================
!    Evaluate elastic parameter by interpolation
!===================================================================================================
!
        integer(kind=8) :: iNode, node
!
        E = 0.d0
        nu = 0.d0
!
        do iNode = 1, nbNodes
            node = cplMap%face2cellSl(iNode)
            E = E+cplData%E(node)*BSEvalFaceSl(iNode)
            nu = nu+cplData%nu(node)*BSEvalFaceSl(iNode)
        end do
!
        lambda = E*nu/(1.d0+nu)/(1.d0-2.d0*nu)
        mu = 0.5d0*E/(1.d0+nu)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplEvalMatrHB(grad, lambda, mu, matHB)
!
        implicit none
!
        real(kind=8), intent(in) :: grad(3), lambda, mu
        real(kind=8), intent(out) :: matHB(6, 3)
!
!===================================================================================================
!    Compute {H} : {Eps} with H the Hooke matrix for one basis function
!    Sigma = 2*mu*eps + lambda*trace(eps)*Id
!    Component [XX, YY, ZZ, XY, XZ, YZ]
!===================================================================================================
!
!
!       sigma_xx(ux, uy, uz)
        matHB(1, 1) = (2.d0*mu+lambda)*grad(1)
        matHB(1, 2) = lambda*grad(2)
        matHB(1, 3) = lambda*grad(3)
!
!       sigma_yy(ux, uy, uz)
        matHB(2, 1) = lambda*grad(1)
        matHB(2, 2) = (2.d0*mu+lambda)*grad(2)
        matHB(2, 3) = lambda*grad(3)
!
!       sigma_zz(ux, uy, uz)
        matHB(3, 1) = lambda*grad(1)
        matHB(3, 2) = lambda*grad(2)
        matHB(3, 3) = (2.d0*mu+lambda)*grad(3)
!
!       sigma_xy(ux, uy, uz)
        matHB(4, 1) = mu*grad(2)
        matHB(4, 2) = mu*grad(1)
        matHB(4, 3) = 0.d0
!
!       sigma_xz(ux, uy, uz)
        matHB(5, 1) = mu*grad(3)
        matHB(5, 2) = 0.d0
        matHB(5, 3) = mu*grad(1)
!
!       sigma_yz(ux, uy, uz)
        matHB(6, 1) = 0.d0
        matHB(6, 2) = mu*grad(3)
        matHB(6, 3) = mu*grad(2)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplEvalMatrHBn(matHB, normal, sigma_n)
!
        implicit none
!
        real(kind=8), intent(in) :: matHB(6, 3), normal(3)
        real(kind=8), intent(out) :: sigma_n(3, 3)
!
!===================================================================================================
!    Compute sigma.n
!===================================================================================================
!
!       sigma_xx(ux).nx+sigma_xy(ux).ny+sigma_xz(ux).nz
        sigma_n(1, 1) = matHB(1, 1)*normal(1)+matHB(4, 1)*normal(2)+matHB(5, 1)*normal(3)
!       sigma_xx(uy).nx+sigma_xy(uy).ny+sigma_xz(uy).nz
        sigma_n(1, 2) = matHB(1, 2)*normal(1)+matHB(4, 2)*normal(2)+matHB(5, 2)*normal(3)
!       sigma_xx(uz).nx+sigma_xy(uz).ny+sigma_xz(uz).nz
        sigma_n(1, 3) = matHB(1, 3)*normal(1)+matHB(4, 3)*normal(2)+matHB(5, 3)*normal(3)
!
!       sigma_yx(ux).nx+sigma_yy(ux).ny+sigma_yz(ux).nz
        sigma_n(2, 1) = matHB(4, 1)*normal(1)+matHB(2, 1)*normal(2)+matHB(6, 1)*normal(3)
!       sigma_yx(uy).nx+sigma_yy(uy).ny+sigma_yz(uy).nz
        sigma_n(2, 2) = matHB(4, 2)*normal(1)+matHB(2, 2)*normal(2)+matHB(6, 2)*normal(3)
!       sigma_yx(uz).nx+sigma_yy(uz).ny+sigma_yz(uz).nz
        sigma_n(2, 3) = matHB(4, 3)*normal(1)+matHB(2, 3)*normal(2)+matHB(6, 3)*normal(3)
!
!       sigma_zx(ux).nx+sigma_zy(ux).ny+sigma_zz(ux).nz
        sigma_n(3, 1) = matHB(5, 1)*normal(1)+matHB(6, 1)*normal(2)+matHB(3, 1)*normal(3)
!       sigma_zx(uy).nx+sigma_zy(uy).ny+sigma_zz(uy).nz
        sigma_n(3, 2) = matHB(5, 2)*normal(1)+matHB(6, 2)*normal(2)+matHB(3, 2)*normal(3)
!       sigma_zx(uz).nx+sigma_zy(uz).ny+sigma_zz(uz).nz
        sigma_n(3, 3) = matHB(5, 3)*normal(1)+matHB(6, 3)*normal(2)+matHB(3, 3)*normal(3)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
end module
