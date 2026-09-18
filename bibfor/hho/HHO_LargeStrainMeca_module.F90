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

module HHO_LargeStrainMeca_module
!
    use Behaviour_module
    use Behaviour_type
    use FE_algebra_module
    use HHO_algebra_module
    use HHO_basis_module
    use HHO_compor_module
    use HHO_eval_module
    use HHO_matrix_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
    use HHO_utils_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/desymt46.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/lagmodtonommod.h"
#include "asterfort/lcdetf.h"
#include "asterfort/matinv.h"
#include "asterfort/nmcomp.h"
#include "asterfort/pk2sig.h"
#include "asterfort/pk2topk1.h"
#include "asterfort/poslog.h"
#include "asterfort/prelog.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - mechanics
!
! Module for large deformations with hho
!
! --------------------------------------------------------------------------------------------------
    public :: hhoLargeStrainLCMeca, hhoCalculF, hhoCalculGreenLagrange
    public :: hhoComputeLhsLarge, hhoComputeRhsLarge, hhoCalculFinger, hhoCalculEpsDef
    public :: hhoAssembleLhsLarge
    private :: hhoComputeAgphi, hhoMatrixToVoigth
    private :: select_behavior, gdeflog, nbsigm_cmp, greenlagr
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoLargeStrainLCMeca(hhoCell, hhoData, hhoQuadCellRigi, hhoCS, gradrec, &
                                    time_prev, time_curr, depl_prev, depl_curr, lhs, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuadCellRigi
        type(HHO_Compor_State), intent(inout) :: hhoCS
        type(HHO_matrix), intent(in) :: gradrec
        real(kind=8), intent(in) :: time_prev
        real(kind=8), intent(in) :: time_curr
        real(kind=8), intent(in) :: depl_prev(MSIZE_TDOFS_VEC)
        real(kind=8), intent(in) :: depl_curr(MSIZE_TDOFS_VEC)
        type(HHO_matrix), intent(inout) :: lhs
        real(kind=8), intent(inout) :: rhs(MSIZE_TDOFS_VEC)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the local contribution for mechanics in large deformations
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   In hhoQuadCellRigi : quadrature rules from the rigidity family
!   InOut hhoCS        : hho compor state
!   In gradrec      : local gradient reconstruction
!   In time_prev    : previous time T-
!   In time_curr    : current time T+
!   In depl_prev    : displacement at T-
!   In depl_curr    : displacement at T+
!   Out lhs         : local contribution (lhs)
!   Out rhs         : local contribution (rhs)
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), parameter :: ksp = 1
        type(HHO_basis_cell) :: hhoBasisCell
        real(kind=8), dimension(MSIZE_CELL_MAT) :: bT, G_prev_coeff, G_curr_coeff
        real(kind=8) :: module_tang(3, 3, 3, 3), G_prev(3, 3), G_curr(3, 3)
        real(kind=8) :: F_prev(3, 3), F_curr(3, 3), Pk1_curr(3, 3)
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
        type(HHO_matrix) :: AT
        real(kind=8) :: jac_prev, jac_curr, coorpg(3), weight
        integer(kind=8) :: cbs, fbs, total_dofs, faces_dofs, gbs, kpg, gbs_cmp, gbs_sym
        integer(kind=8) :: cod(MSIZE_QP_CELL), nbsig, cbs_cmp, gbs_axis
        aster_logical :: l_gdeflog, l_green_lagr, l_lhs, l_rhs
!
! --------------------------------------------------------------------------------------------------
!
        cod = 0
!
! ------ number of dofs
!
        call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, &
                           gbs, gbs_sym, gbs_axis)
        faces_dofs = total_dofs-cbs
        gbs_cmp = (gbs-gbs_axis)/(hhoCell%ndim*hhoCell%ndim)
        cbs_cmp = cbs/hhoCell%ndim
!
        nbsig = nbsigm_cmp(hhoCell%ndim)
        ASSERT(nbsig == hhoCS%nbsigm)
        bT = 0.d0
        G_prev_coeff = 0.d0
        G_curr_coeff = 0.d0
        !print*, "GT", hhoNorm2Mat(gradrec(1:gbs,1:total_dofs))

! ----- Type of behavior
        call select_behavior(hhoCS%compor, l_gdeflog, l_green_lagr)

! ----- Prepare external state variables (geometry)
        call behaviourPrepESVAGeomHHO(hhoCell, hhoQuadCellRigi, hhoCS%BEHInteg)

! ----- Vector and/or matrix
        l_lhs = L_MATR(hhoCS%option)
        l_rhs = L_VECT(hhoCS%option)
!
        if (hhoCS%c_plan) then
            ASSERT(ASTER_FALSE)
        end if
!
        if (l_lhs) then
            call AT%initialize(gbs, gbs, 0.d0)
        end if

! ----- init basis
        call hhoBasisCell%initialize(hhoCell)

! ----- compute G_prev = gradrec * depl_prev
        call gradrec%dot(depl_prev, G_prev_coeff)

! ----- compute G_curr = gradrec * depl_curr
        call gradrec%dot(depl_curr, G_curr_coeff)

! ----- Loop on quadrature points
        do kpg = 1, hhoQuadCellRigi%nbQuadPoints
            coorpg(1:3) = hhoQuadCellRigi%points(1:3, kpg)
            weight = hhoQuadCellRigi%weights(kpg)

! --------- Eval basis function at the quadrature point
            call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)

! --------- Eval gradient at T- and T+
            G_prev = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, G_prev_coeff)
            G_curr = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, G_curr_coeff)

! --------- Eval gradient of the deformation at T- and T+
            call hhoCalculF(G_prev, F_prev)
            call hhoCalculF(G_curr, F_curr)

! --------- Check the jacobian jac >= r8prem
! --------- be carrefull with c_plan, I don't know the result
            call lcdetf(hhoCell%ndim, F_prev, jac_prev)
            cod(kpg) = merge(1, 0, jac_prev .le. r8prem())
            if (cod(kpg) .ne. 0) goto 999
!
            call lcdetf(hhoCell%ndim, F_curr, jac_curr)
            cod(kpg) = merge(1, 0, jac_curr .le. r8prem())
            if (cod(kpg) .ne. 0) goto 999

! --------- Set main parameters for behaviour (on point)
            call behaviourSetParaPoin(kpg, ksp, hhoCS%BEHInteg)

! --------- Integrate
            if (l_gdeflog) then
                call gdeflog(hhoCS, hhoCell%ndim, kpg, &
                             time_prev, time_curr, &
                             F_prev, F_curr, Pk1_curr, module_tang, cod(kpg))
            else if (l_green_lagr) then
                call greenlagr(hhoCS, hhoCell%ndim, kpg, &
                               time_prev, time_curr, F_prev, F_curr, &
                               Pk1_curr, module_tang, cod(kpg))
            else
                ASSERT(ASTER_FALSE)
            end if
!
! -------- Test the code of the LDC
!
            if (cod(kpg) .eq. 1) goto 999
!
! ------- Compute rhs
!
            if (l_rhs) then
                call hhoComputeRhsLarge(hhoCell, Pk1_curr, weight, BSCEval, gbs, bT)
            end if
!
! ------- Compute lhs
!
            if (l_lhs) then
                call hhoComputeLhsLarge(hhoCell, module_tang, weight, BSCEval, gbs, AT)
            end if
!
!     print*,"vi_prev", vi_prev(1:lgpg, kpg)
!     print*,"vi_curr", vi_curr(1:lgpg, kpg)
! print*,"sig_prev", sig_prev(1:nbsig, kpg)
! print*,"sig_curr", sig_curr(1:nbsig, kpg)
! print*,"Fp"
! call hhoPrintMat(F_curr)
! print*,"dPK1dF"
! call hhoPrintTensor4(module_tang)
! print*,"dPK1dF"
! call hhoPrintTensor4Mangle(module_tang)
! print*,"PK1p"
! call hhoPrintMat(Pk1_curr)
! print*,"module tangent"
! call hhoPrintTensor4(module_tang)
        end do
!
! ----- compute rhs += Gradrec**T * bT
!
        if (l_rhs) then
            call hho_dgemv_T(1.d0, gradrec, bT, 1.d0, rhs)
        end if
!
! ----- compute lhs += gradrec**T * AT * gradrec
! ----- step1: TMP = AT * gradrec
!
        if (l_lhs) then
            call hhoAssembleLhsLarge(gradrec, AT, lhs)
        end if
! print*, "KT", hhoNorm2Mat(lhs(1:total_dofs,1:total_dofs))
! print*, "fT", norm2(rhs)
!
999     continue
!
! - SYNTHESE DES CODES RETOURS
!
        call codere(cod, hhoQuadCellRigi%nbQuadPoints, hhoCS%codret)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeRhsLarge(hhoCell, stress, weight, BSCEval, gbs, bT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: stress(3, 3)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gbs
        real(kind=8), intent(inout) :: bT(MSIZE_CELL_MAT)
!
! ------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product bT += (stress, gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In stress       : stress tensor
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs          : number of rows of bT
!   Out bT          : contribution of bt
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_stress(3, 3)
        integer(kind=8) :: i, j, gbs_cmp, deca, axis_dim
! --------------------------------------------------------------------------------------------------
!
        axis_dim = 0
        if (hhoCell%l_axis) then
            axis_dim = 1
        end if
        gbs_cmp = gbs/(hhoCell%ndim*hhoCell%ndim+axis_dim)
        qp_stress = weight*stress
! -------- Compute scalar_product of (stress, gphi)_T
! -------- (RAPPEL: the composents of the gradient are saved by rows)
        deca = 0
        do i = 1, hhoCell%ndim
            do j = 1, hhoCell%ndim
                call daxpy_1(gbs_cmp, qp_stress(i, j), BSCEval, bT(deca+1))
                deca = deca+gbs_cmp
            end do
        end do
!
        if (hhoCell%l_axis) then
            call daxpy_1(gbs_cmp, qp_stress(3, 3), BSCEval, bT(deca+1))
            deca = deca+gbs_cmp
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsLarge(hhoCell, module_tang, weight, BSCEval, gbs, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: module_tang(3, 3, 3, 3)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer(kind=8), intent(in) :: gbs
        type(HHO_matrix), intent(inout) :: AT
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (gphi, module_tang:gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In module_tang  : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs          : number of rows of bT
!   Out AT          : contribution of bt
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_Agphi(MSIZE_CELL_MAT, 3, 3)
        integer(kind=8) :: gbs_cmp, d1, d2, ig, row, axis_dim
! --------------------------------------------------------------------------------------------------
!
        axis_dim = 0
        if (hhoCell%l_axis) then
            axis_dim = 1
        end if
        gbs_cmp = gbs/(hhoCell%ndim*hhoCell%ndim+axis_dim)
! --------- Eval (A : gphi)_T
        call hhoComputeAgphi(hhoCell, module_tang, BSCEval, gbs_cmp, weight, qp_Agphi)
!
! -------- Compute scalar_product of (gphi, A:gphi)_T
!
! ----- Gradient is saved by component with row-major - [XX, XY, YX, YY]
        row = 0
        do d1 = 1, hhoCell%ndim
            do d2 = 1, hhoCell%ndim
                do ig = 1, gbs_cmp
                    row = row+1
                    ! only component G_phi(d1,d2) .ne. 0 and G_phi(d1,d2) = BSCEval(ig)
                    call daxpy_1(gbs, BSCEval(ig), qp_Agphi(:, d1, d2), AT%m(row, :))
                end do
            end do
        end do

        if (hhoCell%l_axis) then
            do ig = 1, gbs_cmp
                row = row+1
                ! only component G_phi(3,3) .ne. 0 and G_phi(3,3) = BSCEval(ig)
                ! (gphi, A:Pi^k_T)_T
                call daxpy_1(gbs, BSCEval(ig), qp_Agphi(:, 3, 3), AT%m(row, :))
            end do
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoAssembleLhsLarge(gradrec, AT, lhs)
!
        implicit none
!
        type(HHO_matrix), intent(in) :: gradrec
        type(HHO_matrix), intent(inout) :: AT
        type(HHO_matrix), intent(inout) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics - assemble LHS
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_matrix) :: TMP
        integer(kind=8) :: gbs, total_dofs
! --------------------------------------------------------------------------------------------------
!
!
        gbs = gradrec%nrows
        total_dofs = gradrec%ncols
!
! ----- compute lhs += gradrec**T * AT * gradrec
!
        call TMP%initialize(gbs, total_dofs, 0.d0)
!
! ----- step1: TMP = AT * gradrec
        call hho_dgemm_NN(1.d0, AT, gradrec, 0.d0, TMP)
!
! ----- step2: lhs += gradrec**T * TMP
        call hho_dgemm_TN(1.d0, gradrec, TMP, 1.d0, lhs)
!
        call TMP%free()
        call AT%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function nbsigm_cmp(ndim)
!
        implicit none
!
        integer(kind=8), intent(in) :: ndim
        integer(kind=8) :: nbsigm_cmp
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Return the number of componant of the Cauchy stress tensor
!   In ndim         : the current HHO Cell
! --------------------------------------------------------------------------------------------------
!
        if (ndim == 2) then
            nbsigm_cmp = 4
        else if (ndim == 3) then
            nbsigm_cmp = 6
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeAgphi(hhoCell, module_tang, BSCEval, gbs_cmp, weight, Agphi)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        integer(kind=8), intent(in) :: gbs_cmp
        real(kind=8), intent(in) :: module_tang(3, 3, 3, 3)
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(out) :: Agphi(MSIZE_CELL_MAT, 3, 3)
!
! -----------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product (module_tang, gphi)_T
!
!   In hhoCell      : the current HHO Cell
!   In module_tang  : elasto_plastic moduli
!
!   In BSCEval      : Basis of one composant gphi
!   In gbs          : number of cols of Aphi
!   In weight       : quadrature weight
!   Out Agphi       : matriw of scalar product
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_module_tang(3, 3, 3, 3)
        integer(kind=8) :: col, jg, d1, d2
! --------------------------------------------------------------------------------------------------
!
        Agphi = 0.d0
        qp_module_tang = weight*module_tang
!
! ----- Gradient is saved by component with row-major - [XX, XY, YX, YY]
        col = 0
        do d1 = 1, hhoCell%ndim
            do d2 = 1, hhoCell%ndim
                do jg = 1, gbs_cmp
                    col = col+1
                    ! only component G_phi(d1,d2) .ne. 0 and G_phi(d1,d2) = BSCEval(jg)
                    Agphi(col, 1:3, 1:3) = qp_module_tang(1:3, 1:3, d1, d2)*BSCEval(jg)
                end do
            end do
        end do
!
        if (hhoCell%l_axis) then
            do jg = 1, gbs_cmp
                col = col+1
                ! only component G_phi(d1,d2) .ne. 0 and G_phi(d1,d2) = BSCEval(jg)
                Agphi(col, 1:3, 1:3) = qp_module_tang(1:3, 1:3, 3, 3)*BSCEval(jg)
            end do
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCalculF(G, F)
!
        implicit none
!
        real(kind=8), intent(in) :: G(3, 3)
        real(kind=8), intent(out) :: F(3, 3)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the matrix F = G +I
!   In ndim         : dimension
!   In G            : gradient
!   In F            : gradient of the deformation
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: idim
! --------------------------------------------------------------------------------------------------
!
        F = G
!
        do idim = 1, 3
            F(idim, idim) = F(idim, idim)+1.d0
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoMatrixToVoigth(ndim, mat) result(vect)
!
        implicit none
!
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(in) :: mat(3, 3)
        real(kind=8) :: vect(6)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!  Voigt Notation (XX, YY, ZZ, XY*rac2, XZ*rac2, YZ*rac2)
! --------------------------------------------------------------------------------------------------
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
! --------------------------------------------------------------------------------------------------
!
        if (ndim == 2) then
            vect = [mat(1, 1), mat(2, 2), mat(3, 3), mat(1, 2)*rac2, 0.d0, 0.d0]
        else if (ndim == 3) then
            vect = [mat(1, 1), mat(2, 2), mat(3, 3), mat(1, 2)*rac2, mat(1, 3)*rac2, &
                    mat(2, 3)*rac2]
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCalculGreenLagrange(ndim, F, GLvec)
!
        implicit none
!
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(in) :: F(3, 3)
        real(kind=8), intent(out) :: GLvec(6)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the matrix E = 0.5*(F^T.F-I)
!   In ndim         : dimension
!   In F            : deformation gradient
!   In GLvec        : Green-Lagrange using Voigt Notation (XX, YY, ZZ, XY*rac2, XZ*rac2, YZ*rac2)
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i, j, k
        real(kind=8) :: GL_(3, 3)
! --------------------------------------------------------------------------------------------------
!
! ---- Compute C/2 = F^T*F/2
!
        GL_ = 0.d0
        do j = 1, 3
            do i = 1, j
                do k = 1, 3
                    GL_(i, j) = GL_(i, j)+F(k, i)*F(k, j)
                end do
                GL_(i, j) = 0.5d0*GL_(i, j)
                GL_(j, i) = GL_(i, j)
            end do
        end do
!
        do k = 1, 3
            GL_(k, k) = GL_(k, k)-0.5d0
        end do
!
        GLvec = hhoMatrixToVoigth(ndim, GL_)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCalculFinger(ndim, F, FingerVec)
!
        implicit none
!
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(in) :: F(3, 3)
        real(kind=8), intent(out) :: FingerVec(6)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the matrix Finger = 0.5*(F.F^T-I)
!   In ndim         : dimension
!   In F            : deformation gradient
!   In Fingervec        : Finger using Voigt Notation (XX, YY, ZZ, XY*rac2, XZ*rac2, YZ*rac2)
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i, j, k
        real(kind=8) :: FG_(3, 3)
! --------------------------------------------------------------------------------------------------
!
! ---- Compute B/2 = F^T*F/2
!
        FG_ = 0.d0
        do j = 1, 3
            do i = 1, j
                do k = 1, 3
                    FG_(i, j) = FG_(i, j)+F(i, k)*F(j, k)
                end do
                FG_(i, j) = 0.5d0*FG_(i, j)
                FG_(j, i) = FG_(i, j)
            end do
        end do
!
        do k = 1, 3
            FG_(k, k) = FG_(k, k)-0.5d0
        end do
!
        FingerVec = hhoMatrixToVoigth(ndim, FG_)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCalculEpsDef(ndim, F, G_incr, EDVec)
!
        implicit none
!
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(in) :: F(3, 3), G_incr(3, 3)
        real(kind=8), intent(out) :: EDVec(6)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the matrix EpsDef = 0.5*(grad_x + grad_x^T) with grad_x = G_incr * F^-1
!   In ndim         : dimension
!   In F            : deformation gradient (for the displacement u_x)
!   In G_inr        : gradient in the reference configuration for incremental dispacement
!                     u_incr = u - u_x
!   In EDvec        : EP using Voigt Notation (XX, YY, ZZ, XY*rac2, XZ*rac2, YZ*rac2)
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: ED_(3, 3), Finv(3, 3), detF, grad_x(3, 3)
! --------------------------------------------------------------------------------------------------
!
        call matinv("S", 3, F, Finv, detF)
        grad_x = matmul(G_incr, Finv)
!
        ED_ = 0.5d0*(grad_x+transpose(grad_x))
!
        EDVec = hhoMatrixToVoigth(ndim, ED_)
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine select_behavior(behavior, l_gdeflog, l_green_lagr)
!
        implicit none
!
        character(len=16), intent(in) :: behavior(*)
        aster_logical, intent(out) :: l_gdeflog
        aster_logical, intent(out) :: l_green_lagr
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Select the appropriate behavior
!   In behavior     : type of behavior
!   Out l_gdeflog   : use GDEF_LOG ?
!   Out l_green_lagr: use GREEN_LAGRANGE ?
! --------------------------------------------------------------------------------------------------
!
        l_gdeflog = ASTER_FALSE
        l_green_lagr = ASTER_FALSE
!
        select case (behavior(DEFO))
        case ('GDEF_LOG')
            l_gdeflog = ASTER_TRUE
        case ('GREEN_LAGRANGE')
            l_green_lagr = ASTER_TRUE
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine gdeflog(hhoCS, ndim, kpg, &
                       time_prev, time_curr, &
                       F_prev, F_curr, &
                       PK1_curr, module_tang, cod)
!
        implicit none
!
        type(HHO_Compor_State), intent(inout) :: hhoCS
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: kpg
        real(kind=8), intent(in) :: time_prev
        real(kind=8), intent(in) :: time_curr
        real(kind=8), intent(in) :: F_prev(3, 3)
        real(kind=8), intent(in) :: F_curr(3, 3)
        real(kind=8), intent(out) :: PK1_curr(3, 3)
        real(kind=8), intent(out) :: module_tang(3, 3, 3, 3)
        integer(kind=8), intent(out) :: cod
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the behavior laws for GDEF_LOF
!   IO BEHInteg     : integration informations
!   In ndim         : dimension of the problem
!   In kpg          : i-th quadrature point
!   In time_prev    : previous time T-
!   In time_curr    : current time T+
!   In F_prev       : previous deformation gradient at T-
!   In F_curr       : curr deformation gradient at T+
!   Out Pk1_curr    : piola-kirschooff 1 at T+
!   Out module_tang : tangent modulus dPK1/dF(Fp)
!   Out cod         : info on integration of the LDC
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: gn(3, 3), lamb(3), logl(3), epslPrev(6), epslIncr(6)
        real(kind=8) :: tlogPrev(6), tlogCurr(6)
        real(kind=8) :: dtde(6, 6), PK2_prev(6), PK2_curr(6), sig(6)
        real(kind=8) :: dpk2dc(6, 6), me(3, 3, 3, 3)
        aster_logical :: lCorr, lMatr, lSigm, lVari
!
        lCorr = L_CORR(hhoCS%option)
        lMatr = L_MATR(hhoCS%option)
        lSigm = L_SIGM(hhoCS%option)
        lVari = L_VARI(hhoCS%option)
!
! ----- Compute pre-processing Elog
!
        call prelog(ndim, hhoCS%lgpg, hhoCS%vari_prev((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), &
                    gn, lamb, &
                    logl, F_prev, F_curr, epslPrev, epslIncr, &
                    tlogPrev, lCorr, cod)
        if (cod .ne. 0) then
            goto 999
        end if

! ----- Compute Stress and module_tangent
        dtde = 0.d0
        tlogCurr = 0.d0
        call nmcomp(hhoCS%BEHInteg, &
                    ndim, hhoCS%option, hhoCS%typmod, &
                    time_prev, time_curr, &
                    hhoCS%compor, hhoCS%carcri, hhoCS%multComp, &
                    6, epslPrev, epslIncr, &
                    6, tlogPrev, &
                    hhoCS%vari_prev((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), &
                    tlogCurr, &
                    hhoCS%vari_curr((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), &
                    36, dtde, cod)

! ----- Test the code of the LDC
        if (cod .eq. 1) goto 999

! ----- Compute post-processing Elog
        call poslog(hhoCS%BEHInteg, &
                    lCorr, lMatr, lSigm, lVari, &
                    tlogPrev, tlogCurr, F_prev, &
                    hhoCS%lgpg, &
                    hhoCS%vari_curr((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), ndim, &
                    F_curr, dtde, hhoCS%sig_prev((kpg-1)*hhoCS%nbsigm+1:kpg*hhoCS%nbsigm), &
                    hhoCS%c_plan, time_curr, gn, &
                    lamb, logl, sig, dpk2dc, PK2_prev, &
                    PK2_curr, cod)

! ----- Test the code of the LDC
        if (cod .ne. 0) goto 999
!
        if (.not. lCorr) then
            PK2_curr = PK2_prev
        end if
!
        if (lSigm) then
            hhoCS%sig_curr((kpg-1)*hhoCS%nbsigm+1:kpg*hhoCS%nbsigm) = sig(1:hhoCS%nbsigm)
        end if

! ----- Compute PK1
        call pk2topk1(ndim, PK2_curr, F_curr, Pk1_curr)
!
        module_tang = 0.d0
        if (lMatr) then

! --------- Unpack lagrangian tangent modulus
            call desymt46(dpk2dc, me)

! --------- Compute nominal tangent modulus
            call lagmodtonommod(me, PK2_curr, F_curr, module_tang)
        end if
!
999     continue
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine greenlagr(hhoCS, ndim, kpg, &
                         time_prev, time_curr, F_prev, F_curr, &
                         PK1_curr, module_tang, cod)
!
        implicit none
!
        type(HHO_Compor_State), intent(inout) :: hhoCS
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: kpg
        real(kind=8), intent(in) :: time_prev
        real(kind=8), intent(in) :: time_curr
        real(kind=8), intent(in) :: F_prev(3, 3)
        real(kind=8), intent(in) :: F_curr(3, 3)
        real(kind=8), intent(out) :: PK1_curr(3, 3)
        real(kind=8), intent(out) :: module_tang(3, 3, 3, 3)
        integer(kind=8), intent(out) :: cod
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the behavior laws for GREEN_LAGRANGE
!   In ndim         : dimension of the problem
!   In kpg          : i-th quadrature point
!   In time_prev    : previous time T-
!   In time_curr    : current time T+
!   In F_prev       : previous deformation gradient at T-
!   In F_curr       : curr deformation gradient at T+
!   Out Pk1_curr    : piola-kirschooff 1 at T+
!   Out module_tang : tangent modulus dPK1/dF(Fp)
!   Out cod         : info on integration of the LDC
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: GL_prev(6), GL_curr(6), GL_incr(6), dpk2dc(6, 6), me(3, 3, 3, 3)
        real(kind=8) :: detF_prev, F_incr(3, 3)
        real(kind=8) :: PK2_prev(6), PK2_curr(6), detF_curr, sig(6)
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
        aster_logical :: lMFront
!
        ASSERT(.not. hhoCS%axis)
        lMFront = nint(hhoCS%carcri(EXTE_TYPE)) == 1 .or. nint(hhoCS%carcri(EXTE_TYPE)) == 2
!
! ----- Compute PK2 at T-
!
        PK2_prev = 0.d0
        sig = 0.d0
        sig(1:hhoCS%nbsigm) = hhoCS%sig_prev((kpg-1)*hhoCS%nbsigm+1:kpg*hhoCS%nbsigm)
        call lcdetf(ndim, F_prev, detF_prev)
        call pk2sig(ndim, F_prev, detF_prev, PK2_prev, sig, -1)
        PK2_prev(4:6) = PK2_prev(4:6)*rac2
!
! ----- Compute behaviour
!
        PK2_curr = 0.d0
        dpk2dc = 0.d0
        module_tang = 0.d0
!
        if (lMFront) then
!
! --------- Compute pre-processing F
!
            F_incr = F_curr-F_prev
!
! --------- Compute behaviour
!
            call nmcomp(hhoCS%BEHInteg, &
                        ndim, hhoCS%option, hhoCS%typmod, &
                        time_prev, time_curr, &
                        hhoCS%compor, hhoCS%carcri, hhoCS%multComp, &
                        9, F_prev, F_incr, &
                        6, PK2_prev, &
                        hhoCS%vari_prev((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), &
                        PK2_curr, &
                        hhoCS%vari_curr((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), &
                        36, dpk2dc, cod)
        else
!
! --------- Compute pre-processing E (Green-Lagrange)
!
            call hhoCalculGreenLagrange(ndim, F_prev, GLvec=GL_prev)
            call hhoCalculGreenLagrange(ndim, F_curr, GLvec=GL_curr)
            GL_incr = GL_curr-GL_prev
!
! --------- Compute behaviour
!
            call nmcomp(hhoCS%BEHInteg, &
                        ndim, hhoCS%option, hhoCS%typmod, &
                        time_prev, time_curr, &
                        hhoCS%compor, hhoCS%carcri, hhoCS%multComp, &
                        6, GL_prev, GL_incr, &
                        6, PK2_prev, &
                        hhoCS%vari_prev((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), &
                        PK2_curr, &
                        hhoCS%vari_curr((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), &
                        36, dpk2dc, cod)
        end if
!
! ----- Test the code of the LDC
!
        if (cod .eq. 1) goto 999
!
        if (.not. L_CORR(hhoCS%option)) then
            PK2_curr = PK2_prev
        end if
!
! ----- Compute Cauchy stress and save them
!
        if (L_SIGM(hhoCS%option)) then
            call lcdetf(ndim, F_curr, detF_curr)
            call pk2sig(ndim, F_curr, detF_curr, PK2_curr, &
                        hhoCS%sig_curr((kpg-1)*hhoCS%nbsigm+1:kpg*hhoCS%nbsigm), 1)
        end if
!
! ----- Compute PK1
!
        call pk2topk1(ndim, PK2_curr, F_curr, PK1_curr)
!
        if (L_MATR(hhoCS%option)) then
!
! ----- Unpack lagrangian tangent modulus
            call desymt46(dpk2dc, me)
!
! ----- Compute nominal tangent modulus
!
            call lagmodtonommod(me, PK2_curr, F_curr, module_tang)
        end if
!
999     continue
!
    end subroutine
!
end module
