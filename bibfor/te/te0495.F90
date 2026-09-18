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
subroutine te0495(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use HHO_basis_module
    use HHO_compor_module
    use HHO_eval_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_LargeStrainMeca_module
    use HHO_matrix_module
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
    use HHO_utils_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/pi0000.h"
#include "asterfort/pidefo.h"
#include "asterfort/readVector.h"
!
    character(len=16), intent(in) :: nomte, option
!
! --------------------------------------------------------------------------------------------------
!
! HHO
! Mechanics - STAT_NON_LINE - Pilotage
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI", typmod2 = "HHO"
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_Meca_State) :: hhoMecaState
    type(HHO_Quadrature) :: hhoQuadCellRigi
    type(HHO_Compor_State) :: hhoCS
    type(HHO_basis_cell) :: hhoBasisCell
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
    real(kind=8), dimension(MSIZE_TDOFS_VEC) :: depl_0r, depl_1r
    real(kind=8), dimension(MSIZE_CELL_MAT) :: E_prev_coeff, E_incr_coeff
    real(kind=8), dimension(MSIZE_CELL_MAT) :: E_0r_coeff, E_1r_coeff
    real(kind=8) :: tau, etamin, etamax
    real(kind=8) :: coorpg(3), E_prev(6), E_pilo(6), E_0(6), E_incr(6)
    real(kind=8) :: copilo(5), sigma(6), E_cst(6)
    real(kind=8) :: G_prev(3, 3), G_incr(3, 3), G_1(3, 3), G_0(3, 3)
    real(kind=8) :: F_prev(3, 3)
    integer(kind=8) :: cbs, fbs, total_dofs, gbs, gbs_sym, gbs_axis
    integer(kind=8) :: kpg, npg, k
    integer(kind=8) :: iborne, ictau, itype
    character(len=16) :: pilo
    aster_logical :: l_pred_elas
    real(kind=8), pointer :: v_copilo(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    if (option /= "PILO_PRED_DEFO" .and. option /= "PILO_PRED_ELAS") then
        ASSERT(ASTER_FALSE)
    end if
!
    call jevech('PTYPEPI', 'L', itype)
    pilo = zk16(itype)
    l_pred_elas = pilo == "PILO_PRED_ELAS"
!
! - Get element parameters
    call elrefe_info(fami=fami, npg=npg)

! - Get HHO data on the modelisation
    call hhoInfoInitCell(hhoCell, hhoData)

! - Number of dofs
    call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs, gbs_sym, gbs_axis)
    ASSERT(total_dofs <= MSIZE_TDOFS_VEC)

! - Initialize quadrature for the rigidity
    call hhoQuadCellRigi%initCell(hhoCell, npg)

! - Type of finite element
    call hhoCS%initialize(fami, option, hhoCell%ndim, hhoCell%barycenter, typmod2)
    call hhoMecaState%initialize(hhoCell, hhoData, hhoCS)
    call hhoBasisCell%initialize(hhoCell)

! - Compute Operators
    hhoCS%l_largestrain = hhoCS%l_largestrain .and. (.not. l_pred_elas)
    if (hhoData%precompute()) then
        call hhoReloadPreCalcMeca(hhoCell, hhoData, hhoCS%l_largestrain, hhoMecaState%grad)
    else
        call hhoCalcOpMeca(hhoCell, hhoData, hhoCS%l_largestrain, hhoMecaState%grad)
    end if
!
    call readVector('PDEPL0R', total_dofs, depl_0r)
    call readVector('PDEPL1R', total_dofs, depl_1r)
!
    call jevech('PCOPILO', 'E', vr=v_copilo)
!
    call jevech('PCDTAU', 'L', ictau)
    tau = zr(ictau)
    if (pilo .eq. 'PRED_ELAS') then
        call jevech('PBORNPI', 'L', iborne)
        etamin = zr(iborne+1)
        etamax = zr(iborne)
    end if

! - Prepare external state variables (geometry)
    call behaviourPrepESVAGeomHHO(hhoCell, hhoQuadCellRigi, hhoCS%BEHInteg)

! - compute E_prev = gradrec_sym * depl_prev
    call hhoMecaState%grad%dot(hhoMecaState%depl_prev, E_prev_coeff)
    call hhoMecaState%grad%dot(hhoMecaState%depl_incr, E_incr_coeff)
    call hhoMecaState%grad%dot(depl_0r, E_0r_coeff)
    call hhoMecaState%grad%dot(depl_1r, E_1r_coeff)
!
    call matfpe(-1)
    copilo = r8vide()
    sigma = 0.d0
!
    do kpg = 1, hhoQuadCellRigi%nbQuadPoints
        coorpg(1:3) = hhoQuadCellRigi%points(1:3, kpg)
! ----- Eval basis function at the quadrature point
        call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)
!
        if (hhoCS%l_largestrain) then
            G_prev = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, E_prev_coeff)
            call hhoCalculF(G_prev, F_prev)
            call hhoCalculFinger(hhoCell%ndim, F_prev, E_prev)
            G_incr = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, E_incr_coeff)
            call hhoCalculEpsDef(hhoCell%ndim, F_prev, G_incr, E_incr)
            G_0 = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, E_0r_coeff)
            call hhoCalculEpsDef(hhoCell%ndim, F_prev, G_0, E_0)
            E_cst = E_incr+E_0
            G_1 = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, E_1r_coeff)
            call hhoCalculEpsDef(hhoCell%ndim, F_prev, G_1, E_pilo)
        else
! --------- Eval deformations
            E_prev = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_prev_coeff)
            E_incr = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_incr_coeff)
            E_0 = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_0r_coeff)
            E_cst = E_incr+E_0
            E_pilo = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_1r_coeff)
        end if
!
        if (pilo .eq. 'DEFORMATION') then
            call pidefo(E_prev, E_cst, E_pilo, tau, copilo)
        else if (pilo .eq. 'PRED_ELAS') then
            sigma(1:hhoCS%nbsigm) = hhoCS%sig_prev((kpg-1)*hhoCS%nbsigm+1:kpg*hhoCS%nbsigm)
            do k = 4, hhoCS%nbsigm
                sigma(k) = sigma(k)*rac2
            end do
            call pi0000(hhoCS%BEHInteg, hhoCS%compor, hhoCS%typmod, hhoCell%ndim, &
                        E_prev(1:hhoCS%nbsigm), E_cst(1:hhoCS%nbsigm), E_pilo(1:hhoCS%nbsigm), &
                        sigma(1:hhoCS%nbsigm), &
                        hhoCS%vari_prev((kpg-1)*hhoCS%lgpg+1:kpg*hhoCS%lgpg), &
                        tau, etamin, etamax, copilo)
        else
            ASSERT(ASTER_FALSE)
        end if
!
        v_copilo((kpg-1)*5+1:kpg*5) = copilo(1:5)
    end do
!
    call matfpe(1)
!
end subroutine
