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
subroutine te0425(option, nomte)
!
    use HHO_basis_module
    use HHO_compor_module
    use HHO_eval_module
    use HHO_GV_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_LargeStrainMeca_module
    use HHO_matrix_module
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
    use HHO_utils_module
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
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
! Mechanics - STAT_NON_LINE - Pilotage - GRAD_HHO
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI", typmod2 = "GRADVARI"
    type(HHO_Data) :: hhoDataMk, hhoDataGv
    type(HHO_Cell) :: hhoCell
    type(HHO_Meca_State) :: hhoMecaState
    type(HHO_Quadrature) :: hhoQuadCellRigi
    type(HHO_Compor_State) :: hhoCS
    type(HHO_basis_cell) :: hhoBasisCell
    type(HHO_GV_State):: hhoGVState
    integer(kind=8), parameter :: nmax = 11
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
    real(kind=8), dimension(MSIZE_TDOFS_MIX) :: tmp_0, tmp_pilo
    real(kind=8), dimension(MSIZE_TDOFS_VEC) :: depl_0, depl_pilo
    real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: vari_0, vari_pilo
    real(kind=8), dimension(MSIZE_CELL_SCAL) :: lagv_0, lagv_pilo
    real(kind=8), dimension(MSIZE_CELL_VEC) :: GV_prev_coeff, GV_incr_coeff
    real(kind=8), dimension(MSIZE_CELL_VEC) :: GV_0_coeff, GV_pilo_coeff
    real(kind=8), dimension(MSIZE_CELL_MAT) :: E_prev_coeff, E_incr_coeff
    real(kind=8), dimension(MSIZE_CELL_MAT) :: E_0_coeff, E_pilo_coeff
    real(kind=8) :: tau, etamin, etamax, coorpg(3)
    real(kind=8) :: E_prev(nmax), E_cste(nmax), E_pilo(nmax), E_0(nmax), E_incr(nmax)
    real(kind=8) :: copilo(5), sigma(nmax)
    real(kind=8) :: G_prev(3, 3), G_incr(3, 3), G_pilo(3, 3), G_0(3, 3)
    real(kind=8) :: F_prev(3, 3)
    real(kind=8) :: GV_prev(3), GV_incr(3), GV_0(3), GV_pilo(3), GV_cste(3)
    real(kind=8) :: var_prev, var_incr, var_0, var_pilo, var_cste
    real(kind=8) :: lag_prev, lag_incr, lag_0, lag_pilo, lag_cste
    integer(kind=8) :: mk_cbs, mk_fbs, mk_total_dofs, mk_gbs, mk_gbs_sym, mkb_gbs_axis
    integer(kind=8) :: gv_cbs, gv_fbs, gv_total_dofs, gv_gbs, total_dofs
    integer(kind=8) :: ipg, npg, k, neps, nmk, gv_faces_dofs, gv_cell_offset
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
    call hhoInfoInitCell(hhoCell, hhoDataMk)
    call hhoDataGVInit(hhoDataGv)

! - Number of dofs
    call hhoMecaNLDofs(hhoCell, hhoDataMk, mk_cbs, mk_fbs, mk_total_dofs, &
                       mk_gbs, mk_gbs_sym, mkb_gbs_axis)
    call hhoTherNLDofs(hhoCell, hhoDataGv, gv_cbs, gv_fbs, gv_total_dofs, gv_gbs)
    total_dofs = mk_total_dofs+gv_total_dofs
    gv_faces_dofs = gv_total_dofs-gv_cbs
    gv_cell_offset = gv_faces_dofs+1

! - Initialize quadrature for the rigidity
    call hhoQuadCellRigi%initCell(hhoCell, npg)

! - Type of finite element
    call hhoCS%initialize(fami, option, hhoCell%ndim, hhoCell%barycenter, typmod2)
    hhoCS%typmod(2) = 'GRADVARI'
    call hhoMecaState%initialize(hhoCell, hhoDataMk, hhoCS, hhoDataGv)
    call hhoGVState%initialize(hhoCell, hhoDataMk, hhoDataGv, hhoCS)
    call hhoBasisCell%initialize(hhoCell)

! - Compute Operators
    hhoCS%l_largestrain = hhoCS%l_largestrain .and. (.not. l_pred_elas)
    call hhoCalcOpGv(hhoCell, hhoDataMk, hhoDataGv, hhoCS%l_largestrain, &
                     hhoMecaState, hhoGvState)
!
    call readVector('PDEPL0R', total_dofs, tmp_0)
    call readVector('PDEPL1R', total_dofs, tmp_pilo)
!
    call hhoExtrField(hhoCell, hhoDataMk, hhoDataGv, &
                      tmp_0, depl_0, vari_0, lagv_0)
    call hhoExtrField(hhoCell, hhoDataMk, hhoDataGv, &
                      tmp_pilo, depl_pilo, vari_pilo, lagv_pilo)
!
    call jevech('PCDTAU', 'L', ictau)
    tau = zr(ictau)
    if (pilo .eq. 'PRED_ELAS') then
        call jevech('PBORNPI', 'L', iborne)
        etamin = zr(iborne+1)
        etamax = zr(iborne)
    end if

! - Compute E_prev = gradrec_sym * depl_prev
    call hhoMecaState%grad%dot(hhoMecaState%depl_prev, E_prev_coeff)
    call hhoMecaState%grad%dot(hhoMecaState%depl_incr, E_incr_coeff)
    call hhoMecaState%grad%dot(depl_0, E_0_coeff)
    call hhoMecaState%grad%dot(depl_pilo, E_pilo_coeff)
!
    call hhoGVState%grad%dot(hhoGVState%vari_prev, GV_prev_coeff)
    call hhoGVState%grad%dot(hhoGVState%vari_incr, GV_incr_coeff)
    call hhoGVState%grad%dot(vari_0, GV_0_coeff)
    call hhoGVState%grad%dot(vari_pilo, GV_pilo_coeff)

!
    call jevech('PCOPILO', 'E', vr=v_copilo)
    copilo = r8vide()
!
    nmk = 2*hhoCell%ndim
    neps = nmk+hhoCell%ndim+2
    sigma = 0.d0
!
    do ipg = 1, hhoQuadCellRigi%nbQuadPoints
        coorpg(1:3) = hhoQuadCellRigi%points(1:3, ipg)
! --------- Eval basis function at the quadrature point
        call hhoBasisCell%BSEval(coorpg(1:3), 0, &
                                 max(hhoDataMk%grad_degree(), hhoDataMk%cell_degree(), &
                                     hhoDataGv%grad_degree(), hhoDataGv%cell_degree()), &
                                 BSCEval)
!
        if (hhoCS%l_largestrain) then
            G_prev = hhoEvalMatCell(hhoCell%ndim, mk_gbs, BSCEval, E_prev_coeff)
            call hhoCalculF(G_prev, F_prev)
            call hhoCalculFinger(hhoCell%ndim, F_prev, E_prev(1:6))
!
            G_incr = hhoEvalMatCell(hhoCell%ndim, mk_gbs, BSCEval, E_incr_coeff)
            call hhoCalculEpsDef(hhoCell%ndim, F_prev, G_incr, E_incr(1:6))
!
            G_0 = hhoEvalMatCell(hhoCell%ndim, mk_gbs, BSCEval, E_0_coeff)
            call hhoCalculEpsDef(hhoCell%ndim, F_prev, G_0, E_0(1:6))
!
            E_cste(1:6) = E_incr(1:6)+E_0(1:6)
!
            G_pilo = hhoEvalMatCell(hhoCell%ndim, mk_gbs, BSCEval, E_pilo_coeff)
            call hhoCalculEpsDef(hhoCell%ndim, F_prev, G_pilo, E_pilo(1:6))
        else
! --------- Eval deformations
            E_prev(1:6) = hhoEvalSymMatCell(hhoCell%ndim, mk_gbs_sym, BSCEval, E_prev_coeff)
            E_incr(1:6) = hhoEvalSymMatCell(hhoCell%ndim, mk_gbs_sym, BSCEval, E_incr_coeff)
            E_0(1:6) = hhoEvalSymMatCell(hhoCell%ndim, mk_gbs_sym, BSCEval, E_0_coeff)
            E_cste(1:6) = E_incr(1:6)+E_0(1:6)
            E_pilo(1:6) = hhoEvalSymMatCell(hhoCell%ndim, mk_gbs_sym, BSCEval, E_pilo_coeff)
        end if
!
        GV_prev = hhoEvalVecCell(hhoCell%ndim, gv_gbs, BSCEval, GV_prev_coeff)
        GV_incr = hhoEvalVecCell(hhoCell%ndim, gv_gbs, BSCEval, GV_incr_coeff)
        GV_0 = hhoEvalVecCell(hhoCell%ndim, gv_gbs, BSCEval, GV_0_coeff)
        GV_cste = GV_incr+GV_0
        GV_pilo = hhoEvalVecCell(hhoCell%ndim, gv_gbs, BSCEval, GV_pilo_coeff)
!
        var_prev = hhoEvalScalCell(gv_cbs, BSCEval, hhoGVState%vari_prev(gv_cell_offset:))
        var_incr = hhoEvalScalCell(gv_cbs, BSCEval, hhoGVState%vari_incr(gv_cell_offset:))
        var_0 = hhoEvalScalCell(gv_cbs, BSCEval, vari_0(gv_cell_offset:))
        var_cste = var_incr+var_0
        var_pilo = hhoEvalScalCell(gv_cbs, BSCEval, vari_pilo(gv_cell_offset:))
!
        lag_prev = hhoEvalScalCell(gv_cbs, BSCEval, hhoGVState%lagv_prev)
        lag_incr = hhoEvalScalCell(gv_cbs, BSCEval, hhoGVState%lagv_incr)
        lag_0 = hhoEvalScalCell(gv_cbs, BSCEval, lagv_0)
        lag_cste = lag_incr+lag_0
        lag_pilo = hhoEvalScalCell(gv_cbs, BSCEval, lagv_pilo)
!
        E_prev(nmk+1) = var_prev
        E_prev(nmk+2) = lag_prev
        E_prev(nmk+2+1:neps) = GV_prev(1:hhoCell%ndim)
!
        E_cste(nmk+1) = var_cste
        E_cste(nmk+2) = lag_cste
        E_cste(nmk+2+1:neps) = GV_cste(1:hhoCell%ndim)
!
        E_pilo(nmk+1) = var_pilo
        E_prev(nmk+2) = lag_pilo
        E_pilo(nmk+2+1:neps) = GV_pilo(1:hhoCell%ndim)
!
! --- PILOTAGE PAR L'INCREMENT DE DEFORMATION
!
        if (option == "PILO_PRED_ELAS") then
            sigma(1:neps) = hhoCS%sig_prev((ipg-1)*neps+1:ipg*neps)
            do k = 4, nmk
                sigma(k) = sigma(k)*rac2
            end do
            call pi0000(hhoCS%BEHInteg, hhoCS%compor, hhoCS%typmod, hhoCell%ndim, &
                        E_prev(1:nmk), E_cste(1:nmk), E_pilo(1:nmk), &
                        sigma(1:nmk), hhoCS%vari_prev((ipg-1)*hhoCS%lgpg+1:ipg*hhoCS%lgpg), &
                        tau, etamin, etamax, copilo)

        else if (option == "PILO_PRED_DEFO") then
            call pidefo(E_prev(1:nmk), E_cste(1:nmk), E_pilo(1:nmk), tau, copilo)

        else
            ASSERT(ASTER_FALSE)
        end if
!
        v_copilo((ipg-1)*5+1:ipg*5) = copilo(1:5)
    end do
!
end subroutine
