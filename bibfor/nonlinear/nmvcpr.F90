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
subroutine nmvcpr(modelZ, caraElemZ, hvalIncr, &
                  ds_material, ds_constitutive)
!
    use NonLin_Datastructure_type
    use ExternalStateVariableComp_module, only: varcCompElem
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/nmchex.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: modelZ, caraElemZ
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19), intent(in) :: hvalIncr(*)
!
! --------------------------------------------------------------------------------------------------
!
! Nonlinear mechanics (algorithm)
!
! Command variables - Second member for prediction
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  caraElem         : name of elementary characteristics (field)
! In  hvalIncr         : hat-variable for incremental values
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: jvBase = "V"
    aster_logical, parameter :: lSuperElem = ASTER_FALSE
    integer(kind=8), parameter :: numeHarm = 0
    integer(kind=8) :: ifm, niv
    character(len=8) :: vevcprCurr, vevcprPrev
    character(len=19) :: varcRefe, compor, varcPrev, varcCurr
    character(len=24) :: materField, materCode, sigmPrev, variPrev
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_14')
    end if

! - Initializations
    materField = ds_material%mater
    materCode = ds_material%mateco
    varcRefe = ds_material%varc_refe(1:19)
    compor = ds_constitutive%compor(1:19)
    vevcprPrev = ds_material%vevcprPrev
    vevcprCurr = ds_material%vevcprCurr

! - Get fields from hat-variables
    call nmchex(hvalIncr, 'VALINC', 'SIGMOI', sigmPrev)
    call nmchex(hvalIncr, 'VALINC', 'VARMOI', variPrev)
    call nmchex(hvalIncr, 'VALINC', 'COMMOI', varcPrev)
    call nmchex(hvalIncr, 'VALINC', 'COMPLU', varcCurr)

! - Compute elementary vectors - Previous
    call varcCompElem(lSuperElem, '-', &
                      numeHarm, modelZ, caraElemZ, materField, materCode, &
                      varcRefe, varcPrev, &
                      jvBase, vevcprPrev, &
                      variZ_=variPrev, sigmZ_=sigmPrev, comporZ_=compor)

! - Compute elementary vectors - Current
    call varcCompElem(lSuperElem, '+', &
                      numeHarm, modelZ, caraElemZ, materField, materCode, &
                      varcRefe, varcPrev, &
                      jvBase, vevcprCurr, &
                      varcCurrZ_=varcCurr, variZ_=variPrev, sigmZ_=sigmPrev, comporZ_=compor)
!
end subroutine
