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
subroutine nonlinDSMaterialInit(model, materField, materCode, caraElem, &
                                compor, hvalIncr, &
                                nume_dof, timeInit, &
                                ds_material)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/detrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/nmvcre.h"
#include "asterfort/nmvcfo.h"
#include "asterfort/assvec.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmvcle.h"
!
    character(len=24), intent(in) :: model, materCode, caraElem, materField
    character(len=24), intent(in) :: compor
    character(len=19), intent(in) :: hvalIncr(*)
    character(len=24), intent(in) :: nume_dof
    real(kind=8), intent(in) :: timeInit
    type(NL_DS_Material), intent(inout) :: ds_material
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Material parameters management
!
! Initializations for material parameters management
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  materField       : name of material characteristics (field)
! In  materCode        : name of coded material
! In  caraElem         : name of elementary characteristics (field)
! In  compor           : name of comportment definition (field)
! In  hvalIncr         : hat-variable for incremental values
! In  nume_dof         : name of numbering object (NUME_DDL)
! IO  ds_material      : datastructure for material parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: varcPrev
    character(len=24), parameter :: vectElem = '&&VARCINIT_ELEM'
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_8')
    end if

! - Set material fields in datastructure
    ds_material%mateco = materCode
    ds_material%mater = materField

! - Create external state variables for reference state
    call nmvcre(model, materField, caraElem, ds_material%varc_refe)

! - Get name of field for external state variables
    call nmchex(hvalIncr, 'VALINC', 'COMMOI', varcPrev)

! - Create external state variables at initial time step
    call nmvcle(model, materField, caraElem, timeInit, varcPrev)

! - Compute CHAR_MECA_*_R at initial time step
    call nmvcfo('-', model, materField, materCode, caraElem, compor, &
                ds_material%varc_refe, hvalIncr, vectElem)
    call assvec('V', ds_material%fvarc_init, 1, vectElem, [1.d0], nume_dof)
    call detrsd('RESUELEM', vectElem)
!
end subroutine
