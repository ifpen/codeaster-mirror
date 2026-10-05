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
subroutine nonlinDSMaterialTimeStep(model, ds_material, caraElem, &
                                    ds_constitutive, hvalIncr, &
                                    nume_dof, sddisc, numeInst)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/diinst.h"
#include "asterfort/detrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/nmvcfo.h"
#include "asterfort/assvec.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmvcle.h"
!
    character(len=24), intent(in) :: model, caraElem
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19), intent(in) :: hvalIncr(*)
    character(len=24), intent(in) :: nume_dof
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: numeInst
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Material parameters management
!
! Update material parameters for new time step
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  caraElem         : name of elementary characteristics (field)
! In  hvalIncr         : hat-variable for incremental values
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  sddisc           : datastructure for time discretization
! In  numeInst         : index of current time step
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: varcCurr
    real(kind=8) :: timeCurr
    character(len=24), parameter :: vectElem = '&&VARCINIT_ELEM'
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_32')
    end if

! - Current time
    timeCurr = diinst(sddisc, numeInst)

! - Get name of fields for external state variables
    call nmchex(hvalIncr, 'VALINC', 'COMPLU', varcCurr)

! - Create external state variables for current time
    call nmvcle(model, ds_material%mater, caraElem, timeCurr, varcCurr)

! - Compute CHAR_MECA_*_R for convergence criteria
    call nmvcfo('+', model, ds_material%mater, ds_material%mateco, &
                caraElem, ds_constitutive%compor, &
                ds_material%varc_refe, hvalIncr, vectElem)
    call assvec('V', ds_material%fvarc_curr, 1, vectElem, [1.d0], nume_dof)
    call detrsd('RESUELEM', vectElem)
!
end subroutine
