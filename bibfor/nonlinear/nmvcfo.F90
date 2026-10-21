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
subroutine nmvcfo(poum, model, materField, materCode, caraElem, compor, &
                  varcRefe, hvalIncr, vectElem)
!
    use ExternalStateVariableComp_module, only: varcCompElem
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/nmchex.h"
#include "asterfort/reajre.h"
#include "asterfort/vemare.h"
!
    character(len=1), intent(in) :: poum
    character(len=24), intent(in) :: model, materField, materCode, caraElem, compor
    character(len=24), intent(in) :: varcRefe
    character(len=19), intent(in) :: hvalIncr(*)
    character(len=19), intent(in) :: vectElem
!
! --------------------------------------------------------------------------------------------------
!
! Nonlinear mechanics (algorithm)
!
! Command variables - Vector for reference (residual evaluation)
!
! --------------------------------------------------------------------------------------------------
!
! In  poum           : type of computation
!                      '-' - Previous step
!                      '+' - Current step
! In  model          : name of model
! In  materField     : name of material characteristics (field)
! In  materCode      : name of coded material
! In  caraElem       : name of elementary characteristics (field)
! In  varcRefe       : name of reference command variables vector
! In  compor         : name of comportment definition (field)
! In  hvalIncr       : hat-variable for incremental values
! In  vect_elem      : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: jvBase = "V"
    aster_logical, parameter :: lSuperElem = ASTER_FALSE, lCumul = ASTER_FALSE
    integer(kind=8), parameter :: numeHarm = 0
    character(len=19) :: sigmPrev, variPrev, varcPrev, varcCurr
    integer(kind=8) :: iret
!
! --------------------------------------------------------------------------------------------------
!

! - Get fields from hat-variables
    call nmchex(hvalIncr, 'VALINC', 'SIGMOI', sigmPrev)
    call nmchex(hvalIncr, 'VALINC', 'VARMOI', variPrev)
    call nmchex(hvalIncr, 'VALINC', 'COMMOI', varcPrev)
    call nmchex(hvalIncr, 'VALINC', 'COMPLU', varcCurr)

! - Prepare elementary vectors
    call jeexin(vectElem(1:19)//'.RELR', iret)
    if (iret .eq. 0) then
        call vemare('V', vectElem, model)
    end if
    call jedetr(vectElem(1:19)//'.RELR')
    call reajre(vectElem, ' ', 'V')

! - Computation of elementary vectors
    call varcCompElem(lSuperElem, poum, &
                      numeHarm, model, caraElem, materField, materCode, &
                      varcRefe, varcPrev, &
                      jvBase, vectElem, &
                      lCumul_=lCumul, varcCurrZ_=varcCurr, variZ_=variPrev, sigmZ_=sigmPrev, &
                      comporZ_=compor)
!
end subroutine
