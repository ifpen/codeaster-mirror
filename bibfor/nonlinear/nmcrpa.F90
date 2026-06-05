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
subroutine nmcrpa(factorKeywordZ, iFactorKeyword, stepSlctListJv, jvBase, &
                  nbStepSlct, stepSlctMini)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmcrpm.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: factorKeywordZ
    integer(kind=8), intent(in) :: iFactorKeyword
    character(len=24), intent(in) :: stepSlctListJv
    character(len=1), intent(in) :: jvBase
    integer(kind=8), intent(out) :: nbStepSlct
    real(kind=8), intent(out) :: stepSlctMini
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Time selector management
!
! Get list of time step
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword to read
! In  iFactorKeyword   : index of factor keyword
! In  stepSlctListJv   : name of object for list of time step
! In  jvBase           : JEVEUX base to create objects
! Out nbStepSlct       : number of step in list of time step
! Out stepSlctMini     : minimum length of time step in list of time step
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: n2, n3, iret
    character(len=19) :: list
    character(len=16) :: factorKeyword
    real(kind=8), pointer :: vale(:) => null()
    real(kind=8), pointer :: stepSlctList(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    nbStepSlct = 0
    factorKeyword = factorKeywordZ
    stepSlctMini = 0.d0

! - CREATION ET INITIALISATION SD
    call getvid(factorKeyword, 'LIST_INST', iocc=iFactorKeyword, scal=list, nbret=n2)
    call getvr8(factorKeyword, 'INST', iocc=iFactorKeyword, nbval=0, nbret=n3)
    n3 = -n3

! - Get number of time step
    if ((n2 .ge. 1) .and. (n3 .ge. 1)) then
        ASSERT(ASTER_FALSE)
    end if
    if (n3 .ge. 1) then
        nbStepSlct = n3
    else if (n2 .ge. 1) then
        call jelira(list//'.VALE', 'LONMAX', ival=nbStepSlct)
    else
        nbStepSlct = 0
        goto 99
    end if

! - Create object for list
    call wkvect(stepSlctListJv, jvBase//' V R', nbStepSlct, vr=stepSlctList)

! - Set list
    if (n3 .ge. 1) then
        call getvr8(factorKeyword, 'INST', iocc=iFactorKeyword, nbval=nbStepSlct, &
                    vect=stepSlctList, nbret=iret)
    else
        call jeveuo(list//'.VALE', 'L', vr=vale)
        stepSlctList(1:nbStepSlct) = vale(1:nbStepSlct)
    end if

! - Get minimum length of time step in list
    call nmcrpm(stepSlctList, nbStepSlct, stepSlctMini)
!
99  continue
!
    call jedema()
!
end subroutine
