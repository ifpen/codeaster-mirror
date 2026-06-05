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
subroutine nmcrpx(factorKeywordZ, stepKeywordZ, iFactorKeyword, stepSlct, jvBase)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmcrpa.h"
#include "asterfort/nmcrpp.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: factorKeywordZ, stepKeywordZ
    integer(kind=8), intent(in) :: iFactorKeyword
    character(len=19), intent(in) :: stepSlct
    character(len=1), intent(in) :: jvBase
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Time selector management
!
! Read parameters from user and create datastructure for time selector
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword to read
! In  stepKeyword      : keyword for step
! In  iFactorKeyword   : index of factor keyword
! In  stepSlct         : name of object to time selector
! In  jvBase           : JEVEUX base to create object
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: factorKeyword, stepKeyword
    real(kind=8) :: stepSlctMini, stepSlctTole
    integer(kind=8) :: nbStepSlct, n1, stepSlctFreq
    character(len=24) :: stepSlctListJv, stepSlctInflJv
    real(kind=8), pointer :: stepSlctInfl(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - INITIALISATIONS
    factorKeyword = factorKeywordZ
    stepKeyword = stepKeywordZ
    nbStepSlct = 0
    stepSlctFreq = 0

! - Names of datastructure: list of time step and object for parameters of list management
    stepSlctListJv = stepSlct(1:19)//'.LIST'
    stepSlctInflJv = stepSlct(1:19)//'.INFL'
    call wkvect(stepSlctInflJv, jvBase//' V R', 4, vr=stepSlctInfl)

! - Get parameters from user
    if (factorKeyword .eq. ' ') then
        stepSlctFreq = 1
        stepSlctTole = 0.d0
        nbStepSlct = 0
        stepSlctMini = 0.d0
    else
! ----- Get tolerance to select one time step
        call nmcrpp(factorKeyword, iFactorKeyword, stepSlctTole)

! ----- Get list of time step to select
        call nmcrpa(factorKeyword, iFactorKeyword, stepSlctListJv, jvBase, &
                    nbStepSlct, stepSlctMini)

! ----- Get frequency step
        n1 = 0
        if (nbStepSlct .eq. 0) then
            call getvis(factorKeyword, stepKeyword, iocc=iFactorKeyword, &
                        scal=stepSlctFreq, nbret=n1)
            if (n1 .ne. 0) then
                ASSERT(stepSlctFreq .ge. 0)
            end if
        end if

! ----- AUCUN MOT-CLE : PAS  = 1
        if (n1+nbStepSlct .eq. 0) then
            stepSlctFreq = 1
        end if

    end if

    stepSlctInfl(1) = stepSlctFreq
    stepSlctInfl(2) = stepSlctTole
    stepSlctInfl(3) = nbStepSlct
    stepSlctInfl(4) = stepSlctMini
!
    call jedema()
!
end subroutine
