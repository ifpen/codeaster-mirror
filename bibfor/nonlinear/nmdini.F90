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
subroutine nmdini(factorKeyword, listInstJv, tole, nbInst, numeInstInit)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utacli.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: factorKeyword
    character(len=19), intent(in) :: listInstJv
    real(kind=8), intent(in) :: tole
    integer(kind=8), intent(in) :: nbInst
    integer(kind=8), intent(out) :: numeInstInit
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Time discretization datastructure
!
! Index of initial time in list of times
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword
! In  listInstJv       : name of JEVEUX object for list of times from INCREMENT/LIST_INST
! In  tole             : tolerance to search time
! In  nbInst           : number of time steps in list
! Out numeInstInit     : index of initial time in list of times
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: instInit
    integer(kind=8) :: n1, n2, iInst
    real(kind=8), pointer :: listInst(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    numeInstInit = 0

! - Acces to list of time steps
    call jeveuo(listInstJv, 'L', vr=listInst)

! - Get keywords
    call getvis(factorKeyword, 'NUME_INST_INIT', iocc=1, scal=numeInstInit, nbret=n1)
    call getvr8(factorKeyword, 'INST_INIT', iocc=1, scal=instInit, nbret=n2)

! - Find index of initial time step
    if (n1 .eq. 0) then
        if (n2 .eq. 0) then
            numeInstInit = 0
        else
            call utacli(instInit, listInst, nbInst, tole, numeInstInit)
            if (numeInstInit .lt. 0) then
                do iInst = 1, nbInst
                    if (listInst(iInst) .ge. instInit) then
                        numeInstInit = iInst
                        exit
                    end if
                end do
                call utmess('I', 'DISCRETISATION_90', sr=listInst(numeInstInit))
                if (numeInstInit .lt. 2) then
                    call utmess('F', 'DISCRETISATION_91')
                else
                    numeInstInit = numeInstInit-2
                end if
            end if
        end if
    end if

! - Checks
    ASSERT(numeInstInit .ge. 0)
    ASSERT(numeInstInit .le. nbInst)
!
end subroutine
