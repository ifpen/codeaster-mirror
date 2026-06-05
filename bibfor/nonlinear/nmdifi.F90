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
subroutine nmdifi(factorKeyword, listInstJv, tole, nbInst, numeInstEnd)
!
    implicit none
!
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
    integer(kind=8), intent(out) :: numeInstEnd
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Time discretization datastructure
!
! Index of final time in list of times
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword
! In  listInstJv       : name of JEVEUX object for list of times from INCREMENT/LIST_INST
! In  tole             : tolerance to search time
! In  nbInst           : number of time steps in list
! Out numeInstEnd      : index of final time in list of times
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: n1, n2
    real(kind=8) :: instEnd
    real(kind=8), pointer :: listInst(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    numeInstEnd = 0

! - Acces to list of time steps
    call jeveuo(listInstJv, 'L', vr=listInst)

! - Get keywords
    call getvis(factorKeyword, 'NUME_INST_FIN', iocc=1, scal=numeInstEnd, nbret=n1)
    call getvr8(factorKeyword, 'INST_FIN', iocc=1, scal=instEnd, nbret=n2)

! - Find index of final time step
    if (n1+n2 .eq. 0) then
        numeInstEnd = nbInst-1
    else if (n1 .eq. 0) then
        call utacli(instEnd, listInst, nbInst, tole, numeInstEnd)
    end if

! - Checks
    if (numeInstEnd .lt. 0 .or. numeInstEnd .gt. (nbInst-1)) then
        call utmess('F', 'DISCRETISATION_94')
    end if
!
end subroutine
