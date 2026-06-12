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
subroutine storeSaveLast(sdarchZ, numeStore, timeCurr)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeexin.h"
!
    character(len=*), intent(in) :: sdarchZ
    integer(kind=8), intent(in) :: numeStore
    real(kind=8), intent(in) :: timeCurr
!
! --------------------------------------------------------------------------------------------------
!
! Storing management
!
! Save current storing
!
! --------------------------------------------------------------------------------------------------
!
! In  sdarch           : datastructure for storing
! In  numeStore        : current storing index
! In  timeCurr         : current time step
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=19) :: sdarch
    character(len=24) :: sdarchLastJv
    real(kind=8), pointer :: sdarchLast(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sdarch = sdarchZ
    sdarchLastJv = sdarch(1:19)//'.LAST'
    call jeexin(sdarchLastJv, iret)
    if (iret .ne. 0) then
        call jeveuo(sdarchLastJv, 'E', vr=sdarchLast)
        if (numeStore .ge. 0) then
            sdarchLast(1) = numeStore
            sdarchLast(2) = timeCurr
        else
            sdarchLast(1:2) = 1.d0
        end if
    end if
!
end subroutine
