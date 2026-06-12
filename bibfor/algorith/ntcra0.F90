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
subroutine ntcra0(sdarchZ)
!
    implicit none
!
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmcrpx.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: sdarchZ
!
! --------------------------------------------------------------------------------------------------
!
! THER_* - Storing management
!
! Create storing objects for transient computation
!
! --------------------------------------------------------------------------------------------------
!
! In  sdarch           : datastructure for storing
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: iocc = 0
    character(len=16), parameter :: factorKeyword = ' ', stepKeyword = ' '
    character(len=1), parameter :: jvBase = 'V'
    character(len=19) :: sdarch
    character(len=24) :: sdarchAinfJv, sdarchAexcJv, sdarchLastJv
    integer(kind=8), pointer :: sdarchAinf(:) => null()
    character(len=16), pointer :: sdarchAexc(:) => null()
    real(kind=8), pointer :: sdarchLast(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    sdarch = sdarchZ

! - Create datastructures
    sdarchAinfJv = sdarch(1:19)//'.AINF'
    sdarchAexcJv = sdarch(1:19)//'.AEXC'
    sdarchLastJv = sdarch(1:19)//'.LAST'
    call wkvect(sdarchAexcJv, 'V V K16', 1, vk16=sdarchAexc)
    call wkvect(sdarchAinfJv, 'V V I', 4, vi=sdarchAinf)
    call wkvect(sdarchLastJv, 'V V R', 2, vr=sdarchLast)

! - Get parameters from ARCHIVAGE
    call nmcrpx(factorKeyword, stepKeyword, iocc, sdarch, jvBase)

! - Save
    sdarchLast(1) = -1.d0
    sdarchLast(2) = -1.d0
!
    call jedema()
!
end subroutine
