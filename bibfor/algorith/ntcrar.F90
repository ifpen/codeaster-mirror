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
subroutine ntcrar(resultZ, sddiscZ, sdarchZ, lReuse)
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmarex.h"
#include "asterfort/nmarnr.h"
#include "asterfort/nmarpr.h"
#include "asterfort/nmcrpx.h"
#include "asterfort/nmdide.h"
#include "asterfort/wkvect.h"
!
    character(len=*), intent(in) :: resultZ, sddiscZ, sdarchZ
    aster_logical, intent(in) :: lReuse
!
! --------------------------------------------------------------------------------------------------
!
! THER_NON_LINE - Storing management
!
! Create storing objects
!
! --------------------------------------------------------------------------------------------------
!
! In  result           : name of datastructure for results
! In  sddisc           : datastructure for time discretization
! In  sdarch           : datastructure for storing
! In  l_reuse          : .true. if reuse results datastructure
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8), parameter :: iocc = 1
    character(len=16), parameter :: factorKeyword = 'ARCHIVAGE', keywStep = 'PAS_ARCH'
    character(len=1), parameter :: jvBase = 'V'
    character(len=19) :: sddisc, sdarch
    character(len=8) :: result
    integer(kind=8) :: nbFactorKeyword
    integer(kind=8) :: lastIndex, numeReuseCalc, numeReuse
    integer(kind=8) :: numeStoring
    character(len=24) :: sdarchAinfJv, sdarchLastJv
    integer(kind=8), pointer :: sdarchAinf(:) => null()
    real(kind=8), pointer :: sdarchLast(:) => null()
    real(kind=8) :: lastTime
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<THERNONLINE> ... CREATION SD ARCHIVAGE'
    end if

! - Initializations
    result = resultZ
    sdarch = sdarchZ
    sddisc = sddiscZ
    numeStoring = -1
    numeReuse = -1
    numeReuseCalc = -1
    call getfac(factorKeyword, nbFactorKeyword)
    ASSERT(nbFactorKeyword .le. 1)

! - Get last time in result datastructure if initial state given
    call nmdide(lReuse, result, lastIndex, lastTime)

! - Get parameters from ARCHIVAGE
    call nmcrpx(factorKeyword, keywStep, iocc, sdarch, jvBase)

! - List of CHAM_EXCLU
    call nmarex(factorKeyword, sdarch)

! - Get index to save first time to store
    call nmarpr(result, sddisc, lReuse, lastIndex, lastTime, numeStoring)

! - Get reuse index from TABLE OBSERVATION
    call nmarnr(result, 'OBSERVATION', numeReuse)

! - Create datastructures
    sdarchAinfJv = sdarch(1:19)//'.AINF'
    sdarchLastJv = sdarch(1:19)//'.LAST'
    call wkvect(sdarchAinfJv, 'V V I', 4, vi=sdarchAinf)
    call wkvect(sdarchLastJv, 'V V R', 2, vr=sdarchLast)

! - Save
    ASSERT(numeStoring .ge. 0)
    ASSERT(numeReuse .ge. 0)
    sdarchAinf(1) = numeStoring
    sdarchAinf(2) = numeReuse
    sdarchAinf(3) = numeReuseCalc
    sdarchAinf(4) = -1
    sdarchLast(1) = -1.d0
    sdarchLast(2) = -1.d0
!
    call jedema()
!
end subroutine
