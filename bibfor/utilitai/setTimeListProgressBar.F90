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
subroutine setTimeListProgressBar(sddisc, numeInst, final_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
#include "asterfort/getTimeListBounds.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: numeInst
    aster_logical, optional, intent(in) :: final_
!
! --------------------------------------------------------------------------------------------------
!
! Display progress bar for time step list
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  numeInst        : index of current inst step
! In  final            : flag for final time step
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: perc, numeStore
    real(kind=8) :: timeCurr, timeInit, timeEnd, timeStore
    character(len=19) :: sdarch
    character(len=24) :: sdarchAinfJv, sdArchLastJv
    integer(kind=8), pointer :: sdarchAinf(:) => null()
    real(kind=8), pointer :: sdArchLast(:) => null()
    integer(kind=8) :: stateStoring
    aster_logical :: lStoreInitState, lStoreState
!
! --------------------------------------------------------------------------------------------------
!
    timeCurr = diinst(sddisc, numeInst)

! - Access to store datastructures
    sdarch = sddisc(1:14)//'.ARCH'
    sdarchAinfJv = sdarch(1:19)//'.AINF'
    sdarchLastJv = sdarch(1:19)//'.LAST'
    call jeveuo(sdarchAinfJv, 'L', vi=sdarchAinf)
    call jeveuo(sdarchLastJv, 'L', vr=sdarchLast)

! - Compute percentage
    call getTimeListBounds(sddisc, timeInit, timeEnd)
    if (present(final_)) then
        ASSERT(final_)
        perc = 100.
    else
        perc = int(100.d0*(timeCurr-timeInit)/(timeEnd-timeInit))
    end if

! - Get status
    lStoreInitState = ASTER_FALSE
    lStoreState = ASTER_FALSE
    stateStoring = sdarchAinf(4)
    if (stateStoring .eq. -1) then
        lStoreState = ASTER_FALSE
        lStoreInitState = ASTER_FALSE
    elseif (stateStoring .eq. 0) then
        lStoreState = ASTER_FALSE
        lStoreInitState = ASTER_TRUE
    else
        lStoreState = ASTER_TRUE
        lStoreInitState = ASTER_FALSE
    end if

! - Cas où l'on force l'archivage du dernier instant
    if (present(final_)) then
        ASSERT(final_)
        lStoreState = ASTER_TRUE
    end if

! - Get parameters of storage
    numeStore = nint(sdarchLast(1))
    timeStore = sdarchLast(2)

! - Print bar
    if (lStoreInitState) then
        call utmess('I', 'PROGRESS_4', si=perc, sr=timeCurr)
    elseif (lStoreState) then
        call utmess('I', 'PROGRESS_1', ni=2, vali=[perc, numeStore], &
                    nr=2, valr=[timeCurr, timeStore])
    else
        call utmess('I', 'PROGRESS_3', si=perc, sr=timeCurr)
    end if
!
end subroutine
