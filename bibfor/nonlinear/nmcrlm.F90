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
subroutine nmcrlm(listRealJv, sddisc, listInstWorkJv)
!
    implicit none
!
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/jedup1.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "event_def.h"
!
    character(len=19), intent(in) :: listRealJv, sddisc, listInstWorkJv
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Time discretization datastructure
!
! Create list of times and information vector from LISTR8_SDASTER
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  listRealJv       : name of object for list of reals
! In  listInstWorkJv   : name of working list of time
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: timeListMethod
    integer(kind=8) :: nbReal, iReal
    real(kind=8) :: timeIncrMini, timeIncr
    character(len=24) :: sddiscLinfJv
    real(kind=8), pointer :: sddiscLinf(:) => null()
    real(kind=8), pointer :: listRealVale(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sddiscLinfJv = sddisc(1:19)//'.LINF'
    timeIncrMini = r8maem()

! - Access to object
    call jeveuo(listRealJv//'.VALE', 'L', vr=listRealVale)
    call jelira(listRealJv//'.VALE', 'LONMAX', nbReal)

! - At least one step
    if (nbReal .lt. 2) then
        call utmess('F', 'DISCRETISATION_95')
    end if

! - Minimum time between two steps
    do iReal = 1, nbReal-1
        timeIncr = listRealVale(iReal+1)-listRealVale(iReal)
        timeIncrMini = min(timeIncr, timeIncrMini)
    end do

! - List must increase
    if (timeIncrMini .le. r8prem()) then
        call utmess('F', 'DISCRETISATION_87')
    end if

! - Copy list of reals in list of times
    call jedup1(listRealJv(1:19)//'.VALE', 'V', listInstWorkJv)

! - Create information vector
    call wkvect(sddiscLinfJv, 'V V R', SIZE_LLINR, vr=sddiscLinf)

! - Update information vector
    timeListMethod = 'MANUEL'
    call utdidt('E', sddisc, 'LIST', 'METHODE', valk_=timeListMethod)
    call utdidt('E', sddisc, 'LIST', 'DTMIN', valr_=timeIncrMini)
    call utdidt('E', sddisc, 'LIST', 'NBINST', vali_=nbReal)
!
end subroutine
