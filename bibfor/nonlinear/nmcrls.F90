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
subroutine nmcrls(sddisc, listInstJv, numeInstInit, numeInstEnd, &
                  nbInstNew, timeIncrMini)
!
    implicit none
!
#include "asterc/r8maem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utdidt.h"
#include "asterfort/wkvect.h"
!
    character(len=19), intent(in) :: sddisc, listInstJv
    integer(kind=8), intent(in) :: numeInstInit, numeInstEnd
    integer(kind=8), intent(out) :: nbInstNew
    real(kind=8), intent(out) :: timeIncrMini
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Time discretization datastructure
!
! Resize list of times
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  listInstJv       : name of JEVEUX object for list of times from INCREMENT/LIST_INST
! In  numeInstInit     : index of initial time in list of times
! In  numeInstEnd      : index of final time in list of times
! Out nbInstNew        : number of time steps in list after resize
! Out timeIncrMini     : minimum time between two steps
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: pos, iInst, nbInst
    real(kind=8) :: deltat
    real(kind=8), pointer :: listInst(:) => null()
    character(len=24) :: sddiscDitrJv
    real(kind=8), pointer :: sddiscDitr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call utdidt('L', sddisc, 'LIST', 'NBINST', vali_=nbInst)

! - Final number of time steps
    nbInstNew = (numeInstEnd-numeInstInit)+1
    ASSERT(nbInstNew .le. nbInst)

! - Acces to list of times
    call jeveuo(listInstJv, 'L', vr=listInst)

! - Create new list of time
    sddiscDitrJv = sddisc(1:19)//'.DITR'
    call wkvect(sddiscDitrJv, 'V V R', nbInstNew, vr=sddiscDitr)

! - Update new list of time
    pos = 1
    do iInst = numeInstInit, numeInstEnd
        sddiscDitr(pos) = listInst(iInst+1)
        pos = pos+1
    end do
    ASSERT(pos-1 .eq. nbInstNew)

! - New minimum time between two steps
    timeIncrMini = r8maem()
    do iInst = 1, nbInstNew-1
        deltat = sddiscDitr(iInst+1)-sddiscDitr(iInst)
        timeIncrMini = min(deltat, timeIncrMini)
    end do
!
end subroutine
