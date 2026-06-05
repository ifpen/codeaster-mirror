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
subroutine dinuar(result, sddisc, numeInst, lForceStore, &
                  numeStore, numeReuseCalc_, lStoreInitState_)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/diinst.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmcrpo.h"
#include "asterfort/rsadpa.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: result
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: numeInst
    aster_logical, intent(in) :: lForceStore
    integer(kind=8), intent(out) :: numeStore
    integer(kind=8), optional, intent(out) :: numeReuseCalc_
    aster_logical, intent(in), optional :: lStoreInitState_
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output datastructure
!
! Get storing index
!
! --------------------------------------------------------------------------------------------------
!
! In  result           : name of datastructure for results
! In  sddisc           : datastructure for time discretization
! In  numeInst         : index of current time step
! In  lForceStore      : to force storage (ex.: error)
! Out numeStore        : index to store in results
! Out numeReuseCalc    : index for reuse rsults datastructure
! In  lStoreInitState  : flag to store initial state
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdarchAinfJv
    integer(kind=8), pointer :: sdarchAinf(:) => null()
    integer(kind=8) :: numeReuseCalc, jvPara
    real(kind=8) :: timeCurr, timePrev
    aster_logical :: l_store, lStoreInitState
    character(len=19) :: sdarch
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    l_store = ASTER_FALSE
    lStoreInitState = ASTER_FALSE
    if (present(lStoreInitState_)) then
        lStoreInitState = lStoreInitState_
    end if
    numeStore = -1
    numeReuseCalc = -1

! - Acces to storing objects
    sdarch = sddisc(1:14)//'.ARCH'
    sdarchAinfJv = sdarch(1:19)//'.AINF'
    call jeveuo(sdarchAinfJv, 'E', vi=sdarchAinf)

! - Current time step
    timeCurr = 0.d0
    if (numeInst .ne. 0) then
        timeCurr = diinst(sddisc, numeInst)
    end if

! - Store or not ?
    if (lForceStore) then
        l_store = ASTER_TRUE

    else
        if (numeInst .eq. 0) then
! --------- Initial state: always
            l_store = ASTER_TRUE
        else
! --------- Other: depends on ARCHIVAGE keywords
            call nmcrpo(sdarch, numeInst, timeCurr, l_store)
        end if

    end if

! - Get storing index
    if (l_store) then
        numeStore = sdarchAinf(1)
    else
        numeStore = -1
    end if

! - REUSE for PARA_CALC table
    numeReuseCalc = sdarchAinf(3)

! - Already stored ?
    if (numeStore .ge. 2) then
        call rsadpa(result, 'L', 1, 'INST', numeStore-1, 0, sjv=jvPara)
        timePrev = zr(jvPara)
        if (timeCurr .le. timePrev) then
            numeStore = -1
            l_store = ASTER_FALSE
        end if
    end if

! - Increase storing index
    if (l_store) then
        sdarchAinf(1) = sdarchAinf(1)+1
        if (lStoreInitState) then
            sdarchAinf(4) = 0
        else
            sdarchAinf(4) = sdarchAinf(4)+1
        end if
    end if
!
    if (present(numeReuseCalc_)) then
        numeReuseCalc_ = numeReuseCalc
    end if
!
    call jedema()
!
end subroutine
