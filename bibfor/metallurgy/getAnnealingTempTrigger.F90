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
subroutine getAnnealingTempTrigger(temp, T1, T2, &
                                   lHardIsot, lHardKine, lHardMixed, &
                                   epsq, epsqMini, &
                                   xcin, xcinMini, &
                                   l_anneal, l_end_anneal)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterf_types.h"
!
    real(kind=8), intent(in) :: temp, T1, T2
    aster_logical, intent(in) :: lHardIsot, lHardKine, lHardMixed
    real(kind=8), intent(in) :: epsq, epsqMini
    real(kind=8), intent(in) :: xcin, xcinMini
    aster_logical, intent(out) :: l_anneal, l_end_anneal
!
! --------------------------------------------------------------------------------------------------
!
! - Protect agains division by zero
    real(kind=8), parameter :: toleTemp = 1.d0
    aster_logical :: lHardMini
!
! --------------------------------------------------------------------------------------------------
!

! - Detect minimum hardening
    lHardMini = ASTER_FALSE
    if (lHardIsot) then
        lHardMini = (abs(epsq) .gt. epsqMini)
    elseif (lHardKine) then
        lHardMini = (abs(xcin) .gt. xcinMini)
    elseif (lHardMixed) then
        lHardMini = (abs(epsq) .gt. epsqMini)
    else
        ASSERT(ASTER_FALSE)
    end if

    l_anneal = ((temp .gt. T1) .and. lHardMini)
! Je suis très proche de la borne fin de la restauration
! Pour éviter la division par zéro je décide
    l_end_anneal = ASTER_FALSE
    if ((abs(temp-T2) .le. toleTemp) .or. ((temp-T2) .gt. toleTemp)) then
        l_end_anneal = ASTER_TRUE
    end if
!
end subroutine
