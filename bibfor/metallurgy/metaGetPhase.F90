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
subroutine metaGetPhase(fami, poum, kpg, ksp, &
                        metaType, nbPhase, &
                        phase_, zSumCold_, zHot_, tolePhase_)
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/rcvarc.h"
!
    character(len=*), intent(in) :: fami
    character(len=1), intent(in) :: poum
    integer(kind=8), intent(in) :: kpg, ksp
    integer(kind=8), intent(in) :: metaType, nbPhase
    real(kind=8), optional, intent(out) :: phase_(nbPhase)
    real(kind=8), optional, intent(out) :: zSumCold_
    real(kind=8), optional, intent(out) :: zHot_
    real(kind=8), optional, intent(in) :: tolePhase_
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Get phases
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  poum         : '-' or '+' for parameters evaluation (previous or current temperature)
! In  kpg          : current point gauss
! In  ksp          : current "sous-point" gauss
! In  metaType     : type of metallurgy
! In  nbPhase      : total number of phase (cold and hot)
! Out phase        : phases
! Out zSumCold     : sum of cold phases
! Out zHot         : hot phase
! In  tolePhase    : tolerance to project phase proportion on boundary
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: steel(NBPHASESTEEL) = (/'FERRITE', 'PERLITE', &
                                                           'BAINITE', 'MARTENS', &
                                                           'AUSTENI', 'COLDSUM'/)
    character(len=8), parameter :: steelTemper(NBPHASESTEELR) = (/'FERRITE', 'PERLITE', &
                                                                  'BAINITB', 'BAINITR', &
                                                                  'MARTENB', 'MARTENR', &
                                                                  'AUSTENI', 'COLDSUM'/)
    character(len=8), parameter :: zirc(3) = (/'ALPHPUR ', 'ALPHBETA', &
                                               'BETA    '/)
    integer(kind=8) :: iPhaseCold, iPhase, iret
    real(kind=8) :: zSumCold, zHot, phase(META_MECA_NBPHASE_MAXI), tolePhase
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbPhase .le. META_MECA_NBPHASE_MAXI)
    ASSERT(nbPhase .gt. 0)
    phase = 0.d0
    if (present(tolePhase_)) then
        tolePhase = tolePhase_
    else
        tolePhase = r8prem()
    end if

! - Get phases
    do iPhase = 1, nbPhase
        if (metaType .eq. META_STEEL) then
            call rcvarc('F', steel(iPhase), poum, fami, kpg, &
                        ksp, phase(iPhase), iret)
            if (iret .eq. 1) then
                phase(iPhase) = 0.d0
            end if

        elseif (metaType .eq. META_ZIRC) then
            call rcvarc('F', zirc(iPhase), poum, fami, kpg, &
                        ksp, phase(iPhase), iret)
            if (iret .eq. 1) then
                phase(iPhase) = 0.d0
            end if

        else if (metaType .eq. META_STEEL_TEMPER) then
            call rcvarc('F', steelTemper(iPhase), poum, fami, kpg, &
                        ksp, phase(iPhase), iret)
            if (iret .eq. 1) then
                phase(iPhase) = 0.d0
            end if

        else
            ASSERT(ASTER_FALSE)
        end if
    end do

! - Get sum of cold phases
    zSumCold = 0.d0
    if (metaType .eq. META_STEEL) then
        call rcvarc('F', steel(NBPHASESTEEL), poum, fami, kpg, &
                    ksp, zSumCold, iret)
        if (iret .eq. 1) then
            zSumCold = 0.d0
        end if

    elseif (metaType .eq. META_ZIRC) then
        do iPhaseCold = 1, nbPhase-1
            zSumCold = zSumCold+phase(iPhaseCold)
        end do
        if (zSumCold .le. tolePhase) then
            zSumCold = 0.d0
        end if
        if (zSumCold .ge. 1.d0) then
            zSumCold = 1.d0
        end if

    else if (metaType .eq. META_STEEL_TEMPER) then
        call rcvarc('F', steelTemper(NBPHASESTEELR), poum, fami, kpg, &
                    ksp, zSumCold, iret)
        if (iret .eq. 1) then
            zSumCold = 0.d0
        end if
    end if

! - Outputs
    zHot = phase(nbPhase)
    if (present(phase_)) then
        phase_(1:nbPhase) = phase(1:nbPhase)
    end if
    if (present(zSumCold_)) then
        zSumCold_ = zSumCold
    end if
    if (present(zHot_)) then
        zHot_ = zHot
    end if
!
end subroutine
