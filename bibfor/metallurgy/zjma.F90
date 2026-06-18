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
subroutine zjma(metaSteelPara, &
                nbVari, nbVariTemper, &
                deltaTime12, &
                infoTemper, metaIn, metaOut)
!
    use Metallurgy_type
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/Metallurgy_type.h"
!
    type(META_SteelParameters), intent(in) :: metaSteelPara
    integer(kind=8), intent(in) :: nbVari, nbVariTemper
    real(kind=8), intent(in) :: deltaTime12
    real(kind=8), intent(in) :: infoTemper(NB_PARAIN_TEMPER), metaIn(nbVari)
    real(kind=8), intent(out) :: metaOut(nbVariTemper)
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY -  Tempering law for steel (bainite and martensite) phase computing
!
! --------------------------------------------------------------------------------------------------
!
! In  metaSteelPara       : parameters for metallurgy of steel
! In  nbVari              : number of internal state variables without tempering
! In  nbVariTemper        : number of internal state variables with tempering
! In  deltaTime12         : increment of time [N, N+1]
! In  infoTemper          : value parameters for tempering
! In  metaIn              : value of internal state variable without tempering
! Out metaOut             : value of internal state variable with tempering
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: kelvin = 273.0d0
    real(kind=8) :: tempPgPrev, tempPgCurr, tempPgIncr, tempTempering, tempHold
    integer(kind=8) :: cyclTherPrev, cyclTherCurr
    real(kind=8) :: ZTildeMartPrev, ZTildeBainPrev
    real(kind=8) :: ZTildeMartCurr, ZTildeBainCurr
    real(kind=8) :: tau_0_bain, tau_0_mart
    real(kind=8) :: ZBainBrut, ZMartBrut
    real(kind=8) :: ZBainRevePrev, ZMartRevePrev
    real(kind=8) :: ZBainBaseCurr, ZMartBaseCurr
    real(kind=8) :: ZBainReveCurr, ZMartReveCurr
    real(kind=8) :: ZBainTotaPrev, ZMartTotaPrev
    real(kind=8) :: bainCoefB, bainCoefN
    real(kind=8) :: martCoefB, martCoefN
    real(kind=8) :: R, delta_H, C0, deltaTimeEqui, Pa
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbVariTemper .eq. NBVARISTEELR)

! - Get temperatures
    tempPgPrev = infoTemper(1)
    tempPgCurr = metaIn(STEEL_TEMP)
    tempPgIncr = abs(tempPgCurr-tempPgPrev)

! - Get internal state variables without tempering
    ZBainBrut = metaIn(PBAINITE)
    ZMartBrut = metaIn(PMARTENS)

! - Get parameters
    cyclTherPrev = nint(infoTemper(2))
    ZBainRevePrev = infoTemper(3)
    ZMartRevePrev = infoTemper(4)

! - Compute ratio
    ZBainTotaPrev = ZBainRevePrev+ZBainBrut
    ZMartTotaPrev = ZMartRevePrev+ZMartBrut
    ZTildeMartPrev = 0.d0
    if (abs(ZMartTotaPrev) .ge. r8prem()) then
        ZTildeMartPrev = ZMartRevePrev/ZMartBrut
    end if
    ZTildeBainPrev = 0.d0
    if (abs(ZBainTotaPrev) .ge. r8prem()) then
        ZTildeBainPrev = ZBainRevePrev/ZBainBrut
    end if

! - Get material parameters
    tempTempering = metaSteelPara%temper%temp
    tempHold = metaSteelPara%temper%tempHold
    bainCoefB = metaSteelPara%temper%bainite_b
    bainCoefN = metaSteelPara%temper%bainite_n
    martCoefB = metaSteelPara%temper%martensite_b
    martCoefN = metaSteelPara%temper%martensite_n
    R = 2.0
    delta_H = 100000.0
    C0 = (log(10.d0)*R)/delta_H

! - Computes
    cyclTherCurr = cyclTherPrev
    if (abs(deltaTime12) .le. r8prem()) then
        ZBainBaseCurr = ZBainBrut
        ZMartBaseCurr = ZMartBrut
        ZBainReveCurr = ZBainRevePrev
        ZMartReveCurr = ZMartRevePrev
        cyclTherCurr = cyclTherPrev
    else
        ZTildeMartCurr = 0.d0
        ZTildeBainCurr = 0.d0
        if (tempPgCurr .gt. metaSteelPara%ac3) then
            cyclTherCurr = 0
            ZTildeMartCurr = 0.0
            ZTildeBainCurr = 0.0
        else if (tempPgCurr .gt. tempTempering) then
            if (cyclTherPrev .eq. 2) then
! ------------- Bainite
                tau_0_bain = (-1.0*log(1-ZBainRevePrev)/bainCoefB)**(1.0/bainCoefN)
                Pa = 1.0/(1.0/(tempPgCurr+kelvin)-C0*log(deltaTime12))
                deltaTimeEqui = exp((1.0/C0)*((1.0/(tempHold+kelvin))-(1.0/Pa)))
                ZTildeBainCurr = &
                    1.0-exp(-1.0*bainCoefB*(tau_0_bain+deltaTimeEqui)**bainCoefN)
! ------------- Martensite
                tau_0_mart = (-1.0*log(1-ZMartRevePrev)/martCoefB)**(1.0/martCoefN)
                Pa = 1.0/(1.0/(tempPgCurr+kelvin)-C0*log(deltaTime12))
                deltaTimeEqui = exp((1.0/C0)*((1.0/(tempHold+kelvin))-(1.0/Pa)))
                ZTildeMartCurr = &
                    1.0-exp(-1.0*martCoefB*(tau_0_mart+deltaTimeEqui)**martCoefN)
            end if
        else
            if (cyclTherPrev .eq. 2) then
                ZTildeMartCurr = 0.d0
                ZTildeBainCurr = 0.d0
                if (ZMartBrut .ge. r8prem()) then
                    ZTildeMartCurr = ZMartRevePrev/ZMartBrut
                end if
                if (ZBainBrut .ge. r8prem()) then
                    ZTildeBainCurr = ZBainRevePrev/ZBainBrut
                end if
                cyclTherCurr = cyclTherPrev
            elseif (cyclTherPrev .eq. 1) then
                ZTildeMartCurr = 0.d0
                ZTildeBainCurr = 0.d0
                if (tempPgIncr .ge. r8prem()) then
                    cyclTherCurr = 2
                end if
            else if (cyclTherPrev .eq. 0) then
                ZTildeMartCurr = 0.d0
                ZTildeBainCurr = 0.d0
                if (ZBainBrut .ge. r8prem() .or. ZMartBrut .ge. r8prem()) then
                    cyclTherCurr = 1
                end if
            else
                WRITE (6, *) "cyclTherPrev: ", cyclTherPrev
                ASSERT(ASTER_FALSE)
            end if
        end if
        ZMartReveCurr = ZMartBrut*ZTildeMartCurr
        ZBainReveCurr = ZBainBrut*ZTildeBainCurr
        ZMartBaseCurr = ZMartBrut-ZMartReveCurr
        ZBainBaseCurr = ZBainBrut-ZBainReveCurr
    end if

! - Update internal state variables
    metaOut(PRFERRITE) = metaIn(PFERRITE)
    metaOut(PRPERLITE) = metaIn(PPERLITE)
    metaOut(PRBAINITEB) = ZBainBaseCurr
    metaOut(PRMARTENSB) = ZMartBaseCurr
    metaOut(PRBAINITER) = ZBainReveCurr
    metaOut(PRMARTENSR) = ZMartReveCurr
    metaOut(PRAUSTENITE) = metaIn(PAUSTENITE)
    metaOut(PRSUMCOLD) = 1.d0-metaIn(PAUSTENITE)
    metaOut(SIZE_GRAINR) = metaIn(SIZE_GRAIN)
    metaOut(STEEL_TEMPR) = metaIn(STEEL_TEMP)
    metaOut(TEMP_MARTENSITER) = metaIn(TEMP_MARTENSITE)
    metaOut(THER_CYCL) = cyclTherCurr
!
end subroutine
