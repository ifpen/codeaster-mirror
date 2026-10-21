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
subroutine te0064(option, nomte)
!
    use Metallurgy_type
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/nzcomp_prep.h"
#include "asterfort/nzcomp.h"
#include "asterfort/nzcompTemper.h"
#include "asterfort/tecach.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: THERMIQUE - 3D
! Option: META_ELNO
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: deltaTime01, deltaTime12, time2
    real(kind=8) :: temp1, tempInit, temp2
    integer(kind=8) :: nbNode, iNode, iMeta
    integer(kind=8) :: jvMater, itempe, itempi, jvTempInit, jvTime
    character(len=16) :: metaName, metaNameTemper
    integer(kind=8) :: numeComp, numeCompTemper
    integer(kind=8) :: nbVari, nbPhase
    integer(kind=8) :: nbVariIn, nbVariOut
    integer(kind=8) :: nbVariTemper, nbPhaseTemper, nbVariPrev
    integer(kind=8) :: nbMetaIn, nbMetaOut
    integer(kind=8) :: jvPhaseIn, jvPhaseOut, jvPhasePrev
    integer(kind=8) :: jvMaterCode
    integer(kind=8) :: itab(7), iret, jvComporMetaTemper, jvComporMeta
    aster_logical :: hasTemperLaw, hasMetaLaw, prevMetaIsTemper
    type(META_MaterialParameters) :: metaPara, metaParaTemper
    real(kind=8) :: infoTemper(NB_PARAIN_TEMPER)
    real(kind=8) :: metaIn(META_META_NBPHASE_MAXI), metaOut(META_META_NBPHASE_MAXI)
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nbNode)
    ASSERT(nbNode .le. MT_NNOMAX3D)

! - Input fields
    call jevech('PMATERC', 'L', jvMater)
    jvMaterCode = zi(jvMater)
    call jevech('PTEMPER', 'L', itempe)
    call jevech('PTEMPIR', 'L', itempi)
    call jevech('PTIMMTR', 'L', jvTime)

! - Get input and output fields for META_ELNO
    call tecach("ONN", "PPHASIN", 'L', iret, nval=7, itab=itab)
    ASSERT(iret .eq. 0)
    jvPhaseIn = itab(1)
    nbVariIn = itab(6)
    ASSERT(nbVariIn .le. META_META_NBPHASE_MAXI)
    ASSERT(nbVariIn .gt. 0)
    call tecach("ONN", "PPHASOUT", 'E', iret, nval=7, itab=itab)
    ASSERT(iret .eq. 0)
    jvPhaseOut = itab(1)
    nbVariOut = itab(6)
    ASSERT(nbVariOut .le. META_META_NBPHASE_MAXI)
    ASSERT(nbVariOut .gt. 0)

! - Get parameters from metallurgy behaviour (without tempering)
    call tecach("ONN", "PCOMPME", 'L', iret, nval=7, itab=itab)
    hasMetaLaw = ASTER_FALSE
    jvComporMeta = 0
    metaName = " "
    numeComp = 0
    nbPhase = 0
    nbVari = 0
    if (iret .eq. 0) then
        hasMetaLaw = ASTER_TRUE
        jvComporMeta = itab(1)
        metaName = zk16(jvComporMeta-1+ZMETATYPE)
        read (zk16(jvComporMeta-1+ZNUMECOMP), '(I16)') numeComp
        read (zk16(jvComporMeta-1+ZNBPHASE), '(I16)') nbPhase
        read (zk16(jvComporMeta-1+ZNBVARI), '(I16)') nbVari
        ASSERT(nbVari .le. META_META_NBPHASE_MAXI)
    end if

! - Get parameters from metallurgy behaviour (with tempering)
    call tecach("ONN", "PCOMPMT", 'L', iret, nval=7, itab=itab)
    hasTemperLaw = ASTER_FALSE
    jvComporMetaTemper = 0
    metaNameTemper = " "
    numeCompTemper = 0
    nbPhaseTemper = 0
    nbVariTemper = 0
    if (iret .eq. 0) then
        hasTemperLaw = ASTER_TRUE
        jvComporMetaTemper = itab(1)
        metaNameTemper = zk16(jvComporMetaTemper-1+ZMETATYPE)
        read (zk16(jvComporMetaTemper-1+ZNUMECOMP), '(I16)') numeCompTemper
        read (zk16(jvComporMetaTemper-1+ZNBPHASE), '(I16)') nbPhaseTemper
        read (zk16(jvComporMetaTemper-1+ZNBVARI), '(I16)') nbVariTemper
        ASSERT(metaNameTemper .eq. 'ACIER_REVENU')
        ASSERT(nbVariTemper .le. NBVARISTEELR)
        nbVari = NBVARISTEEL
    end if

! - Specific input/output fields
    jvTempInit = 0
    if (hasMetaLaw) then
        call jevech('PTEMPAR', 'L', jvTempInit)
    end if
    jvPhasePrev = 0
    nbVariPrev = 0
    prevMetaIsTemper = ASTER_FALSE
    if (hasTemperLaw) then
        call tecach("ONN", "PPHASEP", 'L', iret, nval=7, itab=itab)
        ASSERT(iret .eq. 0)
        jvPhasePrev = itab(1)
        nbVariPrev = itab(6)
        ASSERT(nbVariPrev .le. META_META_NBPHASE_MAXI)
        ASSERT(nbVariOut .eq. nbVariTemper)
        if (nbVariPrev .eq. NBVARISTEELR) then
            prevMetaIsTemper = ASTER_TRUE
        elseif (nbVariPrev .eq. NBVARISTEEL) then
            prevMetaIsTemper = ASTER_FALSE
        else
            ASSERT(ASTER_FALSE)
        end if
    end if

! - Get material and TRC parameters
    if (hasMetaLaw) then
        call nzcomp_prep(jvMaterCode, metaName, metaPara)
    end if
    if (hasTemperLaw) then
        call nzcomp_prep(jvMaterCode, metaNameTemper, metaParaTemper)
    end if

! - Time parameters
    deltaTime01 = zr(jvTime+1)
    deltaTime12 = zr(jvTime+2)
    time2 = zr(jvTime)+deltaTime12

! - Select sizes
    nbMetaIn = nbVari
    if (hasTemperLaw) then
        nbMetaOut = nbVariTemper
    else
        nbMetaOut = nbVari
    end if

! - Loop on nodes
    do iNode = 1, nbNode
! ----- Temperatures: 1 - 2
        temp1 = zr(itempe+iNode-1)
        temp2 = zr(itempi+iNode-1)

! ----- General switches
        if (hasMetaLaw) then
            metaIn = 0.d0
            do iMeta = 1, nbMetaIn
                metaIn(iMeta) = zr(jvPhaseIn+nbVariIn*(iNode-1)+iMeta-1)
            end do
            metaOut = 0.d0
            tempInit = zr(jvTempInit+iNode-1)
            call nzcomp(jvMaterCode, metaPara, numeComp, &
                        nbPhase, nbVari, &
                        deltaTime01, deltaTime12, time2, &
                        tempInit, temp1, temp2, &
                        metaIn, metaOut)
        end if

        if (hasTemperLaw) then
            metaIn = 0.d0
            do iMeta = 1, nbMetaIn
                metaIn(iMeta) = zr(jvPhaseIn+nbVariIn*(iNode-1)+iMeta-1)
            end do
            metaOut = 0.d0
            infoTemper = 0.d0
            if (prevMetaIsTemper) then
                infoTemper(1) = zr(jvPhasePrev+nbVariPrev*(iNode-1)-1+STEEL_TEMPR)
                infoTemper(2) = zr(jvPhasePrev+nbVariPrev*(iNode-1)-1+THER_CYCL)
                infoTemper(3) = zr(jvPhasePrev+nbVariPrev*(iNode-1)-1+PRBAINITER)
                infoTemper(4) = zr(jvPhasePrev+nbVariPrev*(iNode-1)-1+PRMARTENSR)
                infoTemper(5) = zr(jvPhasePrev+nbVariPrev*(iNode-1)-1+PRBAINITEB)
                infoTemper(6) = zr(jvPhasePrev+nbVariPrev*(iNode-1)-1+PRMARTENSB)
            else
                infoTemper(1) = zr(jvPhasePrev+nbVariPrev*(iNode-1)-1+STEEL_TEMP)
            end if
            call nzcompTemper(metaParaTemper, numeCompTemper, &
                              nbVari, nbVariTemper, &
                              deltaTime12, &
                              infoTemper, metaIn, metaOut)
        end if
        do iMeta = 1, nbMetaOut
            zr(jvPhaseOut+nbVariOut*(iNode-1)+iMeta-1) = metaOut(iMeta)
        end do

    end do
!
end subroutine
