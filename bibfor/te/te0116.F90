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
subroutine te0116(option, nomte)
!
    use MetallurgyMeca_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/getAnnealingPara.h"
#include "asterfort/getAnnealingParaOnCell.h"
#include "asterfort/getAnnealingTempTrigger.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvarc.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Computing the option REST_ECRO
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8), parameter :: kspg = 1, nbVariMaxi = 30
    integer(kind=8) :: kpg, iVari, iVariEcro
    integer(kind=8) :: npg, nbVari, nbVariEcro, nbVariAnneal
    integer(kind=8) :: jvMater, jvVariOut, jvVariIn, jvTimePrev, jvTimeCurr
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: relaComp, postIncr, defoComp
    real(kind=8) :: timeCurr, timePrev
    real(kind=8) :: tauInf, x0, alpha, epsqMini, xcinMini
    real(kind=8) :: epsq, xcin
    real(kind=8) :: coefPrag, pragerTempEcroIni, prager
    real(kind=8) :: T1, T2, temp, k, dt
    real(kind=8) :: p_in, p_out, var_ecro_in, var_ecro_out, valExp
    integer(kind=8) :: iret, idxgdef
    integer(kind=8) :: variEcroIndx(nbVariMaxi), variSaveIndx(nbVariMaxi), indxEpseq
    aster_logical :: l_anneal, l_end_anneal, lTrigger
    aster_logical :: lHardIsot, lHardKine, lHardMixed
! - Protect against large overflow for exponential
    real(kind=8), parameter :: maxExp = 500.d0
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'REST_ECRO')

! - Get informations on current element
    call elrefe_info(fami=fami, npg=npg)

! - Get input fields
    call jevech('PMATERC', 'L', jvMater)
    call jevech('PVARIMR', 'L', jvVariIn)

! - Get time
    call jevech('PINSTMR', 'L', jvTimePrev)
    call jevech('PINSTPR', 'L', jvTimeCurr)
    timePrev = zr(jvTimePrev)
    timeCurr = zr(jvTimeCurr)
    dt = timeCurr-timePrev

! - Get parameters of behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)
    relaComp = compor(RELA_NAME)
    defoComp = compor(DEFO)
    read (compor(NVAR), '(I16)') nbVari
    postIncr = compor(POSTINCR)
    ASSERT(nbVari .le. nbVariMaxi)

! - Get output field
    call jevech('PVARIPR', 'E', jvVariOut)
!
    if (postIncr .eq. 'REST_ECRO') then
! ----- Type of hardening to apply annealing
        call metaAnnealGetType(relaComp, lHardIsot, lHardKine, lHardMixed, nbVariAnneal)

! ----- Get annealing parameters on cell
        call getAnnealingParaOnCell(fami, zi(jvMater), &
                                    T1, T2, &
                                    epsqMini, xcinMini)

! ----- Identify index of variables to anneal
        idxgdef = merge(6, 0, defoComp == 'GDEF_LOG')
        nbVariEcro = 0
        variEcroIndx = 0
        variSaveIndx = 0
        indxEpseq = -1
        if (relaComp .eq. 'VMIS_ISOT_LINE' .or. &
            relaComp .eq. 'VMIS_ISOT_TRAC') then
            nbVariEcro = 1
            variEcroIndx(1) = 1
            variSaveIndx(1) = nbVari-nbVariAnneal+1-idxgdef
            indxEpseq = 1
        elseif (relaComp .eq. 'VMIS_ISOT_NL') then
            nbVariEcro = 1
            variEcroIndx(1) = 1
            variSaveIndx(1) = nbVari-nbVariAnneal+1-idxgdef
            indxEpseq = 1
        elseif (relaComp .eq. 'VMIS_CINE_LINE') then
            nbVariEcro = 6
            do iVariEcro = 1, nbVariEcro
                variEcroIndx(iVariEcro) = iVariEcro
                variSaveIndx(iVariEcro) = nbVari-nbVariAnneal+iVariEcro-idxgdef
            end do
            indxEpseq = -1
        elseif (relaComp .eq. 'VMIS_ECMI_LINE') then
            nbVariEcro = 7
            variEcroIndx(1) = 1
            variSaveIndx(1) = nbVari-nbVariAnneal+1-idxgdef
            ! décalage due a l'indicateur de plasticité placé en position 2
            do iVariEcro = 2, nbVariEcro
                variEcroIndx(iVariEcro) = iVariEcro+1
                variSaveIndx(iVariEcro) = nbVari-nbVariAnneal+iVariEcro-idxgdef
            end do
            indxEpseq = 1
        elseif (relaComp .eq. 'VMIS_CIN1_CHAB') then
            nbVariEcro = 1
            variEcroIndx(1) = 1
            variSaveIndx(1) = nbVari-nbVariAnneal+1-idxgdef
            indxEpseq = 1
        elseif (relaComp .eq. 'VMIS_CIN2_CHAB') then
            nbVariEcro = 1
            variEcroIndx(1) = 1
            variSaveIndx(1) = nbVari-nbVariAnneal+1-idxgdef
            indxEpseq = 1
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Modify internal state variables
        do kpg = 1, npg
! --------- Get current temperature
            call rcvarc(' ', 'TEMP', '+', fami, kpg, kspg, temp, iret)
            ASSERT(iret .eq. 0)

! --------- Get annealing parameters
            call getAnnealingPara(fami, zi(jvMater), kpg, kspg, &
                                  T1, temp, epsqMini, &
                                  alpha, tauInf, &
                                  lHardMixed, prager, pragerTempEcroIni)

! --------- Copy internal state variables
            do iVari = 1, nbVari
                zr(jvVariOut-1+nbVari*(kpg-1)+iVari) = zr(jvVariIn-1+nbVari*(kpg-1)+iVari)
            end do

            do iVariEcro = 1, nbVariEcro
! ------------- Get hardening parameters
                epsq = 0.d0
                if (indxEpseq .ne. -1) then
                    ASSERT(lHardMixed .or. lHardIsot)
                    epsq = zr(jvVariIn-1+nbVari*(kpg-1)+indxEpseq)
                end if
                xcin = 0.d0
                if (lHardKine) then
                    ASSERT(nbVariEcro .eq. 6)
                    xcin = zr(jvVariIn-1+nbVari*(kpg-1)+variEcroIndx(iVariEcro))
                end if
                if (lHardMixed .and. iVariEcro .ge. 2) then
                    ASSERT(nbVariEcro .eq. 7)
                    xcin = zr(jvVariIn-1+nbVari*(kpg-1)+variEcroIndx(iVariEcro))
                end if

! ------------- Detect beginning of annealing
                call getAnnealingTempTrigger(temp, T1, T2, &
                                             lHardIsot, lHardKine, lHardMixed, &
                                             epsq, epsqMini, &
                                             xcin, xcinMini, &
                                             l_anneal, l_end_anneal)
! ------------- Apply annealing (or not)
                if (l_anneal) then
! ----------------- Copy internal state variables for annealing
                    p_in = zr(jvVariIn-1+nbVari*(kpg-1)+variEcroIndx(iVariEcro))
                    var_ecro_in = zr(jvVariIn-1+nbVari*(kpg-1)+variSaveIndx(iVariEcro))
                    p_out = p_in
                    var_ecro_out = var_ecro_in

! ----------------- Evaluate parameters for annealing
                    if (l_end_anneal) then
                        p_out = 0.d0
                    else
! --------------------- Get prager coefficient (for mixed hardening)
                        coefPrag = 1.d0
                        if (lHardMixed) then
                            if (iVariEcro .gt. 1) then
                                coefPrag = prager/pragerTempEcroIni
                            end if
                        end if

! --------------------- Detect trigger for annealing
                        x0 = var_ecro_in
                        if (iVariEcro .eq. 1) then
                            lTrigger = p_in .gt. abs(coefPrag*x0*tauInf)

                        else
                            lTrigger = abs(p_in) .gt. abs(coefPrag*x0*tauInf)
                        end if

                        if (lTrigger) then
                            valExp = alpha*(temp-T1)/(T2-temp)
                            if (valExp .ge. maxExp) then
                                p_out = 0.
                            else
                                k = exp(valExp)-1.d0
                                p_out = p_in*(1.d0/(1.d0+k*dt))+ &
                                        (1.d0-1.d0/(1.d0+k*dt))*(coefPrag*x0*tauInf)
                                p_out = max(p_out, 0.d0)
                            end if
                        end if
                    end if
                    zr(jvVariOut-1+nbVari*(kpg-1)+variEcroIndx(iVariEcro)) = p_out
                    zr(jvVariOut-1+nbVari*(kpg-1)+variSaveIndx(iVariEcro)) = var_ecro_out

                else
                    p_in = zr(jvVariIn-1+nbVari*(kpg-1)+variEcroIndx(iVariEcro))
                    p_out = p_in
                    var_ecro_out = p_in
                    zr(jvVariOut-1+nbVari*(kpg-1)+variEcroIndx(iVariEcro)) = p_out
                    zr(jvVariOut-1+nbVari*(kpg-1)+variSaveIndx(iVariEcro)) = var_ecro_out

                end if
            end do
        end do
    else
! ---- Internal variables OUT = Internal variables IN (no REST_ECRO)
        do kpg = 1, npg
            do iVari = 1, nbVari
                zr(jvVariOut-1+nbVari*(kpg-1)+iVari) = zr(jvVariIn-1+nbVari*(kpg-1)+iVari)
            end do
        end do
    end if
!
end subroutine
