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
! aslint: disable=C1505,W1504
!
subroutine nzedga(option, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, carcri, &
                  timePrev, timeCurr, &
                  neps, epsm, deps, &
                  nsig, sigm, &
                  nvi, vim, &
                  sigp, vip, ndsde, dsidep, &
                  codret)
!
    use MetallurgyMeca_module, only: metaGetType, metaGetParaElas, metaGetMechanism, &
                                     metaGetParaMixture, metaGetParaAnneal, metaGetParaVisc, &
                                     metaGetParaPlasTransf, metaGetParaHardLine, metaGetParaHardTrac
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/nzcalc.h"
#include "asterfort/rcvarc.h"
#include "asterfort/verift.h"
!
    character(len=16), intent(in) :: option
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, ndim, jvMaterCode
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: timePrev, timeCurr
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps), deps(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    real(kind=8), intent(out) :: sigp(nsig)
    real(kind=8), intent(out) :: vip(nvi)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(out) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                                        merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Comportment
!
! META_P_I* / META_V_I* for small strains and zirconium metallurgy
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  IMAT    : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  carcri  : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  EPSM    : DEFORMATIONS A L'INSTANT DU CALCUL PRECEDENT
! IN  DEPS    : INCREMENT DE DEFORMATION
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MATRICE CARREE
!     IRET    : CODE RETOUR DE LA RESOLUTION DE L'EQUATION SCALAIRE
!               (NZCALC)
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ECHEC
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX YY ZZ XY XZ YZ
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(kind=8), parameter :: hardCoef = 1.d0
    integer(kind=8), parameter :: viIndxIPlas = 5, viIndxSigy = 4
    integer(kind=8) :: metaType
    integer(kind=8) :: nbPhase, nbPhaseCold, iPhaseCold, iPhaseHot, iPhase
    real(kind=8) :: phaseCurr(PZIRC_NB), phasePrev(PZIRC_NB)
    real(kind=8) :: phaseIncr(PZIRC_NB), zSumCold
    real(kind=8) :: dz1(PZIRC_NB-1), dz2(PZIRC_NB-1)
    real(kind=8) :: epseq(PZIRC_NB), epseqPrev(PZIRC_NB)
    real(kind=8) :: epseqCurr(PZIRC_NB)
    real(kind=8) :: sy(PZIRC_NB)
    real(kind=8) :: h(PZIRC_NB), hplus(PZIRC_NB)
    real(kind=8) :: r(PZIRC_NB)
    real(kind=8) :: annealTheta(2*(PZIRC_NB-1))
    real(kind=8) :: eta(PZIRC_NB), n(PZIRC_NB)
    real(kind=8) :: unsurn(PZIRC_NB), c(PZIRC_NB)
    real(kind=8) :: m(PZIRC_NB)
    real(kind=8) :: kpt(PZIRC_NB-1), fpt(PZIRC_NB-1)
    real(kind=8) :: n0(PZIRC_NB)
    real(kind=8) :: sigyPrev, sigyCurr
    integer(kind=8) :: iValeMaxi, nbValeMaxi
    integer(kind=8) :: ndimsi, i, j, mode, iretTemp
    real(kind=8) :: temp, timeIncr
    real(kind=8) :: epsth, e, deuxmu, deumum, troisk
    real(kind=8) :: fmix, hmoy, rmoy
    real(kind=8) :: cmoy, mmoy, cr
    real(kind=8) :: annealFunc, vimoy, ds
    real(kind=8) :: trans
    real(kind=8) :: trepsm, trdeps, trsigm, trsigp
    real(kind=8) :: dvdeps(6), dvsigm(6), dvsigp(6)
    real(kind=8) :: sigel(6), sig0(6), sieleq, sigeps
    real(kind=8) :: dp, seuil
    real(kind=8) :: coef1, coef2, coef3, dv, b
    real(kind=8) :: precr
    character(len=1) :: poum
    integer(kind=8) :: test
    real(kind=8) :: indicPlasPrev, indicPlasCurr, indicPlas
    aster_logical :: lMatr, lVari, l_temp
    aster_logical :: l_visc, l_anneal, lPlasTran, l_hard_isotline, l_hard_isotnlin
    character(len=16) :: metaRela, metaGlob
!
! --------------------------------------------------------------------------------------------------
!
    ndimsi = 2*ndim
    codret = 0
    lMatr = L_MATR(option)
    lVari = L_VARI(option)
    timeIncr = timeCurr-timePrev
    precr = r8prem()

! - Behaviour in kit
    metaRela = compor(META_RELA)
    metaGlob = compor(META_GLOB)

! - Get metallurgy type
    call metaGetType(metaType, nbPhase)
    ASSERT(metaType .eq. META_ZIRC)
    ASSERT(nbPhase .eq. PZIRC_NB)
    nbPhaseCold = nbPhase-1
    iPhaseHot = nbPhase

! - Get phases
    if (lVari) then
        poum = '+'
        call metaGetPhase(fami, '+', kpg, ksp, metaType, &
                          nbPhase, phaseCurr, zSumCold_=zSumCold)
        call metaGetPhase(fami, '-', kpg, ksp, metaType, &
                          nbPhase, phasePrev)
    else
        poum = '-'
        call metaGetPhase(fami, '-', kpg, ksp, metaType, &
                          nbPhase, phaseCurr, zSumCold_=zSumCold)
    end if
    phaseIncr(1:nbPhaseCold) = phaseCurr(1:nbPhaseCold)-phasePrev(1:nbPhaseCold)

! - Compute thermic strain
    call verift(fami, kpg, ksp, poum, jvMaterCode, epsth_meta_=epsth)
    call rcvarc(' ', 'TEMP', poum, fami, kpg, ksp, temp, iretTemp)
    l_temp = iretTemp .eq. 0

! - Get active mechanisms in behaviour
    call metaGetMechanism(metaRela, metaGlob, &
                          lVisc_=l_visc, &
                          lAnneal_=l_anneal, &
                          lPlasTran_=lPlasTran, &
                          lHardIsotLine_=l_hard_isotline, &
                          lHardIsotNLin_=l_hard_isotnlin)

! - Get elastic parameters
    call metaGetParaElas(poum, fami, kpg, ksp, jvMaterCode, &
                         e_=e, deuxmu_=deuxmu, troisk_=troisk, &
                         deuxmum_=deumum)

! - Get internal state variables
    epseqPrev(1:nbPhase) = vim(1:nbPhase)
    indicPlasPrev = vim(viIndxIPlas)
    sigyPrev = vim(viIndxSigy)

! - Mixture law (yield limit)
    call metaGetParaMixture(poum, fami, kpg, ksp, jvMaterCode, &
                            l_visc, metaType, nbPhase, zSumCold, &
                            fmix, sy)

! - Others parameters
    trans = 0.d0
    indicPlas = indicPlasPrev

! - Integration of behaviour
    if (lVari) then
! ----- Get parameters for annealing (*_RE)
        annealTheta = 1.d0
        if (l_anneal) then
            call metaGetParaAnneal(poum, fami, kpg, ksp, jvMaterCode, &
                                   metaType, nbPhaseCold, &
                                   annealTheta)
        end if

! ----- Parameters for viscosity
        if (l_visc) then
            call metaGetParaVisc(poum, fami, kpg, ksp, jvMaterCode, &
                                 metaType, nbPhase, &
                                 eta, n, unsurn, &
                                 c, m)
        else
            eta = 0.d0
            n = 20.d0
            unsurn = 1.d0
            c = 0.d0
            m = 20.d0
        end if

! ----- Get positive part of +increment of phases (dz1) and -increment of phases (dz2)
        do iPhaseCold = 1, nbPhaseCold
            if (phaseIncr(iPhaseCold) .ge. 0.d0) then
                dz1(iPhaseCold) = phaseIncr(iPhaseCold)
                dz2(iPhaseCold) = 0.d0
            else
                dz1(iPhaseCold) = 0.d0
                dz2(iPhaseCold) = -phaseIncr(iPhaseCold)
            end if
        end do

! ----- Metallurgical annealing for hot phase
        if (phaseCurr(iPhaseHot) .gt. 0.d0) then
            annealFunc = 0.d0
            do iPhaseCold = 1, nbPhaseCold
                annealFunc = annealFunc+ &
                             dz2(iPhaseCold)* &
                             (annealTheta(nbPhaseCold+iPhaseCold)*epseqPrev(iPhaseCold)- &
                              epseqPrev(iPhaseHot))/phaseCurr(iPhaseHot)
            end do
            epseq(iPhaseHot) = epseqPrev(iPhaseHot)+annealFunc
            vimoy = phaseCurr(iPhaseHot)*epseq(iPhaseHot)
        else
            epseq(iPhaseHot) = 0.d0
            vimoy = 0.d0
        end if

! ----- Metallurgical annealing for cold phases
        do iPhaseCold = 1, nbPhaseCold
            if (phaseCurr(iPhaseCold) .gt. 0.d0) then
                annealFunc = dz1(iPhaseCold)* &
                             (annealTheta(iPhaseCold)*epseqPrev(iPhaseHot)- &
                              epseqPrev(iPhaseCold))/phaseCurr(iPhaseCold)
                epseq(iPhaseCold) = epseqPrev(iPhaseCold)+annealFunc
                vimoy = vimoy+phaseCurr(iPhaseCold)*epseq(iPhaseCold)
            else
                epseq(iPhaseCold) = 0.d0
            end if
        end do

! ----- Viscuous annealing
        cmoy = 0.d0
        mmoy = 0.d0
        do iPhase = 1, nbPhase
            cmoy = cmoy+phaseCurr(iPhase)*c(iPhase)
            mmoy = mmoy+phaseCurr(iPhase)*m(iPhase)
        end do
        cr = cmoy*vimoy
        if (cr .le. 0.d0) then
            ds = 0.d0
        else
            ds = timeIncr*(cr**mmoy)
        end if

        do iPhase = 1, nbPhase
            if (phaseCurr(iPhase) .gt. 0.d0) then
                epseq(iPhase) = epseq(iPhase)-ds
                if (epseq(iPhase) .le. 0.d0) then
                    epseq(iPhase) = 0.d0
                end if
            end if
        end do

! ----- Parameters for plasticity of tranformation
        trans = 0.d0
        if (lPlasTran) then
            call metaGetParaPlasTransf('+', fami, 1, 1, jvMaterCode, &
                                       metaType, nbPhase, nbPhaseCold, &
                                       phaseIncr, zSumCold, &
                                       kpt, fpt)
            do iPhaseCold = 1, nbPhaseCold
                if (phaseIncr(iPhaseCold) .gt. 0.d0) then
                    trans = trans+kpt(iPhaseCold)*fpt(iPhaseCold)*phaseIncr(iPhaseCold)
                end if
            end do
        end if

    else
        epseq(1:nbPhase) = epseqPrev(1:nbPhase)

    end if

! - Get current elasticity yield and current hardening slope (linear hardening)
    nbValeMaxi = 0
    if (l_hard_isotline) then
        call metaGetParaHardLine(poum, fami, kpg, ksp, jvMaterCode, &
                                 metaType, nbPhase, &
                                 e, hardCoef, h)
        do iPhase = 1, nbPhase
            r(iPhase) = h(iPhase)*epseq(iPhase)+sy(iPhase)
        end do
    end if

! - Get current elasticity yield and current hardening slope (non-linear hardening)
    if (l_hard_isotnlin) then
        call metaGetParaHardTrac(jvMaterCode, metaType, nbPhase, &
                                 l_temp, temp, &
                                 epseq, h, r, nbValeMaxi)
        do iPhase = 1, nbPhase
            r(iPhase) = r(iPhase)+sy(iPhase)
        end do
    end if

! - Apply mixture law on hardening parameters
    rmoy = 0.d0
    hmoy = 0.d0
    if (zSumCold .gt. 0.d0) then
        do iPhaseCold = 1, nbPhaseCold
            rmoy = rmoy+phaseCurr(iPhaseCold)*r(iPhaseCold)
            hmoy = hmoy+phaseCurr(iPhaseCold)*h(iPhaseCold)
        end do
        rmoy = rmoy/zSumCold
        hmoy = hmoy/zSumCold
    end if
    rmoy = (1.d0-fmix)*r(iPhaseHot)+fmix*rmoy
    hmoy = (1.d0-fmix)*h(iPhaseHot)+fmix*hmoy

! - Current strain state
    trdeps = (deps(1)+deps(2)+deps(3))/3.d0
    trepsm = (epsm(1)+epsm(2)+epsm(3))/3.d0

! - Current stress state
    trsigm = (sigm(1)+sigm(2)+sigm(3))/3.d0
    trsigp = troisk*(trepsm+trdeps)-troisk*epsth
    dvsigm = 0.d0
    dvdeps = 0.d0
    do i = 1, ndimsi
        dvdeps(i) = deps(i)-trdeps*kron(i)
        dvsigm(i) = sigm(i)-trsigm*kron(i)
    end do
    sieleq = 0.d0
    sigel = 0.d0
    do i = 1, ndimsi
        sigel(i) = deuxmu*dvsigm(i)/deumum+deuxmu*dvdeps(i)
        sieleq = sieleq+sigel(i)**2
    end do
    sieleq = sqrt(1.5d0*sieleq)
    if (sieleq .gt. 0.d0) then
        do i = 1, ndimsi
            sig0(i) = sigel(i)/sieleq
        end do
    else
        do i = 1, ndimsi
            sig0(i) = 0.d0
        end do
    end if

! - Integration of behaviour
    indicPlasCurr = indicPlasPrev
    dp = 0.d0
    if (lVari) then
        indicPlasCurr = 0.d0
        seuil = sieleq-(1.5d0*deuxmu*trans+1.d0)*rmoy

! ----- Compute DP (plastic multiplier)
        if (seuil .lt. 0.d0) then
            indicPlasCurr = 0.d0
            dp = 0.d0
        else
            indicPlasCurr = 1.d0
            call nzcalc(carcri, nbPhase, phaseCurr, zSumCold, &
                        fmix, seuil, timeIncr, trans, &
                        hmoy, deuxmu, eta, unsurn, &
                        dp, codret)
            if (codret .eq. 1) goto 999
!
! DANS LE CAS NON LINEAIRE
! VERIFICATION QU ON EST DANS LE BON INTERVALLE
!
            if (l_hard_isotnlin) then

                do iValeMaxi = 1, nbValeMaxi
                    test = 0

! ----------------- Update cumulated plastic strain and hardening slope
                    epseqCurr(1:nbPhase) = epseq(1:nbPhase)+dp
                    hplus(1:nbPhase) = h(1:nbPhase)

! ----------------- Get point on all traction curves
                    call metaGetParaHardTrac(jvMaterCode, metaType, nbPhase, &
                                             l_temp, temp, &
                                             epseqCurr, h, r)

                    do iPhase = 1, nbPhase
                        if (phaseCurr(iPhase) .gt. 0.d0) then
                            r(iPhase) = r(iPhase)+sy(iPhase)
                            if (abs(h(iPhase)-hplus(iPhase)) .gt. precr) then
                                test = 1
                            end if
                        end if
                    end do
                    if (test .eq. 0) exit
                    hmoy = 0.d0
                    rmoy = 0.d0
                    if (zSumCold .gt. 0.d0) then
                        do iPhaseCold = 1, nbPhaseCold
                            if (phaseCurr(iPhaseCold) .gt. 0.d0) then
                                rmoy = rmoy+phaseCurr(iPhaseCold)*(r(iPhaseCold)-h(iPhaseCold)*dp)
                                hmoy = hmoy+phaseCurr(iPhaseCold)*h(iPhaseCold)
                            end if
                        end do
                        rmoy = fmix*rmoy/zSumCold
                        hmoy = fmix*hmoy/zSumCold
                    end if
                    if (phaseCurr(iPhaseHot) .gt. 0.d0) then
                        rmoy = (1.d0-fmix)*(r(iPhaseHot)-h(iPhaseHot)*dp)+rmoy
                        hmoy = (1.d0-fmix)*h(iPhaseHot)+hmoy
                    end if
                    seuil = sieleq-(1.5d0*deuxmu*trans+1.d0)*rmoy
                    call nzcalc(carcri, nbPhase, phaseCurr, zSumCold, &
                                fmix, seuil, timeIncr, trans, &
                                hmoy, deuxmu, eta, unsurn, &
                                dp, codret)
                    if (codret .eq. 1) goto 999
                end do
            end if
        end if

! ----- Compute stresses
        sigp(1:2*ndim) = 0.d0
        do i = 1, ndimsi
            dvsigp(i) = sigel(i)-1.5d0*deuxmu*dp*sig0(i)
            dvsigp(i) = dvsigp(i)/(1.5d0*deuxmu*trans+1.d0)
            sigp(i) = dvsigp(i)+trsigp*kron(i)
        end do

! ----- Update cumumated plastic strain
        do iPhase = 1, nbPhase
            if (phaseCurr(iPhase) .gt. 0.d0) then
                epseqCurr(iPhase) = epseq(iPhase)+dp
            else
                epseqCurr(iPhase) = 0.d0
            end if
        end do

! ----- Update elasticity yield
        sigyCurr = 0.d0
        if (phaseCurr(iPhaseHot) .gt. 0.d0) then
            if (l_hard_isotline) then
                sigyCurr = sigyCurr+(1-fmix)*h(iPhaseHot)*epseqCurr(iPhaseHot)
            end if
            if (l_hard_isotnlin) then
                sigyCurr = sigyCurr+(1-fmix)*(r(iPhaseHot)-sy(iPhaseHot))
            end if
        end if
!
        if (zSumCold .gt. 0.d0) then
            do iPhaseCold = 1, nbPhaseCold
                if (l_hard_isotline) then
                    sigyCurr = sigyCurr+ &
                               fmix*phaseCurr(iPhaseCold)*h(iPhaseCold)*epseqCurr(iPhaseCold)/ &
                               zSumCold
                end if
                if (l_hard_isotnlin) then
                    sigyCurr = sigyCurr+ &
                               fmix*phaseCurr(iPhaseCold)*(r(iPhaseCold)-sy(iPhaseCold))/zSumCold
                end if
            end do
        end if

        indicPlas = indicPlasCurr
    end if

! - Compute tangent matrix
    if (lMatr) then
        mode = 2
        if (l_visc) mode = 1
        dsidep(1:ndimsi, 1:ndimsi) = 0.d0
        do i = 1, ndimsi
            dsidep(i, i) = 1.d0
        end do
        do i = 1, 3
            do j = 1, 3
                dsidep(i, j) = dsidep(i, j)-1.d0/3.d0
            end do
        end do
        if (option(1:9) .eq. 'FULL_MECA') then
            coef1 = (1.5d0*deuxmu*trans+1.d0)
        else
            coef1 = 1.d0
        end if
        do i = 1, ndimsi
            do j = 1, ndimsi
                dsidep(i, j) = dsidep(i, j)*deuxmu/coef1
            end do
        end do

! ----- Plastic part
        b = 1.d0
        coef2 = 0.d0
        coef3 = 0.d0
        if (indicPlas .ge. 0.5d0) then
            if (option(1:9) .eq. 'FULL_MECA') then
                sigeps = 0.d0
                do i = 1, ndimsi
                    sigeps = sigeps+dvsigp(i)*dvdeps(i)
                end do
                if ((mode .eq. 1) .or. ((mode .eq. 2) .and. (sigeps .ge. 0.d0))) then
                    b = 1.d0-(1.5d0*deuxmu*dp/sieleq)
                    dv = 0.d0
                    if (mode .eq. 1) then
                        do iPhase = 1, nbPhase
                            n0(iPhase) = (1-n(iPhase))/n(iPhase)
                        end do
                        dv = (1-fmix)*phaseCurr(nbPhase)*(eta(nbPhase)/n(nbPhase)/timeIncr)* &
                             ((dp/timeIncr)**n0(nbPhase))
                        if (zSumCold .gt. 0.d0) then
                            do iPhaseCold = 1, nbPhaseCold
                                if (phaseCurr(iPhaseCold) .gt. 0.d0) then
                                    dv = dv+ &
                                         fmix*(phaseCurr(iPhaseCold)/zSumCold)* &
                                         (eta(iPhaseCold)/n(iPhaseCold)/timeIncr)* &
                                         ((dp/timeIncr)**n0(iPhaseCold))
                                end if
                            end do
                        end if
                    end if
                    coef2 = hmoy+dv
                    coef2 = (1.5d0*deuxmu*trans+1.d0)*coef2
                    coef2 = (1.5d0*deuxmu)+coef2
                    coef2 = 1/coef2-dp/sieleq
                    coef2 = ((1.5d0*deuxmu)**2)*coef2
                end if
            end if
            if (option(1:14) .eq. 'RIGI_MECA_TANG') then
                if (mode .eq. 2) coef2 = ((1.5d0*deuxmu)**2)/(1.5d0*deuxmu+hmoy)
            end if
            coef3 = coef2/coef1
        end if
        do i = 1, ndimsi
            do j = 1, ndimsi
                dsidep(i, j) = dsidep(i, j)*b
            end do
        end do
        do i = 1, 3
            do j = 1, 3
                dsidep(i, j) = dsidep(i, j)+troisk/3.d0
            end do
        end do
        do i = 1, ndimsi
            do j = 1, ndimsi
                dsidep(i, j) = dsidep(i, j)-coef3*sig0(i)*sig0(j)
            end do
        end do
    end if

! - Update internal state variables
    if (lVari) then
        vip(1:nbPhase) = epseqCurr(1:nbPhase)
        vip(viIndxIPlas) = indicPlasCurr
        vip(viIndxSigy) = sigyCurr
    end if
!
999 continue
!
end subroutine
