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
subroutine nzcifw(option, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, carcri, &
                  timePrev, timeCurr, &
                  neps, epsm, deps, &
                  nsig, sigm, &
                  nvi, vim, &
                  sigp, vip, ndsde, dsidep, &
                  codret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/metaGetMechanism.h"
#include "asterfort/metaGetParaAnneal.h"
#include "asterfort/metaGetParaElas.h"
#include "asterfort/metaGetParaHardLine.h"
#include "asterfort/metaGetParaMixture.h"
#include "asterfort/metaGetParaPlasTransf.h"
#include "asterfort/metaGetParaVisc.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/metaGetType.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/nzcalc.h"
#include "asterfort/verift.h"
!
    character(len=16), intent(in) :: option
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, ndim, jvMaterCode
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: timePrev, timeCurr
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps)
    real(kind=8), intent(in) :: deps(neps)
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
! META_P_C* / META_V_C* for small strains and steel metallurgy
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  IMAT    : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT : RELCOM ET DEFORM
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
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
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8), parameter :: hardCoef = 2.d0/3.d0
    integer(kind=8) :: metaType
    integer(kind=8) :: nbPhase, nbPhaseCold, iPhaseCold, iPhaseHot, iPhase
    real(kind=8) :: phaseCurr(PSTEEL_NB), phasePrev(PSTEEL_NB), phaseIncr(PSTEEL_NB), zalpha
    real(kind=8) :: xcinColdPrev(6*PCSTEEL_NB), xcinColdCurr(6*PCSTEEL_NB)
    real(kind=8) :: xcinCold(6*PCSTEEL_NB)
    real(kind=8) :: xcinColdSave(6*PCSTEEL_NB)
    real(kind=8) :: xcinHotPrev(6), xcinHot(6), xcinHotSave(6), xcinMoy(6), xcinHotCurr(6)
    integer(kind=8) :: ndimsi, i, j, mode
    real(kind=8) :: timeIncr
    real(kind=8) :: epsth, e, deuxmu, deumum, troisk
    real(kind=8) :: fmix, sy(5), symoy, h(5), hmoy, rprim
    real(kind=8) :: annealTheta(8)
    real(kind=8) :: eta(5), n(5), unsurn(5), c(5), m(5), cmoy, mmoy, cr
    real(kind=8) :: dz1(4), dz2(4), annealFunc
    real(kind=8) :: xmoy(6), ds(6), xmoyeq
    real(kind=8) :: trans, kpt(4), fpt(4)
    real(kind=8) :: trepsm, trdeps, trsigm, trsigp
    real(kind=8) :: dvdeps(6), dvsigm(6), dvsigp(6)
    real(kind=8) :: sigel(6), sigel2(6), sig0(6), sieleq, sigeps
    real(kind=8) :: dp, seuil
    real(kind=8) :: coef1, coef2, coef3, dv, n0(5), b
    character(len=1) :: poum
    real(kind=8) :: indicPlasPrev, indicPlasCurr, indicPlas
    aster_logical :: lMatr, lVari
    aster_logical :: l_visc, l_plas, l_anneal, l_plas_tran, l_hard_isotline, l_hard_isotnlin
    character(len=16) :: metaRela, metaGlob
!
! --------------------------------------------------------------------------------------------------
!
    ndimsi = 2*ndim
    codret = 0
    lMatr = L_MATR(option)
    lVari = L_VARI(option)
    timeIncr = timeCurr-timePrev

! - Behaviour in kit
    metaRela = compor(META_RELA)
    metaGlob = compor(META_GLOB)

! - Get metallurgy type
    call metaGetType(metaType, nbPhase)
    ASSERT(metaType .eq. META_STEEL)
    ASSERT(nbPhase .eq. PSTEEL_NB)
    nbPhaseCold = nbPhase-1
    iPhaseHot = nbPhase
    ASSERT(nbPhaseCold .eq. PCSTEEL_NB)

! - Get phases
    if (lVari) then
        poum = '+'
        call metaGetPhase(fami, '+', kpg, ksp, metaType, &
                          nbPhase, phaseCurr, zcold_=zalpha)
        call metaGetPhase(fami, '-', kpg, ksp, metaType, &
                          nbPhase, phasePrev)
    else
        poum = '-'
        call metaGetPhase(fami, '-', kpg, ksp, metaType, &
                          nbPhase, phaseCurr, zcold_=zalpha)
    end if
    phaseIncr(1:nbPhaseCold) = phaseCurr(1:nbPhaseCold)-phasePrev(1:nbPhaseCold)

! - Compute thermic strain
    call verift(fami, kpg, ksp, poum, jvMaterCode, &
                epsth_meta_=epsth)

! - Get active mechanisms in behaviour
    call metaGetMechanism(metaRela, metaGlob, &
                          l_plas=l_plas, &
                          l_visc=l_visc, &
                          l_anneal=l_anneal, &
                          l_plas_tran=l_plas_tran, &
                          l_hard_isotline=l_hard_isotline, &
                          l_hard_isotnlin=l_hard_isotnlin)

! - Get elastic parameters
    call metaGetParaElas(poum, fami, kpg, ksp, jvMaterCode, &
                         e_=e, deuxmu_=deuxmu, troisk_=troisk, &
                         deuxmum_=deumum)

! - Get internal state variables
    xcinColdPrev(1:6*nbPhaseCold) = vim(1:6*nbPhaseCold)
    do i = 1, ndimsi
        xcinHotPrev(i) = vim(6*nbPhaseCold+i)
    end do
    indicPlasPrev = vim(IDX_C_IPLAS)

! - Mixture law (yield limit)
    call metaGetParaMixture(poum, fami, kpg, ksp, jvMaterCode, &
                            l_visc, metaType, nbPhase, zalpha, fmix, &
                            sy)

! - Get hardening slope (linear)
    call metaGetParaHardLine(poum, fami, kpg, ksp, jvMaterCode, &
                             metaType, nbPhase, &
                             e, hardCoef, h)
    hmoy = 0.d0
    do iPhase = 1, nbPhase
        hmoy = hmoy+phaseCurr(iPhase)*h(iPhase)
    end do

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
                                 metaType, nbPhase, eta, n, unsurn, &
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
            do i = 1, ndimsi
                annealFunc = 0.d0
                do iPhaseCold = 1, nbPhaseCold
                    annealFunc = annealFunc+ &
                                 dz2(iPhaseCold)* &
                                 (annealTheta(4+iPhaseCold)*xcinColdPrev(6*(iPhaseCold-1)+i)- &
                                  xcinHotPrev(i))/phaseCurr(iPhaseHot)
                end do
                xcinHot(i) = xcinHotPrev(i)+annealFunc
                if ((xcinHot(i)*xcinHotPrev(i)) .lt. 0.d0) then
                    xcinHot(i) = 0.d0
                end if
            end do
        else
            xcinHot = 0.d0
        end if

! ----- Metallurgical annealing for cold phases
        do iPhaseCold = 1, nbPhaseCold
            do i = 1, ndimsi
                if (phaseCurr(iPhaseCold) .gt. 0.d0) then
                    annealFunc = dz1(iPhaseCold)* &
                                 (annealTheta(iPhaseCold)*xcinHotPrev(i)- &
                                  xcinColdPrev(6*(iPhaseCold-1)+i))/phaseCurr(iPhaseCold)
                    xcinCold(6*(iPhaseCold-1)+i) = xcinColdPrev(6*(iPhaseCold-1)+i)+annealFunc
                    if ((xcinCold(6*(iPhaseCold-1)+i)* &
                         xcinColdPrev(6*(iPhaseCold-1)+i)) .lt. 0.d0) then
                        xcinCold(6*(iPhaseCold-1)+i) = 0.d0
                    end if
                else
                    xcinCold(6*(iPhaseCold-1)+i) = 0.d0
                end if
            end do
        end do
        do i = 4, ndimsi
            do iPhaseCold = 1, nbPhaseCold
                xcinCold(6*(iPhaseCold-1)+i) = xcinCold(6*(iPhaseCold-1)+i)*rac2
            end do
            xcinHot(i) = xcinHot(i)*rac2
        end do

! ----- Mean hardening tensor
        do i = 1, ndimsi
            xmoy(i) = 0.d0
            do iPhaseCold = 1, nbPhaseCold
                xmoy(i) = xmoy(i)+ &
                          phaseCurr(iPhaseCold)*h(iPhaseCold)*xcinCold(6*(iPhaseCold-1)+i)
            end do
            xmoy(i) = xmoy(i)+phaseCurr(iPhaseHot)*h(iPhaseHot)*xcinHot(i)
        end do
        xmoyeq = 0.d0
        do i = 1, ndimsi
            xmoyeq = xmoyeq+xmoy(i)**2.d0
        end do
        xmoyeq = sqrt(1.5d0*xmoyeq)

! ----- Viscuous annealing
        cmoy = 0.d0
        mmoy = 0.d0
        do iPhase = 1, nbPhase
            cmoy = cmoy+phaseCurr(iPhase)*c(iPhase)
            mmoy = mmoy+phaseCurr(iPhase)*m(iPhase)
        end do
        cr = cmoy*xmoyeq
        if (xmoyeq .gt. 0.d0) then
            do i = 1, ndimsi
                ds(i) = 3.d0*timeIncr*(cr**mmoy)*xmoy(i)/(2.d0*xmoyeq)
            end do
        else
            do i = 1, ndimsi
                ds(i) = 0.d0
            end do
        end if

        do i = 1, ndimsi
            do iPhaseCold = 1, nbPhaseCold
                if (phaseCurr(iPhaseCold) .gt. 0.d0) then
                    xcinColdSave(6*(iPhaseCold-1)+i) = xcinCold(6*(iPhaseCold-1)+i)
                    xcinCold(6*(iPhaseCold-1)+i) = xcinCold(6*(iPhaseCold-1)+i)-ds(i)
                    if ((xcinCold(6*(iPhaseCold-1)+i)* &
                         xcinColdSave(6*(iPhaseCold-1)+i)) .lt. 0.d0) then
                        xcinCold(6*(iPhaseCold-1)+i) = 0.d0
                    end if
                end if
            end do
            if (phaseCurr(iPhaseHot) .gt. 0.d0) then
                xcinHotSave(i) = xcinHot(i)
                xcinHot(i) = xcinHot(i)-ds(i)
                if ((xcinHot(i)*xcinHotSave(i)) .lt. 0.d0) then
                    xcinHot(i) = 0.d0
                end if
            end if
        end do

! ----- Parameters for plasticity of tranformation
        trans = 0.d0
        if (l_plas_tran) then
            call metaGetParaPlasTransf('+', fami, 1, 1, jvMaterCode, &
                                       metaType, nbPhase, phaseIncr, zalpha, &
                                       kpt, fpt)
            do iPhaseCold = 1, nbPhaseCold
                if (phaseIncr(iPhaseCold) .gt. 0.d0) then
                    trans = trans+kpt(iPhaseCold)*fpt(iPhaseCold)*phaseIncr(iPhaseCold)
                end if
            end do
        end if

    else
        do i = 1, ndimsi
            do iPhaseCold = 1, nbPhaseCold
                xcinCold(6*(iPhaseCold-1)+i) = xcinColdPrev(6*(iPhaseCold-1)+i)
                if (i .ge. 4) then
                    xcinCold(6*(iPhaseCold-1)+i) = xcinCold(6*(iPhaseCold-1)+i)*rac2
                end if
            end do
            xcinHot(i) = xcinHotPrev(i)
            if (i .ge. 4) then
                xcinHot(i) = xcinHot(i)*rac2
            end if
        end do
        trans = 0.d0
        do i = 1, ndimsi
            xmoy(i) = 0.d0
            do iPhaseCold = 1, nbPhaseCold
                xmoy(i) = xmoy(i)+ &
                          phaseCurr(iPhaseCold)*h(iPhaseCold)*xcinCold(6*(iPhaseCold-1)+i)
            end do
            xcinHot(i) = xmoy(i)+ &
                         phaseCurr(iPhaseHot)*h(iPhaseHot)*xcinHotPrev(i)
        end do
    end if

! - Apply mixture law on hardening parameters
    if (zalpha .gt. 0.d0) then
        symoy = phaseCurr(1)*sy(1)+phaseCurr(2)*sy(2)+phaseCurr(3)*sy(3)+phaseCurr(4)*sy(4)
        symoy = symoy/zalpha
    else
        symoy = 0.d0
    end if
    symoy = (1.d0-fmix)*sy(nbPhase)+fmix*symoy

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
    sigel2 = 0.d0
    do i = 1, ndimsi
        sigel(i) = deuxmu*dvsigm(i)/deumum+deuxmu*dvdeps(i)
        sigel2(i) = sigel(i)-(1.5d0*deuxmu*trans+1.d0)*xmoy(i)
        sieleq = sieleq+sigel2(i)**2
    end do
    sieleq = sqrt(1.5d0*sieleq)
    if (sieleq .gt. 0.d0) then
        do i = 1, ndimsi
            sig0(i) = sigel2(i)/sieleq
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
        seuil = sieleq-(1.5d0*deuxmu*trans+1.d0)*symoy

! ----- Compute DP (plastic multiplier)
        if (seuil .lt. 0.d0) then
            indicPlasCurr = 0.d0
            dp = 0.d0
        else
            indicPlasCurr = 1.d0
            rprim = 3.d0*hmoy/2.d0
            if (l_plas) then
                dp = seuil/(1.5d0*deuxmu+(1.5d0*deuxmu*trans+1.d0)*rprim)
            else
                call nzcalc(carcri, nbPhase, phaseCurr, zalpha, &
                            fmix, seuil, timeIncr, trans, &
                            rprim, deuxmu, eta, unsurn, &
                            dp, codret)
                if (codret .eq. 1) goto 999
            end if
        end if

! ----- Compute stresses
        sigp(1:2*ndim) = 0.d0
        do i = 1, ndimsi
            dvsigp(i) = sigel(i)-1.5d0*deuxmu*dp*sig0(i)
            dvsigp(i) = dvsigp(i)/(1.5d0*deuxmu*trans+1.d0)
            sigp(i) = dvsigp(i)+trsigp*kron(i)
        end do

! ----- Update internal state variables
        do i = 1, ndimsi
            do iPhaseCold = 1, nbPhaseCold
                if (phaseCurr(iPhaseCold) .gt. 0.d0) then
                    xcinColdCurr(6*(iPhaseCold-1)+i) = xcinCold(6*(iPhaseCold-1)+i)+ &
                                                       3.d0*dp*sig0(i)/2.d0
                    if (i .ge. 4) then
                        xcinColdCurr(6*(iPhaseCold-1)+i) = xcinColdCurr(6*(iPhaseCold-1)+i)/rac2
                    end if
                else
                    xcinColdCurr(6*(iPhaseCold-1)+i) = 0.d0
                end if
            end do
            if (phaseCurr(iPhaseHot) .gt. 0.d0) then
                xcinHotCurr(i) = xcinHot(i)+3.d0*dp*sig0(i)/2.d0
                if (i .ge. 4) then
                    xcinHotCurr(i) = xcinHotCurr(i)/rac2
                end if
            else
                xcinHotCurr(i) = 0.d0
            end if
        end do

! ----- Mean value of kinematic hardening tensor
        xcinMoy = 0.d0
        do i = 1, ndimsi
            xcinMoy(i) = xmoy(i)+3.d0*hmoy*dp*sig0(i)/2.d0
            if (i .ge. 4) then
                xcinMoy(i) = xcinMoy(i)/rac2
            end if
        end do

        indicPlas = indicPlasCurr
    end if

! - Compute tangent matrix
    if (lMatr) then
        mode = 2
        if (l_visc) mode = 1
        dsidep(1:6, 1:6) = 0.d0
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
                    dvsigp(i) = dvsigp(i)-xmoy(i)
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
                        if (zalpha .gt. 0.d0) then
                            do iPhaseCold = 1, nbPhaseCold
                                if (phaseCurr(iPhaseCold) .gt. 0.d0) then
                                    dv = dv+ &
                                         fmix*(phaseCurr(iPhaseCold)/zalpha)* &
                                         (eta(iPhaseCold)/n(iPhaseCold)/timeIncr)* &
                                         ((dp/timeIncr)**n0(iPhaseCold))
                                end if
                            end do
                        end if
                    end if
                    coef2 = 1.5d0*hmoy+dv
                    coef2 = (1.5d0*deuxmu*trans+1.d0)*coef2
                    coef2 = (1.5d0*deuxmu)+coef2
                    coef2 = 1/coef2-dp/sieleq
                    coef2 = ((1.5d0*deuxmu)**2)*coef2
                end if
            end if
            if (option(1:14) .eq. 'RIGI_MECA_TANG') then
                if (mode .eq. 2) coef2 = ((1.5d0*deuxmu)**2)/(1.5d0*deuxmu+1.5d0*hmoy)
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
        vip(1:6*nbPhaseCold) = xcinColdCurr(1:6*nbPhaseCold)
        do i = 1, ndimsi
            vip(6*nbPhaseCold+i) = xcinHotCurr(i)
            vip(6*nbPhaseCold+6+i) = xcinMoy(i)
        end do
        vip(IDX_C_IPLAS) = indicPlas
    end if
!
999 continue
!
end subroutine
