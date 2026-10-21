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
! aslint: disable=W1501
!
subroutine lcgdpm(option, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, carcri, &
                  timePrev, timeCurr, &
                  fm, df, &
                  nsig, sigm, &
                  nvi, vim, &
                  sigp, vip, dsigdf, &
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
#include "asterfort/zerop3.h"
!
    character(len=16), intent(in) :: option
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, ndim, jvMaterCode
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: timePrev, timeCurr
    real(kind=8), intent(in) :: fm(3, 3), df(3, 3)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    real(kind=8), intent(out) :: sigp(nsig)
    real(kind=8), intent(out) :: vip(nvi)
    real(kind=8), intent(out) :: dsigdf(6, 3, 3)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Comportment
!
! META_P_I* / META_V_I* for SIMO_MIEHE (-> META2P_I* / META2V_I*) strains and steel metallurgy
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  DF      : INCREMENT DU GRADIENT DE LA TRANSFORMATION
! IN  FM      : GRADIENT DE LA TRANSFORMATION A L INSTANT PRECEDENT
! IN  SIGM    : CONTRAINTES DE CAUCHY A L INSTANT PTECEDENT
! IN  VIM     : VARIABLES INTERNES A L INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! OUT SIGP    : CONTRAINTES DE CAUCHY A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIGDF  : DERIVEE DE SIGMA PAR RAPPORT A F
!
!          ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!          L'ORDRE :  XX,YY,ZZ,XY,XZ,YZ
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: kr(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(kind=8), parameter :: pdtsca(6) = (/1.d0, 1.d0, 1.d0, 2.d0, 2.d0, 2.d0/)
    real(kind=8), parameter :: hardCoef = 1.d0
    integer(kind=8) :: viIndxSigy, viIndxIPlas, viIndxSimo
    integer(kind=8) :: metaType
    integer(kind=8) :: nbPhase, nbPhaseCold, iPhaseCold, iPhaseHot, iPhase
    real(kind=8) :: phaseCurr(META_MECA_NBPHASE_MAXI), phasePrev(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: phaseIncr(META_MECA_NBPHASE_MAXI), zSumCold
    real(kind=8) :: dz1(META_MECA_NBPHASE_MAXI-1), dz2(META_MECA_NBPHASE_MAXI-1)
    real(kind=8) :: epseq(META_MECA_NBPHASE_MAXI), epseqPrev(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: epseqCurr(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: sy(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: h(META_MECA_NBPHASE_MAXI), hplus(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: r(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: annealTheta(2*(META_MECA_NBPHASE_MAXI-1))
    real(kind=8) :: eta(META_MECA_NBPHASE_MAXI), n(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: unsurn(META_MECA_NBPHASE_MAXI), c(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: m(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: kpt(META_MECA_NBPHASE_MAXI-1), fpt(META_MECA_NBPHASE_MAXI-1)
    real(kind=8) :: n0(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: sigyPrev, sigyCurr
    integer(kind=8) :: iValeMaxi, nbValeMaxi
    integer(kind=8) :: i, j, k, l, mode, iretTemp
    integer(kind=8) :: ind(3, 3), nbr
    real(kind=8) :: temp, timeIncr
    real(kind=8) :: epsth, e, mu, mum, troisk
    real(kind=8) :: fmix, hmoy, rmoy
    real(kind=8) :: cmoy, mmoy, cr
    real(kind=8) :: annealFunc, vimoy, ds
    real(kind=8) :: trans
    real(kind=8) :: jm, jp, dj, dfb(3, 3)
    real(kind=8) :: taum(6), dvtaum(6), trtaum, eqtaum
    real(kind=8) :: taup(6), dvtaup(6), trtaup
    real(kind=8) :: tau(6), dvtau(6), trtau, eqtau
    real(kind=8) :: bem(6), bel(6), dvbel(6), trbel
    real(kind=8) :: dvtel(6), eqtel
    real(kind=8) :: dp, seuil, mutild
    real(kind=8) :: je2, je3, xm, xp, sol(3)
    real(kind=8) :: coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7
    real(kind=8) :: coeff8, coeff9, dv, rb
    real(kind=8) :: mat0(3, 3), mat1(3, 3), mat2(6, 3, 3), mat3(3, 3)
    real(kind=8) :: precr
    real(kind=8) :: trElasPrev, trElasCurr
    character(len=1) :: poum
    integer(kind=8) :: test
    real(kind=8) :: indicPlasPrev, indicPlasCurr, indicPlas
    aster_logical :: lMatr, lVari, l_temp, lSigm
    aster_logical :: l_visc, l_anneal, lPlasTran, l_hard_isotline, l_hard_isotnlin
    character(len=16) :: metaRela, metaGlob
    data ind/1, 4, 5, &
        4, 2, 6, &
        5, 6, 3/
!  FB(I,J)=(J**(-1/3))*F(I,J)
!  GPM(I,J) : DEFORMATION PLASTIQUE LAGRANGIENNE A L'INSTANT PRECEDENT
!  B(I,J)=FB(I,K)*GPM(K,L)*FB(J,L)
!  BETR(1)=B(1,1), BETR(2)=B(2,2),...,BETR(6)=B(2,3)=B(3,2)
!
! TENSEUR DE DEFORMATION ELASTIQUE EULERIEN : BE
! TENSEUR DE DEFORMATION PLASTIQUE LAGRANGIEN : GP
! TAU : TENSEUR DE KIRCHHOFF
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    lMatr = L_MATR(option)
    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    timeIncr = timeCurr-timePrev
    precr = r8prem()

    if (lVari) then
        vip(1:nvi) = 0.d0
    end if
    if (lMatr) then
        dsigdf(1:6, 1:3, 1:3) = 0.d0
    end if
    if (lSigm) then
        sigp(1:2*ndim) = 0.d0
    end if

! - Behaviour in kit
    metaRela = compor(META_RELA)
    metaGlob = compor(META_GLOB)

! - Get metallurgy type
    call metaGetType(metaType, nbPhase)
    ASSERT(metaType .eq. META_STEEL .or. metaType .eq. META_STEEL_TEMPER)
    ASSERT(nbPhase .le. META_MECA_NBPHASE_MAXI)
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
                         e_=e, mu_=mu, troisk_=troisk, &
                         mum_=mum)

! - Get internal state variables
    viIndxSigy = nbPhase+1
    viIndxIPlas = nbPhase+2
    viIndxSimo = nbPhase+3
    epseqPrev(1:nbPhase) = vim(1:nbPhase)
    indicPlasPrev = vim(viIndxIPlas)
    sigyPrev = vim(viIndxSigy)
    trElasPrev = vim(viIndxSimo)

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

! - Determinant
    jm = fm(1, 1)*(fm(2, 2)*fm(3, 3)-fm(2, 3)*fm(3, 2))- &
         fm(2, 1)*(fm(1, 2)*fm(3, 3)-fm(1, 3)*fm(3, 2))+ &
         fm(3, 1)*(fm(1, 2)*fm(2, 3)-fm(1, 3)*fm(2, 2))
    dj = df(1, 1)*(df(2, 2)*df(3, 3)-df(2, 3)*df(3, 2))- &
         df(2, 1)*(df(1, 2)*df(3, 3)-df(1, 3)*df(3, 2))+ &
         df(3, 1)*(df(1, 2)*df(2, 3)-df(1, 3)*df(2, 2))
    jp = jm*dj
    do i = 1, 3
        do j = 1, 3
            if (dj .le. 0.d0) then
                codret = 1
                goto 999
            else
                dfb(i, j) = (dj**(-1.d0/3.d0))*df(i, j)
            end if
        end do
    end do

! - CONTRAINTES DE KIRCHHOFF A L INSTANT PRECEDENT
    taum(5) = 0.d0
    taum(6) = 0.d0
    do i = 1, 2*ndim
        taum(i) = jm*sigm(i)
    end do
!
    trtaum = taum(1)+taum(2)+taum(3)
    eqtaum = 0.d0
    do i = 1, 6
        dvtaum(i) = taum(i)-kr(i)*trtaum/3.d0
        eqtaum = eqtaum+pdtsca(i)*(dvtaum(i)**2.d0)
    end do
    eqtaum = sqrt(1.5d0*eqtaum)

! 3.3 - DEFORMATIONS ELASTIQUES A L INSTANT PRECEDENT :
! BEM=DVTAUM/MUM+KR*TRBEM/3.D0
    xm = (jm**(-2.d0/3.d0))*(1.d0-2.d0*trElasPrev/3.d0)
    do i = 1, 6
        bem(i) = dvtaum(i)/mum+kr(i)*xm
    end do

! 3.4 - BEL(I,J)=DFB(I,K)*BEM(K,L)*DFB(J,L)
    do i = 1, 3
        do j = 1, 3
            bel(ind(i, j)) = 0.d0
            do k = 1, 3
                do l = 1, 3
                    bel(ind(i, j)) = bel(ind(i, j))+dfb(i, k)*bem(ind(k, l))*dfb(j, l)
                end do
            end do
        end do
    end do

! 3.5 - TRACE ET PARTIE DEVIATORIQUE DU TENSEUR BEL
    trbel = bel(1)+bel(2)+bel(3)
    do i = 1, 6
        dvbel(i) = bel(i)-kr(i)*trbel/3.d0
    end do

! 3.6 - CONTRAINTE ELASTIQUE (PARTIE DEVIATORIQUE)
    do i = 1, 6
        dvtel(i) = mu*dvbel(i)
    end do

! 3.7 - CONTRAINTE EQUIVALENTE ELASTIQUE ET SEUIL
    eqtel = 0.d0
    do i = 1, 6
        eqtel = eqtel+pdtsca(i)*dvtel(i)*dvtel(i)
    end do
    eqtel = sqrt(1.5d0*eqtel)

! 3.8 - TRACE DU TENSEUR DE KIRCHHOFF (CONNUE CAR NE DEPEND QUE DE J)
    trtaup = (troisk*((jp*jp)-1.d0)/6.d0)-(troisk*epsth*(jp+(1.d0/jp))/2.d0)

! - Integration of behaviour
    indicPlasCurr = indicPlasPrev
    dp = 0.d0
    if (lVari) then
        indicPlasCurr = 0.d0
        seuil = eqtel-(1.d0+mu*trans*trbel)*rmoy

! ----- Compute DP (plastic multiplier)
        if (seuil .lt. 0.d0) then
            indicPlasCurr = 0.d0
            dp = 0.d0
        else
            indicPlasCurr = 1.d0
            mutild = 2.d0*mu*trbel/3.d0
            call nzcalc(carcri, nbPhase, phaseCurr, zSumCold, &
                        fmix, seuil, timeIncr, trans, &
                        hmoy, mutild, eta, unsurn, &
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
                        do iPhase = 1, nbPhase-1
                            if (phaseCurr(iPhase) .gt. 0.d0) then
                                rmoy = rmoy+phaseCurr(iPhase)*(r(iPhase)-h(iPhase)*dp)
                                hmoy = hmoy+phaseCurr(iPhase)*h(iPhase)
                            end if
                        end do
                        rmoy = fmix*rmoy/zSumCold
                        hmoy = fmix*hmoy/zSumCold
                    end if
                    if (phaseCurr(iPhaseHot) .gt. 0.d0) then
                        rmoy = (1.d0-fmix)*(r(iPhaseHot)-h(iPhaseHot)*dp)+rmoy
                        hmoy = (1.d0-fmix)*h(iPhaseHot)+hmoy
                    end if
                    seuil = eqtel-(1.d0+mu*trans*trbel)*rmoy
                    call nzcalc(carcri, nbPhase, phaseCurr, zSumCold, &
                                fmix, seuil, timeIncr, trans, &
                                hmoy, mutild, eta, unsurn, &
                                dp, codret)
                    if (codret .eq. 1) goto 999
                end do
            end if
        end if

! ----- Compute stresses
        do i = 1, 6
            if (eqtel .gt. 0.d0) then
                dvtaup(i) = dvtel(i)-mu*dp*trbel*dvtel(i)/eqtel
                dvtaup(i) = dvtaup(i)/(1.d0+mu*trans*trbel)
            else
                dvtaup(i) = dvtel(i)/(1.d0+mu*trans*trbel)
            end if
            taup(i) = dvtaup(i)+kr(i)*trtaup
        end do
        do i = 1, 2*ndim
            sigp(i) = taup(i)/jp
        end do

! ----- Update cumulated plastic strain
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
        if (option(1:14) .eq. 'RIGI_MECA_TANG') then
            do i = 1, 6
                tau(i) = taum(i)
            end do
        else
            do i = 1, 6
                tau(i) = taup(i)
            end do
        end if
        mat0(1, 1) = df(2, 2)*df(3, 3)-df(2, 3)*df(3, 2)
        mat0(2, 2) = df(1, 1)*df(3, 3)-df(1, 3)*df(3, 1)
        mat0(3, 3) = df(1, 1)*df(2, 2)-df(1, 2)*df(2, 1)
        mat0(1, 2) = df(3, 1)*df(2, 3)-df(2, 1)*df(3, 3)
        mat0(2, 1) = df(1, 3)*df(3, 2)-df(1, 2)*df(3, 3)
        mat0(1, 3) = df(2, 1)*df(3, 2)-df(3, 1)*df(2, 2)
        mat0(3, 1) = df(1, 2)*df(2, 3)-df(1, 3)*df(2, 2)
        mat0(2, 3) = df(3, 1)*df(1, 2)-df(1, 1)*df(3, 2)
        mat0(3, 2) = df(2, 1)*df(1, 3)-df(1, 1)*df(2, 3)
        do i = 1, 3
            do j = 1, 3
                mat1(i, j) = 0.d0
                do k = 1, 3
                    mat1(i, j) = mat1(i, j)+dfb(i, k)*bem(ind(k, j))
                end do
            end do
        end do
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    do l = 1, 3
                        mat2(ind(i, j), k, l) = kr(ind(i, k))*mat1(j, l)+kr(ind(j, k))*mat1(i, l)
                    end do
                end do
            end do
        end do
        eqtau = 0.d0
        trtau = (tau(1)+tau(2)+tau(3))/3.d0
        do i = 1, 6
            dvtau(i) = tau(i)-trtau*kr(i)
            eqtau = eqtau+pdtsca(i)*(dvtau(i)**2.d0)
        end do
        eqtau = sqrt(1.5d0*eqtau)
        if (eqtau .gt. 0.d0) then
            coeff1 = 1.d0+mu*trans*trbel+mu*dp*trbel/eqtau
        else
            coeff1 = 1.d0+mu*trans*trbel
        end if
        coeff2 = (dj**(-1.d0/3.d0))*mu/(coeff1*jp)
        coeff3 = -2.d0*mu/(3.d0*coeff1*jp*dj)
        coeff4 = (troisk*jp/3.d0)-troisk*epsth*(1.d0-(jp**(-2.d0)))/2.d0
        coeff4 = coeff4/dj
        do i = 1, 6
            do j = 1, 3
                do k = 1, 3
                    dsigdf(i, j, k) = mat2(i, j, k)-2.d0*kr(i)*mat1(j, k)/3.d0
                    dsigdf(i, j, k) = coeff2*dsigdf(i, j, k)+coeff3* &
                                      dvbel(i)*mat0(j, k)+coeff4*kr(i)*mat0(j, k) &
                                      -tau(i)*mat0(j, k)/(jp*dj)
                end do
            end do
        end do
!
        if (abs(indicPlas) .lt. 0.5d0) then
            coeff5 = -2.d0*trans*coeff2
            coeff6 = 2.d0*trans*mu*trbel/(3.d0*jp*dj*coeff1)
            do i = 1, 6
                do j = 1, 3
                    do k = 1, 3
                        dsigdf(i, j, k) = dsigdf(i, j, k)+coeff5*dvtau(i)*mat1(j, k)+ &
                                          coeff6*dvtau(i)*mat0(j, k)
                    end do
                end do
            end do
        else
            mode = 2
            if (l_visc) mode = 1
            if (mode .eq. 1) then
                if (dp .gt. 0.d0) then
                    do iPhase = 1, nbPhase
                        n0(iPhase) = (1-n(iPhase))/n(iPhase)
                    end do
                    dv = (1-fmix)*phaseCurr(iPhaseHot)*(eta(iPhaseHot)/n(iPhaseHot)/timeIncr)* &
                         ((dp/timeIncr)**n0(iPhaseHot))
                    if (zSumCold .gt. 0.d0) then
                        do iPhaseCold = 1, nbPhaseCold
                            if (phaseCurr(iPhaseCold) .gt. 0.d0) then
                                dv = dv+fmix* &
                                     (phaseCurr(iPhaseCold)/zSumCold)* &
                                     (eta(iPhaseCold)/n(iPhaseCold)/timeIncr)* &
                                     ((dp/timeIncr)**n0(iPhaseCold))
                            end if
                        end do
                    end if
                else
                    dv = 0.d0
                end if
                if (option(1:9) .eq. 'FULL_MECA') then
                    rb = hmoy+dv
                else
                    rb = 0.d0
                end if
            else
                rb = hmoy
            end if
            if ((option(1:9) .eq. 'FULL_MECA') .or. &
                ((option(1:14) .eq. 'RIGI_MECA_TANG') .and. (mode .eq. 2))) then
                coeff5 = mu*trbel+rb*(1.d0+mu*trans*trbel)
                coeff6 = -3.d0*mu*trbel*(eqtau-rb*dp)/((eqtau**3.d0)*coeff5)
                coeff6 = coeff6*coeff2
                coeff7 = -2.d0*coeff1*rb*(eqtau*trans+dp)/(eqtau*coeff5)
                coeff7 = coeff7*coeff2
                coeff8 = 0.d0
                do i = 1, 6
                    coeff8 = coeff8+pdtsca(i)*dvtau(i)*dvbel(i)
                end do
                coeff9 = coeff1*rb*(eqtau*trans+dp)/eqtau+ &
                         3.d0*mu*coeff8*(eqtau-rb*dp)/(2.d0*(eqtau**3.d0))
                coeff9 = -coeff9*coeff3*trbel/coeff5
                do i = 1, 3
                    do j = 1, 3
                        mat3(i, j) = 0.d0
                        do k = 1, 3
                            mat3(i, j) = mat3(i, j)+dvtau(ind(i, k))*mat1(k, j)
                        end do
                    end do
                end do
                do i = 1, 6
                    do j = 1, 3
                        do k = 1, 3
                            dsigdf(i, j, k) = dsigdf(i, j, k)+ &
                                              coeff6*dvtau(i)*mat3(j, k)+ &
                                              coeff7*dvtau(i)*mat1(j, k)+ &
                                              coeff9*dvtau(i)*mat0(j, k)
                        end do
                    end do
                end do
            end if
        end if
    end if

! - CORRECTION SUR TRBE
    if (lVari) then
        trtau = taup(1)+taup(2)+taup(3)
        eqtau = 0.d0
        do i = 1, 6
            dvtau(i) = taup(i)-kr(i)*trtau/3.d0
            eqtau = eqtau+pdtsca(i)*(dvtau(i)**2.d0)
        end do
        eqtau = sqrt(1.5d0*eqtau)
        je2 = (eqtau**2.d0)/(3.d0*(mu**2.d0))
        je3 = dvtau(1)*(dvtau(2)*dvtau(3)-dvtau(6)*dvtau(6))- &
              dvtau(4)*(dvtau(4)*dvtau(3)-dvtau(5)*dvtau(6))+ &
              dvtau(5)*(dvtau(4)*dvtau(6)-dvtau(5)*dvtau(2))
        je3 = je3/(mu**3.d0)
        call zerop3(0.d0, -je2, je3-1.d0, sol, nbr)
        if (nbr .le. 1) then
            xp = sol(1)
        else
            xp = sol(1)
            do i = 2, nbr
                if ((abs(sol(i)-xm)) .lt. (abs(sol(i-1)-xm))) xp = sol(i)
            end do
        end if
        trElasCurr = 3.d0*(1.d0-(jp**(2.d0/3.d0))*xp)/2.d0
    end if

! - Update internal state variables
    if (lVari) then
        vip(1:nbPhase) = epseqCurr(1:nbPhase)
        vip(viIndxIPlas) = indicPlasCurr
        vip(viIndxSigy) = sigyCurr
        vip(viIndxSimo) = trElasCurr
    end if
!
999 continue
!
end subroutine
