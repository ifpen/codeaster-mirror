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
subroutine lcgdpm(fami, kpg, ksp, ndim, imat, &
                  compor, carcri, instam, instap, fm, &
                  df, sigm, vim, option, sigp, &
                  vip, dsigdf, iret)
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/metaGetMechanism.h"
#include "asterfort/metaGetParaAnneal.h"
#include "asterfort/metaGetParaElas.h"
#include "asterfort/metaGetParaHardLine.h"
#include "asterfort/metaGetParaHardTrac.h"
#include "asterfort/metaGetParaMixture.h"
#include "asterfort/metaGetParaPlasTransf.h"
#include "asterfort/metaGetParaVisc.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/metaGetType.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/nzcalc.h"
#include "asterfort/rcvarc.h"
#include "asterfort/verift.h"
#include "asterfort/zerop3.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: imat
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam
    real(kind=8), intent(in) :: instap
    real(kind=8), intent(in) :: fm(3, 3)
    real(kind=8), intent(in) :: df(3, 3)
    real(kind=8), intent(in) :: sigm(*)
    real(kind=8), intent(in) :: vim(*)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(*)
    real(kind=8), intent(out) :: vip(*)
    real(kind=8), intent(out) :: dsigdf(6, 3, 3)
    integer(kind=8), intent(out) :: iret
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
    integer(kind=8), parameter :: vixIndiPlas = 7, vixSigy = 6, vixSimo = 8
    integer(kind=8) :: maxval, nbPhase, metaType
    integer(kind=8) :: i, j, k, l, mode, iret2
    integer(kind=8) :: ind(3, 3), nbr
    real(kind=8) :: phaseCurr(5), phasePrev(5), zalpha, deltaz(4)
    real(kind=8) :: temp, dt, coef_hard
    real(kind=8) :: epsth, e, mu, mum, troisk
    real(kind=8) :: fmel, sy(5), h(5), hmoy, hplus(5), r(5), rmoy
    real(kind=8) :: theta(8)
    real(kind=8) :: eta(5), n(5), unsurn(5), c(5), m(5), cmoy, mmoy, cr
    real(kind=8) :: dz(4), dz1(4), dz2(4), vi(5), dvin, vimoy, ds
    real(kind=8) :: trans, kpt(4), fpt(4)
    real(kind=8) :: jm, jp, dj, dfb(3, 3)
    real(kind=8) :: taum(6), dvtaum(6), trtaum, eqtaum
    real(kind=8) :: taup(6), dvtaup(6), trtaup
    real(kind=8) :: tau(6), dvtau(6), trtau, eqtau
    real(kind=8) :: bem(6), bel(6), dvbel(6), trbel
    real(kind=8) :: dvtel(6), eqtel
    real(kind=8) :: plasti, dp, seuil, mutild
    real(kind=8) :: je2, je3, xm, xp, sol(3)
    real(kind=8) :: coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7
    real(kind=8) :: coeff8, coeff9, dv, rb, n0(5)
    real(kind=8) :: mat0(3, 3), mat1(3, 3), mat2(6, 3, 3), mat3(3, 3)
    character(len=1) :: poum
    integer(kind=8) :: test
    real(kind=8), parameter :: kr(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
    real(kind=8), parameter :: pdtsca(6) = (/1.d0, 1.d0, 1.d0, 2.d0, 2.d0, 2.d0/)
    aster_logical :: resi, rigi, l_temp
    aster_logical :: l_visc, l_anneal, l_plas_tran, l_hard_isotline, l_hard_isotnlin
    aster_logical :: lMatr, lVari, lSigm
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
! TAU : TENSEUR DE KIRSHHOFF
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
    resi = option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL'
    rigi = option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL'
    dt = instap-instam
    lVari = L_VARI(option)
    lSigm = L_SIGM(option)
    lMatr = L_MATR(option)
    if (lVari) then
        vip(1:8) = 0.d0
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
    ASSERT(metaType .eq. META_STEEL)
    ASSERT(nbPhase .eq. 5)

! - Get phases
    if (resi) then
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
    do k = 1, nbPhase-1
        deltaz(k) = phaseCurr(k)-phasePrev(k)
    end do
!
! - Compute thermic strain
!
    call verift(fami, kpg, ksp, poum, imat, &
                epsth_meta_=epsth)
    call rcvarc(' ', 'TEMP', poum, fami, kpg, &
                ksp, temp, iret2)
    l_temp = iret2 .eq. 0

! - Mechanisms of comportment law
    call metaGetMechanism(metaRela, metaGlob, &
                          l_visc=l_visc, &
                          l_anneal=l_anneal, &
                          l_plas_tran=l_plas_tran, &
                          l_hard_isotline=l_hard_isotline, &
                          l_hard_isotnlin=l_hard_isotnlin)
!
! - Get elastic parameters
!
    call metaGetParaElas(poum, fami, kpg, ksp, imat, &
                         e_=e, mu_=mu, troisk_=troisk, &
                         mum_=mum)
    plasti = vim(vixIndiPlas)

! - Mixture law (yield limit)
    call metaGetParaMixture(poum, fami, kpg, ksp, imat, &
                            l_visc, metaType, nbPhase, zalpha, fmel, &
                            sy)

    if (resi) then
! ----- Parameters for annealing
        if (l_anneal) then
            call metaGetParaAnneal(poum, fami, kpg, ksp, imat, &
                                   metaType, nbPhase, &
                                   theta)
        else
            do i = 1, 8
                theta(i) = 1.d0
            end do
        end if
! ----- Parameters for viscosity
        if (l_visc) then
            call metaGetParaVisc(poum, fami, kpg, ksp, imat, &
                                 metaType, nbPhase, eta, n, unsurn, &
                                 c, m)
        else
            eta(:) = 0.d0
            n(:) = 20.d0
            unsurn(:) = 1.d0
            c(:) = 0.d0
            m(:) = 20.d0
        end if
!
! 2.6 - CALCUL DE VIM+DG-DS
!
        do k = 1, nbPhase-1
            dz(k) = phaseCurr(k)-phasePrev(k)
            if (dz(k) .ge. 0.d0) then
                dz1(k) = dz(k)
                dz2(k) = 0.d0
            else
                dz1(k) = 0.d0
                dz2(k) = -dz(k)
            end if
        end do
        if (phaseCurr(nbPhase) .gt. 0.d0) then
            dvin = 0.d0
            do k = 1, nbPhase-1
                dvin = dvin+dz2(k)*(theta(4+k)*vim(k)-vim(nbPhase))/ &
                       phaseCurr(nbPhase)
            end do
            vi(nbPhase) = vim(nbPhase)+dvin
            vimoy = phaseCurr(nbPhase)*vi(nbPhase)
        else
            vi(nbPhase) = 0.d0
            vimoy = 0.d0
        end if
        do k = 1, nbPhase-1
            if (phaseCurr(k) .gt. 0.d0) then
                dvin = dz1(k)*(theta(k)*vim(nbPhase)-vim(k))/phaseCurr(k)
                vi(k) = vim(k)+dvin
                vimoy = vimoy+phaseCurr(k)*vi(k)
            else
                vi(k) = 0.d0
            end if
        end do
!
! 2.7 - RESTAURATION D ORIGINE VISQUEUSE
!
        cmoy = 0.d0
        mmoy = 0.d0
        do k = 1, nbPhase
            cmoy = cmoy+phaseCurr(k)*c(k)
            mmoy = mmoy+phaseCurr(k)*m(k)
        end do
        cr = cmoy*vimoy
        if (cr .le. 0.d0) then
            ds = 0.d0
        else
            ds = dt*(cr**mmoy)
        end if
        do k = 1, nbPhase
            if (phaseCurr(k) .gt. 0.d0) then
                vi(k) = vi(k)-ds
                if (vi(k) .le. 0.d0) vi(k) = 0.d0
            end if
        end do
! ----- Parameters for plasticity of tranformation
        trans = 0.d0
        if (l_plas_tran) then
            call metaGetParaPlasTransf('+', fami, 1, 1, imat, &
                                       metaType, nbPhase, deltaz, zalpha, &
                                       kpt, fpt)
            do k = 1, nbPhase-1
                if (deltaz(k) .gt. 0.d0) then
                    trans = trans+kpt(k)*fpt(k)*deltaz(k)
                end if
            end do
        end if
    else
        trans = 0.d0
        do k = 1, nbPhase
            vi(k) = vim(k)
        end do
    end if
!
! 2.9 - CALCUL DE HMOY ET RMOY (ON INCLUE LE SIGY)
!
    if (l_hard_isotline) then
! ----- Get hardening slope (linear)
        coef_hard = (1.d0)
        call metaGetParaHardLine(poum, fami, kpg, ksp, imat, &
                                 metaType, nbPhase, &
                                 e, coef_hard, h)
        do k = 1, nbPhase
            r(k) = h(k)*vi(k)+sy(k)
        end do
    end if
    if (l_hard_isotnlin) then
! ----- Get hardening slope (non-linear)
        call metaGetParaHardTrac(imat, metaType, nbPhase, &
                                 l_temp, temp, &
                                 vi, h, r, maxval)
        do k = 1, nbPhase
            r(k) = r(k)+sy(k)
        end do
    end if
    if (zalpha .gt. 0.d0) then
        rmoy = phaseCurr(1)*r(1)+phaseCurr(2)*r(2)+phaseCurr(3)*r(3)+phaseCurr(4)*r(4)
        rmoy = rmoy/zalpha
        hmoy = phaseCurr(1)*h(1)+phaseCurr(2)*h(2)+phaseCurr(3)*h(3)+phaseCurr(4)*h(4)
        hmoy = hmoy/zalpha
    else
        rmoy = 0.d0
        hmoy = 0.d0
    end if
    rmoy = (1.d0-fmel)*r(nbPhase)+fmel*rmoy
    hmoy = (1.d0-fmel)*h(nbPhase)+fmel*hmoy
!
! ********************************
! 3 - DEBUT DE L ALGORITHME
! ********************************
!
! 3.1 - JM=DET(FM),DJ=DET(DF),J=JM*DJ ET DFB
!
    jm = fm(1, 1)*(fm(2, 2)*fm(3, 3)-fm(2, 3)*fm(3, 2)) &
         -fm(2, 1)*(fm(1, 2)*fm(3, 3)-fm(1, 3)*fm(3, 2)) &
         +fm(3, 1)*(fm(1, 2)*fm(2, 3)-fm(1, 3)*fm(2, 2))
    dj = df(1, 1)*(df(2, 2)*df(3, 3)-df(2, 3)*df(3, 2)) &
         -df(2, 1)*(df(1, 2)*df(3, 3)-df(1, 3)*df(3, 2)) &
         +df(3, 1)*(df(1, 2)*df(2, 3)-df(1, 3)*df(2, 2))
    jp = jm*dj
    do i = 1, 3
        do j = 1, 3
            if (dj .le. 0.d0) then
                iret = 1
                goto 999
            else
                dfb(i, j) = (dj**(-1.d0/3.d0))*df(i, j)
            end if
        end do
    end do
!
! 3.2 - CONTRAINTES DE KIRSHHOFF A L INSTANT PRECEDENT
!
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
!
! 3.3 - DEFORMATIONS ELASTIQUES A L INSTANT PRECEDENT :
! BEM=DVTAUM/MUM+KR*TRBEM/3.D0
!
    xm = (jm**(-2.d0/3.d0))*(1.d0-2.d0*vim(vixSimo)/3.d0)
    do i = 1, 6
        bem(i) = dvtaum(i)/mum+kr(i)*xm
    end do
!
! 3.4 - BEL(I,J)=DFB(I,K)*BEM(K,L)*DFB(J,L)
!
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
!
! 3.5 - TRACE ET PARTIE DEVIATORIQUE DU TENSEUR BEL
!
    trbel = bel(1)+bel(2)+bel(3)
    do i = 1, 6
        dvbel(i) = bel(i)-kr(i)*trbel/3.d0
    end do
!
! 3.6 - CONTRAINTE ELASTIQUE (PARTIE DEVIATORIQUE)
!
    do i = 1, 6
        dvtel(i) = mu*dvbel(i)
    end do
!
! 3.7 - CONTRAINTE EQUIVALENTE ELASTIQUE ET SEUIL
!
    eqtel = 0.d0
    do i = 1, 6
        eqtel = eqtel+pdtsca(i)*dvtel(i)*dvtel(i)
    end do
    eqtel = sqrt(1.5d0*eqtel)

! 3.8 - TRACE DU TENSEUR DE KIRSHHOFF (CONNUE CAR NE DEPEND QUE DE J)
!
    trtaup = (troisk*((jp*jp)-1.d0)/6.d0)-(troisk*epsth*(jp+(1.d0/jp))/2.d0)
    dp = 0.d0
!
! ************************
! 4 - RESOLUTION
! ************************
!
    if (resi) then
!
! 4.2.1 - CALCUL DE DP
!
        seuil = eqtel-(1.d0+mu*trans*trbel)*rmoy
        if (seuil .lt. 0.d0) then
            vip(vixIndiPlas) = 0.d0
            dp = 0.d0
        else
            vip(vixIndiPlas) = 1.d0
            mutild = 2.d0*mu*trbel/3.d0
            call nzcalc(carcri, nbPhase, phaseCurr, zalpha, &
                        fmel, seuil, dt, trans, &
                        hmoy, mutild, eta, unsurn, &
                        dp, iret)
            if (iret .eq. 1) goto 999
!
! DANS LE CAS NON LINEAIRE
! VERIFICATION QU ON EST DANS LE BON INTERVALLE
!
            if (l_hard_isotnlin) then
                do j = 1, maxval
                    test = 0
                    vip(1:5) = vi(1:5)+dp
                    hplus(1:5) = h(1:5)
                    call metaGetParaHardTrac(imat, metaType, nbPhase, &
                                             l_temp, temp, &
                                             vip, h, r)
                    do k = 1, nbPhase
                        if (phaseCurr(k) .gt. 0.d0) then
                            r(k) = r(k)+sy(k)
                            if (abs(h(k)-hplus(k)) .gt. r8prem()) then
                                test = 1
                            end if
                        end if
                    end do
                    if (test .eq. 0) goto 600
                    hmoy = 0.d0
                    rmoy = 0.d0
                    if (zalpha .gt. 0.d0) then
                        do k = 1, nbPhase-1
                            if (phaseCurr(k) .gt. 0.d0) then
                                rmoy = rmoy+phaseCurr(k)*(r(k)-h(k)*dp)
                                hmoy = hmoy+phaseCurr(k)*h(k)
                            end if
                        end do
                        rmoy = fmel*rmoy/zalpha
                        hmoy = fmel*hmoy/zalpha
                    end if
                    if (phaseCurr(nbPhase) .gt. 0.d0) then
                        rmoy = (1.d0-fmel)*(r(nbPhase)-h(nbPhase)*dp)+rmoy
                        hmoy = (1.d0-fmel)*h(nbPhase)+hmoy
                    end if
                    seuil = eqtel-(1.d0+mu*trans*trbel)*rmoy
                    call nzcalc(carcri, nbPhase, phaseCurr, zalpha, &
                                fmel, seuil, dt, trans, &
                                hmoy, mutild, eta, unsurn, &
                                dp, iret)
                    if (iret .eq. 1) goto 999
                end do
                ASSERT((test .ne. 1) .or. (j .ne. maxval))
600             continue
            end if
        end if
!
! 4.2.2 - CALCUL DE SIGMA
!
        plasti = vip(vixIndiPlas)
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
!
! 4.2.3 - CALCUL DE VIP ET RMOY
        do k = 1, nbPhase
            if (phaseCurr(k) .gt. 0.d0) then
                vip(k) = vi(k)+dp
            else
                vip(k) = 0.d0
            end if
        end do
        vip(vixSigy) = 0.d0
        if (phaseCurr(nbPhase) .gt. 0.d0) then
            if (l_hard_isotline) then
                vip(vixSigy) = vip(vixSigy)+(1-fmel)*h(nbPhase)*vip(nbPhase)
            end if
            if (l_hard_isotnlin) then
                vip(vixSigy) = vip(vixSigy)+(1-fmel)*(r(nbPhase)-sy(nbPhase))
            end if
        end if
        if (zalpha .gt. 0.d0) then
            do k = 1, nbPhase-1
                if (l_hard_isotline) then
                    vip(vixSigy) = vip(vixSigy)+fmel*phaseCurr(k)*h(k)*vip(k)/zalpha
                end if
                if (l_hard_isotnlin) then
                    vip(vixSigy) = vip(vixSigy)+fmel*phaseCurr(k)*(r(k)-sy(k))/zalpha
                end if
            end do
        end if
    end if
!
! *******************************
! 5 - MATRICE TANGENTE DSIGDF
! *******************************
!
    if (rigi) then
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
        if (abs(plasti) .lt. 0.5d0) then
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
                    do i = 1, nbPhase
                        n0(i) = (1-n(i))/n(i)
                    end do
                    dv = (1-fmel)*phaseCurr(nbPhase)*(eta(nbPhase)/n(nbPhase)/dt)* &
                         ((dp/dt)**n0(nbPhase))
                    if (zalpha .gt. 0.d0) then
                        do i = 1, nbPhase-1
                            if (phaseCurr(i) .gt. 0.d0) then
                                dv = dv+fmel* &
                                     (phaseCurr(i)/zalpha)*(eta(i)/n(i)/dt)* &
                                     ((dp/dt)**n0(i))
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
!
! *******************************
! 6 - CORRECTION SUR TRBE
! *******************************
!
    if (resi) then
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
        vip(vixSimo) = 3.d0*(1.d0-(jp**(2.d0/3.d0))*xp)/2.d0
    end if
!
999 continue
!
end subroutine
