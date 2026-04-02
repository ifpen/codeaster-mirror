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
subroutine te0339(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
#include "asterc/r8prem.h"
!
    character(len=*) :: option, nomte
!     FONCTION REALISEE :
!
!         CALCUL DU TAUX DE CROISSANCE DE CAVITES SELON UNE LOI DE
!         RICE ET TRACEY EN COMPORTEMENT NON-LINEAIRE.
!         ELEMENTS ISOPARAMETRIQUES 3D.
!
!         OPTION : 'RICE_TRACEY'
!
! ENTREE  --->  OPTION : NOM DE L'OPTION DE CALCUL
!         --->  NOMTE  : NOM DU TYPE D'ELEMENT
!
!
!-DEL CHARACTER*32 JEXNUM,JEXNOM,JEXATR,JEXR8
!
    character(len=16) :: optcal(2), rela_name
    character(len=16), pointer :: compor(:) => null()
!
    real(kind=8) :: sig(6), triax, volu, rsr0, numema, depseq
    real(kind=8) :: poids, dvpg, sigm, sigeq, lrsr0m, lrsr0p
    real(kind=8) :: cong(6), varigp, varigm, crois, vk
    integer(kind=8) :: jgano, nno, npg, i, kp, iritra, ndim, iret
    integer(kind=8) :: issopt, ima, iadzi, iazk24, nbvari, ipopp
    integer(kind=8) :: ipoids, ivf, idfde, nnos
    integer(kind=8) :: igeom, icong, ivarpg, ivarmg, isdrmr, isdrpr, jtab(7)
!
!======================== CORPS DU PROGRAMME ===========================
!
!     1. RECUPERATION DES INFOS
!     -------------------------
!     1.1 NOMBRE DE NOEUDS ET DE POINTS DE GAUSS
!     ------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
!     1.2 NUMERO DE LA MAILLE
!     -----------------------
    call tecael(iadzi, iazk24, noms=0)
    ima = zi(iadzi)
    numema = dble(ima)
!
!     1.3 CHAMPS IN
!     -------------
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTPR', 'L', icong)
    call jevech('PVARIMR', 'L', ivarmg)
    call jevech('PVARIPR', 'L', ivarpg)
    call jevech('PSDRMR', 'L', isdrmr)
    call jevech('PSOUSOP', 'L', issopt)
    call tecach('OOO', 'PVARIPR', 'L', iret, nval=7, &
                itab=jtab)
    nbvari = max(jtab(6), 1)*jtab(7)
    call jevech('PCOMPOR', 'L', vk16=compor)
    rela_name = compor(RELA_NAME)
!
    if ((rela_name .eq. 'VMIS_ISOT_TRAC') .or. (rela_name .eq. 'VMIS_ISOT_LINE') .or. &
        (rela_name .eq. 'LEMAITRE') .or. (rela_name .eq. 'VMIS_ECMI_TRAC') .or. &
        (rela_name .eq. 'VMIS_ECMI_LINE') .or. (rela_name .eq. 'VISC_CIN1_CHAB') .or. &
        (rela_name .eq. 'VISC_CIN2_CHAB')) then
        ipopp = 1
    else
        call utmess('F', 'ELEMENTS3_74', sk=rela_name)
    end if
!   /* ========================================================= */
!   /* PVARIMR = DEF PLAST EQ A L'INSTANT PRECEDENT              */
!   /* PVARIPR = DEF PLAST EQ A L'INSTANT COURRANT               */
!   /* PSDRMR  = LOG DU TAUX DE CROISSANCE A L'INSTANT PRECEDENT */
!   /* ========================================================= */
!
!     1.4 CHAMPS OUT
!     --------------
    call jevech('PRICTRA', 'E', iritra)
    call jevech('PSDRPR', 'E', isdrpr)
!   /* ========================================================= */
!   /* PRICTRA = CHAM_ELEM RICE-TRACEY (5 CMP)                   */
!   /* PSDRPR  = LOG DU TAUX DE CROISSANCE A L'INSTANT COURRANT  */
!   /* ========================================================= */
!
!     1.5 OPTIONS DE CALCUL
!     ---------------------
    optcal(1) = zk24(issopt) (1:16)
    optcal(2) = zk24(issopt) (17:19)
!
!     1.6 INITIALISATION
!     ------------------
    poids = 0.d0
    triax = 0.d0
    rsr0 = 0.d0
    volu = 0.d0
    vk = 0.d0
    dvpg = 0.d0
    depseq = 0.d0
    do i = 1, 6, 1
        cong(i) = 0.d0
    end do
    varigm = 0.d0
    varigp = 0.d0
!
!
!     2. BOUCLE SUR POINTS DE GAUSS SUIVANT OPTIONS DE CALCUL
!     -------------------------------------------------------
!   /* ============================================================== */
!   /* CMP ACTIVES DU CHAM_ELEM RICE-TRACEY SUIVANT LES OPTIONS       */
!   /* -------------------------------------------------------------- */
!   /* 1. CALCUL DU TAUX MOYEN : PREPRATION POUR INTEGRATION DE RT    */
!   /*       TAUX DE TRIAXIALITE SUR LA MAILLE               (TRIAX ) */
!   /*       VARIATION DE DEF PLAST EQ                       (DEPSEQ) */
!   /*       VOLUME PRIS EN COMPTE                           (VOLU  ) */
!   /*       NUMERO DE LA MAILLE                             (NUMEMA) */
!   /*    IE. CE QU'IL FAUDRA MOYENNER AVANT D'INTEGRER RT            */
!   /*                                                                */
!   /* 2. CALCUL DU TAUX MAX : INTEGRATION DE RT SUR LA MAILLE        */
!   /*       VOLUME PRIS EN COMPTE                           (VOLU  ) */
!   /*       TAUX DE CROISSANCE SUR LA MAILLE INSTANT COURRANT (RSR0) */
!   /*       NUMERO DE LA MAILLE                             (NUMEMA) */
!   /*    DANS CE CAS, LE PARAM OUT PSDRPR JOUE VRAIMENT SON ROLE     */
!   /* ============================================================== */
!
!     2.1 CHAM_ELEM POUR LE CALCUL DU TAUX MOYEN AVEC CHAMPS IN MOYENNES
!     ------------------------------------------------------------------
    if ((optcal(1) .eq. 'SIGM_ELMOY') .and. (optcal(2) .eq. 'NON')) then
!        2.1.1 INTEGRATION PAR QUADRATURE DES CHAMPS IN
!        ----------------------------------------------
        do kp = 1, npg, 1
            call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids)
            dvpg = poids
            vk = vk+dvpg
            do i = 1, 6, 1
                cong(i) = cong(i)+dvpg*zr(icong+6*kp+i-7)
            end do
            varigm = varigm+dvpg*zr(ivarmg+nbvari*(kp-1)+ipopp-1)
            varigp = varigp+dvpg*zr(ivarpg+nbvari*(kp-1)+ipopp-1)
        end do
!
!        2.1.2 VALEUR MOYENNE DES CHAMPS IN SUR LA MAILLE
!        ------------------------------------------------
        do i = 1, 6, 1
            sig(i) = cong(i)/vk
        end do
        varigm = varigm/vk
        varigp = varigp/vk
!
!        2.1.3 INVARIANTS
!        ----------------
        sigm = (sig(1)+sig(2)+sig(3))/3.d0
        sigeq = sig(4)*sig(4)+sig(5)*sig(5)+sig(6)*sig(6)
        sigeq = sigeq+sigeq
        sigeq = sigeq+( &
                sig(1)-sigm)*(sig(1)-sigm)+(sig(2)-sigm)*(sig(2)-sigm)+(sig(3)-sigm)*(sig(&
                &3)-sigm &
                )
        sigeq = sqrt(1.5d0*sigeq)
!
!        2.1.4 CHAMPS OUT
!        ----------------
        if (sigeq .gt. r8prem()) then
            triax = sigm/sigeq
        else
            triax = 0.d0
        end if
        volu = vk
        depseq = varigp-varigm
        do i = 1, npg, 1
            zr(isdrpr+i-1) = zr(isdrmr+i-1)
        end do
!
!     2.2 CHAM_ELEM POUR CALCUL DU TAUX MOYEN AVEC CHAMPS IN ORIGINAUX
!     ----------------------------------------------------------------
    else if ((optcal(1) .eq. 'SIGM_ELGA') .and. (optcal(2) .eq. 'NON')) &
        then
        do kp = 1, npg, 1
!           2.2.1 RECUPERATION DES CHAMPS IN
!           --------------------------------
            do i = 1, 6, 1
                cong(i) = zr(icong+6*kp+i-7)
            end do
            varigm = zr(ivarmg+nbvari*(kp-1)+ipopp-1)
            varigp = zr(ivarpg+nbvari*(kp-1)+ipopp-1)
!
!           2.2.2 CALCUL DE LA TRIAXIALITE LOCALE
!           -------------------------------------
            sigm = (cong(1)+cong(2)+cong(3))/3.d0
            sigeq = cong(4)*cong(4)+cong(5)*cong(5)+cong(6)*cong(6)
            sigeq = sigeq+sigeq
            sigeq = sigeq+( &
                    cong(1)-sigm)*(cong(1)-sigm)+(cong(2)-sigm)*(cong(2)-sigm)+(cong(3)-si&
                    &gm)*(cong(3)-sigm &
                    )
            sigeq = sqrt(1.5d0*sigeq)
!
!           2.2.3 INTEGRATION PAR QUADRATURE
!           --------------------------------
            call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids)
            dvpg = poids
            vk = vk+dvpg
            if (sigeq .gt. r8prem()) then
                triax = triax+dvpg*(sigm/sigeq)
            end if
            depseq = depseq+dvpg*(varigp-varigm)
        end do
!
!        2.2.4 CHAMPS OUT
!        ----------------
        triax = triax/vk
        volu = vk
        depseq = depseq/vk
        do i = 1, npg, 1
            zr(isdrpr+i-1) = zr(isdrmr+i-1)
        end do
!
!     2.3 CHAM_ELEM POUR LE CALCUL DU TAUX MAX AVEC CHAMPS IN MOYENNES
!     ----------------------------------------------------------------
    else if ((optcal(1) .eq. 'SIGM_ELMOY') .and. (optcal(2) .eq. 'OUI')) &
        then
!        2.3.1 INTEGRATION PAR QUADRATURE DES CHAMPS IN
!        ----------------------------------------------
        do kp = 1, npg, 1
            call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids)
            dvpg = poids
            vk = vk+dvpg
            do i = 1, 6, 1
                cong(i) = cong(i)+dvpg*zr(icong+6*kp+i-7)
            end do
            varigm = varigm+dvpg*zr(ivarmg+nbvari*(kp-1)+ipopp-1)
            varigp = varigp+dvpg*zr(ivarpg+nbvari*(kp-1)+ipopp-1)
        end do
!
!        2.3.2 VALEUR MOYENNE DES CHAMPS IN SUR LA MAILLE
!        ------------------------------------------------
        do i = 1, 6, 1
            sig(i) = cong(i)/vk
        end do
        varigm = varigm/vk
        varigp = varigp/vk
!
!        2.3.3 INVARIANTS
!        ----------------
        sigm = (sig(1)+sig(2)+sig(3))/3.d0
        sigeq = sig(4)*sig(4)+sig(5)*sig(5)+sig(6)*sig(6)
        sigeq = sigeq+sigeq
        sigeq = sigeq+( &
                sig(1)-sigm)*(sig(1)-sigm)+(sig(2)-sigm)*(sig(2)-sigm)+(sig(3)-sigm)*(sig(&
                &3)-sigm &
                )
        sigeq = sqrt(1.5d0*sigeq)
!
!        2.3.4 INTEGRATION DE LA LOI RT
!        ------------------------------
        if (sigeq .gt. r8prem()) then
            triax = sigm/sigeq
        else
            triax = 0.d0
        end if
        volu = vk
        depseq = varigp-varigm
        lrsr0m = zr(isdrmr)
        lrsr0p = lrsr0m+0.283d0*sign(1.0d0, triax)*exp(1.5d0*abs(triax))*depseq
!
!        2.3.5 CHAMPS OUT
!        ----------------
        rsr0 = exp(lrsr0p)
        do i = 1, npg, 1
            zr(isdrpr+i-1) = lrsr0p
        end do
!
!     2.4 CHAM_ELEM POUR LE CALCUL DU TAUX MAX AVEC CHAMPS IN ORIGINAUX
!     -----------------------------------------------------------------
    else if ((optcal(1) .eq. 'SIGM_ELGA') .and. (optcal(2) .eq. 'OUI')) &
        then
        do kp = 1, npg, 1
!           2.4.1 RECUPERATION DES CHAMPS IN
!           --------------------------------
            do i = 1, 6, 1
                cong(i) = zr(icong+6*kp+i-7)
            end do
            varigm = zr(ivarmg+nbvari*(kp-1)+ipopp-1)
            varigp = zr(ivarpg+nbvari*(kp-1)+ipopp-1)
            call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids)
            dvpg = poids
            volu = volu+dvpg
!
!           2.4.2 CALCUL DE LA TRIAXIALITE LOCALE
!           -------------------------------------
            sigm = (cong(1)+cong(2)+cong(3))/3.d0
            sigeq = cong(4)*cong(4)+cong(5)*cong(5)+cong(6)*cong(6)
            sigeq = sigeq+sigeq
            sigeq = sigeq+( &
                    cong(1)-sigm)*(cong(1)-sigm)+(cong(2)-sigm)*(cong(2)-sigm)+(cong(3)-si&
                    &gm)*(cong(3)-sigm &
                    )
            sigeq = sqrt(1.5d0*sigeq)
            if (sigeq .gt. r8prem()) then
                triax = sigm/sigeq
            else
                triax = 0.d0
            end if
!
!           2.4.3 INTEGRATION DE LA LOI RT AU PG COURRANT
!           ---------------------------------------------
            depseq = varigp-varigm
            lrsr0m = zr(isdrmr+kp-1)
            lrsr0p = lrsr0m+0.283d0*sign(1.0d0, triax)*exp(1.5d0*abs(triax))*depseq
            crois = exp(lrsr0p)
!
!           2.4.4 CHAMPS OUT
!           ----------------
            zr(isdrpr+kp-1) = lrsr0p
            if (crois .gt. rsr0) then
                rsr0 = crois
            end if
        end do
!       ON SORT LE VOLUME ASSOCIE A LA "SOUS-MAILLE" (CF DOC R)
!       PLUTOT QU ECRIRE UN POIDS RELATIF DE PT DE GAUSS,
!       VARIABLE SUIVANT LE PG QUI "ACCROCHE" LE MAX,
!       ON SORT LE SOUS-VOLUME MOYEN PAR PG :
        volu = volu/dvpg
!
!     2.5 TRAITEMENT DES OPTIONS INVALIDES
!     ------------------------------------
    else
!       OPTION DE CALCUL NON VALIDE
        ASSERT(.false.)
    end if
!
!
!     3. ECRITURE DES CMP DU CHAM_ELEM OUT DE TYPE RICE-TRACEY
!     --------------------------------------------------------
    zr(iritra) = triax
    zr(iritra+1) = rsr0
    zr(iritra+2) = volu
    zr(iritra+3) = numema
    zr(iritra+4) = depseq
!
end subroutine
