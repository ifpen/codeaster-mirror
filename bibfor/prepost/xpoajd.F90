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
subroutine xpoajd(elrefp, ino, nnop, lsn, lst, &
                  ninter, iainc, ncompa, typma, co, &
                  igeom, jdirno, nfiss, jheavn, ncompn, &
                  he, ndime, ndim, cmp, nbcmp, &
                  nfh, nfe, ddlc, ima, jconx1, &
                  jconx2, jcnsv1, jcnsv2, jcnsl2, nbnoc, &
                  inntot, inn, nnn, contac, lmeca, &
                  pre1, heavno, nlachm, lacthm, jbaslo, &
                  jlsn, jlst, jstno, ka, mu)
! aslint: disable=W1306,W1504
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/elelin.h"
#include "asterfort/iselli.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xdeffe.h"
#include "asterfort/xlacti.h"
#include "asterfort/xmoffc.h"
#include "asterfort/xmofhm.h"
#include "asterfort/xpoffo.h"
#include "MeshTypes_type.h"
!
    integer(kind=8) :: ino, nnop, igeom, ndim, ndime, ddlc, jdirno
    integer(kind=8) :: nbcmp, cmp(*), nfe, ima, jconx1, jconx2, jcnsv1, heavno(20, 3)
    integer(kind=8) :: jcnsv2, jcnsl2, nbnoc, inntot, iainc, contac, ncompa
    integer(kind=8) :: nfiss, he(nfiss), nfh, inn, nnn, ninter(4), jheavn, ncompn
    integer(kind=8) :: nlachm(2), lacthm(16)
    integer(kind=8) :: jbaslo, jlsn, jlst, jstno
    aster_logical :: lmeca
    character(len=8) :: elrefp, typma
    real(kind=8) :: co(3), lsn(nfiss), lst(nfiss), ka, mu
!
!
!     BUT:  CALCUL DES DEPLACEMENTS AUX SOMMENTS DES SOUS-ELEMENTS
!           ET REPORT DES LAGRANGES SI CONTACT
!
!   IN
!     ELREFP : ÉLÉMENT DE RÉFÉRENCE PARENT
!     INO   : NUMÉRO DU NOEUD OU DU POINT D'INTERSECTION
!     NNOP   : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
!     LSN    : LEVEL SETS NORMALES EN INO
!     LST    : LEVEL SET TANGENTE EN INO
!     NINTER : NOMBRE D'ARETES INTERSECTÉS DE L'ELEMENT PARENT
!     IAINC  : ADRESSE DE TOPOFAC.AI DE L'ELEMENT PARENT
!     TYPMA  : TYPE DE LA MAILLE PARENTE
!     CO     : COORDONNÉES INITIALES DE INO
!     IGEOM  : COORDONNÉES DES NOEUDS DE L'ÉLÉMENT PARENT
!     JDIRNO : ADRESSE DU TABLEAU DIRNO LOCAL
!     NFISS  : NOMBRE DE FISSURES "VUES" PAR L'ÉLÉMENT PARENT
!     HE     : VALEURS DE(S) FONCTION(S) HEAVISIDE SUR LE SOUS ÉLÉMENT
!     NDIME  : DIMENSION TOPOLOGIQUE DE LA MAILLE PARENT
!     NDIM   : DIMENSION DU MAILLAGE
!     CMP    : POSITION DES DDLS DE DEPL X-FEM DANS LE CHAMP_NO DE DEPL1
!     NBCMP  : NOMBRE DE COMPOSANTES DU CHAMP_NO DE DEPL1
!     NFH    : NOMBRE DE FONCTIONS HEAVISIDE (PAR NOEUD)
!     NFE    : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT (1 A 4)
!     DDLC   : NOMBRE DE DDL DE CONTACT DE L'ÉLÉMENT PARENT
!     IMA    : NUMERO DE MAILLE COURANTE PARENT
!     JCONX1 : ADRESSE DE LA CONNECTIVITE DU MAILLAGE SAIN
!              (CONNECTIVITE QUADRATIQUE SI LAGRANGES DE CONTACT
!              AUX ARETES)
!     JCONX2 : LONGUEUR CUMULEE DE LA CONNECTIVITE DU MAILLAGE SAIN
!              (CONNECTIVITE QUADRATIQUE SI LAGRANGES DE CONTACT
!              AUX ARETES)
!     JCNSV1 : ADRESSE DU CNSV DU CHAM_NO DE DEPLACEMENT 1
!     NBNOC  : NOMBRE DE NOEUDS CLASSIQUES DU MAILLAGE FISSURE
!     INN    : COMPTEUR LOCAL DU NOMBRE DE NOUVEAU NOEUDS CREES
!     NNN    : NOMBRE DE NOUVEAU NOEUDS A CREER SUR LA MAILLE PARENT
!     LMECA  : VRAI DANS LE CAS MECANIQUE (SINON CAS THERMIQUE)
!
!   OUT
!     JCNSV2 : ADRESSE DU CNSV DU CHAM_NO DE DEPLACEMENT 2
!     JCNSL2 : ADRESSE DU CNSL DU CHAM_NO DE DEPLACEMENT 2
!     INNTOT : COMPTEUR TOTAL DU NOMBRE DE NOUVEAU NOEUDS CREES
!      INN   : COMPTEUR LOCAL DU NOMBRE DE NOUVEAU NOEUDS CREES
!
!
!
!
    character(len=8) :: elrefc, elref2
    integer(kind=8) :: nnops, iaindec, nptint
    real(kind=8) :: ff(nnop), ffc(16), fe(4), crilsn, minlsn
    real(kind=8) :: r, theta, chpri(3), lagrs(9), laghm(3), lagrc(3*ndim)
    real(kind=8) :: ff2(8), press
    real(kind=8) :: fk(MT_NNOMAX3D, 3, 3)
    integer(kind=8) :: i, j, iad, ipos, ig, ino2, ndimc, idecv2, idecl2
    integer(kind=8) :: nnol, ngl(8), ibid, ifiss, fiss, npr(8), nlag
    integer(kind=8) :: lact(8), nlact, hea_se, alp
    aster_logical :: lpint, lcont, pre1
    parameter(crilsn=1.d-6)
!
!     ------------------------------------------------------------------
    call jemarq()
!
! --- LPINT EST VRAI SI LE NOEUD DU MAILLAGE X-FEM EST SUR LA FISSURE
! --- SI LA MAILLE PARENTE POSSEDE DES DDLS DE CONTACT, ON CALCULERA
! --- ALORS LES LAGRANGES DE CONTACT FROTTEMENT POUR CE NOEUDS
! --- ATTENTION, IL SERA VRAI SEULEMENT SI ON EST DU COTÉ ESCAVE.
    if (ino .lt. 1000) then
        lpint = .false.
        do ifiss = 1, nfiss
            if (lsn(ifiss) .eq. 0.d0) lpint = .true.
        end do
    else if (ino .gt. 1000 .and. ino .lt. 2000) then
        lpint = .true.
    else if (ino .gt. 2000) then
        lpint = .false.
        do ifiss = 1, nfiss
            if (abs(lsn(ifiss)) .lt. crilsn) lpint = .true.
        end do
    end if
!
    fiss = 1
    if (lpint) then
        minlsn = r8maem()
        do ifiss = 1, nfiss
!     ON DETECTE LA FISSURE CORESPONDANTE AU POINT D'INTERSECTION
            if (abs(lsn(ifiss)) .lt. minlsn .and. he(ifiss) .ne. 0) then
                minlsn = abs(lsn(ifiss))
                fiss = ifiss
            end if
        end do
        if (he(fiss) .eq. 1) then
            lpint = .false.
        end if
    end if
    if (ddlc .gt. 0) then
        ASSERT(lmeca)
    end if
    nptint = ninter(1)
    if (pre1) nptint = ninter(fiss)
    lcont = (ddlc .gt. 0) .and. (ndime .eq. ndim) .and. (nptint .gt. 0)
!
    inn = inn+1
    inntot = inntot+1
    ASSERT(inn .le. nnn)
!
    zi(jdirno-1+(2+nfiss)*(inn-1)+1) = ino
    zi(jdirno-1+(2+nfiss)*(inn-1)+2) = nbnoc+inntot
    do ifiss = 1, nfiss
        zi(jdirno-1+(2+nfiss)*(inn-1)+2+ifiss) = he(ifiss)
    end do
!
!     FF : FONCTIONS DE FORMES AU NOEUD SOMMET OU D'INTERSECTION
    call xpoffo(ndim, ndime, elrefp, nnop, igeom, &
                co, ff)
!
    if (pre1 .or. (nfe .gt. 0 .and. lmeca .and. .not. iselli(elrefp))) then
!       ON RECUPERE L'ELEMENT LINEAIRE ASSOCIE A L'ELEMENT PARENT
!       QUADRATIQUE ET LE NOMBRE DE NOEUDS SOMMETS
        call elelin(3, elrefp, elref2, ibid, nnops)
!
!       FF2 : FONCTIONS DE FORME AUX NOEUDS SOMMETS POUR INTERPOLER LE
!       CHAMP PRIMAL CORRESPONDANT À LA PRESSION POUR LE CAS HM-XFEM
        ff2(:) = 0.d0
        call xpoffo(ndim, ndime, elref2, nnops, igeom, &
                    co, ff2)
    end if
!
!     RQ : "NDIMC" CORRESPOND AU NOMBRE DE COMPOSANTE VECTORIELLE DU
!     CHAMP PRIMAL (DEPL EN MECA -> NDIM CMP / TEMP EN THERMIQUE
!     SOIT 1 CMP)
    if (lmeca) then
        ndimc = ndim
    else
        ndimc = 1
    end if
    chpri(:) = 0.d0
    fk(:, :, :) = 0.
    fe(:) = 0.
!
!   FK : FONCTIONS D'ENRICHISSEMENT
    if (nfe .gt. 0 .and. ndimc .eq. 1) then
        r = sqrt(lsn(1)**2+lst(1)**2)
        if (r .gt. r8prem()) then
!         LE POINT N'EST PAS SUR LE FOND DE FISSURE
            theta = he(1)*abs(atan2(lsn(1), lst(1)))
        else
!         LE POINT EST SUR LE FOND DE FISSURE :
!         L'ANGLE N'EST PAS DÉFINI, ON LE MET À ZÉRO
            theta = 0.d0
        end if
!
        call xdeffe(r, theta, fe)
    else if (nfe .gt. 0 .and. lmeca) then
        if (abs(lsn(1)) .lt. crilsn .and. lst(1) .le. 0) then
            if (he(1) .gt. 0) then
                call xcalfev_wrap(ndim, nnop, zr(jbaslo), zi(jstno), real(he(1), 8), &
                                  zr(jlsn), zr(jlst), zr(igeom), ka, mu, &
                                  ff, fk, face='MAIT', elref=elrefp, nnop2=nnops, &
                                  ff2=ff2)
            else
                call xcalfev_wrap(ndim, nnop, zr(jbaslo), zi(jstno), real(he(1), 8), &
                                  zr(jlsn), zr(jlst), zr(igeom), ka, mu, &
                                  ff, fk, face='ESCL', elref=elrefp, nnop2=nnops, &
                                  ff2=ff2)
            end if
        else
            call xcalfev_wrap(ndim, nnop, zr(jbaslo), zi(jstno), real(he(1), 8), &
                              zr(jlsn), zr(jlst), zr(igeom), ka, mu, &
                              ff, fk, elref=elrefp, nnop2=nnops, ff2=ff2)
        end if
    end if
!
!   CALCUL DE L IDENTIFIANT SOUS ELEMENT
    hea_se = xcalc_code(nfiss, he)
!
!     CALCUL DE L'APPROXIMATION DU CHAMP PRIMAL "CHPRI" (DEPLACEMENT
!     EN MECA / TEMPERATURE EN THERMIQUE)
    if (pre1) then
        do j = 1, nnop
!
!         ADRESSE DE LA 1ERE CMP DU CHAMP PRIMAL DU NOEUD INO
            iad = jcnsv1-1+nbcmp*(zi(jconx1-1+zi(jconx2+ima-1)+j-1)-1)
!
            ipos = 0
!
!         DDLS CLASSIQUES POUR LES DEPLACEMENTS
            do i = 1, ndimc
                ipos = ipos+1
                chpri(i) = chpri(i)+ff(j)*zr(iad+cmp(ipos))
            end do
!
!         ON ZAPPE LA POSITION DES DDLS DE PRESSION QUELQUE SOIT LA
!         NATURE DU NOEUD
            ipos = ipos+1
!
!         DDLS HEAVISIDE POUR LES DEPLACEMENTS
            do ig = 1, nfh
                do i = 1, ndimc
                    ipos = ipos+1
                    chpri(i) = chpri(i)+xcalc_heav(zi(jheavn-1+(j-1)*ncompn+ig), hea_se, zi(jhea&
                               &vn-1+(j-1)*ncompn+ncompn))*ff(j)*zr(iad+cmp(ipos))
                end do
!         ON ZAPPE LA POSITION DES DDLS HEAVISIDE DE PRESSION QUELQUE SOIT LA
!         NATURE DU NOEUD
                ipos = ipos+1
!
            end do
!
        end do
!         ON TRAITE ICI LES DDLS DE PRESSION (NOEUDS SOMMETS UNIQUEMENT)
        do i = 1, nnops
            npr(i) = zi(jconx1-1+zi(jconx2+ima-1)+i-1)
        end do
!
        press = 0.d0
!         DDLS CLASSIQUES POUR LA PRESSION DANS LE MASSIF
        do i = 1, nnops
            press = press+ff2(i)*zr(jcnsv1-1+nbcmp*(npr(i)-1)+cmp(ndim+1))
!         DDLS HEAVISIDE POUR LA PRESSION DANS LE MASSIF
            do ig = 1, nfh
                press = press+xcalc_heav( &
                        zi(jheavn-1+(i-1)*ncompn+ig), hea_se, &
                        zi( &
                        jheavn-1+(i-1)*ncompn+ncompn))*zr(jcnsv1-1+nbcmp*(npr(i)-1)+cmp(ndim+1+ig&
                                                                                       &*(ndim+1)) &
                                                          )*ff2(i &
                                                                )
            end do
        end do
    else
        do j = 1, nnop
!
!         ADRESSE DE LA 1ERE CMP DU CHAMP PRIMAL DU NOEUD INO
            iad = jcnsv1-1+nbcmp*(zi(jconx1-1+zi(jconx2+ima-1)+j-1)-1)
!
            ipos = 0
!
!         DDLS CLASSIQUES
            do i = 1, ndimc
                ipos = ipos+1
                chpri(i) = chpri(i)+ff(j)*zr(iad+cmp(ipos))
            end do
!
!         DDLS HEAVISIDE
            do ig = 1, nfh
                do i = 1, ndimc
                    ipos = ipos+1
                    chpri(i) = chpri(i)+xcalc_heav(zi(jheavn-1+(j-1)*ncompn+ig), hea_se, zi(jhea&
                               &vn-1+(j-1)*ncompn+ncompn))*ff(j)*zr(iad+cmp(ipos))
                end do
            end do
!
!         DDL ENRICHIS EN FOND DE FISSURE
            do ig = 1, nfe
                if (ndimc .eq. 1) then
                    ipos = ipos+1
                    chpri(i) = chpri(i)+fe(1)*ff(j)*zr(iad+cmp(ipos))
                else if (ndimc .gt. 1) then
                    do alp = 1, ndimc
                        ipos = ipos+1
                        do i = 1, ndimc
                            chpri(i) = chpri(i)+fk(j, alp, i)*zr(iad+cmp(ipos))
                        end do
                    end do
                end if
            end do
        end do
    end if
    if (pre1) then
!
!      CALCUL DES LAGRANGES DE CONTACT HM-XFEM
!      SEULEMENT POUR LES POINTS D'INTERSECTION
!
        laghm(:) = 0.d0
        lagrc(:) = 0.d0
        if (lpint .and. lcont) then
!
!      NOEUD(S) GLOBAUX PORTANT LE(S) LAMBDA(S)
            if (contac .ge. 2) then
! ---  FONCTIONS DE FORMES LINEAIRES POUR LE P2P1
                call elelin(contac, elrefp, elrefc, ibid, nnol)
                call xpoffo(ndim, ndime, elrefc, nnol, igeom, &
                            co, ff)
            else
                ASSERT(.false.)
            end if
! ---  FONCTIONS DE FORMES MODIFIÉES
            iaindec = iainc+ncompa*(fiss-1)
            nptint = ninter(fiss)
            call xlacti(typma, nptint, iaindec, lact, nlact)
            nlachm(1) = nlact
            lacthm(1:8) = lact(1:8)
            ASSERT(nlachm(1) .gt. 0 .or. nlachm(2) .gt. 0)
            call xmofhm(lacthm, nlachm, nnol, ff, ffc)
            do j = 1, nnol
                ngl(j) = zi(jconx1-1+zi(jconx2+ima-1)+j-1)
            end do
!
            if (contac .eq. 3) then
                do i = 1, 3
! ---  CALCUL DE PRE_FLU, LAG_FLI ET LAG_FLS
                    do j = 1, nnol
                        laghm(i) = laghm(i)+zr(jcnsv1-1+nbcmp*(ngl(j)-1) &
                                               +cmp(ndim+1+nfh*(ndim+1)+ &
                                                    (3+ndim)*(heavno(j, fiss)-1)+i)) &
                                   *ffc(j)
                    end do
                end do
!
                do i = 1, ndim
! ---  CALCUL DE LAGS_C, LAGS_F1 ET LAGS_F2 (DDLS LIES A LA PARTIE COHESIVE)
                    do j = 1, nnol
                        lagrc(i) = lagrc(i)+zr(jcnsv1-1+nbcmp*(ngl(j)-1) &
                                               +cmp(ndim+1+nfh*(ndim+1)+ &
                                                    (3+ndim)*(heavno(j, fiss)-1)+3+i)) &
                                   *ffc(j)
                    end do
                end do
            else if (contac .eq. 2) then
                do i = 1, 3
! ---  CALCUL DE PRE_FLU, LAG_FLI ET LAG_FLS
                    do j = 1, nnol
                        laghm(i) = laghm(i)+zr(jcnsv1-1+nbcmp*(ngl(j)-1) &
                                               +cmp(ndim+1+nfh*(ndim+1) &
                                                    +(3+3*ndim)*(heavno(j, fiss)-1)+i))*ffc(j)
                    end do
                end do
!
                do i = 1, ndim
! ---  CALCUL DE LAGS_C, LAGS_F1 ET LAGS_F2 (DDLS LIES A LA PARTIE COHESIVE MORTAR)
                    do j = 1, nnol
                        lagrc(i) = lagrc(i)+zr(jcnsv1-1+nbcmp*(ngl(j)-1) &
                                               +cmp(ndim+1+nfh*(ndim+1) &
                                                    +(3+3*ndim)*(heavno(j, fiss)-1)+3+i))*ffc(j)
                    end do
                end do
!
                do i = ndim+1, 2*ndim
! ---  CALCUL DE JUPS_C, JUPS_F1 ET JUPS_F2 (DDLS LIES A LA PARTIE COHESIVE MORTAR)
                    do j = 1, nnol
                        lagrc(i) = lagrc(i)+zr(jcnsv1-1+nbcmp*(ngl(j)-1) &
                                               +cmp(ndim+1+nfh*(ndim+1) &
                                                    +(3+3*ndim)*(heavno(j, fiss)-1)+3+i))*ffc(j)
                    end do
                end do
!
                do i = 2*ndim+1, 3*ndim
! ---  CALCUL DE MUS_C, MUS_F1 ET MUS_F2 (DDLS LIES A LA PARTIE COHESIVE MORTAR)
                    do j = 1, nnol
                        lagrc(i) = lagrc(i)+zr(jcnsv1-1+nbcmp*(ngl(j)-1) &
                                               +cmp(ndim+1+nfh*(ndim+1) &
                                                    +(3+3*ndim)*(heavno(j, fiss)-1)+3+i))*ffc(j)
                    end do
                end do
            end if
        end if
    else
!
!     CALCUL DES LAGRANGES DE CONTACT FROTTEMENT
!     SEULEMENT POUR LES POINTS D'INTERSECTION
!
        lagrs(:) = 0.d0
        if (lpint .and. lcont .and. nfiss .eq. 1) then
!
!       NOEUD(S) GLOBAUX PORTANT LE(S) LAMBDA(S)
            nptint = ninter(1)
            iaindec = iainc
            call xlacti(typma, nptint, iaindec, lact, nlact)
            ASSERT(nlact .gt. 0)
            if (contac .eq. 1 .or. contac .eq. 2) then
                nnol = nnop
            else if (contac .eq. 3) then
! --- FONCTIONS DE FORMES LINEAIRES POUR LE P2P1
                call elelin(contac, elrefp, elrefc, ibid, nnol)
                call xpoffo(ndim, ndime, elrefc, nnol, igeom, &
                            co, ff)
            end if
! --- FONCTIONS DE FORMES MODIFIÉES
            call xmoffc(lact, nlact, nnol, ff, ffc)
            do j = 1, nnol
                ngl(j) = zi(jconx1-1+zi(jconx2+ima-1)+j-1)
            end do
!
            do i = 1, ddlc
! --- CALCUL AVEC LES FF DE CONTACT FFC, LINÉAIRES ET MODIFIÉES
                do j = 1, nnol
                    lagrs(i) = lagrs(i)+zr(jcnsv1-1+nbcmp*(ngl(j)-1) &
                                           +cmp((1+nfh+nfe)*ndim+i))*ffc(j)
                end do
            end do
        end if
    end if
!
!       ECRITURE DANS LE .VALE2 POUR LE NOEUD INO2
    ino2 = nbnoc+inntot
    if (pre1) then
        idecv2 = jcnsv2-1+(4*ndimc+4)*(ino2-1)
        idecl2 = jcnsl2-1+(4*ndimc+4)*(ino2-1)
        do i = 1, ndimc
            zr(idecv2+i) = chpri(i)
            zl(idecl2+i) = .true.
        end do
!       ADRESSAGE DANS LE TABLEAU ZR DES DDLS DE PRESSIONS
        zr(idecv2+ndimc+1) = press
        zl(idecl2+ndimc+1) = .true.
        if (lpint .and. lcont) then
!       ADRESSAGE POUR LES LAGRANGES DE CONTACT HM-XFEM
            do i = 1, 3
                zr(idecv2+ndimc+1+i) = laghm(i)
                zl(idecl2+ndimc+1+i) = .true.
            end do
            do i = 1, ndim
                zr(idecv2+ndimc+1+3+i) = lagrc(i)
                zl(idecl2+ndimc+1+3+i) = .true.
            end do
            if (contac .eq. 2) then
                do i = ndim+1, 3*ndim
                    zr(idecv2+ndimc+1+3+i) = lagrc(i)
                    zl(idecl2+ndimc+1+3+i) = .true.
                end do
            end if
        end if
    else
        if (lmeca) then
!         POUR LA MECA
            idecv2 = jcnsv2-1+4*ndimc*(ino2-1)
            idecl2 = jcnsl2-1+4*ndimc*(ino2-1)
        else
!         POUR LA THERMIQUE
            idecv2 = jcnsv2-1+ndimc*(ino2-1)
            idecl2 = jcnsl2-1+ndimc*(ino2-1)
        end if
        do i = 1, ndimc
            zr(idecv2+i) = chpri(i)
            zl(idecl2+i) = .true.
            if (lpint .and. lcont .and. nfiss .eq. 1) then
!             POUR LES LAGRANGES DE CONTACT EN MECA
                nlag = ddlc/ndimc
                do j = 1, nlag
                    zr(idecv2+ndimc*j+i) = lagrs(ndimc*(j-1)+i)
                    zl(idecl2+ndimc*j+i) = .true.
                end do
            end if
        end do
    end if
!
    call jedema()
!
end subroutine
