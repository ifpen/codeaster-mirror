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

subroutine vppara(modes, typcon, knega, lraide, lmasse, &
                  lamor, mxresf, neq, nfreq, omecor, &
                  dlagr, dbloq, vectr, vectc, nbpari, &
                  nbparr, nbpark, nopara, mod45, resui, &
                  resur, resuk, ktyp, lcomod, icom1, &
                  icom2, typres, nfreqg)
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lxlgut.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexis.h"
#include "asterfort/vpermc.h"
#include "asterfort/vpermo.h"
#include "asterfort/vpnorx.h"
#include "asterfort/vppcom.h"
#include "asterfort/vppfac.h"
#include "asterfort/vppgec.h"
#include "asterfort/vppgen.h"
#include "asterfort/vpstor.h"
#include "asterfort/wpermo.h"
#include "asterfort/wpnorx.h"
#include "asterfort/wppgen.h"
    character(len=4) :: mod45
    character(len=8) :: modes, knega
    character(len=1) :: ktyp
    character(len=16) :: typcon, typres
    character(len=*) :: resuk(*), nopara(*)
    integer(kind=8) :: lraide, lmasse, lamor, mxresf, neq, nfreq, dlagr(*), dbloq(*)
    integer(kind=8) :: resui(*)
    integer(kind=8) :: nbpari, nbparr, nbpark, icom1, icom2, nfreqg
    real(kind=8) :: vectr(*), resur(*), omecor
    complex(kind=8) :: vectc(*)
    aster_logical :: lcomod
!
!     CALCUL DES PARAMETRES MODAUX, DES NORMES D'ERREUR ET STOCKAGE DES
!     INFORMATIONS
!     ------------------------------------------------------------------
! IN  MODES    : K8  : NOM UTILISATEUR DU CONCEPT MODAL PRODUIT
! IN  KTYP     : K1  : TYPE DE LA MATRICE DE RAIDEUR
! IN  TYPCON   : K16 : TYPE DE LA STRUCTURE DE DONNEES PRODUITE
! IN  KNEGA    : K8  : VALEUR DU MOT-CLE NUME_MODE_NEGA
! IN  LRAIDE   : IS  : DESCRIPTEUR DE LA MATRICE DE "RAIDEUR"
! IN  LMASSE   : IS  : DESCRIPTEUR DE LA MATRICE DE "MASSE"
! IN  LAMOR    : IS  : DESCRIPTEUR DE LA MATRICE DE "AMORTISSEMENT"
! IN  MXRESF   : IS  : PARAMETRE DE DIMENSIONNEMENT DE RESUR
! IN  NEQ      : IS  : NPMBRE DE DDL
! IN/OUT  NFREQ    : IS  : NOMBRE DE FREQUENCES (OUT SI LCOMOD=.TRUE.)
! IN  FCORIG   : IS  : VALEUR SEUIL EN FREQUNCE
! IN DLAGR     : IS  : POSITION DES DDL DE LAGRANGE
! IN DBLOQ     : IS  : POSITION DES DDL BLOQUES
! IN/OUT VECTR : R   : VECTEURS PROPRES REELS
! IN/OUT VECTC : C   : VECTEURS PROPRES COMPLEXES
! OUT RESUR    : R   : STRUTURE DE DONNEES RESULTAT CONTENANT TOUS LES
!                      PARAMETRES MODAUX
! IN LCOMOD   : LOG  : CALC_MODES SUR PLUSIEURS SOUS-BANDES // OR NOT
! IN ICOM1/ICOM2 :  IN   : PARAMETRES // ASSOCIES A LCOMOD=.TRUE.
! IN TYPRES    :  IN  : TYPE DE RESULTATS (DYNAMIQUE OU FLAMB)
! IN NFREQG    : IN   : NBRE TOTAL DE MODES
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ineg, iprec, iret, ilgcon, nrscr, ifreq
    real(kind=8) :: rbid, fr2, fr3, rmin, raux
    complex(kind=8) :: zbid, caux
    aster_logical :: lns
!     ------------------------------------------------------------------
!
    call jemarq()
    zbid = (0.d0, 0.d0)
    rbid = 0.d0
    rmin = r8prem()
!
!     --- PRISE EN COMPTE DES MODES NEGATIFS ?
    ineg = +1
    if (knega .eq. 'OUI') ineg = -1
!
!     --- PREPARATION AU STOCKAGE DANS LA STRUCTURE DE DONNEES ---
    ilgcon = lxlgut(typcon)
    if (typcon(ilgcon-1:ilgcon) .eq. '_C') ilgcon = ilgcon-2
    call rsexis(modes, iret)
!
    if (lcomod) then
        nrscr = nfreqg
    else
        nrscr = nfreq
    end if
    if (iret .eq. 0) call rscrsd('G', modes, typcon(:ilgcon), nrscr)
    iprec = 0
!
!     --- MATRICE K ET/OU M NON SYMETRIQUE(S)
    if (zi(lraide+4)*zi(lmasse+4) .eq. 0) then
        lns = .true.
    else
        lns = .false.
    end if
!
!
!     --- NORMALISATION ET CALCUL DES PARAMETRES MODAUX ---
!     -----------------------------------------------------
!
    if ((lamor .eq. 0) .and. ((ktyp .eq. 'R') .and. (.not. lns))) then
!
! --- GENERALISE MODES REELS
!
        if (mod45 .ne. 'STAB') then
!
!        - NORMALISATION A LA + GRANDE DES COMPOSANTES /= LAGRANGE --
            call vpnorx(nfreq, neq, dlagr, vectr, resuk)
!
!
!        - CALCUL DES PARAMETRES GENERALISES ---
            call vppgen(lmasse, lamor, lraide, resur(4*mxresf+1), resur(6*mxresf+1), &
                        resur(5*mxresf+1), vectr, neq, nfreq, dbloq)
!
!
!        CALCUL DES FACTEURS DE PARTICIPATIONS ET DES MASSES EFFECTIVES
            call vppfac(lmasse, resur(4*mxresf+1), vectr, neq, nfreq, &
                        mxresf, resur(7*mxresf+1), resur(10*mxresf+1), &
                        resur(13*mxresf+1), resur(16*mxresf+1))
!
!
!        - CALCUL DE LA NORME D'ERREUR SUR LE MODE ---
            call vpermo(lmasse, lraide, nfreq, vectr, resur(mxresf+1), &
                        dbloq, omecor, resur(3*mxresf+1))
!
        end if
!        - SI PARALLELISME CALC_MODES SUR PLUSIEURS SOUS-BANDES, COM DES DONNEES:
!        - VECTEURS PROPRES ET DES RESUI/R/K
        call vppcom(lcomod, icom1, icom2, resui, resur, &
                    resuk, nbpari, nbparr, nbpark, mxresf, &
                    vectr, nfreq, neq, typres)
!
!        - STOCKAGE DES VECTEURS PROPRES ---
        call vpstor(ineg, 'R', modes, nfreq, neq, &
                    vectr, [zbid], mxresf, nbpari, nbparr, &
                    nbpark, nopara, mod45, resui, resur, &
                    resuk, iprec)
!
    else if ((lamor .eq. 0) .and. ((ktyp .eq. 'C') .or. lns)) then
        if (lcomod) then
            ASSERT(.false.)
        end if
! --- GENERALISE MODES COMPLEXES
!        - NORMALISATION A LA + GRANDE DES COMPOSANTES /= LAGRANGE --
        call wpnorx(nfreq, neq, dlagr, vectc, resuk)
!
!        - CALCUL DES PARAMETRES GENERALISES ---
        call vppgec(lmasse, lamor, lraide, resur(4*mxresf+1), resur(6*mxresf+1), &
                    resur(5*mxresf+1), vectc, neq, nfreq, dbloq)
!
!        CALCUL DES FACTEURS DE PARTICIPATIONS ET DES MASSES EFFECTIVES
!         CALL VPPFAC(LMASSE,RESUR(4*MXRESF+1),VECTR,NEQ,NFREQ,MXRESF,
!     &               RESUR(7*MXRESF+1),RESUR(10*MXRESF+1))
!
!        - CALCUL DE LA NORME D'ERREUR SUR LE MODE ---
        call vpermc(lmasse, lraide, nfreq, vectc, resur(mxresf+1), &
                    resur(2*mxresf+1), dbloq, omecor, resur(3*mxresf+1))
!
        if (typres(1:9) .eq. 'DYNAMIQUE') then
! Dans le cas GEP a modes complexes + typres=DYNAMIQUE , on stocke jusqu'ici
! avec om=(2*pi*f)**2, via rectfc:
! resufr(ivec, 2)=dble(om) et resufr(ivec,3)=(dimag(om)/dble(om))/2.d0
! Ceci afin detre coherent avec le calcul du residu de vpermc.
! On fige ici la bonne valeur pour l'impression et le sockage ds la SD_RESULTAT
            do ifreq = 1, nfreq
                fr2 = resur(mxresf+ifreq)
                fr3 = resur(2*mxresf+ifreq)
                caux = dcmplx(fr2, 2*fr2*fr3)
                caux = sqrt(caux)/r8depi()
                raux = abs(caux)
                if (raux .lt. rmin) then
                    raux = rmin
                end if
! valeur 'FREQUENCE' pour GEP a mode complexe (abs(lambda)): resufr(ifreq, 1)
                resur(ifreq) = raux
! valeur 'AMORTISSEMENT' pour GEP a mode complexe (Im(lambda)/abs(lambda)): resufr(ifreq,3)
                resur(2*mxresf+ifreq) = dimag(caux)/raux
            end do
        end if
!        - STOCKAGE DES VECTEURS PROPRES ---
        call vpstor(ineg, 'C', modes, nfreq, neq, &
                    [rbid], vectc, mxresf, nbpari, nbparr, &
                    nbpark, nopara, '    ', resui, resur, &
                    resuk, iprec)
!
    else if ((lamor .ne. 0) .and. (ktyp .eq. 'R')) then
        if (lcomod) then
            ASSERT(.false.)
        end if
! --- QUADRATIQUE MODES COMPLEXES AVEC K REELLE
!        - NORMALISATION A LA + GRANDE DES COMPOSANTES /= LAGRANGE --
        call wpnorx(nfreq, neq, dlagr, vectc, resuk)
!
!        - CALCUL DES PARAMETRES GENERALISES ---
        call wppgen(lmasse, lamor, lraide, resur(4*mxresf+1), resur(6*mxresf+1), &
                    resur(5*mxresf+1), vectc, neq, nfreq, dbloq)
!
!        - CALCUL DE LA NORME D'ERREUR SUR LE MODE ---
        call wpermo(lmasse, lraide, lamor, nfreq, vectc, &
                    resur(mxresf+1), resur(2*mxresf+1), dbloq, omecor, resur(3*mxresf+1))
!
!        - STOCKAGE DES VECTEURS PROPRES ---
        call vpstor(ineg, 'C', modes, nfreq, neq, &
                    [rbid], vectc, mxresf, nbpari, nbparr, &
                    nbpark, nopara, '    ', resui, resur, &
                    resuk, iprec)
    else if ((lamor .ne. 0) .and. (ktyp .eq. 'C')) then
        if (lcomod) then
            ASSERT(.false.)
        end if
        if (lns) then
            ASSERT(.false.)
        end if
! --- QUADRATIQUE MODES COMPLEXES AVEC K COMPLEXE
!        - NORMALISATION A LA + GRANDE DES COMPOSANTES /= LAGRANGE --
        call wpnorx(nfreq, neq, dlagr, vectc, resuk)
!
!        - CALCUL DES PARAMETRES GENERALISES ---
        call wppgen(lmasse, lamor, lraide, resur(4*mxresf+1), resur(6*mxresf+1), &
                    resur(5*mxresf+1), vectc, neq, nfreq, dbloq)
!
!        - CALCUL DE LA NORME D'ERREUR SUR LE MODE ---
        call wpermo(lmasse, lraide, lamor, nfreq, vectc, &
                    resur(mxresf+1), resur(2*mxresf+1), dbloq, omecor, resur(3*mxresf+1))
!
!        - STOCKAGE DES VECTEURS PROPRES ---
        call vpstor(ineg, 'C', modes, nfreq, neq, &
                    [rbid], vectc, mxresf, nbpari, nbparr, &
                    nbpark, nopara, '    ', resui, resur, &
                    resuk, iprec)
    end if
!     ------------------------------------------------------------------
    call jedema()
end subroutine
