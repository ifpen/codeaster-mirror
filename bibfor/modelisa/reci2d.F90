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
! aslint: disable=W0413

subroutine reci2d(lirela, mesh, nnoeca, noebe, nbcnx, &
                  cxma, normal, itria, xbar, iproj, &
                  excent)
    implicit none
!  DESCRIPTION : DETERMINATION DES RELATIONS CINEMATIQUES ENTRE LES DDLS
!  -----------   D'UN NOEUD DU CABLE ET LES DDLS DES NOEUDS VOISINS DE
!                LA STRUCTURE BETON
!                CAS OU LA STRUCTURE BETON EST MODELISEE PAR DES
!                ELEMENTS 2D
!                APPELANT : PROJCA
!
!  IN     : LIRELA : CHARACTER*19 , SCALAIRE
!                    NOM DE LA SD DE TYPE LISTE_DE_RELATIONS
!  IN     : MAILLA : CHARACTER*8 , SCALAIRE
!                    NOM DU CONCEPT MAILLAGE ASSOCIE A L'ETUDE
!  IN     : NNOECA : CHARACTER*8 , SCALAIRE
!                    NOM DU NOEUD DU CABLE
!  IN     : NOEBE  : INTEGER , SCALAIRE
!                    NUMERO DU NOEUD VOISIN DE LA STRUCTURE BETON LE
!                    PLUS PROCHE DU NOEUD DU CABLE
!  IN     : NBCNX  : INTEGER , SCALAIRE
!                    NOMBRE DE NOEUDS DE LA MAILLE VOISINE DE LA
!                    STRUCTURE BETON
!  IN     : CXMA   : INTEGER , VECTEUR DE DIMENSION AU PLUS NNOMAX
!                    CONTIENT LES NUMEROS DES NOEUDS DE LA MAILLE
!                    VOISINE DE LA STRUCTURE BETON
!                    (TABLE DE CONNECTIVITE)
!  IN     : NORMAL : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DANS LE REPERE GLOBAL DU VECTEUR NORMAL
!                    AU PLAN MOYEN DE LA MAILLE VOISINE DE LA STRUCTURE
!                    BETON
!  IN     : ITRIA  : INTEGER , SCALAIRE
!                    INDICATEUR DU SOUS-DOMAINE AUQUEL APPARTIENT LE
!                    POINT PROJETE :
!                    ITRIA = 1 : TRIANGLE 1-2-3
!                    ITRIA = 2 : TRIANGLE 3-4-1
!  IN     : XBAR   : REAL*8 , VECTEUR DE DIMENSION 3
!                    SI IPROJ.NE.2 : COORDONNEES BARYCENTRIQUES DU POINT
!                    PROJETE (BARYCENTRE DES SOMMETS DU TRIANGLE 1-2-3
!                    OU 3-4-1)
!  IN     : IPROJ  : INTEGER , SCALAIRE
!                    INDICE DE PROJECTION
!                    IPROJ =  0  LE POINT PROJETE EST A L'INTERIEUR
!                                DE LA MAILLE VOISINE
!                    IPROJ =  1X LE POINT PROJETE EST SUR UNE FRONTIERE
!                                DE LA MAILLE VOISINE
!                    IPROJ =  2  LE POINT PROJETE COINCIDE AVEC UN DES
!                                NOEUDS DE LA MAILLE VOISINE
!  IN     : EXCENT : REAL*8 , SCALAIRE
!                    EXCENTRICITE DU NOEUD DU CABLE PAR RAPPORT A LA
!                    MAILLE VOISINE DE LA STRUCTURE BETON
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/afrela.h"
#include "asterfort/ante2d.h"
#include "asterfort/elrfvf.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeexin.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
!
! ARGUMENTS
! ---------
    character(len=19) :: lirela
    character(len=8) :: mesh, nnoeca
    integer(kind=8) :: noebe, nbcnx, cxma(*), itria, iproj
    real(kind=8) :: normal(*), xbar(*), excent
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: i1, i2, i3, ibloc, icnx, iterm, ind, ier
    integer(kind=8) :: nbbloc, nbsom, nbterm, nbtmax, nnomax, noeca
    real(kind=8) :: ksi1, ksi2, zero
    complex(kind=8) :: cbid
    character(len=8) :: k8b
    aster_logical :: notlin, l_excent, lcolleNd, lcolleEl
!
    real(kind=8) :: ffel2d, x(2), ff(12)
    real(kind=8), pointer :: coemur(:) => null()
    integer(kind=8), pointer :: dimens(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    character(len=8), pointer :: nomddl(:) => null()
    character(len=8), pointer :: nomnoe(:) => null()
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call jemarq()
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   CREATION DES OBJETS DE TRAVAIL - INITIALISATIONS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    nnomax = 9
    nbtmax = 1+2*nnomax
    AS_ALLOCATE(vr=coemur, size=nbtmax)
    AS_ALLOCATE(vk8=nomddl, size=nbtmax)
    AS_ALLOCATE(vk8=nomnoe, size=nbtmax)
    AS_ALLOCATE(vi=dimens, size=nbtmax)
    AS_ALLOCATE(vr=direct, size=3*nbtmax)
    lcolleNd = .false.
    call jeexin(mesh(1:8)//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolleNd = .true.
    end if
    lcolleEl = .false.
    call jeexin(mesh(1:8)//'.NOMMAI', ier)
    if (ier .ne. 0) then
        lcolleEl = .true.
    end if
!
    notlin = (nbcnx .gt. 4)
    if ((nbcnx .eq. 3) .or. (nbcnx .eq. 6)) then
        nbsom = 3
    else
        nbsom = 4
    end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   DETERMINATION DE L'ANTECEDENT DU POINT PROJETE DANS L'ELEMENT
!     DE REFERENCE ASSOCIE A L'ELEMENT REEL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    if (iproj .ne. 2) call ante2d(itria, xbar(1), ksi1, ksi2)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 3   DETERMINATION DES RELATIONS CINEMATIQUES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    zero = 0.0d0
!
    nomnoe(1) = nnoeca
    nomddl(1) = 'DEPL'
    coemur(1) = 1.0d0

    l_excent = .true.
    if (excent .eq. 0.0d0) l_excent = .false.

!
    if (iproj .eq. 2) then
!
!       pas de liaisons si les noeuds sont les mêmes
        noeca = char8_to_int(nnoeca, lcolleNd, mesh, "NOEUD")
        if (noeca .eq. noebe) goto 110
!
        nbterm = 2
        nomnoe(1+1) = int_to_char8(noebe, lcolleNd, mesh, "NOEUD")
        nomddl(1+1) = 'DEPL'
        coemur(1+1) = -1.0d0
!
        if (l_excent) then
            nbterm = 3
            nomnoe(1+2) = nomnoe(1+1)
            nomddl(1+2) = 'ROTA'
            coemur(1+2) = -excent
        end if
!
    else
        if (nbcnx .eq. 3) then
            ff(1) = 0.5d0*(1.0d0+ksi2)
            ff(2) = -0.5d0*(ksi1+ksi2)
            ff(3) = 0.5d0*(1.0d0+ksi1)
        else if (nbcnx .eq. 6) then
            ff(1) = 0.5d0*(1.0d0+ksi2)*ksi2
            ff(2) = 0.5d0*(ksi1+ksi2)*(ksi1+ksi2+1.0d0)
            ff(3) = 0.5d0*(1.0d0+ksi1)*ksi1
            ff(4) = -1.0d0*(1.0d0+ksi2)*(ksi1+ksi2)
            ff(5) = -1.0d0*(1.0d0+ksi1)*(ksi1+ksi2)
            ff(6) = (1.0d0+ksi1)*(1.0d0+ksi2)
        else
            x(1) = ksi1
            x(2) = ksi2
            if (nbcnx .eq. 4) then
                call elrfvf('QU4', x, ff)
            else if (nbcnx .eq. 8) then
                call elrfvf('QU8', x, ff)
            else if (nbcnx .eq. 9) then
                call elrfvf('QU9', x, ff)
            else if (nbcnx .eq. 12) then
                call elrfvf('Q12', x, ff)
            end if
            ffel2d = ff(1)
            ff(1) = ff(4)
            ff(4) = ff(3)
            ff(3) = ff(2)
            ff(2) = ffel2d
            if (nbcnx .ge. 8) then
                ffel2d = ff(5)
                ff(5) = ff(8)
                ff(8) = ff(7)
                ff(7) = ff(6)
                ff(6) = ffel2d
            end if
        end if
!
        if (iproj .gt. 10) then
!
            nbterm = 3
            i1 = iproj-10
            i2 = i1+1
            if (i2 .gt. nbsom) i2 = 1
            if (l_excent) then
                ind = 3
            else
                ind = 2
            end if
            nomnoe(1+1) = int_to_char8(cxma(i1), lcolleNd, mesh, "NOEUD")
            nomnoe(1+ind) = int_to_char8(cxma(i2), lcolleNd, mesh, "NOEUD")
            nomddl(1+1) = 'DEPL'
            nomddl(1+ind) = 'DEPL'
!
            coemur(1+1) = -ff(i1)
            coemur(1+ind) = -ff(i2)
!
            if (l_excent) then
                nbterm = 5
                nomnoe(1+2) = nomnoe(1+1)
                nomnoe(1+4) = nomnoe(1+ind)
                nomddl(1+2) = 'ROTA'
                nomddl(1+4) = 'ROTA'
                coemur(1+2) = excent*coemur(1+1)
                coemur(1+4) = excent*coemur(1+ind)
            end if
!
            if (notlin) then
                nbterm = nbterm+1
                i3 = i1+nbsom
                ind = 3
                if (l_excent) ind = 5
                nomnoe(1+ind) = int_to_char8(cxma(i3), lcolleNd, mesh, "NOEUD")
                nomddl(1+ind) = 'DEPL'
!
                coemur(1+ind) = -ff(i3)
!
                if (l_excent) then
                    nbterm = 7
                    nomnoe(1+6) = nomnoe(1+ind)
                    nomddl(1+6) = 'ROTA'
                    coemur(1+6) = excent*coemur(1+ind)
                end if
            end if
!
        else
!
            nbterm = 1+nbcnx
            if (l_excent) nbterm = 1+2*nbcnx
            do icnx = 1, nbcnx
                ind = icnx
                if (l_excent) ind = 2*icnx-1
                nomnoe(1+ind) = int_to_char8(cxma(icnx), lcolleNd, mesh, "NOEUD")
                nomddl(1+ind) = 'DEPL'
                coemur(1+ind) = -ff(icnx)
!
                if (l_excent) then
                    nomnoe(1+ind+1) = nomnoe(1+ind)
                    nomddl(1+ind+1) = 'ROTA'
                    coemur(1+ind+1) = excent*coemur(1+ind)
                end if
            end do
!
        end if
    end if
!
!....... UNE RELATION PAR DDL DE TRANSLATION DU NOEUD DU CABLE
!        .....................................................
!
!....... LE VECTEUR ZI(JDIME) DOIT ETRE REINITIALISE AFIN DE PRENDRE
!....... EN COMPTE LES DIFFERENTS COEFFICIENTS PAR DIRECTION DEFINIS
!....... DANS LE VECTEUR ZR(JDIREC)
!
    do iterm = 1, nbterm
        dimens(iterm) = 3
    end do
!
    nbbloc = nbterm-1
    ind = 3
    if (l_excent) then
        ind = 6
        nbbloc = (nbterm-1)/2
    end if
!
!....... COEFFICIENTS PAR DIRECTIONS POUR LA PREMIERE RELATION (DDL DX)
!....... PUIS AFFECTATION
!
    direct(1) = 1.0d0
    direct(1+1) = 0.0d0
    direct(1+2) = 0.0d0
    do ibloc = 1, nbbloc
        direct(1+3+ind*(ibloc-1)) = 1.0d0
        direct(1+3+ind*(ibloc-1)+1) = 0.0d0
        direct(1+3+ind*(ibloc-1)+2) = 0.0d0
        if (l_excent) then
            direct(1+3+6*(ibloc-1)+3) = 0.0d0
            direct(1+3+6*(ibloc-1)+4) = normal(3)
            direct(1+3+6*(ibloc-1)+5) = -normal(2)
        end if
    end do
!
    call afrela(coemur, [cbid], nomddl, nomnoe, dimens, &
                direct, nbterm, zero, cbid, k8b, &
                'REEL', 'REEL', 0.d0, lirela)
!
!....... COEFFICIENTS PAR DIRECTIONS POUR LA DEUXIEME RELATION (DDL DY)
!....... PUIS AFFECTATION
!
    direct(1) = 0.0d0
    direct(1+1) = 1.0d0
    direct(1+2) = 0.0d0
    do ibloc = 1, nbbloc
        direct(1+3+ind*(ibloc-1)) = 0.0d0
        direct(1+3+ind*(ibloc-1)+1) = 1.0d0
        direct(1+3+ind*(ibloc-1)+2) = 0.0d0
        if (l_excent) then
            direct(1+3+6*(ibloc-1)+3) = -normal(3)
            direct(1+3+6*(ibloc-1)+4) = 0.0d0
            direct(1+3+6*(ibloc-1)+5) = normal(1)
        end if
    end do
!
    call afrela(coemur, [cbid], nomddl, nomnoe, dimens, &
                direct, nbterm, zero, cbid, k8b, &
                'REEL', 'REEL', 0.d0, lirela)
!
!....... COEFFICIENTS PAR DIRECTIONS POUR LA TROISIEME RELATION (DDL DZ)
!....... PUIS AFFECTATION
!
    direct(1) = 0.0d0
    direct(1+1) = 0.0d0
    direct(1+2) = 1.0d0
    do ibloc = 1, nbbloc
        direct(1+3+ind*(ibloc-1)) = 0.0d0
        direct(1+3+ind*(ibloc-1)+1) = 0.0d0
        direct(1+3+ind*(ibloc-1)+2) = 1.0d0
        if (l_excent) then
            direct(1+3+6*(ibloc-1)+3) = normal(2)
            direct(1+3+6*(ibloc-1)+4) = -normal(1)
            direct(1+3+6*(ibloc-1)+5) = 0.0d0
        end if
    end do
!
    call afrela(coemur, [cbid], nomddl, nomnoe, dimens, &
                direct, nbterm, zero, cbid, k8b, &
                'REEL', 'REEL', 0.d0, lirela)
!
110 continue
!
! --- MENAGE
    AS_DEALLOCATE(vr=coemur)
    AS_DEALLOCATE(vk8=nomddl)
    AS_DEALLOCATE(vk8=nomnoe)
    AS_DEALLOCATE(vi=dimens)
    AS_DEALLOCATE(vr=direct)
!
    call jedema()
!
! --- FIN DE RECI2D.
end subroutine
