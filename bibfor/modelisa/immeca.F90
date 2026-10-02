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
!
subroutine immeca(tablca, lirela, mesh, nbnobe, nunobe, &
                  icabl, nbnoca, xnoca, ynoca, znoca, &
                  ncncin, nmabet, gromai)
    implicit none
!  DESCRIPTION : IMMERSION DES NOEUDS D'UN CABLE DANS LE MAILLAGE BETON
!  -----------   ET DETERMINATION DES RELATIONS CINEMATIQUES ENTRE LES
!                DDLS DES NOEUDS DU CABLE ET LES DDLS DES NOEUDS VOISINS
!                DE LA STRUCTURE BETON
!                APPELANT : OP0180 , OPERATEUR DEFI_CABLE_BP
!
!                EN SORTIE ON AJOUTE DES LIGNES DANS LA TABLE RESULTAT
!                LES CASES RENSEIGNEES CORRESPONDENT AUX PARAMETRES
!                <MAILLE_BETON_VOISINE>, <NOEUD_BETON_VOISIN> ET
!                <INDICE_IMMERSION>
!                LA SD DE TYPE LISTE_DE_RELATIONS EST MISE A JOUR
!
!  IN     : TABLCA : CHARACTER*19
!                    NOM DE LA TABLE DECRIVANT LES CABLES
!  IN     : LIRELA : CHARACTER*19 , SCALAIRE
!                    NOM DE LA SD DE TYPE LISTE_DE_RELATIONS
!  IN     : MAILLA : CHARACTER*8 , SCALAIRE
!                    NOM DU CONCEPT MAILLAGE ASSOCIE A L'ETUDE
!  IN     : NBNOBE : INTEGER , SCALAIRE
!                    NOMBRE DE NOEUDS APPARTENANT A LA STRUCTURE BETON
!  IN     : NUNOBE : CHARACTER*19 , SCALAIRE
!                    NOM D'UN VECTEUR D'ENTIERS POUR STOCKAGE DES
!                    NUMEROS DES NOEUDS APPARTENANT A LA STRUCTURE BETON
!  IN     : ICABL  : INTEGER , SCALAIRE
!                    NUMERO DU CABLE
!  IN     : NBNOCA : INTEGER , VECTEUR DE DIMENSION NBCABL
!                    CONTIENT LES NOMBRES DE NOEUDS DE CHAQUE CABLE
!  IN     : XNOCA  : CHARACTER*19 , SCALAIRE
!                    NOM D'UN VECTEUR DE REELS POUR STOCKAGE DES
!                    ABSCISSES X DES NOEUDS APPARTENANT AUX CABLES
!  IN     : YNOCA  : CHARACTER*19 , SCALAIRE
!                    NOM D'UN VECTEUR DE REELS POUR STOCKAGE DES
!                    ORDONNEES Y DES NOEUDS APPARTENANT AUX CABLES
!  IN     : ZNOCA  : CHARACTER*19 , SCALAIRE
!                    NOM D'UN VECTEUR DE REELS POUR STOCKAGE DES
!                    COTES Z DES NOEUDS APPARTENANT AUX CABLES
!  IN     : NCNCIN : CHARACTER*24 ,
!                    OBJET CONNECTIVITE INVERSE POUR LES MAILLES BETON
!  IN     : NMABET : CHARACTER*24 ,
!                    OBJET CONTENANT LES MAILLES BETON
!  IN     : GROMAI   CHARACTER*24
!                    NOM DU VECTEUR CONTENANT LES PLUS GRANDS
!                    DIAMETRES DES MAILLES DE BETON SELON LES
!                    DIRECTIONS X, Y ET Z
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
!
! ARGUMENTS
! ---------
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/immeno.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reci3d.h"
#include "asterfort/jeexin.h"
#include "asterfort/tbajli.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: mesh
    character(len=19) :: lirela, nunobe, xnoca, ynoca, znoca, tablca
    integer(kind=8) :: nbnobe, icabl, nbnoca(*)
    character(len=24) :: ncncin, nmabet, gromai
    character(len=24) :: valk(2)
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: nselec, nbnob2
    parameter(nselec=5)
    integer(kind=8) :: ideca, immer, inob1, inob2, inobe, inoca, ipara, itetra, jcoor
    integer(kind=8) :: jnoca, jnunob, jxca
    integer(kind=8) :: jyca, jzca, nbcnx, nblign, nbno, nbpara, nnomax, noe
    integer(kind=8) :: noebe(nselec), numail, nbval, nbval2, iret, ibid, noebec
    real(kind=8) :: d2, d2min(nselec), dx, dy, dz, rbid, x3dca(3), d2_max
    real(kind=8) :: x3dca2(3), axe(3), xnorm, xnorm2, zero, xbar(4)
    real(kind=8) :: rayon
    real(kind=8) :: long, longcy, longca, d2minc
    integer(kind=8) :: ifm, niv, ier
    character(len=8) :: nnoec2, k8b, presen(2), k8vide, noancr(2)
    complex(kind=8) :: cbid
    character(len=3) :: k3b
    character(len=8) :: nnoeca, voisin(2)
    character(len=24) :: coorno, nonoca, nogrna(2)
    integer(kind=8) :: n1, ibe, jbe, jgmai
    aster_logical :: lcolleNd, lcolleEl
!
    character(len=24) :: param(3), parcr
    integer(kind=8), pointer :: cnx_maille(:) => null()
    real(kind=8), pointer :: d2_min_max(:) => null()
    integer(kind=8), pointer :: no_min_max(:) => null()
    real(kind=8), pointer :: xyz_noemai(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
    character(len=24), pointer :: tblp(:) => null()
    data param/'MAILLE_BETON_VOISINE    ',&
     &                     'NOEUD_BETON_VOISIN      ',&
     &                     'INDICE_IMMERSION        '/
    data parcr/'NOEUD_CABLE             '/
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call infmaj()
    call jemarq()
    call infniv(ifm, niv)
!
    cbid = (0.d0, 0.d0)
    rbid = 0.d0
    zero = 0.0d0
    longcy = zero
    longca = zero
    k8vide = '        '
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   ACCES AUX DONNEES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 1.1 OBJETS DU MAILLAGE
! ---

    call jeveuo(gromai, 'L', jgmai)
    dx = zr(jgmai)*.55
    dy = zr(jgmai+1)*.55
    dz = zr(jgmai+2)*.55
    d2_max = dx**2+dy**2+dz**2

    coorno = mesh//'.COORDO    .VALE'
    call jeveuo(coorno, 'L', jcoor)

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
!
! RECUPERATION DES MOTS-CLES
!
!     TRAITEMENT DU MOT-CLE 'CONE'
    call getvr8('CONE', 'RAYON', iocc=1, scal=rayon, nbret=nbval)
    call getvr8('CONE', 'LONGUEUR', iocc=1, scal=long, nbret=nbval2)
    if (nbval .eq. 0) then
        rayon = zero
    end if
    if (nbval2 .eq. 0) then
        long = zero
    end if
    presen(1) = k8vide
    presen(2) = k8vide
    call getvtx('CONE', 'PRESENT', iocc=1, nbval=2, vect=presen, &
                nbret=n1)
!
!
!     TRAITEMENT DU MOT-CLE 'NOEUD_ANCRAGE'
    noancr(1) = k8vide
    noancr(2) = k8vide
    call getvtx('DEFI_CABLE', 'NOEUD_ANCRAGE', iocc=icabl, nbval=2, vect=noancr, &
                nbret=n1)
!
!     TRAITEMENT DU MOT-CLE 'GROUP_NO_ANCRAGE'
    if (n1 .eq. 0) then
        call getvem(mesh, 'GROUP_NO', 'DEFI_CABLE', 'GROUP_NO_ANCRAGE', icabl, &
                    2, nogrna(1), ibid)
!
        call utnono(' ', mesh, 'NOEUD', nogrna(1), k8b, &
                    iret)
        if (iret .eq. 10) then
            call utmess('F', 'ELEMENTS_67', sk=nogrna(1))
        else if (iret .eq. 1) then
            valk(1) = nogrna(1)
            valk(2) = k8b
            call utmess('A', 'SOUSTRUC_87', nk=2, valk=valk)
        end if
        noancr(1) = k8b
!
        call utnono(' ', mesh, 'NOEUD', nogrna(2), k8b, &
                    iret)
        if (iret .eq. 10) then
            call utmess('F', 'ELEMENTS_67', sk=nogrna(2))
        else if (iret .eq. 1) then
            valk(1) = nogrna(2)
            valk(2) = k8b
            call utmess('A', 'SOUSTRUC_87', nk=2, valk=valk)
        end if
        noancr(2) = k8b
    end if
!
!
! 1.2 DONNEES RELATIVES AU CABLE
! ---
!.... NOMBRE DE NOEUDS
!
    nbno = nbnoca(icabl)
!
!.... NOMS DES NOEUDS
!
    call jeveuo(tablca//'.TBNP', 'L', vi=tbnp)
    nbpara = tbnp(1)
    nblign = tbnp(2)
    ideca = nblign-nbno
    call jeveuo(tablca//'.TBLP', 'L', vk24=tblp)
    do ipara = 1, nbpara
        if (tblp(1+4*(ipara-1)) .eq. parcr) then
            nonoca = tblp(1+4*(ipara-1)+2)
            call jeveuo(nonoca, 'L', jnoca)
            goto 11
        end if
    end do
11  continue
!
!.... COORDONNEES DES NOEUDS
!
    call jeveuo(xnoca, 'L', jxca)
    call jeveuo(ynoca, 'L', jyca)
    call jeveuo(znoca, 'L', jzca)
!
! 1.3 NUMEROS DES NOEUDS DE LA STRUCTURE BETON
! ---
    call jeveuo(nunobe, 'L', jnunob)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   IMMERSION DES NOEUDS DU CABLE DANS LA STRUCTURE BETON
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! 2.1 CREATION D'OBJETS DE TRAVAIL
! ---
!.... LES MAILLES APPARTENANT A LA STRUCTURE BETON SONT DES MAILLES
!.... TETRA4, TETRA10, PYRAM5, PYRAM13, PENTA6, PENTA15,
!.... HEXA8, HEXA20 OU HEXA27
!.... LA VERIFICATION A ETE EFFECTUEE EN AMONT PAR LA ROUTINE TOMABE
!.... LE NOMBRE DE NOEUDS MAXIMAL SUR UNE MAILLE VAUT DONC 27
!
    nnomax = 27
    AS_ALLOCATE(vr=xyz_noemai, size=3*nnomax)
    AS_ALLOCATE(vi=cnx_maille, size=nnomax)
!
    AS_ALLOCATE(vr=d2_min_max, size=nbnobe)
    AS_ALLOCATE(vi=no_min_max, size=nbnobe)
!
!.... CALCUL DE LA LONGUEUR TOTALE DU CABLE
!
    do inoca = 1, (nbno-1)
!
        nnoeca = zk8(jnoca+ideca+inoca-1)
!
!       DANS LE CAS OU LE NOEUD INITIAL CORRESPOND AU NOEUD FINAL DONNE
!       PAR L'UTILISATEUR ON DOIT INVERSER LES DIRECTIVES 'PRESENT'
        if (inoca .eq. 1) then
            if (nnoeca .eq. noancr(2)) then
                k8b = presen(2)
                presen(2) = presen(1)
                presen(1) = k8b
            end if
        end if
!
        x3dca(1) = zr(jxca+ideca+inoca-1)
        x3dca(2) = zr(jyca+ideca+inoca-1)
        x3dca(3) = zr(jzca+ideca+inoca-1)
!
        nnoec2 = zk8(jnoca+ideca+inoca-1+1)
        x3dca2(1) = zr(jxca+ideca+inoca-1+1)
        x3dca2(2) = zr(jyca+ideca+inoca-1+1)
        x3dca2(3) = zr(jzca+ideca+inoca-1+1)
!
        axe(1) = (x3dca2(1)-x3dca(1))
        axe(2) = (x3dca2(2)-x3dca(2))
        axe(3) = (x3dca2(3)-x3dca(3))
        xnorm2 = axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3)
!
        if (xnorm2 .eq. zero) then
            call utmess('F', 'MODELISA4_70')
        end if
!
        xnorm = sqrt(xnorm2)
        longca = longca+xnorm
!
    end do
!
    if (niv .eq. 2) then
        write (ifm, *) '------------------------------------------'
        write (ifm, *) ' DEFINITION DES RELATIONS CINEMATIQUES'
        if (rayon .eq. zero) then
            write (ifm, *) '  CONE : PAS DE CONE'
        else
            write (ifm, *) '  RAYON DU CONE : ', rayon
            write (ifm, *) '  LONGUEUR DU CONE : ', long
        end if
        write (ifm, *) '  LONGUEUR DU CABLE : ', longca
        write (ifm, *) ' '
    end if
!
!
!
! 2.2 BOUCLE SUR LE NOMBRE DE NOEUDS DU CABLE
! ---
    do inoca = 1, nbno
!
!
        nnoeca = zk8(jnoca+ideca+inoca-1)
        x3dca(1) = zr(jxca+ideca+inoca-1)
        x3dca(2) = zr(jyca+ideca+inoca-1)
        x3dca(3) = zr(jzca+ideca+inoca-1)
!
        nnoec2 = zk8(jnoca+ideca+inoca-1+1)
        x3dca2(1) = zr(jxca+ideca+inoca-1+1)
        x3dca2(2) = zr(jyca+ideca+inoca-1+1)
        x3dca2(3) = zr(jzca+ideca+inoca-1+1)
!
        if (niv .eq. 2) then
            write (ifm, *) ' '
            write (ifm, *) ' '
            if (inoca .lt. nbno) then
                write (ifm, *) 'NOEUDS CABLE : ', nnoeca, ' - ', nnoec2
            else
                write (ifm, *) 'NOEUD CABLE : ', nnoeca
            end if
        end if
!
!
! 2.2.0  CREATION DU VECTEUR AXE, RELIANT DEUX NOEUDS CABLES CONSECUTIFS
! .....  POUR LE CALCUL DES DISTANCES AU CYLINDRE
!
        if (inoca .ne. nbno) then
            axe(1) = (x3dca2(1)-x3dca(1))
            axe(2) = (x3dca2(2)-x3dca(2))
            axe(3) = (x3dca2(3)-x3dca(3))
            xnorm2 = axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3)
            xnorm = sqrt(xnorm2)
        else
            xnorm = 0.d0
        end if
!
! ... CHOIX DU TRAITEMENT :
!
!  SI LA LONGUEUR (OU LE RAYON) DU CONE EST NULLE (PAS DE CONE)
        if ((long .eq. zero) .or. (rayon .eq. zero)) goto 112
!
!  TESTE SI ON EST EN DEHORS DES ZONES DE DEFINITIONS DES TUNNELS
        if ((longcy .gt. long) .and. (longcy .lt. (longca-long))) goto 112
!
!  SINON TESTE SI ON A DEMANDE DES TUNNELS
!
        if ((longcy .lt. long) .and. (presen(1) (1:3) .eq. 'NON')) goto 112
!
        if ((longcy .gt. (longca-long)) .and. (presen(2) (1:3) .eq. 'NON')) goto 112
!
!  SINON ON DEFINI LE CONE
!
!     --------------------------------------------------------
!     CAS 1 : ON DEFINIT LE CONE POUR ATTACHER LES NOEUDS
!     --------------------------------------------------------
!      Note : on ne fait rien au niveau du fortran
!      (voir la macrocommande DEFI_CABLE_BP)
!
!
        if (niv .eq. 2) then
            write (ifm, *) '-> ON DEFINIT LE CYLINDRE D''AXE ', nnoeca, &
                ' - ', nnoec2
        end if
!
        longcy = longcy+xnorm
        goto 100
!
!
!     ---------------------------------------------------------
!     CAS 2 : ON ATTACHE UN SEUL NOEUD DU BETON AU NOEUD CABLE
!     ---------------------------------------------------------
!
112     continue
        if (niv .eq. 2) then
            write (ifm, *) '-> ON ATTACHE LE NOEUD OU LA MAILLE BETON '//&
     &             'LA PLUS PROCHE'
        end if
!
        longcy = longcy+xnorm
!
!
! 2.2.1  DETERMINATION DU NOEUD DE LA STRUCTURE BETON LE PLUS PROCHE
! .....  DU NOEUD CABLE COURANT
!
!
! FILTRE DES NOEUDS TROP LOIN
!
        nbnob2 = 0
        do inobe = 1, nbnobe
            noe = zi(jnunob+inobe-1)
            dx = x3dca(1)-zr(jcoor+3*(noe-1))
            dy = x3dca(2)-zr(jcoor+3*(noe-1)+1)
            dz = x3dca(3)-zr(jcoor+3*(noe-1)+2)
            d2 = dx*dx+dy*dy+dz*dz
            if (d2 .lt. d2_max) then
                nbnob2 = nbnob2+1
                d2_min_max(nbnob2) = d2
                no_min_max(nbnob2) = noe
            end if
        end do
!
! ON DETERMINE LES NSELEC NOEUDS LES PLUS PROCHES
!
        do ibe = 1, nselec
            d2min(ibe) = d2_max
            noebe(ibe) = 0
        end do
        do inobe = 1, nbnob2
            d2 = d2_min_max(inobe)
            if (d2 .gt. d2min(nselec)) goto 113
            do ibe = 1, nselec
                if (d2 .lt. d2min(ibe)) then
                    do jbe = 0, nselec-ibe-1
                        d2min(nselec-jbe) = d2min(nselec-jbe-1)
                        noebe(nselec-jbe) = noebe(nselec-jbe-1)
                    end do
                    d2min(ibe) = d2
                    noebe(ibe) = no_min_max(inobe)
                    goto 113
                end if
            end do
113         continue
        end do
!
        if (niv .eq. 2) then
            write (ifm, *) '   INFOS : DISTANCE MINIMALE : ', sqrt(d2min(1))
        end if
!
! 2.2.2  TENTATIVE D'IMMERSION DU NOEUD CABLE DANS LES MAILLES
! .....  AUXQUELLES APPARTIENT LE NOEUD BETON LE PLUS PROCHE
!
        noebec = 0
        immer = -1
        do ibe = 1, nselec
!          ATTENTION IL PEUT Y AVOIR MOINS QUE NSELEC NOEUDS
!          DE BETON
            if (noebe(ibe) .eq. 0) goto 116
!
            call immeno(ncncin, nmabet, mesh, x3dca(1), noebe(ibe), &
                        numail, nbcnx, cnx_maille, xyz_noemai, itetra, &
                        xbar(1), immer)
            if (immer .ge. 0) then
                noebec = noebe(ibe)
                goto 116
            end if
        end do
116     continue
!
! -  TEST DE COINCIDENCE GEOGRAPHIQUEMENT AVEC LE NOEUD BETON LE PLUS PROCHE
!
        if (immer > 0 .and. sqrt(d2min(1)) < 1d2*r8prem()*xnorm) then
            immer = 2
        end if
!
! 2.2.3  EN CAS D'ECHEC DE LA TENTATIVE PRECEDENTE
! .....
        if (immer .lt. 0) then
!
!.......... ON CREE UNE LISTE ORDONNEE DES NOEUDS DE LA STRUCTURE BETON
!.......... DU PLUS PROCHE AU PLUS ELOIGNE DU NOEUD CABLE CONSIDERE
!
            do inob1 = 1, nbnob2-1
                d2minc = d2_min_max(inob1)
                noebec = no_min_max(inob1)
                inobe = inob1
                do inob2 = inob1+1, nbnob2
                    if (d2_min_max(inob2) .lt. d2minc) then
                        d2minc = d2_min_max(inob2)
                        noebec = no_min_max(inob2)
                        inobe = inob2
                    end if
                end do
                if (inobe .gt. inob1) then
                    d2 = d2_min_max(inob1)
                    noe = no_min_max(inob1)
                    d2_min_max(inob1) = d2minc
                    no_min_max(inob1) = noebec
                    d2_min_max(inobe) = d2
                    no_min_max(inobe) = noe
                end if
            end do
!
            if (niv .eq. 2) then
                write (ifm, *) '   INFOS : DISTANCE MINIMALE : ', sqrt(d2)
            end if
!
!
!.......... LA TENTATIVE D'IMMERSION DANS LES MAILLES AUXQUELLES
!.......... APPARTIENT LE NOEUD BETON LE PLUS PROCHE A DEJA ETE
!.......... EFFECTUEE, SANS SUCCES
!.......... ON EFFECTUE DE NOUVELLES TENTATIVES EN UTILISANT LES NOEUDS
!.......... DE LA LISTE ORDONNEE PRECEDENTE, DU SECOND JUSQU'AU DERNIER
!.......... REPETER
            do inobe = nselec, nbnob2
                noebec = no_min_max(inobe)
!............. TENTATIVE D'IMMERSION DU NOEUD CABLE DANS LES MAILLES
!............. AUXQUELLES APPARTIENT LE NOEUD BETON COURANT
                call immeno(ncncin, nmabet, mesh, x3dca(1), noebec, &
                            numail, nbcnx, cnx_maille, xyz_noemai, itetra, &
                            xbar(1), immer)
!............. SORTIE DU BLOC REPETER EN CAS DE SUCCES
                if (immer .ge. 0) goto 131
            end do
131         continue
!
        end if
!
!
!
! 2.2.4  SORTIE EN ERREUR FATALE SI ECHEC PERSISTANT
! .....
        if (immer .lt. 0) then
            write (k3b, '(I3)') icabl
            valk(1) = k3b
            valk(2) = nnoeca
            call utmess('F', 'MODELISA4_71', nk=2, valk=valk)
        end if
!
! 2.2.5  DETERMINATION DES RELATIONS CINEMATIQUES
! .....
        call reci3d(lirela, mesh, nnoeca, noebec, nbcnx, &
                    cnx_maille, itetra, xbar(1), immer)
!
! 2.2.6  MISE A JOUR DE LA SD TABLE
! .....
        voisin(1) = int_to_char8(numail, lcolleEl, mesh, "MAILLE")
        ASSERT(noebec .ne. 0)
        voisin(2) = int_to_char8(noebec, lcolleNd, mesh, "NOEUD")
        call tbajli(tablca, 3, param, [immer], [rbid], &
                    [cbid], voisin(1), ideca+inoca)
!
100     continue
    end do
!
!  FIN BOUCLE SUR NBNO
!
!
! --- MENAGE
!
    AS_DEALLOCATE(vr=xyz_noemai)
    AS_DEALLOCATE(vi=cnx_maille)
    AS_DEALLOCATE(vr=d2_min_max)
    AS_DEALLOCATE(vi=no_min_max)
!
    call jedema()
!
! --- FIN DE IMMECA.
end subroutine
