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
subroutine irmpga(nofimd, chanom, nochmd, typech, nomtyp, &
                  nbimpr, caimpi, caimpk, modnum, nuanom, &
                  sdcarm, carael, lfichUniq, field_type, codret)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterc/utflsh.h"
#include "asterfort/infniv.h"
#include "asterfort/irmase.h"
#include "asterfort/irmpg1.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/teattr.h"
#include "asterfort/uteref.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: nbimpr
    integer(kind=8) :: caimpi(10, nbimpr)
    integer(kind=8) :: modnum(MT_NTYMAX), nuanom(MT_NTYMAX, MT_NNOMAX)
    character(len=8) :: nomtyp(*)
    character(len=8) :: typech, sdcarm, carael
    character(len=16) :: tuyau, coque, grille, typmod2
    character(len=19) :: chanom
    character(len=80) :: caimpk(3, nbimpr)
    character(len=*) :: nofimd
    character(len=64) :: nochmd
    aster_logical :: lfichUniq
    character(len=16) :: field_type
    integer(kind=8) :: codret

!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE AU FORMAT MED - LOCALISATION POINTS DE GAUSS
!
! --------------------------------------------------------------------------------------------------
!
!     ENTREES :
!       NOFIMD : NOM DU FICHIER MED
!       CHANOM : NOM ASTER DU CHAMP
!       TYPECH : TYPE DU CHAMP ('ELEM', 'ELGA')
!       NOMTYP : NOM DES TYPES DE MAILLES ASTER
!       NBIMPR : NOMBRE D'IMPRESSIONS
!         CAIMPI : ENTIERS POUR CHAQUE IMPRESSION
!                  CAIMPI(1,I) = TYPE D'EF / MAILLE ASTER (0, SI NOEUD)
!                  CAIMPI(2,I) = NOMBRE DE POINTS (GAUSS OU NOEUDS)
!                  CAIMPI(3,I) = NOMBRE DE SOUS-POINTS
!                  CAIMPI(4,I) = NOMBRE DE COUCHES
!                  CAIMPI(5,I) = NOMBRE DE SECTEURS
!                  CAIMPI(6,I) = NOMBRE DE FIBTRES
!                  CAIMPI(7,I) = NOMBRE DE MAILLES A ECRIRE
!                  CAIMPI(8,I) = TYPE DE MAILLES ASTER (0, SI NOEUD)
!                  CAIMPI(9,I) = TYPE GEOMETRIQUE AU SENS MED
!                  CAIMPI(10,I) = NOMBRE TOTAL DE MAILLES IDENTIQUES
!       MODNUM : INDICATEUR SI LA SPECIFICATION DE NUMEROTATION DES
!                NOEUDS DES MAILLES EST DIFFERENTES ENTRE ASTER ET MED:
!                     MODNUM = 0 : NUMEROTATION IDENTIQUE
!                     MODNUM = 1 : NUMEROTATION DIFFERENTE
!       NUANOM : TABLEAU DE CORRESPONDANCE DES NOEUDS MED/ASTER.
!                NUANOM(ITYP,J): NUMERO DANS ASTER DU J IEME NOEUD DE LA
!                MAILLE DE TYPE ITYP DANS MED.
!     SORTIES :
!         CAIMPK : CARACTERES POUR CHAQUE IMPRESSION
!                  CAIMPK(1,I) = NOM DE LA LOCALISATION ASSOCIEE
!                  CAIMPK(2,I) = NOM DU PROFIL AU SENS MED
!                  CAIMPK(3,I) = NOM DE L'ELEMENT DE STRUCTURE
!       CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6), parameter :: nompro = 'IRMPGA'
    integer(kind=8) :: ifm, niv, nbcouc, nbsect, nummai
    integer(kind=8) :: iaux, jaux, kaux
    integer(kind=8) :: nbrepg, nbnoso, nbnoto, ndim
    integer(kind=8) :: ntypef, tygeom, tymast
    integer(kind=8) :: nbpg, nbsp
    integer(kind=8) :: nrimpr, iret
    integer(kind=8), parameter :: lgmax = 1000
    real(kind=8) :: refcoo(3*lgmax), gscoo(3*lgmax), wg(lgmax)
    real(kind=8) :: raux1(3*lgmax), raux2(3*lgmax), raux3(lgmax)
    aster_logical :: okgr, okcq, oktu, okpf
    character(len=4)  :: chnbco, chnbse
    character(len=10) :: nonuma
    character(len=16) :: nomtef, nomfpg, typsec
    character(len=16) :: valk(2)
    character(len=64) :: nolopg, nomasu
    real(kind=8) :: start, end
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    refcoo = 0.d0
    gscoo = 0.d0
    wg = 0.d0
    raux1 = 0.d0
    raux2 = 0.d0
    raux3 = 0.d0
    start = 0.d0
    end = 0.d0
!
    call infniv(ifm, niv)
!
    if (niv .gt. 1) then
        call cpu_time(start)
        write (ifm, 100) 'DEBUT DE '//nompro
        call utmess('I', 'MED_74')
        call utflsh(codret)
    end if
100 format(/, 4x, 10('='), a, 10('='),/)
!
!====
! 2. BOUCLAGE SUR LES IMPRESSIONS DEMANDEES
!    ON NE S'INTERESSE QU'A CELLES QUI COMPORTENT PLUS DE 1 POINT DE
!    GAUSS ET/OU PLUS DE 1 SOUS_POINT
!====
!
    do nrimpr = 1, nbimpr
        nbpg = caimpi(2, nrimpr)
        nbsp = caimpi(3, nrimpr)
!
        nomasu = ' '
!
        tymast = caimpi(8, nrimpr)
        tygeom = caimpi(9, nrimpr)
!
! 2.1.  CARACTERISATIONS DE L'ELEMENT FINI
!       ON DOIT RECUPERER LES COORDONNEES SOUS LA FORME :
!           ELEMENT 1D : X1 X2 ... ... XN
!           ELEMENT 2D : X1 Y1 X2 Y2 ... ... XN YN
!           ELEMENT 3D : X1 Y1 Z1 X2 Y2 Z2 ... ... XN YN ZN
!           C'EST CE QUE MED APPELLE LE MODE ENTRELACE. ON DOIT RECUPERER LES POIDS SOUS LA FORME :
!               WG1 WG2 ... ... WGN
!
        if (codret .eq. 0) then
            if (niv .gt. 1) then
                if (typech(1:4) .eq. 'ELGA') then
                    write (ifm, 210) typech, nbpg, nbsp
                else
                    write (ifm, 211) typech, nbpg, nbsp
                end if
            end if
210         format('CHAMP DE TYPE ', a, ', AVEC,', i5, ' POINTS DE GAUSS ET', i5, ' SOUS_POINTS')
211         format('CHAMP DE TYPE ', a, ', AVEC,', i5, ' POINTS ET', i5, ' SOUS_POINTS')
!
! 2.1.1.    CARACTERISATIONS DE L'ELEMENT FINI QUAND C'EST UN CHAMP
!           AUX POINTS DE GAUSS AVEC PLUS DE 1 POINT DE GAUSS
            if (typech(1:4) .eq. 'ELGA') then
!
! 2.1.1.1.      FORMULE GENERALE. CODE DE RETOUR :
!               CODRET = 0 : PAS DE PB
!               CODRET = 1 : LE CHAMP N'EST PAS DEFINI SUR CE TYPE D'ELEMENT
                ntypef = caimpi(1, nrimpr)
                call jenuno(jexnum('&CATA.TE.NOMTE', ntypef), nomtef)
                call teattr('C', 'TUYAU', tuyau, iret, typel=nomtef)
                call teattr('C', 'COQUE', coque, iret, typel=nomtef)
                call teattr('C', 'GRILLE', grille, iret, typel=nomtef)
                call teattr('C', 'TYPMOD2', typmod2, iret, typel=nomtef)
                if (tuyau .eq. 'OUI' .or. coque .eq. 'OUI' .or. grille .eq. 'OUI' &
                    .or. typmod2 .eq. 'PMF') then
                    if (nbsp .gt. 1 .and. carael .eq. ' ') then
                        valk(1) = field_type
                        valk(2) = nomtef
                        call utmess('F', 'MED2_14', nk=2, valk=valk)
                    end if
                end if
!
                call uteref(chanom, typech, ntypef, nomtef, lfichUniq, &
                            nomfpg, nbnoso, nbnoto, nbrepg, ndim, &
                            refcoo, gscoo, wg, nochmd, codret)
                if (codret .eq. 1) then
                    codret = 0
                    goto 20
                end if
!
                nbcouc = caimpi(4, nrimpr)
                nbsect = caimpi(5, nrimpr)
                nummai = caimpi(6, nrimpr)
                typsec = ' '
                ! CAS GRILLE, COQUE, TUYAU, PMF
                okgr = (nummai .eq. 0) .and. (nbcouc .eq. 1) .and. (nbsect .eq. 0) &
                       .and. (nbsp .eq. 1)
                okcq = (nummai .eq. 0) .and. (nbcouc .ge. 1) .and. (nbsect .eq. 0) &
                       .and. (nbsp .eq. 3*nbcouc)
                oktu = (nummai .eq. 0) .and. (nbcouc .ge. 1) .and. (nbsect .ge. 1)
                okpf = (nummai .ne. 0) .and. (nbcouc .eq. 0) .and. (nbsect .eq. 0)
                !
                if (okgr .or. okcq .or. oktu .or. okpf) then
                    if (okgr) then
                        typsec = 'GRILLE'
                        write (chnbco, '(I4)') nbcouc
                        write (chnbse, '(I4)') nbsp
                        nomasu = 'SECT_SHELL'//chnbco//'C '//chnbse//'P'
                    else if (okcq) then
                        typsec = 'COQUE'
                        write (chnbco, '(I4)') nbcouc
                        write (chnbse, '(I4)') nbsp
                        nomasu = 'SECT_SHELL'//chnbco//'C '//chnbse//'P'
                    else if (oktu) then
                        typsec = 'TUYAU'
                        write (chnbco, '(I4)') nbcouc
                        write (chnbse, '(I4)') nbsect
                        nomasu = 'SECT_PIPE '//chnbco//'C '//chnbse//'S'
                    else if (okpf) then
                        typsec = 'PMF'
                        write (nonuma, '(I10)') nummai
                        nomasu = 'SECT_BEAM '//nonuma
                    end if
                    call irmase(nofimd, typsec, nbcouc, nbsect, nummai, sdcarm, nomasu)
                    nbrepg = nbpg
                else
! 2.1.1.3           Extension avec des sous_points :
!                       on reproduit la meme description dans chaque 'couche'.
!                       c'est une solution temporaire, dans l'attente de l'evolution med
                    if (nbsp .gt. 1) then
                        do iaux = nbpg, 1, -1
                            do kaux = 1, nbsp
                                do jaux = ndim, 1, -1
                                    gscoo(ndim*(nbsp*(iaux-1)+kaux-1)+jaux) = &
                                        gscoo(ndim*(iaux-1)+jaux)
                                end do
                                wg(nbsp*(iaux-1)+kaux) = wg(iaux)
                            end do
                        end do
                    end if
                    nbrepg = nbpg*nbsp
                end if
!
! 2.1.2.    CARACTERISATIONS DE L'ELEMENT FINI QUAND C'EST UN CHAMP AUX POINTS DE GAUSS AVEC
!           UN SEUL POINT DE GAUSS MAIS PLUSIEURS SOUS_POINTS
!           ON DEFINIT UNE PSEUDO-LOCALISATION AUX POINTS DE GAUSS
!           ON DEDUIT LA DIMENSION DU CODAGE MED DU TYPE DE MAILLE SOUS-JACENTE
!
            else if ((typech(1:4) .eq. 'ELEM')) then
                iaux = mod(caimpi(9, nrimpr), 100)
                ndim = (caimpi(9, nrimpr)-iaux)/100
!
                ntypef = 0
                jaux = 3*lgmax
                do iaux = 1, jaux
                    refcoo(iaux) = 0.d0
                    gscoo(iaux) = 0.d0
                end do
                do iaux = 1, lgmax
                    wg(iaux) = 0.d0
                end do
!
                nomfpg(1:8) = nomtyp(tymast)
                do iaux = 1, 8
                    if (nomfpg(iaux:iaux) .eq. ' ') then
                        nomfpg(iaux:iaux) = '_'
                    end if
                end do
                nomfpg(9:16) = typech
!
                nbnoto = mod(tygeom, 100)
                nbrepg = nbsp
            end if
!
! 2.2.      ON ECRIT LA LOCALISATION
            if (typech(1:4) .eq. 'ELEM' .and. nbpg .eq. 1) then
                caimpk(1, nrimpr) = ' '
                caimpk(3, nrimpr) = ' '
            else if (ndim .gt. 0) then
                call irmpg1(nofimd, nomfpg, nbnoto, nbrepg, nbsp, &
                            ndim, tygeom, refcoo, gscoo, wg, &
                            raux1, raux2, raux3, nolopg, nomasu, &
                            lfichUniq, codret)
!
                caimpk(1, nrimpr) = nolopg
                caimpk(3, nrimpr) = nomasu
            end if
        end if
20      continue
    end do
!
! 3. LA FIN
    if (niv .gt. 1) then
        call cpu_time(end)
        write (ifm, *) '    ========== FIN DE '//nompro//' EN', end-start, 'sec ==========='
    end if
!
end subroutine
