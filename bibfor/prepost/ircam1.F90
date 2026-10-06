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
! aslint: disable=W1504
!
subroutine ircam1(nofimd, nochmd, existc, ncmprf, numpt, &
                  instan, numord, adsd, adsv, adsl, &
                  adsk, partie, indcmp, ncmpve, ntlcmp, &
                  ntncmp, ntucmp, ntproa, nbimpr, caimpi, &
                  caimpk, typech, nomamd, nomtyp, modnum, &
                  nuanom, lfichUniq, nosdfu, dschmd, nosdpr, &
                  codret)
!
    use as_med_module, only: as_med_open
    implicit none
!
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterc/utflsh.h"
#include "asterfort/as_mfdfin.h"
#include "asterfort/as_mficlo.h"
#include "asterfort/as_mfiodw_field.h"
#include "asterfort/infniv.h"
#include "asterfort/ircmcc.h"
#include "asterfort/ircmec.h"
#include "asterfort/ircmva.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: nbimpr
    integer(kind=8) :: caimpi(10, nbimpr)
    integer(kind=8) :: numpt, numord
    integer(kind=8) :: adsd, adsv, adsl, adsk
    integer(kind=8) :: existc, ncmprf
    integer(kind=8) :: ncmpve
    integer(kind=8) :: typent, tygeom
    integer(kind=8) :: modnum(MT_NTYMAX), nuanom(MT_NTYMAX, *)
    character(len=8) :: typech
    character(len=8) :: nomtyp(*)
    character(len=24) :: ntlcmp, ntncmp, ntucmp, ntproa, indcmp
    character(len=*) :: nofimd, partie
    character(len=*) :: nomamd, nosdpr
    character(len=*) :: caimpk(3, nbimpr)
    character(len=64) :: nochmd, dschmd
    real(kind=8) :: instan
    aster_logical :: lfichUniq
    character(len=8) :: nosdfu
    integer(kind=8) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE D'UN CHAMP - FORMAT MED - PHASE 1
!
! --------------------------------------------------------------------------------------------------
!
!     ENTREES :
!       NOFIMD : NOM DU FICHIER MED
!       NOCHMD : NOM MED DU CHAMP A ECRIRE
!       NCMPRF : NOMBRE DE COMPOSANTES DU CHAMP DE REFERENCE
!       NUMPT  : NUMERO DE PAS DE TEMPS
!       INSTAN : VALEUR DE L'INSTANT A ARCHIVER
!       NUMORD : NUMERO D'ORDRE DU CHAMP
!       ADSK, D, ... : ADRESSES DES TABLEAUX DES CHAMPS SIMPLIFIES
!       PARTIE: IMPRESSION DE LA PARTIE IMAGINAIRE OU REELLE POUR
!               UN CHAMP COMPLEXE
!       NCMPVE : NOMBRE DE COMPOSANTES VALIDES EN ECRITURE
!       NTPROA : PROFIL ASTER. C'EST LA LISTE DES NUMEROS ASTER DES
!                ELEMENTS/NOEUDS POUR LESQUELS LE CHAMP EST DEFINI
!       NBIMPR : NOMBRE D'IMPRESSIONS
!         NCAIMI : ENTIERS POUR CHAQUE IMPRESSION
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
!         NCAIMK : CARACTERES POUR CHAQUE IMPRESSION
!                  CAIMPK(1,I) = NOM DE LA LOCALISATION ASSOCIEE
!                  CAIMPK(2,I) = NOM DU PROFIL AU SENS MED
!                  CAIMPK(3,I) = NOM DE L'ELEMENT DE STRUCTURE
!       NOMAMD : NOM DU MAILLAGE MED
!       NOMTYP : NOM DES TYPES POUR CHAQUE MAILLE
!       MODNUM : INDICATEUR SI LA SPECIFICATION DE NUMEROTATION DES
!                NOEUDS DES MAILLES EST DIFFERENTES ENTRE ASTER ET MED:
!                     MODNUM = 0 : NUMEROTATION IDENTIQUE
!                     MODNUM = 1 : NUMEROTATION DIFFERENTE
!       NUANOM : TABLEAU DE CORRESPONDANCE DES NOEUDS MED/ASTER.
!                NUANOM(ITYP,J): NUMERO DANS ASTER DU J IEME NOEUD DE LA
!                MAILLE DE TYPE ITYP DANS MED.
!     SORTIES:
!       CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6), parameter :: nompro = 'IRCAM1'
    integer(kind=8) :: edleaj
    integer(kind=8), parameter :: ednoeu = 3
    integer(kind=8), parameter :: edmail = 0
    integer(kind=8), parameter :: ednoma = 4
    integer(kind=8), parameter :: typnoe = 0
    integer(kind=8), parameter :: ednopg = 1
    integer(kind=8), parameter :: edelst = 5
    character(len=8) :: saux08
    character(len=24) :: ntvale
    character(len=64) :: nomprf, nolopg, nomam2
    integer(kind=8) :: nbrepg, nbpt, iret
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nbenty, nvalec, nbpg, nbsp, nbco, nbse, nbfi
    integer(kind=8) :: tymast
    integer(kind=8) :: advale
    integer(kind=8) :: adproa, adnucm
    integer(kind=8) :: nrimpr, codre2, retsav
    integer(kind=8) :: ideb, ifin
    med_idt :: idfimd
    integer(kind=8) :: iaux
    aster_logical :: ficexi, lnvalec
    character(len=16), pointer :: cname(:) => null()
    character(len=16), pointer :: cunit(:) => null()
    real(kind=8) :: rbid(1)
    real(kind=8) :: start_time, end_time, start2, end2
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv > 1) then
        call cpu_time(start_time)
        write (ifm, *) "    ========== DEBUT DE IRCAM1 =========="
    end if
!
! 1.2. ==> NOMS DES TABLEAUX DE TRAVAIL
!
    ntvale = '&&'//nompro//'.VALEURS        '
!
! 1.3. ==> ADRESSES
!
    call jeveuo(ntproa, 'L', adproa)
    call jeveuo(ntlcmp, 'L', adnucm)
!
!====
! 2. OUVERTURE FICHIER MED EN MODE 'LECTURE_AJOUT'
!    CELA SIGNIFIE QUE LE FICHIER EST ENRICHI MAIS ON NE PEUT PAS
!    ECRASER UNE DONNEE DEJA PRESENTE.
!    (LE FICHIER EXISTE DEJA CAR SOIT ON A TROUVE UN MAILLAGE DEDANS,
!     SOIT ON EST PASSE PAR IRMAIL/IRMHDF)
!====
!
    call cpu_time(start2)
    inquire (file=nofimd, exist=ficexi)
    if (ficexi) then
        edleaj = 1
        if (lfichUniq) then
            call as_med_open(idfimd, nofimd, edleaj, codret, .true._1)
        else
            call as_med_open(idfimd, nofimd, edleaj, codret, .false._1)
        end if
        if (codret .ne. 0) then
            edleaj = 3
            if (lfichUniq) then
                call as_med_open(idfimd, nofimd, edleaj, codret, .true._1)
            else
                call as_med_open(idfimd, nofimd, edleaj, codret, .false._1)
            end if
        end if
    else
        edleaj = 3
        if (lfichUniq) then
            call as_med_open(idfimd, nofimd, edleaj, codret, .true._1)
        else
            call as_med_open(idfimd, nofimd, edleaj, codret, .false._1)
        end if
    end if
    if (codret .ne. 0) then
        saux08 = 'mfiope'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
    call cpu_time(end2)
    if (niv > 1) then
        call cpu_time(end_time)
        write (ifm, *) "    ========== IRCAM1 : OUVERTURE DU FICHIER EN ", &
            end2-start2, "sec ============"
    end if
!
!====
! 3. CREATION DU CHAMP
!====
!
! 3.1. ==> CREATION DU TABLEAUX DES COMPOSANTES
!
    nbpt = 0
    if (lfichUniq) then
        AS_ALLOCATE(vk16=cname, size=ncmpve)
        AS_ALLOCATE(vk16=cunit, size=ncmpve)
    else
        AS_ALLOCATE(vk16=cname, size=ncmprf)
        AS_ALLOCATE(vk16=cunit, size=ncmprf)
    end if
    nomam2 = ' '
    iret = 0
    call as_mfdfin(idfimd, nochmd, nomam2, nbpt, cunit(1), &
                   cname(1), iret)
    if (iret .eq. 0 .and. nbpt .ne. 0 .and. nomam2 .ne. nomamd) then
        call utmess('F', 'MED_94')
    end if
!
! 3.2. ==> CREATION DU CHAMP DANS LE FICHIER
!
    call ircmcc(idfimd, nomamd, nochmd, existc, ncmpve, &
                ntncmp, ntucmp, codret)
!
! 3.3. ==> DESCRIPTION DU CHAMP
!
    call as_mfiodw_field(idfimd, nochmd, dschmd, codret)
!
    if (codret .ne. 0) then
        saux08 = 'mfiodw'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
!====
! 4. ECRITURE POUR CHAQUE IMPRESSION SELECTIONNEE
!====
!
    ifin = 0
    retsav = 0
    do nrimpr = 1, nbimpr
        if (codret .eq. 0) then
            nvalec = caimpi(7, nrimpr)
            lnvalec = .false._1
            if (lfichUniq) then
                lnvalec = .true._1
            else
                if (nvalec .gt. 0) lnvalec = .true._1
            end if
            if (lnvalec) then
!
! 4.1. ==> ON DOIT ECRIRE DES VALEURS CORRESPONDANTS A NVALEC SUPPORTS
!          DU TYPE EN COURS.
!
                if (niv .gt. 1) then
                    write (ifm, 400)
                    call utflsh(codret)
                end if
!
                tygeom = caimpi(9, nrimpr)
                tymast = caimpi(8, nrimpr)
                ideb = ifin+1
                ifin = ideb+nvalec-1
!
                if (tygeom .eq. typnoe) then
                    typent = ednoeu
                else
                    if (typech .eq. 'ELNO') then
                        typent = ednoma
                    else
                        typent = edmail
                    end if
                end if
!
                if (tygeom .eq. typnoe) then
                    nbpg = 1
                    nbsp = 1
                    nbco = 0
                    nbse = 0
                    nbfi = 0
                    if (niv .gt. 1) then
                        write (ifm, 401)
                    end if
                else
                    nbpg = caimpi(2, nrimpr)
                    nbsp = caimpi(3, nrimpr)
                    nbco = caimpi(4, nrimpr)
                    nbse = caimpi(5, nrimpr)
                    nbfi = caimpi(6, nrimpr)
                    if (niv .gt. 1) then
                        write (ifm, 402) nomtyp(tymast), tygeom
                    end if
                end if
!
                nbenty = caimpi(10, nrimpr)
                nolopg = caimpk(1, nrimpr)
                nomprf = caimpk(2, nrimpr)
!
! 4.2. ==> CREATION DES TABLEAUX DE VALEURS A ECRIRE
!
                if (codret .eq. 0 .and. nvalec .ne. 0) then
!
                    iaux = ncmpve*nbsp*nbpg*nvalec
                    call wkvect(ntvale, 'V V R', iaux, advale)
!
                    call ircmva(zi(adnucm), indcmp, ncmpve, ncmprf, nvalec, &
                                nbpg, nbsp, adsv, adsd, adsl, &
                                adsk, partie, tymast, modnum, nuanom, &
                                typech, zr(advale), zi(adproa), ideb, ifin, &
                                codre2)
                    if (codre2 .ne. 0) retsav = 100
!
                end if
!
! 4.4. ==> ECRITURE VRAIE
!
                if (codret .eq. 0) then
!
                    nbrepg = ednopg
                    if ((tygeom .ne. typnoe) .and. (nbpg*nbsp .ne. 1)) then
                        nbrepg = nbpg
                        if (nbco .eq. 0 .and. nbse .eq. 0 .and. nbfi .eq. 0) then
                            nbrepg = nbpg*nbsp
                        else
                            typent = edelst
                        end if
                    end if
!
                    if (nvalec .ne. 0) then
                        call ircmec(idfimd, nochmd, nomprf, nolopg, numpt, &
                                    instan, numord, zr(advale), ncmpve, nbenty, &
                                    nbrepg, nvalec, typent, tygeom, nosdfu, &
                                    tymast, nrimpr, nosdpr, codret)
                    else
                        call ircmec(idfimd, nochmd, nomprf, nolopg, numpt, &
                                    instan, numord, rbid, ncmpve, nbenty, &
                                    nbrepg, nvalec, typent, tygeom, nosdfu, &
                                    tymast, nrimpr, nosdpr, codret)
                    end if
!
                    call jedetr(ntvale)
!
                end if
!
            end if
!
        end if
!
    end do
!
    if (niv .gt. 1) then
        call utflsh(codret)
        write (ifm, 400)
    end if
!
400 format(/, 80('-'),/)
401 format('  * VALEURS AUX NOEUDS',/)
402 format(/, '  * VALEURS SUR LES MAILLES DE TYPE ASTER ', a, ' ET MED', i4)
!
!====
! 5. FERMETURE DU FICHIER MED
!====
!
    call as_mficlo(idfimd, codret)
    if (codret .ne. 0) then
        saux08 = 'mficlo'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
    if (retsav .eq. 100) codret = 100
!
!====
! 6. LA FIN
!====
!
    AS_DEALLOCATE(vk16=cname)
    AS_DEALLOCATE(vk16=cunit)
    if (niv > 1) then
        call cpu_time(end_time)
        write (ifm, *) "    ========== FIN DE IRCAM1 EN ", end_time-start_time, "sec ============"
    end if
!
end subroutine
