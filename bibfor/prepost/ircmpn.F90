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

subroutine ircmpn(nofimd, ncmprf, ncmpve, numcmp, exicmp, &
                  nbvato, nbnoec, linoec, adsl, caimpi, &
                  caimpk, profas, innoce, nosdfu)
!
!_______________________________________________________________________
!  ECRITURE D'UN CHAMP - FORMAT MED - PROFIL POUR LES NOEUDS
!     -  -       -              -     -               -
!_______________________________________________________________________
!     ENTREES :
!       NOFIMD : NOM DU FICHIER MED
!       NCMPRF : NOMBRE DE COMPOSANTES DU CHAMP DE REFERENCE
!       NCMPVE : NOMBRE DE COMPOSANTES VALIDES EN ECRITURE
!       NUMCMP : NUMEROS DES COMPOSANTES VALIDES
!       EXICMP : EXISTENCE DES COMPOSANTES PAR MAILLES
!       NBVATO : NOMBRE DE VALEURS TOTALES
!       NBNOEC : NOMBRE D'ENTITES A ECRIRE (O, SI TOUTES)
!       LINOEC : LISTE DES ENTITES A ECRIRE SI EXTRAIT
!       ADSK, D, ... : ADRESSES DES TABLEAUX DES CHAMPS SIMPLIFIES
!       INNOCE : TABLEAU INDICATEUR DE NOEUD CENTRE
!                INNOCE(INO)=1 SI LE NOEUD INO EST UN NOEUD CENTRE
!                D'UNE MAILLE TRIA7,QUAD9,PENTA18 OU HEXA27
!
!     SORTIES :
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
!         CAIMPK : CARACTERES POUR CHAQUE IMPRESSION
!                  CAIMPK(1) = NOM DE LA LOCALISATION ASSOCIEE
!                  CAIMPK(2) = NOM DU PROFIL AU SENS MED
!                  CAIMPK(3) = NOM DE L'ELEMENT DE STRUCTURE
!       PROFAS : PROFIL ASTER. C'EST LA LISTE DES NUMEROS ASTER DES
!                NOEUDS POUR LESQUELS LE CHAMP EST DEFINI
!
!  COMMENTAIRE : C'EST AVEC L'USAGE DE ZL(ADSL) QU'IL FAUT FILTRER LES
!                COMPOSANTES. SEUL MOYEN FIABLE AVEC UN CHAMELEM
!                SIMPLIFIE.
!_______________________________________________________________________
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/infniv.h"
#include "asterfort/ircmpf.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: nbvato, ncmprf, ncmpve
    integer(kind=8) :: numcmp(ncmprf), innoce(nbvato)
    integer(kind=8) :: nbnoec
    integer(kind=8) :: linoec(*)
    integer(kind=8) :: adsl
    integer(kind=8) :: caimpi(10)
    integer(kind=8) :: profas(nbvato)
!
    character(len=*) :: nofimd
    character(len=80) :: caimpk(3)
!
    aster_logical :: exicmp(nbvato)
    character(len=8) :: nosdfu
!
! 0.2. ==> COMMUNS
!
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='IRCMPN')
!
    character(len=80) :: ednopf
    parameter(ednopf=' ')
    character(len=80) :: ednoga
    parameter(ednoga=' ')
!                         12345678901234567890123456789012
!
    integer(kind=8) :: typnoe
    parameter(typnoe=0)
    integer(kind=8) :: ednopg
    parameter(ednopg=1)
!
    character(len=64) :: noprof
!
    integer(kind=8) :: ifm, nivinf
!
    integer(kind=8) :: iaux, jaux
    integer(kind=8) :: nrcmp, nbtot, nbloc, start
    integer(kind=8) :: nval, rang, nbproc, jnbno, jno, nbnov, nbnoect(1)
    mpi_int :: mrank, msize
    aster_logical :: lficUniq, lnbnol, lnoec
    real(kind=8) :: start_time, end_time
!
!====
! 1. PREALABLES
!====
!
    call infniv(ifm, nivinf)
!
    if (nivinf .gt. 1) then
        call cpu_time(start_time)
        write (ifm, 1001) 'DEBUT DE '//nompro
    end if
1001 format(/, 4x, 10('='), a, 10('='),/)
!
!====
! 2. ON REMPLIT UN PREMIER TABLEAU PAR NOEUD :
!    VRAI DES QU'UNE DES COMPOSANTES DU CHAMP EST PRESENTE SUR LE NOEUD
!    FAUX SINON
!    REMARQUES: 1- ON EXAMINE LES NCMPVE COMPOSANTES QUI SONT DEMANDEES,
!    MAIS IL FAUT BIEN TENIR COMPTE DE NCMPRF, NOMBRE DE COMPOSANTES DE
!    REFERENCE, POUR L'ADRESSAGE DANS LE TABLEAU ADSL
!               2- SI LE NOEUD EST UN NOEUD CENTRE, ON L'OUBLIE
!====
!
    do iaux = 0, nbvato-1
!
        if (innoce(iaux+1) .eq. 1) then
            exicmp(iaux+1) = .false.
            goto 21
        end if
!
        jaux = adsl-1+iaux*ncmprf
        do nrcmp = 1, ncmpve
            if (zl(jaux+numcmp(nrcmp))) then
                exicmp(iaux+1) = .true.
                goto 21
            end if
        end do
!
21  end do
!
!====
! 3. PROFAS : LISTE DES NOEUDS POUR LESQUELS ON AURA IMPRESSION
!    UN NOEUD EN FAIT PARTIE SI ET SEULEMENT SI AU MOINS UNE COMPOSANTE
!    Y EST DEFINIE ET S'IL FAIT PARTIE DU FILTRAGE DEMANDE
!====
!
    nval = 0
    lficUniq = .false._1
    if (nosdfu .ne. ' ') then
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
!
        call jeveuo(nosdfu//'.NBNO', 'L', jnbno)
        nbnov = zi(jnbno+1)
        call jeveuo(nosdfu//'.NOEU', 'L', jno)
        lficUniq = .true._1
    else
        nbnov = nbvato
    end if
!
! 3.1. ==> SANS FILTRAGE : C'EST LA LISTE DES NOEUDS AVEC UNE COMPOSANTE
!          VALIDE
!
    if (lficUniq) then
        nbnoect(1) = nbnoec
        call asmpi_comm_vect('MPI_SUM', 'I', nbval=1, vi=nbnoect)
        lnoec = nbnoect(1) .eq. 0
    else
        lnoec = nbnoec .eq. 0
    end if
    if (lnoec) then
!
        do iaux = 1, nbvato
            if (exicmp(iaux)) then
                if (lficUniq) then
                    if (zi(jno+iaux-1) .gt. 0) then
                        nval = nval+1
                        profas(nval) = iaux
                    end if
                else
                    nval = nval+1
                    profas(nval) = iaux
                end if
            end if
        end do
!
! 3.2. ==> AVEC FILTRAGE
!
    else
!
        do jaux = 1, nbnoec
            iaux = linoec(jaux)
            if (exicmp(iaux)) then
                if (lficUniq) then
                    if (zi(jno+iaux-1) .gt. 0) then
                        nval = nval+1
                        profas(nval) = iaux
                    end if
                else
                    nval = nval+1
                    profas(nval) = iaux
                end if
            end if
        end do
!
    end if
!
!====
! 4. CARACTERISATIONS DES IMPRESSIONS
!====
!
!                  CAIMPI(1,I) = TYPE D'EF / MAILLE ASTER (0, SI NOEUD)
    caimpi(1) = 0
!                  CAIMPI(2,I) = NOMBRE DE POINTS DE GAUSS
    caimpi(2) = ednopg
!                  CAIMPI(3,I) = NOMBRE DE SOUS-POINTS
    caimpi(3) = ednopg
!                  CAIMPI(4,I) = NOMBRE DE MAILLES A ECRIRE
    caimpi(4) = ednopg
!                  CAIMPI(4,I) = NOMBRE DE MAILLES A ECRIRE
    caimpi(5) = ednopg
!                  CAIMPI(4,I) = NOMBRE DE MAILLES A ECRIRE
    caimpi(6) = ednopg
!                  CAIMPI(4,I) = NOMBRE DE MAILLES A ECRIRE
    caimpi(7) = nval
!                  CAIMPI(5,I) = TYPE DE MAILLES ASTER (0, SI NOEUD)
    caimpi(8) = 0
!                  CAIMPI(6,I) = TYPE GEOMETRIQUE AU SENS MED
    caimpi(9) = typnoe
!                  CAIMPI(7,I) = NOMBRE DE MAILLES IDENTIQUES
    caimpi(10) = nbvato
!                  CAIMPK(1) = NOM DE LA LOCALISATION ASSOCIEE
    caimpk(1) = ednoga
!                  CAIMPK(2) = NOM DU PROFIL AU SENS MED
    caimpk(2) = ednopf
!                  CAIMPK(3) = NOM DE L'ELEMENT DE STRUCTURE AU SENS MED
    caimpk(3) = ednopf
!
!GN      WRITE(IFM,*) 'A LA FIN DE 4, CAIMPI :'
!GN      WRITE(IFM,4000) (CAIMPI(IAUX,1),IAUX = 1 , 7)
!GN 4000 FORMAT('TYPN =',I4,', NB PG =',I4,', NB SPT =',I4,
!GN     >       ', NBVAL ECR =',I8,', MA ASTER =',I4,
!GN     >       ', TYPE GEO MED =',I4,', NBVAL TOT =',I8)
!
    if (nivinf .gt. 1) then
        write (ifm, 3301) nompro, ' : NOMBRE TOTAL DE VALEURS    : ',&
     &                   nbvato
        write (ifm, 3301) nompro, ' : NOMBRE DE VALEURS A ECRIRE : ', &
            nval
    end if
3301 format(4x, a6, a, i8)
!
!====
! 5. STOCKAGE DU PROFIL DANS LE FICHIER MED
!    REMARQUE : DANS LE CAS DES NOEUDS, IL Y A IDENTITE ENTRE LES
!               NUMEROTATIONS ASTER ET MED DES NOEUDS (CF IRMMNO)
!====
!
    iaux = 0
    lnbnol = .false._1
    if (lficUniq) then
        if (nval .ne. nbnov .and. nval .ne. 0) then
            iaux = 1
        end if
        call asmpi_comm_vect('MPI_MAX', 'I', 1, 0, sci=iaux)
        if (iaux .eq. 1) lnbnol = .true._1
    else
        if (nval .ne. nbnov .and. nval .ne. 0) lnbnol = .true._1
    end if
    if (lnbnol) then
!
        call ircmpf(nofimd, nval, profas, noprof, nosdfu, &
                    0, 0, nbtot, nbloc, start)
!
        caimpk(2) = noprof
!
    end if
!
!====
! 6. LA FIN
!====
!
    if (nivinf .gt. 1) then
        write (ifm, 1001) 'FIN DE '//nompro
        call cpu_time(end_time)
        write (ifm, *) "=========== EN ", end_time-start_time, "sec"
    end if
!
end subroutine
