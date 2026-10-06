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
subroutine ircmec(idfimd, nochmd, nomprf, nolopg, numpt, &
                  instan, numord, val, ncmpve, nbenty, &
                  nbrepg, nvalec, typent, typgeo, nosdfu, &
                  tymast, nrimpr, nosdpr, codret)
!_______________________________________________________________________
!     ECRITURE D'UN CHAMP -  FORMAT MED - ECRITURE
!        -  -       -               -     --
!_______________________________________________________________________
!     ENTREES :
!       IDFIMD : IDENTIFIANT DU FICHIER MED
!       NOMAM2 : NOM DU MAILLAGE MED
!       NOCHMD : NOM MED DU CHAMP A ECRIRE
!       NOMPRF : NOM MED DU PROFIL ASSOCIE AU CHAMP
!       NOLOPG : NOM MED LOCALISATION DES PTS DE GAUSS ASSOCIEE AU CHAMP
!       NUMPT  : NUMERO DE PAS DE TEMPS
!       INSTAN : VALEUR DE L'INSTANT A ARCHIVER
!       NUMORD : NUMERO D'ORDRE DU CHAMP
!       VAL    : VALEURS EN MODE ENTRELACE
!       NCMPVE : NOMBRE DE COMPOSANTES VALIDES EN ECRITURE
!       NBENTY : NOMBRE D'ENTITES DU TYPE CONSIDERE
!       NBREPG : NOMBRE DE POINTS DE GAUSS
!       NVALEC : NOMBRE DE VALEURS A ECRIRE EFFECTIVEMENT
!       TYPENT : TYPE D'ENTITE MED DU CHAMP A ECRIRE
!       TYPGEO : TYPE GEOMETRIQUE MED DU CHAMP A ECRIRE
!     SORTIES:
!       CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!_______________________________________________________________________
!
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "jeveux.h"
#include "asterfort/as_mfdraw.h"
#include "asterfort/as_mfdrpw.h"
#include "asterfort/as_mfrall.h"
#include "asterfort/as_mfrblc.h"
#include "asterfort/as_mfrdea.h"
#include "asterfort/infniv.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    character(len=*) :: nochmd, nomprf, nolopg, nosdpr
    character(len=8) :: nosdfu
!
    med_idt :: idfimd
    integer(kind=8) :: numpt, numord
    integer(kind=8) :: ncmpve, nbenty, nbrepg, nvalec
    integer(kind=8) :: typent, typgeo, tymast, nrimpr
!
    real(kind=8) :: instan
    real(kind=8) :: val(*)
!
    integer(kind=8) :: codret
!
! 0.2. ==> COMMUNS
!
!
! 0.3. ==> VARIABLES LOCALES
!
!
    character(len=6) :: nompro
    parameter(nompro='IRCMEC')
!
    character(len=32) :: ednopf
!                         12345678901234567890123456789012
    parameter(ednopf='                                ')
    integer(kind=8) :: edfuin
    parameter(edfuin=0)
    integer(kind=8) :: edall
    parameter(edall=0)
    integer(kind=8) :: ednopt
    parameter(ednopt=-1)
    integer(kind=8) :: ednopg
    parameter(ednopg=1)
    integer(kind=8) :: edcomp
    parameter(edcomp=2)
    character(len=32) :: ednoga
    parameter(ednoga='                                ')
!
    integer(kind=8) :: ifm, nivinf
    integer(kind=8) :: iaux, jnbno, nbentl, nbentt, start, filter(1)
    integer(kind=8) :: jnbma, nbbloc
    integer(kind=8), pointer :: v_nb_ma(:) => null()
    real(kind=8) :: start_time, end_time
!
    character(len=8) :: saux08
    character(len=14) :: saux14
    character(len=35) :: saux35
!
!====
! 1. PREALABLES
!====
!
! 1.1. ==> RECUPERATION DU NIVEAU D'IMPRESSION
!
    call infniv(ifm, nivinf)
!
! 1.2. ==> INFORMATION
!
    if (nivinf .gt. 1) then
        call cpu_time(start_time)
        write (ifm, 1001) 'DEBUT DE '//nompro
1001    format(/, 4x, 10('='), a, 10('='),/)
        call utmess('I', 'MED_49')
        write (ifm, 13001) nbrepg, typent, typgeo
        do 13, iaux = 1, ncmpve
            write (ifm, 13002)&
         &    '. PREMIERE ET DERNIERE VALEURS A ECRIRE POUR LA COMPOSANTE',&
         &    iaux, ' : ', val(iaux), val((nvalec*nbrepg-1)*ncmpve+iaux)
13          continue
            end if
13001       format(2x, '. NBREPG =', i4, ', TYPENT =', i4, ', TYPGEO =', i4)
13002       format(2x, a, i3, a3, 5g16.6)
!
!====
! 2. ECRITURE DES VALEURS
!    LE TABLEAU DE VALEURS EST UTILISE AINSI :
!        TV(NCMPVE,NBSP,NBPG,NVALEC)
!    TV(1,1,1,1), TV(2,1,1,1), ..., TV(NCMPVE,1,1,1),
!    TV(1,2,1,1), TV(2,2,1,1), ..., TV(NCMPVE,2,1,1),
!            ...     ...     ...
!    TV(1,NBSP,NBPG,NVALEC), TV(2,NBSP,NBPG,NVALEC), ... ,
!                                      TV(NCMPVE,NBSP,NBPG,NVALEC)
!    C'EST CE QUE MED APPELLE LE MODE ENTRELACE
!    REMARQUE : LE 6-EME ARGUMENT DE as_mfdrpw EST LE NOMBRE DE VALEURS
!               C'EST LE PRODUIT DU NOMBRE TOTAL D'ENTITES DU TYPE EN
!               COURS PAR LE PRODUIT DES NOMBRES DE POINTS DE GAUSS
!               ET DE SOUS-POINT.
!               ATTENTION, CE N'EST DONC PAS LE NOMBRE DE VALEURS
!               REELLEMENT ECRITES MAIS PLUTOT LE NOMBRE MAXIMUM QU'ON
!               POURRAIT ECRIRE.
!====
!
! 2.1. ==> MESSAGES
!
!GN      PRINT *,'TABLEAU REELLEMENT ECRIT'
!GN      PRINT 1789,(VAL(IAUX),
!GN     >  IAUX=1,NVALEC*NBREPG*NCMPVE-1)
!GN 1789  FORMAT(10G12.5)
!
            if (nivinf .gt. 1) then
!                  12345678901235
                saux14 = '. ECRITURE DES'
!                  12345678901234567890123456789012345
                saux35 = ' VALEURS POUR LE NUMERO D''ORDRE : '
!
                if (nbrepg .eq. ednopg) then
                    write (ifm, 20001) saux14, ncmpve, nvalec, saux35, numord
                else
                    write (ifm, 20002) saux14, ncmpve, nbrepg, nvalec, saux35, &
                        numord
                end if
                if (numpt .ne. ednopt) then
                    write (ifm, 20003) numpt, instan
                end if
                if (nomprf .eq. ednopf) then
                    write (ifm, 20004)
                else
                    write (ifm, 20005) nomprf
                end if
                if (nolopg .eq. ednoga) then
                    write (ifm, 20006)
                else
                    write (ifm, 20007) nolopg
                end if
            end if
!
20001       format(2x, a14, i3, ' * ', i8, a35, i5)
20002       format(2x, a14, 2(i3, ' * '), i8, a35, i5)
20003       format(5x, '( PAS DE TEMPS NUMERO :', i5, ', T = ', g13.5, ' )')
20004       format(2x, '. PAS DE PROFIL')
20005       format(2x, '. NOM DU PROFIL : ', a)
20006       format(2x, '. PAS DE LOCALISATION DE POINTS DE GAUSS')
20007       format(2x, '. NOM DE LA LOCALISATION DES POINTS DE GAUSS : ', a)
!
! 2.2. ==> NOMBRE DE VALEURS
!
            iaux = nbenty
!
! 2.3. ==> ECRITURE VRAIE
!
            if (nosdfu .ne. ' ') then
!
                if (typent .eq. 0 .or. typent .eq. 4) then
                    if (nomprf .eq. ' ') then
                        call jeveuo(nosdfu//'.MATY', 'L', jnbma)
                        start = zi(jnbma+3*(tymast-1))
                        nbentl = zi(jnbma+3*(tymast-1)+1)
                        nbentt = zi(jnbma+3*(tymast-1)+2)
                    else
                        call jeveuo(nosdpr, 'L', vi=v_nb_ma)
                        nbentt = v_nb_ma(3*nrimpr-2)
                        nbentl = v_nb_ma(3*nrimpr-1)
                        start = v_nb_ma(3*nrimpr)
                    end if
                else
                    if (nomprf .eq. ' ') then
                        call jeveuo(nosdfu//'.NBNO', 'L', jnbno)
                    else
                        call jeveuo(nosdfu//'.NBNOP', 'L', jnbno)
                    end if
                    start = zi(jnbno)
                    nbentl = zi(jnbno+1)
                    nbentt = zi(jnbno+2)
                end if
                nbbloc = 1
                if (nbentl .eq. 0) nbbloc = 0
                call as_mfrall(1, filter, codret)
                call as_mfrblc(idfimd, nbentt, nbrepg, ncmpve, 0, &
                               edfuin, 2, nomprf, start, nbentl, &
                               nbbloc, nbentl, 0, filter(1), codret)
!
                if (codret .ne. 0) then
                    saux08 = 'mfrblc'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
!
                call as_mfdraw(idfimd, nochmd, filter(1), val, nolopg, &
                               typent, typgeo, numpt, instan, numord, &
                               codret)
!
                if (codret .ne. 0) then
                    saux08 = 'mfdraw'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
!
                call as_mfrdea(1, filter, codret)
                if (codret .ne. 0) then
                    saux08 = 'mfrdea'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
            else
!
                call as_mfdrpw(idfimd, nochmd, val, edfuin, iaux, &
                               nolopg, edall, nomprf, edcomp, typent, &
                               typgeo, numpt, instan, numord, codret)
!
                if (codret .ne. 0) then
                    saux08 = 'mfdrpw'
                    call utmess('F', 'DVP_97', sk=saux08, si=codret)
                end if
!
            end if
!
            if (nivinf .gt. 1) then
                call cpu_time(end_time)
                write (ifm, *) '    ==========FIN DE '//nompro//' EN ', &
                    end_time-start_time, " sec=========="
            end if
!
            end subroutine
