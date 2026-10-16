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

subroutine ircmpf(nofimd, nvalty, profil, noprof, nosdfu, &
                  typent, typmai, nbtot, nbloc, start)
!
!_______________________________________________________________________
!  ECRITURE D'UN CHAMP - FORMAT MED - PROFIL
!     -  -       -              -     -  -
!_______________________________________________________________________
!     ENTREES :
!       NOFIMD : NOM DU FICHIER MED
!       NVALTY : NOMBRE DE VALEURS DU TYPE
!       PROFIL : PROFIL ENTIER
!     SORTIES :
!       NOPROF : NOM DU PROFIL AU SENS MED
!_______________________________________________________________________
!
    use as_med_module, only: as_med_open
    implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "jeveux.h"
#include "asterfort/as_mficlo.h"
#include "asterfort/as_mpfnpf.h"
#include "asterfort/as_mpfpfi.h"
#include "asterfort/as_mpfprr.h"
#include "asterfort/as_mpfprw.h"
#include "asterfort/asmpi_barrier.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/infniv.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterc/asmpi_bcast_char80.h"
#include "asterc/asmpi_comm.h"
    integer(kind=8) :: nvalty, profil(nvalty), typent, typmai, nbtot, nbloc, start
!
    character(len=*) :: nofimd
    character(len=*) :: noprof
    character(len=8) :: nosdfu
!
! 0.2. ==> COMMUNS
!
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='IRCMPF')
    integer(kind=8) :: edlect
    parameter(edlect=0)
    integer(kind=8) :: edleaj
    parameter(edleaj=1)
!
    integer(kind=8) :: ifm, nivinf
!
    med_idt :: idfimd
    integer(kind=8) :: nbprof, lgprof, adprof, adnopf, nrprty, jprof
    integer(kind=8) :: iaux, jaux, jent, jproc
    integer(kind=8) :: codret, decal, nvalty2, jdecal
    integer(kind=8), pointer :: v_decal(:) => null()
!
    character(len=8) :: saux08
    character(len=24) :: ntprof, ntnopf
    character(len=64) :: nopr64
    character(len=80) :: comprf(1)
    integer(kind=8) :: rang, nbproc
    mpi_int :: mrank, msize, world, proc, taille
    real(kind=8) :: start_time, end_time
!
!====
! 1. PREALABLES
!====
! 1.1. ==> RECUPERATION DU NIVEAU D'IMPRESSION
!
    call infniv(ifm, nivinf)
!
    if (nivinf .gt. 1) then
        call cpu_time(start_time)
        write (ifm, 1001) 'DEBUT DE '//nompro
    end if
1001 format(/, 4x, 10('='), a, 10('='),/)
!
! 1.2. ==> NOMS DES TABLEAUX
!               12   345678   9012345678901234
    ntprof = '&&'//nompro//'.PROFIL_MED_LU  '
    ntnopf = '&&'//nompro//'.NOM_PROFIL_MED '
!
    nvalty2 = 0
    decal = 0
    if (nosdfu .ne. ' ') then
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
        call jedetr('&&IRCMPF.PROCS')
        call wkvect('&&IRCMPF.PROCS', 'V V I', nbproc+1, jproc)
        do iaux = rang+1, nbproc
            zi(jproc+iaux) = nvalty
        end do
        call asmpi_comm_vect('MPI_SUM', 'I', nbval=nbproc+1, vi=zi(jproc))
        if (typent .eq. 0) then
            call jeveuo(nosdfu//'.NOEU', 'L', vi=v_decal)
        else
            call jeveuo(nosdfu//'.MAIL', 'L', vi=v_decal)
        end if
        nvalty2 = zi(jproc+nbproc)
        call jedetr('&&IRCMPF.PROFIL')
        call wkvect('&&IRCMPF.PROFIL', 'V V I', nvalty2, jprof)
        decal = zi(jproc+rang)
        do iaux = 1, nvalty
            zi(jprof+iaux+decal-1) = v_decal(profil(iaux))
        end do
        call asmpi_comm_vect('MPI_SUM', 'I', nbval=nvalty2, vi=zi(jprof))
        if (typent .eq. 0) then
            call jeveuo(nosdfu//'.NBNOP', 'E', jent)
            zi(jent) = decal+1
            zi(jent+1) = nvalty
            zi(jent+2) = nvalty2
        end if
        if (rang .ne. 0) goto 500
    else
        nvalty2 = nvalty
    end if
!
! 1.3. ==> A PRIORI, PAS DE PROFIL DANS LE FICHIER
!
    nrprty = 0
!
!====
! 2. ON OUVRE LE FICHIER EN LECTURE
!====
!
    call as_med_open(idfimd, nofimd, edlect, codret)
    if (codret .ne. 0) then
        saux08 = 'mfiope'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
!====
! 3. REPERAGE DU NOMBRE DE PROFILS DEJA ENREGISTRES
!    S'IL Y EN A, ON ALLOUE UN TABLEAU POUR STOCKER LEURS NOMS
!====
!
    nbprof = 0
    call as_mpfnpf(idfimd, nbprof, codret)
    if (codret .ne. 0) then
        saux08 = 'mpfnpf'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
    if (nivinf .gt. 1) then
        write (ifm, *) '   NOMBRE DE PROFILS DANS LE FICHIER : ', nbprof
    end if
!
    if (nbprof .ne. 0) then
        call wkvect(ntnopf, 'V V K80', nbprof, adnopf)
    end if
!
!====
! 4. LECTURE DE CHACUN DES PROFILS ET COMPARAISON AVEC CELUI RETENU
!====
!
    do iaux = 1, nbprof
!
! 4.1. ==> NOM ET NOMBRE DE VALEURS DU IAUX-EME PROFIL
!
        call as_mpfpfi(idfimd, iaux, nopr64, lgprof, codret)
        noprof = nopr64
        if (codret .ne. 0) then
            saux08 = 'mpfpfi'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
        if (nivinf .gt. 1) then
            write (ifm, 4101) iaux, noprof, lgprof
        end if
4101    format(5x, 'LECTURE DU PROFIL NUMERO', i8, &
               /, 5x, '...NOM       : ', a, &
               /, 5x, '... LONGUEUR : ', i8)
!
        zk80(adnopf+iaux-1) = noprof
!
! 4.2. ==> SI LA LONGUEUR EST LA MEME ET QU'ON N'A TOUJOURS PAS TROUVE
!          UN PPOFIL IDENTIQUE, ON LIT LE PROFIL COURANT ET ON COMPARE
!
        if (nvalty2 .eq. lgprof .and. nrprty .eq. 0) then
!
! 4.2.1. ==> LECTURE DES VALEURS DU PROFIL
!
            call wkvect(ntprof, 'V V I', lgprof, adprof)
!
            call as_mpfprr(idfimd, zi(adprof), lgprof, noprof, codret)
            if (codret .ne. 0) then
                saux08 = 'mpfprr'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            end if
!
            if (nivinf .gt. 1) then
                write (ifm, 4201) zi(adprof), zi(adprof+lgprof-1)
            end if
4201        format(5x, '... 1ERE ET DERNIERE VALEURS : ', 2i8)
!
! 4.2.2. ==> ON COMPARE TERME A TERME.
!            DES QU'UNE VALEUR DIFFERE, ON PASSE AU PROFIL SUIVANT.
!            SI TOUS LES TERMES SONT EGAUX, C'EST LE BON PROFIL !
!            ON PEUT SORTIR DE LA RECHERCHE DES PROFILS.
!
            if (nosdfu .eq. ' ') then
                do jaux = 1, lgprof
                    if (profil(jaux) .ne. zi(adprof+jaux-1)) then
                        goto 423
                    end if
                end do
            else
                do jaux = 1, lgprof
                    if (zi(jprof+jaux-1) .ne. zi(adprof+jaux-1)) then
                        goto 423
                    end if
                end do
            end if
!
            nrprty = iaux
!
            if (nivinf .gt. 1) then
                write (ifm, 4202)
            end if
4202        format('...... CE PROFIL EST IDENTIQUE A CELUI VOULU')
!
            goto 51
!
! 4.2.3. ==> MENAGE AVANT D'EXAMINER UN NOUVEAU PROFIL
!
423         continue
!
            call jedetr(ntprof)
!
        end if
!
    end do
!
!====
! 5. FERMETURE DU FICHIER
!====
!
51  continue
!
    call as_mficlo(idfimd, codret)
    if (codret .ne. 0) then
        saux08 = 'mficlo'
        call utmess('F', 'DVP_97', sk=saux08, si=codret)
    end if
!
!====
! 6. SI AUCUN PROFIL N'A ETE TROUVE, ON ECRIT LE NOTRE DANS LE FICHIER
!====
!
    if (nrprty .eq. 0) then
!
! 6.1. ==> OUVERTURE FICHIER MED EN MODE 'LECTURE_AJOUT'
!    CELA SIGNIFIE QUE LE FICHIER EST ENRICHI MAIS ON NE PEUT PAS
!    ECRASER UNE DONNEE EXISTANTE
!
        call as_med_open(idfimd, nofimd, edleaj, codret)
        if (codret .ne. 0) then
            saux08 = 'mfiope'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
! 6.2. ==> ELABORATION D'UN NOM DE PROFIL
!          IL SERA DU TYPE 'PROFIL_00000000N'
!          POUR ETRE SUR DE FINIR PAR TROUVER UN NOM QUI N'A PAS SERVI,
!          IL SUFFIT D'ESSAYER PLUS DE NUMEROS QUE DE PROFILS DEJA
!          ENREGISTRES. QUELLE RUSE DIABOLIQUE !
!
        noprof(17:64) = '                                                '
!
        do iaux = 1, nbprof+1
!
            call codent(iaux, 'D0', saux08)
!                         12345678
            noprof(1:16) = 'PROFIL__'//saux08
!
            do jaux = 0, nbprof-1
                if (noprof .eq. zk80(adnopf+jaux) (1:64)) then
                    goto 62
                end if
            end do
            goto 622
62      end do
!
622     continue
!
! 6.3. ==> ECRITURE DU PROFIL
!
        if (nivinf .gt. 1) then
            write (ifm, 6301) noprof, nvalty2, profil(1), profil(nvalty2)
6301        format(4x, 'PROFIL A CREER :', &
                   /, 4x, '. NOM                      = ', a, &
                   /, 4x, '. LONGUEUR                 = ', i8, &
                   /, 4x, '. 1ERE ET DERNIERE VALEURS = ', 2i8)
        end if
        if (nosdfu .eq. ' ') then
            call as_mpfprw(idfimd, profil, nvalty2, noprof, codret)
        else
            call as_mpfprw(idfimd, zi(jprof), nvalty2, noprof, codret)
        end if
        if (codret .ne. 0) then
            saux08 = 'mpfprw'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
! 6.4. ==> FERMETURE FICHIER MED
!
        call as_mficlo(idfimd, codret)
        if (codret .ne. 0) then
            saux08 = 'mficlo'
            call utmess('F', 'DVP_97', sk=saux08, si=codret)
        end if
!
    end if

500 continue
    if (nosdfu .ne. ' ') then
        call asmpi_comm('GET', world)
        comprf(1) = noprof
        taille = 1
        proc = 0
        call asmpi_bcast_char80(comprf, taille, proc, world)
        noprof = comprf(1)
        call asmpi_barrier()
    end if
!
!====
! 7. LA FIN
!====
!
    call jedetr(ntnopf)
    call jedetr(ntprof)
    nbtot = nvalty2
    nbloc = nvalty
    start = decal+1
!
    if (nivinf .gt. 1) then
        write (ifm, 1001) 'FIN DE '//nompro
        call cpu_time(end_time)
        print *, "=========== EN ", end_time-start_time, "sec"
    end if
!
end subroutine
