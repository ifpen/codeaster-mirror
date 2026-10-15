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

subroutine amumpi(option, lquali, ldist, kxmps, type, lmhpc, lbloc)
!
!
    implicit none
!--------------------------------------------------------------
! BUT : ROUTINE DE PARAMETRAGE MUMPS POUR AMUMPS/C/D/Z
!
! IN  OPTION:   IN   : OPTION D'UTILISATION.
! IN  LQUALI:  LOG   : LOGICAL EN CAS DE CRITERE DE QUALITE
! IN  LDIST :  LOG   : LOGICAL MUMPS DISTRIBUE OR NOT
! IN  KXMPS :   IN   : INDICE DE L'INSTANCE MUMPS DANS DMPS
! IN  TYPE  :   K1   : TYPE DU POINTEUR R OU C
! IN  LBLOC :  LOG   : LOGIQUE PRECISANT SI ON EFFECTUE L ANALYSE PAR BLOCS
!---------------------------------------------------------------
!
#include "asterc/asmpi_comm.h"
#include "asterc/r4maem.h"
#include "asterf_config.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeexin.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "threading_interfaces.h"
    integer(kind=4) :: option
    integer(kind=8) :: kxmps
    aster_logical :: lquali, ldist, lmhpc, lbloc
    character(len=1) :: type
!
#ifdef ASTER_HAVE_MUMPS
#include "asterf_mumps.h"
    mpi_int :: mpicou, mpimum
    integer(kind=8) :: nicntl, ncntl
    parameter(nicntl=60, ncntl=15)
    type(smumps_struc), pointer :: smpsk => null()
    type(cmumps_struc), pointer :: cmpsk => null()
    type(dmumps_struc), pointer :: dmpsk => null()
    type(zmumps_struc), pointer :: zmpsk => null()
    integer(kind=8) :: ifm, niv, i, isymm, isymv, isym, nbproc, i1, i2
    integer(kind=8) :: nprec, nbomp, nbrhs, iret1, iret2, iret3
    integer(kind=8) :: k371, k268, iaux, redmpi
    mumps_int :: i4, icntl(nicntl)
    real(kind=8) :: cntl(ncntl), rr4max, blreps
    aster_logical :: lverbose
    character(len=4) :: typm, etam
    character(len=8) :: kacmum
    character(len=14) :: nonu
    character(len=19) :: nomat, nosolv
    character(len=24) :: kmumps_sparse(3)
    character(len=24), pointer :: refa(:) => null()
    real(kind=8), pointer :: slvr(:) => null()
    character(len=24), pointer :: slvk(:) => null()
    integer(kind=8), pointer :: slvi(:) => null()
    call jemarq()
! --- COMMUNICATEUR MPI DE TRAVAIL
    call asmpi_comm('GET', mpicou)
    call infniv(ifm, niv)
! pour forcer le mode VERBOSE de MUMPS
    lverbose = .false.
!
!       ------------------------------------------------
!        INITS
!       ------------------------------------------------
    rr4max = r4maem()
!
    if (type .eq. 'S') then
        smpsk => smps(kxmps)
    else if (type .eq. 'C') then
        cmpsk => cmps(kxmps)
    else if (type .eq. 'D') then
        dmpsk => dmps(kxmps)
    else if (type .eq. 'Z') then
        zmpsk => zmps(kxmps)
    else
        ASSERT(.false.)
    end if
!
    nomat = nomats(kxmps)
    nosolv = nosols(kxmps)
    nonu = nonus(kxmps)
    etam = etams(kxmps)
    call jeveuo(nomat//'.REFA', 'L', vk24=refa)
    call jeveuo(nosolv//'.SLVK', 'L', vk24=slvk)
    call jeveuo(nosolv//'.SLVI', 'L', vi=slvi)
    nprec = slvi(1)
    call jeveuo(nosolv//'.SLVR', 'L', vr=slvr)

    kacmum = 'XXXXXXXX'
    kacmum = trim(adjustl(slvk(5)))
    blreps = slvr(4)
    redmpi = slvi(7)
    nbrhs = slvi(8)
!
!       -----------------------------------------------------
!        INITIALISATION SYM, PAR ET JOB POUR MUMPS (CREATION)
!       -----------------------------------------------------
    if (option .eq. 0) then
!
        if (type .eq. 'S') then
            smpsk%comm = mpicou
        else if (type .eq. 'C') then
            cmpsk%comm = mpicou
        else if (type .eq. 'D') then
            dmpsk%comm = mpicou
        else if (type .eq. 'Z') then
            zmpsk%comm = mpicou
        else
            ASSERT(.false.)
        end if
!
! ---     ISYM = 0 => NON-SYMETRIQUE
! ---     ISYM = 1 => SYMETRIQUE DEFINIE POSITIVE
! ---     ISYM = 2 => SYMETRIQUE  GENERAL
! ---     ISYMM DEDUIT DE LA MATRICE : NONSYM OU SYMGEN
        typm = refa(9) (1:4)
        if (typm .eq. 'MR') then
            isymm = 0
        else if (typm .eq. 'MS') then
            isymm = 2
        else
            ASSERT(.false.)
        end if
!
! ---     PRISE EN COMPTE DE LA VOLONTE DE L'UTILISATEUR
! ---     => ISYMV
        if (slvk(3) .eq. 'NONSYM') then
            isymv = 0
        else if (slvk(3) .eq. 'SYMDEF') then
            isymv = 1
        else if (slvk(3) .eq. 'SYMGEN') then
            isymv = 2
        else if (slvk(3) .eq. 'AUTO') then
            isymv = -1
        else
            ASSERT(.false.)
        end if
!
! ---     STRATEGIE PRUDENTE ET CONSERVATIVE
! ---     SI AUTO: NONSYM OU SYMGEN SUIVANT LA STRUCTURE DE LA MATRICE
! ---     SINON, ON APPLIQUE LE CHOIX DE L'UTILISATEUR
        if (isymv .eq. -1) then
            isym = isymm
        else if (isymv .eq. 0) then
            isym = isymv
        else
            if (isymm .eq. 0) then
                call utmess('F', 'FACTOR_56', sk=slvk(3))
            else
                isym = isymv
            end if
        end if
!
! ---     PARAMETRES D'INITIALISATION DE L'OCCURENCE MUMPS KXMPS
        i4 = to_mumps_int(isym)
        if (type .eq. 'S') then
            smpsk%sym = i4
            smpsk%par = 1
            smpsk%job = -1
        else if (type .eq. 'C') then
            cmpsk%sym = i4
            cmpsk%par = 1
            cmpsk%job = -1
        else if (type .eq. 'D') then
            dmpsk%sym = i4
            dmpsk%par = 1
            dmpsk%job = -1
        else if (type .eq. 'Z') then
            zmpsk%sym = i4
            zmpsk%par = 1
            zmpsk%job = -1
        else
            ASSERT(.false.)
        end if
!
!       ------------------------------------------------------
!        INITIALISATION ICNTL/CNTL POUR MUMPS (ANALYSE +FACTO)
!       ------------------------------------------------------
    else if (option .eq. 2) then
! ECRIRE LA MATRICE ET LE RHS VUS PAR MUMPS
!        dmpsk%WRITE_PROBLEM='chemin complet/nomfichier'
        if (type .eq. 'S') then
            nbproc = smpsk%nprocs
        else if (type .eq. 'C') then
            nbproc = cmpsk%nprocs
        else if (type .eq. 'D') then
            nbproc = dmpsk%nprocs
        else if (type .eq. 'Z') then
            nbproc = zmpsk%nprocs
        else
            ASSERT(.false.)
        end if
! SI 1 SEUL MPI PAS D'OPTION REDUCTION_MPI
        if (nbproc < 2) then
            redmpi = -9999
        end if
!
! ---     INIT
        do i = 1, nicntl
            icntl(i) = 0
        end do
        do i = 1, ncntl
            cntl(i) = 0.d0
        end do

!
! ---     OPTIONS AVANCEES (ACCELERATIONS)
! ------     TEST DE COMPATIBILITE ACCELERATION/VERSIONS
#ifndef ASTER_MUMPS_REDUCMPI
        if (redmpi > 1) then
! version non compatible avec option REDUCTION_MPI
            call utmess('I', 'FACTOR_49')
            redmpi = -9999
        end if
#endif

#ifdef ASTER_MUMPS_CONSORTIUM
        select case (kacmum)
        case ('FR', 'FR+', 'FR++', 'LR', 'LR+', 'LR++')
            !ok
        case ('AUTO')
            kacmum = 'FR+'
        case default
            ASSERT(.false.)
        end select
#else
        select case (kacmum)
        case ('FR', 'LR')
            !ok
        case ('FR+', 'FR++')
            kacmum = 'FR'
            call utmess('I', 'FACTOR_48', sk=kacmum)
        case ('LR+', 'LR++')
            kacmum = 'LR'
            call utmess('I', 'FACTOR_48', sk=kacmum)
        case ('AUTO')
            kacmum = 'FR'
        case default
            ASSERT(.false.)
        end select
#endif

! ------     API MUMPS ACCELERATION

        if (lbloc) then
            if ((kacmum .eq. 'FR') .or. (kacmum .eq. 'LR')) then
                ! parametrage par defaut de MUMPS
            else
                icntl(15) = 1
            end if
        else
            icntl(15) = 0
        end if
!
! ---    REDISTRIBUTION DE PARALLELISME: MPI_TO_OPENMP FEATURE
! ---    IMPACT SUR //ISME OPENMP L0 POTENTIEL
        icntl(17) = 0
        if (redmpi > 1) then
            i1 = nbproc/redmpi
            i2 = nbproc-(redmpi*i1)
            if (i2 .ne. 0) then
                call utmess('I', 'FACTOR_47')
            else
                icntl(17) = to_mumps_int(redmpi)
            end if
        end if
!
#ifdef ASTER_HAVE_OPENMP
        nbomp = omp_get_max_threads()
#else
        nbomp = 1
#endif
! POUR PRENDRE POTENTIEL IMPACT REDUCTION_MPI
        if (icntl(17) > 1) then
            nbomp = nbomp*icntl(17)
        end if
        k268 = 0
        k371 = 0
        icntl(35) = 0
        icntl(36) = 0
        icntl(37) = 0
        cntl(7) = -999.d0
        select case (kacmum)
        case ('FR')
! FR std
        case ('FR+', 'FR++')
! FR +  aggressive optimizations + MUMPS feature in advance if ++ (consortium version)
            k371 = 1
            if (kacmum(3:4) .eq. '++') then
                k268 = -2

#ifdef ASTER_MUMPS_CONSORTIUM
                icntl(58) = 2
#endif
                if (nbomp > 1) then
                    icntl(48) = 1
                end if
            end if
        case ('LR', 'LR+', 'LR++')
! BLR std, UFSC, ou BLR UCFS (si blreps<0)
            icntl(35) = 1
            if (blreps >= 0) then
                cntl(7) = blreps
            else
                icntl(36) = 1
                cntl(7) = -blreps
            end if
            if (kacmum(3:3) .eq. '+') then
! BLR+: priorite temps (consortium version)
#ifdef ASTER_MUMPS_CONSORTIUM
                icntl(58) = 2
#endif
                k268 = -2
                k371 = 1
                if (nbomp > 1) then
                    icntl(48) = 1
                end if
                if (kacmum(4:4) .eq. '+') then
! BLR+: priorite RAM (consortium version)
                    icntl(37) = 1
                    icntl(40) = 1
                end if
            end if
        case default
            ASSERT(.false.)
        end select

        if (type .eq. 'S') then
            smpsk%keep(268) = to_mumps_int(k268)
            smpsk%keep(371) = to_mumps_int(k371)
        else if (type .eq. 'C') then
            cmpsk%keep(268) = to_mumps_int(k268)
            cmpsk%keep(371) = to_mumps_int(k371)
        else if (type .eq. 'D') then
            dmpsk%keep(268) = to_mumps_int(k268)
            dmpsk%keep(371) = to_mumps_int(k371)
        else if (type .eq. 'Z') then
            zmpsk%keep(268) = to_mumps_int(k268)
            zmpsk%keep(371) = to_mumps_int(k371)
        end if
!
! ---     MESSAGES/ALERTES MUMPS
        icntl(1) = -1
        icntl(2) = -1
        icntl(3) = -1
        icntl(4) = 0
        if ((niv .ge. 2) .or. (lverbose)) then
            icntl(1) = to_mumps_int(ifm)
            icntl(2) = to_mumps_int(ifm)
            icntl(3) = to_mumps_int(ifm)
            icntl(4) = 2
        end if
! ---     FORMAT MATRICE
        icntl(5) = 0
! ---     PRETRAITEMENTS (SCALING/PERMUTATION)
        select case (slvk(2) (1:4))
        case ('SANS')
            icntl(6) = 0
            icntl(8) = 0
            icntl(12) = 1
        case ('AUTO')
            icntl(6) = 7
            icntl(8) = 77
            icntl(12) = 0
        case default
            ASSERT(.false.)
        end select
!
! ---     Ordering phase
! automatic choice of sequential or parallel analysis
        icntl(28) = 0
! automatic choice of the sequential ordering
        icntl(7) = 7
! automatic choice of the parallel ordering
        icntl(29) = 0

        if (slvk(4) .eq. 'AMD') then
            icntl(28) = 1
            icntl(7) = 0
        else if (slvk(4) .eq. 'AMF') then
            icntl(28) = 1
            icntl(7) = 2
        else if (slvk(4) .eq. 'SCOTCH') then
            icntl(28) = 1
            icntl(7) = 3
        else if (slvk(4) .eq. 'PTSCOTCH') then
#ifdef ASTER_HAVE_MPI
            icntl(28) = 2
            icntl(29) = 1
#else
            icntl(28) = 1
            icntl(7) = 3
            call utmess('I', 'FACTOR_89')
#endif
        else if (slvk(4) .eq. 'PORD') then
            icntl(28) = 1
            icntl(7) = 4
        else if (slvk(4) .eq. 'METIS') then
            icntl(28) = 1
            icntl(7) = 5
        else if (slvk(4) .eq. 'PARMETIS') then
#ifdef ASTER_HAVE_MPI
            icntl(28) = 2
            icntl(29) = 2
#else
            icntl(28) = 1
            icntl(7) = 5
            call utmess('I', 'FACTOR_90')
#endif
        else if (slvk(4) .eq. 'QAMD') then
            icntl(28) = 1
            icntl(7) = 6
        else if (slvk(4) .eq. 'AUTO') then
! choosen par default
        else if (slvk(4) .eq. 'SANS') then
        else
            ASSERT(.false.)
        end if
!
! ---     INITIALISATION EN DUR (EN DOUBLONS VS CALL DMUMPS JOB=-1)
! ---     MAIS ON NE SAIT JAMAIS AVEC LES EVOLUTIONS DES INITS DU PACKAG
! ---     ET CELA PERMET DE SURCHARGER PLUS RAPIDEMENT POUR TESTER
!
! ---     TYPE DE RESOLUTION: A OU AT
        icntl(9) = 1
!
! ---     RAFFINEMENT ITERATIF ET ANALYSE QUALITE SOLUTION
! ---     PARAMETRES ACTIVES JUSTE AVANT SOLVE VIA AMUMPI OPTION=3
        icntl(10) = 0
        cntl(2) = 0.d0
        icntl(11) = 0
!
! ---     PARALLELISME INDUIT PAR SCALAPACK (VOIR NPREC PLUS BAS)
        icntl(13) = 0
!
! ---     MEMOIRE SUPPL. POUR PIVOTAGE (DEFAUT:20)
        icntl(14) = to_mumps_int(slvi(2))
!
! ---     PAS UTILISE
        icntl(16) = 0
!
! --      DETECTION DE SINGULARITE/NOYAU
        icntl(25) = 0
        icntl(56) = 0
        if (nprec .ge. 0) then
            icntl(13) = 1
            icntl(24) = 1
            icntl(56) = 0
            cntl(3) = -10.d0**(-(nprec))
            cntl(5) = 1.d+20
            iaux = nprec-100
! -- AU CAS OU, POUR TESTER LE NOUVEAU PARAMETRAGE AVEC NPREC=108
! -- AU LIEU DE 8 PAR EXEMPLE
            if ((iaux .gt. 0) .and. (iaux .lt. 100)) then
                icntl(56) = 1
                cntl(3) = 10.d0**(-(iaux))
            end if
        else
            icntl(24) = 0
            cntl(3) = 0.d0
            cntl(5) = 0.d0
        end if
!
! ---     PIVOTAGE STATIQUE DESACTIVE
        cntl(4) = -1.d0
!
! ---     PARALLELISME/DISTRIBUTION SECOND MEMBRE/SOLUTION
        if (ldist .or. lmhpc) then
            icntl(18) = 3
        else
            icntl(18) = 0
        end if
        icntl(20) = 0
        icntl(21) = 0
!
! ---     GESTION MEMOIRE MUMPS
! ---     PARAMETRES ACTIVES APRES L'ANALYSE VIA AMUMPU OPTION=1
!
        icntl(22) = -999
        icntl(23) = -999
        if (type .eq. 'S') then
            smpsk%ooc_tmpdir = 'XXXX'
        else if (type .eq. 'C') then
            cmpsk%ooc_tmpdir = 'XXXX'
        else if (type .eq. 'D') then
            dmpsk%ooc_tmpdir = 'XXXX'
        else if (type .eq. 'Z') then
            zmpsk%ooc_tmpdir = 'XXXX'
        else
            ASSERT(.false.)
        end if
!
! ---     COMPLEMENT DE SCHUR
        icntl(19) = 0
        icntl(26) = 0
!
! ---     MULTIPLE RHS
        icntl(27) = -32
!
! ---     PAS DE CALCUL DE TERMES DE A-1
        icntl(30) = 0
!
! ---     ON GARDER LA FACTO EN MEMOIRE POUR LE SOLVE
        icntl(31) = 0
!
! ---     NON UTILISE
        icntl(32) = 0
!
! ---     PAS DE CALCUL DU DETERMINANT
        icntl(33) = 0
!
! ---   REMPLISSAGE DE DIFFERENTS OBJETS SUIVANT LE TYPE DU POINTEUR
! ---   DE MUMPS: DMUMPS_STRUC OU ZMUMPS_STRUC
        if (type .eq. 'S') then
            do i = 1, nicntl
                smpsk%icntl(i) = icntl(i)
            end do
            do i = 2, ncntl
                if (abs(cntl(i)) .gt. rr4max) then
                    ASSERT(.false.)
                end if
                smpsk%cntl(i) = real(cntl(i), kind=4)
            end do
        else if (type .eq. 'C') then
            do i = 1, nicntl
                cmpsk%icntl(i) = icntl(i)
            end do
            do i = 2, ncntl
                if (abs(cntl(i)) .gt. rr4max) then
                    ASSERT(.false.)
                end if
                cmpsk%cntl(i) = real(cntl(i), kind=4)
            end do
        else if (type .eq. 'D') then
            do i = 1, nicntl
                dmpsk%icntl(i) = icntl(i)
            end do
            do i = 2, ncntl
                dmpsk%cntl(i) = cntl(i)
            end do
        else if (type .eq. 'Z') then
            do i = 1, nicntl
                zmpsk%icntl(i) = icntl(i)
            end do
            do i = 2, ncntl
                zmpsk%cntl(i) = cntl(i)
            end do
        else
            ASSERT(.false.)
        end if
!
!       ------------------------------------------------------
!        INITIALISATION ICNTL/CNTL POUR MUMPS (SOLVE)
!       ------------------------------------------------------
    else if (option .eq. 3) then
! ---   POUR CMD ECLATEE RESOUDRE PRINCIPALEMENT
! ---   TEST DU COMMUNICATEUR COURANT AU CAS OU (ERREUR PROGRAMMEUR).
! ---   IL DOIT ETRE IDENTIQUE A CELUI PARAMETRE DS L'OCCURENCE MUMPS
        if (type .eq. 'S') then
            mpimum = smpsk%comm
        else if (type .eq. 'C') then
            mpimum = cmpsk%comm
        else if (type .eq. 'D') then
            mpimum = dmpsk%comm
        else if (type .eq. 'Z') then
            mpimum = zmpsk%comm
        else
            ASSERT(.false.)
        end if
        if (mpimum .ne. mpicou) then
            ASSERT(.false.)
        end if
!
! ---     MESSAGE/ALERTES MUMPS
        icntl(1) = to_mumps_int(ifm)
        icntl(2) = 0
        icntl(3) = 0
        icntl(4) = 1
        if ((niv .ge. 2) .or. (lverbose)) then
! ---     ICNTL(4) = 1/ERROR MESSAGES ONLY 2/ERRORS, WARNINGS, 3 PUIS 4
            icntl(3) = to_mumps_int(ifm)
            icntl(4) = 2
        end if
!
! ---     RAFFINEMENT ITERATIF ET ETUDE DE LA QUALITE
        icntl(10) = 0
        icntl(11) = 0
        cntl(2) = 0.d0
        slvk(11) = trim(adjustl(slvk(11)))
        select case (slvk(11))
        case ('MINI')
            icntl(10) = -2
            if (lquali) icntl(11) = 2
        case default
            if (lquali) then
                if (slvk(11) .eq. 'SANS') then
                    icntl(11) = 2
                else if (slvk(11) .eq. 'AUTO') then
                    icntl(10) = 4
                    cntl(2) = 1.d-14
                    icntl(11) = 2
                else if (slvk(11) .eq. 'FORCE') then
                    icntl(10) = 10
                    cntl(2) = 10.d-50
                    icntl(11) = 1
                    if (type .eq. 'S' .or. type .eq. 'C') then
                        cntl(2) = 1.d-38
                    end if
                end if
            end if
        end select

! --     AFFECTATIONS A LA SD MUMPS
! ---     PARAMETRE POUR RESOLUTIONS SIMULTANEES
        iret1 = 0
        iret2 = 0
        iret3 = 0
        kmumps_sparse(1) = '&&MUMPS.SPARSE1'
        kmumps_sparse(2) = '&&MUMPS.SPARSE2'
        kmumps_sparse(3) = '&&MUMPS.SPARSE3'
        call jeexin(kmumps_sparse(1), iret1)
        call jeexin(kmumps_sparse(2), iret2)
        call jeexin(kmumps_sparse(3), iret3)
        if (type .eq. 'S') then
            smpsk%icntl(1) = icntl(1)
            smpsk%icntl(2) = icntl(2)
            smpsk%icntl(3) = icntl(3)
            smpsk%icntl(4) = icntl(4)
            smpsk%icntl(10) = icntl(10)
            smpsk%icntl(11) = icntl(11)
            smpsk%cntl(2) = real(cntl(2), kind=4)
        else if (type .eq. 'C') then
            cmpsk%icntl(1) = icntl(1)
            cmpsk%icntl(2) = icntl(2)
            cmpsk%icntl(3) = icntl(3)
            cmpsk%icntl(4) = icntl(4)
            cmpsk%icntl(10) = icntl(10)
            cmpsk%icntl(11) = icntl(11)
            cmpsk%cntl(2) = real(cntl(2), kind=4)
        else if (type .eq. 'D') then
            dmpsk%icntl(1) = icntl(1)
            dmpsk%icntl(2) = icntl(2)
            dmpsk%icntl(3) = icntl(3)
            dmpsk%icntl(4) = icntl(4)
            dmpsk%icntl(10) = icntl(10)
            dmpsk%icntl(11) = icntl(11)
            dmpsk%cntl(2) = cntl(2)
!
! cas particulier des rhs par blocs, denses (nbrhs>1) ou sparses (nbrhs quelconque, meme 1)
! pour conserver les post-traitements qualite avec nbrhs=1
            if (nbrhs .ne. 1) then
                if ((nbrhs .lt. 0) .and. ((iret1*iret2*iret3) .ne. 0)) dmpsk%icntl(20) = 3
                dmpsk%icntl(27) = int(-abs(nbrhs), 4)
            end if
        else if (type .eq. 'Z') then
            zmpsk%icntl(1) = icntl(1)
            zmpsk%icntl(2) = icntl(2)
            zmpsk%icntl(3) = icntl(3)
            zmpsk%icntl(4) = icntl(4)
            zmpsk%icntl(10) = icntl(10)
            zmpsk%icntl(11) = icntl(11)
            zmpsk%cntl(2) = cntl(2)
        else
            ASSERT(.false.)
        end if
!
!       ------------------------------------------------
!        MAUVAISE OPTION
!       ------------------------------------------------
    else
        ASSERT(.false.)
    end if
    call jedema()
#endif
end subroutine
