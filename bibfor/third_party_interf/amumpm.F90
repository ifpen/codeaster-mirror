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

subroutine amumpm(ldist, kxmps, kmonit, impr, ifmump, &
                  klag2, type, lmd, epsmat, ktypr, &
                  lpreco, lmhpc, lbloc)
!
!
    implicit none
!--------------------------------------------------------------
! BUT : REMPLIR LES OBJETS F90 DE MUMPS REPRESENTANT LA MATRICE A PARTIR
!       DE CELLE DE CODE_ASTER.
!
! IN  LDIST  :  LOG   : LOGICAL MUMPS DISTRIBUE OR NOT
! IN  KXMPS  :   IN   : INDICE DE L'INSTANCE MUMPS DANS DMPS
! IN  KMONIT :  K24   : VECTEUR DE NOMS DES OBJ JEVEUX
! IN  IMPR   :  K14   : FLAG POUR IMPRESSION MATRICE
! IN  IFMUMP :   IN   : UNITE LOGIQUE POUR IMPRESSION FICHIER
! IN  KLAG2  :   K5   : PARAMETRE DE SOLVEUR/ELIM_LAGR
! IN  TYPE   :   K1   : TYPE DU POINTEUR R OU C
! IN  LMD    :  LOG   : LOGIQUE PRECISANT SI ON EST EN MATR_DISTRIBUEE
! IN  EPSMAT :   R8   : SEUIL DE FILTRAGE DES TERMES DE LA MATRICE
! IN  KTYPR  :   K8   : TYPE DE RESOLUTION MUMPS (SYMDEF...)
! IN  LPRECO :  LOG   : MUMPS EST-IL UTILISE COMME PRECONDITIONNEUR ?
! IN  LMHPC  :  LOG   : LOGIQUE PRECISANT SI ON EST EN MODE DISTRIBUE ASTERXX
! IN  LBLOC  :  LOG   : LOGIQUE PRECISANT SI ON EFFECTUE L ANALYSE PAR BLOCS
!---------------------------------------------------------------
! aslint: disable=W1501
!
#include "asterf.h"
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r4maem.h"
#include "asterc/r4miem.h"
#include "asterc/r8maem.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/amumpm_hpc.h"
#include "asterfort/asmpi_comm_jev.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibd.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/vecint.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: kxmps, ifmump
    aster_logical :: ldist, lmd, lpreco, lmhpc, lbloc
    real(kind=8) :: epsmat
    character(len=1) :: type
    character(len=5) :: klag2
    character(len=8) :: ktypr
    character(len=14) :: impr
    character(len=24) :: kmonit(12)
!
#ifdef ASTER_HAVE_MUMPS
#include "asterf_mumps.h"
!
!
    type(smumps_struc), pointer :: smpsk => null()
    type(cmumps_struc), pointer :: cmpsk => null()
    type(dmumps_struc), pointer :: dmpsk => null()
    type(zmumps_struc), pointer :: zmpsk => null()
    integer(kind=8) :: nsmdi, jsmhc, nsmhc, jdelg, n, n1, nz, nvale, jvale
    integer(kind=8) :: nlong, jvale2, nzloc, kterm, iterm, ifm, niv, k, maxnz2
    integer(kind=8) :: sym, iret, jcoll, iligl, jnulogl, ltot, iok, iok2, coltmp
    integer(kind=8) :: kzero, iad, ifiltr, vali(2), nbproc, nfilt1, nfilt2, nblk
    integer(kind=8) :: nfilt3, isizemu, nsizemu, rang, esizemu, jdeeq, iblock
    mumps_int :: nbeq, nz2, iligg, jcolg
    character(len=4) :: etam
    character(len=8) :: k8bid
    character(len=14) :: nonu
    character(len=16) :: k16bid, nomcmd
    character(len=19) :: nomat, nosolv
    character(len=24) :: kfiltr, kpiv, kpiv2, ksizemu, kblock
    real(kind=8) :: raux, rfiltr, epsmac, rmax, rmin, rtest
    complex(kind=8) :: caux
    aster_logical :: lmnsy, ltypr, lnn, lfiltr, lspd, eli2lg, lsimpl, lcmde
    aster_logical ::  lvbloc
    integer(kind=8), pointer :: smdi(:) => null()
    integer(kind=8), pointer :: nequ(:) => null()
    character(len=24), pointer :: refn(:) => null()
!
!-----------------------------------------------------------------------
    call jemarq()
    call infdbg('SOLVEUR', ifm, niv)
    if (lmhpc) then
        call amumpm_hpc(kxmps, kmonit, impr, ifmump, klag2, &
                        type, epsmat, ktypr, lpreco, lbloc)
        goto 999
    end if
! pour forcer le mode VERBOSE de l'option bloc de MUMPS
    lvbloc = .true.
    lvbloc = .false.
!
!       ------------------------------------------------
!        INITS
!       ------------------------------------------------
!
    epsmac = r8prem()
    nomat = nomats(kxmps)
    nosolv = nosols(kxmps)
    nonu = nonus(kxmps)
    etam = etams(kxmps)
    call jeveuo(nonu//'.NUME.REFN', 'L', vk24=refn)
    if (refn(4) .eq. "SIMPLE_LAGRANGE" .and. klag2 .ne. "NON") then
        call utmess('I', 'ALGELINE_10')
        klag2 = "NON"
    end if

! --- REMPLISSAGE DE DIFFERENTS OBJETS SUIVANT LE TYPE DU POINTEUR
! --- DE MUMPS: DMUMPS_STRUC OU ZMUMPS_STRUC
    if (type .eq. 'S') then
        smpsk => smps(kxmps)
        ltypr = .true.
        sym = smpsk%sym
        rmax = r4maem()*real(0.5d0, 4)
        rmin = r4miem()*real(2.0d0, 4)
        nbproc = smpsk%nprocs
        rang = smpsk%myid
        esizemu = 4
    else if (type .eq. 'C') then
        cmpsk => cmps(kxmps)
        ltypr = .false.
        sym = cmpsk%sym
        rmax = r4maem()*real(0.5d0, 4)
        rmin = r4miem()*real(2.0d0, 4)
        nbproc = cmpsk%nprocs
        rang = cmpsk%myid
        esizemu = 8
    else if (type .eq. 'D') then
        dmpsk => dmps(kxmps)
        ltypr = .true.
        sym = dmpsk%sym
        rmax = r8maem()*0.5d0
        rmin = r8miem()*2.0d0
        nbproc = dmpsk%nprocs
        rang = dmpsk%myid
        esizemu = 8
    else if (type .eq. 'Z') then
        zmpsk => zmps(kxmps)
        ltypr = .false.
        sym = zmpsk%sym
        rmax = r8maem()*0.5d0
        rmin = r8miem()*2.0d0
        nbproc = zmpsk%nprocs
        rang = zmpsk%myid
        esizemu = 16
    else
        ASSERT(.false.)
    end if
!
    if (lmd) then
        if (lpreco) then
            call jeveuo(nonu//'.NUML.NLGP', 'L', jnulogl)
        else
            call jeveuo(nonu//'.NUML.NULG', 'L', jnulogl)
        end if
    else
        jnulogl = 1
    end if
    if (ktypr(1:6) .eq. 'SYMDEF') then
        if ((type .eq. 'C') .or. (type .eq. 'Z')) then
            call utmess('F', 'FACTOR_80')
        end if
        lspd = .true.
    else
        lspd = .false.
    end if
!
    lsimpl = (type .eq. 'S') .or. (type .eq. 'C')
!
!     ------------------------------------------
!     DETERMINATION DE LA SYMETRIE DE LA MATRICE
!     ------------------------------------------
    call jelira(nomat//'.VALM', 'NMAXOC', nvale)
! --- LMNSY EST INDEPENDANT DE XMPSK%SYM POUR POUVOIR TRAITER
! --- DES CAS SYMETRIQUES EN MUMPS NON SYMETRIQUE
    if (nvale .eq. 1) then
        lmnsy = .false.
    else if (nvale .eq. 2) then
        lmnsy = .true.
    else
        ASSERT(.false.)
    end if
!
!
!       ------------------------------------------------
!        lecture d'adresses et de parametres preliminaires
!       ------------------------------------------------
    if (((rang .eq. 0) .and. (.not. ldist)) .or. (ldist)) then
        call jeveuo(nonu//'.SMOS.SMDI', 'L', vi=smdi)
        call jelira(nonu//'.SMOS.SMDI', 'LONMAX', nsmdi)
        call jeveuo(nonu//'.SMOS.SMHC', 'L', jsmhc)
        call jelira(nonu//'.SMOS.SMHC', 'LONMAX', nsmhc)
        if (lmd) then
            call jeveuo(nonu//'.NUML.DELG', 'L', jdelg)
            call jelira(nonu//'.NUML.DELG', 'LONMAX', n1)
        else
            call jeveuo(nonu//'.NUME.DELG', 'L', jdelg)
            call jelira(nonu//'.NUME.DELG', 'LONMAX', n1)
        end if
        call jeveuo(nonu//'.NUME.NEQU', 'L', vi=nequ)
        nbeq = to_mumps_int(nequ(1))
        ASSERT(n1 .eq. nsmdi)
! --- CALCUL DE N
        n = nsmdi
!
! --- GESTION DES BLOCS POUR ANALYSE (A PARTIR DE 5.4Consortium)
        if ((lbloc) .and. (rang .eq. 0)) then
            call jeveuo(nonu//'.NUME.DEEQ', 'L', jdeeq)
            if (lvbloc) then
! pour monitoring bloc
                do k = 1, n
                    write (ifm, *) 'k=', k, ' deeq=', &
                        zi(jdeeq+2*(k-1)), zi(jdeeq+2*(k-1)+1)
                end do
            end if
            kblock = '&&AMUMPM.BLOCKS'
            call wkvect(kblock, 'V V I', n+1, iblock)
            call vecint(n, -1_8, zi(iblock))
            nblk = 1
! premier bloc
            zi(iblock+nblk-1) = 1
            do k = 2, n
                if (zi(jdeeq+2*(k-1)) .gt. 0) then
                    if (zi(jdeeq+2*(k-1)) .eq. zi(jdeeq+2*(k-2))) then
! on ne fait rien, element du meme bloc, ddl physique ou Lagrange
                    else
! nouveau bloc
                        nblk = nblk+1
                        zi(iblock+nblk-1) = k
                    end if
                else
! nouveau bloc: lagrange pour CL egalite
                    nblk = nblk+1
                    zi(iblock+nblk-1) = k
                end if
            end do
! pour marquer la fin des blocs
            zi(iblock+nblk) = n+1

            if (lvbloc) then
! pour monitoring bloc
                write (ifm, *) 'n/nblock=', n, nblk
                do k = 1, nblk+1
                    write (ifm, *) 'block=', k, ' valeur=', zi(iblock+k-1)
                end do
            end if
            if (type .eq. 'S') then
                allocate (smpsk%blkptr(nblk+1))
                smpsk%nblk = to_mumps_int(nblk)
                do k = 1, nblk+1
                    smpsk%blkptr(k) = to_mumps_int(zi(iblock+k-1))
                end do
            else if (type .eq. 'C') then
                allocate (cmpsk%blkptr(nblk+1))
                cmpsk%nblk = to_mumps_int(nblk)
                do k = 1, nblk+1
                    cmpsk%blkptr(k) = to_mumps_int(zi(iblock+k-1))
                end do
            else if (type .eq. 'D') then
                allocate (dmpsk%blkptr(nblk+1))
                dmpsk%nblk = to_mumps_int(nblk)
                do k = 1, nblk+1
                    dmpsk%blkptr(k) = to_mumps_int(zi(iblock+k-1))
                end do
            else if (type .eq. 'Z') then
                allocate (zmpsk%blkptr(nblk+1))
                zmpsk%nblk = to_mumps_int(nblk)
                do k = 1, nblk+1
                    zmpsk%blkptr(k) = to_mumps_int(zi(iblock+k-1))
                end do
            end if
            call jedetr(kblock)
        end if

! --- GESTION ELIM_LAGR='LAGR2'
! --- on a essaye une basule automatique ELIM_LAGR='LAGR2'/'NON'
! --- en fonction de la proportion de lagranges. en fait, 'LAGR2'
! --- offre la plupart du temps le meilleur compromis:
! ---     cpu x memoire x qualite --> on ne programme pas de
! --- bascule, on laisse 'LAGR2' par defaut (pour l'instant)
        select case (klag2(1:5))
        case ('LAGR2')
            eli2lg = .true.
        case ('OUI', 'NON', 'XXXX')
            eli2lg = .false.
        case default
            ASSERT(.false.)
        end select
!
! --- CALCUL DE NZ2
        nz = smdi(n)
        ASSERT(nz .le. nsmhc)
        nz2 = to_mumps_int(nz)
        if (sym .eq. 0) nz2 = to_mumps_int(2*nz-n)
!
        call jeveuo(jexnum(nomat//'.VALM', 1_8), 'L', jvale)
        call jelira(jexnum(nomat//'.VALM', 1_8), 'LONMAX', nlong)
        ASSERT(nlong .eq. nz)
        if (lmnsy) then
            call jeveuo(jexnum(nomat//'.VALM', 2_8), 'L', jvale2)
            call jelira(jexnum(nomat//'.VALM', 2_8), 'LONMAX', nlong)
            ASSERT(nlong .eq. nz)
        end if
!
! ---  DETERMINATION DES TERMES DE REFERENCE POUR LE FILTRAGE
! ---  NFILT1 : TERMES RIGOUREUSEMENT NULS
! ---  NFILT2 : DEPASSEMENT DE CAPACITE PAR VALEUR MAX
! ---  NFILT3 : DEPASSEMENT DE CAPACITE PAR VALEUR MIN
        lfiltr = .false.
        kfiltr = '&&AMUMPM.FILTRAGE'
        nfilt1 = 0
        nfilt2 = 0
        nfilt3 = 0
        if (epsmat .gt. 0.d0) then
            lfiltr = .true.
            call jeexin(kfiltr, iret)
            if (iret .ne. 0) call jedetr(kfiltr)
            call wkvect(kfiltr, 'V V R', n, ifiltr)
            if (ltypr) then
                do k = 1, n
                    zr(ifiltr-1+k) = epsmat*abs(zr(jvale-1+smdi(k)))
                end do
            else
                do k = 1, n
                    zr(ifiltr-1+k) = epsmat*abs(zc(jvale-1+smdi(k)))
                end do
            end if
! --- SEUILLAGE DES TERMES DE FILTRAGE POUR EVITER LES VALEURS ABBERANTE
            do k = 1, n
                rtest = zr(ifiltr-1+k)
                if (rtest .lt. rmin) zr(ifiltr-1+k) = 0.d0
                if (rtest .gt. rmax) zr(ifiltr-1+k) = rmax
            end do
        else
            lfiltr = .false.
        end if
    end if
!
!       ------------------------------------------------
!       DETERMINATION DU NBRE LOCAL DE TERMES PAR PROC: NZLOC
!       EN CENTRALISE (LDIST=.FALSE.): PROC 0 GERE TOUS LES TERMES POUR
!         LES FOURNIR A MUMPS. CE DERNIER ENSUITE LES REDISPATCHE PAR
!         PAQUETS SUR TOUS LES AUTRES PROCS.
!       EN DISTRIBUE (LDIST=.TRUE.): CHAQUE PROC FOURNIT LES TERMES
!         DONT IL A LA RESPONSABILITE.
!       ------------------------------------------------
!
    if (((rang .eq. 0) .and. (.not. ldist)) .or. (ldist)) then
!
! --- VECTEURS AUXILIAIRES POUR FILTRER LES TERMES MATRICIELS
! --- KPIV: PROFIL TRIANGULAIRE INF
! --- KPIV2: IDEM SUP
        kpiv = '&&AMUMPM.TERMEOK'
        kpiv2 = '&&AMUMPM.TERMEOK2'
        call jeexin(kpiv, iad)
        if (iad .ne. 0) then
            ASSERT(.false.)
        else
            call wkvect(kpiv, 'V V S', nz, iok)
            do k = 1, nz
                zi4(iok+k-1) = 0
            end do
        end if
        if (sym .eq. 0) then
            call jeexin(kpiv2, iad)
            if (iad .ne. 0) then
                ASSERT(.false.)
            else
                call wkvect(kpiv2, 'V V S', nz, iok2)
                do k = 1, nz
                    zi4(iok2+k-1) = 0
                end do
            end if
        end if
!
        nzloc = 0
        jcoll = 1
        do kterm = 1, nz
!
! --- PREPARATION DES DONNES (NUM DE COLONNE, TERME DE FILTRAGE...)
            if (smdi(jcoll) .lt. kterm) jcoll = jcoll+1
            iligl = zi4(jsmhc-1+kterm)
            if (lfiltr) then
                rfiltr = zr(ifiltr-1+iligl)+zr(ifiltr-1+jcoll)
            else
                rfiltr = -1.d0
            end if
!
!
! --- PARTIE TRIANGULAIRE INF. SI REEL
            if (ltypr) then
                raux = zr(jvale-1+kterm)
                if (raux .ne. 0.d0) then
                    rtest = abs(raux)
                    if (rtest .gt. rfiltr) then
                        if (rtest .lt. rmin) then
                            nfilt3 = nfilt3+1
                            zi4(iok+kterm-1) = -2
                        else if (rtest .gt. rmax) then
                            nfilt2 = nfilt2+1
                            zi4(iok+kterm-1) = -1
                        else
                            zi4(iok+kterm-1) = 1
                        end if
                        nzloc = nzloc+1
                    end if
                else
! ---   TERME RIGOUREUSEMENT NUL
                    nfilt1 = nfilt1+1
                end if
!
            else
!
! --- PARTIE TRIANGULAIRE INF. SI COMPLEXE
                caux = zc(jvale-1+kterm)
                if (caux .ne. (0.d0, 0.d0)) then
                    rtest = abs(caux)
                    if (rtest .gt. rfiltr) then
                        if (rtest .lt. rmin) then
                            nfilt3 = nfilt3+1
                            zi4(iok+kterm-1) = -2
                        else if (rtest .gt. rmax) then
                            nfilt2 = nfilt2+1
                            zi4(iok+kterm-1) = -1
                        else
                            zi4(iok+kterm-1) = 1
                        end if
                        nzloc = nzloc+1
                    end if
                else
! ---   TERME RIGOUREUSEMENT NUL
                    nfilt1 = nfilt1+1
                end if
            end if

!
! --- PARTIE TRIANGULAIRE SUP. SI REEL
            if ((sym .eq. 0) .and. (iligl .ne. jcoll)) then
!
                if (ltypr) then
                    if (lmnsy) raux = zr(jvale2-1+kterm)
                    if (raux .ne. 0.d0) then
                        rtest = abs(raux)
                        if (rtest .gt. rfiltr) then
                            if (rtest .lt. rmin) then
                                nfilt3 = nfilt3+1
                                zi4(iok2+kterm-1) = -2
                            else if (rtest .gt. rmax) then
                                nfilt2 = nfilt2+1
                                zi4(iok2+kterm-1) = -1
                            else
                                zi4(iok2+kterm-1) = 1
                            end if
                            nzloc = nzloc+1
                        end if
                    else
! ---   TERME RIGOUREUSEMENT NUL
                        nfilt1 = nfilt1+1
                    end if
!
                else
!
! --- PARTIE TRIANGULAIRE SUP. SI COMPLEXE
                    if (lmnsy) caux = zc(jvale2-1+kterm)
                    if (caux .ne. (0.d0, 0.d0)) then
                        rtest = abs(caux)
                        if (rtest .gt. rfiltr) then
                            if (rtest .lt. rmin) then
                                nfilt3 = nfilt3+1
                                zi4(iok2+kterm-1) = -2
                            else if (rtest .gt. rmax) then
                                nfilt2 = nfilt2+1
                                zi4(iok2+kterm-1) = -1
                            else
                                zi4(iok2+kterm-1) = 1
                            end if
                            nzloc = nzloc+1
                        end if
                    else
! ---   TERME RIGOUREUSEMENT NUL
                        nfilt1 = nfilt1+1
                    end if
!
                end if
            end if
        end do
        nz2 = to_mumps_int(nzloc)

!
! --- GARDE-FOU POUR NE PAS TRANSMETTRE A MUMPS DES SYSTEMES VIDES
        maxnz2 = nzloc
        if (ldist) then
            call asmpi_comm_vect('MPI_MAX', 'I', nbval=1_8, sci=maxnz2)
        end if
        if ((maxnz2 .le. 0) .or. (n .le. 0)) then
            call utmess('F', 'FACTOR_41')
        end if
        if (niv .ge. 2) then
            write (ifm, *) '<AMUMPM>     NZLOC: ', nzloc
            write (ifm, *) '       TERMES NULS: ', nfilt1
            write (ifm, *) '   UNDER/OVERFLOWS: ', nfilt3, '/', nfilt2
        end if
! if (ldist...
    end if
!
!       ------------------------------------------------
!       ALLOCATION DES OBJETS MUMPS F90
!       ------------------------------------------------
!
! ---   OBJET JEVEUX STOCKANT, PAR PROC, LA TAILLE DES OBJETS ALLOUES
! ---   POUR MUMPS. UTILE A AMUMPU. ON SUPPOSE UN SEUL RHS.
! ---   EN FIN DE ROUTINE, TOUS LES PROCS CONNAISSENT CE VECTEUR
! ----  VIA UN MPI_ALLREDUCE + SUM.
    ksizemu = '&&TAILLE_OBJ_MUMPS'
    call jeexin(ksizemu, iret)
    if (iret .eq. 0) then
        call wkvect(ksizemu, 'V V I', nbproc, isizemu)
    else
        call jeveuo(ksizemu, 'E', isizemu)
    end if
    do k = 1, nbproc
        zi(isizemu+k-1) = 0
    end do
    nsizemu = 0
    if (((rang .eq. 0) .and. (.not. ldist)) .or. (ldist)) then
    if (ldist) then
    if (type .eq. 'S') then
        if (lmd) then
            smpsk%n = nbeq
        else
            smpsk%n = to_mumps_int(n)
        end if
        smpsk%nnz_loc = nz2
        allocate (smpsk%irn_loc(nz2))
        allocate (smpsk%jcn_loc(nz2))
        allocate (smpsk%a_loc(nz2))
    else if (type .eq. 'C') then
        cmpsk%n = to_mumps_int(n)
        cmpsk%nnz_loc = nz2
        allocate (cmpsk%irn_loc(nz2))
        allocate (cmpsk%jcn_loc(nz2))
        allocate (cmpsk%a_loc(nz2))
    else if (type .eq. 'D') then
        if (lmd) then
            dmpsk%n = nbeq
        else
            dmpsk%n = to_mumps_int(n)
        end if
        dmpsk%nnz_loc = nz2
        allocate (dmpsk%irn_loc(nz2))
        allocate (dmpsk%jcn_loc(nz2))
        allocate (dmpsk%a_loc(nz2))
    else if (type .eq. 'Z') then
        zmpsk%n = to_mumps_int(n)
        zmpsk%nnz_loc = nz2
        allocate (zmpsk%irn_loc(nz2))
        allocate (zmpsk%jcn_loc(nz2))
        allocate (zmpsk%a_loc(nz2))
    else
        ASSERT(.false.)
    end if
    if (lmd) then
        nsizemu = nz2*(4+4+esizemu)+esizemu*nbeq
    else
        nsizemu = nz2*(4+4+esizemu)+esizemu*n
    end if
    else
    if (type .eq. 'S') then
        smpsk%n = to_mumps_int(n)
        smpsk%nnz = nz2
        allocate (smpsk%irn(nz2))
        allocate (smpsk%jcn(nz2))
        allocate (smpsk%a(nz2))
    else if (type .eq. 'C') then
        cmpsk%n = to_mumps_int(n)
        cmpsk%nnz = nz2
        allocate (cmpsk%irn(nz2))
        allocate (cmpsk%jcn(nz2))
        allocate (cmpsk%a(nz2))
    else if (type .eq. 'D') then
        dmpsk%n = to_mumps_int(n)
        dmpsk%nnz = nz2
        allocate (dmpsk%irn(nz2))
        allocate (dmpsk%jcn(nz2))
        allocate (dmpsk%a(nz2))
    else if (type .eq. 'Z') then
        zmpsk%n = to_mumps_int(n)
        zmpsk%nnz = nz2
        allocate (zmpsk%irn(nz2))
        allocate (zmpsk%jcn(nz2))
        allocate (zmpsk%a(nz2))
    else
        ASSERT(.false.)
    end if
    nsizemu = nz2*(4+4+esizemu)+esizemu*n
    end if
    zi(isizemu+rang) = 1+(nsizemu/(1024*1024))
!
!       ------------------------------------------------
!       INTERPRETATION DES PBS RENCONTRES LORS DU FILTRAGE
!       REMPLISSAGE DE LA MATRICE
!       EN CAS D'OVERFLOW:
!         SI SOLVEUR DIRECT: UTMESS_F
!         SI SIMPLE PRECISION : ON SEUILLE LES TERMES IMPACTES
!                               LORS DU REMPLISSAGE EFFECTIF.
!       EN CAS D'UNDERFLOW: ON FIXE A ZERO.
!       ------------------------------------------------
!
    if (.not. lsimpl) then
        vali(1) = nfilt3
        vali(2) = nfilt2
        do kterm = 1, nz
            if (zi4(iok+kterm-1) .eq. -1) then
                call utmess('F', 'FACTOR_78', ni=2_8, vali=vali)
            end if
        end do
        if (sym .eq. 0) then
            do kterm = 1, nz
                if (zi4(iok2+kterm-1) .eq. -1) then
                    call utmess('F', 'FACTOR_78', ni=2_8, vali=vali)
                end if
            end do
        end if
    end if
!
! --- REMPLISSAGE EFFECTIF DES TERMES DE LA MATRICE
    jcoll = 1
    iterm = 0
    do kterm = 1, nz
!
        if (smdi(jcoll) .lt. kterm) jcoll = jcoll+1
        iligl = zi4(jsmhc-1+kterm)
        lnn = .false.
! --- PARTIE TRIANGULAIRE INF. SI REEL
        if (ltypr) then
            if (zi4(iok+kterm-1) .eq. 1) then
                lnn = .true.
                raux = zr(jvale-1+kterm)
            else if (zi4(iok+kterm-1) .eq. -1) then
                lnn = .true.
                raux = zr(jvale-1+kterm)
                raux = rmax*sign(1.d0, raux)
            end if
        else
! --- PARTIE TRIANGULAIRE INF. SI COMPLEXE
            if (zi4(iok+kterm-1) .eq. 1) then
                lnn = .true.
                caux = zc(jvale-1+kterm)
            else if (zi4(iok+kterm-1) .eq. -1) then
                lnn = .true.
                caux = zc(jvale-1+kterm)
                caux = rmax*dcmplx(1.d0*sign(1.d0, dble(caux)), 1.d0*sign(1.d0, imag(caux)))
            end if
        end if
        if (lmd) then
            iligg = to_mumps_int(zi(jnulogl+iligl-1))
            jcolg = to_mumps_int(zi(jnulogl+jcoll-1))
        else
            iligg = to_mumps_int(iligl)
            jcolg = to_mumps_int(jcoll)
        end if
        if ((sym .ne. 0) .and. (iligg .ge. jcolg)) then
            coltmp = jcolg
            jcolg = to_mumps_int(iligg)
            iligg = to_mumps_int(coltmp)
        end if
!
! ---- PARTIE TRIANGULAIRE INF. TERME NON NUL
        if (lnn) then
            iterm = iterm+1
            if (ldist) then
                if (type .eq. 'S') then
                    smpsk%irn_loc(iterm) = iligg
                    smpsk%jcn_loc(iterm) = jcolg
                    smpsk%a_loc(iterm) = real(raux, kind=4)
                else if (type .eq. 'C') then
                    cmpsk%irn_loc(iterm) = iligg
                    cmpsk%jcn_loc(iterm) = jcolg
                    cmpsk%a_loc(iterm) = cmplx(caux, kind=4)
                else if (type .eq. 'D') then
                    dmpsk%irn_loc(iterm) = iligg
                    dmpsk%jcn_loc(iterm) = jcolg
                    dmpsk%a_loc(iterm) = raux
                else if (type .eq. 'Z') then
                    zmpsk%irn_loc(iterm) = iligg
                    zmpsk%jcn_loc(iterm) = jcolg
                    zmpsk%a_loc(iterm) = caux
                else
                    ASSERT(.false.)
                end if
            else
                if (type .eq. 'S') then
                    smpsk%irn(iterm) = iligg
                    smpsk%jcn(iterm) = jcolg
                    smpsk%a(iterm) = real(raux, kind=4)
                else if (type .eq. 'C') then
                    cmpsk%irn(iterm) = iligg
                    cmpsk%jcn(iterm) = jcolg
                    cmpsk%a(iterm) = cmplx(caux, kind=4)
                else if (type .eq. 'D') then
                    dmpsk%irn(iterm) = iligg
                    dmpsk%jcn(iterm) = jcolg
                    dmpsk%a(iterm) = raux
                else if (type .eq. 'Z') then
                    zmpsk%irn(iterm) = iligg
                    zmpsk%jcn(iterm) = jcolg
                    zmpsk%a(iterm) = caux
                else
                    ASSERT(.false.)
                end if
            end if
            kzero = 0
            if (eli2lg) then
! ------      ON ELIMINE LE DERNIER TERME DE A/A_loc SI LAG2
                if (zi(jdelg-1+iligl) .eq. -1) then
                    if (jcoll .eq. iligl) kzero = 1
                    if (zi(jdelg-1+jcoll) .eq. -2) kzero = 1
                end if
                if (zi(jdelg-1+iligl) .eq. -2) then
                    if (jcoll .ne. iligl) kzero = 1
                end if
                if (zi(jdelg-1+jcoll) .eq. -2) then
                    if (jcoll .ne. iligl) kzero = 1
                end if
                if (kzero .eq. 1) iterm = iterm-1
            end if
! --- SI RESOLUTION SPD DEMANDEE ET TERME NEGATIF OU NUL ON S'ARRETE
            if ((lspd) .and. (kzero .eq. 0)) then
                if ((iligg .eq. jcolg) .and. (raux .lt. epsmac)) then
                    call utmess('F', 'FACTOR_84')
                end if
            end if
        end if

!
! --- PARTIE TRIANGULAIRE SUP. SI REEL
        if ((sym .eq. 0) .and. (iligl .ne. jcoll)) then
!
            lnn = .false.
            if (ltypr) then
                if (zi4(iok2+kterm-1) .eq. 1) then
                    lnn = .true.
                    if (lmnsy) raux = zr(jvale2-1+kterm)
                else if (zi4(iok2+kterm-1) .eq. -1) then
                    lnn = .true.
                    raux = zr(jvale-1+kterm)
                    raux = rmax*sign(1.d0, raux)
                end if
            else
! --- PARTIE TRIANGULAIRE SUP. SI COMPLEXE
                if (zi4(iok2+kterm-1) .eq. 1) then
                    lnn = .true.
                    if (lmnsy) caux = zc(jvale2-1+kterm)
                else if (zi4(iok2+kterm-1) .eq. -1) then
                    lnn = .true.
                    caux = zc(jvale-1+kterm)
                    caux = rmax*dcmplx(1.d0*sign(1.d0, dble(caux)), 1.d0*sign(1.d0, imag(caux)))
                end if
            end if
!
! ---- PARTIE TRIANGULAIRE SUP. TERME NON NUL
            if (lnn) then
                iterm = iterm+1
                if (lmd) then
                    iligg = to_mumps_int(zi(jnulogl+iligl-1))
                    jcolg = to_mumps_int(zi(jnulogl+jcoll-1))
                else
                    iligg = to_mumps_int(iligl)
                    jcolg = to_mumps_int(jcoll)
                end if
                if (ldist) then
                    if (type .eq. 'S') then
                        smpsk%irn_loc(iterm) = jcolg
                        smpsk%jcn_loc(iterm) = iligg
                        smpsk%a_loc(iterm) = real(raux, kind=4)
                    else if (type .eq. 'C') then
                        cmpsk%irn_loc(iterm) = jcolg
                        cmpsk%jcn_loc(iterm) = iligg
                        cmpsk%a_loc(iterm) = cmplx(caux, kind=4)
                    else if (type .eq. 'D') then
                        dmpsk%irn_loc(iterm) = jcolg
                        dmpsk%jcn_loc(iterm) = iligg
                        dmpsk%a_loc(iterm) = raux
                    else if (type .eq. 'Z') then
                        zmpsk%irn_loc(iterm) = jcolg
                        zmpsk%jcn_loc(iterm) = iligg
                        zmpsk%a_loc(iterm) = caux
                    else
                        ASSERT(.false.)
                    end if
                else
                    if (type .eq. 'S') then
                        smpsk%irn(iterm) = jcolg
                        smpsk%jcn(iterm) = iligg
                        smpsk%a(iterm) = real(raux, kind=4)
                    else if (type .eq. 'C') then
                        cmpsk%irn(iterm) = jcolg
                        cmpsk%jcn(iterm) = iligg
                        cmpsk%a(iterm) = cmplx(caux, kind=4)
                    else if (type .eq. 'D') then
                        dmpsk%irn(iterm) = jcolg
                        dmpsk%jcn(iterm) = iligg
                        dmpsk%a(iterm) = raux
                    else if (type .eq. 'Z') then
                        zmpsk%irn(iterm) = jcolg
                        zmpsk%jcn(iterm) = iligg
                        zmpsk%a(iterm) = caux
                    else
                        ASSERT(.false.)
                    end if
                end if
                if (eli2lg) then
                    if (kzero .eq. 1) iterm = iterm-1
                end if
! --- SI RESOLUTION SPD DEMANDEE ET TERME NEGATIF OU NUL ON S'ARRETE
                if ((lspd) .and. (kzero .eq. 0)) then
                    if ((iligg .eq. jcolg) .and. (raux .lt. epsmac)) then
                        call utmess('F', 'FACTOR_84')
                    end if
                end if
            end if
        end if
! ---FIN DE LA BOUCLE SUR NZ
    end do
!
    ASSERT(iterm .le. nz2)
    nz2 = to_mumps_int(iterm)
    if (ldist) then
    if (type .eq. 'S') then
        smpsk%nnz_loc = nz2
    else if (type .eq. 'C') then
        cmpsk%nnz_loc = nz2
    else if (type .eq. 'D') then
        dmpsk%nnz_loc = nz2
    else if (type .eq. 'Z') then
        zmpsk%nnz_loc = nz2
    else
        ASSERT(.false.)
    end if
    else
    if (type .eq. 'S') then
        smpsk%nnz = nz2
    else if (type .eq. 'C') then
        cmpsk%nnz = nz2
    else if (type .eq. 'D') then
        dmpsk%nnz = nz2
    else if (type .eq. 'Z') then
        zmpsk%nnz = nz2
    else
        ASSERT(.false.)
    end if
    end if
    ASSERT(iligl .eq. n)
    ASSERT(jcoll .eq. n)
    call jelibe(nonu//'.SMOS.SMDI')
    call jelibe(nonu//'.SMOS.SMHC')
    call jelibe(jexnum(nomat//'.VALM', 1_8))
    if (lmnsy) call jelibe(jexnum(nomat//'.VALM', 2_8))
    call jelibe(nonu//'.NUME.DELG')
    if (lmd) then
        call jelibe(nonu//'.NUML.DELG')
    end if
!
    if (niv .ge. 2) then
! --- TEST POUR EVITER LE MONITORING DES CMDES ECLATEES ET DE LDLT_SP
!     LES OBJETS TEMPORAIRES DE MONITORING SONT EFFACES A CHAQUE
!     FIN DE COMMANDE (NUM_DDL/FACTORISER/RESOUDRE)
        call getres(k8bid, k16bid, nomcmd)
        if ((nomcmd(1:8) .eq. 'NUME_DDL') .or. (nomcmd(1:10) .eq. 'FACTORISER') .or. &
            (nomcmd(1:8) .eq. 'RESOUDRE')) then
            lcmde = .true.
        else
            lcmde = .false.
        end if
        if ((.not. lcmde) .and. (.not. lpreco)) then
            call jeexin(kmonit(1), iret)
            if (iret .ne. 0) then
! --- CAS CMDE STD AVEC MUMPS SOLVEUR DIRECT
                call jeveuo(kmonit(1), 'E', iad)
                zi(iad+rang) = nz2
            else
! --- L'OBJET KMONIT(1) DEVRAIT EXISTER
                ASSERT(.false.)
            end if
        end if
    end if
!
!
!       ------------------------------------------------
!       IMPRESSION DE LA MATRICE (SI DEMANDEE) :
!       ------------------------------------------------
    if (impr(1:3) .eq. 'OUI') then
        write (ifmump, *) sym, '   : SYM', rang, '   : RANG'
        write (ifmump, *) n, '   : N'
        if (ldist) then
            write (ifmump, *) nz2, '   : NZ_loc'
        else
            write (ifmump, *) nz2, '   : NZ'
        end if
        if (type .eq. 'S') then
            do k = 1, nz2
                if (ldist) then
                    write (ifmump, *) smpsk%irn_loc(k), smpsk%jcn_loc(k), smpsk%a_loc(k)
                else
                    write (ifmump, *) smpsk%irn(k), smpsk%jcn(k), smpsk%a(k)
                end if
            end do
        else if (type .eq. 'C') then
            do k = 1, nz2
                if (ldist) then
                    write (ifmump, *) cmpsk%irn_loc(k), cmpsk%jcn_loc(k), cmpsk%a_loc(k)
                else
                    write (ifmump, *) cmpsk%irn(k), cmpsk%jcn(k), cmpsk%a(k)
                end if
            end do
        else if (type .eq. 'D') then
            do k = 1, nz2
                if (ldist) then
                    write (ifmump+rang, *) dmpsk%irn_loc(k), dmpsk%jcn_loc(k), dmpsk%a_loc(k)
                else
                    write (ifmump, *) dmpsk%irn(k), dmpsk%jcn(k), dmpsk%a(k)
                end if
            end do
        else if (type .eq. 'Z') then
            do k = 1, nz2
                if (ldist) then
                    write (ifmump, *) zmpsk%irn_loc(k), zmpsk%jcn_loc(k), zmpsk%a_loc(k)
                else
                    write (ifmump, *) zmpsk%irn(k), zmpsk%jcn(k), zmpsk%a(k)
                end if
            end do
        else
            ASSERT(.false.)
        end if
        if (ldist) then
            write (ifmump, *) 'MUMPS FIN A_loc'
        else
            write (ifmump, *) 'MUMPS FIN A'
        end if

!  -------   VIDANGE DES BUFFERS D'IMPRESSION
        flush (ifmump+rang)

    end if
! FIN DU IF LDIST
    end if
!
! --- COMMUNICATION DU VECTEUR KSIZEMU A TOUS LES PROCS
    call asmpi_comm_jev('MPI_SUM', ksizemu)
!
! --- NETTOYAGE VECTEURS TEMPORAIRES LOCAUX
    if (((rang .eq. 0) .and. (.not. ldist)) .or. (ldist)) then
        call jeexin(kfiltr, iret)
        if (iret .ne. 0) call jedetr(kfiltr)
        call jedetr(kpiv)
        if (sym .eq. 0) call jedetr(kpiv2)
    end if
!
! --- DECHARGEMENT CIBLE D'OBJETS JEVEUX
    call jelibd(nonu//'.SMOS.SMDI', ltot)
    call jelibd(nonu//'.SMOS.SMHC', ltot)
    call jelibd(nonu//'.NUME.DEEQ', ltot)
    call jelibd(nonu//'.NUME.NUEQ', ltot)
    call jelibd(nonu//'.NUME.LILI', ltot)
    call jelibd(jexnum(nomat//'.VALM', 1_8), ltot)
    if (lmnsy) call jelibd(jexnum(nomat//'.VALM', 2_8), ltot)
    if (lmd) then
        call jelibd(nonu//'.NUML.NULG', ltot)
        call jelibd(nonu//'.NUML.DELG', ltot)
    else
        call jelibd(nonu//'.NUME.DELG', ltot)
    end if
!
999 continue
    call jedema()
#endif
end subroutine
