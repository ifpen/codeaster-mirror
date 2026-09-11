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
subroutine amumpp(option, nbsol, kxmps, ldist, type, &
                  impr, ifmump, eli2lg, rsolu, csolu, &
                  vcine, prepos, lpreco, lmhpc)
!
!
    use iso_c_binding, only: c_ptr, c_f_pointer
    implicit none
!-----------------------------------------------------------------------
! BUT : ROUTINE DE PRE/POST-TRAITEMENT DE LA SOLUTION ET DU
!       SECOND MEMBRE POUR AMUMPS/C/D/Z
!
! IN  OPTION :   IN   : OPTION D'UTILISATION.
! IN  NBSOL  :   IN   : NBRE DE SYSTEMES A RESOUDRE
! IN  KXMPS  :   IN   : INDICE DE L'INSTANCE MUMPS DANS DMPS
! IN  LDIST  :  LOG   : LOGICAL MUMPS DISTRIBUE OR NOT
! IN  TYPE   :   K1   : TYPE DU POINTEUR R OU C
! IN  IMPR   :  K14   : FLAG POUR IMPRESSION MATRICE
! IN  IFMUMP :   IN   : UNITE LOGIQUE POUR IMPRESSION FICHIER
! IN  ELI2LG :  LOG   : LOGICAL POUR NE LAISSER QU'1 LAGRANGE ACTIF
! I/O RSOLU  :    R   : EN ENTREE : SECONDS MEMBRES REELS
!                     : EN SORTIE : SOLUTIONS
! I/O CSOLU  :    C   : EN ENTREE : SECONDS MEMBRES COMPLEXES
!                     : EN SORTIE : SOLUTIONS
! IN  VCINE  :  K19   : NOM DU CHAM_NO DE CHARGEMENT CINEMATIQUE
! IN  PREPOS :  LOG   : SI .TRUE. ON FAIT LES PRE ET POSTTRAITEMENTS DE
!           MISE A L'ECHELLE DU RHS ET DE LA SOLUTION (MRCONL) ET DE LA
!           PRISE EN COMPTE DES AFFE_CHAR_CINE (CSMBGG).
!           SI .FALSE. ON NE LES FAIT PAS (PAR EXEMPLE EN MODAL).
! IN  LPRECO :  LOG   : MUMPS EST-IL UTILISE COMME PRECONDITIONNEUR ?
!-----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "asterf.h"
#include "jeveux.h"
#include "asterc/r4maem.h"
#include "asterc/r4miem.h"
#include "asterc/r8maem.h"
#include "asterc/r8miem.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/csmbgg.h"
#include "asterfort/filter_rhs_c.h"
#include "asterfort/filter_rhs.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/jgetptc.h"
#include "asterfort/mcconl.h"
#include "asterfort/mrconl.h"
#include "asterfort/mtdscr.h"
#include "asterfort/nudlg2.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
#include "blas/zcopy.h"
    integer(kind=4) :: option
    integer(kind=8) :: nbsol, kxmps, ifmump
    aster_logical :: ldist, eli2lg, prepos, lpreco, lmhpc
    character(len=1) :: type
    character(len=14) :: impr
    character(len=19) :: vcine
    real(kind=8) :: rsolu(*)
    complex(kind=8) :: csolu(*)
!
#ifdef ASTER_HAVE_MUMPS
#include "asterf_mumps.h"
    type(smumps_struc), pointer :: smpsk => null()
    type(cmumps_struc), pointer :: cmpsk => null()
    type(dmumps_struc), pointer :: dmpsk => null()
    type(zmumps_struc), pointer :: zmpsk => null()
    integer(kind=4) :: n, nnbsol
    integer(kind=8) :: rang, lmat, i, ierd, idvalc, k, ifm, niv
    integer(kind=8) :: jj, nbeql, nuglo, jrhs, j, dime, ntot, nec, nlili, ili
    integer(kind=8) :: jrefn, jmlogl, jdeeq, nuno, nucmp, lonmax, iret
    integer(kind=8) :: nbno_prno, ino, ieq, nbcmp, jdee2, idprn1, idprn2, jconl
    integer(kind=8) :: iret1, iret2, iret3, js1, js2, js3, izrhs, i1, i2, nzrhs
    character(len=1) :: rouc
    character(len=4) :: etam
    character(len=8) :: mesh, k8bid
    character(len=14) :: nonu
    character(len=19) :: nomat, nosolv
    character(len=24) :: vcival, nonulg, kmumps_sparse(3), conl
    aster_logical :: ltypr, lverif, lrhs_sparse, lverbose
    aster_logical, parameter :: l_debug = ASTER_FALSE
    real(kind=8) :: rr4max, raux, rmin, rmax, rtest, valr(2)
    complex(kind=8) :: cbid, caux
    integer(kind=8), pointer :: delg(:) => null()
    integer(kind=8), pointer :: dlg2(:) => null()
    integer(kind=8), pointer :: nequ(:) => null()
    integer(kind=8), pointer :: pddl(:) => null()
    integer(kind=8), pointer :: nulg(:) => null()
    real(kind=8), pointer :: rsolu2(:) => null()
    complex(kind=8), pointer :: csolu2(:) => null()
    type(c_ptr) :: pteur_c
    blas_int :: b_incx, b_incy, b_n
    cbid = dcmplx(0.d0, 0.d0)
#define zzprno(ili,nunoel,l) zi(idprn1-1+zi(idprn2+ili-1)+ (nunoel-1)* (nec+2)+l-1)
!
!-----------------------------------------------------------------------
    call jemarq()
    call infdbg('SOLVEUR', ifm, niv)
! pour verifier la norme de la solution produite
    lverif = .true.
    lverif = .false.
! pour forcer le mode VERBOSE en multiple RHS creux
    lverbose = .true.
    lverbose = .false.
!
!     ------------------------------------------------
!     INITS
!     ------------------------------------------------
    rr4max = r4maem()
    if (type .eq. 'S') then
        smpsk => smps(kxmps)
        rang = smpsk%myid
        n = smpsk%n
        smpsk%nrhs = to_mumps_int(nbsol)
        smpsk%lrhs = to_mumps_int(n)
        ltypr = .true.
        rmax = r4maem()*real(0.5d0, 4)
        rmin = r4miem()*real(2.0d0, 4)
    else if (type .eq. 'C') then
        cmpsk => cmps(kxmps)
        rang = cmpsk%myid
        n = cmpsk%n
        cmpsk%nrhs = to_mumps_int(nbsol)
        cmpsk%lrhs = to_mumps_int(n)
        ltypr = .false.
        rmax = r4maem()*real(0.5d0, 4)
        rmin = r4miem()*real(2.0d0, 4)
    else if (type .eq. 'D') then
        dmpsk => dmps(kxmps)
        rang = dmpsk%myid
        n = dmpsk%n
        dmpsk%nrhs = to_mumps_int(nbsol)
        dmpsk%lrhs = to_mumps_int(n)
        ltypr = .true.
        rmax = r8maem()*0.5d0
        rmin = r8miem()*2.0d0
    else if (type .eq. 'Z') then
        zmpsk => zmps(kxmps)
        rang = zmpsk%myid
        n = zmpsk%n
        zmpsk%nrhs = to_mumps_int(nbsol)
        zmpsk%lrhs = to_mumps_int(n)
        ltypr = .false.
        rmax = r8maem()*0.5d0
        rmin = r8miem()*2.0d0
    else
        ASSERT(.false.)
    end if
    nnbsol = int(n*nbsol, 4)
    nomat = nomats(kxmps)
    nosolv = nosols(kxmps)
    nonu = nonus(kxmps)
    etam = etams(kxmps)
!
    vcival = vcine//'.VALE'
    call jeexin(vcival, ierd)
    call mtdscr(nomat)
    call jeveuo(nomat//'.&INT', 'E', lmat)
!
    if (lmhpc) then
        ASSERT(nbsol .eq. 1)
        call jeveuo(nonu//'.NUME.PDDL', 'L', vi=pddl)
        call jeveuo(nonu//'.NUME.NULG', 'L', vi=nulg)
        call jeveuo(nonu//'.NUME.NEQU', 'L', vi=nequ)
        nbeql = nequ(1)
    else
        nbeql = n
    end if
!
! Pour appel a MUMPS solve en mode sparse: initialisations
    kmumps_sparse(1) = '&&MUMPS.SPARSE1'
    kmumps_sparse(2) = '&&MUMPS.SPARSE2'
    kmumps_sparse(3) = '&&MUMPS.SPARSE3'
    lrhs_sparse = .false.
    iret1 = 0
    iret2 = 0
    iret3 = 0
    call jeexin(kmumps_sparse(1), iret1)
    call jeexin(kmumps_sparse(2), iret2)
    call jeexin(kmumps_sparse(3), iret3)
    if (((iret1*iret2*iret3) .ne. 0) .and. (.not. lmhpc) .and. (ierd .eq. 0)) then
        ASSERT(type .eq. 'D')
        ASSERT(.not. lpreco)
        ASSERT(prepos)
        lrhs_sparse = .true.
    else
        lrhs_sparse = .false.
    end if
!
!   Adresses needed to get the solution wrt nodes and dof numbers (see below)
    if (l_debug) then
        call jeveuo(nonu//'.NUME.REFN', 'L', jrefn)
        call jeveuo(nonu//'.NUME.DEEQ', 'L', jdeeq)
        call wkvect(nonu//'.NUME.DEE2', 'V V I', int(2*nnbsol, 8), jdee2)
        mesh = zk24(jrefn) (1:8)
        if (lmhpc) then
            nonulg = mesh//'.NUNOLG'
            call jeveuo(nonulg, 'L', jmlogl)
        end if
!
        call jeveuo(mesh//'.DIME', 'L', dime)
        call jelira(jexnum(nonu//'.NUME.PRNO', 1_8), 'LONMAX', ntot, k8bid)
        call jeveuo(nonu//'.NUME.PRNO', 'L', idprn1)
        call jeveuo(jexatr(nonu//'.NUME.PRNO', 'LONCUM'), 'L', idprn2)
        nec = ntot/zi(dime)-2
        call jelira(nonu//'.NUME.PRNO', 'NMAXOC', nlili)
        do ili = 1, nlili
            call jelira(jexnum(nonu//'.NUME.PRNO', ili), 'LONMAX', lonmax)
            nbno_prno = lonmax/(nec+2)
            do ino = 1, nbno_prno
                ieq = zzprno(ili, ino, 1)
                nbcmp = zzprno(ili, ino, 2)
                do j = 1, nbcmp
                    if (zi(jdeeq+2*(ieq-1)) .eq. 0) then
                        if (ili .eq. 1) then
                            zi(jdee2+2*(ieq-1)) = ino
                            zi(jdee2+2*(ieq-1)+1) = -j
                        else
                            zi(jdee2+2*(ieq-1)) = -ino
                            zi(jdee2+2*(ieq-1)+1) = j
                        end if
                    else
                        zi(jdee2+2*(ieq-1)) = zi(jdeeq+2*(ieq-1))
                        zi(jdee2+2*(ieq-1)+1) = zi(jdeeq+2*(ieq-1)+1)
                    end if
                    ieq = ieq+1
                end do
            end do
        end do
        jdeeq = jdee2
    end if
!
!
!! Pour appel a MUMPS solve en mode sparse: preparation des donnees
    if ((option .eq. 0) .and. (lrhs_sparse)) then
        call jeveuo(kmumps_sparse(1), 'L', js1)
        js1 = js1-1
        call jeveuo(kmumps_sparse(2), 'L', js2)
        js2 = js2-1
        call jeveuo(kmumps_sparse(3), 'L', js3)
        izrhs = zi(js3)
        i1 = zi(js3+izrhs)
        i2 = zi(js3+izrhs+nbsol)
        nzrhs = i2-i1
! Traitements simplifies
        if (rang .eq. 0) then
            conl = zk24(zi(lmat+1)) (1:19)//'.CONL'
            call jeexin(conl, iret)
            ASSERT(iret .ne. 0)
            call jeveuo(conl, 'L', jconl)
            jconl = jconl-1
            allocate (dmpsk%rhs(nnbsol))
            dmpsk%nz_rhs = int(nzrhs, 4)
            allocate (dmpsk%rhs_sparse(nzrhs))
            allocate (dmpsk%irhs_sparse(nzrhs))
            do i = 0, dmpsk%nz_rhs-1
                dmpsk%rhs_sparse(i+1) = zr(js1+i1+i)*zr(jconl+zi(js2+i1+i))
                dmpsk%irhs_sparse(i+1) = int(zi(js2+i1+i), 4)
            end do
            allocate (dmpsk%irhs_ptr(nbsol+1))
            dmpsk%irhs_ptr(1) = 1
            do i = 0, nbsol
                dmpsk%irhs_ptr(i+1) = int(zi(js3+izrhs+i)-(i1-1), 4)
            end do
!
            if (lverbose) then
                write (6, *) '<amumpp> ibatch/batch_size/nz_rhs/i1/i2=', izrhs, nbsol, nzrhs, i1, i2
                do i = 1, dmpsk%nz_rhs
                    write (6, *) i, dmpsk%rhs_sparse(i), dmpsk%irhs_sparse(i)
                end do
                do i = 1, nbsol+1
                    write (6, *) i, dmpsk%irhs_ptr(i)
                end do
            end if
        end if
    else if ((option .eq. 0) .and. (.not. lrhs_sparse)) then
!
        if (rang .eq. 0) then
            if (type .eq. 'S') then
                allocate (smpsk%rhs(nnbsol))
            else if (type .eq. 'C') then
                allocate (cmpsk%rhs(nnbsol))
            else if (type .eq. 'D') then
                allocate (dmpsk%rhs(nnbsol))
            else if (type .eq. 'Z') then
                allocate (zmpsk%rhs(nnbsol))
            else
                ASSERT(.false.)
            end if
        end if
!
!       ------------------------------------------------
!        PRETRAITEMENTS ASTER DU/DES SECONDS MEMBRES :
!       ------------------------------------------------
!
!        --- PAS DE PRETRAITEMENT SI NON DEMANDE
        if (.not. lpreco .and. prepos) then
!
            if (rang .eq. 0 .or. lmhpc) then
!           --- DANS LE CAS DES LAGRANGE EN HPC, CERTAINES CONTRIBUTIONS
!           --- AU SECOND MEMBRE PEUVENT VENIR D'UN AUTRE PROCESSEUR
!           --- CAR LES MAILLES TARDIVES NE SONT PAS PARTAGEES DONC
!           --- CALCULEES SEULEMENT SUR UN PROC
                if (lmhpc) then
                    if (.not. ltypr) then
                        call filter_rhs_c(csolu, nonu//".NUME")
                    else
                        call filter_rhs(rsolu, nonu//".NUME")
                    end if
                end if

!           --- MISE A L'ECHELLE DES LAGRANGES DANS LE SECOND MEMBRE
!           --- RANG 0 UNIQUEMENT
                if (ltypr) then
                    call mrconl('MULT', lmat, nbeql, 'R', rsolu, &
                                nbsol)
                else
                    call mcconl('MULT', lmat, nbeql, 'C', csolu, &
                                nbsol)
                end if
            end if
!
!           --- PRISE EN COMPTE DES CHARGES CINEMATIQUES :
            if (ierd .ne. 0) then
!              --- ON RAZ RSOLU SUR LES RANGS > 0 POUR NE CUMULER QUE
!              --- LA CONTRIBUTION DES CHARGES CINEMATIQUES EN DISTRIBUE
                if (ldist .and. rang .ne. 0) then
                    do i = 1, nnbsol
                        if (ltypr) then
                            rsolu(i) = 0.d0
                        else
                            csolu(i) = dcmplx(0.d0, 0.d0)
                        end if
                    end do
                end if
                call jeveuo(vcival, 'L', idvalc)
                call jelira(vcival, 'TYPE', cval=rouc)
                if (ltypr) then
                    ASSERT(rouc .eq. 'R')
                    do i = 1, nbsol
                        call csmbgg(lmat, rsolu(nbeql*(i-1)+1), zr(idvalc), [cbid], [cbid], &
                                    'R')
                    end do
                else
                    ASSERT(rouc .eq. 'C')
                    do i = 1, nbsol
                        call csmbgg(lmat, [0.d0], [0.d0], csolu(nbeql*(i-1)+1), zc(idvalc), &
                                    'C')
                    end do
                end if
!
!         --- REDUCTION DU SECOND MEMBRE AU PROC MAITRE EN DISTRIBUE
!         --- POUR ETRE COHERENT AVEC LA MATRICE QUI CONTIENT DES N
!         --- SUR LA DIAGONALE
                if (ldist) then
                    if (ltypr) then
                        call asmpi_comm_vect('REDUCE', 'R', nbval=int(nnbsol, 8), vr=rsolu)
                    else
                        call asmpi_comm_vect('REDUCE', 'C', nbval=int(nnbsol, 8), vc=csolu)
                    end if
                end if
!
            end if
!
        end if
!
        if (lmhpc) then
            ASSERT(nbsol .eq. 1)
            if (ltypr) then
                call wkvect('&&AMUMPP.RHS', 'V V R', int(nnbsol, 8), jrhs)
            else
                call wkvect('&&AMUMPP.RHS', 'V V C', int(nnbsol, 8), jrhs)
            end if
            if (.not. lpreco) then
                do j = 1, nbeql
                    if (pddl(j) .eq. rang) then
                        nuglo = nulg(j)
                        if (ltypr) then
                            if (l_debug) then
                                nuno = zi(jdeeq+2*(j-1))
                                if (nuno .gt. 0) nuno = zi(jmlogl+nuno-1)+1
                                nucmp = zi(jdeeq+2*(j-1)+1)
!                    numéro noeud global, num comp du noeud, rhs
                                write (101+rang, *) nuno, nucmp, nuglo, rsolu(j)
                            end if
                            zr(jrhs+nuglo) = rsolu(j)
                        else
                            zc(jrhs+nuglo) = csolu(j)
                        end if
                    end if
                    if (l_debug) then
                        nuno = zi(jdeeq+2*(j-1))
                        if (nuno .gt. 0) nuno = zi(jmlogl+nuno-1)+1
                        nucmp = zi(jdeeq+2*(j-1)+1)
!                    numéro noeud global, num comp du noeud, rhs
                        write (201+rang, *) nuno, nucmp, rsolu(j), pddl(j), j, nulg(j)
                    end if
                end do
                if (l_debug) flush (101+rang)
                if (l_debug) flush (201+rang)
                if (ltypr) then
                    call asmpi_comm_vect('REDUCE', 'R', nbval=int(nnbsol, 8), vr=zr(jrhs))
                    call jgetptc(jrhs, pteur_c, vr=zr(1))
                    call c_f_pointer(pteur_c, rsolu2, [nnbsol])
                else
                    call asmpi_comm_vect('REDUCE', 'C', nbval=int(nnbsol, 8), vc=zc(jrhs))
                    call jgetptc(jrhs, pteur_c, vc=zc(1))
                    call c_f_pointer(pteur_c, csolu2, [nnbsol])
                end if
            else
!               if mumps is used as a preconditionner in HPC mode (asterxx),
!               each process knows the RHS and the numberings are the same
                do j = 1, nnbsol
                    if (ltypr) then
                        zr(jrhs+j-1) = rsolu(j)
                    else
                        zc(jrhs+j-1) = csolu(j)
                    end if
                end do
                if (ltypr) then
                    call jgetptc(jrhs, pteur_c, vr=zr(1))
                    call c_f_pointer(pteur_c, rsolu2, [nnbsol])
                else
                    call jgetptc(jrhs, pteur_c, vc=zc(1))
                    call c_f_pointer(pteur_c, csolu2, [nnbsol])
                end if
            end if
            call jedetr("&&AMUMPP.LAGR")
        else
            if (ltypr) then
                if (l_debug) then
                    do j = 1, nnbsol
                        nuno = zi(jdeeq+2*(j-1))
                        nucmp = zi(jdeeq+2*(j-1)+1)
!                numéro noeud, num comp du noeud, solution
                        write (49, *) nuno, nucmp, j, rsolu(j)
                    end do
                    flush (49)
                end if
                call jgetptc(1_8, pteur_c, vr=rsolu(1))
                call c_f_pointer(pteur_c, rsolu2, [nnbsol])
            else
                call jgetptc(1_8, pteur_c, vc=csolu(1))
                call c_f_pointer(pteur_c, csolu2, [nnbsol])
            end if
        end if
!
!        --- COPIE DE RSOLU DANS %RHS:
        if (rang .eq. 0) then
            if (type .eq. 'S') then
                do i = 1, nnbsol
                    raux = rsolu2(i)
                    rtest = abs(raux)
                    if (rtest .lt. rmin) then
                        raux = 0.d0
                    else if (rtest .gt. rmax) then
                        raux = rmax*sign(1.d0, raux)
                    end if
                    smpsk%rhs(i) = real(raux, kind=4)
                end do
            else if (type .eq. 'C') then
                do i = 1, nnbsol
                    caux = csolu2(i)
                    rtest = abs(caux)
                    if (rtest .lt. rmin) then
                        caux = dcmplx(0.d0, 0.d0)
                    else if (rtest .gt. rmax) then
                        caux = dcmplx(rmax*sign(1.d0, dble(caux)), 0.d0)
                        caux = rmax*dcmplx( &
                               1.d0*sign(1.d0, dble(caux)), 1.d0*sign(1.d0, imag(caux)))
                    end if
                    cmpsk%rhs(i) = cmplx(caux, kind=4)
                end do
            else if (type .eq. 'D') then
                do i = 1, nnbsol
                    raux = rsolu2(i)
                    rtest = abs(raux)
                    if (rtest .lt. rmin) then
                        raux = 0.d0
                    else if (rtest .gt. rmax) then
                        valr(1) = rtest
                        valr(2) = rmax
                        call utmess('F', 'FACTOR_79', si=i, nr=2_8, valr=valr)
                    end if
                    dmpsk%rhs(i) = raux
                end do
            else if (type .eq. 'Z') then
                do i = 1, nnbsol
                    caux = csolu2(i)
                    rtest = abs(caux)
                    if (rtest .lt. rmin) then
                        caux = dcmplx(0.d0, 0.d0)
                    else if (rtest .gt. rmax) then
                        valr(1) = rtest
                        valr(2) = rmax
                        call utmess('F', 'FACTOR_79', si=i, nr=2_8, valr=valr)
                    end if
                    zmpsk%rhs(i) = caux
                end do
            else
                ASSERT(.false.)
            end if
        end if
        if (lmhpc) call jedetr('&&AMUMPP.RHS')
!
!         -- IMPRESSION DU/DES SECONDS MEMBRES (SI DEMANDE) :
        if (impr(1:3) .eq. 'OUI') then
            if (rang .eq. 0) then
                if (type .eq. 'S') then
                    do k = 1, nnbsol
                        write (ifmump, *) k, smpsk%rhs(k)
                    end do
                else if (type .eq. 'C') then
                    do k = 1, nnbsol
                        write (ifmump, *) k, cmpsk%rhs(k)
                    end do
                else if (type .eq. 'D') then
                    do k = 1, nnbsol
                        write (ifmump, *) k, dmpsk%rhs(k)
                    end do
                else if (type .eq. 'Z') then
                    do k = 1, nnbsol
                        write (ifmump, *) k, zmpsk%rhs(k)
                    end do
                else
                    ASSERT(.false.)
                end if
                write (ifmump, *) 'MUMPS FIN RHS'
            end if
            if (impr(1:11) .eq. 'OUI_NOSOLVE') then
                call utmess('F', 'FACTOR_71', si=ifmump)
            end if
        end if
!
!
    else if (option .eq. 2) then
!
        if (lmhpc) then
            if (ltypr) then
                call wkvect('&&AMUMPP.RHS', 'V V R', int(nnbsol, 8), jrhs)
                call jgetptc(jrhs, pteur_c, vr=zr(1))
                call c_f_pointer(pteur_c, rsolu2, [nnbsol])
            else
                call wkvect('&&AMUMPP.RHS', 'V V C', int(nnbsol, 8), jrhs)
                call jgetptc(jrhs, pteur_c, vc=zc(1))
                call c_f_pointer(pteur_c, csolu2, [nnbsol])
            end if
        else
            if (ltypr) then
                call jgetptc(1_8, pteur_c, vr=rsolu(1))
                call c_f_pointer(pteur_c, rsolu2, [nnbsol])
            else
                call jgetptc(1_8, pteur_c, vc=csolu(1))
                call c_f_pointer(pteur_c, csolu2, [nnbsol])
            end if
        end if
!
!       ------------------------------------------------
!        POST-TRAITEMENTS ASTER DE LA SOLUTION :
!       ------------------------------------------------
        if (rang .eq. 0 .or. lmhpc) then
            if (rang .eq. 0) then
                if (type .eq. 'S') then
                    do i = 1, nnbsol
                        rsolu2(i) = smpsk%rhs(i)
                    end do
                    deallocate (smpsk%rhs)
                else if (type .eq. 'C') then
                    do i = 1, nnbsol
                        csolu2(i) = cmpsk%rhs(i)
                    end do
                    deallocate (cmpsk%rhs)
                else if (type .eq. 'D') then
                    b_n = to_blas_int(nnbsol)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, dmpsk%rhs, b_incx, rsolu2, b_incy)
                    deallocate (dmpsk%rhs)
                else if (type .eq. 'Z') then
                    b_n = to_blas_int(nnbsol)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call zcopy(b_n, zmpsk%rhs, b_incx, csolu2, b_incy)
                    deallocate (zmpsk%rhs)
                else
                    ASSERT(.false.)
                end if
            end if
!
            if (lmhpc) then
                if (ltypr) then
                    call asmpi_comm_vect('BCAST', 'R', nbval=int(nnbsol, 8), bcrank=0_8, vr=rsolu2)
                    if (l_debug) then
                        do j = 1, nbeql
                            nuno = zi(jdeeq+2*(j-1))
                            if (nuno .gt. 0) nuno = zi(jmlogl+nuno-1)+1
                            nucmp = zi(jdeeq+2*(j-1)+1)
!                    numero ddl local, numéro noeud local, numéro noeud global, num comp du noeud,
!                                num ddl global, num proc proprio, solution
! write(51+rang,*) j, zi(jdeeq+2*(j-1)), nuno, nucmp,  &
!                      nulg(j), pddl(j), rsolu(j)
                            write (51+rang, *) nuno, nucmp, rsolu2(j)
                        end do
                        flush (51+rang)
                    end if
                    if (.not. lpreco) then
                        do j = 1, nbeql
                            rsolu(j) = rsolu2(nulg(j)+1)
                        end do
                    else
!               if mumps is used as a preconditionner in HPC mode (asterxx),
!               each process knows the SOL and the numberings are the same
                        do j = 1, nnbsol
                            rsolu(j) = rsolu2(j)
                        end do
                    end if
                else
                    call asmpi_comm_vect('BCAST', 'C', nbval=int(nnbsol, 8), bcrank=0_8, vc=csolu2)
                    if (.not. lpreco) then
                        do j = 1, nbeql
                            csolu(j) = csolu2(nulg(j)+1)
                        end do
                    else
                        do j = 1, nnbsol
                            csolu(j) = csolu2(j)
                        end do
                    end if
                end if
            end if
!
            if (eli2lg) then
!           -- PRISE EN COMPTE DES LAGRANGES "2" :
!           -- EN SORTIE DE RESOLUTION AVEC ELIM_LAGR='LAGR2' ON A :
!           -- LAGR1 = LAGR1 + LAGR2, ON DOIT RECTIFIER CELA :
!           -- LAGR1 = LAGR1/2 PUIS LAGR2 = LAGR1
!           -- VALIDE QUE SUR PROC 0, MAIS C'EST OK CAR ON
!           -- BROADCAST LA SOLUTION APRES
                call jeveuo(nonu//'.NUME.DELG', 'L', vi=delg)
                call nudlg2(nonu)
                call jeveuo(nonu//'.NUME.DLG2', 'L', vi=dlg2)
                if (ltypr) then
                    do i = 1, nbsol
                        do k = 1, nbeql
                            if (delg(k) .eq. -1) then
                                rsolu((i-1)*nbeql+k) = 0.5d0*rsolu((i-1)*nbeql+k)
                                jj = dlg2(k)
                                if (jj .ne. 0) then
                                    rsolu((i-1)*nbeql+jj) = rsolu((i-1)*nbeql+k)
                                end if
                            end if
                        end do
                    end do
                else
                    do i = 1, nbsol
                        do k = 1, nbeql
                            if (delg(k) .eq. -1) then
                                csolu((i-1)*nbeql+k) = 0.5d0*csolu((i-1)*nbeql+k)
                                jj = dlg2(k)
                                if (jj .ne. 0) then
                                    csolu((i-1)*nbeql+jj) = csolu((i-1)*nbeql+k)
                                end if
                            end if
                        end do
                    end do
                end if
            end if
!
!         --- MISE A L'ECHELLE DES LAGRANGES DANS LA SOLUTION :
!         ON NE LE FAIT PAS SI NON DEMANDE
            if (.not. lpreco .and. prepos) then
                if (ltypr) then
                    call mrconl('MULT', lmat, nbeql, 'R', rsolu, &
                                nbsol)
                else
                    call mcconl('MULT', lmat, nbeql, 'C', csolu, &
                                nbsol)
                end if
            end if
        end if
!
!       -- BROADCAST DE SOLU A TOUS LES PROC
        if (.not. lmhpc) then
            if (ltypr) then
                call asmpi_comm_vect('BCAST', 'R', nbval=int(nnbsol, 8), bcrank=0_8, vr=rsolu)
            else
                call asmpi_comm_vect('BCAST', 'C', nbval=int(nnbsol, 8), bcrank=0_8, vc=csolu)
            end if
        end if
!
!       -- IMPRESSION DU/DES SOLUTIONS (SI DEMANDE) :
        if ((lverif) .and. (rang .eq. 0)) then
            raux = 0.d0
            if (ltypr) then
                do k = 1, nnbsol
                    raux = raux+abs(rsolu2(k))
                end do
            else
                do k = 1, nnbsol
                    raux = raux+abs(csolu2(k))
                end do
            end if
            write (ifm, *) 'NORME L1 MUMPS SOLUTION=', raux
        end if
        if (impr(1:9) .eq. 'OUI_SOLVE') then
            if (rang .eq. 0) then
                if (ltypr) then
                    do k = 1, nnbsol
                        write (ifmump, *) k, rsolu2(k)
                    end do
                else
                    do k = 1, nnbsol
                        write (ifmump, *) k, csolu2(k)
                    end do
                end if
                write (ifmump, *) 'MUMPS FIN SOLUTION'
            end if
        end if
!
        if (l_debug .and. .not. lmhpc) then
            do j = 1, nnbsol
                nuno = zi(jdeeq+2*(j-1))
                nucmp = zi(jdeeq+2*(j-1)+1)
!                numéro noeud, num comp du noeud, solution
                write (50, *) nuno, nucmp, rsolu(j)
            end do
            flush (50)
        end if
        if (lmhpc) call jedetr('&&AMUMPP.RHS')
    else
!       ------------------------------------------------
!        MAUVAISE OPTION
!       ------------------------------------------------
        ASSERT(.false.)
    end if
    call jedetr(nonu//'.NUME.DEE2')
    call jedema()
#endif
end subroutine
