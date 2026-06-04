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

subroutine sscgma(ma, nbgmp, nbgmin)
!
    implicit none
!
#include "asterc/getres.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cgmaap.h"
#include "asterfort/cgmaba.h"
#include "asterfort/cgmacy.h"
#include "asterfort/cgmafn.h"
#include "asterfort/cgmasp.h"
#include "asterfort/cgmaxf.h"
#include "asterfort/cgmftm.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jeexin.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utlisi.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/addGrpMa.h"
#include "jeveux.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: ma
    integer(kind=8), intent(in)          :: nbgmp
    integer(kind=8), intent(in)          :: nbgmin

!     BUT: TRAITER LE MOT CLEF CREA_GROUP_MA
!          DE L'OPERATEUR: DEFI_GROUP
!
!     IN:
!          MA    : NOM DU MAILLAGE
!          NBGMP : NOMBRE DE GROUP_MA A CREER
!     ------------------------------------------------------------------
!
    character(len=8) :: noma, kbid, kpos, nom1
    character(len=8) :: alarm, tyma
    character(len=16) :: concep, cmd, option
    character(len=24) :: lisma, nogma, nogma2
    character(len=24) :: valk(2)
    character(len=132) :: card
    integer(kind=8) :: i, iagm1, iagm2, ialii1, ialii2, iret
    integer(kind=8) :: idlima, ier, ierr, ifm, igm, igm1
    integer(kind=8) :: igm2, ii, iii, ili1, ili2, im1
    integer(kind=8) :: ima, ind1, ind2, iocc, ireste, jjj
    integer(kind=8) :: jlisma, jmail, kkk, maxcol, n, n1
    integer(kind=8) :: n2, n3, n4, n5, n6, n6a, n6b
    integer(kind=8) :: n7, n8, nalar, nb, nbcol
    integer(kind=8) :: nbgnaj, nbgrmn, nbid, nbis, nbk8, nbline, nbma
    integer(kind=8) :: nbmat, niv, ntrou, ntyp, num, nverif
    aster_logical :: l_parallel_mesh, l_added_grpma, lcolle
    character(len=24), pointer :: lik8(:) => null()
    character(len=8), pointer :: l_maille(:) => null()
    integer(kind=8), pointer :: maille2(:) => null()
    character(len=24), pointer :: lik8_2(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
!     RECUPERATION DU NIVEAU D'IMPRESSION
!     -----------------------------------
    call infniv(ifm, niv)
!
    call getres(kbid, concep, cmd)
    l_parallel_mesh = isParallelMesh(ma)
    lisma = '&&SSCGMA.LISTE_MAILLES'
    if (l_parallel_mesh) then
        call jelira(ma//'.PAR_GRPMAI', 'NOMMAX', nbgrmn)
    else
        call jelira(ma//'.GROUPEMA', 'NMAXOC', nbgrmn)
    end if
    nbis = nbgrmn
    nbk8 = nbgrmn
    call wkvect('&&SSCGMA.LII1', 'V V I', nbis, ialii1)
    call wkvect('&&SSCGMA.LII2', 'V V I', nbis, ialii2)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbmat)
!
    call getvtx(' ', 'ALARME', scal=alarm, nbret=nalar)
!
    nbgnaj = 0
    do iocc = 1, nbgmp
!
        call getvtx('CREA_GROUP_MA', 'NOM', iocc=iocc, scal=nogma, nbret=n1)
!
        call jenonu(jexnom(ma//'.GROUPEMA', nogma), iret)
        if (iret .gt. 0) then
            call utmess('F', 'ALGELINE3_7', sk=nogma)
        end if
!
        call getvem(ma, 'MAILLE', 'CREA_GROUP_MA', 'MAILLE', iocc, &
                    0, kbid, n2)
        call getvtx('CREA_GROUP_MA', 'INTERSEC', iocc=iocc, nbval=0, nbret=n3)
        call getvtx('CREA_GROUP_MA', 'UNION', iocc=iocc, nbval=0, nbret=n4)
        call getvtx('CREA_GROUP_MA', 'DIFFE   ', iocc=iocc, nbval=0, nbret=n5)
        call getvem(ma, 'GROUP_MA', 'CREA_GROUP_MA', 'GROUP_MA', iocc, &
                    0, kbid, n6)
        call getvtx('CREA_GROUP_MA', 'OPTION', iocc=iocc, nbval=0, nbret=n7)
        call getvtx('CREA_GROUP_MA', 'TOUT', iocc=iocc, nbval=0, nbret=n8)
        n2 = -n2
        n3 = -n3
        n4 = -n4
        n5 = -n5
        n6 = -n6
        n7 = -n7
        n8 = -n8
        nbma = 0

!
!
!       -- MOT CLEF TOUT:
!       -------------------
        if (n8 .gt. 0) then
            nbma = nbmat
            call wkvect(lisma, 'V V I', nbma, jlisma)
            do ima = 1, nbmat
                zi(jlisma-1+ima) = ima
            end do
            goto 219
        end if
!
!
!       -- MOT CLEF MAILLE:
!       -------------------
        if (n2 .gt. 0) then
            if (l_parallel_mesh) then
                call utmess('F', 'MODELISA7_86')
            end if
            AS_ALLOCATE(vk8=l_maille, size=n2)
            call getvem(ma, 'MAILLE', 'CREA_GROUP_MA', 'MAILLE', iocc, &
                        n2, l_maille, n1)
            call wkvect('&&SSCGMA.MAILLE', 'V V I', n2, jmail)
            call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbmat)
            AS_ALLOCATE(vi=maille2, size=nbmat)
            nbma = 0
            lcolle = .false.
            call jeexin(ma//'.NOMMAI', ier)
            if (ier .ne. 0) then
                lcolle = .true.
            end if
            ier = 0
            do im1 = 1, n2
                nom1 = l_maille(im1)
                num = char8_to_int(nom1, lcolle, ma, "MAILLE")
                if (num .eq. 0) then
                    ier = ier+1
                    call utmess('E', 'SOUSTRUC_31', sk=nom1)
                    goto 20
                end if
                maille2(num) = maille2(num)+1
                if (maille2(num) .eq. 2) then
                    valk(1) = nom1
                    valk(2) = nogma
                    call utmess('A', 'SOUSTRUC_32', nk=2, valk=valk)
                    goto 20
                end if
                nbma = nbma+1
                zi(jmail+nbma-1) = num
20              continue
            end do
            if (ier .ne. 0) then
                ASSERT(.false.)
            end if
            call wkvect(lisma, 'V V I', nbma, jlisma)
            do ima = 0, nbma-1
                zi(jlisma+ima) = zi(jmail+ima)
            end do
            call jedetr('&&SSCGMA.MAILLE')
            AS_DEALLOCATE(vi=maille2)
            AS_DEALLOCATE(vk8=l_maille)
            goto 219
        end if
!
!
!       -- MOT CLEF GROUP_MA:
!       ---------------------
        if (n6 .gt. 0) then
            call getvem(ma, 'GROUP_MA', 'CREA_GROUP_MA', 'GROUP_MA', iocc, &
                        1, nogma2, nbid)
            call getvtx('CREA_GROUP_MA', 'POSITION', iocc=iocc, nbval=0, nbret=n6b)
            if (nbid .ne. 0) then
                call jenonu(jexnom(ma//'.GROUPEMA', nogma2), igm2)
                call jelira(jexnum(ma//'.GROUPEMA', igm2), 'LONUTI', ili2)
                call jeveuo(jexnum(ma//'.GROUPEMA', igm2), 'L', iagm2)
                ind1 = 0
                ind2 = 0
                if (n6b .eq. 0) then
                    call getvis('CREA_GROUP_MA', 'NUME_INIT', iocc=iocc, scal=ind1, nbret=n6a)
                    if (n6a .eq. 0) ind1 = 1
                    call getvis('CREA_GROUP_MA', 'NUME_FIN', iocc=iocc, scal=ind2, nbret=n6a)
                    if (n6a .eq. 0) ind2 = ili2
                    if (ind2 .lt. ind1) then
                        call utmess('F', 'SOUSTRUC_33')
                    end if
                    if (ili2 .lt. ind2) then
                        call utmess('F', 'SOUSTRUC_34')
                    end if
                    n6a = ind2-ind1+1
                else
                    n6a = 1
                end if
                call wkvect(lisma, 'V V I', n6a, jlisma)
                nbma = n6a
                if (n6b .eq. 0) then
                    n = ind2-ind1+1
                    do ii = 1, n
                        zi(jlisma-1+ii) = zi(iagm2-2+ind1+ii)
                    end do
                    goto 219
                end if
                call getvtx('CREA_GROUP_MA', 'POSITION', iocc=iocc, scal=kpos, nbret=n6b)
                if (kpos .eq. 'INIT') then
                    zi(jlisma) = zi(iagm2)
                else if (kpos .eq. 'FIN') then
                    ii = ili2
                    zi(jlisma) = zi(iagm2+ii-1)
                else if (kpos .eq. 'MILIEU') then
                    ii = (ili2+1)/2
                    zi(jlisma) = zi(iagm2+ii-1)
                end if
            end if
            goto 219
        end if
!
!
!       -- MOT CLEF INTER:
!       -------------------
        if (n3 .gt. 0) then
            AS_ALLOCATE(vk24=lik8, size=n3)
            call getvem(ma, 'GROUP_MA', 'CREA_GROUP_MA', 'INTERSEC', iocc, &
                        n3, lik8, nbid)
            n3 = nbid
            do igm = 1, n3
                call jenonu(jexnom(ma//'.GROUPEMA', lik8(igm)), igm2)
                if (igm2 .eq. 0) then
                    call utmess('F', 'SOUSTRUC_35', sk=lik8(igm))
                end if
            end do
!
            if (n3 .ne. 0) then
                call jenonu(jexnom(ma//'.GROUPEMA', lik8(1)), igm1)
                call jelira(jexnum(ma//'.GROUPEMA', igm1), 'LONUTI', ili1)
                call jeveuo(jexnum(ma//'.GROUPEMA', igm1), 'L', iagm1)
                if (ili1 .gt. nbis) then
                    nbis = 2*ili1
                    call jedetr('&&SSCGMA.LII1')
                    call jedetr('&&SSCGMA.LII2')
                    call wkvect('&&SSCGMA.LII1', 'V V I', nbis, ialii1)
                    call wkvect('&&SSCGMA.LII2', 'V V I', nbis, ialii2)
                end if
                n = ili1
                do ii = 1, n
                    zi(ialii1-1+ii) = zi(iagm1-1+ii)
                end do
!
                do igm = 2, n3
                    call jenonu(jexnom(ma//'.GROUPEMA', lik8(igm)), igm2)
                    call jelira(jexnum(ma//'.GROUPEMA', igm2), 'LONUTI', ili2)
                    call jeveuo(jexnum(ma//'.GROUPEMA', igm2), 'L', iagm2)
                    call utlisi('INTER', zi(ialii1), n, zi(iagm2), ili2, &
                                zi(ialii2), nbis, ntrou)
                    n = ntrou
                    do ii = 1, n
                        zi(ialii1-1+ii) = zi(ialii2-1+ii)
                    end do
                end do
                AS_DEALLOCATE(vk24=lik8)
!
                if (n .eq. 0) then
                    if (alarm .eq. 'OUI') then
                        call utmess('A', 'SOUSTRUC_36', sk=nogma)
                    end if
                else
                    call wkvect(lisma, 'V V I', n, jlisma)
                    nbma = n
                    do ii = 1, n
                        zi(jlisma-1+ii) = zi(ialii1-1+ii)
                    end do
                end if
            else
                AS_DEALLOCATE(vk24=lik8)
            end if
            goto 219
        end if
!
!
!       -- MOT CLEF UNION:
!       -------------------
        if (n4 .gt. 0) then
            AS_ALLOCATE(vk24=lik8, size=n4)
            call getvem(ma, 'GROUP_MA', 'CREA_GROUP_MA', 'UNION', iocc, &
                        n4, lik8, nbid)
            n4 = nbid
            do igm = 1, n4
                call jenonu(jexnom(ma//'.GROUPEMA', lik8(igm)), igm2)
                if (igm2 .eq. 0) then
                    call utmess('F', 'SOUSTRUC_35', sk=lik8(igm))
                end if
            end do
!
            if (n4 .ne. 0) then
                call jenonu(jexnom(ma//'.GROUPEMA', lik8(1)), igm1)
                call jelira(jexnum(ma//'.GROUPEMA', igm1), 'LONUTI', ili1)
                call jeveuo(jexnum(ma//'.GROUPEMA', igm1), 'L', iagm1)
                if (ili1 .gt. nbis) then
                    nbis = 2*ili1
                    call jedetr('&&SSCGMA.LII1')
                    call jedetr('&&SSCGMA.LII2')
                    call wkvect('&&SSCGMA.LII1', 'V V I', nbis, ialii1)
                    call wkvect('&&SSCGMA.LII2', 'V V I', nbis, ialii2)
                end if
                n = ili1
                do ii = 1, n
                    zi(ialii1-1+ii) = zi(iagm1-1+ii)
                end do
!
                do igm = 2, n4
                    call jenonu(jexnom(ma//'.GROUPEMA', lik8(igm)), igm2)
                    call jelira(jexnum(ma//'.GROUPEMA', igm2), 'LONUTI', ili2)
                    call jeveuo(jexnum(ma//'.GROUPEMA', igm2), 'L', iagm2)
                    call utlisi('UNION', zi(ialii1), n, zi(iagm2), ili2, &
                                zi(ialii2), nbis, ntrou)
!
                    if (ntrou .lt. 0) then
                        nbis = -2*ntrou
                        call jedetr('&&SSCGMA.LII2')
                        call wkvect('&&SSCGMA.LII2', 'V V I', nbis, ialii2)
                        call utlisi('UNION', zi(ialii1), n, zi(iagm2), ili2, &
                                    zi(ialii2), nbis, ntrou)
                        call jedetr('&&SSCGMA.LII1')
                        call wkvect('&&SSCGMA.LII1', 'V V I', nbis, ialii1)
                    end if
                    n = ntrou
                    do ii = 1, n
                        zi(ialii1-1+ii) = zi(ialii2-1+ii)
                    end do
                end do
                AS_DEALLOCATE(vk24=lik8)
!
                if (n .eq. 0) then
                    if (alarm .eq. 'OUI') then
                        call utmess('A', 'SOUSTRUC_36', sk=nogma)
                    end if
                else
                    call wkvect(lisma, 'V V I', n, jlisma)
                    nbma = n
                    do ii = 1, n
                        zi(jlisma-1+ii) = zi(ialii1-1+ii)
                    end do
                end if
            else
                AS_DEALLOCATE(vk24=lik8)
            end if
            goto 219
        end if
!
!
!       -- MOT CLEF DIFFE:
!       -------------------
        if (n5 .gt. 0) then
            AS_ALLOCATE(vk24=lik8_2, size=n5)
            call getvtx('CREA_GROUP_MA', 'DIFFE', iocc=iocc, nbval=n5, vect=lik8_2, nbret=nverif)
            AS_ALLOCATE(vk24=lik8, size=n5)
            call getvem(ma, 'GROUP_MA', 'CREA_GROUP_MA', 'DIFFE', iocc, &
                        n5, lik8, nbid)
            n5 = nbid
            if (lik8_2(1) .ne. lik8(1)) then
                n5 = 0
            end if
            if (nverif .ne. 0) then
                AS_DEALLOCATE(vk24=lik8_2)
            end if
            do igm = 1, n5
                call jenonu(jexnom(ma//'.GROUPEMA', lik8(igm)), igm2)
                if (igm2 .eq. 0) then
                    call utmess('F', 'SOUSTRUC_35', sk=lik8(igm))
                end if
            end do
!
            if (n5 .ne. 0) then
                call jenonu(jexnom(ma//'.GROUPEMA', lik8(1)), igm1)
                call jelira(jexnum(ma//'.GROUPEMA', igm1), 'LONUTI', ili1)
                call jeveuo(jexnum(ma//'.GROUPEMA', igm1), 'L', iagm1)
                if (ili1 .gt. nbis) then
                    nbis = 2*ili1
                    call jedetr('&&SSCGMA.LII1')
                    call jedetr('&&SSCGMA.LII2')
                    call wkvect('&&SSCGMA.LII1', 'V V I', nbis, ialii1)
                    call wkvect('&&SSCGMA.LII2', 'V V I', nbis, ialii2)
                end if
                n = ili1
                do ii = 1, n
                    zi(ialii1-1+ii) = zi(iagm1-1+ii)
                end do
!
                do igm = 2, n5
                    call jenonu(jexnom(ma//'.GROUPEMA', lik8(igm)), igm2)
                    call jelira(jexnum(ma//'.GROUPEMA', igm2), 'LONUTI', ili2)
                    call jeveuo(jexnum(ma//'.GROUPEMA', igm2), 'L', iagm2)
                    call utlisi('DIFFE', zi(ialii1), n, zi(iagm2), ili2, &
                                zi(ialii2), nbis, ntrou)
                    n = ntrou
                    do ii = 1, n
                        zi(ialii1-1+ii) = zi(ialii2-1+ii)
                    end do
                end do
                AS_DEALLOCATE(vk24=lik8)
!
                if (n .eq. 0) then
                    if (alarm .eq. 'OUI') then
                        call utmess('A', 'SOUSTRUC_36', sk=nogma)
                    end if
                else
                    call wkvect(lisma, 'V V I', n, jlisma)
                    nbma = n
                    do ii = 1, n
                        zi(jlisma-1+ii) = zi(ialii1-1+ii)
                    end do
                end if
            else
                AS_DEALLOCATE(vk24=lik8)
            end if
            goto 219
        end if
!
!
!       -- MOT CLEF OPTION:
!       -------------------
        if (n7 .gt. 0) then
!
            call getvtx('CREA_GROUP_MA', 'OPTION', iocc=iocc, scal=option, nbret=nb)
!
!            -- TRAITEMENT DE L'OPTION FACE_NORMALE :
!               -----------------------------------
            if (option(1:12) .eq. 'FACE_NORMALE') then
                call cgmafn('CREA_GROUP_MA', iocc, ma, lisma, nbma)
!
!            -- TRAITEMENT DE L'OPTION SPHERE :
!               -----------------------------
            else if (option(1:6) .eq. 'SPHERE') then
                call cgmasp('CREA_GROUP_MA', iocc, ma, lisma, nbma)
!
!            -- TRAITEMENT DE L'OPTION CYLINDRE :
!               -------------------------------
            else if (option(1:8) .eq. 'CYLINDRE') then
                call cgmacy('CREA_GROUP_MA', iocc, ma, lisma, nbma)
!
!            -- TRAITEMENT DE L'OPTION BANDE :
!               ----------------------------
            else if (option(1:5) .eq. 'BANDE') then
                call cgmaba('CREA_GROUP_MA', iocc, ma, lisma, nbma)
!
!            -- TRAITEMENT DE L'OPTION APPUI_STRICT :
!               ----------------------------------
            else if (option(1:5) .eq. 'APPUI') then
                call cgmaap('CREA_GROUP_MA', iocc, ma, lisma, nbma)
!
!            -- TRAITEMENT DE L'OPTION FISS_XFEM :
!               ----------------------------------
            else if (option(1:9) .eq. 'FISS_XFEM') then
                call cgmaxf('CREA_GROUP_MA', iocc, ma, lisma, nbma)
            end if
        end if
!
!
!       -- ON FILTRE LES TYPES DE MAILLES :
!       -----------------------------------
219     continue
!
        if (nbma .gt. 0) then
            call getvtx('CREA_GROUP_MA', 'TYPE_MAILLE', iocc=iocc, scal=tyma, nbret=ntyp)
            if (tyma(1:4) .ne. 'TOUT') then
                call cgmftm(tyma, ma, lisma, nbma, ierr)
                if (ierr .ne. 0) then
                    if (alarm .eq. 'OUI') then
                        call utmess('A', 'SOUSTRUC_36', sk=nogma)
                    end if
                end if
            end if
        end if

!
!       -- CREATION ET AFFECTATION DU GROUP_MA :
!       ----------------------------------
        call jeexin(lisma, iret)
        if (iret .ne. 0) then
            call jeveuo(lisma, 'L', idlima)
            call addGrpMa(ma, nogma, zi(idlima), nbma, l_added_grpma)

            if (l_added_grpma) then
                nbgnaj = nbgnaj+1
            end if
        end if
!
        call jedetr(lisma)
!
    end do
!
!     IMPRESSIONS NIVEAUX 1 ET 2
!     --------------------------
    if (niv .ge. 1 .and. nbgnaj .ne. 0) then
        write (ifm, '(/,/,A,I6,/,45(''=''))')&
     &    'NOMBRE  DE GROUPES DE MAILLES CREES : ', nbgnaj
!
        write (ifm, '(/,15X,54(''-''),2(/,15X,A),/,15X,54(''-''))')&
     &    '!         NOM DU GROUPE         ! NBRE DE MAILLES DU !',&
     &    '!            MAILLES            !     GROUPE_MA      !'
!
        do i = 1, nbgnaj
            ii = nbgmin+i
            call jenuno(jexnum(ma//'.GROUPEMA', ii), nogma)
            call jelira(jexnum(ma//'.GROUPEMA', ii), 'LONUTI', nbma)
            write (ifm, '(15X,A,2X,A24,5X,A,2X,I8,10X,A)') '!', nogma, '!',&
     &      nbma, '!'
        end do
        write (ifm, '(15X,54(''-''),/)')
    end if
!
!     IMPRESSIONS NIVEAU 2
!     --------------------
    if (niv .eq. 2 .and. nbgnaj .ne. 0) then
        maxcol = 8
        do i = 1, nbgnaj
            ii = nbgmin+i
            call jeveuo(jexnum(ma//'.GROUPEMA', ii), 'L', jlisma)
            call jenuno(jexnum(ma//'.GROUPEMA', ii), nogma)
            call jelira(jexnum(ma//'.GROUPEMA', ii), 'LONUTI', nbma)
            write (ifm, '(/,3A,/,28(''-''))') 'MAILLES DU GROUPE ', &
                nogma, ' :'
            nbline = nbma/maxcol
            ireste = mod(nbma, maxcol)
            if (ireste .ne. 0) nbline = nbline+1
            nbcol = maxcol
            kkk = 0
            lcolle = .false.
            call jeexin(ma//'.NOMMAI', ier)
            if (ier .ne. 0) then
                lcolle = .true.
            end if
            do jjj = 1, nbline
                if (ireste .ne. 0 .and. jjj .eq. nbline) nbcol = ireste
                do iii = 1, nbcol
                    kkk = kkk+1
                    noma = int_to_char8(zi(jlisma-1+kkk), lcolle, ma, "MAILLE")
                    card((iii-1)*10+1:) = ' '//noma//' '
                end do
                write (ifm, '(A)') card(:10*nbcol)
            end do
        end do
        write (ifm, '(/,/)')
    end if
!
!
! --- MENAGE
    call jedetr(lisma)
    call jedetr('&&SSCGMA.LII1')
    call jedetr('&&SSCGMA.LII2')
!
    call jedema()
!
end subroutine
