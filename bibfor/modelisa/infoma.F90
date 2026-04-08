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
subroutine infoma(nomu, niv_)
    implicit none
#include "jeveux.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/dismoi.h"
!
    character(len=8) :: nomu
    integer(kind=8), optional :: niv_
!
!     IMPRESSION DES INFOS (1 OU 2)
!
!
    character(len=32) :: lisnoe, lismai, lisgrn, lisgrm
    character(len=32) :: comnoe, commai, comgrn, comgrm
    character(len=24) :: conxv, grpnov, grpmav, titre, cooval
    character(len=24) :: nom
    character(len=8) :: type
    integer(kind=8) :: niv, ifm, nn, nbno, j, idec, iad1, nbcoor, nbma
    integer(kind=8) :: nbltit, iad, i, nbnoeu, nbmail, nbgrno, nbgrma
    integer(kind=8) :: nbmmai, n1, nbmmax, ityp, nbnoeuL, nbmailL
    integer(kind=8) :: rang, nbproc
    parameter(nbmmax=100)
    integer(kind=8) :: dimmai(nbmmax), iret, iret2
    character(len=8) :: mclmai(nbmmax)
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    integer(kind=8), pointer :: noex(:) => null()
    integer(kind=8), pointer :: maex(:) => null()
    aster_logical :: lPMesh
    mpi_int :: mrank, msize
!
!
!
    data lisnoe/'LISTE DES NOEUDS                '/
    data lismai/'LISTE DES MAILLES               '/
    data lisgrn/'LISTE DES GROUPES DE NOEUDS     '/
    data lisgrm/'LISTE DES GROUPES DE MAILLES    '/
    data comnoe/'NOMBRE DE NOEUDS                '/
    data commai/'NOMBRE DE MAILLES               '/
    data comgrn/'NOMBRE DE GROUPES DE NOEUDS     '/
    data comgrm/'NOMBRE DE GROUPES DE MAILLES    '/
!
    call jemarq()
!
    call jeexin(nomu//'.NOEX', iret)
    call jeexin(nomu//'.MAEX', iret2)
    lPMesh = .false.
    if (iret .ne. 0 .and. iret2 .ne. 0) then
        lPMesh = .true.
    end if
!
    conxv = nomu//'.CONNEX'
    grpnov = nomu//'.GROUPENO'
    grpmav = nomu//'.GROUPEMA'
    titre = nomu//'           .TITR'
    cooval = nomu//'.COORDO    .VALE'
!
    call dismoi('NB_MA_MAILLA', nomu, 'MAILLAGE', repi=nbmail)
    call dismoi('NB_NO_MAILLA', nomu, 'MAILLAGE', repi=nbnoeu)
!
    if (lPMesh) then
        call jeveuo(nomu//'.NOEX', 'L', vi=noex)
        call jeveuo(nomu//'.MAEX', 'L', vi=maex)
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
        nbmailL = nbmail
        nbnoeuL = 0
        do i = 1, nbnoeu
            if (noex(i) .eq. rang) then
                nbnoeuL = nbnoeuL+1
            end if
        end do
        call asmpi_comm_vect('MPI_SUM', 'I', sci=nbnoeuL)
    else
        nbmailL = nbmail
        nbnoeuL = nbnoeu
    end if
!
    if (present(niv_)) then
        ifm = 6
        niv = niv_
    else
        call infniv(ifm, niv)
    end if
!
    call jeexin(grpmav, iret)
    if (iret .gt. 0) then
        call jelira(grpmav, 'NMAXOC', nbgrma)
    else
        nbgrma = 0
    end if
    call jeexin(grpnov, iret)
    if (iret .gt. 0) then
        call jelira(grpnov, 'NMAXOC', nbgrno)
    else
        nbgrno = 0
    end if
!
!
    call jeexin(titre, iret)
    if (iret .gt. 0) then
        call jelira(titre, 'LONMAX', nbltit)
    else
        nbltit = 0
    end if
    call jeveuo(nomu//'.DIME', 'L', vi=dime)
    call jeveuo(nomu//'.TYPMAIL', 'L', vi=typmail)
    nbcoor = dime(6)
!
!
    call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbmmai)
    do i = 1, nbmmai
        dimmai(i) = 0
        call jenuno(jexnum('&CATA.TM.NOMTM', i), mclmai(i))
    end do
    if (lPMesh) then
        nbmailL = 0
        do i = 1, nbmail
            ityp = typmail(i)
            ASSERT((ityp .gt. 0) .and. (ityp .lt. 100))
            if (maex(i) .eq. rang) then
                dimmai(ityp) = dimmai(ityp)+1
                nbmailL = nbmailL+1
            end if
        end do
        call asmpi_comm_vect('MPI_SUM', 'I', 100, vi=dimmai)
        call asmpi_comm_vect('MPI_SUM', 'I', sci=nbmailL)
    else
        do i = 1, nbmail
            ityp = typmail(i)
            ASSERT((ityp .gt. 0) .and. (ityp .lt. 100))
            dimmai(ityp) = dimmai(ityp)+1
        end do
    end if
!
!
!
!
! -     ECRITURE DE L EN TETE
! ----------------------------------
    if (niv .ge. 1) then
        write (ifm, 802) nomu, niv
        if (nbltit .gt. 0) then
            call jeveuo(titre, 'L', iad)
            do i = 1, nbltit
                write (ifm, 801) zk80(iad+i-1)
            end do
        end if
        write (ifm, 804) comnoe, nbnoeuL
        write (ifm, 804) commai, nbmailL
        do i = 1, nbmmai
            if (dimmai(i) .ne. 0) write (ifm, 806) mclmai(i), dimmai(i)
        end do
!
        if (niv .lt. 2) then
            if (nbgrno .ne. 0) then
                write (ifm, 804) comgrn, nbgrno
            end if
            !
            if (nbgrma .ne. 0) then
                write (ifm, 804) comgrm, nbgrma
            end if
        end if
    end if
!
!
    if (niv .ge. 2) then
!
        if (nbgrno .ne. 0) then
            write (ifm, 804) comgrn, nbgrno
            do i = 1, nbgrno
                call jeexin(jexnum(grpnov, i), iret)
                if (iret .eq. 0) goto 60
                call jenuno(jexnum(grpnov, i), nom)
                call jelira(jexnum(grpnov, i), 'LONUTI', n1)
                write (ifm, 808) nom, n1
60              continue
            end do
        end if
!
        if (nbgrma .ne. 0) then
            write (ifm, 804) comgrm, nbgrma
            do i = 1, nbgrma
                call jeexin(jexnum(grpmav, i), iret)
                if (iret .eq. 0) goto 70
                call jenuno(jexnum(grpmav, i), nom)
                call jelira(jexnum(grpmav, i), 'LONUTI', n1)
                write (ifm, 808) nom, n1
70              continue
            end do
        end if

        write (ifm, 803) lisnoe
        call jeveuo(cooval, 'L', iad)
        do i = 1, nbnoeu
            nom = int_to_char8(i)
            idec = iad+(i-1)*3
            write (ifm, 701) i, nom, (zr(idec+j-1), j=1, nbcoor)
        end do
!
        write (ifm, 803) lismai
        do i = 1, nbmail
            nom = int_to_char8(i)
            call jeveuo(jexnum(conxv, i), 'L', iad1)
            call jelira(jexnum(conxv, i), 'LONMAX', nbno)
            ityp = typmail(i)
            call jenuno(jexnum('&CATA.TM.NOMTM', ityp), type)
            if (nbno .le. 5) then
                write (ifm, 702) i, nom, type, (zi(iad1+j-1), j=1, nbno)
            else
                write (ifm, 702) i, nom, type, (zi(iad1+j-1), j=1, 5)
                write (ifm, 703) (zi(iad1+j-1), j=6, nbno)
            end if
        end do
!
        if (nbgrno .ne. 0) then
            write (ifm, 803) lisgrn
            do i = 1, nbgrno
                call jeexin(jexnum(grpnov, i), iret)
                if (iret .eq. 0) goto 100
                call jenuno(jexnum(grpnov, i), nom)
                call jeveuo(jexnum(grpnov, i), 'L', iad)
                call jelira(jexnum(grpnov, i), 'LONUTI', nbno)
                nn = nbno
                if (nn .le. 5) then
                    write (ifm, 704) i, nom, nbno, (zi(iad+j-1), j=1, nn)
                else
                    write (ifm, 704) i, nom, nbno, (zi(iad+j-1), j=1, 5)
                    write (ifm, 703) (zi(iad+j-1), j=6, nn)
                end if
100             continue
            end do
        end if
!
        if (nbgrma .ne. 0) then
            write (ifm, 803) lisgrm
            do i = 1, nbgrma
                call jeexin(jexnum(grpmav, i), iret)
                if (iret .eq. 0) goto 110
                call jenuno(jexnum(grpmav, i), nom)
                call jeveuo(jexnum(grpmav, i), 'L', iad)
                call jelira(jexnum(grpmav, i), 'LONUTI', nbma)
                nn = nbma
                if (nbma .le. 5) then
                    write (ifm, 704) i, nom, nbma, (zi(iad+j-1), j=1, nn)
                else
                    write (ifm, 704) i, nom, nbma, (zi(iad+j-1), j=1, 5)
                    write (ifm, 703) (zi(iad+j-1), j=6, nn)
                end if
110             continue
            end do
        end if
    end if
    write (ifm, 809)
!
    call jedema()
!
!
701 format(2x, i8, 2x, a24, 10x, 3(d14.5, 2x))
702 format(2x, i8, 2x, a24, 2x, a8, 5(2x, i8))
703 format(100(30x, 5(2x, i8),/))
704 format(2x, i8, 2x, a24, 2x, i8, 5(2x, i8))
801 format(a80)
802 format(/, '------------ MAILLAGE ', a8,&
     &       ' - IMPRESSIONS NIVEAU ', i2, ' ------------',/)
803 format(/, 15x, '------  ', a32, '  ------',/)
804 format(/, a32, i12)
806 format(30x, a8, 5x, i12)
808 format(30x, a24, 2x, i12)
809 format(/, 80('-'),/)
!
end subroutine
