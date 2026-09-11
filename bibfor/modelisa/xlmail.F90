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
subroutine xlmail(fiss, nmaen1, nmaen2, nmaen3, nmafon, &
                  jmaen1, jmaen2, jmaen3, jmafon, nfon, &
                  jfon, jnofaf, nbfond, jbas, jtail, &
                  jfonmu, ndim, goinop)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: fiss
    integer(kind=8) :: nmaen1, nmaen2, nmaen3, nmafon
    integer(kind=8) :: jmaen1, jmaen2, jmaen3, jmafon
    integer(kind=8) :: nfon
    integer(kind=8) :: jfon, jnofaf, jbas, jtail, jfonmu
    integer(kind=8) :: nbfond, ndim
    aster_logical :: goinop
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (CREATION DES SD)
!
! CREATION SD SUR MAILLAGE XFEM
!
! ----------------------------------------------------------------------
!
!
! IN  FISS   : NOM DE LA FISSURE
! IN  NOMA   : NOM DU MAILLAGE
! IN  NMAEN1 : NOMBRE DE MAILLES 'HEAVISIDE'
! IN  NMAEN2 : NOMBRE DE MAILLES 'CRACKTIP'
! IN  NMAEN3 : NOMBRE DE MAILLES 'HEAVISIDE-CRACKTIP'
! IN  NMAFON : NOMBRE DE MAILLES CONTENANT LE FOND DE FISSURE
! IN  JMAEN1 : POINTEUR SUR MAILLES 'HEAVISIDE'
! IN  JMAEN2 : POINTEUR SUR MAILLES 'CRACKTIP'
! IN  JMAEN3 : POINTEUR SUR MAILLES 'HEAVISIDE-CRACKTIP'
! IN  JMAFON : POINTEUR SUR MAILLES CONTENANT LE FOND DE FISSURE
! IN  NFON   : NOMBRE DE POINTS DE FOND DE FISSURE
! IN  JFON   : POINTEUR SUR POINTS DE FOND DE FISSURE
! IN  JNOFAF : POINTEUR SUR NUMERO DES NOEUDS DES FACES DES ELEMENTS
!              PARENTS QUI CONTIENNENT LES POINTS DU FOND DE FISSURE
! IN  JBAS   : POINTEUR SUR DIRECTION DE PROPAGATION
! IN  JTAIL  : POINTEUR SUR TAILLES MAXIMALES DE MAILLES
! IN  NBFOND : NOMBRE DE FONDS DE FISSURES DETECTES
! IN  JFONMU : POINTEUR SUR DEBUT ET ARRIVEE DU FOND DE FISSURE
! IN  GOINOP : .TRUE. SI PASSAGE DANS OPOO10 (UPWIND-SIMPLEXE/GRILLE/3D)
!              .FALSE. SINON
!
    character(len=24) :: xheav, xctip, xhect, xmafon, xfonfi, xbasfo, xfonmu
    character(len=24) :: xtailr, xnofaf
    character(len=24) :: xfonfg
!
    integer(kind=8) :: jma1, jma2, jma3, jma4, jfo, jfomu, jba, jta, jnf
    integer(kind=8) :: i, k
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- ACCES AUX OBJETS
!
    if (.not. goinop) then
        xheav = fiss(1:8)//'.MAILFISS.HEAV'
        xctip = fiss(1:8)//'.MAILFISS.CTIP'
        xhect = fiss(1:8)//'.MAILFISS.HECT'
        xmafon = fiss(1:8)//'.MAILFISS.MAFOND'
        xfonfi = fiss(1:8)//'.FONDFISS'
        xnofaf = fiss(1:8)//'.NOFACPTFON'
        xbasfo = fiss(1:8)//'.BASEFOND'
        xtailr = fiss(1:8)//'.FOND.TAILLE_R'
        xfonmu = fiss(1:8)//'.FONDMULT'
!
! --- ENREGISTREMENT DES GROUP_MA 'HEAVISIDE'
!
        if (nmaen1 .ne. 0) then
            call wkvect(xheav, 'G V I', nmaen1, jma1)
            do i = 1, nmaen1
                zi(jma1-1+i) = zi(jmaen1-1+i)
            end do
        end if
!
! --- ENREGISTREMENT DES GROUP_MA 'CRACKTIP'
!
        if (nmaen2 .ne. 0) then
            call wkvect(xctip, 'G V I', nmaen2, jma2)
            do i = 1, nmaen2
                zi(jma2-1+i) = zi(jmaen2-1+i)
            end do
        end if
!
! --- ENREGISTREMENT DES GROUP_MA ''HEAVISIDE-CRACKTIP'
!
        if (nmaen3 .ne. 0) then
            call wkvect(xhect, 'G V I', nmaen3, jma3)
            do i = 1, nmaen3
                zi(jma3-1+i) = zi(jmaen3-1+i)
            end do
        end if
!
! --- ENREGISTREMENT DES MAILLES CONTENANT LE FOND DE FISSURE
!
        if (nmafon .ne. 0) then
            call wkvect(xmafon, 'G V I', nmafon, jma4)
            do i = 1, nmafon
                zi(jma4-1+i) = zi(jmafon-1+i)
            end do
        end if
!
! --- ENREGISTREMENT DES COORD ET DES ABS CURV
!
        if (nfon .ne. 0) then
            call wkvect(xfonfi, 'G V R', 4*nfon, jfo)
            do i = 1, nfon
                do k = 1, 4
                    zr(jfo-1+4*(i-1)+k) = zr(jfon-1+4*(i-1)+k)
                end do
            end do
!
            call wkvect(xnofaf, 'G V I', 4*nfon, jnf)
            do i = 1, nfon
                do k = 1, 4
                    zi(jnf-1+4*(i-1)+k) = zi(jnofaf-1+4*(i-1)+k)
                end do
            end do
!
            ! consistent with '&&XENRCH.BASFO'
            call wkvect(xbasfo, 'G V R', 2*ndim*(nfon+1), jba)
            do i = 1, nfon
                do k = 1, ndim
                    zr(jba-1+2*ndim*(i-1)+k) = zr(jbas-1+2*ndim*(i-1)+k)
                    zr(jba-1+2*ndim*(i-1)+k+ndim) = zr(jbas-1+2*ndim*(i-1)+k+ndim)
                end do
            end do
!
            call wkvect(xtailr, 'G V R', nfon, jta)
            do i = 1, nfon
                zr(jta-1+i) = zr(jtail-1+i)
            end do
!
        end if
!
! --- ENREGISTREMENT DES FONDS MULTIPLES
!
        if (nbfond .ne. 0) then
            call wkvect(xfonmu, 'G V I', 2*nbfond, jfomu)
            do i = 1, nbfond
                zi(jfomu-1+2*(i-1)+1) = zi(jfonmu-1+2*(i-1)+1)
                zi(jfomu-1+2*(i-1)+2) = zi(jfonmu-1+2*(i-1)+2)
            end do
        end if
!
!
    else if (goinop) then
!     CREATION DU FOND DE FISSURE SUR LA GRILLE FONDFISG
!     ET DES POINTS VIRTUELS FONDFISV
!     (UNIQUEMENT SI PROPAGATION UPWIND/SIMPLEXE AVEC GRILLE EN 3D)
        xfonfg = fiss(1:8)//'.FONDFISG'
!
! --- ENREGISTREMENT DES COORD SUR LA GRILLE
!
        call wkvect(xfonfg, 'G V R', 4*nfon, jfo)
        do i = 1, nfon
            do k = 1, 4
                zr(jfo-1+4*(i-1)+k) = zr(jfon-1+4*(i-1)+k)
            end do
        end do
    end if
!
    call jedema()
end subroutine
