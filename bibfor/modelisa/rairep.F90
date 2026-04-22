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
subroutine rairep(noma, ioc, km, rigiRep, nbgr, &
                  ligrma, zjdlm, nbno, tabnoe, rignoe, rirot, ndim, jdme)
!
!
    implicit none
    integer(kind=8) :: ioc, nbgr, nbno, ndim, jdme
    integer(kind=8) :: zjdlm(*)
    character(len=8) :: noma, tabnoe(*), km
    character(len=24) :: ligrma(nbgr)
    real(kind=8) :: rignoe(*), rirot(3)
!
! --------------------------------------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/in_liste_entier.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/ddot.h"
#include "asterfort/int_to_char8.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ii, ij, posi, in, inoe, iret, nb_ma_surf, NbNoeud, nbparno
    integer(kind=8) :: ldgm, ldgn, ldnm, ltyp, nb, ncg
    integer(kind=8) :: nfg, ngn, nm, nn, ntopo
    integer(kind=8) :: num_maille, NbMaille
    integer(kind=8) :: appui
!
    real(kind=8) :: x(9), y(9), z(9), rigiRep(6)
    real(kind=8) :: a(3), b(3), c(3), u(3)
    real(kind=8) :: coef, dist, xyzc(3)
    real(kind=8) :: r1, r2, r3, r4, r5, r6, rig3, rig4, rig5, rig6
    real(kind=8) :: surf, surtot
    real(kind=8) :: xx, yy, zz, xyzg(3)
!
    character(len=8) :: k8b, nomnoe, typm, nommai
    character(len=24) :: nomgr, magrno, magrma, manoma, matyma
!
    aster_logical :: lfonc, trans
!
    integer(kind=8), pointer :: parno(:) => null()
    integer(kind=8), pointer :: mailles_surf(:) => null()
    real(kind=8), pointer   :: coegro(:) => null()
    real(kind=8), pointer   :: coeno(:) => null()
    real(kind=8), pointer   :: surmai(:) => null()
    real(kind=8), pointer   :: coord(:) => null()
!
    character(len=8), pointer :: fongro(:) => null()
    character(len=16) :: ktyelm
!
    blas_int :: b_1, b_2, b_3
!
! --------------------------------------------------------------------------------------------------
    call jemarq()
    lfonc = .false.
!
    magrno = noma//'.GROUPENO'
    magrma = noma//'.GROUPEMA'
    manoma = noma//'.CONNEX'
    matyma = noma//'.TYPMAIL'

!
!   Coordonnées des noeuds
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=coord)
!
!   Récupération du centre
    call getvr8('RIGI_PARASOL', 'COOR_CENTRE', iocc=ioc, nbval=0, nbret=ncg)
    call getvtx('RIGI_PARASOL', 'GROUP_NO_CENTRE', iocc=ioc, nbval=0, vect=k8b, nbret=ngn)
    xyzg(:) = 0.0
    if (ncg .ne. 0) then
        ASSERT(ngn .eq. 0)
        call getvr8('RIGI_PARASOL', 'COOR_CENTRE', iocc=ioc, nbval=3, vect=xyzg, nbret=ncg)
    else if (ngn .ne. 0) then
        ASSERT(ncg .eq. 0)
        call getvtx('RIGI_PARASOL', 'GROUP_NO_CENTRE', iocc=ioc, nbval=1, vect=nomgr, nbret=ngn)
        call jeveuo(jexnom(magrno, nomgr), 'L', ldgn)
        inoe = zi(ldgn)
        nomnoe = int_to_char8(inoe)
        xyzg(1) = coord(3*(inoe-1)+1)
        xyzg(2) = coord(3*(inoe-1)+2)
        xyzg(3) = coord(3*(inoe-1)+3)
    else
        ASSERT(ASTER_FALSE)
    end if
!
!   Récuperation des coefs ou fonctions de groupe
    call getvr8('RIGI_PARASOL', 'COEF_GROUP', iocc=ioc, nbval=0, nbret=ncg)
    if (ncg .ne. 0) then
        AS_ALLOCATE(vr=coegro, size=nbgr)
        call getvr8('RIGI_PARASOL', 'COEF_GROUP', iocc=ioc, nbval=nbgr, vect=coegro, nbret=ncg)
    else
        AS_ALLOCATE(vk8=fongro, size=nbgr)
        lfonc = .true.
        call getvid('RIGI_PARASOL', 'FONC_GROUP', iocc=ioc, nbval=nbgr, vect=fongro, nbret=nfg)
    end if
!
!   La dimension de l'appui n'est pas encore determinée
    appui = -1
    if (ndim .eq. 2) then
        appui = 1
    end if
!
    NbMaille = 0
    NbNoeud = 0
    call jeveuo(matyma, 'L', ltyp)
    do ii = 1, nbgr
        call jelira(jexnom(magrma, ligrma(ii)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, ligrma(ii)), 'L', ldgm)
        do in = 0, nb-1
            num_maille = zi(ldgm+in)
            call jenuno(jexnum('&CATA.TE.NOMTE', zi(jdme-1+num_maille)), ktyelm)
            if ((ktyelm .eq. "MEC3TR7H") .or. (ktyelm .eq. "MEC3QU9H")) then
                call utmess('F', 'MODELISA6_39')
            end if
            if (num_maille .le. 0) then
                nommai = '????'
                call utmess('F', 'AFFECARAELEM_25', si=ioc, nk=2, valk=[ligrma(ii), nommai])
            else if (zjdlm(num_maille) .eq. 0) then
                nommai = int_to_char8(num_maille)
                call utmess('F', 'AFFECARAELEM_25', si=ioc, nk=2, valk=[ligrma(ii), nommai])
            end if
            NbMaille = NbMaille+1
            call jelira(jexnum(manoma, num_maille), 'LONMAX', nm)
            call jenuno(jexnum('&CATA.TM.NOMTM', zi(ltyp-1+num_maille)), typm)
            call dismoi('DIM_TOPO', typm, 'TYPE_MAILLE', repi=ntopo)
!
            if (appui .eq. -1) then
!               la dimension de la première maille définit l'appui
                appui = ntopo
            else if ((appui .eq. 1) .or. (appui .eq. 2)) then
                if (appui .ne. ntopo) then
                    call utmess('F', 'MODELISA6_35')
                end if
            else
                call utmess('F', 'MODELISA6_29')
            end if
            NbNoeud = NbNoeud+nm
        end do
    end do
    ASSERT(appui .ne. -1)
    ASSERT(NbMaille .ne. 0)
!
    b_1 = to_blas_int(1)
    b_2 = to_blas_int(2)
    b_3 = to_blas_int(3)
!
!   Coefficients des noeuds de l interface
    AS_ALLOCATE(vr=coeno, size=NbNoeud)
!   Participation des noeuds de l interface
    AS_ALLOCATE(vi=parno, size=NbNoeud)
!   Surfaces élémentaires de la maille
    AS_ALLOCATE(vr=surmai, size=NbMaille)
!   Numéro des mailles
    AS_ALLOCATE(vi=mailles_surf, size=NbMaille)
!
    nb_ma_surf = 0; nbparno = 0
    surtot = 0.0d0
    do ii = 1, nbgr
        call jelira(jexnom(magrma, ligrma(ii)), 'LONUTI', nb)
        call jeveuo(jexnom(magrma, ligrma(ii)), 'L', ldgm)
        cymaille: do in = 0, nb-1
            num_maille = zi(ldgm+in)
            if (.not. in_liste_entier(num_maille, mailles_surf(1:nb_ma_surf), posi)) then
                nb_ma_surf = nb_ma_surf+1
                ASSERT(nb_ma_surf .le. NbMaille)
                mailles_surf(nb_ma_surf) = num_maille
                posi = nb_ma_surf
            else
!               MESSAGE <A> si une maille est en double
                nommai = int_to_char8(num_maille)
                call utmess('A', 'AFFECARAELEM_24', si=ioc, nk=2, valk=[ligrma(ii), nommai])
                cycle cymaille
            end if
            call jelira(jexnum(manoma, num_maille), 'LONMAX', nm)
            call jeveuo(jexnum(manoma, num_maille), 'L', ldnm)
            xyzc(:) = 0.0d0
            do nn = 1, nm
                inoe = zi(ldnm+nn-1)
!               On enregistre le numéro du noeud dans parno, s'il n'y est pas déjà
                if (.not. in_liste_entier(inoe, parno(1:nbparno))) then
                    nbparno = nbparno+1
                    ASSERT(nbparno .le. NbNoeud)
                    parno(nbparno) = inoe
                end if
                x(nn) = coord(3*(inoe-1)+1)
                y(nn) = coord(3*(inoe-1)+2)
                z(nn) = coord(3*(inoe-1)+3)
                xyzc(1) = xyzc(1)+x(nn)
                xyzc(2) = xyzc(2)+y(nn)
                xyzc(3) = xyzc(3)+z(nn)
            end do
            xyzc(:) = xyzc(:)/nm
!
            if (appui .eq. 1) then
                a(1) = x(2)-x(1)
                a(2) = y(2)-y(1)
                a(3) = z(2)-z(1)
                surf = ddot(b_2, a, b_1, a, b_1)
                surmai(posi) = sqrt(surf)
            else if (appui .eq. 2) then
                a(1) = x(3)-x(1)
                a(2) = y(3)-y(1)
                a(3) = z(3)-z(1)
                if (nm .eq. 3 .or. nm .eq. 6 .or. nm .eq. 7) then
                    b(1) = x(2)-x(1)
                    b(2) = y(2)-y(1)
                    b(3) = z(2)-z(1)
                else if (nm .eq. 4 .or. nm .eq. 8 .or. nm .eq. 9) then
                    b(1) = x(4)-x(2)
                    b(2) = y(4)-y(2)
                    b(3) = z(4)-z(2)
                else
                    ASSERT(.false.)
                end if
                call provec(a, b, c)
                surf = ddot(b_3, c, b_1, c, b_1)
                surmai(posi) = sqrt(surf)*0.5d0
            else
                ASSERT(.false.)
            end if
            if (lfonc) then
                u(1:3) = xyzg(1:3)-xyzc(1:3)
                dist = ddot(b_3, u, b_1, u, b_1)
                dist = sqrt(dist)
                call fointe('F ', fongro(ii), 1, ['X'], [dist], coef, iret)
                surmai(posi) = surmai(posi)*coef
            else
                surmai(posi) = surmai(posi)*coegro(ii)
            end if
            surtot = surtot+surmai(posi)
!           Surface de la maille affectée à chacun des noeuds
            surmai(posi) = surmai(posi)/nm
        end do cymaille
    end do
!   Ceinture et bretelle
    ASSERT(nbparno .le. nbno)
    nbno = nbparno
!
!   Calcul des pondérations élémentaires
    do ii = 1, nb_ma_surf
        num_maille = mailles_surf(ii)
        call jelira(jexnum(manoma, num_maille), 'LONMAX', nm)
        call jeveuo(jexnum(manoma, num_maille), 'L', ldnm)
        do nn = 1, nm
            inoe = zi(ldnm+nn-1)
!           Le noeud doit être dans parno
            if (in_liste_entier(inoe, parno(1:nbparno), posi)) then
                coeno(posi) = coeno(posi)+surmai(ii)/surtot
            else
                ASSERT(.false.)
            end if
        end do
    end do
!
!   Calcul des raideurs de rotation
    rig3 = 0.0d0
    rig4 = 0.0d0
    rig5 = 0.0d0
    rig6 = 0.0d0
    do ij = 1, nbparno
        inoe = parno(ij)
        xx = coord(3*(inoe-1)+1)-xyzg(1)
        yy = coord(3*(inoe-1)+2)-xyzg(2)
        zz = coord(3*(inoe-1)+3)-xyzg(3)
        if (ndim .eq. 3) then
            rig3 = 0.d0
            rig4 = rig4+(rigiRep(2)*zz**2+rigiRep(3)*yy**2)*coeno(ij)
            rig5 = rig5+(rigiRep(1)*zz**2+rigiRep(3)*xx**2)*coeno(ij)
            rig6 = rig6+(rigiRep(2)*xx**2+rigiRep(1)*yy**2)*coeno(ij)
        else
            rig3 = rig3+(rigiRep(2)*xx**2+rigiRep(1)*yy**2)*coeno(ij)
        end if
    end do
!
    trans = (km(1:7) .eq. 'K_T_D_N') .or. (km(1:7) .eq. 'K_T_D_L') .or. &
            (km(1:7) .eq. 'A_T_D_N') .or. (km(1:7) .eq. 'A_T_D_L')
!
    if (trans) then
!       Pas de raideur en rotation sur les discrets
        if (ndim .eq. 2) then
            rigiRep(3) = 0.0d0
            rig3 = 0.0d0
        end if
        rigiRep(4) = 0.0d0
        rigiRep(5) = 0.0d0
        rigiRep(6) = 0.0d0
        rig4 = 0.0d0
        rig5 = 0.0d0
        rig6 = 0.0d0
        rirot(1) = 0.0d0
        rirot(2) = 0.0d0
        rirot(3) = 0.0d0
    else
        rig3 = rigiRep(3)-rig3
        rig4 = rigiRep(4)-rig4
        rig5 = rigiRep(5)-rig5
        rig6 = rigiRep(6)-rig6
        if (ndim .eq. 3) then
            rirot(1) = rig4
            rirot(2) = rig5
            rirot(3) = rig6
        else
            rirot(1) = rig3
        end if
    end if
!
    do ij = 1, nbparno
        inoe = parno(ij)
        r1 = rigiRep(1)*coeno(ij)
        r2 = rigiRep(2)*coeno(ij)
        if (ndim .eq. 3) then
            r3 = rigiRep(3)*coeno(ij)
            r4 = rig4*coeno(ij)
            r5 = rig5*coeno(ij)
            r6 = rig6*coeno(ij)
        else
            r3 = rig3*coeno(ij)
            r4 = 0.0d0
            r5 = 0.0d0
            r6 = 0.0d0
        end if
        nomnoe = int_to_char8(inoe)
        rignoe(6*(ij-1)+1) = r1
        rignoe(6*(ij-1)+2) = r2
        rignoe(6*(ij-1)+3) = r3
        rignoe(6*(ij-1)+4) = r4
        rignoe(6*(ij-1)+5) = r5
        rignoe(6*(ij-1)+6) = r6
        tabnoe(ij) = nomnoe
    end do
!
    AS_DEALLOCATE(vr=coegro)
    AS_DEALLOCATE(vk8=fongro)
    AS_DEALLOCATE(vr=coeno)
    AS_DEALLOCATE(vi=parno)
    AS_DEALLOCATE(vr=surmai)
    AS_DEALLOCATE(vi=mailles_surf)
!
    call jedema()
end subroutine
