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

subroutine uteref2(ntypef, imolo, renum, &
                   elrefe, famil, nomtypmail, nnos, nno, nbpg, ndim, dimrefe, &
                   refcoo, gscoo, wg, codret)
!
!-----------------------------------------------------------------------
!     UTILITAIRE - ELEMENT DE REFERENCE
!     --           -          ---
!-----------------------------------------------------------------------
!     ENTREES :
!       NTYPEF : TYPE D'ELEMENT ASTER
!       IMOLO  : MODE LOCAL
!       RENUM  : TRUE POUR RENUMEROTER REFCOO AU FORMAT MED
!     SORTIES :
!       ELREFE : NOM DE L'ELEMENT DE REFERENCE
!       FAMIL : NOM DE LA FAMILLE DE POINTS DE GAUSS
!       NOMFPG : NOM DU TYPE D'ELEMENT
!       NNOS   : NOMBRE DE NOEUDS SOMMETS
!       NNO    : NOMBRE DE NOEUDS TOTAL
!       NBPG   : NOMBRE DE POINTS DE GAUSS
!       NDIM   : DIMENSION DE L'ELEMENT
!       DIMREFE: DIMENSION DE L'ELEMENT DE REFERENCE
!       REFCOO : COORDONNEES DES NNO NOEUDS
!       GSCOO  : COORDONNEES DES POINTS DE GAUSS, SI CHAMP ELGA
!       WG     : POIDS DES POINTS DE GAUSS, SI CHAMP ELGA
!       CODRET : CODE DE RETOUR
!                0 : PAS DE PB
!                1 : CORRECTION POUR ELEMENTS JOINT NON PROGRAMMEE
!                2 : FAMILLE NON GEREE AU FORMAT MED
!                3 : FAMILLE LISTE
!     REMARQUE :
!     ON DOIT RETOURNER LES COORDONNEES SOUS LA FORME :
!     . ELEMENT 1D : X1 X2 ... ... XN
!     . ELEMENT 2D : X1 Y1 X2 Y2 ... ... XN YN
!     . ELEMENT 3D : X1 Y1 Z1 X2 Y2 Z2 ... ... XN YN ZN
!     C'EST CE QUE MED APPELLE LE MODE ENTRELACE
!     ON DOIT RETOURNER LES POIDS SOUS LA FORME :
!     WG1 WG2 ... ... WGN
!-----------------------------------------------------------------------
!
    implicit none

#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/gnomsd_xx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/dismoi.h"
#include "asterfort/elref2.h"
#include "asterfort/jexnum.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnom.h"
#include "asterfort/jenonu.h"
#include "asterfort/utmess.h"
#include "asterfort/elraca.h"
#include "asterfort/elraga.h"
#include "asterfort/lrmtyp.h"
#include "MeshTypes_type.h"

    integer(kind=8), intent(in) ::  ntypef, imolo
    aster_logical, intent(in) :: renum
    character(len=8), intent(out) :: elrefe, famil, nomtypmail
    integer(kind=8), intent(out) :: nnos, nno, nbpg, ndim, dimrefe
    real(kind=8), intent(out) :: refcoo(*), gscoo(*), wg(*)
    integer(kind=8), intent(out) ::  codret

    integer(kind=8) :: i, nec, kfpg, nbelr, jmolo, shift
    character(len=8) :: lielrf(MT_NBFAMX), nomgd, elrefe_joint
    character(len=16) :: nomte, nomfpg
    aster_logical :: ljoint, lpenta

    integer(kind=8) :: nbpg00(MT_NBFAMX), modnum(MT_NTYMAX), nbfpg, nbfpg2, ifam, ipg
    integer(kind=8) :: idime, dimtopo, typmail, nuanom(MT_NTYMAX, MT_NNOMAX)
!
    integer(kind=8), parameter:: lgmax = 1000
!
    real(kind=8) :: gscoo2(3*lgmax), refcoo2(3*lgmax)

!
    character(len=8) :: fapg(MT_NBFAMX), fapg2(MT_NBFAMX)
!
    real(kind=8) :: vol

    elrefe = ' '
    ljoint = ASTER_FALSE
    lpenta = ASTER_FALSE
    famil = ' '
    codret = 0

    call jemarq()

    if ((ntypef .gt. 0) .and. (imolo .gt. 0)) then

        !     DETERMINATION DE ELREFE

        call jenuno(jexnum('&CATA.TE.NOMTE', ntypef), nomte)
        call elref2(nomte, MT_NBFAMX, lielrf, nbelr)
        ASSERT(nbelr .gt. 0)
        elrefe = lielrf(1)

        call dismoi('DIM_TOPO', nomte, 'TYPE_ELEM', repi=dimtopo)
        call dismoi('NOM_TYPMAIL', nomte, 'TYPE_ELEM', repk=nomtypmail)
        call dismoi('NUM_TYPMAIL', nomtypmail, 'TYPE_MAILLE', repi=typmail)

        ! DETERMINATION DE LA FAMILLE DE POINTS DE GAUSS

        call jeveuo(jexnum('&CATA.TE.MODELOC', imolo), 'L', jmolo)

        ! Assert ELGA
        ASSERT(zi(jmolo-1+1) .eq. 3)

        ! Get NOM_GRANDEUR
        call jenuno(jexnum('&CATA.GD.NOMGD', zi(jmolo-1+2)), nomgd)
        call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nec)

        kfpg = zi(jmolo-1+4+nec+1)
        ! Assert famille ordinaire, voir manopg.F90 pour gestion familles liste
        if (kfpg .lt. 0) then
            codret = 3
        end if
    end if

    if (codret .eq. 0) then
        !
        call jenuno(jexnum('&CATA.TM.NOFPG', kfpg), nomfpg)
        ASSERT(elrefe .eq. nomfpg(1:8))
        famil = nomfpg(9:16)
    end if

    if (codret .eq. 0) then
        call elraca(elrefe, &
                    nbfpg, fapg, nbpg00, &
                    ndim, nno, nnos, &
                    refcoo2, vol)
!
        dimrefe = ndim
!       Glute pour les joints
        elrefe_joint = elrefe
        if (dimtopo .ne. ndim) then
            select case (nomtypmail)
            case ('HEXA8')
                elrefe_joint = 'HE8'
            case ('QUAD8')
                elrefe_joint = 'QU8'
            case ('PENTA6')
                elrefe_joint = 'PE6'
                lpenta = ASTER_TRUE
            case ('PENTA15')
                elrefe_joint = 'P15'
                lpenta = ASTER_TRUE
            case ('HEXA20')
                elrefe_joint = 'H20'
            case default
                codret = 1
            end select
            if (codret .eq. 0) then
                call elraca(elrefe_joint, &
                            nbfpg2, fapg2, nbpg00, &
                            dimtopo, nno, nnos, &
                            refcoo2, vol)
                ljoint = ASTER_TRUE
            end if
        end if
    end if
!
    if (codret .eq. 0) then
        ASSERT(nbfpg .le. MT_NBFAMX)
        ASSERT(nno .le. MT_NNOMAX)
        ifam = indik8(fapg, famil, 1, nbfpg)
        if (ifam .le. 0) then
            codret = 2
        end if
    end if

    if (codret .eq. 0) then
        if (ljoint) then
            call elraga(elrefe, famil, ndim, nbpg, gscoo2, wg)
            if (lpenta) then
                do ipg = 1, nbpg
                    gscoo((ipg-1)*dimtopo+1) = 0.0d0
                    do idime = 2, dimtopo
                        gscoo((ipg-1)*dimtopo+idime) = gscoo2((ipg-1)*ndim+idime-1)
                    end do
                end do
            else
                do ipg = 1, nbpg
                    do idime = 1, ndim
                        gscoo((ipg-1)*dimtopo+idime) = gscoo2((ipg-1)*ndim+idime)
                    end do
                    do idime = ndim+1, dimtopo
                        gscoo((ipg-1)*dimtopo+idime) = 0.0d0
                    end do
                end do
            end if
            ndim = dimtopo
        else
            call elraga(elrefe, famil, ndim, nbpg, gscoo, wg)
        end if
!
        ASSERT(nbpg .le. 125)

    end if

    if (codret .eq. 0 .and. renum) then
        call lrmtyp(modnum=modnum, nuanom=nuanom)
    end if

    if (codret .eq. 0 .and. renum .and. modnum(typmail) .eq. 1) then
        ! RENUMEROTATION DES NOEUDS SI NECESSAIRE
        do i = 1, nno
            shift = ndim*(nuanom(typmail, i)-1)
            refcoo(ndim*(i-1)+1:ndim*i) = refcoo2(shift+1:shift+ndim)
        end do
    else
        refcoo(1:ndim*nno) = refcoo2(1:ndim*nno)
    end if
    call jedema()

end subroutine
