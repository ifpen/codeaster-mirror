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
subroutine oriintpou(mesh, listCellNume, nbCell, norien, vxorie)
!
    use mesh_module, only: checkCellsAreVolume
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfvf.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "MeshTypes_type.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: nbCell
    integer(kind=8), pointer :: listCellNume(:)
    integer(kind=8), intent(out) :: norien
    real(kind=8), intent(in) :: vxorie(3)
!.======================================================================
!
!   ORHEXA  --  LE BUT EST QUE POUR TOUTES LES MAILLES DE LA LISTE,
!               L'AXE X DU REPERE LOCAL SOIT ORIENTE SUIVANT LE VECTEUR.
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NOMA           IN    K8      NOM DU MAILLAGE
!    LISTMA         IN    I       LISTE DES MAILLES A REORIENTER
!    NBMAIL         IN    I       NB DE MAILLES DE LA LISTE
!    NORIEN        VAR            NOMBRE DE MAILLES REORIENTEES
!    VXORIE         IN    R       VECTEUR DIRECTEUR
!.========================= DEBUT DES DECLARATIONS ====================
! -----  VARIABLES LOCALES
    integer(kind=8) :: iCell, cellNume
    integer(kind=8) :: jcoor, p1, p2, ifm, niv
    integer(kind=8) :: jdesm1
    integer(kind=8), parameter :: nbnds = 27
    aster_logical :: hasVolume, hasVoluNotHexaBiQ, onlyHexaBiQuad
    character(len=8), pointer :: ori5(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    !
    integer(kind=8), parameter :: dp = kind(0.0d0)
    integer(kind=8) :: i, ii, ic, ideb
    real(dp) :: vec(3), projs(3), vx(3), vy(3), vz(3)
    real(dp) :: vxn(3), vyn(3), vzn(3), nvxn, nvyn, nvzn, nvec
    real(dp) :: ni(3, 2), rot_mat(3, 3), ml(3, MT_NNOMAX3D), ff(nbnds)
    real(dp) :: coorn(3)
    integer(kind=8) :: connn(nbnds), nrange(nbnds), couples(2, 3)
    integer(kind=8) :: imax, psign, cxn(2), ndnum, ndpos
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
    norien = 0
    if (nbCell .eq. 0) goto 999
!
    call infniv(ifm, niv)

! - Initialization
    onlyHexaBiQuad = ASTER_TRUE
    ! node couples in the center of opposite faces (in parametric space)
    couples = reshape((/23, 25, 22, 24, 21, 26/), shape(couples))
    ! vector basis of parametric space
    vx = (/1.0_dp, 0.0_dp, 0.0_dp/)
    vy = (/0.0_dp, 1.0_dp, 0.0_dp/)
    vz = (/0.0_dp, 0.0_dp, 1.0_dp/)
    call elrfno('H27', nodeCoor=ml)
    nrange = [(ii, ii=1, nbnds)]

! - Options

! - Access to mesh datastructures
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
    call jeveuo(mesh//'.COORDO    .VALE', 'L', jcoor)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', p2)
    call jeveuo(mesh//'.CONNEX', 'E', p1)
    ! p1 points to the start of the table of connectivity of the whole mesh.
    !  It contains the tables of connectivity of all the elements one after the
    !  other in the order of the elements
    ! p2 is primarily used to get the values in p1. It is the pointer to the
    !  cumulated lengths of the table of connectivity. So it contains, at the
    !  index of the element number, the position in p1 of the start of the table
    !  of connectivity of the element

! - Working vectors
    AS_ALLOCATE(vk8=ori5, size=nbCell)

! - Check type of cells (only bi-quadratic hexa)
    call checkCellsAreVolume(mesh, &
                             nbCell, listCellNume, &
                             onlyHexaBiQuad, &
                             ori5, &
                             hasVolume, hasVoluNotHexaBiQ)

! - Iterate over mesh elements
    do iCell = 1, nbCell
        !
        cellNume = listCellNume(iCell)
        jdesm1 = zi(p2+cellNume-1)-1

        ! Projection of new x-vector on the 3 vectors that join face centers
        projs = 0.0_dp
        do i = 1, 3
            ni = 0.0_dp
            do ii = 1, 2
                ! node number in mesh of the node positioned at couples(ii, i)
                !  in the original element connectivity
                ndnum = zi(p1+jdesm1-1+couples(ii, i))
                ideb = jcoor-1+3*(ndnum-1)
                do ic = 1, 3
                    ! coordinates of that node in 3D-physical space
                    ni(ic, ii) = zr(ideb+ic)
                end do
            end do
            ! vector that joins two face centers
            vec = ni(:, 2)-ni(:, 1)
            call normev(vec, nvec)
            projs(i) = dot_product(vxorie, vec)
        end do

        ! Get nodes of x-vector in new orientation
        imax = maxloc(abs(projs), 1)
        psign = int(sign(1.0_dp, projs(imax)))
        if (psign > 0) then
            cxn = couples(:, imax)
        else
            cxn = (/couples(2, imax), couples(1, imax)/)
        end if

        ! reorientation is not needed
        if (all(cxn .eq. (/25, 23/))) goto 200

        ! Get new basis in parametric space
        ! new x-vector
        vxn = ml(:, cxn(2))-ml(:, cxn(1))
        if (all(cxn .eq. (/23, 25/))) then
            ! arbitrarily the opposite of the initial y-vector is chosen as the
            !  new y-vector
            vyn = ml(:, 22)-ml(:, 24)
        else
            ! arbitrarily the initial x-vector (in parametric space) is chosen
            !  as the new y-vector
            vyn = ml(:, 23)-ml(:, 25)
        end if
        ! cross product of vxn and vyn
        call provec(vxn, vyn, vzn)
        ! normalize basis
        call normev(vxn, nvxn)
        call normev(vyn, nvyn)
        call normev(vzn, nvzn)

        ! Define the rotation matrix between bases
        rot_mat = reshape((/vxn, vyn, vzn/), shape(rot_mat))

        ! Get new table of connectivity by identifying node positions in
        !  connectivity from rotated coordinates in parametric space
        do i = 1, nbnds
            ! apply the transformation on element node coordinates in parametric
            !  space ==> new coordinates
            coorn = matmul(ml(:, i), rot_mat)
            ! evaluate shape functions at new coordinates
            call elrfvf('H27', coorn, ff)
            ! identify node position from coordinates
            ndpos = nint(dot_product(nrange, ff))
            ! new table of connectivity
            connn(ndpos) = zi(p1+jdesm1-1+i)
        end do

        ! Update connectivity
        do i = 1, nbnds
            zi(p1+jdesm1-1+i) = connn(i)
        end do

        ! Inform on reorientation
        norien = norien+1

200     continue
!
    end do
!
    AS_DEALLOCATE(vk8=ori5)
!
999 continue
    call jedema()
end subroutine
