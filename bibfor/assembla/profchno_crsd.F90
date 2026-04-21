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

subroutine profchno_crsd(nume_equaz, base, nb_equa, meshz, nb_ligrz, &
                         nb_ecz, gran_namez, prno_lengthz, l_coll_const)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jedetr.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecreo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jenonu.h"
#include "asterfort/wkvect.h"
#include "asterfort/isParallelMesh.h"
!
!
    character(len=*), intent(in) :: nume_equaz
    character(len=1), intent(in) :: base
    integer(kind=8), intent(in) :: nb_equa
    character(len=*), optional, intent(in) :: meshz
    character(len=*), optional, intent(in) :: gran_namez
    integer(kind=8), optional, intent(in) :: nb_ecz
    integer(kind=8), optional, intent(in) :: nb_ligrz
    integer(kind=8), optional, intent(in) :: prno_lengthz
    aster_logical, optional, intent(in) :: l_coll_const
!
! --------------------------------------------------------------------------------------------------
!
! NUME_EQUA
!
! Create object
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_equa    : name of NUME_EQUA
! In  base         : JEVEUX base to create NUME_EQUA
! In  nb_equa      : number of equations
! In  nb_ligr      : number of LIGREL in .LILI object
!                   if not present => only mesh (nb_ligr=1)
! In  mesh         : name of mesh
! In  gran_name    : name of GRANDEUR
! In  nb_ec        : number of coding integers for GRANDEUR
!                   if not present => get from gran_name
! In  prno_length  : length of first PRNO object (on mesh)
! In  l_coll_const : .true. if PRNO colelction is CONSTANT (not variable) for cnscno
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: nume_equa
    integer(kind=8) :: i_equa, i_ligr_mesh, nb_node_mesh, nb_ec, nb_ligr, prno_length
    integer(kind=8), pointer :: prchno_nueq(:) => null()
    integer(kind=8), pointer :: prchno_deeq(:) => null()
    aster_logical :: l_pmesh
!
! --------------------------------------------------------------------------------------------------
!
    prno_length = 0
    nb_ec = 0
    nume_equa = nume_equaz
    call jedetr(nume_equa//'.DEEQ')
    call jedetr(nume_equa//'.LILI')
    call jedetr(nume_equa//'.NUEQ')
    call jedetr(nume_equa//'.PRNO')
    if (present(nb_ligrz)) then
        nb_ligr = nb_ligrz
    else
        nb_ligr = 1
    end if
    if (present(nb_ecz)) then
        nb_ec = nb_ecz
    else
        ASSERT(present(gran_namez) .or. present(prno_lengthz))
        if (present(gran_namez)) then
            call dismoi('NB_EC', gran_namez, 'GRANDEUR', repi=nb_ec)
        end if
    end if
    if (present(meshz)) then
        ASSERT(.not. present(prno_lengthz))
        call dismoi('NB_NO_MAILLA', meshz, 'MAILLAGE', repi=nb_node_mesh)
        l_pmesh = isParallelMesh(meshz)
        if (.not. l_pmesh) then
            ASSERT(nb_equa > 0)
        end if
    end if
!
! - Create object NUEQ
!
    call wkvect(nume_equa//'.NUEQ', base//' V I', max(1, nb_equa), vi=prchno_nueq)
    call jeecra(nume_equa//'.NUEQ', "LONUTI", nb_equa)
!
! - Set to identity
!
    do i_equa = 1, nb_equa
        prchno_nueq(i_equa) = i_equa
    end do
!
! - Create object DEEQ
!
    call wkvect(nume_equa//'.DEEQ', base//' V I', max(1, 2*nb_equa), vi=prchno_deeq)
    call jeecra(nume_equa//'.DEEQ', "LONUTI", 2*nb_equa)
!
! - Create object LILI (name repertory)
!
    call jecreo(nume_equa//'.LILI', base//' N K24')
    call jeecra(nume_equa//'.LILI', 'NOMMAX', nb_ligr)
!
! - Create &MAILLA object in LILI
!
    call jecroc(jexnom(nume_equa(1:19)//'.LILI', '&MAILLA'))
    call jenonu(jexnom(nume_equa//'.LILI', '&MAILLA'), i_ligr_mesh)
    ASSERT(i_ligr_mesh .eq. 1)
!
! - Length of first PRNO object (on mesh)
!
    if (present(meshz)) then
        prno_length = (2+nb_ec)*nb_node_mesh
    else
        prno_length = prno_lengthz
    end if
!
! - Create object PRNO (collection)
!
    if (present(l_coll_const)) then
        if (l_coll_const) then
            call jecrec(nume_equa//'.PRNO', base//' V I', 'NU', 'CONTIG', 'CONSTANT', nb_ligr)
        else
            call jecrec(nume_equa//'.PRNO', base//' V I', 'NU', 'CONTIG', 'VARIABLE', nb_ligr)
        end if
    else
        call jecrec(nume_equa//'.PRNO', base//' V I', 'NU', 'CONTIG', 'VARIABLE', nb_ligr)
    end if
!
! - Length of &MAILLA object in PRNO
!
    call jeecra(jexnum(nume_equa//'.PRNO', i_ligr_mesh), 'LONMAX', max(prno_length, 1))
    call jeecra(jexnum(nume_equa//'.PRNO', i_ligr_mesh), 'LONUTI', prno_length)
!
end subroutine
