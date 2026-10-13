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

subroutine get_med_connectivity(mesh, med_connect)

    implicit none

#include "asterfort/gnomsd_xx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lrmtyp.h"
#include "asterfort/jedupo.h"
#include "asterfort/jexatr.h"
#include "asterfort/assert.h"

    character(len=8), intent(in) :: mesh
    character(len=24), intent(inout) :: med_connect

! ------------------------------------------------------------------------------
!
! Return the connectivity of 'mesh' using Med conventions
!
! In  mesh          : name of mesh
! InOut med_connect : collection identical to 'mesh.CONNECT' but with the
!                     Med conventions
!    If it is empty as input (" "), a default temporary name is used.
!
! ------------------------------------------------------------------------------

    integer(kind=8) :: modnum(MT_NTYMAX)
    integer(kind=8) :: nnotyp(MT_NTYMAX)
    integer(kind=8) :: nuanom(MT_NTYMAX, MT_NNOMAX)
    character(len=24) :: connex
    integer(kind=8) :: i, j, ityp
    integer(kind=8) :: nbcell, nbnode, offset
    integer(kind=8), pointer :: dime(:) => null()
    integer(kind=8), pointer :: med_values(:) => null()
    integer(kind=8), pointer :: typmast(:) => null()
    integer(kind=8), pointer :: loncum(:) => null()

    integer(kind=8) :: cell_connect(MT_NNOMAX)

    call jemarq()

    if (med_connect .eq. ' ') then
        call gnomsd_xx(med_connect)
    end if

!   get conversion arrays code_aster <---> med
    call lrmtyp(nnotyp=nnotyp, modnum=modnum, nuanom=nuanom)

    connex = mesh//'.CONNEX'
    call jedupo(connex, 'V', med_connect, ASTER_TRUE)
    call jeveuo(med_connect, 'E', vi=med_values)
    call jeveuo(jexatr(med_connect, 'LONCUM'), 'L', vi=loncum)

    call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmast)
    call jeveuo(mesh//'.DIME', 'L', vi=dime)
    nbcell = dime(3)

    do i = 1, nbcell
        ityp = typmast(i)

        ! Nothing to do for identical numbering.
        if (modnum(ityp) == 0) cycle

        ! LONCUM contains 1-based offsets
        offset = loncum(i)
        nbnode = loncum(i+1)-offset
        ASSERT(nbnode .eq. nnotyp(ityp))
        ! Copy
        cell_connect(1:nbnode) = med_values(offset:offset+nbnode-1)

        ! Permutation
        do j = 1, nbnode
            med_values(offset+j-1) = cell_connect(nuanom(ityp, j))
        end do
    end do
    call jedema()

end subroutine
