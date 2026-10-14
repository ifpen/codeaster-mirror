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

subroutine get_med_map(nbtyp, aster, med)

    implicit none

#include "asterfort/gnomsd_xx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lrmtyp.h"
#include "asterfort/wkvect.h"

    integer(kind=8), intent(out) :: nbtyp, aster(*), med(*)

! ------------------------------------------------------------------------------
!
! Return all the possible MED cell types and their corresponding aster type
!
! ------------------------------------------------------------------------------

    character(len=8) :: nomtyp(MT_NTYMAX)
    integer(kind=8) :: i
    integer(kind=8) :: nnotyp(MT_NTYMAX), typgeo(MT_NTYMAX)
    integer(kind=8) :: renumd(MT_NTYMAX), modnum(MT_NTYMAX)
    integer(kind=8) :: nuanom(MT_NTYMAX, MT_NNOMAX)
    integer(kind=8) :: numnoa(MT_NTYMAX, MT_NNOMAX)

    call lrmtyp(nbtyp, nomtyp, nnotyp, typgeo, renumd, modnum, nuanom, numnoa)

    do i = 1, nbtyp
        aster(i) = renumd(i)
        med(i) = typgeo(renumd(i))
    end do

end subroutine
