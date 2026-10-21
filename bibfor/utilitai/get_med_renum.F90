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

subroutine get_med_renum(ityp, nno, noa, nom)

    implicit none

#include "asterfort/gnomsd_xx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lrmtyp.h"
#include "asterfort/wkvect.h"

    integer(kind=8), intent(out) :: ityp, nno, noa(*), nom(*)

! ------------------------------------------------------------------------------
!
! Return the renum map for a specific aster ityp
!
! ------------------------------------------------------------------------------

    character(len=8) :: nomtyp(MT_NTYMAX)
    integer(kind=8) :: i, nbtyp
    integer(kind=8) :: nnotyp(MT_NTYMAX), typgeo(MT_NTYMAX)
    integer(kind=8) :: renumd(MT_NTYMAX), modnum(MT_NTYMAX)
    integer(kind=8) :: nuanom(MT_NTYMAX, MT_NNOMAX)
    integer(kind=8) :: numnoa(MT_NTYMAX, MT_NNOMAX)

    call lrmtyp(nbtyp, nomtyp, nnotyp, typgeo, renumd, modnum, nuanom, numnoa)
    nno = nnotyp(ityp)
    do i = 1, nno
        if (modnum(ityp) .eq. 1) then
            noa(i) = numnoa(ityp, i)-1
            nom(i) = nuanom(ityp, i)-1
        else
            noa(i) = i-1
            nom(i) = i-1
        end if
    end do

end subroutine
