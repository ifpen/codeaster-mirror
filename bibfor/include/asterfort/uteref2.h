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

#include "asterf_types.h"

interface
    subroutine uteref2(ntypef, imolo, renum, &
                       elrefe, famil, nomtypmail, nnos, nno, nbpg, ndim, dimtopo, &
                       refcoo, gscoo, wg, codret)

        integer(kind=8), intent(in) ::  ntypef, imolo
        aster_logical, intent(in) :: renum
        character(len=8), intent(out) :: elrefe, famil, nomtypmail
        integer(kind=8), intent(out) :: nnos, nno, nbpg, ndim, dimtopo
        real(kind=8), intent(out) :: refcoo(*), gscoo(*), wg(*)
        integer(kind=8), intent(out) ::  codret
    end subroutine uteref2
end interface
