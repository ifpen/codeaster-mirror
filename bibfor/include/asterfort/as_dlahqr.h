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
!
#include "asterf_types.h"
!
interface
    subroutine as_dlahqr(wantt, wantz, n, ilo, ihi,&
                      h, ldh, wr, wi, iloz,&
                      ihiz, z, ldz, info)
        integer(kind=8) :: ldz
        integer(kind=8) :: ldh
        aster_logical :: wantt
        aster_logical :: wantz
        integer(kind=8) :: n
        integer(kind=8) :: ilo
        integer(kind=8) :: ihi
        real(kind=8) :: h(ldh, *)
        real(kind=8) :: wr(*)
        real(kind=8) :: wi(*)
        integer(kind=8) :: iloz
        integer(kind=8) :: ihiz
        real(kind=8) :: z(ldz, *)
        integer(kind=8) :: info
    end subroutine as_dlahqr
end interface
