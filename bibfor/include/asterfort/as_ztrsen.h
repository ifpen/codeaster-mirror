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
    subroutine as_ztrsen(select, n, t, ldt, q,&
                      ldq, w, m, info)
        integer(kind=8) :: ldq
        integer(kind=8) :: ldt
        aster_logical :: select(*)
        integer(kind=8) :: n
        complex(kind=8) :: t(ldt, *)
        complex(kind=8) :: q(ldq, *)
        complex(kind=8) :: w(*)
        integer(kind=8) :: m
        integer(kind=8) :: info
    end subroutine as_ztrsen
end interface
