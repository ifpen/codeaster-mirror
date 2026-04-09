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
    subroutine as_zgetv0(ido, bmat, initv, n, j,&
                      v, ldv, resid, rnorm, ipntr,&
                      workd, ierr, alpha)
        integer(kind=8) :: ldv
        integer(kind=8) :: j
        integer(kind=8) :: n
        integer(kind=8) :: ido
        character(len=1) :: bmat
        aster_logical :: initv
        complex(kind=8) :: v(ldv, j)
        complex(kind=8) :: resid(n)
        real(kind=8) :: rnorm
        integer(kind=8) :: ipntr(3)
        complex(kind=8) :: workd(2*n)
        integer(kind=8) :: ierr
        real(kind=8) :: alpha
    end subroutine as_zgetv0
end interface
