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
    subroutine as_dlasy2(ltranl, ltranr, isgn, n1, n2,&
                      tl, ldtl, tr, ldtr, b,&
                      ldb, scale, x, ldx, xnorm,&
                      info)
        integer(kind=8) :: ldx
        integer(kind=8) :: ldb
        integer(kind=8) :: ldtr
        integer(kind=8) :: ldtl
        aster_logical :: ltranl
        aster_logical :: ltranr
        integer(kind=8) :: isgn
        integer(kind=8) :: n1
        integer(kind=8) :: n2
        real(kind=8) :: tl(ldtl, *)
        real(kind=8) :: tr(ldtr, *)
        real(kind=8) :: b(ldb, *)
        real(kind=8) :: scale
        real(kind=8) :: x(ldx, *)
        real(kind=8) :: xnorm
        integer(kind=8) :: info
    end subroutine as_dlasy2
end interface
