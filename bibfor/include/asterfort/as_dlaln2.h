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
    subroutine as_dlaln2(ltrans, na, nw, smin, ca,&
                      a, lda, d1, d2, b,&
                      ldb, wr, wi, x, ldx,&
                      scale, xnorm, info)
        integer(kind=8) :: ldx
        integer(kind=8) :: ldb
        integer(kind=8) :: lda
        aster_logical :: ltrans
        integer(kind=8) :: na
        integer(kind=8) :: nw
        real(kind=8) :: smin
        real(kind=8) :: ca
        real(kind=8) :: a(lda, *)
        real(kind=8) :: d1
        real(kind=8) :: d2
        real(kind=8) :: b(ldb, *)
        real(kind=8) :: wr
        real(kind=8) :: wi
        real(kind=8) :: x(ldx, *)
        real(kind=8) :: scale
        real(kind=8) :: xnorm
        integer(kind=8) :: info
    end subroutine as_dlaln2
end interface
