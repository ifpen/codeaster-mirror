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
interface
    subroutine as_dneigh(rnorm, n, h, ldh, ritzr,&
                      ritzi, bounds, q, ldq, workl,&
                      ierr)
        integer(kind=8) :: ldq
        integer(kind=8) :: ldh
        integer(kind=8) :: n
        real(kind=8) :: rnorm
        real(kind=8) :: h(ldh, n)
        real(kind=8) :: ritzr(n)
        real(kind=8) :: ritzi(n)
        real(kind=8) :: bounds(n)
        real(kind=8) :: q(ldq, n)
        real(kind=8) :: workl(n*(n+3))
        integer(kind=8) :: ierr
    end subroutine as_dneigh
end interface
