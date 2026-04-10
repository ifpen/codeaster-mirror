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
    subroutine as_dlaexc(wantq, n, t, ldt, q,&
                      ldq, j1, n1, n2, work,&
                      info)
        integer(kind=8) :: ldq
        integer(kind=8) :: ldt
        aster_logical :: wantq
        integer(kind=8) :: n
        real(kind=8) :: t(ldt, *)
        real(kind=8) :: q(ldq, *)
        integer(kind=8) :: j1
        integer(kind=8) :: n1
        integer(kind=8) :: n2
        real(kind=8) :: work(*)
        integer(kind=8) :: info
    end subroutine as_dlaexc
end interface
