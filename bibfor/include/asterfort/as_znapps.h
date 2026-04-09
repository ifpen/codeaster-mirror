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
    subroutine as_znapps(n, kev, np, shift, v,&
                      ldv, h, ldh, resid, q,&
                      ldq, workl, workd)
        integer(kind=8) :: ldq
        integer(kind=8) :: ldh
        integer(kind=8) :: ldv
        integer(kind=8) :: np
        integer(kind=8) :: kev
        integer(kind=8) :: n
        complex(kind=8) :: shift(np)
        complex(kind=8) :: v(ldv, kev+np)
        complex(kind=8) :: h(ldh, kev+np)
        complex(kind=8) :: resid(n)
        complex(kind=8) :: q(ldq, kev+np)
        complex(kind=8) :: workl(kev+np)
        complex(kind=8) :: workd(2*n)
    end subroutine as_znapps
end interface
