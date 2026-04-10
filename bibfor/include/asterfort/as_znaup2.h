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
    subroutine as_znaup2(ido, bmat, n, which, nev,&
                      np, tol, resid, ishift, mxiter,&
                      v, ldv, h, ldh, ritz,&
                      bounds, q, ldq, workl, ipntr,&
                      workd, rwork, info, neqact, alpha)
        integer(kind=8) :: ldq
        integer(kind=8) :: ldh
        integer(kind=8) :: ldv
        integer(kind=8) :: np
        integer(kind=8) :: nev
        integer(kind=8) :: n
        integer(kind=8) :: ido
        character(len=1) :: bmat
        character(len=2) :: which
        real(kind=8) :: tol
        complex(kind=8) :: resid(n)
        integer(kind=8) :: ishift
        integer(kind=8) :: mxiter
        complex(kind=8) :: v(ldv, nev+np)
        complex(kind=8) :: h(ldh, nev+np)
        complex(kind=8) :: ritz(nev+np)
        complex(kind=8) :: bounds(nev+np)
        complex(kind=8) :: q(ldq, nev+np)
        complex(kind=8) :: workl((nev+np)*(nev+np+3))
        integer(kind=8) :: ipntr(13)
        complex(kind=8) :: workd(3*n)
        real(kind=8) :: rwork(nev+np)
        integer(kind=8) :: info
        integer(kind=8) :: neqact
        real(kind=8) :: alpha
    end subroutine as_znaup2
end interface
