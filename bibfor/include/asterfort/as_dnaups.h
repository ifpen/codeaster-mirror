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
    subroutine as_dnaups(ido, bmat, n, which, nev,&
                      tol, resid, ncv, v, ldv,&
                      iparam, ipntr, workd, workl, lworkl,&
                      info, neqact, alpha, nsta, ddlsta,&
                      vstab, csta, ldynfa, ddlexc, redem)
        integer(kind=8) :: lworkl
        integer(kind=8) :: ldv
        integer(kind=8) :: ncv
        integer(kind=8) :: n
        integer(kind=8) :: ido
        character(len=1) :: bmat
        character(len=2) :: which
        integer(kind=8) :: nev
        real(kind=8) :: tol
        real(kind=8) :: resid(n)
        real(kind=8) :: v(ldv, ncv)
        integer(kind=8) :: iparam(11)
        integer(kind=8) :: ipntr(14)
        real(kind=8) :: workd(3*n)
        real(kind=8) :: workl(lworkl)
        integer(kind=8) :: info
        integer(kind=8) :: neqact
        real(kind=8) :: alpha
        integer(kind=8) :: nsta
        integer(kind=8) :: ddlsta(n)
        real(kind=8) :: vstab(n)
        real(kind=8) :: csta
        integer(kind=8) :: ldynfa
        integer(kind=8) :: ddlexc(n)
        integer(kind=8) :: redem
    end subroutine as_dnaups
end interface
