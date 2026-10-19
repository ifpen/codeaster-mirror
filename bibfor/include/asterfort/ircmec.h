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
    subroutine ircmec(idfimd, nochmd, nomprf, nolopg, numpt,&
                      instan, numord, val, ncmpve, nbenty,&
                      nbrepg, nvalec, typent, typgeo, nosdfu,&
                      tymast, nrimpr, nosdpr, codret)
        med_idt :: idfimd
        character(len=*) :: nochmd
        character(len=*) :: nomprf
        character(len=*) :: nolopg
        integer(kind=8) :: numpt
        real(kind=8) :: instan
        integer(kind=8) :: numord
        real(kind=8) :: val(*)
        integer(kind=8) :: ncmpve
        integer(kind=8) :: nbenty
        integer(kind=8) :: nbrepg
        integer(kind=8) :: nvalec
        integer(kind=8) :: typent
        integer(kind=8) :: typgeo
        character(len=8) :: nosdfu
        integer(kind=8) :: tymast
        integer(kind=8) :: nrimpr
        character(len=*) :: nosdpr
        integer(kind=8) :: codret
    end subroutine ircmec
end interface
