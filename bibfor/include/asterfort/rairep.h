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
    subroutine rairep(noma, ioc, km, rigi, nbgr, &
                      ligrma, zjdlm, nbno, tabnoe, rignoe, rirot, ndim, jdme)
        integer(kind=8)   :: nbgr
        character(len=8)  :: noma
        integer(kind=8)   :: ioc
        character(len=8)  :: km
        real(kind=8)      :: rigi(6)
        character(len=24) :: ligrma(nbgr)
        integer(kind=8)   :: zjdlm(*)
        integer(kind=8)   :: nbno
        character(len=8)  :: tabnoe(*)
        real(kind=8)      :: rignoe(*)
        real(kind=8)      :: rirot(3)
        integer(kind=8)   :: ndim
        integer(kind=8)   :: jdme
    end subroutine rairep
end interface
