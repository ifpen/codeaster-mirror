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
interface
    subroutine nueffe(nbLigr, listLigr, base, numeDofZ, renumZ, &
                      modelZ, modeLocZ_, idenRelaZ_, lMacrElemZ_)
        integer(kind=8), intent(in) :: nbLigr
        character(len=24), pointer :: listLigr(:)
        character(len=2), intent(in) :: base
        character(len=*), intent(in) :: numeDofZ, renumZ, modelZ
        character(len=*), optional, intent(in) :: modeLocZ_, idenRelaZ_
        aster_logical, optional, intent(in) :: lMacrElemZ_
    end subroutine nueffe
end interface
