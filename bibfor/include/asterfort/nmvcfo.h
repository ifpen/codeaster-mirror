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
    subroutine nmvcfo(poum, model, materField, materCode, caraElem, compor, &
                      varcRefe, hvalIncr, vectElem)
        character(len=1), intent(in) :: poum
        character(len=24), intent(in) :: model
        character(len=24), intent(in) :: materField, materCode, caraElem
        character(len=24), intent(in) :: varcRefe
        character(len=24), intent(in) :: compor
        character(len=19), intent(in) :: hvalIncr(*)
        character(len=19), intent(in) :: vectElem
    end subroutine nmvcfo
end interface
