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
#include "asterf_types.h"
interface
    subroutine projOrthoNewton(cellCode, cellNbNode, cellDime, cellCoor, poinCoor, &
                               newtIterMaxi, newtToleMaxi, &
                               ksi1, ksi2, &
                               tang_1, tang_2, &
                               projError, lLineSearch_)
        character(len=8), intent(in) :: cellCode
        integer(kind=8), intent(in) :: cellNbNode, cellDime
        real(kind=8), intent(in) :: cellCoor(27), poinCoor(3)
        integer(kind=8), intent(in) :: newtIterMaxi
        real(kind=8), intent(in) :: newtToleMaxi
        real(kind=8), intent(out) :: ksi1, ksi2
        real(kind=8), intent(out) :: tang_1(3), tang_2(3)
        integer(kind=8), intent(out) :: projError
        aster_logical, intent(in), optional :: lLineSearch_
    end subroutine projOrthoNewton
end interface
