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
    subroutine apnorm(elem_nbnode, elem_code, elem_dime, elem_coor,&
                      ksi1       , ksi2     , elem_norm, elem_tau1, elem_tau2, elem_noor)
#include "MeshTypes_type.h"
        integer(kind=8), intent(in) :: elem_nbnode
        character(len=8), intent(in) :: elem_code
        integer(kind=8), intent(in) :: elem_dime
        real(kind=8), intent(in) :: elem_coor(3,MT_NNOMAX2D)
        real(kind=8), intent(in) :: ksi1
        real(kind=8), intent(in) :: ksi2
        real(kind=8), intent(out) :: elem_norm(3)
        real(kind=8), intent(out), optional :: elem_tau1(3), elem_tau2(3), elem_noor
    end subroutine apnorm
end interface
