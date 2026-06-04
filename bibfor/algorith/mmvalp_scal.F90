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
subroutine mmvalp_scal(cellCode, cellNbNode, ksi1, ksi2, &
                       valeCell, valePoin)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/mmnonf.h"
!
    character(len=8), intent(in) :: cellCode
    integer(kind=8), intent(in) :: cellNbNode
    real(kind=8), intent(in) :: ksi1
    real(kind=8), intent(in) :: ksi2
    real(kind=8), intent(in) :: valeCell(*)
    real(kind=8), intent(out) :: valePoin
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Utility
!
! Continue method - Interpolate ONE component at point in given element
!
! --------------------------------------------------------------------------------------------------
!
! In  cellCode         : type of element
! In  cellNbNode       : number of nodes
! In  ksi1             : first parametric coordinate of the point
! In  ksi2             : second parametric coordinate of the point
! In  valeCell        : value of components at nodes
! Out valePoin        : value of components at point
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: shape_func(9)
    integer(kind=8) :: i_node
!
! --------------------------------------------------------------------------------------------------
!
    valePoin = 0.d0
    ASSERT(cellNbNode .le. 9)

! - Shape functions
    call mmnonf(cellCode, ksi1, ksi2, shape_func)

! - Compute
    do i_node = 1, cellNbNode
        valePoin = shape_func(i_node)*valeCell(i_node)+valePoin
    end do
!
end subroutine
