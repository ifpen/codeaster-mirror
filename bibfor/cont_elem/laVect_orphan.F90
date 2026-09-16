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

subroutine laVect_orphan(parameters, geom, vect_cont, vect_fric)
!
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "contact_module.h"
#include "asterfort/assert.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(in) :: geom
    real(kind=8), intent(inout) :: vect_cont(MAX_LAGA_DOFS), vect_fric(MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Handle the orphan nodes in the tangent matrix = nodes that do not belong to any virtual cell
!
! --------------------------------------------------------------------------------------------------
!
! In  parameters       : numerical parameters
! In  geom             : geometrical information
! IO  matr_cont        : matrix (only upper part)
! IO  matr_fric        : matrix
!
! --------------------------------------------------------------------------------------------------
!
    if (geom%elem_mast_code == "LAGR") then
        vect_cont(geom%elem_dime+1) = geom%lagc_slav_curr(1)
        if (parameters%l_fric) then
            vect_fric(geom%elem_dime+2) = geom%lagf_slav_curr(1, 1)
            if (geom%elem_dime == 3) then
                vect_fric(geom%elem_dime+3) = geom%lagf_slav_curr(2, 1)
            end if
        end if
    else
        if (geom%elem_mast_code .ne. "NOLAGR") then
            ASSERT(ASTER_FALSE)
        end if
    end if
!
end subroutine
