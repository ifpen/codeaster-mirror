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

subroutine te0524(nomopt, nomte)
!
    use FE_algebra_module
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/readVector.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/writeVector.h"
!
    character(len=16) :: nomte, nomopt
!
!===================================================================================================
!
!
! Mechanics - Coupling FEM/FEM with augmented Lagrangian - Set LAGR to zero.
!
!===================================================================================================
!
    integer(kind=8) :: ndim, i, nbDoFs
    real(kind=8) :: rhs(6), disp_prev(6), disp_curr(6)
    real(kind=8) :: lhs(6, 6)
!
    if (nomte == "CL_POI2D") then
        ndim = 2
    elseif (nomte == "CL_POI3D") then
        ndim = 3
    else
        ASSERT(ASTER_FALSE)
    end if
    nbDoFs = 2*ndim
!
! - Computation
!
    if (nomopt == "CHAR_MECA_CPL") then
!
! --- Compute coupling residual
!
        call readVector('PDEPLMR', nbDoFs, disp_prev)
        ! this is incremental displacement
        call readVector('PDEPLPR', nbDoFs, disp_curr)
        call daxpy_1(nbDoFs, 1.d0, disp_prev, disp_curr)
        rhs = 0.d0
        do i = ndim+1, nbDoFs
            rhs(i) = disp_curr(i)
        end do
!
! --- Write vector
!
        call writeVector('PVECTUR', nbDoFs, rhs)
!
    elseif (nomopt == "RIGI_CPL" .or. nomopt == "RIGI_ELAS_CPL") then
!
! --- Compute coupling matrix
!
        lhs = 0.d0
!
        do i = ndim+1, nbDoFs
            lhs(i, i) = 1.d0
        end do
!
! - Write matrix
!
        call writeMatrix("PMATUUR", nbDoFs, nbDoFs, ASTER_TRUE, lhs)
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
