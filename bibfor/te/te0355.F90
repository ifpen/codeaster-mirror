! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

subroutine te0355(nomopt, nomte)
!
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/laelem.h"
#include "asterfort/laMatr.h"
#include "asterfort/laMatr_diff.h"
#include "asterfort/laParam.h"
#include "asterfort/laVect.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/writeVector.h"
#include "contact_module.h"
#include "jeveux.h"
!
    character(len=16) :: nomte, nomopt
!
! --------------------------------------------------------------------------------------------------
!
!   CHAR_MECA_CONT and RIGI_CONT for augmented Lagrangian method (Mortar)
!
! --------------------------------------------------------------------------------------------------
!
    type(ContactParameters) :: parameters
    type(ContactGeom) :: geom
    real(kind=8) :: vect_cont(MAX_LAGA_DOFS), vect_fric(MAX_LAGA_DOFS)
    real(kind=8), pointer :: matr_cont(:, :) => null()
    real(kind=8), pointer :: matr_fric(:, :) => null()
    aster_logical :: diff_num, l_vari
!
! --------------------------------------------------------------------------------------------------
!

! - Large arrays allocated on the heap rather than on the stack
    allocate (matr_cont(MAX_LAGA_DOFS, MAX_LAGA_DOFS))
    allocate (matr_fric(MAX_LAGA_DOFS, MAX_LAGA_DOFS))

! - Informations about finite element
    call laelem(nomte, geom, parameters)

! - Get Parameters
    call laParam(parameters)

    if ((parameters%vari_cont == CONT_VARI_ROBU) .and. &
        (parameters%l_fric)) then
        ASSERT(parameters%coef_cont(1) == parameters%coef_fric(1))
    end if

! - Use finite difference Jacobian
    diff_num = (parameters%jac_type .gt. 0)

! - Verif
    ASSERT(parameters%algo_cont == CONT_ALGO_LAGR)
    l_vari = ((parameters%vari_cont == CONT_VARI_ROBU) .or. &
              (parameters%vari_cont == CONT_VARI_CLAS) .or. &
              (parameters%vari_cont == CONT_VARI_NONE))
    ASSERT(l_vari)
!
! - Computation
!
    if (nomopt == "CHAR_MECA_CONT") then
!
! --- Compute contact residual
!
        call laVect(parameters, geom, vect_cont, vect_fric)
!
! --- Write vector
!
        call writeVector('PVECTCR', geom%nb_dofs, vect_cont)
!
        if (parameters%l_fric) then
            call writeVector('PVECTFR', geom%nb_dofs, vect_fric)
        end if
!
    elseif (nomopt == "RIGI_CONT") then
!
! --- Compute contact matrix
!
        if (diff_num) then
            call laMatr_diff(parameters, geom, matr_cont, matr_fric)
        else
            call laMatr(parameters, geom, matr_cont, matr_fric)
        end if
!
! --- Write matrix
!
        if (.not. parameters%l_fric) then
            ! Contact only
            call writeMatrix('PMATUNS', geom%nb_dofs, geom%nb_dofs, ASTER_FALSE, matr_cont)
        else
            ! Contact with friction
            call writeMatrix('PMATUNS', geom%nb_dofs, geom%nb_dofs, ASTER_FALSE, matr_cont)
!
            call writeMatrix('PMATUNS', geom%nb_dofs, geom%nb_dofs, ASTER_FALSE, matr_fric)
        end if
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
    deallocate (matr_cont)
    deallocate (matr_fric)
!
end subroutine
