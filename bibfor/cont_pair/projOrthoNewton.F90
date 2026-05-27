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
subroutine projOrthoNewton(cellCode, cellNbNode, cellDime, cellCoor, poinCoor, &
                           newtIterMaxi, newtToleMaxi, &
                           ksi1, ksi2, &
                           tang_1, tang_2, &
                           projError, lLineSearch_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/mmnewt.h"
#include "asterfort/elrfno.h"
!
    character(len=8), intent(in) :: cellCode
    integer(kind=8), intent(in) :: cellNbNode, cellDime
    real(kind=8), intent(in) :: cellCoor(27), poinCoor(3)
    integer(kind=8), intent(in) :: newtIterMaxi
    real(kind=8), intent(in) :: newtToleMaxi
    real(kind=8), intent(out) :: ksi1, ksi2
    real(kind=8), intent(out) :: tang_1(3), tang_2(3)
    integer(kind=8), intent(out) :: projError
    aster_logical, intent(in), optional :: lLineSearch_
!
! --------------------------------------------------------------------------------------------------
!
! Contact (all methods)
!
! Projection of point on element (Newton algorithm) - Minimum distance
!
! --------------------------------------------------------------------------------------------------
!
! In  cellCode : element type
! In  cellNbNode   : number of nodes of element
! In  cellDime    : dimension of element (2 or 3)
! In  cellCoor : coordinates of nodes of the element
! In  poinCoor   : coordinates of poitn to project
! In  newtIterMaxi : Newton algorithm - Maximum number of iterations
! In  newtToleMaxi : Newton algorithm - Tolerance
! Out ksi1      : first parametric coordinate of projection of point on element
! Out ksi2      : second parametric coordinate of projection of point on element
! Out tang_1    : first tangent of local basis for the projection of point on element
! Out tang_2    : second tangent of local basis for the projection of point on element
! Out projError     : projError code
!                  0  OK
!                  1  NON-OK
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: cellLineNbNode, nbNode, nbNodeS
    aster_logical :: lLineSearch, lPrintDbg, lReproject, lCurvature
    character(len=8) :: cellLineCode
!
! --------------------------------------------------------------------------------------------------
!
    projError = 0
    lLineSearch = ASTER_TRUE
    if (present(lLineSearch_)) then
        lLineSearch = lLineSearch_
    end if
    call elrfno(cellCode, nbNode, nbNodeS)
    lReproject = nbNode .ne. nbNodeS

    lPrintDbg = ASTER_TRUE
    if (lReproject) then
        lPrintDbg = ASTER_FALSE
    end if
    call mmnewt(cellCode, cellNbNode, cellDime, cellCoor, poinCoor, &
                newtIterMaxi, newtToleMaxi, &
                ksi1, ksi2, &
                tang_1, tang_2, &
                projError, lLineSearch, lPrintDbg)
    if (lReproject .and. projError .ne. 0) then
        lPrintDbg = ASTER_TRUE
        lCurvature = ASTER_FALSE
        call mmnewt(cellCode, cellNbNode, cellDime, cellCoor, poinCoor, &
                    newtIterMaxi, newtToleMaxi, &
                    ksi1, ksi2, &
                    tang_1, tang_2, &
                    projError, lLineSearch, lPrintDbg, lCurvature)
    end if
!
end subroutine
