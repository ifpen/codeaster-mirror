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
subroutine mmnewt(cellCode, cellNbNode, cellDime, cellCoor, poinCoor, &
                  newtIterMaxi, newtToleMaxi, &
                  ksi1, ksi2, &
                  tang_1, tang_2, &
                  projError, lLineSearch_, lPrintDbg_, lCurvature_)
!
    implicit none
!
#include "asterc/r8gaem.h"
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/mmfonf.h"
#include "asterfort/mmreli.h"
#include "asterfort/mmtang.h"
!
    character(len=8), intent(in) :: cellCode
    integer(kind=8), intent(in) :: cellNbNode, cellDime
    real(kind=8), intent(in) :: cellCoor(27)
    real(kind=8), intent(in) :: poinCoor(3)
    integer(kind=8), intent(in) :: newtIterMaxi
    real(kind=8), intent(in) :: newtToleMaxi
    real(kind=8), intent(out) :: ksi1, ksi2
    real(kind=8), intent(out) :: tang_1(3), tang_2(3)
    integer(kind=8), intent(out) :: projError
    aster_logical, intent(in), optional :: lLineSearch_, lPrintDbg_, lCurvature_
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
! In  cellDime    : dimension of space (2 or 3)
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
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8) :: ino, idim, iter
    real(kind=8) :: ff(9), dff(2, 9), ddff(3, 9)
    real(kind=8) :: vect_posi(3)
    real(kind=8) :: matrix(2, 2), par11(3), par12(3), par22(3)
    real(kind=8) :: residu(2)
    real(kind=8) :: dksi1, dksi2
    real(kind=8) :: det, test, refe
    real(kind=8) :: alpha
    real(kind=8) :: tole_rela, tole_abso, tole_newt
    real(kind=8) :: dist, dist_mini, ksi1Min, ksi2Min
    aster_logical :: lLineSearch, lPrintDbg, lCurvature
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(cellNbNode .le. 9)
    ASSERT(cellDime .le. 3)
    ASSERT(cellDime .ge. 2)
    lLineSearch = ASTER_TRUE
    if (present(lLineSearch_)) then
        lLineSearch = lLineSearch_
    end if
    lPrintDbg = ASTER_FALSE
    if (present(lPrintDbg_)) then
        lPrintDbg = lPrintDbg_
    end if
    if (present(lCurvature_)) then
        lCurvature = lCurvature_
    else
        lCurvature = cellCode .ne. 'QU4' .and. cellCode .ne. 'QU8' .and. cellCode .ne. 'QU9'
        lCurvature = ASTER_TRUE
    end if

! - Initializations
    projError = 0
    ksi1 = zero
    ksi2 = zero
    iter = 0
    tole_abso = newtToleMaxi/100.d0
    tole_rela = newtToleMaxi
    alpha = 1.d0
    dist_mini = r8gaem()

! - Newton loop
20  continue
!
    vect_posi = zero
    tang_1 = zero
    tang_2 = zero
    matrix = zero
    residu = zero
    dksi1 = zero
    dksi2 = zero

! - Shape functions (and derivates) at current point
    call mmfonf(cellDime, cellNbNode, cellCode, ksi1, ksi2, &
                ff, dff, ddff)

! - Position vector of current point
    do idim = 1, 3
        do ino = 1, cellNbNode
            vect_posi(idim) = cellCoor(3*(ino-1)+idim)*ff(ino)+vect_posi(idim)
        end do
    end do

! - Local base
    call mmtang(cellDime, cellNbNode, cellCoor, dff, tang_1, &
                tang_2)

! - Quantity to minimize
    do idim = 1, 3
        vect_posi(idim) = poinCoor(idim)-vect_posi(idim)
    end do
    dist = sqrt(vect_posi(1)*vect_posi(1)+vect_posi(2)*vect_posi(2)+vect_posi(3)*vect_posi(3))

! - Newton residual
    residu(1) = vect_posi(1)*tang_1(1)+vect_posi(2)*tang_1(2)+vect_posi(3)*tang_1(3)
    if (cellDime .eq. 3) then
        residu(2) = vect_posi(1)*tang_2(1)+vect_posi(2)*tang_2(2)+vect_posi(3)*tang_2(3)
    end if

! - Local curvatures
    par11 = zero
    par12 = zero
    par22 = zero
    if (lCurvature) then
        do idim = 1, cellDime
            do ino = 1, cellNbNode
                par11(idim) = cellCoor(3*(ino-1)+idim)*ddff(1, ino)+par11(idim)
                if (cellDime .eq. 3) then
                    par22(idim) = cellCoor(3*(ino-1)+idim)*ddff(2, ino)+par22(idim)
                    par12(idim) = cellCoor(3*(ino-1)+idim)*ddff(3, ino)+par12(idim)
                end if
            end do
        end do
    end if

! - Tangent matrix (Newton)
    do idim = 1, 3
        matrix(1, 1) = -tang_1(idim)*tang_1(idim)+par11(idim)*vect_posi(idim)+matrix(1, 1)
        if (cellDime .eq. 3) then
            matrix(1, 2) = -tang_2(idim)*tang_1(idim)+par12(idim)*vect_posi(idim)+matrix(1, 2)
            matrix(2, 1) = -tang_1(idim)*tang_2(idim)+par12(idim)*vect_posi(idim)+matrix(2, 1)
            matrix(2, 2) = -tang_2(idim)*tang_2(idim)+par22(idim)*vect_posi(idim)+matrix(2, 2)
        end if
    end do

! - System determinant
    if (cellDime .eq. 2) then
        det = matrix(1, 1)
    else if (cellDime .eq. 3) then
        det = matrix(1, 1)*matrix(2, 2)-matrix(1, 2)*matrix(2, 1)
    end if
!
    if (abs(det) .le. r8prem()) then
        projError = 1
        goto 999
    end if

! - Solve system
    if (cellDime .eq. 2) then
        dksi1 = -residu(1)/matrix(1, 1)
        dksi2 = 0.d0
    else if (cellDime .eq. 3) then
        dksi1 = (matrix(2, 2)*(-residu(1))-matrix(1, 2)*(-residu(2)))/det
        dksi2 = (matrix(1, 1)*(-residu(2))-matrix(2, 1)*(-residu(1)))/det
    else
        ASSERT(ASTER_FALSE)
    end if

! - Line search
    if (lLineSearch) then
        call mmreli(cellCode, cellNbNode, cellDime, cellCoor, poinCoor, &
                    ksi1, ksi2, dksi1, dksi2, alpha)
    else
        alpha = 1.d0
    end if

! - Update
    ksi1 = ksi1+alpha*dksi1
    ksi2 = ksi2+alpha*dksi2

! - Save values if Newton avoids
    if (dist .le. dist_mini) then
        dist_mini = dist
        ksi1Min = ksi1
        ksi2Min = ksi2
    end if

! - Convergence
    refe = (ksi1*ksi1+ksi2*ksi2)
    if (refe .le. tole_rela) then
        tole_newt = tole_abso
        test = sqrt(dksi1*dksi1+dksi2*dksi2)
    else
        tole_newt = tole_rela
        test = sqrt(dksi1*dksi1+dksi2*dksi2)/sqrt(refe)
    end if

! - Continue or not ?
    if ((test .gt. tole_newt) .and. (iter .lt. newtIterMaxi)) then
        iter = iter+1
        goto 20
    else if ((iter .ge. newtIterMaxi) .and. (test .gt. tole_newt)) then
        ksi1 = ksi1Min
        ksi2 = ksi2Min
        call mmfonf(cellDime, cellNbNode, cellCode, ksi1, ksi2, &
                    ff, dff, ddff)
        call mmtang(cellDime, cellNbNode, cellCoor, dff, tang_1, &
                    tang_2)
        projError = 1
    end if

! - End of loop
999 continue
!
    if (projError .eq. 1 .and. lPrintDbg) then
        WRITE (6, *) "Courbure ?", lCurvature
        write (6, *) 'POINT A PROJETER : ', poinCoor(1), poinCoor(2), poinCoor(3)
        write (6, *) 'MAILLE             ', cellCode, cellNbNode
        do ino = 1, cellNbNode
            write (6, *) '  NOEUD ', ino
            write (6, *) '   (X,Y,Z)', cellCoor(3*(ino-1)+1), &
                cellCoor(3*(ino-1)+2), &
                cellCoor(3*(ino-1)+3)
        end do
        write (6, *) 'KSI   : ', ksi1, ksi2
        write (6, *) 'ALPHA : ', alpha
    end if
!
end subroutine
