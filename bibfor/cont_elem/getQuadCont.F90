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
subroutine getQuadCont(elem_dime, &
                       elem_slav_code, elem_mast_code, &
                       nbPoinInte, poinInteSlav, &
                       nb_qp, coor_qp, &
                       l_axis_, nb_node_slav_, elem_slav_coor_, &
                       weight_qp_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/latrco.h"
#include "asterfort/lcptga.h"
#include "asterfort/mesh_pairing_type.h"
#include "asterfort/mmdonf.h"
#include "asterfort/mmmjac.h"
#include "asterfort/mmnonf.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: elem_dime
    character(len=8), intent(in) :: elem_slav_code, elem_mast_code
    integer(kind=8), intent(in) :: nbPoinInte
    real(kind=8), intent(in) :: poinInteSlav(2, MAX_NB_INTE)
    real(kind=8), intent(out) :: coor_qp(2, MAX_NB_QUAD)
    integer(kind=8), intent(out) :: nb_qp
    integer(kind=8), optional, intent(in) :: nb_node_slav_
    real(kind=8), optional, intent(in) :: elem_slav_coor_(3, 9)
    aster_logical, optional, intent(in) :: l_axis_
    real(kind=8), optional, intent(out) :: weight_qp_(MAX_NB_QUAD)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Quadrature
!
! Compute quadrature in slave slide
!
! --------------------------------------------------------------------------------------------------
!
! In  modelDime        : dimension of model
! In  elem_slav_code   : code element for slave side from contact element
! In  elem_slav_coor   : coordinates from slave side of contact element
! In  nbPoinInte       : number of intersection points
! In  poinInteSlav     : coordinates of intersection points (in slave parametric space)
! Out nb_qp            : number of quadrature points
! Out coor_qp          : coordinates of quadrature points (parametric slave space)
! Out weight_qp        : weight of quadrature points
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_slav_line, l_mast_line
    integer(kind=8) :: iTria, iGauss
    integer(kind=8) :: nbTria, nbGauss
    real(kind=8) :: triaCoorSlav(2, 3)
    real(kind=8) :: gausWeightSlav(12), gausCoorSlav(2, 12)
    real(kind=8) :: shape_func(9), shape_dfunc(2, 9), jacobian_sl
    character(len=8) :: elga_fami
!
! --------------------------------------------------------------------------------------------------
!
    nb_qp = 0
    coor_qp = 0
    if (present(weight_qp_)) then
        weight_qp_ = 0.d0
    end if
    l_slav_line = elem_slav_code == "SE2" .or. elem_slav_code == "TR3" .or. elem_slav_code == "QU4"
    l_mast_line = elem_mast_code == "SE2" .or. elem_mast_code == "TR3" .or. elem_mast_code == "QU4"

! - Triangulation of convex polygon defined by intersection points
    if (elem_dime .eq. 3) then
        if (nbPoinInte == 3) then
            nbTria = 1
        else
            nbTria = nbPoinInte
        end if
        if (l_slav_line .and. l_mast_line) then
            ! order 3 by triangle
            elga_fami = 'FPG4'
        else
            ! order 5 by triangle
            elga_fami = 'FPG7'
        end if
    elseif (elem_dime .eq. 2) then
        nbTria = 1
        if (l_slav_line .and. l_mast_line) then
            ! order 5
            elga_fami = 'FPG3'
        else
            ! order 7
            elga_fami = 'FPG4'
        end if
    else
        ASSERT(ASTER_FALSE)
    end if

! - Loop on triangles
    do iTria = 1, nbTria
! ----- Coordinates of current triangle (slave)
        triaCoorSlav = 0.d0
        if (elem_dime .eq. 3) then
            call latrco(iTria, nbPoinInte, poinInteSlav, triaCoorSlav)
        elseif (elem_dime .eq. 2) then
            ASSERT(nbPoinInte .eq. 2)
            triaCoorSlav(1:2, 1:2) = poinInteSlav(1:2, 1:2)
        end if

! ----- Get integration points for slave element
        call lcptga(elem_dime, triaCoorSlav, elga_fami, &
                    nbGauss, gausCoorSlav, gausWeightSlav)

! ----- Loop on integration points
        do iGauss = 1, nbGauss
            nb_qp = nb_qp+1
            ASSERT(nb_qp <= MAX_NB_QUAD)

! --------- Get current integration point (slave)
            coor_qp(1:2, nb_qp) = gausCoorSlav(1:2, iGauss)
            ! WRITE (6, *) "coor_qp: ", coor_qp(:, nb_qp)

            if (present(weight_qp_)) then
! ------------- Get shape functions and first derivative only (for perf)
                call mmnonf(elem_slav_code, &
                            gausCoorSlav(1, iGauss), gausCoorSlav(2, iGauss), &
                            shape_func)
                call mmdonf(elem_slav_code, &
                            gausCoorSlav(1, iGauss), gausCoorSlav(2, iGauss), &
                            shape_dfunc)

! ------------- Compute jacobian
                call mmmjac(l_axis_, nb_node_slav_, elem_dime, &
                            elem_slav_code, elem_slav_coor_, &
                            shape_func, shape_dfunc, jacobian_sl)

                weight_qp_(nb_qp) = jacobian_sl*gausWeightSlav(iGauss)
            end if
        end do
    end do
!
end subroutine
