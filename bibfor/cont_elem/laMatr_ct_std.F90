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
!
subroutine laMatr_ct_std(parameters, geom, matr_cont, matr_fric)
!
    use contact_type
    use contact_module
    use contact_algebra_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getInterCont.h"
#include "asterfort/getQuadCont.h"
#include "asterfort/laElemCont.h"
#include "blas/dger.h"
#include "contact_module.h"
#include "asterfort/laMatr_orphan.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(in) :: geom
    real(kind=8), intent(inout) :: matr_cont(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
    real(kind=8), intent(inout) :: matr_fric(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute matrices
!
! WARNING
! Note that we build the transpose of the Jacobian - it allows the computer expressions to be
! closer to the mathematical ones. The matrix is then transposed at the end of the subroutine.
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
    aster_logical :: l_cont_qp, l_fric_qp
    integer(kind=8) :: i_qp, nb_qp, nbPoinInte
    real(kind=8) :: weight_sl_qp, coeff, hF
    real(kind=8) :: coor_qp_sl(2)
    real(kind=8) :: coor_qp(2, MAX_NB_QUAD), weight_qp(MAX_NB_QUAD)
    real(kind=8) :: gap, gamma_c, projRmVal
    real(kind=8) :: gamma_f
    real(kind=8) :: dGap(MAX_LAGA_DOFS)
    real(kind=8) :: lagr_c, mu_c(MAX_LAGA_DOFS)
    real(kind=8) :: d2Gap(MAX_LAGA_DOFS, MAX_LAGA_DOFS)
    real(kind=8) :: poinInteSlav(2, MAX_NB_INTE)
    blas_int :: b_nb_dofs
    blas_int :: b_1, b_3, b_MAX_LAGA_DOFS
!
! --------------------------------------------------------------------------------------------------
!
    b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
    b_1 = 1
    b_3 = 3
!
    matr_cont = 0.d0
    matr_fric = 0.d0
!
! - Slave node is not paired -> Special treatment
!
    if (geom%elem_slav_code == "PO1") then
        call laMatr_orphan(parameters, geom, matr_cont, matr_fric)
!
        go to 999
    end if
!
! - Get intersection points
!
    call getInterCont(nbPoinInte, poinInteSlav)
!
! - Get quadrature (slave side)
!
    call getQuadCont(geom%elem_dime, &
                     geom%elem_slav_code, geom%elem_mast_code, &
                     nbPoinInte, poinInteSlav, &
                     nb_qp, coor_qp, &
                     geom%l_axis, geom%nb_node_slav, geom%coor_slav_init, &
                     weight_qp)
!
! - Diameter of slave side
!
    hF = diameter(geom%nb_node_slav, geom%coor_slav_init)
!
! - Loop on quadrature points
!
    do i_qp = 1, nb_qp
!
! ----- Get current quadrature point (slave side)
!
        coor_qp_sl(1:2) = coor_qp(1:2, i_qp)
        weight_sl_qp = weight_qp(i_qp)
!
! ----- Compute contact quantities
!
        call laElemCont(parameters, geom, coor_qp_sl, hF, &
                        gap, gamma_c, projRmVal, l_cont_qp, &
                        gamma_f, l_fric_qp, lagr_c=lagr_c, &
                        dGap=dGap, d2Gap=d2Gap, mu_c=mu_c)
!
!
! ------ CONTACT PART
!

        if (l_cont_qp) then
!
! ------ Compute displacement / displacement (slave and master side)
!        term: (gamma_c*H*D(gap(u))[v], D(gap(u))[du])
!
            coeff = weight_sl_qp*gamma_c
            b_nb_dofs = to_blas_int(geom%nb_dofs)
            call dger(b_nb_dofs, b_nb_dofs, coeff, dGap, b_1, &
                      dGap, b_1, matr_cont, b_MAX_LAGA_DOFS)

!
! ------ Compute displacement / displacement (slave and master side)
!        term: (H*[lagr_c + gamma_c * gap(u)]_R-, D2(gap(u))[v, du])
!
            coeff = weight_sl_qp*projRmVal
            matr_cont(1:geom%nb_dofs, 1:geom%nb_dofs) = matr_cont(1:geom%nb_dofs, 1:geom%nb_dofs)&
                                                        &+coeff*d2Gap(1:geom%nb_dofs, 1:geom%nb_d&
                                                                  &ofs)
!
! ------ Compute displacement / Lagrange and Lagrange / displacement
!        term: (H * D(gap(u))[v], dlagr_c) + (H * mu_c,  D(gap(u))[du])
!
            coeff = weight_sl_qp
            b_nb_dofs = to_blas_int(geom%nb_dofs)
            call dger(b_nb_dofs, b_nb_dofs, coeff, dGap, b_1, &
                      mu_c, b_1, matr_cont, b_MAX_LAGA_DOFS)
            call dger(b_nb_dofs, b_nb_dofs, coeff, mu_c, b_1, &
                      dGap, b_1, matr_cont, b_MAX_LAGA_DOFS)
        else
!
! ------ Compute Lagrange / Lagrange (slave side)
!        term: ((H-1) / gamma_c * mu_c, dlagr_c) = (- mu_c / gamma_c, dlagr_c) since H = 0
!
            coeff = -weight_sl_qp/gamma_c
            b_nb_dofs = to_blas_int(geom%nb_dofs)
            call dger(b_nb_dofs, b_nb_dofs, coeff, mu_c, b_1, &
                      mu_c, b_1, matr_cont, b_MAX_LAGA_DOFS)
!
        end if
    end do
!
999 continue
!
end subroutine
