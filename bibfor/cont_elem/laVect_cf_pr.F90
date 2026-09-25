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
subroutine laVect_cf_pr(parameters, geom, vect_cont, vect_fric)
!
    use contact_module
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getInterCont.h"
#include "asterfort/getQuadCont.h"
#include "asterfort/laElemCont.h"
#include "blas/dgemv.h"
#include "contact_module.h"
#include "asterfort/laVect_orphan.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(in) :: geom
    real(kind=8), intent(inout) :: vect_cont(MAX_LAGA_DOFS), vect_fric(MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute vector
!
! -----
!
! The implementation follows Poulios' & Renard's paper : An unconstrained integral approximation
! of large sliding frictional contact between deformable solids. Computers & Structures, 2015
! The Lagrange multipliers are expressed in *global coordinates*!
!
! --------------------------------------------------------------------------------------------------
!
! In  parameters       : numerical parameters
! In  geom             : geometrical information
! IO  vect_cont        : vector for contact
! IO  vect_fric        : vector for friction
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_cont_qp, l_fric_qp
    integer(kind=8) :: i_qp, nb_qp, nbPoinInte
    real(kind=8) :: weight_sl_qp, coeff, hF
    real(kind=8) :: coor_qp_sl(2), norm_slav(3)
    real(kind=8) :: coor_qp(2, MAX_NB_QUAD), weight_qp(MAX_NB_QUAD)
    real(kind=8) :: gap, gamma_c, projRmVal, lagr_g(3), thres
    real(kind=8) :: gamma_f, projBsVal(3)
    real(kind=8) :: jump_v(MAX_LAGA_DOFS, 3), mu_g(MAX_LAGA_DOFS, 3)
    real(kind=8) :: poinInteSlav(2, MAX_NB_INTE)
    blas_int :: b_1, b_nb_dofs, b_dime
    blas_int :: b_MAX_LAGA_DOFS
!
! --------------------------------------------------------------------------------------------------
!
    b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
    b_1 = 1
!
    vect_cont = 0.d0
    vect_fric = 0.d0
!
! - Slave node is not paired -> Special treatment
!
    if (geom%elem_slav_code == "PO1") then
        call laVect_orphan(parameters, geom, vect_cont, vect_fric)
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
                        gamma_f, l_fric_qp, &
                        projBsVal=projBsVal, jump_v=jump_v, &
                        norm_slav=norm_slav, lagr_g=lagr_g, &
                        mu_g=mu_g, thres=thres)
        ! if (present(k_diff)) write (6, *) '*projBsVal*', k_diff(1:3), '*', projBsVal(1:3)
!
! ------ CONTACT PART (always computed)
!
        if (l_cont_qp) then
!
! ------ Compute displacement (slave and master side)
!        term: -(lagr_g, -v^m + v^s)
!
            coeff = -weight_sl_qp
            b_nb_dofs = to_blas_int(geom%nb_dofs)
            b_dime = to_blas_int(geom%elem_dime)
            call dgemv('N', b_nb_dofs, b_dime, coeff, jump_v, &
                       b_MAX_LAGA_DOFS, lagr_g, b_1, 1.d0, vect_cont, &
                       b_1)
        end if
!
! ------ Compute Lagrange (slave side)
!        term: (([lagr_g.n^s + gamma_c * gap(u)]_R- n^s - lagr_g) / gamma_c, mu_g)
!
        coeff = weight_sl_qp/gamma_c
        b_nb_dofs = to_blas_int(geom%nb_dofs)
        b_dime = to_blas_int(geom%elem_dime)
        ! write (6, *) 'gap=', gap
        ! write (6, *) 'lagr_g=', lagr_g
        ! write (6, *) 'projRmVal=', projRmVal
        call dgemv('N', b_nb_dofs, b_dime, coeff, mu_g, &
                   b_MAX_LAGA_DOFS, projRmVal*norm_slav-lagr_g, b_1, 1.d0, vect_cont, &
                   b_1)
!
! ------ FRICTION PART (computed only if friction)
!
        if (parameters%l_fric) then
!
! ------ Compute Lagrange (slave side)
!        term: ([lagr_g - gamma_f * v(u)]_Bs / gamma_f, mu_g)
!
            coeff = weight_sl_qp/gamma_f
            b_nb_dofs = to_blas_int(geom%nb_dofs)
            b_dime = to_blas_int(geom%elem_dime)
            call dgemv('N', b_nb_dofs, b_dime, coeff, mu_g, &
                       b_MAX_LAGA_DOFS, projBsVal, b_1, 1.d0, vect_fric, &
                       b_1)
        end if
    end do
!
999 continue
!
!
end subroutine
