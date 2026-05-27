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
! aslint: disable=W0413
!
module contact_module
!
    use contact_type
    implicit none
!
    private
!
#include "asterf_types.h"
#include "asterfort/apnorm.h"
#include "asterfort/assert.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/mm2onf.h"
#include "asterfort/mmdonf.h"
#include "asterfort/mmnewd.h"
#include "asterfort/mmnonf.h"
#include "asterfort/projInsideCell.h"
#include "asterfort/reereg.h"
#include "asterfort/reerel.h"
#include "contact_module.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! Contact - generic
!
! Generic method for contact
!
! --------------------------------------------------------------------------------------------------
!
    public :: projQpSl2Ma, shapeFuncDisp, shapeFuncLagr, evalPoly
    public :: diameter, testLagrC, gapEval, testLagrF, barycenter
    public :: speedEval, thresEval, shapeFuncDispVolu, projQpSl2Vo
    public :: gradFuncDispVolu, evalDtestM
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine projQpSl2Ma(geom, coor_qp_sl, proj_tole, coor_qp_ma, gap, &
                           tau_slav, norm_slav, tau_mast, norm_mast)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: coor_qp_sl(2), proj_tole
        real(kind=8), intent(out) :: coor_qp_ma(2)
        real(kind=8), intent(out) :: gap
        real(kind=8), intent(out) :: tau_slav(3, 2), tau_mast(3, 2)
        real(kind=8), intent(out) :: norm_slav(3), norm_mast(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Project slave point to master space and compute quantity
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: iret, iret1, iret2, elem_mast_line_nbnode
        real(kind=8) :: coor_qp_sl_re(3), tau1_mast(3), tau2_mast(3), coor_qp_ma_re(3)
        real(kind=8) :: coor_qp_sl_re_aux(3)
        real(kind=8) :: ksi_line(2)
        character(len=8) :: elem_mast_line_code
        aster_logical :: debug
!
!
        norm_slav = 0.d0
        iret = 0
        iret1 = 0
        iret2 = 0
        norm_mast = 0.d0
        tau_slav = 0.d0
        tau_mast = 0.d0
        gap = 0.d0
        ksi_line = 0.d0
        debug = ASTER_FALSE
!
! ------ Compute outward slave normal (pairing configuration)
!
        call apnorm(geom%nb_node_slav, geom%elem_slav_code, geom%elem_dime, geom%coor_slav_pair, &
                    coor_qp_sl(1), coor_qp_sl(2), norm_slav, tau_slav(1:3, 1), tau_slav(1:3, 2))
!
! ----- Return in real slave space (pairing configuration)
!
        coor_qp_sl_re = 0.d0
        call reerel(geom%elem_slav_code, geom%nb_node_slav, 3, geom%coor_slav_pair, coor_qp_sl, &
                    coor_qp_sl_re)
! ----- Projection of node on master cell (master parametric space)
!
        if (geom%elem_mast_code(1:2) == "SE") then
            elem_mast_line_code = "SE2"
            elem_mast_line_nbnode = 2
        else if (geom%elem_mast_code(1:2) == "TR") then
            elem_mast_line_code = "TR3"
            elem_mast_line_nbnode = 3
        else if (geom%elem_mast_code(1:2) == "QU") then
            elem_mast_line_code = "QU4"
            elem_mast_line_nbnode = 4
        else
            ASSERT(ASTER_FALSE)
        end if
!
!
! ----- Projection on master element
        call mmnewd(geom%elem_mast_code, geom%nb_node_mast, geom%elem_dime, geom%coor_mast_pair, &
                    coor_qp_sl_re, 100, proj_tole, norm_slav, ksi_line(1), &
                    ksi_line(2), tau1_mast, tau2_mast, iret)
        if (iret == 1) then
!
! ----- Try with linearization
!
            call mmnewd(elem_mast_line_code, elem_mast_line_nbnode, geom%elem_dime, &
                        geom%coor_mast_pair, coor_qp_sl_re, 75, proj_tole, norm_slav, &
                        ksi_line(1), ksi_line(2), tau1_mast, tau2_mast, iret1)
!
            call reerel(elem_mast_line_code, elem_mast_line_nbnode, 3, geom%coor_mast_pair, &
                        ksi_line, coor_qp_sl_re_aux)
!
            call mmnewd(geom%elem_mast_code, geom%nb_node_mast, geom%elem_dime, &
                        geom%coor_mast_pair, coor_qp_sl_re_aux, 75, proj_tole, norm_slav, &
                        ksi_line(1), ksi_line(2), tau1_mast, tau2_mast, iret1)
        end if
!
!
!
! ----- Check that projected node is inside cell
!
        call projInsideCell(1e-3, geom%elem_dime, geom%elem_mast_code, ksi_line, iret2)
!
        coor_qp_ma(:) = ksi_line(:)
!
        if ((iret .eq. 0) .and. (iret1 .eq. 0) .and. (iret2 .eq. 0)) then
!print*, 'mast_coor_pair', geom%coor_mast_pair(1:2,1:3)
!print*, 'mast_coor_curr', geom%coor_mast_curr(1:2,1:3)
!print*, 'slav_coor_pair', geom%coor_slav_pair(1:2,1:3)
!print*, 'salv_coor_curr', geom%coor_slav_curr(1:2,1:3)
!print*, "normslav", norm_slav
!print*, "iret", iret
!print*, "iret1", iret1
!print*, "iret2", iret2
            if (debug) then
!print*, "1", ksi_line(1)
!print*, "2", ksi_line(2)
!print*, "3", ksi_line(1)+coor_qp_ma(2)
            end if
        end if
!
! ------ Compute outward master normal (pairing configuration)
!
        call apnorm(geom%nb_node_mast, geom%elem_mast_code, geom%elem_dime, geom%coor_mast_pair, &
                    coor_qp_ma(1), coor_qp_ma(2), norm_mast, tau_mast(1:3, 1), tau_mast(1:3, 2))
!
! ------ Check
!
        if (dot_product(norm_slav, norm_mast) > 0.1d0) then
            if (debug) then
                print *, "Normals have the same direction: ", dot_product(norm_slav, norm_mast)
                print *, "Slave normal: ", norm_slav
                print *, "Master normal: ", norm_mast
                print *, "Distance: ", norm2(barycenter(geom%nb_node_mast, geom%coor_mast_curr)- &
                                             barycenter(geom%nb_node_slav, geom%coor_slav_curr))
                print *, "Diameter slave / master: ", &
                    diameter(geom%nb_node_slav, geom%coor_slav_curr), &
                    diameter(geom%nb_node_mast, geom%coor_mast_curr)
            end if
!
        end if
!
! ----- Compute gap for raytracing gap = -(x^s - x^m).n^s (current configuration)
!
        coor_qp_sl_re = 0.d0
        call reerel(geom%elem_slav_code, geom%nb_node_slav, 3, geom%coor_slav_curr, coor_qp_sl, &
                    coor_qp_sl_re)
        coor_qp_ma_re = 0.d0
        call reerel(geom%elem_mast_code, geom%nb_node_mast, 3, geom%coor_mast_curr, coor_qp_ma, &
                    coor_qp_ma_re)
        gap = gapEval(coor_qp_sl_re, coor_qp_ma_re, norm_slav)
!print*, "COOR_SL: ", geom%coor_slav_curr(1,1:2)
!print*, "COOR_MA: ", geom%coor_mast_curr(1,1:2)
!print*, "NORM_SL: ", norm_slav
!print*, "NORM_MA: ", norm_mast
!print*, "COOR_QP: ", coor_qp_sl
!print*, "COOR_QP_RE: ", coor_qp_sl_re
!print*, "COOR_PJ: ", coor_qp_ma
!print*, "COOR_PJ_RE: ", coor_qp_ma_re
!print*, "GAP: ", gap
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine projQpSl2Vo(geom, coor_qp_sl, coor_qp_vo)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: coor_qp_sl(2)
        real(kind=8), intent(out) :: coor_qp_vo(2)
!
! --------------------------------------------------------------------------------------------------
!
!   Change parametric coordinates from face to volume
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: iret
        real(kind=8) :: coor_qp_sl_re(3)
!
        coor_qp_vo = 0.d0
!
! ----- Return in real face slave space (current configuration)
!
        coor_qp_sl_re = 0.d0
        call reerel(geom%elem_slav_code, geom%nb_node_slav, 3, geom%coor_slav_curr, coor_qp_sl, &
                    coor_qp_sl_re)
!
! ----- Projection of node on volumic slave cell (volumic parametric space)
!
        call reereg('S', geom%elem_volu_code, geom%nb_node_volu, geom%coor_volu_curr, &
                    coor_qp_sl_re, geom%elem_dime, coor_qp_vo, iret, ndim_coor_=3)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function gapEval(slav_pt, mast_pt, norm_slav)
!
        implicit none
!
        real(kind=8), intent(in) :: slav_pt(3), mast_pt(3), norm_slav(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Compute gap for raytracing gap = -(x^s - x^m).n^s
!
! --------------------------------------------------------------------------------------------------
!
        gapEval = -dot_product(slav_pt-mast_pt, norm_slav)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine thresEval(param, l_cont_qp, projRmVal, thres_qp, l_fric_qp)
!
        implicit none
!
        type(ContactParameters), intent(in) :: param
        aster_logical, intent(in) :: l_cont_qp
        real(kind=8), intent(in) :: projRmVal
        aster_logical, intent(out) :: l_fric_qp
        real(kind=8), intent(out) :: thres_qp
!
! --------------------------------------------------------------------------------------------------
!
!   Compute threshold for friction
!
! --------------------------------------------------------------------------------------------------
!
        thres_qp = 0.d0
        l_fric_qp = ASTER_FALSE
!
! ----- Define threshold for friction
        if (param%l_fric) then
            if (param%type_fric == FRIC_TYPE_TRES) then
                l_fric_qp = ASTER_TRUE
                thres_qp = param%threshold_given
            else if (param%type_fric == FRIC_TYPE_NONE) then
                l_fric_qp = ASTER_FALSE
                thres_qp = 0.d0
            else
                l_fric_qp = l_cont_qp
                if (l_cont_qp) then
                    if (param%type_fric == FRIC_TYPE_COUL) then
                        thres_qp = -param%threshold_given*projRmVal
                    else
                        thres_qp = THRES_STICK
                    end if
                    if (thres_qp == 0.d0) l_fric_qp = ASTER_FALSE
                else
                    thres_qp = 0.d0
                end if
            end if
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function speedEval(geom, coor_qp_slav, coor_qp_mast, gap)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: coor_qp_slav(2), coor_qp_mast(2), gap
        real(kind=8) :: speedEval(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Compute speed for raytracing v = -(x^s(t_prev) - x^m(t_prev) + gap * n^s(t_prev))
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: coor_qp_sl_prev(3), coor_qp_ma_prev(3), norm_slav_prev(3)
!
! ----- Return in real slave space
!
        coor_qp_sl_prev = 0.d0
        call reerel(geom%elem_slav_code, geom%nb_node_slav, 3, geom%coor_slav_prev, coor_qp_slav, &
                    coor_qp_sl_prev)
!
! ----- Return in real master space
!
        coor_qp_ma_prev = 0.d0
        call reerel(geom%elem_mast_code, geom%nb_node_mast, 3, geom%coor_mast_prev, coor_qp_mast, &
                    coor_qp_ma_prev)
!
! ------ Compute outward slave normal
!
        call apnorm(geom%nb_node_slav, geom%elem_slav_code, geom%elem_dime, geom%coor_slav_prev, &
                    coor_qp_slav(1), coor_qp_slav(2), norm_slav_prev)
!
        speedEval = -(coor_qp_sl_prev-coor_qp_ma_prev+gap*norm_slav_prev)
        ! write (6, *) 'coor_qp_slav=', coor_qp_slav
        ! write (6, *) 'coor_qp_sl_prev=', coor_qp_sl_prev
        ! write (6, *) 'coor_qp_ma_prev=', coor_qp_ma_prev
        ! write (6, *) 'gap=', gap
        ! write (6, *) 'norm_slav_prev=', norm_slav_prev
        ! write (6, *) 'geom%time_curr=', geom%time_curr
        ! write (6, *) 'geom%time_prev=', geom%time_prev
        ! write (6, *) 'speedEval=', speedEval
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine shapeFuncDisp(elem_dime, elem_nbnode, elem_code, coor_qp, shape_, &
                             dshape_, ddshape_)
!
        implicit none
!
        integer(kind=8), intent(in) :: elem_dime
        integer(kind=8), intent(in) :: elem_nbnode
        character(len=8), intent(in) :: elem_code
        real(kind=8), intent(in) :: coor_qp(2)
        real(kind=8), intent(out), optional :: shape_(9)
        real(kind=8), intent(out), optional :: dshape_(2, 9)
        real(kind=8), intent(out), optional :: ddshape_(3, 9)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate shape function and derivative for displacement
!
! --------------------------------------------------------------------------------------------------
!
        if (present(shape_)) then
            call mmnonf(elem_code, coor_qp(1), coor_qp(2), &
                        shape_)
        end if
        if (present(dshape_)) then
            call mmdonf(elem_code, coor_qp(1), coor_qp(2), &
                        dshape_)
        end if
        if (present(ddshape_)) then
            call mm2onf(elem_dime, elem_nbnode, elem_code, coor_qp(1), coor_qp(2), &
                        ddshape_)
        end if
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine shapeFuncDispVolu(elem_code, coor_qp, shape_, dshape_)
!
        implicit none
!
        character(len=8), intent(in) :: elem_code
        real(kind=8), intent(in) :: coor_qp(3)
        real(kind=8), intent(out), optional :: shape_(27)
        real(kind=8), intent(out), optional :: dshape_(3, 27)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate shape function and derivative for displacement (volumic cell)
!
! --------------------------------------------------------------------------------------------------
!
        if (present(shape_)) then
            shape_ = 0.d0
            call elrfvf(elem_code, coor_qp, shape_)
        end if
!
        if (present(dshape_)) then
            dshape_ = 0.d0
            call elrfdf(elem_code, coor_qp, dshape_)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine gradFuncDispVolu(geom, dshape, gradFunc)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: dshape(3, 27)
        real(kind=8), intent(out) :: gradFunc(3, 27)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate derivative for displacement (volumic cell)
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: g(3, 3), h(3, 3), jac
        integer(kind=8) :: ino, j, k
!
        g = 0.d0
!
        do ino = 1, geom%nb_node_volu
            do j = 1, geom%elem_dime
                g(1, j) = g(1, j)+geom%coor_volu_init(j, ino)*dshape(1, ino)
                g(2, j) = g(2, j)+geom%coor_volu_init(j, ino)*dshape(2, ino)
                g(3, j) = g(3, j)+geom%coor_volu_init(j, ino)*dshape(3, ino)
            end do
        end do
!
        if (geom%elem_dime == 2) g(3, 3) = 1.d0
!
        h(1, 1) = g(2, 2)*g(3, 3)-g(2, 3)*g(3, 2)
        h(2, 1) = g(3, 1)*g(2, 3)-g(2, 1)*g(3, 3)
        h(3, 1) = g(2, 1)*g(3, 2)-g(3, 1)*g(2, 2)
        h(1, 2) = g(1, 3)*g(3, 2)-g(1, 2)*g(3, 3)
        h(2, 2) = g(1, 1)*g(3, 3)-g(1, 3)*g(3, 1)
        h(3, 2) = g(1, 2)*g(3, 1)-g(3, 2)*g(1, 1)
        h(1, 3) = g(1, 2)*g(2, 3)-g(1, 3)*g(2, 2)
        h(2, 3) = g(2, 1)*g(1, 3)-g(2, 3)*g(1, 1)
        h(3, 3) = g(1, 1)*g(2, 2)-g(1, 2)*g(2, 1)
!
        jac = g(1, 1)*h(1, 1)+g(1, 2)*h(2, 1)+g(1, 3)*h(3, 1)
!
        gradFunc = 0.d0
!
        do ino = 1, geom%nb_node_volu
            do j = 1, geom%elem_dime
                do k = 1, geom%elem_dime
                    gradFunc(j, ino) = gradFunc(j, ino)+h(j, k)*dshape(k, ino)/jac
                end do
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine shapeFuncLagr(elem_code, coor_qp, shape_)
!
        implicit none
!
        character(len=8), intent(in) :: elem_code
        real(kind=8), intent(in) :: coor_qp(2)
        real(kind=8), intent(out) :: shape_(4)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate shape function and derivative (P1 for Lagrange)
!
! --------------------------------------------------------------------------------------------------
!
        character(len=8) :: elem_code_lagr
        integer(kind=8) :: elem_nbnode_lagr
        real(kind=8) :: ff(9)
!
        if (elem_code == "SE2") then
            elem_code_lagr = "SE2"
            elem_nbnode_lagr = 2
        else if (elem_code == "SE3") then
            elem_code_lagr = "SE3"
            elem_nbnode_lagr = 3
        else if (elem_code(1:2) == "TR") then
            elem_code_lagr = "TR3"
            elem_nbnode_lagr = 3
        else if (elem_code(1:2) == "QU") then
            elem_code_lagr = "QU4"
            elem_nbnode_lagr = 4
        else
            ASSERT(ASTER_FALSE)
        end if
!
        call mmnonf(elem_code_lagr, coor_qp(1), coor_qp(2), &
                    ff)
        shape_(1:4) = ff(1:4)
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function evalPoly(nb_node, shape, coeff_node)
!
        implicit none
!
        integer(kind=8), intent(in) :: nb_node
        real(kind=8), intent(in) :: coeff_node(*)
        real(kind=8), intent(in) :: shape(*)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate polynome via linear combinaison
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_node
!
        evalPoly = 0.d0
        do i_node = 1, nb_node
            evalPoly = evalPoly+coeff_node(i_node)*shape(i_node)
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function diameter(nb_node, nodes_coor)
!
        implicit none
!
        integer(kind=8), intent(in) :: nb_node
        real(kind=8), intent(in) :: nodes_coor(3, 9)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate diameter of cell
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_node_i, i_node_j
        real(kind=8) :: length
!
        diameter = 0.d0
        do i_node_i = 1, nb_node
            do i_node_j = i_node_i+1, nb_node
                length = norm2(nodes_coor(1:3, i_node_i)-nodes_coor(1:3, i_node_j))
                diameter = max(diameter, length)
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine testLagrC(geom, func_lagr, mu_c)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: func_lagr(4)
        real(kind=8), intent(out) :: mu_c(MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate test Lagrangian function
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_node, index
!
        mu_c = 0.d0
        index = 0
!
! --- Slave side
!
        do i_node = 1, geom%nb_lagr_c
            ASSERT(geom%indi_lagc(i_node) > 0)
            index = index+geom%elem_dime
            mu_c(index+1) = func_lagr(i_node)
            index = index+geom%indi_lagc(i_node)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine testLagrF(geom, func_lagr, mu_f)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: func_lagr(4)
        real(kind=8), intent(out) :: mu_f(MAX_LAGA_DOFS, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate test Lagrangian function
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_node, index
!
        mu_f = 0.d0
        index = 0
!
! --- Slave side
!
        do i_node = 1, geom%nb_lagr_c
            ASSERT(geom%indi_lagc(i_node) > 0)
            index = index+geom%elem_dime
            mu_f(index+2, 1) = func_lagr(i_node)
            if (geom%elem_dime == 3) then
                mu_f(index+3, 2) = func_lagr(i_node)
            end if
            index = index+geom%indi_lagc(i_node)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function evalDtestM(geom, dfunc_ma)
!
        implicit none
!
        type(ContactGeom), intent(in) :: geom
        real(kind=8), intent(in) :: dfunc_ma(2, 9)
        real(kind=8) :: evalDtestM(3, MAX_LAGA_DOFS, 2)
!
! --------------------------------------------------------------------------------------------------
!
!   Evaluate the first derivative of the master test function
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: j_node, index, k_zeta
!
        evalDtestM = 0.d0
        do k_zeta = 1, 2
            ! start at the first master dof (displ only)
            index = geom%nb_node_slav*geom%elem_dime+geom%nb_lagr_c*geom%indi_lagc(1)
            do j_node = 1, geom%nb_node_mast
                evalDtestM(1, index+1, k_zeta) = dfunc_ma(k_zeta, j_node)
                evalDtestM(2, index+2, k_zeta) = dfunc_ma(k_zeta, j_node)
                if (geom%elem_dime == 3) then
                    evalDtestM(3, index+3, k_zeta) = dfunc_ma(k_zeta, j_node)
                end if
                index = index+geom%elem_dime
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function barycenter(nb_nodes, coor_nodes)
!
        implicit none
!
        integer(kind=8), intent(in) :: nb_nodes
        real(kind=8), intent(in) :: coor_nodes(3, nb_nodes)
        real(kind=8) :: barycenter(3)
!
! --------------------------------------------------------------------------------------------------
!
!   Compute barycenter
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i_node
!
        barycenter = 0.d0
!
        do i_node = 1, nb_nodes
            barycenter = barycenter+coor_nodes(:, i_node)
        end do
!
        barycenter = barycenter/real(nb_nodes, kind=8)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
end module
