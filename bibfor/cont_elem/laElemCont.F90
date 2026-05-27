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
subroutine laElemCont(parameters, geom, coor_qp_sl, hF, &
                      gap, gamma_c, projRmVal, l_cont_qp, &
                      gamma_f, l_fric_qp, norm_slav, &
                      dGap, projBsVal, dv, dvT, projBsVal2, lagr_c, lagr_f, lagr_f3, &
                      dThres_du, dThres_dl, lagr_v, thres, jump_v, dNs, &
                      tau_slav, lagr_g, mu_g, dts_ns, speed, dfunc_ma, dZetaM, &
                      mu_c, mu_f, mu_f3, metricTens, d2Gap)
!
    use contact_module
    use contact_type
    use contact_algebra_module
!
! aslint: disable=W1504
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "contact_module.h"
!
    type(ContactParameters), intent(in) :: parameters
    type(ContactGeom), intent(in) :: geom
    real(kind=8), intent(in) :: coor_qp_sl(2), hF
    real(kind=8), intent(out) :: gap, gamma_c, projRmVal
    real(kind=8), intent(out) :: gamma_f
    aster_logical, intent(out) :: l_cont_qp, l_fric_qp
    real(kind=8), intent(out), optional :: dts_ns(MAX_LAGA_DOFS, 2), dGap(MAX_LAGA_DOFS), thres
    real(kind=8), intent(out), optional :: speed(3), dvT(MAX_LAGA_DOFS, 3), mu_c(MAX_LAGA_DOFS)
    real(kind=8), intent(out), optional :: tau_slav(3, 2), dfunc_ma(3, MAX_LAGA_DOFS, 2)
    real(kind=8), intent(out), optional :: mu_f(MAX_LAGA_DOFS, 2), metricTens(2, 2)
    real(kind=8), intent(out), optional :: d2Gap(MAX_LAGA_DOFS, MAX_LAGA_DOFS), lagr_f3(3)
    real(kind=8), intent(out), optional :: lagr_g(3), mu_g(MAX_LAGA_DOFS, 3), lagr_c, lagr_f(2)
    real(kind=8), intent(out), optional :: jump_v(MAX_LAGA_DOFS, 3), dZetaM(MAX_LAGA_DOFS, 2)
    real(kind=8), intent(out), optional :: projBsVal(3), lagr_v(3), dNs(MAX_LAGA_DOFS, 3)
    real(kind=8), intent(out), optional :: mu_f3(MAX_LAGA_DOFS, 3)
    real(kind=8), intent(out), optional :: dv(MAX_LAGA_DOFS, 3), norm_slav(3), projBsVal2(2)
    real(kind=8), intent(out), optional :: dThres_du(MAX_LAGA_DOFS), dThres_dl(MAX_LAGA_DOFS)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (Lagrangian method) - Elementary computations
!
! Compute elementary contact quantities
!
! ------------------------------------l_fric,
! In  elem_dime        : dimension of elements
! In  l_axis           : .true. for axisymmetric element
! In  nb_lagr          : total number of Lagrangian dof on contact element
! In  indi_lagc        : node where Lagrangian dof is present (1) or not (0)
! In  lagr_g           : value of the Lagrange multiplier
! In  l_norm_smooth    : indicator for normals smoothing
! In  nb_node_slav     : number of nodes of for slave side from contact element
! In  elem_slav_code   : code element for slave side from contact element
! In  nb_node_mast     : number of nodes of for master side from contact element
! In  elem_mast_code   : code element for master side from contact element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: shape_func_sl(9), dshape_func_sl(2, 9), ddshape_func_sl(3, 9)
    real(kind=8) :: shape_func_ma(9), dshape_func_ma(2, 9), ddshape_func_ma(3, 9)
    real(kind=8) :: shape_func_lagr(4), lagr_v_(3), mu_f_(MAX_LAGA_DOFS, 2)
    real(kind=8) :: norm_mast(3), tau_slav_(3, 2), tau_mast(3, 2)
    real(kind=8) :: H, coor_qp_ma(2), lagrc_gap, lagr_c_, lagr_f_(2)
    real(kind=8) :: metricTens_(2, 2), invMetricTens(2, 2), lagr_g_(3)
    real(kind=8) :: metricTens_mast(2, 2), invMetricTens_mast(2, 2), H_mast(2, 2)
    real(kind=8) :: thres_qp, speed_(3), projBsVal_(3), mu_c_(MAX_LAGA_DOFS)
    real(kind=8) :: dNs_(MAX_LAGA_DOFS, 3), dGap_(MAX_LAGA_DOFS)
    real(kind=8) :: jump_v_(MAX_LAGA_DOFS, 3), dZetaM_(MAX_LAGA_DOFS, 2)
    real(kind=8) :: dts_ns_(MAX_LAGA_DOFS, 2), norm_slav_(3)
    real(kind=8) :: dTm_nm(MAX_LAGA_DOFS, 2)
    real(kind=8) :: dTs(MAX_LAGA_DOFS, 3, 2)
    real(kind=8) :: d2Ns_nm(MAX_LAGA_DOFS, MAX_LAGA_DOFS), dTs_nm(MAX_LAGA_DOFS, 2)
    aster_logical :: glob_coor
    blas_int :: b_3, b_2, b_MAX_LAGA_DOFS
!
! --------------------------------------------------------------------------------------------------
!
    b_MAX_LAGA_DOFS = to_blas_int(MAX_LAGA_DOFS)
    b_2 = 2
    b_3 = 3
    glob_coor = (parameters%l_fric) .and. (parameters%vari_cont == CONT_VARI_ROBU)
    ! write (6, *) 'glob_coor=', glob_coor
    ! write (6, *) 'parameters%l_fric=', parameters%l_fric
    ! write (6, *) 'parameters%vari_cont=', parameters%vari_cont
    ! write (6, *) 'CONT_VARI_ROBU=', CONT_VARI_ROBU
!
! ----- Project quadrature point (on master side)
!
    call projQpSl2Ma(geom, coor_qp_sl, parameters%proj_tole, &
                     coor_qp_ma, gap, tau_slav_, norm_slav_, tau_mast, norm_mast)
    if (present(norm_slav)) norm_slav = norm_slav_
    if (present(tau_slav)) tau_slav = tau_slav_
    ! if (present(k_diff)) write (6, *) '*tau_slav*', k_diff(1:3), '*', tau_slav_(:, 1)
!
! ----- Evaluate shape function for displacement (slave and master)
!
    call shapeFuncDisp(geom%elem_dime, geom%nb_node_slav, geom%elem_slav_code, coor_qp_sl, &
                       shape_func_sl, dshape_func_sl, ddshape_func_sl)
    call shapeFuncDisp(geom%elem_dime, geom%nb_node_mast, geom%elem_mast_code, coor_qp_ma, &
                       shape_func_ma, dshape_func_ma, ddshape_func_ma)
!
! ----- Store first derivative of the shape functions (master)
!
    if (present(dfunc_ma)) dfunc_ma = evalDtestM(geom, dshape_func_ma)
!
! ----- Evaluate shape function for Lagrange (slave)
!
    call shapeFuncLagr(geom%elem_slav_code, coor_qp_sl, &
                       shape_func_lagr)
!
! ----- Evaluate Lagr_c and gamma_c at quadrature point
!
    lagr_c_ = evalPoly(geom%nb_lagr_c, shape_func_lagr, geom%lagc_slav_curr)
    gamma_c = evalPoly(geom%nb_lagr_c, shape_func_lagr, parameters%coef_cont)/hF
    if (geom%nb_lagr_c == 0) then
        ! gamma_c = parameters%coef_cont(1)/hF
        gamma_c = evalPoly(geom%nb_node_slav, shape_func_sl, parameters%coef_cont)/hF
    end if

    lagr_f_(1) = evalPoly(geom%nb_lagr_c, shape_func_lagr, geom%lagf_slav_curr(1, :))
    lagr_f_(2) = evalPoly(geom%nb_lagr_c, shape_func_lagr, geom%lagf_slav_curr(2, :))
!
    if (present(lagr_c)) lagr_c = lagr_c_
    if (present(lagr_f)) lagr_f = lagr_f_
!
    if (present(lagr_f3)) lagr_f3 = matmul(tau_slav_, lagr_f)
!
! ----- Lagrange multiplier in global coordinate at quadrature point
!
    lagr_g_ = (/lagr_c_, lagr_f_(1), lagr_f_(2)/)
    if (present(lagr_g)) then
        lagr_g = lagr_g_
    end if
!
    if (glob_coor) then
        lagrc_gap = dot_product(lagr_g_, norm_slav_)+gamma_c*gap
    else
        lagrc_gap = lagr_c_+gamma_c*gap
    end if
!
! ----- Contact activate at quadrature point ( H = 0 or 1 )
!
    if (parameters%type_cont == CONT_TYPE_UNIL) then
        if (parameters%cont_init == PAIR_CONT_INTE) then
            H = Heaviside(-lagrc_gap)
            projRmVal = projRm(lagrc_gap)
        elseif (parameters%cont_init == FRIC_ALGO_FALSE) then
            H = 0.d0
            projRmVal = 0.d0
        elseif (parameters%cont_init == FRIC_ALGO_TRUE) then
            H = 1.d0
            projRmVal = lagrc_gap
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        H = 1.d0
        projRmVal = lagrc_gap
    end if
    l_cont_qp = (H > 0.5d0)
!
! ----- Evaluate contact derivative
!
    metricTens_mast = metricTensor(tau_mast)
    invMetricTens_mast = invMetricTensor(geom, metricTens_mast)
    H_mast = secondFundForm(geom%elem_dime, geom%nb_node_mast, geom%coor_mast_pair, &
                            ddshape_func_ma, norm_mast)
!
    jump_v_ = jump(geom, shape_func_sl, shape_func_ma)
    if (present(jump_v)) jump_v = jump_v_
!
    dts_ns_ = dTs_du_ns(geom, dshape_func_sl, norm_slav_)
    if (present(dts_ns)) dts_ns = dts_ns_
!
    dNs_ = dNs_du(geom, tau_slav_, dts_ns_)
    if (present(dNs)) dNs = dNs_
    call dPi_du(geom, norm_slav_, tau_mast, gap, jump_v_, dNs_, dGap_, dZetaM_)
!
    if (present(dZetaM)) dZetaM = dZetaM_
!
! ----- Evaluate Lagr_f and gamma_f at quadrature point
!
    gamma_f = 0.d0
    l_fric_qp = ASTER_FALSE
    ! if (present(projBsVal)) projBsVal = 0.d0
    !
    if (parameters%l_fric) then
!
        call thresEval(parameters, l_cont_qp, projRmVal, thres_qp, l_fric_qp)

        if (present(thres)) thres = thres_qp
!
        metricTens_ = metricTensor(tau_slav_)
        if (present(metricTens)) metricTens = metricTens_
        invMetricTens = invMetricTensor(geom, metricTens_)
!
        if (glob_coor) then
            gamma_f = gamma_c
        else
            gamma_f = evalPoly(geom%nb_lagr_c, shape_func_lagr, parameters%coef_fric)/hF
        end if
!
        speed_ = speedEval(geom, coor_qp_sl, coor_qp_ma, gap)
        if (present(speed)) speed = speed_
!
        if (glob_coor) then
            lagr_v_ = lagr_g_-gamma_f*speed_
        else
            lagr_v_ = (lagr_c_*norm_slav_+matmul(tau_slav_, lagr_f_))-gamma_f*speed_
        end if

        if (present(lagr_v)) lagr_v = lagr_v_
!
        projBsVal_ = projBs(parameters, lagr_v_, thres_qp, norm_slav_)
        if (present(projBsVal)) projBsVal = projBsVal_
!
        if (present(dv)) then
            dv = dv_du(geom, parameters%proj_tole, dGap_, dZetaM_, coor_qp_sl, dshape_func_ma)
        end if
!
        if (present(projBsVal2)) then
            projBsVal2 = cmpTang(invMetricTens, tau_slav_, projBsVal_)
        end if
!
        if (present(dvT)) then
            dvT = dvT_du(geom, speed_, norm_slav_, parameters%proj_tole, &
                         dGap_, invMetricTens, tau_slav_, dZetaM_, dTs_ns_, &
                         coor_qp_sl, dshape_func_ma)
        end if
!
    end if
!
    if (present(dGap)) then
!
! ----- Compute d (gap(u))[v] / du
!
        dGap = dGap_
!
    end if
!
    if (present(d2Gap)) then
!
! ----- Compute d^2 (gap(u))[v, w] / du^2
!
        dTs = dTs_du(geom, dshape_func_sl)
        dTs_nm = dTs_du_nm(geom, dshape_func_sl, norm_mast)
        dTm_nm = dTm_du_nm(geom, dshape_func_ma, norm_mast)
        d2Ns_nm = d2Ns_du2_nm(geom, tau_slav_, norm_mast, dNs_, dTs, dTs_ns_, dTs_nm)
        d2Gap = d2Gap_du2(geom, norm_slav_, norm_mast, gap, H_mast, &
                          dNs_, dGap_, dZetaM_, dTm_nm, d2Ns_nm)

    end if
!
!
! ----- Compute Lagrange multipliers in local coordinates
!
    call testLagrC(geom, shape_func_lagr, mu_c_)
!
    if (parameters%l_fric) then
        call testLagrF(geom, shape_func_lagr, mu_f_)
    else
        mu_f_ = 0.d0
    end if
!
    if (present(mu_c)) mu_c = mu_c_
!
    if (present(mu_f)) mu_f = mu_f_
!
    if (present(mu_g)) then
!
! ----- Compute mu_g in global coordinates
!
        mu_g(:, 1) = mu_c_(:)
        mu_g(:, 2) = mu_f_(:, 1)
        mu_g(:, 3) = mu_f_(:, 2)
!
    end if
!
    if (present(mu_f3)) then
        !
        ! ----- Compute mu_f in global coordinates
        !
        mu_f3 = matmul(mu_f, transpose(tau_slav_))
    end if

!
    if (parameters%l_fric) then
!
! ----- Compute threshold derivatives
!
        if (present(dThres_du)) then
            dThres_du = -H*parameters%threshold_given*(gamma_c*dGap_+ &
                                                       matmul(dNs_, lagr_g_))
        end if

        if (present(dThres_dl)) then
            dThres_dl = -H*parameters%threshold_given*matmul(mu_g, norm_slav)
        end if
    end if

!
end subroutine
