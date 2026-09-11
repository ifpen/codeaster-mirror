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
! aslint: disable=W1504
!
subroutine refthm(ds_thm, jv_mater, ndim, l_axi, fnoevo, &
                  mecani, press1, press2, tempe, second, &
                  nno, nnos, npi, npg, elem_coor, &
                  dt, dimdef, dimcon, dimuel, inte_type, jv_poids, &
                  jv_poids2, jv_func, jv_func2, jv_dfunc, jv_dfunc2, &
                  nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, &
                  nddl_2nd, b, r, vectu)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8miem.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/fnothm.h"
#include "asterfort/r8inir.h"
#include "asterfort/terefe.h"
#include "blas/daxpy.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    integer(kind=8), intent(in) :: jv_mater
    integer(kind=8), intent(in) :: ndim
    aster_logical, intent(in) :: l_axi
    aster_logical, intent(in) :: fnoevo
    integer(kind=8), intent(in) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
    integer(kind=8), intent(in) :: nno, nnos
    integer(kind=8), intent(in) :: npi, npg
    real(kind=8) :: elem_coor(ndim, nno)
    real(kind=8), intent(in) :: dt
    integer(kind=8), intent(in) :: dimuel, dimdef, dimcon
    character(len=3), intent(in) :: inte_type
    integer(kind=8), intent(in) :: jv_poids, jv_poids2
    integer(kind=8), intent(in) :: jv_func, jv_func2, jv_dfunc, jv_dfunc2
    integer(kind=8), intent(in) :: nddls, nddlm
    integer(kind=8), intent(in) :: nddl_meca, nddl_p1, nddl_p2, nddl_2nd
    real(kind=8), intent(out) :: b(dimdef, dimuel)
    real(kind=8), intent(out) :: r(1:dimdef+1)
    real(kind=8), intent(out) :: vectu(dimuel)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute REFE_FORC_NODA at Gauss points on current element
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! In  jv_mater         : coded material address
! In  ndim             : dimension of element (2 ou 3)
! In  l_axi            : flag is axisymmetric model
! In  fnoevo           : .true. if compute in non-linear operator (transient terms)
! In  mecani           : parameters for mechanic
! In  press1           : parameters for hydraulic (first pressure)
! In  press1           : parameters for hydraulic (second pressure)
! In  tempe            : parameters for thermic
! In  second           : parameters for second gradient
! In  nno              : number of nodes (all)
! In  nnos             : number of nodes (not middle ones)
! In  npi              : number of Gauss points for linear
! In  npg              : number of Gauss points
! In  elem_coor        : coordinates of nodes for current element
! In  dt               : time increment
! In  dimdef           : number of generalized strains
! In  dimcon           : dimension of generalized stresses vector
! In  dimuel           : number of dof for element
! In  inte_type        : type of integration - classical, lumped (D), reduced (R)
! In  jv_poids         : JEVEUX adress for weight of Gauss points (linear shape functions)
! In  jv_poids2        : JEVEUX adress for weight of Gauss points (quadratic shape functions)
! In  jv_func          : JEVEUX adress for shape functions (linear shape functions)
! In  jv_func2         : JEVEUX adress for shape functions (quadratic shape functions)
! In  jv_dfunc         : JEVEUX adress for derivative of shape functions (linear shape functions)
! In  jv_dfunc2        : JEVEUX adress for derivative of shape functions (quadratic shape functions)
! In  nddls            : number of dof at nodes (not middle ones)
! In  nddlm            : number of dof at nodes (middle ones)
! In  nddl_meca        : number of dof for mechanical quantity
! In  nddl_p1          : number of dof for first hydraulic quantity
! In  nddl_p2          : number of dof for second hydraulic quantity
! In  nddl_second      : number of dof for second gradient quantity
! Out b                : [B] matrix for generalized strains
! Out r                : stress vector
! Out vectu            : nodal force vector (FORC_NODA)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: indx_vale_refe, kpi, i_dim, k
    integer(kind=8), parameter :: parsig = 27*36
    integer(kind=8), parameter :: partmp = 27*6
    integer(kind=8), parameter :: parbsi = 27*6
    real(kind=8) :: sigtm(parsig), ftemp(partmp), bsigm(parbsi)
    real(kind=8) :: vale_refe, list_vale_refe(6)
    character(len=16) :: type_elem
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call r8inir(dimcon*npi, 0.d0, sigtm(1), 1)
    call r8inir(dimuel, 0.d0, ftemp(1), 1)
    ASSERT(nddls .le. 6)
    ASSERT(nno .le. 27)
    ASSERT(npi .le. 27)
    ASSERT(dimcon .le. 31+5)

    if (ds_thm%ds_elem%l_dof_2nd) then
        type_elem = "THM_DIL"
    else
        type_elem = "THM"
    end if
!
! - Check which *_REFE exist
!
    if (ds_thm%ds_elem%l_dof_meca) then
        call terefe('SIGM_REFE', type_elem, vale_refe)
        indx_vale_refe = 1
        list_vale_refe(indx_vale_refe) = vale_refe
    end if
    if (ds_thm%ds_elem%l_dof_pre1) then
        call terefe('FLUX_HYD1_REFE', type_elem, vale_refe)
        indx_vale_refe = 2
        list_vale_refe(indx_vale_refe) = vale_refe
    end if
    if (ds_thm%ds_elem%l_dof_pre2) then
        call terefe('FLUX_HYD2_REFE', type_elem, vale_refe)
        indx_vale_refe = 3
        list_vale_refe(indx_vale_refe) = vale_refe
    end if
    if (ds_thm%ds_elem%l_dof_ther) then
        call terefe('FLUX_THER_REFE', type_elem, vale_refe)
        indx_vale_refe = 4
        list_vale_refe(indx_vale_refe) = vale_refe
    end if
    if (ds_thm%ds_elem%l_dof_2nd) then
        call terefe('LAGR_REFE', type_elem, vale_refe)
        indx_vale_refe = 5
        list_vale_refe(indx_vale_refe) = vale_refe
        call terefe('EPSI_REFE', type_elem, vale_refe)
        indx_vale_refe = 6
        list_vale_refe(indx_vale_refe) = vale_refe
    end if
!
! - Compute
!
    do kpi = 1, npi
        do i_dim = 1, dimcon
! --------- Get current dof
            if (i_dim .le. mecani(5)) then
                indx_vale_refe = 1
            else if (i_dim .le. (mecani(5)+press1(2)*press1(7))) then
                indx_vale_refe = 2
                if (tempe(5) .gt. 0) then
                    if (i_dim .eq. (mecani(5)+press1(7)) .or. i_dim .eq. &
                        (mecani(5)+press1(2)*press1(7))) then
                        cycle
                    end if
                end if
            else if (i_dim .le. (mecani(5)+press1(2)*press1(7)+press2(2)*press2(7))) then
                indx_vale_refe = 3
                if (tempe(5) .gt. 0) then
                    if (i_dim .eq. (mecani(5)+press1(2)*press1(7)+press2(7)) .or. i_dim &
                        .eq. (mecani(5)+press1(2)*press1(7)+press2(2)*press2(7))) then
                        cycle
                    end if
                end if
            else if (i_dim .le. (mecani(5)+tempe(5))) then
                indx_vale_refe = 4
            else if (i_dim .ge. second(3)) then
                if (i_dim .eq. second(3)) then
                    indx_vale_refe = 5
                else if (i_dim .eq. (second(3)+1)) then
                    indx_vale_refe = 5
                else if ((i_dim .ge. (second(3)+2)) .and. (i_dim .le. (second(3)+1+ndim))) then
                    indx_vale_refe = 1
                else if (i_dim .eq. (second(3)+second(5)-1)) then
                    indx_vale_refe = 6
                end if
            end if
! --------- Get *_REFE value
            vale_refe = list_vale_refe(indx_vale_refe)
! --------- Compute
            if (vale_refe .ne. r8vide()) then
                sigtm(i_dim+dimcon*(kpi-1)) = vale_refe
                call fnothm(ds_thm, jv_mater, ndim, l_axi, fnoevo, &
                            mecani, press1, press2, tempe, second, &
                            nno, nnos, npi, npg, elem_coor, &
                            dt, dimdef, dimcon, dimuel, inte_type, jv_poids, &
                            jv_poids2, jv_func, jv_func2, jv_dfunc, jv_dfunc2, &
                            nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, &
                            nddl_2nd, sigtm, sigtm, b, r, bsigm(1))
                do k = 1, dimuel
                    ftemp(k) = ftemp(k)+abs(bsigm(k))
                end do
                sigtm(i_dim+dimcon*(kpi-1)) = 0.0d0
            end if
        end do
    end do
!
    b_n = to_blas_int(dimuel)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0/npi, ftemp(1), b_incx, vectu(1), &
               b_incy)
!
! - Check
!
    do k = 1, dimuel
        ASSERT(abs(vectu(k)) .gt. r8miem())
    end do
!
end subroutine
