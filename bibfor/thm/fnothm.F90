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
subroutine fnothm(ds_thm, jv_mater, ndim, l_axi, fnoevo, &
                  mecani, press1, press2, tempe, second, &
                  nno, nnos, npi, npg, &
                  elem_coor, deltat, dimdef, dimcon, dimuel, inte_type, &
                  jv_poids, jv_poids2, &
                  jv_func, jv_func2, jv_dfunc, jv_dfunc2, &
                  nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                  congem, congep, b, r, vectu)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cabthm.h"
#include "asterfort/fonoda.h"
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
    real(kind=8), intent(in) :: deltat
    integer(kind=8), intent(in) :: dimuel, dimdef, dimcon
    character(len=3), intent(in) :: inte_type
    integer(kind=8), intent(in) :: jv_poids, jv_poids2
    integer(kind=8), intent(in) :: jv_func, jv_func2, jv_dfunc, jv_dfunc2
    integer(kind=8), intent(in) :: nddls, nddlm
    integer(kind=8), intent(in) :: nddl_meca, nddl_p1, nddl_p2, nddl_2nd
    real(kind=8), intent(inout) :: congem(1:npi*dimcon)
    real(kind=8), intent(inout) :: congep(1:npi*dimcon)
    real(kind=8), intent(inout) :: b(dimdef, dimuel)
    real(kind=8), intent(inout) :: r(1:dimdef+1)
    real(kind=8), intent(out) :: vectu(dimuel)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute nodal force vector (FORC_NODA)
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! In  jv_mater         : coded material address
! In  ndim             : dimension of element (2 ou 3)
! In  l_axi            : flag is axisymmetric model
! In  l_steady         : .true. for steady state
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
! In  deltat           : time increment
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
! In  nddl_2nd         : number of dof for second gradient quantity
! IO  congem           : generalized stresses at the beginning of time step
!                    => output sqrt(2) on SIG_XY, SIG_XZ, SIG_YZ
! IO  congep           : generalized stresses at the end of time step
!                    => output sqrt(2) on SIG_XY, SIG_XZ, SIG_YZ
! IO  b                : [B] matrix for generalized strains
! IO  r                : stress vector
! Out vectu            : nodal force vector (FORC_NODA)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: kpi, i, n
    real(kind=8) :: dfdi(20, 3), dfdi2(20, 3), poids, poids2
    integer(kind=8) :: addeme, addete, addep1, addep2, adde2nd
    real(kind=8) :: r_bar(1:dimdef+1)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nno .le. 20)
    ASSERT(ndim .le. 3)
    ASSERT(npi .le. 27)
    vectu(1:dimuel) = 0.d0
!
! - Get adresses in generalized vectors
!
    addeme = mecani(2)
    addep1 = press1(3)
    addep2 = press2(3)
    addete = tempe(2)
    adde2nd = second(2)
!
! - Loop on Gauss points
!
    do kpi = 1, npi
        r(1:dimdef+1) = 0.d0
! ----- Compute [B] matrix for generalized strains
        call cabthm(ds_thm, l_axi, ndim, &
                    nddls, nddlm, &
                    nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                    nno, nnos, &
                    dimuel, dimdef, kpi, &
                    addeme, addete, addep1, addep2, adde2nd, &
                    elem_coor, &
                    jv_poids, jv_poids2, &
                    jv_func, jv_func2, &
                    jv_dfunc, jv_dfunc2, &
                    dfdi, dfdi2, &
                    poids, poids2, &
                    b)
! ----- Compute stress vector {R}
        call fonoda(ds_thm, jv_mater, ndim, fnoevo, &
                    mecani, press1, press2, tempe, second, &
                    dimdef, dimcon, deltat, congem((kpi-1)*dimcon+1), &
                    congep((kpi-1)*dimcon+1), r)

! ----- Compute stress vector {R_bar} for selective integration
        if (inte_type .eq. 'RED') then
            if (kpi .le. npg) then
                r_bar = r
                if (ds_thm%ds_elem%l_dof_pre1) r_bar(addep1) = 0.d0
                if (ds_thm%ds_elem%l_dof_pre2) r_bar(addep2) = 0.d0
                if (ds_thm%ds_elem%l_dof_ther) r_bar(addete) = 0.d0
            else
                r_bar = 0.d0
                if (ds_thm%ds_elem%l_dof_pre1) r_bar(addep1) = r(addep1)
                if (ds_thm%ds_elem%l_dof_pre2) r_bar(addep2) = r(addep2)
                if (ds_thm%ds_elem%l_dof_ther) r_bar(addete) = r(addete)
            end if
        else
            r_bar = r
        end if

! ----- Compute residual = [B]^T.{R_bar}
        do i = 1, dimuel
            do n = 1, dimdef
                vectu(i) = vectu(i)+b(n, i)*r_bar(n)*poids
            end do
        end do
    end do
!
end subroutine
