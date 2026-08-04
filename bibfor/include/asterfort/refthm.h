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
#include "asterf_types.h"
!
interface
    subroutine refthm(ds_thm   ,&
                      jv_mater , ndim     , l_axi    , fnoevo ,&
                      mecani   , press1   , press2   , tempe    , second, &
                      nno      , nnos     , npi      , npg      ,&
                      elem_coor, dt       , dimdef   , dimcon   , dimuel, inte_type,&
                      jv_poids , jv_poids2,&
                      jv_func  , jv_func2 , jv_dfunc , jv_dfunc2,&
                      nddls    , nddlm    , nddl_meca, nddl_p1  , nddl_p2, nddl_2nd, &
                      b        , r        , vectu )
        use THM_type
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
    end subroutine refthm
end interface
