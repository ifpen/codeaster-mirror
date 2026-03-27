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
subroutine comp_meca_code(rela_comp, defo_comp, type_cpla, kit_comp, &
                          post_iter, regu_visc, post_incr, &
                          comp_code_py)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/lccree.h"
!
    character(len=16), intent(in) :: rela_comp, defo_comp, type_cpla, kit_comp(4)
    character(len=16), intent(in) :: post_iter, regu_visc, post_incr
    character(len=16), intent(out) :: comp_code_py
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Coding composite behaviour
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_comp        : RELATION comportment
! In  defo_comp        : DEFORMATION comportment
! In  type_cpla        : plane stress method
! In  kit_comp         : KIT comportment
! In  post_iter        : type of post_treatment at each Newton iteration
! In  regu_visc        : keyword for viscuous regularization
! In  post_incr        : type of post-treatment at end of time step
! Out comp_code_py     : composite coded comportment (coding in Python)
!
! --------------------------------------------------------------------------------------------------
!
! - Empty kit_comp for KIT_META
    character(len=16), parameter :: NoKitComp(4) = (/'VIDE', 'VIDE', 'VIDE', 'VIDE'/)
    integer(kind=8) :: nb_comp_elem, ikit
    character(len=16) :: comp_elem(20)
!
! --------------------------------------------------------------------------------------------------
!
    nb_comp_elem = 0
    comp_elem(1:20) = 'VIDE'

! - Create composite behaviour
    nb_comp_elem = nb_comp_elem+1
    comp_elem(nb_comp_elem) = rela_comp
    if (rela_comp .eq. 'KIT_META') then
        do ikit = 1, 4
            nb_comp_elem = nb_comp_elem+1
            comp_elem(nb_comp_elem) = NoKitComp(ikit)
        end do
    else
        do ikit = 1, 4
            nb_comp_elem = nb_comp_elem+1
            comp_elem(nb_comp_elem) = kit_comp(ikit)
        end do
    end if
    nb_comp_elem = nb_comp_elem+1
    comp_elem(nb_comp_elem) = regu_visc
    if (post_iter .ne. ' ') then
        nb_comp_elem = nb_comp_elem+1
        comp_elem(nb_comp_elem) = post_iter
    end if
    if (post_incr .ne. ' ') then
        nb_comp_elem = nb_comp_elem+1
        comp_elem(nb_comp_elem) = post_incr
    end if
    nb_comp_elem = nb_comp_elem+1
    comp_elem(nb_comp_elem) = defo_comp
    nb_comp_elem = nb_comp_elem+1
    comp_elem(nb_comp_elem) = type_cpla

! - Coding composite comportment (Python)
    call lccree(nb_comp_elem, comp_elem, comp_code_py)
!
end subroutine
