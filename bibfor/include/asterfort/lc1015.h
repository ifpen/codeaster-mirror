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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine lc1015(fami, kpg, ksp, ndim, imate, &
                      compor, carcri, instam, instap, epsm, &
                      deps, nsig, sigm, nvi, vim, option, &
                      sigp, vip, typmod, &
                      dsidep, codret)
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg
        integer(kind=8), intent(in) :: ksp
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: imate, nvi, nsig
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8), intent(in) :: instam, instap
        real(kind=8), intent(in) :: epsm(*)
        real(kind=8), intent(in) :: deps(*)
        real(kind=8), intent(in) :: sigm(*)
        real(kind=8), intent(in) :: vim(nvi)
        character(len=16), intent(in) :: option
        real(kind=8), intent(out) :: sigp(*)
        real(kind=8), intent(out) :: vip(nvi)
        character(len=8), intent(in) :: typmod(*)
        real(kind=8), intent(out) :: dsidep(*)
        integer(kind=8), intent(out) :: codret
    end subroutine lc1015
end interface
