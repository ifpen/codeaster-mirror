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
#include "asterfort/Metallurgy_type.h"
!
interface
    subroutine metaGetPhase(fami, poum, kpg, ksp, &
                            metaType, nbPhase, &
                            phase_, zSumCold_, zHot_, tolePhase_)
        character(len=*), intent(in) :: fami
        character(len=1), intent(in) :: poum
        integer(kind=8), intent(in) :: kpg, ksp
        integer(kind=8), intent(in) :: metaType, nbPhase
        real(kind=8), optional, intent(out) :: phase_(nbPhase)
        real(kind=8), optional, intent(out) :: zSumCold_
        real(kind=8), optional, intent(out) :: zHot_
        real(kind=8), optional, intent(in) :: tolePhase_
    end subroutine metaGetPhase
end interface
