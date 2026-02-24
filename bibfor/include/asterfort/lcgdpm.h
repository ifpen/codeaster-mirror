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
    subroutine lcgdpm(option, &
                      fami, kpg, ksp, ndim, jvMaterCode, &
                      compor, carcri, &
                      timePrev, timeCurr, &
                      fm, df, &
                      nsig, sigm, &
                      nvi, vim, &
                      sigp, vip, dsigdf, &
                      codret)
        character(len=16), intent(in) :: option
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp, ndim, jvMaterCode
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8), intent(in) :: timePrev, timeCurr
        real(kind=8), intent(in) :: fm(3, 3), df(3, 3)
        integer(kind=8), intent(in) :: nsig
        real(kind=8), intent(in) :: sigm(nsig)
        integer(kind=8), intent(in) :: nvi
        real(kind=8), intent(in) :: vim(nvi)
        real(kind=8), intent(out) :: sigp(nsig)
        real(kind=8), intent(out) :: vip(nvi)
        real(kind=8), intent(out) :: dsigdf(6, 3, 3)
        integer(kind=8), intent(out) :: codret
    end subroutine lcgdpm
end interface
