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
! aslint: disable=W0104,W1504
!
subroutine lc1015(fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  sigp, vip, typmod, &
                  dsidep, codret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lcgdpm.h"
#include "asterfort/nzgdzi.h"
#include "asterfort/postsm.h"
!
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
    real(kind=8), intent(in) :: sigm(nsig)
    real(kind=8), intent(in) :: vim(nvi)
    character(len=16), intent(in) :: option
    real(kind=8), intent(out) :: sigp(nsig)
    real(kind=8), intent(out) :: vip(nvi)
    character(len=8), intent(in) :: typmod(*)
    real(kind=8), intent(out) :: dsidep(*)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour - Special SIMO_MIEHE
!
! Metallurgy META_*
!
! --------------------------------------------------------------------------------------------------
!

    if (compor(META_PHAS) .eq. 'ACIERMECA' .or. &
        compor(META_PHAS) .eq. 'ACIER_REVENUMECA') then
        call lcgdpm(option, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, &
                    instam, instap, &
                    epsm, deps, &
                    nsig, sigm, &
                    nvi, vim, &
                    sigp, vip, dsidep, &
                    codret)
    else if (compor(META_PHAS) .eq. 'ZIRCMECA') then
        call nzgdzi(option, &
                    fami, kpg, ksp, ndim, imate, &
                    compor, carcri, &
                    instam, instap, &
                    epsm, deps, &
                    nsig, sigm, &
                    nvi, vim, &
                    sigp, vip, dsidep, &
                    codret)
    else
        ASSERT(ASTER_FALSE)
    end if
    call postsm(option, epsm, deps, sigm, sigp, dsidep)
!
end subroutine
