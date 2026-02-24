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
! aslint: disable=W1504,W0104,C1505
!
subroutine lc0015(BEHInteg, &
                  option, typmod, &
                  fami, kpg, ksp, ndim, jvMaterCode, &
                  compor, carcri, timePrev, timeCurr, &
                  neps, epsm, deps, &
                  nsig, sigm, &
                  nvi, vim, &
                  sigp, vip, &
                  ndsde, dsidep, codret)
!
    use Behaviour_type
    use MetallurgyMeca_module, only: metaGetMechanism
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/nzcifw.h"
#include "asterfort/nzcizi.h"
#include "asterfort/nzedga.h"
#include "asterfort/nzisfw.h"
!
    type(Behaviour_Integ), intent(in):: BEHInteg
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: typmod(2)
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, ndim, jvMaterCode
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: timePrev, timeCurr
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps), deps(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm(nsig)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    real(kind=8), intent(out) :: sigp(nsig)
    real(kind=8), intent(out)  :: vip(nvi)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(out):: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                                       merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! KIT_META (steel and zircaloy)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: metaRela, metaGlob, metaPhas
    aster_logical :: l_hard_kine
!
! --------------------------------------------------------------------------------------------------
!
    metaPhas = compor(META_PHAS)
    metaRela = compor(META_RELA)
    metaGlob = compor(META_GLOB)
    call metaGetMechanism(metaRela, metaGlob, lHardKine_=l_hard_kine)

    if (l_hard_kine) then
        if (metaPhas .eq. 'ACIERMECA' .or. metaPhas .eq. 'ACIER_REVENUMECA') then
            call nzcifw(option, &
                        fami, kpg, ksp, ndim, jvMaterCode, &
                        compor, carcri, &
                        timePrev, timeCurr, &
                        neps, epsm, deps, &
                        nsig, sigm, &
                        nvi, vim, &
                        sigp, vip, ndsde, dsidep, &
                        codret)
        elseif (metaPhas .eq. 'ZIRCMECA') then
            call nzcizi(option, &
                        fami, kpg, ksp, ndim, jvMaterCode, &
                        compor, carcri, &
                        timePrev, timeCurr, &
                        neps, epsm, deps, &
                        nsig, sigm, &
                        nvi, vim, &
                        sigp, vip, ndsde, dsidep, &
                        codret)
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        if (metaPhas .eq. 'ACIERMECA' .or. metaPhas .eq. 'ACIER_REVENUMECA') then
            call nzisfw(option, &
                        fami, kpg, ksp, ndim, jvMaterCode, &
                        compor, carcri, &
                        timePrev, timeCurr, &
                        neps, epsm, deps, &
                        nsig, sigm, &
                        nvi, vim, &
                        sigp, vip, ndsde, dsidep, &
                        codret)
        elseif (metaPhas .eq. 'ZIRCMECA') then
            call nzedga(option, &
                        fami, kpg, ksp, ndim, jvMaterCode, &
                        compor, carcri, &
                        timePrev, timeCurr, &
                        neps, epsm, deps, &
                        nsig, sigm, &
                        nvi, vim, &
                        sigp, vip, ndsde, dsidep, &
                        codret)
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
end subroutine
