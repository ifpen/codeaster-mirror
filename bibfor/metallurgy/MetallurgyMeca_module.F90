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
module MetallurgyMeca_module
! ==================================================================================================
    use Metallurgy_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: metaAnnealGetType
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Metallurgy_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! metaAnnealGetType
!
! Get type of hardening to apply annealing
!
! In  relaComp         : behaviour (RELATION keyword)
!
! --------------------------------------------------------------------------------------------------
    subroutine metaAnnealGetType(relaComp, lHardIsot, lHardKine, lHardMixed, &
                                 nbVariAnneal_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: relaComp
        aster_logical, intent(out) :: lHardIsot, lHardKine, lHardMixed
        integer(kind=8), optional, intent(out) :: nbVariAnneal_
! ----- Local
        integer(kind=8) :: nbVariAnneal
!   ------------------------------------------------------------------------------------------------
!
        lHardIsot = ASTER_FALSE
        lHardKine = ASTER_FALSE
        lHardMixed = ASTER_FALSE
        nbVariAnneal = 0
        if (relaComp .eq. 'VMIS_ISOT_LINE') then
            lHardIsot = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_ISOT_TRAC') then
            lHardIsot = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_CINE_LINE') then
            lHardKine = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_ECMI_LINE') then
            lHardMixed = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_CIN1_CHAB') then
            lHardIsot = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_CIN2_CHAB') then
            lHardIsot = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_ISOT_NL') then
            lHardIsot = ASTER_TRUE
        else
            ASSERT(ASTER_FALSE)
        end if
        if (lHardIsot) then
            nbVariAnneal = 1
        end if
        if (lHardKine) then
            nbVariAnneal = 6
        end if
        if (lHardMixed) then
            nbVariAnneal = 7
        end if
        if (present(nbVariAnneal_)) then
            nbVariAnneal_ = nbVariAnneal
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module MetallurgyMeca_module
