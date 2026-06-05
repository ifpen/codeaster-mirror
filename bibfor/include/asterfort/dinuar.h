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
    subroutine dinuar(result, sddisc, numeInst, lForceStore, &
                      numeStoring, numeReuseCalc_, lStoreInitState_)
        use NonLin_Datastructure_type
        character(len=8), intent(in) :: result
        character(len=19), intent(in) :: sddisc
        integer(kind=8), intent(in) :: numeInst
        aster_logical, intent(in) :: lForceStore
        integer(kind=8), intent(out) :: numeStoring
        integer(kind=8), optional, intent(out) :: numeReuseCalc_
        aster_logical, intent(in), optional :: lStoreInitState_
    end subroutine dinuar
end interface
