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
subroutine nmvcre(modelZ, materFieldZ, caraElemZ, varcRefe)
!
    implicit none
!
#include "asterfort/detrsd.h"
#include "asterfort/vrcref.h"
!
    character(len=*), intent(in) :: modelZ, materFieldZ, caraElemZ
    character(len=24), intent(in) :: varcRefe
!
! --------------------------------------------------------------------------------------------------
!
! External state variables
!
! Create external state variables for reference state
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  materField       : name of material characteristics (field)
! In  caraElem         : name of elementary characteristics (field)
! In  varcRefe         : name of external state variable field for reference values
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: varcRefeAll
!
! --------------------------------------------------------------------------------------------------
!
    varcRefeAll = varcRefe(1:14)//'.TOUT'
    call detrsd('VARI_COM', varcRefe)
    call vrcref(modelZ, materFieldZ, caraElemZ, varcRefeAll)
!
end subroutine
