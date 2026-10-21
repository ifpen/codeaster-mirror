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

subroutine afvarc(mateField, mesh, model)
!
    use ExternalStateVariable_type
    use ExternalStateVariable_module, only: creaCata, freeCata
    use ExternalStateVariableRead_module, only: readDataFromUser, freeDataFromUser
    use ExternalStateVariable_module, only: shrinkMaps
    use ExternalStateVariable_module, only: creaJvObjects, creaMaps, fillJvObjects, fillMaps
    use ExternalStateVariable_module, only: debugMaps
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/ExternalStateVariable_type.h"
!
    character(len=8), intent(in) :: mateField, mesh, model
!
! --------------------------------------------------------------------------------------------------
!
! Material - External state variables (VARC)
!
! For AFFE_MATERIAU/AFFE_VARC
!
! --------------------------------------------------------------------------------------------------
!
! In  mateField        : name of material field (CHAM_MATER)
! In  mesh             : name of mesh
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: jvBase = "G"
    type(EXTE_VARI_CATA), pointer :: exteVariCata(:) => null()
    type(EXTE_VARI_AFFE) :: exteVariAffe
!
! --------------------------------------------------------------------------------------------------

! - Create catalog of all external state variables
    call creaCata(exteVariCata)

! - Get external state variables affected from command file
    call readDataFromUser(exteVariCata, exteVariAffe)

! - Create
    if (exteVariAffe%nbAffe .ne. 0) then
! ----- Create main objects
        call creaJvObjects(jvBase, mateField, exteVariAffe)

! ----- Create maps
        call creaMaps(jvBase, mesh, mateField, exteVariCata, exteVariAffe)

! ----- Affect values in main objects
        call fillJvObjects(mateField, exteVariCata, exteVariAffe)

! ----- Affect values in maps
        call fillMaps(mateField, mesh, model, &
                      exteVariCata, exteVariAffe)

! ----- Shrink number of components to save memory
        call shrinkMaps(mateField, exteVariAffe)

! ----- Debug if required
        if (EXTEVARI_DBG_READ .eq. 1) then
            call debugMaps(mateField)
        end if
    end if

! - Free objects
    call freeCata(exteVariCata)
    call freeDataFromUser(exteVariAffe)
!
end subroutine
