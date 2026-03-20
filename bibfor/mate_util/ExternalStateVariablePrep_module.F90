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
! ==================================================================================================
!
! Module to prepare External State Variables
!
! ==================================================================================================
!
module ExternalStateVariablePrep_module
! ==================================================================================================
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: checkField
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! checkField
!
! Check conformity of field from user
!
! --------------------------------------------------------------------------------------------------
    subroutine checkField(meshZ, modelZ, physQuanZ, &
                          exteVariNameZ, fieldUserZ, &
                          fieldDisc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ, modelZ, physQuanZ, exteVariNameZ, fieldUserZ
        character(len=8), intent(out) :: fieldDisc
! ----- Locals
        aster_logical :: lFieldELGA, lXFEM
        character(len=24) :: fieldUser, fieldFED, modelFED, valk(3)
        character(len=8) :: mesh, fieldMesh, exteVariName
        character(len=8) :: physQuan, fieldPhysQuan
        character(len=16) :: fieldOption
        integer(kind=8) :: iret, nbVari
        integer(kind=8), pointer :: celd(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        physQuan = physQuanZ
        exteVariName = exteVariNameZ
        fieldUser = fieldUserZ
        fieldDisc = " "

! ----- Get parameters
        call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelFED)
        call exixfe(modelZ, iret)
        lXFEM = iret .ne. 0
        call dismoi('NOM_MAILLA', fieldUser, 'CHAMP', repk=fieldMesh)
        call dismoi('NOM_GD', fieldUser, 'CHAMP', repk=fieldPhysQuan)
        call dismoi('TYPE_CHAMP', fieldUser, 'CHAMP', repk=fieldDisc)
        lFieldELGA = ASTER_FALSE
        if (fieldDisc .eq. "ELGA" .and. .not. lXFEM) then
            lFieldELGA = ASTER_TRUE
            call dismoi('NOM_LIGREL', fieldUser, 'CHAM_ELEM', repk=fieldFED)
            call dismoi('NOM_OPTION', fieldUser, 'CHAM_ELEM', repk=fieldOption)
        end if

        nbVari = 0
        if (fieldDisc .eq. "ELNO") then
            call jeveuo(fieldUser(1:19)//'.CELD', 'L', vi=celd)
            nbVari = CELD(4)
        end if
        WRITE (6, *) "CELD: ", exteVariName, nbVari

! ----- Generic checks
        if (fieldMesh .ne. mesh) then
            call utmess('F', 'VARC1_2')
        end if
        if (physQuan .ne. fieldPhysQuan) then
            valk(1) = exteVariName
            valk(2) = physQuan
            valk(3) = fieldPhysQuan
            call utmess('F', 'VARC1_3', nk=3, valk=valk)
        end if

! ----- Specific checks
        if (lFieldELGA) then
            if (fieldFED .ne. modelFED) then
                call utmess('F', 'VARC1_4', sk=exteVariName)
            end if
            if (fieldOption .ne. 'INI_SP_MATER') then
                valk(1) = exteVariName
                valk(2) = fieldOption
                call utmess('F', 'VARC1_5', nk=2, valk=valk)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module ExternalStateVariablePrep_module
