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
subroutine numero(numeDofZ, base, &
                  nbLigr, listLigr, &
                  numeDofOldZ_, modeLocZ_, &
                  idenRelaZ_, modelZ_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/crnulg.h"
#include "asterfort/dismoi.h"
#include "asterfort/gettco.h"
#include "asterfort/numer2.h"
#include "asterfort/uttcpu.h"
!
    character(len=*), intent(inout) :: numeDofZ
    character(len=2), intent(in) :: base
    integer(kind=8), intent(in) :: nbLigr
    character(len=24), pointer :: listLigr(:)
    character(len=*), optional, intent(in) :: numeDofOldZ_
    character(len=*), optional, intent(in) :: modeLocZ_, idenRelaZ_, modelZ_
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering
!
! --------------------------------------------------------------------------------------------------
!
! IO  numeDof        : name of numbering object (NUME_DDL)
! In  base           : JEVEUX base to create objects
!                      base(1:1) => NUME_EQUA objects
!                      base(2:2) => NUME_DDL objects
! In  nbLigr         : number of LIGREL in list
! Ptr listLigr       : pointer to list of LIGREL
! In  numeDofOld     : name of previous numeDof object
! In  modeLoc        : local mode for GRANDEUR numbering
! In  idenRela       : name of object for identity relations between dof
! In  model          : model
!
! If numeDofOld is present
!   -> try to know if NUME_EQUA in numeDofOld can be reuse
!      In this case numeDof = numeDofOld
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iLigr
    character(len=8) :: mesh, model, modelInList
    character(len=14) :: numeDof
    character(len=16) :: typsd
    character(len=24) :: modeLoc, numeDofOld, idenRela, ligrelName
!
! --------------------------------------------------------------------------------------------------
!
    call uttcpu('CPU.RESO.1', 'DEBUT', ' ')
    call uttcpu('CPU.RESO.2', 'DEBUT', ' ')

    ASSERT(nbligr .ge. 1)

! - Initializations
    idenRela = ' '
    if (present(idenRelaZ_)) then
        idenRela = idenRelaZ_
    end if
    modeLoc = ' '
    if (present(modeLocZ_)) then
        modeLoc = modeLocZ_
    end if
    numeDofOld = ' '
    if (present(numeDofOldZ_)) then
        numeDofOld = numeDofOldZ_
    end if

! ==================================================================================================
!
!  Affreuse glute (tant qu'on stockera le modèle dans NUME_EQUA/REFN)
!
! ==================================================================================================
! - Detect model in list of LIGREL
    modelInList = " "
    do iLigr = 1, nbLigr
        ligrelName = listLigr(iLigr)
        if (ligrelName(10:15) .eq. "MODELE") then
            modelInList = ligrelName(1:8)
        end if
    end do

! - Set model
    model = " "
    if (present(modelZ_)) then
        model = modelZ_
        if (modelInList .ne. ' ') then
            ASSERT(modelInList .eq. model)
        end if
    else
        model = modelInList
    end if
! ==================================================================================================
! ==================================================================================================

! - Create numbering
    call numer2(nbLigr, listLigr, base, numeDofZ, &
                numeDofOld, modeLoc, model, idenRela)

! - Create parallel numbering
    call dismoi('NOM_MAILLA', listLigr(1), 'LIGREL', repk=mesh)
    call gettco(mesh, typsd)
    if (typsd .eq. 'MAILLAGE_P') then
        numeDof = numeDofZ
        call crnulg(numeDof, to_aster_logical(modeLoc .eq. ' '))
    end if
!
    call uttcpu('CPU.RESO.1', 'FIN', ' ')
    call uttcpu('CPU.RESO.2', 'FIN', ' ')
!
end subroutine
