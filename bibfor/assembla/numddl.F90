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
subroutine numddl(numeDofZ, renumZ, base, nbMatrElem, listMatrElem, &
                  lMacrElemZ_)
!
    implicit none
!
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/crnulg.h"
#include "asterfort/dismoi.h"
#include "asterfort/gettco.h"
#include "asterfort/nueffe.h"
#include "asterfort/numoch.h"
#include "asterfort/utmess.h"
!
    character(len=2), intent(in) :: base
    character(len=*), intent(in) :: numeDofZ, renumZ
    character(len=24), pointer :: listMatrElem(:)
    integer(kind=8), intent(in) :: nbMatrElem
    aster_logical, optional, intent(in) :: lMacrElemZ_
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Create numbering from list of elementary matrices
!
! --------------------------------------------------------------------------------------------------
!
! In  numeDof        : name of numeDof object
! In  renum          : method for renumbering equations (SANS/RCMK)
! In  base           : JEVEUX base to create objects
!                      base(1:1) => NUME_EQUA objects
!                      base(2:2) => NUME_DDL objects
! Ptr listMatrElem   : list of elementary matrixes
! In  nbMatrElem     : number of elementary matrixes
! In  lMacrElem      : call from MACR_ELEM_STAT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbLigr
    character(len=24), pointer :: listLigr(:) => null()
    integer(kind=8) :: iMatrElem
    character(len=8) :: mesh, model, modelNew
    character(len=14) :: numeDof
    character(len=16) :: typsd
    aster_logical :: lMacrElem
!
! --------------------------------------------------------------------------------------------------
!
    lMacrElem = ASTER_FALSE
    if (present(lMacrElemZ_)) then
        lMacrElem = lMacrElemZ_
    end if

! - Extract list of LIGREL from elementary matrixes
    call numoch(listMatrElem, nbMatrElem, listLigr, nbLigr)

! - Get model
    ASSERT(nbMatrElem .ge. 1)
    call dismoi('NOM_MODELE', listMatrElem(1), 'MATR_ELEM', repk=model)
    do iMatrElem = 2, nbMatrElem
        call dismoi('NOM_MODELE', listMatrElem(1), 'MATR_ELEM', repk=modelNew)
        if (modelNew .ne. model) then
            call utmess("F", "ASSEMBLA_2")
        end if
    end do

! - Numbering - Create NUME_EQUA objects
    call nueffe(nbLigr, listLigr, base, numeDofZ, renumZ, model, &
                lMacrElemZ_=lMacrElem)
    AS_DEALLOCATE(vk24=listLigr)

! - Numbering - Create parallel object
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call gettco(mesh, typsd)
    if (typsd .eq. 'MAILLAGE_P') then
        numeDof = numeDofZ
        call crnulg(numeDof)
    end if
!
end subroutine
