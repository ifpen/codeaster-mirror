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
subroutine nueffe(nbLigr, listLigr, base, numeDofZ, renumZ, &
                  modelZ, modeLocZ_, idenRelaZ_, lMacrElemZ_)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/cheksd.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/nueffe_lag1.h"
#include "asterfort/nueffe_lag2.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: nbLigr
    character(len=24), pointer :: listLigr(:)
    character(len=2), intent(in) :: base
    character(len=*), intent(in) :: numeDofZ, renumZ, modelZ
    character(len=*), optional, intent(in) :: modeLocZ_, idenRelaZ_
    aster_logical, optional, intent(in) :: lMacrElemZ_
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering - Create NUME_EQUA objects
!
! --------------------------------------------------------------------------------------------------
!
! In  nbLigr         : number of LIGREL in list
! Ptr listLigr       : pointer to list of LIGREL
! In  numeDof        : name of numbering object (NUME_DDL)
! In  base           : JEVEUX base to create objects
!                      base(1:1) => NUME_DDL objects
!                      base(2:2) => NUME_EQUA objects
! In  renum          : method for renumbering equations (SANS/RCMK)
! In  model          : name of model
! In  modeLoc        : local mode for GRANDEUR numbering
! In  idenRela       : name of object for identity relations between dof
! In  lMacrElem      : call from MACR_ELEM_STAT
!
! --------------------------------------------------------------------------------------------------
!
! Attention : ne fait pas jemarq/jedema car nulili
!             recopie des adresses jeveux dans .ADNE et .ADLI
!             Ces objets seront utilises pendant la creation de la sd "stockage" (promor.F90)
!
! --------------------------------------------------------------------------------------------------
!
! Cette routine cree les objets suivants :
!  nume(1:14)//     .ADLI
!                   .ADNE
!              .NUME.DEEQ
!              .NUME.DELG
!              .NUME.LILI
!              .NUME.NEQU
!              .NUME.NUEQ
!              .NUME.PRNO
!              .NUME.REFN
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: debug = ASTER_FALSE
    integer(kind=8) :: iLigr, iret
    character(len=8) :: typeLagr, model, typeLagrC
    character(len=24) :: modeLoc, idenRela
    character(len=19) :: ligrelName
    aster_logical :: lMacrElem
!
! --------------------------------------------------------------------------------------------------
!
!    call jemarq() FORBIDDEN !

    modeLoc = " "
    if (present(modeLocZ_)) then
        modeLoc = modeLocZ_
    end if
    idenRela = " "
    if (present(idenRelaZ_)) then
        idenRela = idenRelaZ_
    end if
    lMacrElem = ASTER_FALSE
    if (present(lMacrElemZ_)) then
        lMacrElem = lMacrElemZ_
    end if
    model = modelZ

! - Only one method for Lagrange
! - The first LIGREL is for model => no Lagrange on it !
! - No double Lagrange on contact virtual cells
    typeLagr = " "
    do iLigr = 2, nbLigr
        ligrelName = listLigr(iLigr) (1:19)
        call exisd('LIGREL', ligrelName, iret)
        if (iret .ne. 0) then
            call dismoi('TYPE_LAGR', ligrelName, 'LIGREL', repk=typeLagrC)
            if (typeLagr .eq. " ") then
                typeLagr = typeLagrC
            end if
            if (typeLagr .ne. typeLagrC) then
                call utmess('F', 'ASSEMBLA_6')
            end if
        end if
    end do

! - Create
    if (typeLagr .eq. 'LAG1') then
        if (lMacrElem) call utmess('F', 'SOUSTRUC_72')
        ASSERT(nbLigr > 1)
! ----- Case with simple Lagrange
        call nueffe_lag1(nbLigr, listLigr, base, numeDofZ, renumZ, &
                         model, modeLoc, idenRela)
    else
! ----- Case without Lagrange or with double Lagrange
        call nueffe_lag2(nbLigr, listLigr, base, numeDofZ, renumZ, &
                         model, modeLoc, idenRela)
    end if

! - Print debug
    if (debug) then
        call cheksd(numeDofZ, 'SD_NUME_DDL', iret)
    end if

!    call jedema() FORBIDDEN !
!
end subroutine
