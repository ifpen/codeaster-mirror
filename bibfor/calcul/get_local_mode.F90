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

subroutine get_local_mode(option, nompar, typef, modeloc)

    implicit none

#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
#include "asterfort/modat2.h"

    character(len=16), intent(in) :: option
    character(len=8), intent(in) ::  nompar
    integer(kind=8), intent(in) :: typef
    integer(kind=8), intent(out) :: modeloc
!-----------------------------------------------------------------------
!   entrees:
!      nompar : nom d'1 parametre de l'option que l'on calcule.
!   sorties:
!      grdeur : grandeur associee au parametre
!-----------------------------------------------------------------------
    integer(kind=8) :: optionNume
!-------------------------------------------------------------------

    call jemarq()

    call jenonu(jexnom('&CATA.OP.NOMOPT', option), optionNume)

    modeloc = modat2(optionNume, typef, nompar)

    call jedema()

end subroutine
