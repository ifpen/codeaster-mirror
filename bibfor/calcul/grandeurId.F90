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

subroutine grandeurId(option, nompar, gd)

    implicit none

#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterc/indik8.h"
#include "asterfort/utmess.h"

    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: nompar
    integer(kind=8), intent(out) :: gd
!-----------------------------------------------------------------------
!   entrees:
!      nompar : nom d'1 parametre de l'option que l'on calcule.
!   sorties:
!      grdeur : grandeur associee au parametre
!-----------------------------------------------------------------------
    integer(kind=8) :: jpar, nbpar, optionNume
    integer(kind=8), pointer :: cataDescopt(:) => null()
    character(len=8), pointer :: cataOptpara(:) => null()
!-------------------------------------------------------------------

    call jemarq()

    call jenonu(jexnom('&CATA.OP.NOMOPT', option), optionNume)
    call jeveuo(jexnum('&CATA.OP.DESCOPT', optionNume), 'L', vi=cataDescopt)
    call jeveuo(jexnum('&CATA.OP.OPTPARA', optionNume), 'L', vk8=cataOptpara)

    nbpar = cataDescopt(2)+cataDescopt(3)+cataDescopt(4)
    jpar = indik8(cataOptpara, nompar, 1, nbpar)
    if (jpar .eq. 0) then
        call utmess('F', 'CALCUL_14', sk=nompar)
    end if
    gd = cataDescopt(4+jpar)

    call jedema()

end subroutine
