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

subroutine nmvcd2(varcNameZ, materFieldZ, exis_varc)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
    character(len=*), intent(in) :: varcNameZ
    character(len=*), intent(in) :: materFieldZ
    aster_logical, intent(out) :: exis_varc
!
! --------------------------------------------------------------------------------------------------
!
! External state variables
!
! Is external state variable exists in material field ?
!
! --------------------------------------------------------------------------------------------------
!
! In  varcName         : name of external state variable
! In  materField       : name of material characteristics (field)
! Out exis_varc        : ASTER_TRUE if state variable exists in material field
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nmax, i, iret
    character(len=8) :: materField
    character(len=8) :: varcName
    character(len=8), pointer :: cvrcvarc(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    materField = materFieldZ
    varcName = varcNameZ
    exis_varc = ASTER_FALSE
!
    call jeexin(materField//'.CVRCVARC', iret)
    if (iret .ne. 0) then
        call jelira(materField//'.CVRCVARC', 'LONMAX', ival=nmax)
        call jeveuo(materField//'.CVRCVARC', 'L', vk8=cvrcvarc)
        do i = 1, nmax
            if (cvrcvarc(i) .eq. varcName) then
                exis_varc = ASTER_TRUE
                exit
            end if
        end do
    end if
!
    call jedema()
end subroutine
