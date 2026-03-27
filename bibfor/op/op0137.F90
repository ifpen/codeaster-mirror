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

subroutine op0137()
    implicit none
!     OPERATEUR :     DEBUG
#include "asterc/jdcset.h"
#include "asterfort/getvtx.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: lundef, idebug
    common/undfje/lundef, idebug
! ----------------------------------------------------------------------
    character(len=3) :: repons
    integer(kind=8) :: l
!
!
!
!     SDVERI=OUI/NON :
!     ----------------
    repons = ' '
    call getvtx(' ', 'SDVERI', scal=repons, nbret=l)
    if (repons .eq. 'OUI') then
        call jdcset('sdveri', 1)
        call utmess('I', 'SUPERVIS_24')
    else if (repons .eq. 'NON') then
        call jdcset('sdveri', 0)
        call utmess('I', 'SUPERVIS_43', sk='SDVERI')
    end if
!
!
!     JEVEUX=OUI/NON :
!     ----------------
    repons = ' '
    call getvtx(' ', 'JEVEUX', scal=repons, nbret=l)
    if (repons .eq. 'OUI') then
        idebug = 1
        call utmess('I', 'SUPERVIS_44', sk='JEVEUX')
        call jdcset('dbgjeveux', 1)
    else if (repons .eq. 'NON') then
        idebug = 0
        call utmess('I', 'SUPERVIS_43', sk='JEVEUX')
        call jdcset('dbgjeveux', 0)
    end if
!
!
!     JXVERI=OUI/NON :
!     ----------------
    repons = ' '
    call getvtx(' ', 'JXVERI', scal=repons, nbret=l)
    if (repons .eq. 'OUI') then
        call utmess('I', 'SUPERVIS_23')
        call jdcset('jxveri', 1)
    else if (repons .eq. 'NON') then
        call utmess('I', 'SUPERVIS_43', sk='JXVERI')
        call jdcset('jxveri', 0)
    end if
!
!
!     IMPR_MACRO=OUI/NON :
!     ---------------------
    repons = ' '
    call getvtx(' ', 'IMPR_MACRO', scal=repons, nbret=l)
    if (repons .eq. 'OUI') then
        call utmess('I', 'SUPERVIS_44', sk='IMPR_MACRO')
        call jdcset('ShowChildCmd', 1)
    else if (repons .eq. 'NON') then
        call utmess('I', 'SUPERVIS_43', sk='IMPR_MACRO')
        call jdcset('ShowChildCmd', 0)
    end if
!
!
end subroutine
