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
subroutine nmcrpp(factorKeywordZ, iFactorKeyword, stepSlctTole)
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: factorKeywordZ
    integer(kind=8), intent(in) :: iFactorKeyword
    real(kind=8), intent(out) :: stepSlctTole
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Time selector management
!
! Get tolerance to select one time step
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword to read
! In  iFactorKeyword   : index of factor keyword
! Out stepSlctTole     : tolerance to select one time step
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: n0, n1, n2
    character(len=16) :: factorKeyword
    real(kind=8) :: predef
    character(len=8) :: criter
    real(kind=8) :: prec
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    prec = 0.d0
    stepSlctTole = 0.d0
    factorKeyword = factorKeywordZ
    predef = 1.d-6

!   CRITERE/PRECISION are only needed if INST or LIST_INST exist
    call getvr8(factorKeyword, 'INST', iocc=iFactorKeyword, nbret=n0)
    if (n0 .eq. 0) then
        call getvid(factorKeyword, 'LIST_INST', iocc=iFactorKeyword, nbret=n0)
    end if

    if (n0 .ne. 0) then
        call getvr8(factorKeyword, 'PRECISION', iocc=iFactorKeyword, scal=prec, nbret=n1)
        call getvtx(factorKeyword, 'CRITERE', iocc=iFactorKeyword, scal=criter, nbret=n2)
        if (criter .eq. 'ABSOLU') then
            if (n1 .eq. 0) then
                call utmess('F', 'LISTINST_1')
            end if
        else if (criter .eq. 'RELATIF') then
            if (n1 .eq. 0) then
                prec = predef
                call utmess('A', 'LISTINST_2', sr=predef)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
        if (prec .le. r8prem()) then
            call utmess('F', 'LISTINST_3')
        end if
        if (criter .eq. 'RELATIF') then
            stepSlctTole = prec
        else if (criter .eq. 'ABSOLU') then
            stepSlctTole = -prec
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
    call jedema()
!
end subroutine
