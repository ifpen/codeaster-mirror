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
subroutine nmarpr(result, sddisc, lReuse, numeInstEnd, timeEnd, &
                  numeStore)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmttch.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: result
    character(len=19), intent(in) :: sddisc
    aster_logical, intent(in) :: lReuse
    integer(kind=8), intent(in) :: numeInstEnd
    real(kind=8), intent(in) :: timeEnd
    integer(kind=8), intent(out) :: numeStore
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE *_NON_LINE (ARCHIVAGE)
!
! PREMIER NUMERO A ARCHIVER
!
! --------------------------------------------------------------------------------------------------
!
! IN  RESULT : NOM DE LA SD RESULTAT
! IN  SDDISC : SD DISCRETISATION
! IN  NUMDER : DERNIER NUMERO ARCHIVE DANS RESULT
!               (OU 0 SI NON REENTRANT)
! IN  INSDER : DERNIER INSTANT ARCHIVE DANS RESULT
!               (R8VIDE SI NON REENTRANT)
! IN  LREUSE : .TRUE. SI CONCEPT REENTRANT
! OUT NUMARC : NUMERO DU PREMIER PAS A ARCHIVER
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sddiscDitrJv
    real(kind=8) :: timeNext
    real(kind=8), pointer :: sddiscDitr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Access to list of time
    sddiscDitrJv = sddisc(1:19)//'.DITR'
    call jeveuo(sddiscDitrJv, 'L', vr=sddiscDitr)
!
    if (lReuse) then
        timeNext = sddiscDitr(2)
        if (timeNext .le. timeEnd) then
            call utmess('I', 'ARCHIVAGE_1', nr=2, valr=[timeEnd, timeNext])
            call nmttch(result, timeNext, numeInstEnd)
            numeStore = numeInstEnd
        else
            numeStore = numeInstEnd+1
        end if
!
    else
        ASSERT(numeInstEnd .eq. 0)
        numeStore = 0
    end if
!
    call jedema()
!
end subroutine
