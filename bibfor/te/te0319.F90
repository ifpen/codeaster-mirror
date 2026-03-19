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
subroutine te0319(option, nomte)
!
    use Metallurgy_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements:
! Option: META_TRAN_ELNO
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbNode, iNode
    integer(kind=8) :: nbVariIn, nbVariOut, iVariIn, iVariOut
    integer(kind=8) :: jvPhaseIn, jvPhaseOut
    integer(kind=8) :: itab(7), iret, indxIn
    real(kind=8) :: metaIn(META_META_NBPHASE_MAXI), metaOut(META_META_NBPHASE_MAXI)
    integer(kind=8), parameter :: indxTransfer(NBVARISTEELR) = (/1, 2, 3, 0, &
                                                                 4, 0, 5, 6, &
                                                                 7, 8, 9, 0/)
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nbNode)

! - Get input and output fields
    call tecach("ONN", "PPHASIN", 'L', iret, nval=7, itab=itab)
    ASSERT(iret .eq. 0)
    jvPhaseIn = itab(1)
    nbVariIn = itab(6)
    ASSERT(nbVariIn .le. META_META_NBPHASE_MAXI)
    ASSERT(nbVariIn .gt. 0)
    call tecach("ONN", "PPHASOUT", 'E', iret, nval=7, itab=itab)
    ASSERT(iret .eq. 0)
    jvPhaseOut = itab(1)
    nbVariOut = itab(6)
    ASSERT(nbVariOut .le. META_META_NBPHASE_MAXI)
    ASSERT(nbVariOut .gt. 0)

    if (nbVariIn .lt. nbVariOut) then
        ASSERT(nbVariIn .eq. NBVARISTEEL)
        ASSERT(nbVariOut .eq. NBVARISTEELR)
        do iNode = 1, nbNode
            metaIn = 0.d0
            do iVariIn = 1, nbVariIn
                metaIn(iVariIn) = zr(jvPhaseIn+nbVariIn*(iNode-1)+iVariIn-1)
            end do
            metaOut = 0.d0
            do iVariOut = 1, nbVariOut
                indxIn = indxTransfer(iVariOut)
                if (indxIn .ne. 0) then
                    metaOut(iVariOut) = metaIn(indxIn)
                end if
            end do
            do iVariOut = 1, nbVariOut
                zr(jvPhaseOut+nbVariOut*(iNode-1)+iVariOut-1) = metaOut(iVariOut)
            end do
        end do
    end if

!
end subroutine
