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
subroutine te0551(option, nomte)
!
    use Metallurgy_type
    use MetallurgySteel_Compute_module
    use MetallurgyMeca_module, only: metaGetType
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/Metallurgy_type.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: THERMIQUE - 3D*, AXIS*, PLAN*
! Option: DURT_ELNO
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iPhase, jvDurt, jvPhaseIn, iret
    integer(kind=8) :: nbNode, jvMaterc, iNode, nbVari, nbPhase, metaType
    real(kind=8) :: phase(PRSTEEL_NB), durtno
    integer(kind=8) :: jtab(6)
    type(META_HardnessParameters) :: metaHardnessPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nbNode)
    ASSERT(nbNode .le. MT_NNOMAX3D)

! - Input/Output fields
    call tecach('OOO', 'PPHASIN', 'L', iret, nval=6, itab=jtab)
    nbVari = jtab(6)
    jvPhaseIn = jtab(1)
    call jevech('PMATERC', 'L', jvMaterc)
    call jevech('PDURT_R', 'E', jvDurt)

! - Type of phases
    if (nbVari .eq. NBVARISTEEL) then
        nbPhase = PSTEEL_NB
        metaType = META_STEEL
    else if (nbVari .eq. NBVARISTEELR) then
        nbPhase = PRSTEEL_NB
        metaType = META_STEEL_TEMPER
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get material properties
    call metaHardnessGetParameters(zi(jvMaterc), metaType, nbPhase, metaHardnessPara)

! - Compute
    do iNode = 1, nbNode
! ----- Get phases
        do iPhase = 1, nbPhase
            phase(iPhase) = zr(jvPhaseIn+nbVari*(iNode-1)+iPhase-1)
        end do

! ----- Compute hardness
        durtno = 0.d0
        do iPhase = 1, nbPhase
            durtno = durtno+phase(iPhase)*metaHardnessPara%hardSteel(iPhase)
        end do
        zr(jvDurt+(iNode-1)) = durtno
    end do
!
end subroutine
