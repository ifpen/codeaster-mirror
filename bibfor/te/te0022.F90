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
subroutine te0022(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/sigvmc.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D*, C_PLAN*, D_PLAN*, AXIS*
!
! Options: SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    real(kind=8), parameter :: nharm = 0.d0
    integer(kind=8) :: ndim, nno, npg, nbsig, i
    integer(kind=8) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    integer(kind=8) :: jvSigm, jvDisp, jvGeom, jvMaterc
    real(kind=8) :: sigm(162), time
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=jvGaussWeight, jvf=jvBaseFunc, jdfde=jvDBaseFunc)

! - Initializations
    nbsig = nbsigm()
    sigm = 0.d0
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Current displacements (nodes)
    call jevech('PDEPLAR', 'L', jvDisp)

! - Get current time
    time = r8vide()

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Compute mechanical stress (without effect of external state variables)
    call sigvmc(materPara, &
                nno, ndim, nbsig, npg, &
                jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                zr(jvGeom), zr(jvDisp), &
                time, nharm, &
                sigm)

! - Final copy of stress
    call jevech('PCONTRR', 'E', jvSigm)
    do i = 1, nbsig*npg
        zr(jvSigm+i-1) = sigm(i)
    end do
!
end subroutine
