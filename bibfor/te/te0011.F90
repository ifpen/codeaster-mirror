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
subroutine te0011(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/bmatmc.h"
#include "asterfort/btdbmc.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/tecach.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
! Option: RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    character(len=8), parameter :: fami = 'RIGI'
    real(kind=8), parameter :: nharm = 0.d0
    integer(kind=8) :: jvGeom, jvMaterc, imatuu, jvInstr
    integer(kind=8) :: i, j, k, kpg, iret
    integer(kind=8) :: nbinco, nbsig
    real(kind=8) :: b(6*MAX_BV_CG), btdb(81, 81), d(36), jacgau
    real(kind=8) :: time
    integer(kind=8) :: ndim, nno, npg
    integer(kind=8) :: ipoids, ivf, idfde
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(ndim .eq. 3)

! - Initializations
    nbsig = nbsigm()
    nbinco = ndim*nno
    btdb = 0.d0

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get current time
    time = r8vide()
    call tecach('ONO', 'PINSTR', 'L', iret, iad=jvInstr)
    if (jvInstr .ne. 0) then
        time = zr(jvInstr)
    end if

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Compute RIGI_MECA
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! ----- Compute matrix [B]: displacement -> strain (first order)
        call bmatmc(kpg, nbsig, zr(jvGeom), ipoids, ivf, &
                    idfde, nno, nharm, jacgau, b)

! ----- Compute Hooke matrix [D]
        call dmatmc(materPara, "+", time, &
                    nbsig, d)

! ----- Compute rigidity matrix [K] = [B]Tx[D]x[B]
        call btdbmc(b, d, jacgau, ndim, nno, &
                    nbsig, materPara%elasID, btdb)

    end do

! - Set matrix in output field
    call jevech('PMATUUR', 'E', imatuu)
    k = 0
    do i = 1, nbinco
        do j = 1, i
            k = k+1
            zr(imatuu+k-1) = btdb(i, j)
        end do
    end do
!
end subroutine
