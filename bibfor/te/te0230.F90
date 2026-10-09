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
subroutine te0230(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/bmatmc.h"
#include "asterfort/btdbmc.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/pmfmats.h"
#include "asterfort/rcvalb.h"
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
! Elements: 2D
! Option: RIGI_MECA_HYST
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    character(len=8), parameter :: fami = 'RIGI'
    real(kind=8), parameter :: nharm = 0.d0
    integer(kind=8), parameter :: nbProp = 2, nbPara = 3
    character(len=8), parameter :: paraName(nbPara) = (/'X', 'Y', 'Z'/)
    real(kind=8) :: propVale(nbProp)
    character(len=16) :: propName(nbProp)
    integer(kind=8) :: propCode(nbProp)
    integer(kind=8) :: i, kpg, jvMaterc, imatuu, j, iret, jvRigiReal(2), rigi
    integer(kind=8) :: k, nbinco, nbsig, ndim, nno, npg
    integer(kind=8) :: jvGeom, ipoids, ivf, idfde
    real(kind=8) :: b(6*MAX_BV_CG), jacgau
    real(kind=8) :: btdbi(81, 81), di(36), eta
    real(kind=8) :: time, coorBary(3)
    character(len=8) :: materParaKpg
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)

! - Initializations
    nbsig = nbsigm()
    nbinco = ndim*nno
    btdbi = 0.d0
    coorBary = 0.d0

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get current time
    time = r8vide()

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Get material for integration point
    call pmfmats(materParaKpg)

! - Compute barycenter of cell from coordinates of nodes
    call compCellBary(ndim, nno, jvGeom, coorBary)

! - Get RIGI_MECA real part
    call tecach('ONO', 'PRIGIEL', 'L', iret, nval=2, itab=jvRigiReal)

! - Compute RIGI_MECA imaginary part
    if (lMaterVisc(materPara)) then
        do kpg = 1, npg
! --------- Initializations of material parameters on current integration point
            call initParaPoin(kpg, ksp, materPara)

! --------- Compute matrix [B]: displacement -> strain (first order)
            call bmatmc(kpg, nbsig, zr(jvGeom), ipoids, ivf, &
                        idfde, nno, nharm, jacgau, b)

! --------- Compute Hooke matrix [D]
            call dmatmc(materPara, "+", time, &
                        nbsig, di_=di)

! --------- Compute rigidity matrix [K] = [B]Tx[D]x[B]
            call btdbmc(b, di, jacgau, ndim, nno, &
                        nbsig, materPara%elasID, btdbi)
        end do
    else
        propName(1) = 'AMOR_HYST'
        propVale(1) = 0.d0
        kpg = 1
        call initParaPoin(kpg, ksp, materPara)
        call rcvalb(fami, kpg, ksp, '+', &
                    zi(jvMaterc), materParaKpg, materPara%elasKeyword, &
                    ndim, paraName, coorBary, &
                    1, propName, propVale, propCode, 0, &
                    nan='NON')
        eta = propVale(1)
    end if

! - Set matrix in output field
    rigi = jvRigiReal(1)
    call jevech('PMATUUC', 'E', imatuu)
    k = 0
    do i = 1, nbinco
        do j = 1, i
            k = k+1
            if (lMaterVisc(materPara)) then
                zc(imatuu+k-1) = dcmplx(zr(rigi+k-1), btdbi(i, j))
            else
                zc(imatuu+k-1) = dcmplx(zr(rigi+k-1), eta*zr(rigi+k-1))
            end if
        end do
    end do
!
end subroutine
