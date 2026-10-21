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
subroutine te0013(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/sigtmc.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 2D et 3D
!
! Option: CHAR_MECA_* for external state variables
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    real(kind=8), parameter :: nharm = 0.d0
    real(kind=8) :: forcVarc(3*MT_NNOMAX), sigmVarc(162)
    real(kind=8) :: time
    integer(kind=8) :: jvTime, jvVect, iret
    integer(kind=8) :: ndim, nno, npg, nbsig, i, indxVarcStrain
    integer(kind=8) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    integer(kind=8) :: jvGeom, jvMaterc
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=jvGaussWeight, jvf=jvBaseFunc, jdfde=jvDBaseFunc)
    nbsig = nbsigm()
    ASSERT(nno .le. MT_NNOMAX)
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)

! - Choice of external state variable
    if (option .eq. 'CHAR_MECA_TEMP_R') then
        indxVarcStrain = VARC_STRAIN_TEMP
    elseif (option .eq. 'CHAR_MECA_HYDR_R') then
        indxVarcStrain = VARC_STRAIN_HYDR
    elseif (option .eq. 'CHAR_MECA_SECH_R') then
        indxVarcStrain = VARC_STRAIN_SECH
    elseif (option .eq. 'CHAR_MECA_PTOT_R') then
        indxVarcStrain = VARC_STRAIN_PTOT
    elseif (option .eq. 'CHAR_MECA_EPSA_R') then
        indxVarcStrain = VARC_STRAIN_EPSA
    else
        ASSERT(ASTER_FALSE)
    end if

! - Compute CHAR_MECA_TEMP_R for metallurgy
    ! if (ndim .eq. 3) then
    !     call metau2(l_meta)
    ! else
    !     call metau1(l_meta)
    ! end if

    ! if (.not. l_meta) then
! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! ----- Material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! ----- Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! ----- Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Get time
    time = r8vide()
    call tecach('ONO', 'PINSTR', 'L', iret, iad=jvTime)
    if (jvTime .ne. 0) then
        time = zr(jvTime)
    end if

! ----- Calcul des contraintes anélastiques
    call sigtmc(materPara, time, &
                nbsig, npg, ndim, &
                indxVarcStrain, sigmVarc)

! - Compute [B]Tx{SIGMVARC}
    call bsigmc(nno, ndim, nbsig, npg, &
                jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                zr(jvGeom), nharm, sigmVarc, &
                forcVarc)

! - Set output vector
    call jevech('PVECTUR', 'E', jvVect)
    do i = 1, ndim*nno
        zr(jvVect+i-1) = forcVarc(i)
    end do
!
    ! end if
end subroutine
