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
subroutine te0592(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/niinit.h"
#include "asterfort/nirmtd.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! FONCTION REALISEE:  CALCUL DE LA RIGIDITE MECANIQUE POUR LES ELEMENTS
!                     INCOMPRESSIBLES A 3 CHAMPS UGP
!                     EN 3D/D_PLAN/AXI
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8) :: ndim, nno1, nno2, nno3, npg, ntrou
    integer(kind=8) :: iw, ivf1, ivf2, ivf3, idf1, idf2, idf3
    integer(kind=8) :: vu(3, MT_NNOMAX3D), vg(MT_NNOMAX3D), vp(MT_NNOMAX3D), vpi(3, MT_NNOMAX3D)
    integer(kind=8) :: jvGeom, jvMaterc, jvMatr
    character(len=8) :: lielrf(10), typmod(2)
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call elref2(nomte, 10, lielrf, ntrou)
    ASSERT(ntrou .ge. 3)
    call elrefe_info(elrefe=lielrf(3), fami=fami, ndim=ndim, nno=nno3, &
                     jvf=ivf3, jdfde=idf3)
    call elrefe_info(elrefe=lielrf(2), fami=fami, ndim=ndim, nno=nno2, &
                     jvf=ivf2, jdfde=idf2)
    call elrefe_info(elrefe=lielrf(1), fami=fami, ndim=ndim, nno=nno1, npg=npg, &
                     jpoids=iw, jvf=ivf1, jdfde=idf1)

! - TYPE DE MODELISATION
    if (ndim .eq. 2 .and. lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS  '
    else if (ndim .eq. 2 .and. lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN  '
    else if (ndim .eq. 3) then
        typmod(1) = '3D'
    else
        ASSERT(ASTER_FALSE)
    end if
    typmod(2) = '        '

! - Get index of dof
    call niinit(typmod, ndim, &
                nno1, nno2, nno3, 0, &
                vu, vg, vp, vpi)

! - Input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PMATERC', 'L', jvMaterc)
    call jevech('PMATUUR', 'E', jvMatr)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno1, jvGeom, materPara%lcsPara)

! - Compute rigidity matrix
    call nirmtd(ndim, nno1, nno2, nno3, npg, &
                iw, zr(ivf2), zr(ivf3), ivf1, idf1, &
                vu, vg, vp, jvGeom, materPara, &
                zr(jvMatr))
!
end subroutine
