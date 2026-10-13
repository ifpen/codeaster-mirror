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
subroutine te0597(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/niinit.h"
#include "asterfort/nurmtd.h"
#include "asterfort/teattr.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! FONCTION REALISEE:  CALCUL DE LA RIGIDITE MECANIQUE POUR LES ELEMENTS
!                     INCOMPRESSIBLES A 2 CHAMPS UP
!                     EN 3D/D_PLAN/AXI
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    aster_logical :: mini
    integer(kind=8) :: ndim, nno1, nno2, npg, ntrou
    integer(kind=8) :: iw, ivf1, ivf2, idf1, idf2
    integer(kind=8) :: vu(3, MT_NNOMAX3D), vg(MT_NNOMAX3D), vp(MT_NNOMAX3D), vpi(3, MT_NNOMAX3D)
    integer(kind=8) :: jvGeom, jvMaterc, jvMatr
    integer(kind=8) :: ibid
    character(len=8) :: lielrf(10), typmod(2), alias8
!
! --------------------------------------------------------------------------------------------------
!

! - MINI ELEMENT ?
    call teattr('S', 'ALIAS8', alias8, ibid)
    if (alias8(6:8) .eq. 'TR3' .or. alias8(6:8) .eq. 'TE4') then
        mini = .true.
    else
        mini = .false.
    end if

! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elref2(nomte, 10, lielrf, ntrou)
    ASSERT(ntrou .ge. 2)
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
                nno1, 0, nno2, 0, &
                vu, vg, vp, vpi)

! - Input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PMATERC', 'L', jvMaterc)
    call jevech('PMATUUR', 'E', jvMatr)

! - Compute rigidity matrix
    call nurmtd(ndim, nno1, nno2, npg, iw, &
                zr(ivf1), zr(ivf2), idf1, vu, &
                vp, typmod, jvGeom, zi(jvMaterc), mini, &
                zr(jvMatr))
!
end subroutine
