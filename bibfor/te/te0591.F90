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
subroutine te0591(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nifnlg.h"
#include "asterfort/nifnpd.h"
#include "asterfort/nifnsm.h"
#include "asterfort/niinit.h"
#include "asterfort/utmess.h"
#include "MeshTypes_type.h"
!
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
! FONCTION REALISEE:  CALCUL DE L'OPTION FORC_NODA POUR LES ELEMENTS
!                     INCOMPRESSIBLES A 3 CHAMPS UGP
!                     EN 3D/D_PLAN/AXI
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno1, nno2, nno3, nnos, npg, jgn, ntrou
    integer(kind=8) :: iw, ivf1, ivf2, ivf3, idf1, idf2, idf3
    integer(kind=8) :: vu(3, MT_NNOMAX3D), vg(MT_NNOMAX3D), vp(MT_NNOMAX3D), vpi(3, MT_NNOMAX3D)
    integer(kind=8) :: igeom, jvSief, jvDisp, imate, ivectu
    character(len=8) :: lielrf(10), typmod(2)
    character(len=16), pointer :: compor(:) => null()
! ----------------------------------------------------------------------
!
! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elref2(nomte, 10, lielrf, ntrou)
    ASSERT(ntrou .ge. 3)
    call elrefe_info(elrefe=lielrf(3), fami='RIGI', ndim=ndim, nno=nno3, nnos=nnos, npg=npg, &
                     jpoids=iw, jvf=ivf3, jdfde=idf3, jgano=jgn)
    call elrefe_info(elrefe=lielrf(2), fami='RIGI', ndim=ndim, nno=nno2, nnos=nnos, npg=npg, &
                     jpoids=iw, jvf=ivf2, jdfde=idf2, jgano=jgn)
    call elrefe_info(elrefe=lielrf(1), fami='RIGI', ndim=ndim, nno=nno1, nnos=nnos, npg=npg, &
                     jpoids=iw, jvf=ivf1, jdfde=idf1, jgano=jgn)
!
! - TYPE DE MODELISATION
    if (ndim .eq. 2 .and. lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS  '
    else if (ndim .eq. 2 .and. lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN  '
    else if (ndim .eq. 3) then
        typmod(1) = '3D'
    else
        call utmess('F', 'ELEMENTS_34', sk=nomte)
    end if
!
! - Get index of dof
!
    call niinit(typmod, ndim, nno1, nno2, nno3, 0, vu, vg, vp, vpi)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PSIEFR', 'L', jvSief)
    call jevech('PDEPLAR', 'L', jvDisp)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PVECTUR', 'E', ivectu)
!
! - CALCUL DES FORCES INTERIEURES
    if (compor(DEFO) (1:6) .eq. 'PETIT ') then
!
        call nifnpd(ndim, nno1, nno2, nno3, npg, &
                    iw, zr(ivf1), zr(ivf2), zr(ivf3), idf1, &
                    vu, vg, vp, typmod, zr(igeom), &
                    zr(jvSief), zr(jvDisp), zr(ivectu))
    else if (compor(DEFO) (1:8) .eq. 'GDEF_LOG') then
!
        call jevech('PMATERC', 'L', imate)
        call nifnlg(ndim, nno1, nno2, nno3, npg, &
                    iw, zr(ivf1), zr(ivf2), zr(ivf3), idf1, &
                    idf2, vu, vg, vp, typmod, &
                    zi(imate), zr(igeom), zr(jvSief), zr(jvDisp), zr(ivectu))

    else if (compor(DEFO) (1:10) .eq. 'SIMO_MIEHE') then
!
        call jevech('PMATERC', 'L', imate)
        call nifnsm(ndim, nno1, nno2, nno3, npg, &
                    iw, zr(ivf1), zr(ivf2), zr(ivf3), idf1, &
                    idf2, vu, vg, vp, typmod, &
                    zi(imate), zr(igeom), zr(jvSief), zr(jvDisp), zr(ivectu))
!
    else
        call utmess('F', 'ELEMENTS3_16', sk=compor(DEFO))
    end if
!
end subroutine
