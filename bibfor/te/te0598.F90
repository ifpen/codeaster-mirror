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
subroutine te0598(option, nomte)
!
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/niinit.h"
#include "asterfort/norfpd.h"
#include "asterfort/nurfgd.h"
#include "asterfort/nurfpd.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
#include "MeshTypes_type.h"
!
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
! FONCTION REALISEE:  CALCUL DE L'OPTION FORC_REFE POUR LES ELEMENTS
!                     INCOMPRESSIBLES A 2 CHAMPS UP
!                     EN 3D/D_PLAN/AXI
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno1, nno2, nnos, npg, jgn, ntrou
    integer(kind=8) :: iw, ivf1, ivf2, idf1, idf2
    integer(kind=8) :: vu(3, MT_NNOMAX3D), vg(MT_NNOMAX3D), vp(MT_NNOMAX3D), vpi(3, MT_NNOMAX3D)
    integer(kind=8) :: igeom, ivectu
    real(kind=8) :: sigref, epsref, piref
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: lielrf(10), typmod(2)
    character(len=24) :: valk
    type(RESI_REFE):: refe
! ----------------------------------------------------------------------
!
!
! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elref2(nomte, 10, lielrf, ntrou)
    ASSERT(ntrou .ge. 2)
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
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PVECTUR', 'E', ivectu)
    call jevech('PCOMPOR', 'L', vk16=compor)

    call refe%Init(nomte)
    sigref = refe%GetRef('SIGM')
    epsref = refe%GetRef('EPSI')
    call refe%Check()
!
! - CALCUL DE REFE_FORC_NODA
    if (compor(DEFO) .eq. 'PETIT ') then
!
        if (lteatt('INCO', 'C2 ')) then
! --------- Get index of dof
            call niinit(typmod, ndim, nno1, 0, nno2, 0, vu, vg, vp, vpi)
!
            call nurfpd(ndim, nno1, nno2, npg, iw, zr(ivf1), zr(ivf2), idf1, vu, vp, &
                        typmod, zr(igeom), sigref, epsref, zr(ivectu))
        elseif (lteatt('INCO', 'C2O')) then
            piref = refe%GetRef('PI')
! --------- Get index of dof
            call niinit(typmod, ndim, nno1, 0, nno2, nno2, vu, vg, vp, vpi)
!
            call norfpd(ndim, nno1, nno2, nno2, npg, iw, zr(ivf1), zr(ivf2), zr(ivf2), idf1, &
                        vu, vp, vpi, typmod, nomte, zr(igeom), sigref, epsref, piref, zr(ivectu))
        else
            valk = compor(DEFO)
            call utmess('F', 'MODELISA10_17', sk=valk)
        end if
    else if (compor(DEFO) .eq. 'GDEF_LOG') then
!
        if (.not. lteatt('INCO', 'C2 ')) then
            valk = compor(DEFO)
            call utmess('F', 'MODELISA10_17', sk=valk)
        end if
! ----- Get index of dof
        call niinit(typmod, ndim, nno1, 0, nno2, 0, vu, vg, vp, vpi)
!
        call nurfgd(ndim, nno1, nno2, npg, iw, zr(ivf1), zr(ivf2), idf1, vu, vp, &
                    typmod, zr(igeom), sigref, epsref, zr(ivectu))
    else
        call utmess('F', 'ELEMENTS3_16', sk=compor(DEFO))
    end if
    call refe%Check()
!
end subroutine
