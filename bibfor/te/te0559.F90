! --------------------------------------------------------------------------------------------------
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
! --------------------------------------------------------------------------------------------------

subroutine te0559(option, nomte)
    use resi_refe_module, only: RESI_REFE
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
    character(len=16) :: option, nomte
! --------------------------------------------------------------------------------------------------
! REALISE LES OPTIONS :
!     REFE_FORC_NODA pour les éléments de BARRE
! --------------------------------------------------------------------------------------------------
! IN OPTION    : K16 :  OPTION DE CALCUL
! IN NOMTE     : K16 : NOM DU TYPE ELEMENT
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: iret, itab(2), jv_vectur, jv_sect, nddl
    real(kind=8):: sigref, aire, forref
    type(RESI_REFE):: refe
! --------------------------------------------------------------------------------------------------

    ! Cross section area
    call jevech('PCAGNBA', 'L', jv_sect)
    aire = zr(jv_sect)

    call refe%Init(nomte)
    sigref = refe%GetRef('SIGM')
    call refe%Check()
    call tecach('OOO', 'PVECTUR', 'E', iret, nval=2, itab=itab)
    jv_vectur = itab(1)
    nddl = itab(2)

    ! pour éviter d'avoir zéro dans les directions perpendiculaires à la barre
    forref = sigref*aire
    zr(jv_vectur:jv_vectur-1+nddl) = forref

end subroutine
