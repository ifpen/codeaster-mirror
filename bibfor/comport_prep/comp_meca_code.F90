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
subroutine comp_meca_code(relaComp, defoCompIn, typeCpla, kitComp, &
                          postIter, reguVisc, postIncr, &
                          compCodePY)
!
    use MetallurgyMeca_module
    implicit none
!
#include "asterc/lccree.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
!
    character(len=16), intent(in) :: relaComp, defoCompIn, typeCpla, kitComp(4)
    character(len=16), intent(in) :: postIter, reguVisc, postIncr
    character(len=16), intent(out) :: compCodePY
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Coding composite behaviour
!
! --------------------------------------------------------------------------------------------------
!
! In  relaComp         : behaviour (RELATION keyword)
! In  defoComp         : model of strain (DEFORMATION keyword)
! In  typeCpla         : plane stress method (analytical or De Borst algorithm)
! In  kitComp          : KIT behaviour
! In  postIter         : type of post_treatment at each Newton iteration (POST_ITER keyword)
! In  reguVisc         : keyword for viscuous regularization (REGU_VISC keyword)
! In  postIncr         : type of post-treatment at end of time step (POST_INCR keyword)
! Out compCodePY       : composite coded comportment (coding in Python)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: NoKitComp(4) = (/'VIDE', 'VIDE', 'VIDE', 'VIDE'/)
    integer(kind=8) :: nbCompElem, iKit
    character(len=16) :: compElem(20), postIncrCode
    aster_logical :: lHardIsot, lHardKine, lHardMixed
    character(len=16) :: defoComp
!
! --------------------------------------------------------------------------------------------------
!
    nbCompElem = 0
    compElem = 'VIDE'

! - Select strain model
    defoComp = defoCompIn
    if (relaComp .eq. 'KIT_META') then
        if (defoCompIn .eq. "SIMO_MIEHE") then
            defoComp = "SIMO_MIEHE_META"
        end if
    end if

! - Create composite behaviour
    nbCompElem = nbCompElem+1
    compElem(nbCompElem) = relaComp
    if (relaComp .eq. 'KIT_META') then
        do iKit = 1, 4
            nbCompElem = nbCompElem+1
            compElem(nbCompElem) = NoKitComp(iKit)
        end do
    else
        do iKit = 1, 4
            nbCompElem = nbCompElem+1
            compElem(nbCompElem) = kitComp(iKit)
        end do
    end if
    nbCompElem = nbCompElem+1
    compElem(nbCompElem) = reguVisc
    if (postIter .ne. ' ') then
        nbCompElem = nbCompElem+1
        compElem(nbCompElem) = postIter
    end if
    if (postIncr .ne. ' ') then
        postIncrCode = postIncr
        if (postIncrCode .eq. "REST_ECRO") then
            call metaAnnealGetType(relaComp, lHardIsot, lHardKine, lHardMixed)
            if (lHardIsot) then
                postIncrCode = "REST_ECRO_ISOT"
            elseif (lHardKine) then
                postIncrCode = "REST_ECRO_CINE"
            elseif (lHardMixed) then
                postIncrCode = "REST_ECRO_ECMI"
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        nbCompElem = nbCompElem+1
        compElem(nbCompElem) = postIncrCode
    end if
    nbCompElem = nbCompElem+1
    compElem(nbCompElem) = defoComp
    nbCompElem = nbCompElem+1
    compElem(nbCompElem) = typeCpla

! - Coding composite comportment (Python)
    call lccree(nbCompElem, compElem, compCodePY)
!
end subroutine
