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
subroutine comp_meca_rkit(factorKeyword, iFactorKeyword, relaComp, kitComp)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/ddi_kit_read.h"
#include "asterfort/getvtx.h"
#include "asterfort/lxlgut.h"
#include "asterfort/thm_kit_read.h"
!
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), intent(in) :: iFactorKeyword
    character(len=16), intent(in) :: relaComp
    character(len=16), intent(out) :: kitComp(4)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Read informations for KIT
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword to read (COMPORTEMENT)
! In  iFactorKeyword   : index of factor keyword
! In  relaComp         : behaviour (RELATION keyword)
! In  kitComp          : KIT behaviour
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbFactorKeyword
    character(len=16) :: rela_thmc, rela_hydr, rela_meca, rela_ther
    character(len=16) :: rela_flua, rela_plas, rela_cpla, rela_coup
    character(len=16) :: rela_cg(2)
    character(len=16) :: metaPhas, metaRela, relaCompMeta, metaGlob, metaPhasUser
    aster_logical :: lIsot, lCine
!
! --------------------------------------------------------------------------------------------------
!
    kitComp = 'VIDE'
!
    if (relaComp .eq. 'KIT_META') then
! ----- Get phase
        metaPhasUser = 'VIDE'
        call getvtx(factorKeyword, 'RELATION_KIT', iocc=iFactorKeyword, &
                    nbval=1, vect=metaPhasUser, nbret=nbFactorKeyword)
        ASSERT(nbFactorKeyword .eq. 1)

! ----- Using list of phases for mechanical  behaviours
        metaPhas = metaPhasUser(1:lxlgut(metaPhasUser))//"MECA"

! ----- Get behaviour
        call getvtx(factorKeyword, 'RELATION', iocc=iFactorKeyword, scal=relaCompMeta)

! ----- Internal state variables (by phase)
        metaRela = 'VIDE'
        lIsot = ASTER_FALSE
        lCine = ASTER_FALSE
        if (relaCompMeta(1:9) .eq. 'META_P_CL') then
            metaRela = 'META_P_CINE_LINE'
            lCine = ASTER_TRUE
        elseif (relaCompMeta(1:9) .eq. 'META_P_IL') then
            metaRela = 'META_P_ISOT_LINE'
            lIsot = ASTER_TRUE
        elseif (relaCompMeta(1:10) .eq. 'META_P_INL') then
            metaRela = 'META_P_ISOT_TRAC'
            lIsot = ASTER_TRUE
        elseif (relaCompMeta(1:9) .eq. 'META_V_CL') then
            metaRela = 'META_V_CINE_LINE'
            lCine = ASTER_TRUE
        elseif (relaCompMeta(1:9) .eq. 'META_V_IL') then
            metaRela = 'META_V_ISOT_LINE'
            lIsot = ASTER_TRUE
        elseif (relaCompMeta(1:10) .eq. 'META_V_INL') then
            metaRela = 'META_V_ISOT_TRAC'
            lIsot = ASTER_TRUE
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Internal state variables (global)
        metaGlob = 'VIDE'
        if (lIsot) then
            metaGlob = 'META_G_ISOT'
        elseif (lCine) then
            metaGlob = 'META_G_CINE'
        else
            ASSERT(ASTER_FALSE)
        end if
        if ((relaCompMeta(11:13) .eq. 'PT ') .or. (relaCompMeta(12:14) .eq. 'PT ')) then
            metaGlob(12:16) = '_PT  '
        end if
        if ((relaCompMeta(11:15) .eq. 'PT_RE') .or. (relaCompMeta(12:16) .eq. 'PT_RE')) then
            metaGlob(12:16) = '_PTRE'
        else if ((relaCompMeta(11:13) .eq. 'RE') .or. (relaCompMeta(12:14) .eq. 'RE')) then
            metaGlob(12:16) = '_RE  '
        end if
        kitComp(1) = metaPhas
        kitComp(2) = metaRela
        kitComp(3) = metaGlob

    else if (relaComp .eq. 'KIT_DDI') then
        call ddi_kit_read(factorKeyword, iFactorKeyword, &
                          rela_flua, rela_plas, rela_cpla, rela_coup)
        kitComp(1) = rela_flua
        kitComp(2) = rela_plas
        kitComp(3) = rela_coup
        kitComp(4) = rela_cpla

    else if (relaComp .eq. 'KIT_CG') then
        call getvtx(factorKeyword, 'RELATION_KIT', iocc=iFactorKeyword, &
                    nbval=2, vect=rela_cg, nbret=nbFactorKeyword)
        ASSERT(nbFactorKeyword .eq. 2)
        kitComp(1) = rela_cg(1)
        kitComp(2) = rela_cg(2)

    elseif ((relaComp(1:5) .eq. 'KIT_H') .or. (relaComp(1:6) .eq. 'KIT_TH')) then
        call thm_kit_read(factorKeyword, iFactorKeyword, &
                          relaComp, rela_thmc, rela_hydr, rela_meca, rela_ther)
        kitComp(1) = rela_meca
        kitComp(2) = rela_hydr
        kitComp(3) = rela_ther
        kitComp(4) = rela_thmc

    else
        ASSERT(ASTER_FALSE)
    end if
end subroutine
