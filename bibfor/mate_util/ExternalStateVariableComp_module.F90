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
! aslint: disable=W1501
!
! ==================================================================================================
!
! Module to compute External State Variables
!
! ==================================================================================================
!
module ExternalStateVariableComp_module
! ==================================================================================================
! ==================================================================================================
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    use BehaviourStrain_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: varcCompElem
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/gcnco2.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/nmvcd2.h"
#include "asterfort/nmvcex.h"
#include "asterfort/reajre.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
#include "asterfort/vemare.h"
#include "asterfort/xajcin.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! varcCompElem
!
! Compute elementary vector for external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine varcCompElem(lSuperElem, poum, &
                            numeHarm, modelZ, caraElemZ, materFieldZ, materCodeZ, &
                            varcRefeZ, varcPrevZ, &
                            jvBase, vectElemZ, &
                            lCumul_, varcCurrZ_, variZ_, sigmZ_, comporZ_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: lSuperElem
        character(len=1), intent(in) :: poum
        integer(kind=8), intent(in) :: numeHarm
        character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ, materCodeZ
        character(len=*), intent(in) :: varcRefeZ, varcPrevZ
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(in) :: vectElemZ
        aster_logical, optional, intent(in) :: lCumul_
        character(len=*), optional, intent(in) :: varcCurrZ_, variZ_, sigmZ_, comporZ_
! ----- Locals
        character(len=24), parameter :: chsith = '&&NMVCPR.CHSITH'
        integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOutMax = 2
        character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
        character(len=24) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
        character(len=8) :: newnom
        character(len=16) :: option
        character(len=24) :: modelLigrel, resuElem
        character(len=24) :: chgeom, chharm, varcCurr, varcPrev, varcRefe
        character(len=24) :: vari, sigm, compor
        character(len=24) :: varcAllRefe, varcAllPrev, varcAllCurr
        character(len=24) :: timePrev, timeCurr
        aster_logical :: lXFEM, lCumul, lExistVarc
        aster_logical :: lExistTemp, lExistMetaZirc, lExistMetaAcier, lCalcMeta
        integer(kind=8) :: nbFieldIn, nbFieldOut, iret, iVarcStrain, iFieldCompor
!   ------------------------------------------------------------------------------------------------
!
        lCumul = ASTER_FALSE
        if (present(lCumul_)) then
            lCumul = lCumul_
        end if
        varcCurr = " "
        if (present(varcCurrZ_)) then
            varcCurr = varcCurrZ_
        end if
        varcPrev = varcPrevZ
        varcRefe = varcRefeZ
        vari = " "
        sigm = " "
        compor = " "
        if (present(variZ_)) then
            vari = variZ_
        end if
        if (present(sigmZ_)) then
            sigm = sigmZ_
        end if
        if (present(comporZ_)) then
            compor = comporZ_
        end if

! ----- Initializations
        call exixfe(modelZ, iret)
        lXFEM = (iret .eq. 1)
        call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)
        lpain = " "
        lpaout = " "
        lchin = " "
        lchout = " "
        newnom = '.VEXXXXX'
        resuElem = vectElemZ(1:8)//'.0000000'
        iFieldCompor = 0

! ----- Get geometry field
        call megeom(modelZ, chgeom)

! ----- Create field for Fourier
        call meharm(modelZ, numeHarm, chharm)

! ----- Command variables affected
        call nmvcd2('TEMP', materFieldZ, lExistTemp)
        call nmvcd2('M_ZIRC', materFieldZ, lExistMetaZirc)
        call nmvcd2('M_ACIER', materFieldZ, lExistMetaAcier)
        lCalcMeta = lExistTemp .and. (lExistMetaZirc .or. lExistMetaAcier) .and. poum .eq. "+"

! ----- Get external state variables
        call nmvcex('TOUT', varcRefe, varcAllRefe)
        call nmvcex('TOUT', varcPrev, varcAllPrev)
        varcAllCurr = " "
        if (poum .eq. "+") then
            call nmvcex('TOUT', varcCurr, varcAllCurr)
        end if

! ----- Get field for time
        call nmvcex('INST', varcPrev, timePrev)
        timeCurr = " "
        if (poum .eq. "+") then
            call nmvcex('INST', varcCurr, timeCurr)
        end if

! ----- Add input fields
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PMATERC'
        lchin(2) = materCodeZ
        lpain(3) = 'PHARMON'
        lchin(3) = chharm
        if (lSuperElem) then
            lpain(4) = 'PVARCRR'
            lchin(4) = varcRefe
        else
            lpain(4) = 'PVARCRR'
            lchin(4) = varcAllRefe
        end if
        lpain(5) = 'PCONTMR'
        lchin(5) = sigm
        lpain(6) = 'PVARIPR'
        lchin(6) = vari
        lpain(7) = 'PCOMPOR'
        lchin(7) = compor
        nbFieldIn = 7

! ----- Add specific input fields
        if (poum .eq. '-') then
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PINSTR'
            lchin(nbFieldIn) = timePrev
            nbFieldIn = nbFieldIn+1
            if (lSuperElem) then
                lpain(nbFieldIn) = 'PVARCPR'
                lchin(nbFieldIn) = varcPrev
            else
                lpain(nbFieldIn) = 'PVARCPR'
                lchin(nbFieldIn) = varcAllPrev
            end if
        else if (poum .eq. '+') then
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PINSTR'
            lchin(nbFieldIn) = timeCurr
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PVARCPR'
            lchin(nbFieldIn) = varcAllCurr
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PVARCMR'
            lchin(nbFieldIn) = varcAllPrev
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Add fields for structural elements
        call setStructFields(caraElemZ, nbFieldInMax, lchin, lpain, nbFieldIn)

! ----- Add fields for orientation
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElemZ)

! ----- Add XFEM input fields
        if (lXFEM) then
            call xajcin(modelZ, 'CHAR_MECA_TEMP_R', nbFieldInMax, lchin, lpain, nbFieldIn)
        end if

! ----- Add HHO fields
        call hhoAddInputField(modelZ, nbFieldInMax, lchin, lpain, nbFieldIn)

! ----- Set output fields
        lpaout(1) = 'PVECTUR'
        lpaout(2) = 'PCONTRT'

! ----- Allocate result
        if (.not. lCumul) then
            call detrsd('VECT_ELEM', vectElemZ)
            call vemare(jvBase, vectElemZ, modelZ)
            call reajre(vectElemZ, ' ', jvBase)
        end if

        do iVarcStrain = 1, VARC_STRAIN_NBMAXI
! --------- Detection
            call nmvcd2(varcStrainName(iVarcStrain), materFieldZ, lExistVarc)

            if (lExistVarc) then
! ------------- Current external state variable
                if (lSuperElem) then
                    if (.not. lExistTemp) then
                        call utmess('F', 'SOUSTRUC_18')
                    end if
                end if

! ------------- Generate name of RESU_ELEM
                call codent(iVarcStrain, 'D0', newnom(4:6))
                resuElem = vectElemZ(1:8)//newnom

! ------------- Set output and option
                call corich('E', resuElem, ichin_=-1)
                lchout(1) = resuElem
                option = varcRHSOption(iVarcStrain)
                nbFieldOut = 1
                if (lXFEM .and. lExistTemp) then
                    call detrsd('CHAM_ELEM', chsith)
                    call alchml(modelLigrel, 'SIEF_ELGA', 'PCONTRR', 'V', chsith, iret, ' ')
                    lchout(2) = chsith(1:19)
                    nbFieldOut = nbFieldOut+1
                end if

! ------------- Compute
                call calcul('C', option, modelLigrel, &
                            nbFieldIn, lchin, lpain, &
                            nbFieldOut, lchout, lpaout, &
                            jvBase, 'OUI')
                call reajre(vectElemZ, lchout(1), jvBase)
            end if

        end do

! ----- Special for mettalurgy
        if (lCalcMeta) then
! --------- Current external state variable
            if (lSuperElem) then
                call utmess('F', 'SOUSTRUC_18')
            end if
! --------- Generate name of RESU_ELEM
            call codent(VARC_STRAIN_NBMAXI+1, 'D0', newnom(4:6))
            resuElem = vectElemZ(1:8)//newnom

! --------- Set output and option
            call corich('E', resuElem, ichin_=-1)
            lchout(1) = resuElem
            option = 'CHAR_MECA_META_Z'
            nbFieldOut = 1

! --------- Compute
            call calcul('S', option, modelLigrel, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        jvBase, 'OUI')
            call reajre(vectElemZ, lchout(1), jvBase)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module ExternalStateVariableComp_module
