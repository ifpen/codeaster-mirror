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
! ==================================================================================================
!
! Module for the management of strains
!
! ==================================================================================================
!
module BehaviourStrain_module
! ==================================================================================================
    use BehaviourStrain_type
    use MaterialPara_module
    use MaterialPara_type
    use MetallurgyMeca_module, only: metaGetType
    use calcul_module, only: ca_nbcvrc_
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: strainDetectVarc, compVarcStrain, getVarcStrain
    private :: detectVarc
    private :: compTherStrain, compTherStrainField
    private :: compSechStrainField, compHydrStrainField, compEpsaStrainField
    public :: compPtotStrainField, getStrainType
! ==================================================================================================
    private
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/get_elasth_para.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! strainDetectVarc
!
! Detect external state variables
!
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  lTHM             : flag for THM
! IO  materPara        : parameters of material
! IO  allVarcStrain    : all external state variables for anelastic strains
! In  indxVarcStrain   : index of external state variable
!
! --------------------------------------------------------------------------------------------------
    subroutine strainDetectVarc(poum, lTHM, materPara, &
                                allVarcStrain, indxVarcStrain_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: poum
        aster_logical, intent(in) :: lTHM
        type(Material_Para), intent(in) :: materPara
        type(All_Varc_Strain), intent(inout) :: allVarcStrain
        integer(kind=8), optional, intent(in) :: indxVarcStrain_
! ----- Local
        integer(kind=8) :: iVarcStrain, indxVarcStrain
!   ------------------------------------------------------------------------------------------------
!
        allVarcStrain%lTHM = lTHM
        if (ca_nbcvrc_ .ne. 0 .or. lTHM) then
            if (present(indxVarcStrain_)) then
                ASSERT(indxVarcStrain_ .gt. 0)
                indxVarcStrain = indxVarcStrain_
            else
                indxVarcStrain = VARC_STRAIN_ALL
            end if
            if (indxVarcStrain .eq. VARC_STRAIN_ALL) then
                do iVarcStrain = 1, VARC_STRAIN_NBMAXI
                    call detectVarc(poum, materPara%schemePara, &
                                    iVarcStrain, allVarcStrain%list(iVarcStrain))
                end do
            else
                ASSERT(indxVarcStrain .le. VARC_STRAIN_NBMAXI)
                call detectVarc(poum, materPara%schemePara, &
                                indxVarcStrain, allVarcStrain%list(indxVarcStrain))
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! detectVarc
!
! Detect external state variables
!
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  schemePara       : parameters of integration scheme
! In  indxVarcStrain   : index of external state variable
! IO  varcStrain       : external state variable for anelastic strains
!
! --------------------------------------------------------------------------------------------------
    subroutine detectVarc(poum, schemePara, &
                          indxVarcStrain, varcStrain)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: poum
        type(Scheme_Para), intent(in) :: schemePara
        integer(kind=8), intent(in) :: indxVarcStrain
        type(Varc_Strain), intent(inout) :: varcStrain
! ----- Local
        integer(kind=8) :: iVarcStrainCmp, strainNbCmp
        integer(kind=8) :: iretPrev, iretCurr, iretRefe
        aster_logical :: exist
        real(kind=8) :: varcPrev, varcCurr, varcRefe
        character(len=16) :: varcName
        character(len=2), parameter :: varcCmpName(6) = (/'XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ'/)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxVarcStrain .gt. 0)
        ASSERT(indxVarcStrain .le. VARC_STRAIN_NBMAXI)
        strainNbCmp = varcStrainNbcmp(indxVarcStrain)
        varcStrain%varcStrainType = indxVarcStrain
        exist = ASTER_FALSE
        do iVarcStrainCmp = 1, strainNbCmp
! --------- Prepare name of external state variable
            if (strainNbCmp .eq. 1) then
                varcName = varcStrainName(indxVarcStrain)
            else
                varcName = varcStrainName(indxVarcStrain) (1:4)//varcCmpName(iVarcStrainCmp)
            end if

! --------- Get previous and current values of external state variable
            varcPrev = r8nnem()
            varcCurr = r8nnem()
            iretPrev = 1
            if (poum .eq. '-' .or. poum .eq. 'T') then
                call rcvarc(' ', varcName, '-', &
                            schemePara%fami, schemePara%kpg, schemePara%ksp, &
                            varcPrev, iretPrev)
            end if
            iretCurr = 1
            if (poum .eq. '+' .or. poum .eq. 'T') then
                call rcvarc(' ', varcName, '+', &
                            schemePara%fami, schemePara%kpg, schemePara%ksp, &
                            varcCurr, iretCurr)
            end if
            if (poum .eq. '-') then
                if (iretPrev .eq. 0) then
                    exist = ASTER_TRUE
                end if
            end if
            if (poum .eq. '+') then
                if (iretCurr .eq. 0) then
                    exist = ASTER_TRUE
                end if
            end if
            if (poum .eq. 'T') then
                if (iretCurr+iretPrev .eq. 0) then
                    exist = ASTER_TRUE
                end if
                if (iretCurr+iretPrev .eq. 1) then
                    call utmess("F", "COMPOR7_2", sk=varcName)
                end if
            end if
            varcStrain%exist = exist
            varcStrain%varcPrev(iVarcStrainCmp) = varcPrev
            varcStrain%varcCurr(iVarcStrainCmp) = varcCurr
            varcStrain%varcIncr(iVarcStrainCmp) = varcCurr-varcPrev

! --------- Get reference value of external state variable (if exist !)
            varcRefe = r8nnem()
            iretRefe = 1
            if (varcStrainHasRefe(indxVarcStrain)) then
                call rcvarc(' ', varcName, 'REF', &
                            schemePara%fami, schemePara%kpg, schemePara%ksp, &
                            varcRefe, iretRefe)
                if (exist .and. iretRefe .eq. 1) then
                    call utmess("F", "COMPOR7_8")
                end if
            end if
            varcStrain%varcRefe = varcRefe
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compVarcStrain
!
! Compute anelastic strains from external state variables
!
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  materPara        : parameters of material
! IO  allVarcStrain    : all external state variables for anelastic strains
!
! --------------------------------------------------------------------------------------------------
    subroutine compVarcStrain(poum, materPara, allVarcStrain)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: poum
        type(Material_Para), intent(in) :: materPara
        type(All_Varc_Strain), intent(inout) :: allVarcStrain
! ----- Local
        integer(kind=8) :: indxVarcStrain
!   ------------------------------------------------------------------------------------------------
!

! ----- For temperature
        indxVarcStrain = VARC_STRAIN_TEMP
        if (allVarcStrain%list(indxVarcStrain)%exist) then
            call compTherStrainField(poum, materPara, &
                                     allVarcStrain%list(indxVarcStrain))
            allVarcStrain%hasInelasticStrains = ASTER_TRUE
        end if

! ----- For drying
        indxVarcStrain = VARC_STRAIN_SECH
        if (allVarcStrain%list(indxVarcStrain)%exist) then
            call compSechStrainField(poum, materPara, &
                                     allVarcStrain%list(indxVarcStrain))
            allVarcStrain%hasInelasticStrains = ASTER_TRUE
        end if

! ----- For hydration
        indxVarcStrain = VARC_STRAIN_HYDR
        if (allVarcStrain%list(indxVarcStrain)%exist) then
            call compHydrStrainField(poum, materPara, &
                                     allVarcStrain%list(indxVarcStrain))
            allVarcStrain%hasInelasticStrains = ASTER_TRUE
        end if

! ----- For pressure
        indxVarcStrain = VARC_STRAIN_PTOT
        if (allVarcStrain%list(indxVarcStrain)%exist) then
            call compPtotStrainField(poum, materPara, &
                                     allVarcStrain%list(indxVarcStrain))
            allVarcStrain%hasInelasticStrains = ASTER_TRUE
        end if

! ----- For any anelastic strains
        indxVarcStrain = VARC_STRAIN_EPSA
        if (allVarcStrain%list(indxVarcStrain)%exist) then
            call compEpsaStrainField(allVarcStrain%list(indxVarcStrain))
            allVarcStrain%hasInelasticStrains = ASTER_TRUE
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getVarcStrain
!
! Get anelastic strains from all external state variables
!
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  indxVarcStrain   : index of external state variable
! In  allVarcStrain    : all external state variables for anelastic strains
! Out epsiVarc         : anelastic strains from all external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine getVarcStrain(poum, indxVarcStrain, allVarcStrain, epsiVarc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: poum
        integer(kind=8), intent(in) :: indxVarcStrain
        type(All_Varc_Strain), intent(in) :: allVarcStrain
        real(kind=8), intent(out) :: epsiVarc(6)
! ----- Local
        integer(kind=8) :: iVarcStrain
!   ------------------------------------------------------------------------------------------------
!
        epsiVarc = 0.d0

        if (indxVarcStrain .eq. VARC_STRAIN_ALL) then
            do iVarcStrain = 1, VARC_STRAIN_NBMAXI
                if (allVarcStrain%list(iVarcStrain)%exist) then
                    if (poum .eq. "+") then
                        epsiVarc = epsiVarc+ &
                                   allVarcStrain%list(iVarcStrain)%fieldCurr
                    elseif (poum .eq. "-") then
                        epsiVarc = epsiVarc+ &
                                   allVarcStrain%list(iVarcStrain)%fieldPrev
                    elseif (poum .eq. "T") then
                        epsiVarc = epsiVarc+ &
                                   allVarcStrain%list(iVarcStrain)%fieldIncr
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                end if
            end do
        else
            ASSERT(indxVarcStrain .gt. 0)
            ASSERT(indxVarcStrain .le. VARC_STRAIN_NBMAXI)
            if (poum .eq. "+") then
                epsiVarc = epsiVarc+ &
                           allVarcStrain%list(indxVarcStrain)%fieldCurr
            elseif (poum .eq. "-") then
                epsiVarc = epsiVarc+ &
                           allVarcStrain%list(indxVarcStrain)%fieldPrev
            elseif (poum .eq. "T") then
                epsiVarc = epsiVarc+ &
                           allVarcStrain%list(indxVarcStrain)%fieldIncr
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compTherStrain
!
! Compute anelastic strains from thermal variable
!
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  materPara        : parameters of material
! In  tempRefe         : reference temperature
! In  temp             : temperature
! Out epsthIsot        : thermal strain (isotropic case)
! Out epsthAnis        : thermal strain (anisotropic case)
! Out epsthMeta        : thermal strain (metallurgical case)
!
! --------------------------------------------------------------------------------------------------
    subroutine compTherStrain(poum, materPara, &
                              tempRefe, temp, &
                              epsthIsot, epsthAnis, epsthMeta)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: poum
        type(Material_Para), intent(in) :: materPara
        real(kind=8), intent(in) :: tempRefe, temp
        real(kind=8), intent(out) :: epsthIsot, epsthAnis(3), epsthMeta
! ----- Local
        real(kind=8) :: alpha(2)
        real(kind=8) :: alphaL, alphaT, alphaN
        integer(kind=8) :: metaType, nbPhase
        real(kind=8) :: zSumCold, zHot
        real(kind=8) :: epsthMetaHot, epsthMetaCold
        real(kind=8) :: z_h_r, deps_ch_tref
!   ------------------------------------------------------------------------------------------------
!
        epsthIsot = 0.d0
        epsthAnis = 0.d0
        epsthMeta = 0.d0

! ----- Get elastic parameters for thermic dilatation
        call get_elasth_para(materPara%schemePara%fami, materPara%jvMaterCode, &
                             poum, materPara%schemePara%kpg, materPara%schemePara%ksp, &
                             materPara%elasID, materPara%elasKeyword, &
                             alpha=alpha, &
                             alpha_l=alphaL, &
                             alpha_t=alphaT, &
                             alpha_n=alphaN)

! ----- Get type of metallurgy and get phases
        metaType = META_NONE
        nbPhase = 0
        zSumCold = 0.d0
        zHot = 0.d0
        if (materPara%elasKeyword .eq. 'ELAS_META') then
            call metaGetType(metaType, nbPhase)
            if (nbPhase .ne. META_NONE) then
                call metaGetPhase(materPara%schemePara%fami, poum, &
                                  materPara%schemePara%kpg, materPara%schemePara%ksp, &
                                  metaType, nbPhase, &
                                  zSumCold_=zSumCold, zhot_=zHot)
            end if
        end if

! ----- Compute thermic strain
        if (materPara%elasID .eq. ELAS_ISOT) then
            if (materPara%elasKeyword .eq. 'ELAS_META') then
                if (.not. materPara%lMetaLemaAni) then
                    epsthMetaHot = zHot*alpha(1)*(temp-tempRefe)
                    epsthMetaCold = zSumCold*alpha(2)*(temp-tempRefe)
                    call get_elasth_para(materPara%schemePara%fami, materPara%jvMaterCode, &
                                         "+", materPara%schemePara%kpg, materPara%schemePara%ksp, &
                                         materPara%elasID, materPara%elasKeyword, &
                                         z_h_r_=z_h_r, deps_ch_tref_=deps_ch_tref)
                    epsthMetaHot = epsthMetaHot+(1-z_h_r)*deps_ch_tref*zHot
                    epsthMetaCold = epsthMetaCold+z_h_r*deps_ch_tref*zSumCold
                    epsthMeta = epsthMetaHot+epsthMetaCold
                end if

            else
                epsthIsot = alpha(1)*(temp-tempRefe)
            end if

        elseif (materPara%elasID .eq. ELAS_ORTH) then
            epsthAnis(1) = alphaL*(temp-tempRefe)
            epsthAnis(2) = alphaT*(temp-tempRefe)
            epsthAnis(3) = alphaN*(temp-tempRefe)

        elseif (materPara%elasID .eq. ELAS_ISTR) then
            epsthAnis(1) = alphaL*(temp-tempRefe)
            epsthAnis(2) = alphaN*(temp-tempRefe)

        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compTherStrainField
!
! Prepare external state variable TEMP for strain
!
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  materPara        : parameters of material
! IO  varcStrainTher   : external state variables parameters for temperature
!
! --------------------------------------------------------------------------------------------------
    subroutine compTherStrainField(poum, materPara, varcStrainTher)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters.
        character(len=*), intent(in) :: poum
        type(Material_Para), intent(in) :: materPara
        type(Varc_Strain), intent(inout) :: varcStrainTher
! ----- Local
        aster_logical :: lElasIsMeta
        real(kind=8) :: tempRefe, tempPrev, tempCurr
        real(kind=8) :: epsthIsotPrev, epsthAnisPrev(3), epsthMetaPrev
        real(kind=8) :: epsthIsotCurr, epsthAnisCurr(3), epsthMetaCurr
!   ------------------------------------------------------------------------------------------------
!
        varcStrainTher%fieldPrev = 0.d0
        varcStrainTher%fieldCurr = 0.d0
        varcStrainTher%fieldIncr = 0.d0
        lElasIsMeta = (materPara%elasKeyword == 'ELAS_META')

! ----- Get external state variable
        tempRefe = varcStrainTher%varcRefe
        tempPrev = varcStrainTher%varcPrev(1)
        tempCurr = varcStrainTher%varcCurr(1)

! ----- Compute thermal strains
        epsthIsotPrev = r8nnem()
        epsthAnisPrev = r8nnem()
        epsthMetaPrev = r8nnem()
        epsthIsotCurr = r8nnem()
        epsthAnisCurr = r8nnem()
        epsthMetaCurr = r8nnem()
        if (poum .eq. '-' .or. poum .eq. 'T') then
            call compTherStrain('-', materPara, &
                                tempRefe, tempPrev, &
                                epsthIsotPrev, epsthAnisPrev, epsthMetaPrev)
        end if
        if (poum .eq. '+' .or. poum .eq. 'T') then
            call compTherStrain('+', materPara, &
                                tempRefe, tempCurr, &
                                epsthIsotCurr, epsthAnisCurr, epsthMetaCurr)
        end if

! ----- Compute thermal strain fields
        if (lElasIsMeta) then
            ASSERT(materPara%elasID .eq. ELAS_ISOT)
            varcStrainTher%fieldPrev(1:3) = epsthMetaPrev
            varcStrainTher%fieldCurr(1:3) = epsthMetaCurr
            if (poum .eq. 'T') then
                varcStrainTher%fieldIncr(1:3) = epsthMetaCurr-epsthMetaPrev
            end if

        else
            if (materPara%elasID == ELAS_ISOT) then
                varcStrainTher%fieldPrev(1:3) = epsthIsotPrev
                varcStrainTher%fieldCurr(1:3) = epsthIsotCurr
                if (poum .eq. 'T') then
                    varcStrainTher%fieldIncr(1:3) = epsthIsotCurr-epsthIsotPrev
                end if

            elseif (materPara%elasID == ELAS_ORTH) then
                varcStrainTher%fieldPrev(1:3) = &
                    epsthAnisPrev(1:3)
                varcStrainTher%fieldCurr(1:3) = &
                    epsthAnisCurr(1:3)
                if (poum .eq. 'T') then
                    varcStrainTher%fieldIncr(1:3) = &
                        epsthAnisCurr(1:3)-epsthAnisPrev(1:3)
                end if

            elseif (materPara%elasID == ELAS_ISTR) then
                varcStrainTher%fieldPrev(1) = epsthAnisPrev(1)
                varcStrainTher%fieldCurr(1) = epsthAnisCurr(1)
                varcStrainTher%fieldPrev(2) = epsthAnisPrev(1)
                varcStrainTher%fieldCurr(2) = epsthAnisCurr(1)
                varcStrainTher%fieldPrev(3) = epsthAnisPrev(2)
                varcStrainTher%fieldCurr(3) = epsthAnisCurr(2)
                if (poum .eq. 'T') then
                    varcStrainTher%fieldIncr = varcStrainTher%fieldCurr-varcStrainTher%fieldPrev
                end if

            else
                ASSERT(ASTER_FALSE)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compSechStrainField
!
! Prepare external state variable SECH for strain
!
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  materPara        : parameters of material
! IO  varcStrainSech   : external state variables parameters for drying
!
! --------------------------------------------------------------------------------------------------
    subroutine compSechStrainField(poum, materPara, varcStrainSech)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters.
        character(len=*), intent(in) :: poum
        type(Material_Para), intent(in) :: materPara
        type(Varc_Strain), intent(inout) :: varcStrainSech
! ----- Local
        integer(kind=8), parameter :: nbProp = 1
        real(kind=8) :: propVale(nbProp)
        character(len=16), parameter :: propName(nbProp) = ('K_DESSIC')
        integer(kind=8) :: codret(nbProp)
        real(kind=8) :: sechRefe, sechPrev, sechCurr
        real(kind=8) :: kdessm, kdessp
!   ------------------------------------------------------------------------------------------------
!
        varcStrainSech%fieldPrev = 0.d0
        varcStrainSech%fieldCurr = 0.d0
        varcStrainSech%fieldIncr = 0.d0

! ----- Get external state variable
        sechRefe = varcStrainSech%varcRefe
        sechPrev = varcStrainSech%varcPrev(1)
        sechCurr = varcStrainSech%varcCurr(1)

! ----- Get elastic parameters
        kdessm = r8nnem()
        kdessp = r8nnem()
        if (poum .eq. '-' .or. poum .eq. 'T') then
            call rcvalb(materPara%schemePara%fami, &
                        materPara%schemePara%kpg, materPara%schemePara%ksp, &
                        '-', materPara%jvMaterCode, ' ', materPara%elasKeyword, &
                        0, ' ', [0.d0], &
                        nbProp, propName, propVale, &
                        codret, 1)
            kdessm = propVale(1)
        end if
        if (poum .eq. '+' .or. poum .eq. 'T') then
            call rcvalb(materPara%schemePara%fami, &
                        materPara%schemePara%kpg, materPara%schemePara%ksp, &
                        '+', materPara%jvMaterCode, ' ', materPara%elasKeyword, &
                        0, ' ', [0.d0], &
                        nbProp, propName, propVale, &
                        codret, 1)
            kdessp = propVale(1)
        end if

! ----- Compute strain fields
        varcStrainSech%fieldPrev(1:3) = -kdessm*(sechRefe-sechPrev)
        varcStrainSech%fieldCurr(1:3) = -kdessp*(sechRefe-sechCurr)
        if (poum .eq. 'T') then
            varcStrainSech%fieldIncr = varcStrainSech%fieldCurr- &
                                       varcStrainSech%fieldPrev
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compHydrStrainField
!
! Prepare external state variable HYDR for strain
!
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  materPara        : parameters of material
! IO  varcStrainHydr   : external state variables parameters for hydratation
!
! --------------------------------------------------------------------------------------------------
    subroutine compHydrStrainField(poum, materPara, varcStrainHydr)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters.
        character(len=*), intent(in) :: poum
        type(Material_Para), intent(in) :: materPara
        type(Varc_Strain), intent(inout) :: varcStrainHydr
! ----- Local
        integer(kind=8), parameter :: nbProp = 1
        real(kind=8) :: propVale(nbProp)
        character(len=16), parameter :: propName(nbProp) = ('B_ENDOGE')
        integer(kind=8) :: codret(nbProp)
        real(kind=8) :: hydrPrev, hydrCurr
        real(kind=8) :: bendom, bendop
!   ------------------------------------------------------------------------------------------------
!
        varcStrainHydr%fieldPrev = 0.d0
        varcStrainHydr%fieldCurr = 0.d0
        varcStrainHydr%fieldIncr = 0.d0

! ----- Get external state variable
        hydrPrev = varcStrainHydr%varcPrev(1)
        hydrCurr = varcStrainHydr%varcCurr(1)

! ----- Get elastic parameters
        bendom = r8nnem()
        bendop = r8nnem()
        if (poum .eq. '-' .or. poum .eq. 'T') then
            call rcvalb(materPara%schemePara%fami, &
                        materPara%schemePara%kpg, materPara%schemePara%ksp, &
                        '-', materPara%jvMaterCode, ' ', materPara%elasKeyword, &
                        0, ' ', [0.d0], &
                        nbProp, propName, propVale, &
                        codret, 1)
            bendom = propVale(1)
        end if
        if (poum .eq. '+' .or. poum .eq. 'T') then
            call rcvalb(materPara%schemePara%fami, &
                        materPara%schemePara%kpg, materPara%schemePara%ksp, &
                        '+', materPara%jvMaterCode, ' ', materPara%elasKeyword, &
                        0, ' ', [0.d0], &
                        nbProp, propName, propVale, &
                        codret, 1)
            bendop = propVale(1)
        end if

! ----- Compute strain fields
        varcStrainHydr%fieldPrev(1:3) = -bendom*hydrPrev
        varcStrainHydr%fieldCurr(1:3) = -bendop*hydrCurr
        if (poum .eq. 'T') then
            varcStrainHydr%fieldIncr = varcStrainHydr%fieldCurr-varcStrainHydr%fieldPrev
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compPtotStrainField
!
! Prepare external state variable PTOT for strain
!
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  materPara        : parameters of material
! IO  varcStrainPtot   : external state variables parameters for PTOT
!
! --------------------------------------------------------------------------------------------------
    subroutine compPtotStrainField(poum, materPara, varcStrainPtot)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters.
        character(len=*), intent(in) :: poum
        type(Material_Para), intent(in) :: materPara
        type(Varc_Strain), intent(inout) :: varcStrainPtot
! ----- Local
        integer(kind=8), parameter :: nbParaBiot = 1
        integer(kind=8)  :: codretBiot(nbParaBiot)
        real(kind=8) :: paraValeBiot(nbParaBiot)
        character(len=16), parameter :: paraNameBiot(nbParaBiot) = (/'BIOT_COEF'/)
        real(kind=8) :: ptotPrev, ptotCurr
        real(kind=8) :: biotp, biotm
        real(kind=8) :: troikp, troikm, nup, num, ep, em
!   ------------------------------------------------------------------------------------------------
!
        varcStrainPtot%fieldPrev = 0.d0
        varcStrainPtot%fieldCurr = 0.d0
        varcStrainPtot%fieldIncr = 0.d0

! ----- Get external state variable
        ptotPrev = varcStrainPtot%varcPrev(1)
        ptotCurr = varcStrainPtot%varcCurr(1)

! ----- Get Biot coefficients
        biotm = r8nnem()
        biotp = r8nnem()
        if (poum .eq. '-' .or. poum .eq. 'T') then
            call rcvalb(materPara%schemePara%fami, &
                        materPara%schemePara%kpg, materPara%schemePara%ksp, &
                        '-', materPara%jvMaterCode, ' ', 'THM_DIFFU', &
                        0, ' ', [0.d0], &
                        nbParaBiot, paraNameBiot, paraValeBiot, &
                        codretBiot, 1)
            if (codretBiot(1) .ne. 0) then
                paraValeBiot(1) = 0.d0
            end if
            biotm = paraValeBiot(1)
        end if
        if (poum .eq. '+' .or. poum .eq. 'T') then
            call rcvalb(materPara%schemePara%fami, &
                        materPara%schemePara%kpg, materPara%schemePara%ksp, &
                        '+', materPara%jvMaterCode, ' ', 'THM_DIFFU', &
                        0, ' ', [0.d0], &
                        nbParaBiot, paraNameBiot, paraValeBiot, &
                        codretBiot, 1)
            if (codretBiot(1) .ne. 0) then
                paraValeBiot(1) = 0.d0
            end if
            biotp = paraValeBiot(1)
        end if

! ----- Get elastic parameters
        if (materPara%elasID .ne. ELAS_ISOT) then
            call utmess("F", "COMPOR7_1")
        end if
        troikm = r8nnem()
        troikp = r8nnem()
        if (poum .eq. '-' .or. poum .eq. 'T') then
            call get_elas_para(materPara%schemePara%fami, &
                               materPara%jvMaterCode, '-', &
                               materPara%schemePara%kpg, materPara%schemePara%ksp, &
                               materPara%elasID, materPara%elasKeyword, &
                               e_=em, nu_=num)
            troikm = em/(1.d0-2.d0*num)
        end if
        if (poum .eq. '+' .or. poum .eq. 'T') then
            call get_elas_para(materPara%schemePara%fami, &
                               materPara%jvMaterCode, '+', &
                               materPara%schemePara%kpg, materPara%schemePara%ksp, &
                               materPara%elasID, materPara%elasKeyword, &
                               e_=ep, nu_=nup)
            troikp = ep/(1.d0-2.d0*nup)
        end if

! ----- Compute strain fields
        varcStrainPtot%fieldPrev(1:3) = (biotm/troikm)*ptotPrev
        varcStrainPtot%fieldCurr(1:3) = (biotp/troikp)*ptotCurr
        if (poum .eq. 'T') then
            varcStrainPtot%fieldIncr = varcStrainPtot%fieldCurr-varcStrainPtot%fieldPrev
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpsaStrainField
!
! Prepare external state variable EPSA for strain
!
! IO  varcStrainEpsa   : external state variables parameters for EPSA
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpsaStrainField(varcStrainEpsa)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters.
        type(Varc_Strain), intent(inout) :: varcStrainEpsa
! ----- Local
        integer(kind=8), parameter :: indxVarcStrain = VARC_STRAIN_EPSA
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
        real(kind=8) :: epsaPrev, epsaCurr, epsaIncr
        integer(kind=8) :: strainNbCmp, iVarcStrainCmp
!   ------------------------------------------------------------------------------------------------
!
        varcStrainEpsa%fieldPrev = 0.d0
        varcStrainEpsa%fieldCurr = 0.d0
        varcStrainEpsa%fieldIncr = 0.d0

! ----- Get external state variables
        strainNbCmp = varcStrainNbcmp(indxVarcStrain)
        do iVarcStrainCmp = 1, strainNbCmp
            epsaPrev = varcStrainEpsa%varcPrev(iVarcStrainCmp)
            epsaCurr = varcStrainEpsa%varcCurr(iVarcStrainCmp)
            epsaIncr = varcStrainEpsa%varcIncr(iVarcStrainCmp)
            if (iVarcStrainCmp .ge. 4) then
                epsaPrev = epsaPrev*rac2
                epsaCurr = epsaCurr*rac2
                epsaIncr = epsaIncr*rac2
            end if
            varcStrainEpsa%fieldPrev(iVarcStrainCmp) = epsaPrev
            varcStrainEpsa%fieldCurr(iVarcStrainCmp) = epsaCurr
            varcStrainEpsa%fieldIncr(iVarcStrainCmp) = epsaIncr
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine

! --------------------------------------------------------------------------------------------------
!
! getStrainType
!
! Get type of strain from option
!
! In  option           : option to compute
! Out strainType       : type of strain (small, Green, log, etc.)
! Out lStrainMeca      : flag to compute mechanical strains
!
! --------------------------------------------------------------------------------------------------
    subroutine getStrainType(optionZ, strainType, lStrainMeca)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters.
        character(len=*), intent(in) :: optionZ
        integer(kind=8), intent(out) :: strainType
        aster_logical, intent(out) :: lStrainMeca
! ----- Local
        character(len=16) :: option
!   ------------------------------------------------------------------------------------------------
!
        option = optionZ
        strainType = STRAIN_TYPE_NONE
        lStrainMeca = option(1:3) .eq. "EPM"
        if (option(4:4) .eq. 'I') then
            strainType = STRAIN_TYPE_SMALL
        elseif (option(4:4) .eq. 'E') then
            strainType = STRAIN_TYPE_SMALL
        elseif (option(4:4) .eq. 'L') then
            strainType = STRAIN_TYPE_LOG
        elseif (option(4:4) .eq. 'G') then
            strainType = STRAIN_TYPE_GREEN
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module BehaviourStrain_module
