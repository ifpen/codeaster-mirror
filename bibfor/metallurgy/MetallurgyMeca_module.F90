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
! aslint: disable=W1306
!
module MetallurgyMeca_module
! ==================================================================================================
    use calcul_module, only: calcul_status
    use Metallurgy_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: metaAnnealGetType, metaGetType, metaGetParaElas, metaGetParaMixture
    public :: metaGetParaAnneal, metaGetParaVisc, metaGetParaPlasTransf, metaGetParaHardLine
    public :: metaGetParaHardTrac, metaGetMechanism
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/rcfonc.h"
#include "asterfort/rctrac.h"
#include "asterfort/rctype.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! metaAnnealGetType
!
! Get type of hardening to apply annealing
!
! --------------------------------------------------------------------------------------------------
    subroutine metaAnnealGetType(relaComp, lHardIsot, lHardKine, lHardMixed, &
                                 nbVariAnneal_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: relaComp
        aster_logical, intent(out) :: lHardIsot, lHardKine, lHardMixed
        integer(kind=8), optional, intent(out) :: nbVariAnneal_
! ----- Local
        integer(kind=8) :: nbVariAnneal
!   ------------------------------------------------------------------------------------------------
!
        lHardIsot = ASTER_FALSE
        lHardKine = ASTER_FALSE
        lHardMixed = ASTER_FALSE
        nbVariAnneal = 0
        if (relaComp .eq. 'VMIS_ISOT_LINE') then
            lHardIsot = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_ISOT_TRAC') then
            lHardIsot = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_CINE_LINE') then
            lHardKine = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_ECMI_LINE') then
            lHardMixed = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_CIN1_CHAB') then
            lHardIsot = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_CIN2_CHAB') then
            lHardIsot = ASTER_TRUE
        elseif (relaComp .eq. 'VMIS_ISOT_NL') then
            lHardIsot = ASTER_TRUE
        else
            ASSERT(ASTER_FALSE)
        end if
        if (lHardIsot) then
            nbVariAnneal = 1
        end if
        if (lHardKine) then
            nbVariAnneal = 6
        end if
        if (lHardMixed) then
            nbVariAnneal = 7
        end if
        if (present(nbVariAnneal_)) then
            nbVariAnneal_ = nbVariAnneal
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaGetType
!
! Get type of metallurgy
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetType(metaType, nbPhase, metaTypeName_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(out) :: metaType, nbPhase
        character(len=16), optional, intent(out) :: metaTypeName_
! ----- Local
        integer(kind=8), parameter :: kpg = 1, ksp = 1
        character(len=8), parameter :: steel = "BAINITE", zirc = "ALPHPUR", steelTemper = "BAINITB"
        character(len=8) :: fami
        character(len=16) :: metaTypeName
        integer(kind=8) :: iret_steel, iret_zirc, iret_steelTemper
        real(kind=8) :: r8dummy
! --------------------------------------------------------------------------------------------------
!
        metaType = META_NONE
        metaTypeName = " "
        nbPhase = 0

! ----- Choice of integration scheme: for CALC_POINT_MAT is PMAT !
        if (calcul_status() .eq. 2) then
            fami = 'PMAT'
        else
            fami = 'RIGI'
        end if

! ----- Detect phases
        call rcvarc(' ', steel, '+', fami, kpg, &
                    ksp, r8dummy, iret_steel)
        if (iret_steel .eq. 0) then
            metaType = META_STEEL
            nbPhase = PSTEEL_NB
            metaTypeName = "ACIER"
        else
            call rcvarc(' ', zirc, '+', fami, kpg, &
                        ksp, r8dummy, iret_zirc)
            if (iret_zirc .eq. 0) then
                metaType = META_ZIRC
                nbPhase = PZIRC_NB
                metaTypeName = "ZIRC"
            else
                call rcvarc(' ', steelTemper, '+', fami, kpg, &
                            ksp, r8dummy, iret_steelTemper)
                if (iret_steelTemper .eq. 0) then
                    metaType = META_STEEL_TEMPER
                    nbPhase = PRSTEEL_NB
                    metaTypeName = "ACIER_REVENU"
                else
                    metaType = META_NONE
                    nbPhase = 0
                end if
            end if
        end if

        if (present(metaTypeName_)) then
            metaTypeName_ = metaTypeName
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaGetMechanism
!
! Get mechanism of behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetMechanism(metaRela, metaGlob, &
                                lPlas_, lVisc_, &
                                lHardIsotLine_, lHardIsotNLin_, &
                                lHardKine_, lHardLine_, lAnneal_, &
                                lPlasTran_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: metaRela, metaGlob
        aster_logical, optional, intent(out) :: lPlas_
        aster_logical, optional, intent(out) :: lVisc_
        aster_logical, optional, intent(out) :: lHardIsotLine_, lHardIsotNLin_
        aster_logical, optional, intent(out) :: lHardKine_
        aster_logical, optional, intent(out) :: lHardLine_
        aster_logical, optional, intent(out) :: lAnneal_
        aster_logical, optional, intent(out) :: lPlasTran_
! --------------------------------------------------------------------------------------------------
!
        if (present(lPlas_)) then
            lPlas_ = ASTER_FALSE
            if (metaRela(6:6) .eq. 'P') then
                lPlas_ = ASTER_TRUE
            end if
        end if
!
        if (present(lVisc_)) then
            lVisc_ = ASTER_FALSE
            if (metaRela(6:6) .eq. 'V') then
                lVisc_ = ASTER_TRUE
            end if
        end if
!
        if (present(lAnneal_)) then
            lAnneal_ = ASTER_FALSE
            if (metaGlob(12:16) .eq. '_RE  ' .or. metaGlob(12:16) .eq. '_PTRE') then
                lAnneal_ = ASTER_TRUE
            end if
        end if
!
        if (present(lPlasTran_)) then
            lPlasTran_ = ASTER_FALSE
            if (metaGlob(12:16) .eq. '_PT  ' .or. metaGlob(12:16) .eq. '_PTRE') then
                lPlasTran_ = ASTER_TRUE
            end if
        end if
!
        if (present(lHardIsotLine_)) then
            lHardIsotLine_ = ASTER_FALSE
            if (metaRela(8:16) .eq. 'ISOT_LINE') then
                lHardIsotLine_ = ASTER_TRUE
            end if
        end if
!
        if (present(lHardIsotNLin_)) then
            lHardIsotNLin_ = ASTER_FALSE
            if (metaRela(8:16) .eq. 'ISOT_TRAC') then
                lHardIsotNLin_ = ASTER_TRUE
            end if
        end if
!
        if (present(lHardKine_)) then
            lHardKine_ = ASTER_FALSE
            if (metaRela(8:16) .eq. 'CINE_LINE') then
                lHardKine_ = ASTER_TRUE
            end if
        end if
!
        if (present(lHardLine_)) then
            lHardLine_ = ASTER_FALSE
            if (metaRela(8:16) .eq. 'ISOT_LINE' .or. metaRela(8:16) .eq. 'CINE_LINE') then
                lHardLine_ = ASTER_TRUE
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaGetParaElas
!
! Get elasticity parameters
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetParaElas(poum, fami, kpg, ksp, jvMaterCode, &
                               e_, deuxmu_, mu_, troisk_, &
                               em_, deuxmum_, mum_, troiskm_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: poum
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp
        integer(kind=8), intent(in) :: jvMaterCode
        real(kind=8), optional, intent(out) :: e_, deuxmu_, mu_, troisk_
        real(kind=8), optional, intent(out) :: em_, deuxmum_, mum_, troiskm_
! ----- Local
        integer(kind=8), parameter :: nbProp = 2
        real(kind=8) :: propVale(nbProp)
        integer(kind=8) :: propCode(nbProp)
        character(len=16), parameter :: propName(nbProp) = (/'E ', 'NU'/)
        real(kind=8) :: em, num, deuxmum, mum, troiskm
        real(kind=8) :: e, nu, deuxmu, mu, troisk
! --------------------------------------------------------------------------------------------------
!
        call rcvalb(fami, kpg, ksp, '-', jvMaterCode, &
                    ' ', 'ELAS_META', 0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 2)
        em = propVale(1)
        num = propVale(2)
        deuxmum = em/(1.d0+num)
        mum = deuxmum/2.d0
        troiskm = em/(1.d0-2.d0*num)
        call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                    ' ', 'ELAS_META', 0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 2)
        e = propVale(1)
        nu = propVale(2)
        deuxmu = e/(1.d0+nu)
        mu = deuxmu/2.d0
        troisk = e/(1.d0-2.d0*nu)
!
        if (present(e_)) then
            e_ = e
        end if
        if (present(deuxmu_)) then
            deuxmu_ = deuxmu
        end if
        if (present(mu_)) then
            mu_ = mu
        end if
        if (present(troisk_)) then
            troisk_ = troisk
        end if
!
        if (present(em_)) then
            em_ = em
        end if
        if (present(deuxmum_)) then
            deuxmum_ = deuxmum
        end if
        if (present(mum_)) then
            mum_ = mum
        end if
        if (present(troiskm_)) then
            troiskm_ = troiskm
        end if
!
! --------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaGetParaMixture
!
! Get parameters for mixing law
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetParaMixture(poum, fami, kpg, ksp, jvMaterCode, &
                                  l_visc, metaType, nbPhase, zColdCurr, &
                                  fmix, sy_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: poum
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode
        integer(kind=8), intent(in) :: metaType
        integer(kind=8), intent(in) :: nbPhase
        aster_logical, intent(in) :: l_visc
        real(kind=8), intent(in) :: zColdCurr
        real(kind=8), intent(out) :: fmix
        real(kind=8), optional, intent(out) :: sy_(nbPhase)
! ----- Local
        real(kind=8) :: propVale(nbPhase)
        integer(kind=8) :: propCode(nbPhase)
        character(len=16) :: propName(nbPhase)
        integer(kind=8) :: nbProp, iProp
! --------------------------------------------------------------------------------------------------
!
        ASSERT(nbPhase .le. META_MECA_NBPHASE_MAXI)
        if (present(sy_)) then
            sy_ = 0.d0
        end if

! ----- Get mixing function
        propName = " "
        nbProp = 1
        propName(1) = 'SY_MELANGE'
        if (l_visc) then
            propName(1) = 'S_VP_MELANGE'
        end if
        fmix = 0.d0
        call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                    ' ', 'ELAS_META', 1, 'META', [zColdCurr], &
                    nbProp, propName, propVale, propCode, 0)
        if (propCode(1) .eq. 0) then
            fmix = propVale(1)
        else
            fmix = zColdCurr
        end if

! ----- Get elasticity yield by phases
        propName = " "
        if (present(sy_)) then
            nbProp = nbPhase
            if (metaType .eq. META_STEEL) then
                propName(1) = 'F1_SY'
                propName(2) = 'F2_SY'
                propName(3) = 'F3_SY'
                propName(4) = 'F4_SY'
                propName(5) = 'C_SY'
                if (l_visc) then
                    propName(1) = 'F1_S_VP'
                    propName(2) = 'F2_S_VP'
                    propName(3) = 'F3_S_VP'
                    propName(4) = 'F4_S_VP'
                    propName(5) = 'C_S_VP'
                end if
            elseif (metaType .eq. META_ZIRC) then
                propName(1) = 'F1_SY'
                propName(2) = 'F2_SY'
                propName(3) = 'C_SY'
                if (l_visc) then
                    propName(1) = 'F1_S_VP'
                    propName(2) = 'F2_S_VP'
                    propName(3) = 'C_S_VP'
                end if
            elseif (metaType .eq. META_STEEL_TEMPER) then
                propName(1) = 'F1_SY'
                propName(2) = 'F2_SY'
                propName(3) = 'F3_SY'
                propName(4) = 'F4_SY'
                propName(5) = 'F5_SY'
                propName(6) = 'F6_SY'
                propName(7) = 'C_SY'
                if (l_visc) then
                    propName(1) = 'F1_S_VP'
                    propName(2) = 'F2_S_VP'
                    propName(3) = 'F3_S_VP'
                    propName(4) = 'F4_S_VP'
                    propName(5) = 'F5_S_VP'
                    propName(6) = 'F6_S_VP'
                    propName(7) = 'C_S_VP'
                end if
            else
                ASSERT(ASTER_FALSE)
            end if

            call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                        ' ', 'ELAS_META', 0, ' ', [0.d0], &
                        nbProp, propName, propVale, propCode, 2)
            do iProp = 1, nbProp
                sy_(iProp) = propVale(iProp)
            end do
        end if
!
! --------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaGetParaAnneal
!
! Get parameters for annealing
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetParaAnneal(poum, fami, kpg, ksp, jvMaterCode, &
                                 metaType, nbPhaseCold, &
                                 annealTheta)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: poum
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode
        integer(kind=8), intent(in) :: metaType, nbPhaseCold
        real(kind=8), intent(out) :: annealTheta(2*nbPhaseCold)
! ----- Local
        integer(kind=8), parameter :: nbPropMaxi = 2*(META_MECA_NBPHASE_MAXI-1)
        real(kind=8) :: propVale(nbPropMaxi)
        integer(kind=8) :: propCode(nbPropMaxi)
        character(len=16) :: propName(nbPropMaxi)
        integer(kind=8) :: nbProp, iProp
! --------------------------------------------------------------------------------------------------
!
        nbProp = 2*(nbPhaseCold-1)
        annealTheta = 0.d0
        ASSERT(nbPhaseCold .le. META_MECA_NBPHASE_MAXI)

! ----- Name of parameters
        propName = " "
        if (metaType .eq. META_STEEL) then
            propName(1) = 'C_F1_THETA'
            propName(2) = 'C_F2_THETA'
            propName(3) = 'C_F3_THETA'
            propName(4) = 'C_F4_THETA'
            propName(5) = 'F1_C_THETA'
            propName(6) = 'F2_C_THETA'
            propName(7) = 'F3_C_THETA'
            propName(8) = 'F4_C_THETA'
        elseif (metaType .eq. META_ZIRC) then
            propName(1) = 'C_F1_THETA'
            propName(2) = 'C_F2_THETA'
            propName(3) = 'F1_C_THETA'
            propName(4) = 'F2_C_THETA'
            annealTheta(1:nbProp) = 0.d0
        elseif (metaType .eq. META_STEEL_TEMPER) then
            propName(1) = 'C_F1_THETA'
            propName(2) = 'C_F2_THETA'
            propName(3) = 'C_F3_THETA'
            propName(4) = 'C_F4_THETA'
            propName(5) = 'C_F5_THETA'
            propName(6) = 'C_F6_THETA'
            propName(7) = 'F1_C_THETA'
            propName(8) = 'F2_C_THETA'
            propName(9) = 'F3_C_THETA'
            propName(10) = 'F4_C_THETA'
            propName(11) = 'F5_C_THETA'
            propName(12) = 'F6_C_THETA'
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Get parameters
        call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                    ' ', 'META_RE', 0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 2)
        do iProp = 1, nbProp
            annealTheta(iProp) = propVale(iProp)
        end do
!
! --------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaGetParaVisc
!
! Get parameters for viscosity
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetParaVisc(poum, fami, kpg, ksp, jvMaterCode, &
                               metaType, nbPhase, &
                               eta_, n_, unsurn_, &
                               c_, m_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: poum
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode
        integer(kind=8), intent(in) :: metaType, nbPhase
        real(kind=8), optional, intent(out) :: eta_(nbPhase)
        real(kind=8), optional, intent(out) :: n_(nbPhase)
        real(kind=8), optional, intent(out) :: unsurn_(nbPhase)
        real(kind=8), optional, intent(out) :: c_(nbPhase)
        real(kind=8), optional, intent(out) :: m_(nbPhase)
! ----- Local
        integer(kind=8), parameter :: nbPropMaxi = META_MECA_NBPHASE_MAXI
        real(kind=8) :: propVale(nbPropMaxi)
        integer(kind=8) :: propCode(nbPropMaxi)
        character(len=16) :: propName(nbPropMaxi)
        integer(kind=8) :: nbProp, iProp
! --------------------------------------------------------------------------------------------------
!
        ASSERT(nbPhase .le. META_MECA_NBPHASE_MAXI)
        nbProp = nbPhase
        if (present(eta_)) eta_ = 0.d0
        if (present(n_)) n_ = 20.d0
        if (present(unsurn_)) unsurn_ = 1.d0
        if (present(c_)) c_ = 0.d0
        if (present(m_)) m_ = 20.d0

! ----- Name of parameters
        propName = " "
        if (present(eta_)) then
            if (metaType .eq. META_STEEL) then
                propName(1) = 'F1_ETA'
                propName(2) = 'F2_ETA'
                propName(3) = 'F3_ETA'
                propName(4) = 'F4_ETA'
                propName(5) = 'C_ETA'
            elseif (metaType .eq. META_ZIRC) then
                propName(1) = 'F1_ETA'
                propName(2) = 'F2_ETA'
                propName(3) = 'C_ETA'
            else if (metaType .eq. META_STEEL_TEMPER) then
                propName(1) = 'F1_ETA'
                propName(2) = 'F2_ETA'
                propName(3) = 'F3_ETA'
                propName(4) = 'F4_ETA'
                propName(5) = 'F5_ETA'
                propName(6) = 'F6_ETA'
                propName(7) = 'C_ETA'
            else
                ASSERT(ASTER_FALSE)
            end if
        end if

! ----- Get parameters
        if (present(eta_)) then
            call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                        ' ', 'META_VISC', 0, ' ', [0.d0], &
                        nbProp, propName, propVale, propCode, 2)
            do iProp = 1, nbProp
                eta_(iProp) = propVale(iProp)
            end do
        end if

! ----- Name of parameters
        propName = " "
        if (present(n_)) then
            if (metaType .eq. META_STEEL) then
                propName(1) = 'F1_N'
                propName(2) = 'F2_N'
                propName(3) = 'F3_N'
                propName(4) = 'F4_N'
                propName(5) = 'C_N'
            elseif (metaType .eq. META_ZIRC) then
                propName(1) = 'F1_N'
                propName(2) = 'F2_N'
                propName(3) = 'C_N'
            else if (metaType .eq. META_STEEL_TEMPER) then
                propName(1) = 'F1_N'
                propName(2) = 'F2_N'
                propName(3) = 'F3_N'
                propName(4) = 'F4_N'
                propName(5) = 'F5_N'
                propName(6) = 'F6_N'
                propName(7) = 'C_N'
            else
                ASSERT(ASTER_FALSE)
            end if
        end if

! ----- Get parameters
        if (present(n_)) then
            call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                        ' ', 'META_VISC', 0, ' ', [0.d0], &
                        nbProp, propName, propVale, propCode, 2)
            do iProp = 1, nbProp
                n_(iProp) = propVale(iProp)
                unsurn_(iProp) = 1.d0/n_(iProp)
            end do
        end if

! ----- Name of parameters
        propName = " "
        if (present(c_)) then
            if (metaType .eq. META_STEEL) then
                propName(1) = 'F1_C'
                propName(2) = 'F2_C'
                propName(3) = 'F3_C'
                propName(4) = 'F4_C'
                propName(5) = 'C_C'
            elseif (metaType .eq. META_ZIRC) then
                propName(1) = 'F1_C'
                propName(2) = 'F2_C'
                propName(3) = 'C_C'
            else if (metaType .eq. META_STEEL_TEMPER) then
                propName(1) = 'F1_C'
                propName(2) = 'F2_C'
                propName(3) = 'F3_C'
                propName(4) = 'F4_C'
                propName(5) = 'F5_C'
                propName(6) = 'F6_C'
                propName(7) = 'C_C'
            else
                ASSERT(ASTER_FALSE)
            end if
        end if

! ----- Get parameters
        if (present(c_)) then
            call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                        ' ', 'META_VISC', 0, ' ', [0.d0], &
                        nbProp, propName, propVale, propCode, 2)
            do iProp = 1, nbProp
                c_(iProp) = propVale(iProp)
            end do
        end if

! ----- Name of parameters
        propName = " "
        if (present(m_)) then
            if (metaType .eq. META_STEEL) then
                propName(1) = 'F1_M'
                propName(2) = 'F2_M'
                propName(3) = 'F3_M'
                propName(4) = 'F4_M'
                propName(5) = 'C_M'
            elseif (metaType .eq. META_ZIRC) then
                propName(1) = 'F1_M'
                propName(2) = 'F2_M'
                propName(3) = 'C_M'
            else if (metaType .eq. META_STEEL_TEMPER) then
                propName(1) = 'F1_M'
                propName(2) = 'F2_M'
                propName(3) = 'F3_M'
                propName(4) = 'F4_M'
                propName(5) = 'F5_M'
                propName(6) = 'F6_M'
                propName(7) = 'C_M'
            else
                ASSERT(ASTER_FALSE)
            end if
        end if

! ----- Get parameters
        if (present(m_)) then
            call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                        ' ', 'META_VISC', 0, ' ', [0.d0], &
                        nbProp, propName, propVale, propCode, 2)
            do iProp = 1, nbProp
                m_(iProp) = propVale(iProp)
            end do
        end if
!
! --------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaGetParaPlasTransf
!
! Get plasticity of transformation
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetParaPlasTransf(poum, fami, kpg, ksp, jvMaterCode, &
                                     metaType, nbPhase, nbPhaseCold, &
                                     phaseIncr, zSumCold, &
                                     kpt, fpt)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: poum
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode
        integer(kind=8), intent(in) :: metaType, nbPhase, nbPhaseCold
        real(kind=8), intent(in) :: phaseIncr(nbPhase), zSumCold
        real(kind=8), intent(out) :: kpt(nbPhaseCold), fpt(nbPhaseCold)
! ----- Local
        integer(kind=8), parameter :: nbPropMaxi = META_MECA_NBPHASE_MAXI-1
        real(kind=8) :: propVale(nbPropMaxi)
        integer(kind=8) :: propCode(nbPropMaxi)
        character(len=16) :: propName(nbPropMaxi)
        integer(kind=8) :: nbProp, iProp
! --------------------------------------------------------------------------------------------------
!
        nbProp = nbPhaseCold
        ASSERT(nbPhaseCold .le. META_MECA_NBPHASE_MAXI-1)
        ASSERT(nbPhase .le. META_MECA_NBPHASE_MAXI)
        kpt = 0.d0
        fpt = 0.d0

! ----- Coefficient K
        propName = " "
        if (metaType .eq. META_STEEL) then
            propName(1) = 'F1_K'
            propName(2) = 'F2_K'
            propName(3) = 'F3_K'
            propName(4) = 'F4_K'
        elseif (metaType .eq. META_ZIRC) then
            propName(1) = 'F1_K'
            propName(2) = 'F2_K'
        elseif (metaType .eq. META_STEEL_TEMPER) then
            propName(1) = 'F1_K'
            propName(2) = 'F2_K'
            propName(3) = 'F3_K'
            propName(4) = 'F4_K'
            propName(5) = 'F5_K'
            propName(6) = 'F6_K'
        else
            ASSERT(ASTER_FALSE)
        end if
        call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                    ' ', 'META_PT', 0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 2)
        do iProp = 1, nbProp
            kpt(iProp) = propVale(iProp)
        end do

! ----- Function F'
        propName = " "
        if (metaType .eq. META_STEEL) then
            propName(1) = 'F1_D_F_META'
            propName(2) = 'F2_D_F_META'
            propName(3) = 'F3_D_F_META'
            propName(4) = 'F4_D_F_META'
        elseif (metaType .eq. META_ZIRC) then
            propName(1) = 'F1_D_F_META'
            propName(2) = 'F2_D_F_META'
        elseif (metaType .eq. META_STEEL_TEMPER) then
            propName(1) = 'F1_D_F_META'
            propName(2) = 'F2_D_F_META'
            propName(3) = 'F3_D_F_META'
            propName(4) = 'F4_D_F_META'
            propName(5) = 'F5_D_F_META'
            propName(6) = 'F6_D_F_META'
        else
            ASSERT(ASTER_FALSE)
        end if

        do iProp = 1, nbProp
            if (phaseIncr(iProp) .gt. 0.d0) then
                call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                            ' ', 'META_PT', 1, 'META', [zSumCold], &
                            1, propName(iProp), propVale(iProp), propCode, 2)
                fpt(iProp) = propVale(iProp)
            else
                fpt(iProp) = 0.d0
            end if
        end do
!
! --------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaGetParaHardLine
!
! Get hardening slope (linear)
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetParaHardLine(poum, fami, kpg, ksp, jvMaterCode, &
                                   metaType, nbPhase, &
                                   young, hardCoef, h)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: poum
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode
        integer(kind=8), intent(in) :: metaType, nbPhase
        real(kind=8), intent(in) :: young, hardCoef
        real(kind=8), intent(out) :: h(nbPhase)
! ----- Local
        integer(kind=8), parameter :: nbPropMaxi = META_MECA_NBPHASE_MAXI
        real(kind=8) :: propVale(nbPropMaxi)
        integer(kind=8) :: propCode(nbPropMaxi)
        character(len=16) :: propName(nbPropMaxi)
        integer(kind=8) :: nbProp, iProp
! --------------------------------------------------------------------------------------------------
!
        nbProp = nbPhase
        ASSERT(nbPhase .le. META_MECA_NBPHASE_MAXI)
        h = 0.d0
!
        propName = " "
        if (metaType .eq. META_STEEL) then
            propName(1) = 'F1_D_SIGM_EPSI'
            propName(2) = 'F2_D_SIGM_EPSI'
            propName(3) = 'F3_D_SIGM_EPSI'
            propName(4) = 'F4_D_SIGM_EPSI'
            propName(5) = 'C_D_SIGM_EPSI'
        elseif (metaType .eq. META_ZIRC) then
            propName(1) = 'F1_D_SIGM_EPSI'
            propName(2) = 'F2_D_SIGM_EPSI'
            propName(3) = 'C_D_SIGM_EPSI'
        elseif (metaType .eq. META_STEEL_TEMPER) then
            propName(1) = 'F1_D_SIGM_EPSI'
            propName(2) = 'F2_D_SIGM_EPSI'
            propName(3) = 'F3_D_SIGM_EPSI'
            propName(4) = 'F4_D_SIGM_EPSI'
            propName(5) = 'F5_D_SIGM_EPSI'
            propName(6) = 'F6_D_SIGM_EPSI'
            propName(7) = 'C_D_SIGM_EPSI'
        else
            ASSERT(ASTER_FALSE)
        end if
!
        call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                    ' ', 'META_ECRO_LINE', 0, ' ', [0.d0], &
                    nbPhase, propName, propVale, propCode, 2)
        do iProp = 1, nbProp
            h(iProp) = propVale(iProp)
            h(iProp) = hardCoef*h(iProp)*young/(young-h(iProp))
        end do
!
! --------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaGetParaHardTrac
!
! Get hardening slope (non-linear)
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetParaHardTrac(jvMaterCode, metaType, nbPhase, &
                                   l_temp, temp, &
                                   epseq, h0, rp_, nbValeMaxi_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvMaterCode
        integer(kind=8), intent(in) :: metaType
        integer(kind=8), intent(in) :: nbPhase
        aster_logical, intent(in) :: l_temp
        real(kind=8), intent(in) :: temp
        real(kind=8), intent(in) :: epseq(nbPhase)
        real(kind=8), intent(out) :: h0(nbPhase)
        real(kind=8), optional, intent(out) :: rp_(nbPhase)
        integer(kind=8), optional, intent(out) :: nbValeMaxi_
! ----- Parameters
        character(len=16), parameter :: tractionName = 'META_TRACTION'
        integer(kind=8), parameter :: nbPropMaxi = META_MECA_NBPHASE_MAXI
        character(len=8) :: propName(nbPropMaxi)
        integer(kind=8) :: iPhase
        integer(kind=8) :: jvProl, jvVale
        integer(kind=8) :: nbVale, nbValeMaxi
        character(len=8) :: paraName
        real(kind=8) :: paraVale, young
! --------------------------------------------------------------------------------------------------
!
        h0 = 0.d0
        ASSERT(nbPhase .le. META_MECA_NBPHASE_MAXI)
        if (present(rp_)) rp_ = 0.d0

! ----- Select name of properties for traction curve
        propName = " "
        if (metaType .eq. META_STEEL) then
            propName(1) = 'SIGM_F1'
            propName(2) = 'SIGM_F2'
            propName(3) = 'SIGM_F3'
            propName(4) = 'SIGM_F4'
            propName(5) = 'SIGM_C'
        elseif (metaType .eq. META_ZIRC) then
            propName(1) = 'SIGM_F1'
            propName(2) = 'SIGM_F2'
            propName(3) = 'SIGM_C'
        elseif (metaType .eq. META_STEEL_TEMPER) then
            propName(1) = 'SIGM_F1'
            propName(2) = 'SIGM_F2'
            propName(3) = 'SIGM_F3'
            propName(4) = 'SIGM_F4'
            propName(5) = 'SIGM_F5'
            propName(6) = 'SIGM_F6'
            propName(7) = 'SIGM_C'
        else
            ASSERT(ASTER_FALSE)
        end if

        nbValeMaxi = -1
        do iPhase = 1, nbPhase
! --------- Get informations (value and type) about parameters for traction curve
            call rctype(jvMaterCode, 1, 'TEMP', [temp], &
                        paraVale, paraName, &
                        keyw_factz=tractionName, keywz=propName(iPhase))
            if ((paraName .eq. 'TEMP') .and. (.not. l_temp)) then
                call utmess('F', 'COMPOR5_5', sk=paraName)
            end if

! --------- Get access of traction curve
            call rctrac(jvMaterCode, 2, propName(iPhase), temp, &
                        jvProl, jvVale, nbVale, young)

            if (present(rp_)) then
                call rcfonc('V', 2, jvProl, jvVale, nbVale, &
                            p=epseq(iPhase), rp=rp_(iPhase), rprim=h0(iPhase))
            else
                call rcfonc('V', 2, jvProl, jvVale, nbVale, &
                            p=epseq(iPhase), rprim=h0(iPhase))
            end if
            if (nbVale .ge. nbValeMaxi) then
                nbValeMaxi = nbVale
            end if
        end do
        if (present(nbValeMaxi_)) then
            nbValeMaxi_ = nbValeMaxi
        end if
!
! --------------------------------------------------------------------------------------------------
    end subroutine
!
end module MetallurgyMeca_module
