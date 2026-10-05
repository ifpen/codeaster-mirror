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
! Types for the management strains
!
! ==================================================================================================
!
module BehaviourStrain_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/ElasticityMaterial_type.h"
! ==================================================================================================
! Global variables
! ==================================================================================================

! - List of external state variables for anelastic strains
    character(len=8), parameter :: varcStrainName(VARC_STRAIN_NBMAXI) = &
                                   (/'TEMP', 'SECH', 'HYDR', 'EPSA', 'PTOT'/)
    integer(kind=8), parameter :: varcStrainNbcmp(VARC_STRAIN_NBMAXI) = &
                                  (/1, 1, 1, VARC_EPSA_NBCMP, 1/)
    aster_logical, parameter :: varcStrainHasRefe(VARC_STRAIN_NBMAXI) = &
        (/.true._1, .true._1, .false._1, .false._1, .false._1/)
    character(len=8), parameter :: varcAnelName(VARC_EPSA_NBCMP) = &
                                   (/'EPSAXX', 'EPSAYY', 'EPSAZZ', 'EPSAXY', 'EPSAXZ', 'EPSAYZ'/)
    character(len=16), parameter :: varcRHSOption(VARC_STRAIN_NBMAXI) = &
                                    (/'CHAR_MECA_TEMP_R', 'CHAR_MECA_SECH_R', &
                                      'CHAR_MECA_HYDR_R', 'CHAR_MECA_EPSA_R', &
                                      'CHAR_MECA_PTOT_R'/)

! ==================================================================================================
! Type: External state variables for anelastic strains
! ==================================================================================================
    type Varc_Strain
! ----- Type of External state variable
        integer(kind=8) :: varcStrainType = VARC_STRAIN_NONE
! ----- Flag
        aster_logical :: exist = ASTER_FALSE
! ----- Values of external state variables
        real(kind=8) :: varcPrev(6) = 0.d0
        real(kind=8) :: varcCurr(6) = 0.d0
        real(kind=8) :: varcIncr(6) = 0.d0
        real(kind=8) :: varcRefe = 0.d0
! ----- Values of field
        real(kind=8) :: fieldPrev(6) = 0.d0
        real(kind=8) :: fieldCurr(6) = 0.d0
        real(kind=8) :: fieldIncr(6) = 0.d0
    end type Varc_Strain
! ==================================================================================================
! Type: All external state variables for inelastic strains
! ==================================================================================================
    type All_Varc_Strain
! ----- Flag for THM
        aster_logical :: lTHM = ASTER_FALSE
! ----- Flag for inelastic strains from external state variables
        aster_logical :: hasInelasticStrains = ASTER_FALSE
! ----- Flag to give time to evaluate
        aster_logical :: hasTime = ASTER_FALSE
        real(kind=8) :: time = 0.d0
! ----- Flag to give specific temperature to evaluate
        aster_logical :: hasTemp = ASTER_FALSE
        real(kind=8) :: temp = 0.d0
! ----- List of all inelastic strains and external states variables
        type(Varc_Strain) :: list(VARC_STRAIN_NBMAXI)
    end type All_Varc_Strain
!===================================================================================================
    public :: All_Varc_Strain, Varc_Strain
    public :: varcStrainName, varcStrainNbcmp, varcStrainHasRefe, varcAnelName, varcRHSOption
contains
!===================================================================================================
end module
