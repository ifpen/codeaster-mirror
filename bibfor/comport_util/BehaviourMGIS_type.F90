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
! Types for the management of integration of behaviour
!
! ==================================================================================================
!
module BehaviourMGIS_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterfort/BehaviourMGIS_type.h"
! ==================================================================================================
! Global variables
! ==================================================================================================
! - Correspondance between MFront and code_aster names
    character(len=64), parameter :: fromAsterToMFront(2, ESVA_EXTE_MGIS_NBMAXI) = &
                                    reshape((/"ELTSIZE1            ", "ElementSize         ", &
                                              "HYDR                ", "ConcreteHydration   ", &
                                              "HYGR                ", "Hygrometry          ", &
                                              "PCAP                ", "CapillaryPressure   ", &
                                              "BAINITE             ", "BainitePhaseRatio   ", &
                                              "FERRITE             ", "FerritePhaseRatio   ", &
                                              "MARTENS             ", "MartensitePhaseRatio", &
                                              "PERLITE             ", "PerlitePhaseRatio   ", &
                                              "SECH                ", "ConcreteDrying      ", &
                                              "TEMP                ", "Temperature         ", &
                                              "TEMPREFE            ", "ReferenceTemperature", &
                                              "TIME                ", "Time                "/), &
                                            (/2, ESVA_EXTE_MGIS_NBMAXI/))

!===================================================================================================
    public :: fromAsterToMFront
contains
!===================================================================================================
end module
