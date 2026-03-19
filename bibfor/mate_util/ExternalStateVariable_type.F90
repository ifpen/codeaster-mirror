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
! Types for External State Variables
!
! ==================================================================================================
!
module ExternalStateVariable_type
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
#include "asterf_types.h"
#include "asterfort/ExternalStateVariable_type.h"
! ==================================================================================================
! ==================================================================================================
!
! Global variables - General
!
! List of all external state variables: name, physical quantity and number of components
!
! ==================================================================================================
    character(len=8), parameter :: listExteStatVari(EXTEVARI_NB_MAXI) = (/"TEMP    ", "GEOM    ", &
                                                                          "CORR    ", "IRRA    ", &
                                                                          "HYDR    ", "SECH    ", &
                                                                          "EPSA    ", "M_ACIER ", &
                                                                          "M_ZIRC  ", "NEUT1   ", &
                                                                          "NEUT2   ", "NEUT3   ", &
                                                                          "PTOT    "/)
    character(len=8), parameter :: listPhysQuantity(EXTEVARI_NB_MAXI) = (/"TEMP_R  ", "GEOM_R  ", &
                                                                          "CORR_R  ", "IRRA_R  ", &
                                                                          "HYDR_R  ", "TEMP_R  ", &
                                                                          "EPSI_R  ", "VARI_R  ", &
                                                                          "VARI_R  ", "NEUT_R  ", &
                                                                          "NEUT_R  ", "NEUT_R  ", &
                                                                          "DEPL_R  "/)
    integer(kind=8), parameter :: listNbCmp(EXTEVARI_NB_MAXI) = (/7, 3, &
                                                                  1, 1, &
                                                                  1, 1, &
                                                                  6, 9, &
                                                                  5, 1, &
                                                                  1, 1, &
                                                                  1/)
! ==================================================================================================
!
! Derivated types
!
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! TYPE - EXTE_VARI_CMP
!
! Components of an external state variable
!
! --------------------------------------------------------------------------------------------------
    type EXTE_VARI_CMP
! ----- Name of component in physical quantity catalog
        character(len=8) :: physQuantityCmp = " "
! ----- Name of component in external state variable
        character(len=8) :: nameCmp = " "
    end type EXTE_VARI_CMP
! --------------------------------------------------------------------------------------------------
!
! TYPE - EXTE_VARI_CATA
!
! External state variables available in catalog
!
! --------------------------------------------------------------------------------------------------
    type EXTE_VARI_CATA
! ----- Name of external state variable
        character(len=8) :: name = " "
! ----- Type of field (DEPL, SECH, META_ELNO, ...)
        character(len=16) :: fieldType = " "
! ----- Physical quantity
        character(len=8) :: physQuantity = " "
! ----- List of components
        integer(kind=8) :: nbCmp = 0
        type(EXTE_VARI_CMP), pointer :: listCmp(:) => null()
    end type EXTE_VARI_CATA

! --------------------------------------------------------------------------------------------------
!
! TYPE - EXTE_VARI_DESC
!
! Descriptor of an external state variable
!
! --------------------------------------------------------------------------------------------------
    type EXTE_VARI_DESC
! ----- Index in catalog
        integer(kind=8) :: cataIndx = 0
! ----- Name of external state variable
        character(len=8) :: exteVariName = " "
! ----- Origin of variable: EVOL, CHAMP
        character(len=8) :: affeType = " "
! ----- Value of *_REFE (TEMP_REFE for instance)
        real(kind=8)  :: valeRefe = 0.d0
! ----- Type of field (DEPL, SECH, META_ELNO, ...)
        character(len=16) :: fieldType = " "
! ----- Name of datastructure from user
        character(len=8) :: dsUser = " "
! ----- PROL_GAUCHE when AFFE_TYPE='EVOL'
        character(len=16) :: funcExtrLeft = " "
! ----- PROL_DROITE when AFFE_TYPE='EVOL'
        character(len=16) :: funcExtrRight = " "
! ----- FONC_INST when AFFE_TYPE='EVOL'
        character(len=8) :: funcResult = " "
    end type EXTE_VARI_DESC
! --------------------------------------------------------------------------------------------------
!
! TYPE - EXTE_VARI_AFFE
!
! External state variables assigned by user
!
! --------------------------------------------------------------------------------------------------
    type EXTE_VARI_AFFE
! ----- Descriptor of external state variables
        type(EXTE_VARI_DESC), pointer :: exteVariList(:) => null()
! ----- Number of external state variables assigned by user
        integer(kind=8) :: nbAffe = 0
! ----- Total number of components for all external state variables
        integer(kind=8) :: nbCmpTotal = 0
    end type EXTE_VARI_AFFE
!
end module ExternalStateVariable_type
