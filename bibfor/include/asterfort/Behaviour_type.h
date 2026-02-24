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
! aslint: disable=C1509
!
!
! --------------------------------------------------------------------------------------------------
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         If this file is modified, bibcxx/include/Behaviour_type.h must also be modified
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! What to compute from option ?
!
! L_VARI: we have internal state variables
! L_SIGM: we have stress and return code error
! L_VECT: compute internal forces vector
! L_MATR: compute tangent matrix
! L_PRED: prediction phase of Newton
! L_CORR: correction phase of Newton
!
! --------------------------------------------------------------------------------------------------
!
#define L_VARI(option) (option(1:9).eq.'RAPH_MECA' .or. option(1:9).eq.'FULL_MECA')
#define L_SIGM(option) (option(1:9).eq.'RAPH_MECA' .or. option(1:9).eq.'FULL_MECA' .or. option .eq. 'RIGI_MECA_TANG')
#define L_VECT(option) (option(1:9).eq.'RAPH_MECA' .or. option(1:9).eq.'FULL_MECA' .or. option .eq. 'RIGI_MECA_TANG')
#define L_MATR(option) (option(1:9).eq.'FULL_MECA' .or. option(1:9).eq.'RIGI_MECA')
#define L_PRED(option) (option .eq. 'RIGI_MECA_TANG')
#define L_CORR(option) (option .ne. 'RIGI_MECA_TANG')
#define L_MATR_PRED(option) (option(1:4) .eq. 'RIGI')
!
! --------------------------------------------------------------------------------------------------
!
! The field <COMPOR>
!   List of strings (K16) for each cell
!   Definition of behaviour (relation, strain model, number of internal state vari, etc.)
!
! --------------------------------------------------------------------------------------------------
!
! Size
!
#define COMPOR_SIZE 25
!
! Slots: general
!
#define RELA_NAME 1
#define NVAR 2
#define DEFO 3
#define INCRELAS 4
#define PLANESTRESS 5
#define NUME 6
#define MULTCOMP 7
#define POSTITER 8
#define DEFO_LDC 21
#define RIGI_GEOM 22
#define REGUVISC 23
#define MGIS_ADDR 24
#define POSTINCR 25
!
! Slots: for KIT
!
#define KIT1_NAME 9
#define KIT2_NAME 10
#define KIT3_NAME 11
#define KIT4_NAME 12
#define KIT1_NUME 13
#define KIT2_NUME 14
#define KIT3_NUME 15
#define KIT4_NUME 16
#define KIT1_NVAR 17
#define KIT2_NVAR 18
#define KIT3_NVAR 19
#define KIT4_NVAR 20
!
! Slots: for KIT_THM
!
#define MECA_NAME 9
#define HYDR_NAME 10
#define THER_NAME 11
#define THMC_NAME 12
#define THMC_NUME 13
#define THER_NUME 14
#define HYDR_NUME 15
#define MECA_NUME 16
#define THMC_NVAR 17
#define THER_NVAR 18
#define HYDR_NVAR 19
#define MECA_NVAR 20
!
! Slots: for KIT_DDI
!
#define CREEP_NAME 9
#define PLAS_NAME 10
#define COUPL_NAME 11
#define CPLA_NAME 12
#define CREEP_NUME 16
#define PLAS_NUME 15
#define CREEP_NVAR 17
#define PLAS_NVAR 18
!
! Slots: for KIT_META
!
#define META_PHAS 9
#define META_RELA 10
#define META_GLOB 11
!
! Slots: for KIT_CG
!
#define CABLE_NAME 9
#define SHEATH_NAME 10
#define CABLE_NUME 13
#define SHEATH_NUME 14
#define CABLE_NVAR 17
#define SHEATH_NVAR 18

! --------------------------------------------------------------------------------------------------
!
! The field <CARCRI>
!   List of real for each cell
!   Parameters to solve behaviour equations
!
! --------------------------------------------------------------------------------------------------
! Size
#define CARCRI_SIZE 22

! Slots
#define ITER_INTE_MAXI 1
#define TYPE_MATR_T 2
#define RESI_INTE 3
#define PARM_THETA 4
#define ITER_INTE_PAS 5
#define ALGO_INTE_R 6
#define VALE_PERT_RELA 7
#define RESI_DEBORST_MAX 8
#define ITER_DEBORST_MAX 9
#define RESI_RADI_RELA 10
#define IVARIEXT1 11
#define PARM_THETA_THM 12
#define IPOSTITER 13
#define CARCRI_MATRSYME 17
#define PARM_ALPHA_THM 18
#define IVARIEXT2 22
! Slots for external solvers (UMAT/MFRONT)
! Pointer to MGISBehaviour (MFront) or function (UMAT)
#define EXTE_PTR 16
! 1 for MFRONT official, 2 for MFRONT proto, 4 for UMAT (default: 0 internal)
#define EXTE_TYPE 15
#define SOLV_BEHAV_ASTER 0
#define SOLV_BEHAV_MGIS_OFFI 1
#define SOLV_BEHAV_MGIS_PROTO 2
#define SOLV_BEHAV_UMAT 4

!       Strain model for (MFRONT only)
#define EXTE_STRAIN 21

! --------------------------------------------------------------------------------------------------
!
! For external state variables
!
! --------------------------------------------------------------------------------------------------

! Number of different type of external states variables
#define VARC_EXTE_NBTYPE 31

! Maximum number of external state variables in external solvers
#define VARC_EXTE_NBMAXI 8

! External state variables to create anelastic strain
#define VARC_STRAIN_NBMAXI 5
#define VARC_STRAIN_NONE 0
#define VARC_STRAIN_TEMP 1
#define VARC_STRAIN_SECH 2
#define VARC_STRAIN_HYDR 3
#define VARC_STRAIN_EPSA 4
#define VARC_STRAIN_PTOT 5
#define VARC_STRAIN_ALL 99

! Number of components for EPSA
#define VARC_EPSA_NBCMP 6

! Maximum number of Gauss points for coordinates
#define VARC_GEOM_NBMAXI 128

! --------------------------------------------------------------------------------------------------
!
! For strains
!
! --------------------------------------------------------------------------------------------------
#define STRAIN_TYPE_NONE 0
#define STRAIN_TYPE_SMALL 1
#define STRAIN_TYPE_GREEN 2
#define STRAIN_TYPE_LOG 3

! Type of external state variables
#define ELTSIZE1 1
#define XXXXXXXX 2
#define GRADVELO 3
#define HYGR 4
#define NEUT1 5
#define NEUT2 6
#define TEMP 7
#define DTX 8
#define DTY 9
#define DTZ 10
#define X 11
#define Y 12
#define Z 13
#define SECH 14
#define HYDR 15
#define CORR 16
#define IRRA 17
#define EPSAXX 18
#define EPSAYY 19
#define EPSAZZ 20
#define EPSAXY 21
#define EPSAXZ 22
#define EPSAYZ 23
#define ZFERRITE 24
#define ZPERLITE 25
#define ZBAINITE 26
#define ZMARTENS 27
#define ZALPHPUR 28
#define ZALPHBET 29
#define TIME 30
#define TEMPREFE 31

! Set to 1 to activate DEBUG (careful, very verbose, at each Gauss point !)
#define LDC_PREP_DEBUG 0

! --------------------------------------------------------------------------------------------------
! Error during integration of behaviour
!
! 0 => No problem
! 1 => convergence default
! 2 => quality problem
! 3 => stress plane algorithm not converged
! 4 => out of bound for validity
!
! --------------------------------------------------------------------------------------------------
#define LDC_ERROR_NONE 0
#define LDC_ERROR_NCVG 1
#define LDC_ERROR_QUAL 2
#define LDC_ERROR_CPLA 3
#define LDC_ERROR_DVAL 4

! --------------------------------------------------------------------------------------------------
! Local time step division
!
! 0 => No time step division
! 1 => Forced time step division
! 2 => Automatic time step division
!
! --------------------------------------------------------------------------------------------------
#define LDC_TIMEDIV_NONE 0
#define LDC_TIMEDIV_FORC 1
#define LDC_TIMEDIV_AUTO 2
