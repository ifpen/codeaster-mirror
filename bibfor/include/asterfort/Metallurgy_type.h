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
! ==================================================================================================
!
! Metallurgical phases
!
! ==================================================================================================
! Type of material
#define META_NONE         0
#define META_STEEL        1
#define META_ZIRC         2

! Number of phases for steel (total and cold)
#define PSTEEL_NB         5
#define PCSTEEL_NB        4

! Index of phases for steel
#define PFERRITE          1
#define PPERLITE          2
#define PBAINITE          3
#define PMARTENS          4
#define PAUSTENITE        5
#define PSUMCOLD          6

! Number of phases for steel with tempering
#define PRSTEEL_NB        7

! Index of phases for steel with tempering
#define PRFERRITE         1
#define PRPERLITE         2
#define PRBAINITEB        3
#define PRBAINITER        4
#define PRMARTENSB        5
#define PRMARTENSR        6
#define PRAUSTENITE       7
#define PRSUMCOLD         8

! Number of phases for zircaloy
#define PZIRC_NB          3

! Index of phases for zircaloy
#define PALPHA1           1
#define PALPHA2           2
#define PBETA             3

! ==================================================================================================
!
! Metallurgical laws (CALC_META)
!
! ==================================================================================================
! Map for metallurgical behaviour <COMPOR_META>: size of map
#define COMPORMETA_SIZE   5

! Map for metallurgical behaviour <COMPOR_META>: slots in map
#define ZMETATYPE         1
#define ZNBVARI           2
#define ZMETALAW          3
#define ZNUMECOMP         4
#define ZNBPHASE          5

! Maximum number of phases from user
#define META_META_NBPHASE_MAXI 12

! ==================================================================================================
! Standard steel
! ==================================================================================================
! Total number of phases (4 cold, 1 hot + 1 sum of cold)
#define NBPHASESTEEL      6

! Index of internal state variable for standard steel
#define SIZE_GRAIN        7
#define STEEL_TEMP        8
#define TEMP_MARTENSITE   9

! Number of internal state variables required for initial state (nb phases + size of grain)
#define PVARIINIT         7

! Number of internal state variables for standard steel
#define NBVARISTEEL       9

! ==================================================================================================
! Steel with tempering
! ==================================================================================================
! Total number of phases (6 cold, 1 hot + 1 sum of cold)
#define NBPHASESTEELR     8

! Index of internal state variable for tempering steel
#define SIZE_GRAINR       9
#define STEEL_TEMPR       10
#define TEMP_MARTENSITER  11
#define THER_CYCL         12

! Number of internal states variables required for initial state
#define PRVARIINIT        9

! Number of internal state variables for tempering steel
#define NBVARISTEELR      12

! ==================================================================================================
! Zircaloy
! ==================================================================================================
! For next ones: add total number of phases to access in internal state variable vector
#define ZIRC_TEMP        1
#define TIME_TRAN        2

! Number of internal state variables for Zircaloy
#define NBVARIZIRC       5

! ==================================================================================================
! Parameters for algorithm
! ==================================================================================================
! Size of vector for input parameters for tempering steel
#define NB_PARAIN_TEMPER 6

! Kinetic
#define COOLING          0
#define HEATING          1

! Parameter for integration : bound to cut
#define STEEL_MIN_CUT 1.d-3
! Parameter for integration : Maximum temperature step
#define STEEL_MAX_TEMP_STEP 5.d0
! Parameter for integration : Maximum number of sub-steps
#define STEEL_MAX_NB_STEP 10

! ==================================================================================================
!
! Mechanical laws with metallurgical phases
!
! ==================================================================================================
! Maximum number of phases from user
#define META_MECA_NBPHASE_MAXI 5

! - Index for internal state variables - Isotropic hardening: updated elasticity yield
#define IDX_I_SIGY       6
#define IDX_I_IPLAS      7

! - Index for internal state variables - Isotropic hardening: updated elasticity yield
#define IDX_C_IPLAS      37
