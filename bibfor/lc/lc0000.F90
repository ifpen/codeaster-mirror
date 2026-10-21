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
! aslint: disable=W1306,W1501,W1504,C1505
!
subroutine lc0000(BEHInteg, &
                  ndim, option, typmod, &
                  instam, instap, &
                  compor, carcri, multComp, &
                  neps, epsm_tot, deps_tot, &
                  nsig, sigm_all, &
                  nvi_all, vim, &
                  sigp, vip, &
                  ndsde, dsidep, codret, &
                  l_epsi_varc, numlc)
!
    use Behaviour_type
    use Behaviour_module
    use MetallurgyMeca_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lc0001.h"
#include "asterfort/lc0002.h"
#include "asterfort/lc0003.h"
#include "asterfort/lc0004.h"
#include "asterfort/lc0007.h"
#include "asterfort/lc0008.h"
#include "asterfort/lc0009.h"
#include "asterfort/lc0015.h"
#include "asterfort/lc0016.h"
#include "asterfort/lc0017.h"
#include "asterfort/lc0018.h"
#include "asterfort/lc0019.h"
#include "asterfort/lc0021.h"
#include "asterfort/lc0022.h"
#include "asterfort/lc0023.h"
#include "asterfort/lc0024.h"
#include "asterfort/lc0025.h"
#include "asterfort/lc0026.h"
#include "asterfort/lc0028.h"
#include "asterfort/lc0029.h"
#include "asterfort/lc0030.h"
#include "asterfort/lc0031.h"
#include "asterfort/lc0032.h"
#include "asterfort/lc0033.h"
#include "asterfort/lc0034.h"
#include "asterfort/lc0035.h"
#include "asterfort/lc0036.h"
#include "asterfort/lc0038.h"
#include "asterfort/lc0040.h"
#include "asterfort/lc0042.h"
#include "asterfort/lc0050.h"
#include "asterfort/lc0054.h"
#include "asterfort/lc0055.h"
#include "asterfort/lc0058.h"
#include "asterfort/lc0059.h"
#include "asterfort/lc0060.h"
#include "asterfort/lc0062.h"
#include "asterfort/lc0075.h"
#include "asterfort/lc0076.h"
#include "asterfort/lc0077.h"
#include "asterfort/lc0078.h"
#include "asterfort/lc0079.h"
#include "asterfort/lc0101.h"
#include "asterfort/lc0102.h"
#include "asterfort/lc0120.h"
#include "asterfort/lc0137.h"
#include "asterfort/lc0145.h"
#include "asterfort/lc0165.h"
#include "asterfort/lc0166.h"
#include "asterfort/lc0167.h"
#include "asterfort/lc0168.h"
#include "asterfort/lc0169.h"
#include "asterfort/lc1002.h"
#include "asterfort/lc1015.h"
#include "asterfort/lc1037.h"
#include "asterfort/lc1137.h"
#include "asterfort/lc2001.h"
#include "asterfort/lc2002.h"
#include "asterfort/lc2036.h"
#include "asterfort/lc4047.h"
#include "asterfort/lc5001.h"
#include "asterfort/lc5002.h"
#include "asterfort/lc5003.h"
#include "asterfort/lc5024.h"
#include "asterfort/lc5038.h"
#include "asterfort/lc5076.h"
#include "asterfort/lc5090.h"
#include "asterfort/lc5091.h"
#include "asterfort/lc5101.h"
#include "asterfort/lc5102.h"
#include "asterfort/lc5152.h"
#include "asterfort/lc6036.h"
#include "asterfort/lc6046.h"
#include "asterfort/lc6057.h"
#include "asterfort/lc6058.h"
#include "asterfort/lc6075.h"
#include "asterfort/lc6076.h"
#include "asterfort/lc7010.h"
#include "asterfort/lc7011.h"
#include "asterfort/lc7013.h"
#include "asterfort/lc7045.h"
#include "asterfort/lc7046.h"
#include "asterfort/lc7047.h"
#include "asterfort/lc7048.h"
#include "asterfort/lc7058.h"
#include "asterfort/lc8028.h"
#include "asterfort/lc8029.h"
#include "asterfort/lc8057.h"
#include "asterfort/lc8146.h"
#include "asterfort/lc8331.h"
#include "asterfort/lc9040.h"
#include "asterfort/lc9041.h"
#include "asterfort/lc9043.h"
#include "asterfort/lc9049.h"
#include "asterfort/lc9051.h"
#include "asterfort/lc9056.h"
#include "asterfort/lc9058.h"
#include "asterfort/lc9077.h"
#include "asterfort/lcvisc.h"
#include "asterfort/utmess.h"
#include "asterfort/lc9078.h"
#include "asterfort/lc9501.h"
#include "asterfort/lc9502.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    integer(kind=8), intent(in) :: ndim
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(in) :: instam, instap
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16), intent(in) :: multComp
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm_tot(neps), deps_tot(neps)
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: sigm_all(nsig)
    integer(kind=8), intent(in) :: nvi_all
    real(kind=8), intent(in) :: vim(nvi_all)

    real(kind=8), intent(inout) :: sigp(nsig)
    real(kind=8), intent(inout) :: vip(nvi_all)
    integer(kind=8), intent(in) :: ndsde
    real(kind=8), intent(inout) :: dsidep(merge(nsig, 6, nsig*neps .eq. ndsde), &
                                          merge(neps, 6, nsig*neps .eq. ndsde))
    integer(kind=8), intent(out):: codret
    aster_logical, intent(in) :: l_epsi_varc
    integer(kind=8), intent(in) :: numlc
!
! --------------------------------------------------------------------------------------------------
!
! Mechanical non-linear behaviours
!
! Main switch to the integration of behaviour laws
!
! --------------------------------------------------------------------------------------------------
!
! IO  BEHInteg         : parameters for integration of behaviour
! In  option           : option to compute
! In  typmod           : type of modeling (3D, 2D, etc.)
! In  instam           : time at beginning of current time step
! In  instap           : time at end of current time step
! In  compor           : description of behaviour
! In  carcri           : parameters for integration of behaviour
! In  multComp         : name of JEVEUX object for multi-behaviour (DEFI_COMPOR)
! In  neps             : size of strain tensor
! In  epsm_inp         : strain tensor at beginning of current time step
! In  deps_inp         : increment of strain tensor from beginning of current time step
! In  nsig             : size of stress tensor
! In  sigm             : stress tensor at beginning of current time step
! In  vim              : internal state variables at beginning of current time step
! IO  sigp             : stress tensor at end of current time step
! IO  vip              : internal state variables at end of current time step
! In  ndsde            : size of jacobian matrix (dSig/dEps)
! IO  dsidep           : jacobian matrix (dSig/dEps)
! IO  codret           : return code from integration of behaviour
!     LDC_ERROR_NONE => No problem
!     LDC_ERROR_NCVG => convergence default
!     LDC_ERROR_QUAL => quality problem
!     LDC_ERROR_CPLA => stress plane algorithm not converged
!     LDC_ERROR_DVAL => out of bound for validity
! In  l_epsi_varc      : flag to compute non-mechanical strains (from external state variables)
! In  num_lc           : index of behaviour to integrate
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8), dimension(6), parameter:: r2 = [1.d0, 1.d0, 1.d0, rac2, rac2, rac2]
    integer(kind=8), parameter :: nvi_regu_visc = 8, nvi_gdef_log = 6
    aster_logical :: lHardIsot, lHardKine, lHardMixed
    character(len=16) :: relaComp
    integer(kind=8) :: nviRestEcro
    integer(kind=8):: nvi, idx_regu_visc, numlcEff, ndimsi
    real(kind=8):: sigm(nsig), epsm(neps), deps(neps)
    integer(kind=8) :: ndt, ndi
    common/tdim/ndt, ndi
    integer(kind=8) :: kpg, ksp, jvMaterCode
    character(len=8) :: fami
!
! --------------------------------------------------------------------------------------------------
!

! - Size of tensors
    ndimsi = BEHInteg%behavPara%ndimsi
    ndt = ndimsi
    ndi = ndim

! - Detect external state variables
    call detectVarc(BEHInteg)

! - Prepare external state variables at Gauss point
    call behaviourPrepESVAPoin(BEHInteg)

! - Prepare input strains for the behaviour law
! - Default: mechanical strains are total strains (no external state variables)
    epsm = epsm_tot
    deps = deps_tot
    call behaviourPrepStrain(BEHInteg, epsm, deps)

! - Prepare external state variables for external solvers (UMAT/MFRONT)
    if (BEHInteg%behavPara%lExteSolver) then
        call behaviourPrepESVAExte(BEHInteg)
    end if

! - How many internal variables really for the constitutive law ?
    nvi = nvi_all
    if (BEHInteg%behavPara%lGdefLog) then
        nvi = nvi-nvi_gdef_log
    end if
    if (BEHInteg%behavPara%lReguVisc) then
        nvi = nvi-nvi_regu_visc
        idx_regu_visc = nvi+1
    end if
    if (BEHinteg%behavPara%lAnnealing) then
        relaComp = compor(RELA_NAME)
        call metaAnnealGetType(relaComp, lHardIsot, lHardKine, lHardMixed, nviRestEcro)
        nvi = nvi-nviRestEcro
    end if
    ASSERT(nvi .ge. 1)

! - What is the stress at t- for the constitutive law ?
    sigm = sigm_all
    if (BEHInteg%behavPara%lReguVisc) then
        sigm(1:ndimsi) = sigm(1:ndimsi)-vim(idx_regu_visc:idx_regu_visc-1+ndimsi)*r2(1:ndimsi)
    end if

! - Initializations of output variables
    codret = 0
    if (BEHInteg%behavPara%lSigm) then
        sigp = 0.d0
    end if
    if (BEHInteg%behavPara%lMatr) then
        dsidep = 0.d0
    end if
    if (BEHInteg%behavPara%lVari .and. BEHInteg%behavPara%lAnnealing) then
        vip(nvi_all) = vim(nvi_all)
    end if

! - Get index of behaviour law
    numlcEff = numlc+BEHInteg%behavPara%lawIndexOffset

! - Get parameters
    fami = BEHInteg%materPara%schemePara%fami
    kpg = BEHInteg%materPara%schemePara%kpg
    ksp = BEHInteg%materPara%schemePara%ksp
    jvMaterCode = BEHInteg%materPara%jvMaterCode

! --------------------------------------------------------------------------------------------------
    select case (numlcEff)
    case (1)
!     ELAS
        call lc0001(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (2)
!     VMIS_ISOT_XXX, VISC_ISOT_XXX
        call lc0002(fami, kpg, ksp, ndim, jvMaterCode, l_epsi_varc, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, typmod, ndsde, &
                    dsidep, codret)

    case (3)
!     VMIS_CINE_LINE, VMIS_ECMI_XXXX
        call lc0003(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, &
                    epsm, deps, sigm, vim, option, &
                    sigp, vip, typmod, nvi, &
                    dsidep, codret)

    case (4)
!     VMIS_CINX_CHAB/MEMO VISC_CINX_CHAB/MEMO,
        call lc0004(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, &
                    deps, sigm, vim, option, &
                    sigp, vip, typmod, &
                    nvi, dsidep, codret)

    case (7)
!     ENDO_ORTH_BETON
        call lc0007(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (8)
!     MAZARS
        call lc0008(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)
    case (9)
!     BETON_REGLE_PR
        call lc0009(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, sigm, vim, option, &
                    sigp, vip, typmod, &
                    nvi, dsidep, codret)

    case (15)
! ----- KIT_META
        call lc0015(BEHInteg, &
                    option, typmod, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, &
                    neps, epsm, deps, &
                    nsig, sigm, &
                    nvi, vim, &
                    sigp, vip, &
                    ndsde, dsidep, codret)

    case (16)
!     DRUCK_PRAGER
        call lc0016(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (17)
!     NORTON_HOFF
        call lc0017(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, sigm, vim, option, &
                    sigp, vip, typmod, &
                    nvi, dsidep, codret)

    case (18)
!     VISC_TAHERI
        call lc0018(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, vim, option, &
                    sigp, vip, typmod, &
                    nvi, dsidep, codret)

    case (19)
!     ELAS_HYPER
        call lc0019(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (21)
! ----- Beton_UMLV
        call lc0021(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (22)
! Cam-Clay
        call lc0022(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)
    case (23)
        call lc0023(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (24)
        call lc0024(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (25)
        call lc0025(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, &
                    epsm, deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    numlcEff, dsidep, codret)

    case (26)
        call lc0026(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (28)
        call lc0028(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (29)
        call lc0029(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (30)
        call lc0030(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, &
                    typmod, dsidep, &
                    codret)

    case (31)
        call lc0031(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (32)
        ! VISCOCHAB
        call lc0032(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (33)
        ! Hoek LAIGLE
        call lc0033(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, &
                    typmod, dsidep, codret)

    case (34)
        ! HUJEUX
        call lc0034(BEHInteg, &
                    fami, kpg, ksp, jvMaterCode, &
                    carcri, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (35)
        ! LETK
        call lc0035(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (36)
!     ENDO_ISOT_BETON
        call lc0036(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, &
                    typmod, ndsde, &
                    dsidep, codret)

    case (38)
!     SANS
        call lc0038(BEHInteg, neps, nsig, nvi, option, sigp, vip, ndsde, dsidep, codret)

    case (40)
!       DRUCKER_PRAGER_NA
        call lc0040(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, &
                    typmod, ndsde, &
                    dsidep, codret)

    case (42)
        call lc0042(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (50)
!     UMAT
        call lc0050(BEHInteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    jvMaterCode, compor, carcri, instam, instap, &
                    neps, epsm, deps, nsig, sigm, &
                    nvi, vim, option, &
                    sigp, vip, dsidep, codret)

    case (54)
        call lc0054(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (55)
        call lc0055(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)
    case (58)
!     MFRONT
        call lc0058(BEHInteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    jvMaterCode, compor, carcri, instam, instap, &
                    neps, epsm, deps, nsig, sigm, &
                    nvi, vim, option, &
                    sigp, vip, ndsde, dsidep, codret)

    case (59)
! - LKR
        call lc0059(BEHInteg, &
                    fami, kpg, ksp, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, &
                    typmod, dsidep, codret)

    case (60)
        call lc0060(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (62)
        call lc0062(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (75)
        call lc0075(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (76)
        call lc0076(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (77)
        call lc0077(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (78)
        call lc0078(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (79)
        call lc0079(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (101)
        call lc0101(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (102)
        call lc0102(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (120)
!     BETON_DOUBLE_DP
        call lc0120(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, l_epsi_varc, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (137)
!     MONOCRISTAL, POLYCRISTAL
        call lc0137(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, multComp, carcri, instam, instap, neps, &
                    epsm, deps, sigm, nvi, vim, option, &
                    sigp, vip, &
                    typmod, &
                    dsidep, codret)

    case (145)
!       BETON_RAG : nouvelle
        call lc0145(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (165)
!     FLUA_PORO_BETON
        call lc0165(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (166)
!     ENDO_PORO_BETON
        call lc0166(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (167)
!     FLUA_ENDO_PORO
        call lc0167(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (168)
!     RGI_BETON
        call lc0168(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (169)
!     RGI_BETON_BA
        call lc0169(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With SIMO_MIEHE
! --------------------------------------------------------------------------------------------------
!
    case (1002)
!     VMIS_ISOT_XXX, VISC_ISOT_XXX
        call lc1002(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, typmod, ndsde, &
                    dsidep, codret)

    case (1015)
! ----- KIT_META
        call lc1015(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (1037)
!     ROUSSELIER
        call lc1037(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (1137)
!     MONOCRISTAL, POLYCRISTAL
        call lc1137(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, multComp, carcri, instam, instap, neps, &
                    epsm, deps, sigm, nvi, vim, option, &
                    sigp, vip, &
                    typmod, &
                    dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With IMPLEX
! --------------------------------------------------------------------------------------------------
!
    case (2001)
!     ELAS
        call lc2001(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    neps, deps, nsig, sigm, option, &
                    sigp, nvi, vip, typmod, ndsde, &
                    dsidep, codret)

    case (2002)
!     VMIS_ISOT_XXX, VISC_ISOT_XXX
        call lc2002(fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, typmod, ndsde, &
                    dsidep, codret)

    case (2036)
!     ENDO_ISOT_BETON
        call lc2036(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, &
                    typmod, ndsde, &
                    dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With GRADSIGM
! --------------------------------------------------------------------------------------------------
!
    case (4047)
!     ENDO_HETEROGENE
        call lc4047(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - 1D constitutive law (5000+numldc)
! --------------------------------------------------------------------------------------------------
!
!
    case (5001)
        call lc5001(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (5002)
        call lc5002(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (5003)
        call lc5003(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (5024)
        call lc5024(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (5038)
!     SANS
        call lc5038(BEHInteg, neps, nsig, nvi, option, sigp, vip, ndsde, dsidep, codret)

    case (5076)
        call lc5076(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (5090)
        call lc5090(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (5091)
        call lc5091(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (5101)
        call lc5101(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (5102)
        call lc5102(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case (5152)
!     CABLE_GAINE
        call lc5152(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, &
                    ndsde, dsidep, codret)

! --------------------------------------------------------------------------------------------------
! - With GRADVARI
! --------------------------------------------------------------------------------------------------
!
    case (6036)
!     ENDO_ISOT_BETON
        call lc6036(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, &
                    typmod, ndsde, &
                    dsidep, codret)

    case (6046)
!     ENDO_SCALAIRE
        call lc6046(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, &
                    typmod, ndsde, &
                    dsidep, codret)
!
    case (6057)
!     ENDO_FISS_EXP
        call lc6057(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)
!
    case (6058)
!     MFRONT
        call lc6058(BEHInteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    jvMaterCode, compor, carcri, instam, instap, &
                    neps, epsm, deps, nsig, sigm, &
                    nvi, vim, option, &
                    sigp, vip, ndsde, dsidep, codret)
!
    case (6075)
!     GTN
        call lc6075(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)
!
    case (6076)
!     VMIS_ISOT_NL
        call lc6076(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With EJ_HYME/ELEMJOIN
! --------------------------------------------------------------------------------------------------
!
    case (7010)
!     CZM_EXP_REG
        call lc7010(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (7011)
!     CZM_LIN_REG
        call lc7011(fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (7013)
!     JOINT_BA
        call lc7013(fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (7045)
        call lc7045(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (7046)
        call lc7046(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (7047)
!     JOINT_MECA_ENDO
        call lc7047(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (7048)
        call lc7048(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, epsm, &
                    deps, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    dsidep, codret)

    case (7058)
!     MFRONT
        call lc7058(BEHInteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    jvMaterCode, compor, carcri, instam, instap, &
                    neps, epsm, deps, nsig, sigm, &
                    nvi, vim, option, &
                    sigp, vip, ndsde, dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - For KIT_DDI
! --------------------------------------------------------------------------------------------------
!
    case (8028)
        call lc8028(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, &
                    typmod, ndsde, dsidep, codret)

    case (8029)
        call lc8029(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, &
                    typmod, ndsde, dsidep, codret)

    case (8057)
        call lc8057(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, &
                    typmod, ndsde, dsidep, codret)
!
    case (8146)
        call lc8146(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, &
                    typmod, ndsde, dsidep, codret)
!
    case (8331)
        call lc8331(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, &
                    epsm, deps, nsig, sigm, nvi, vim, &
                    option, sigp, vip, &
                    typmod, ndsde, dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With INTERFACE
! --------------------------------------------------------------------------------------------------
!
    case (9040)
        call lc9040(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, &
                    ndsde, dsidep, codret)
    case (9041)
        call lc9041(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, &
                    ndsde, dsidep, codret)

    case (9043)
        call lc9043(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, &
                    ndsde, dsidep, codret)

    case (9049)
        call lc9049(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, &
                    ndsde, dsidep, codret)

    case (9051)
        call lc9051(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, &
                    ndsde, dsidep, codret)

    case (9056)
        call lc9056(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, &
                    ndsde, dsidep, codret)

!   MFRONT
    case (9058)
        call lc9058(BEHInteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    jvMaterCode, compor, carcri, instam, instap, &
                    neps, epsm, deps, nsig, sigm, &
                    nvi, vim, option, &
                    sigp, vip, ndsde, dsidep, codret)

    case (9077)
        call lc9077(BEHInteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, &
                    ndsde, dsidep, codret)

    case (9078)
        call lc9078(BEHinteg, &
                    fami, kpg, ksp, ndim, jvMaterCode, &
                    carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)
!
! --------------------------------------------------------------------------------------------------
! - With INSOLPI
! --------------------------------------------------------------------------------------------------
!
    case (9501)
        call lc9501(BEHinteg, fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)
    case (9502)
        call lc9502(BEHinteg, fami, kpg, ksp, ndim, jvMaterCode, &
                    compor, carcri, instam, instap, neps, epsm, &
                    deps, nsig, sigm, nvi, vim, option, &
                    sigp, vip, typmod, &
                    ndsde, dsidep, codret)

    case default
        call utmess('F', 'COMPOR1_43', si=numlcEff)
    end select
! --------------------------------------------------------------------------------------------------

! - For "old" prediction
    if (BEHInteg%behavPara%lPred .and. BEHInteg%behavPara%lSigm .and. &
        .not. BEHInteg%behavPara%lStrainMeca) then
        sigp = sigm
    end if

! - Viscous regularisation
    if (BEHInteg%behavPara%lReguVisc .and. codret .ne. LDC_ERROR_NCVG) then

        call lcvisc(fami, kpg, ksp, ndim, jvMaterCode, &
                    BEHInteg%behavPara%lSigm, BEHInteg%behavPara%lMatr, BEHInteg%behavPara%lVari, &
                    instam, instap, deps(1:ndimsi), &
                    vim(idx_regu_visc:idx_regu_visc+nvi_regu_visc-1), &
                    sigp(1:ndimsi), &
                    vip(idx_regu_visc:idx_regu_visc+nvi_regu_visc-1), &
                    dsidep(1:ndimsi, 1:ndimsi))
    end if
!
    ASSERT(codret .ge. LDC_ERROR_NONE .and. codret .le. LDC_ERROR_QUAL)
!
end subroutine
