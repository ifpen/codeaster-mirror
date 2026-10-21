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
subroutine meta_vpta_coef(metaRela, metaGlob, &
                          lgpg, fami, kpg, jvMaterCode, &
                          l_temp, temp, metaType, nbPhase, phasPrev, &
                          phasCurr, zColdCurr, young, deuxmu, coef, &
                          trans)
!
    use MetallurgyMeca_module, only: metaGetType, metaGetParaElas, metaGetMechanism, &
                                     metaGetParaMixture, metaGetParaAnneal, metaGetParaVisc, &
                                     metaGetParaPlasTransf, metaGetParaHardLine, metaGetParaHardTrac
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: metaRela, metaGlob
    integer(kind=8), intent(in) :: lgpg
    character(len=4), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: jvMaterCode
    aster_logical, intent(in) :: l_temp
    real(kind=8), intent(in) :: temp
    integer(kind=8), intent(in) :: metaType
    integer(kind=8), intent(in) :: nbPhase
    real(kind=8), intent(in) :: phasPrev(META_MECA_NBPHASE_MAXI), phasCurr(META_MECA_NBPHASE_MAXI)
    real(kind=8), intent(in) :: zColdCurr
    real(kind=8), intent(in) :: young, deuxmu
    real(kind=8), intent(out) :: coef
    real(kind=8), intent(out) :: trans
!
! --------------------------------------------------------------------------------------------------
!
! Metallurgy - Comportment
!
! Compute coefficient for second member
! Effect of command variables phases variation on transformation plasticity
!
! --------------------------------------------------------------------------------------------------
!
! In  metaRela      : behaviour for each phase
! In  metaGlob      : global behaviour
! In  lgpg          : length of integration point
! In  fami          : integration point type
! In  kpg           : integration point number
! In  jvMaterCode   : coded material address
! In  l_temp        : .true. if temperature command variable is affected
! In  temp          : temperature
! In  metaType      : type of metallurgy
! In  nbPhase       : total number of phases (cold and hot)
! In  phasPrev      : previous phaszs
! In  phasCurr      : current phases
! In  zColdCurr     : sum of cold phases
! In  young         : Young modulus
! In  deuxmu        : e/(1.d0+nu)
! Out coef          : coefficient for command variable second member
! Out trans         : transformation for command variable second member
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: jvVari
    integer(kind=8) :: iPhase, iPhasesCold, nbPhaseCold
    real(kind=8) :: epseq(META_MECA_NBPHASE_MAXI), h0(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: kpt(META_MECA_NBPHASE_MAXI-1), fpt(META_MECA_NBPHASE_MAXI-1)
    real(kind=8) :: eta(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: n(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: unsurn(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: c(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: m(META_MECA_NBPHASE_MAXI)
    real(kind=8) :: rprim, phaseIncr(META_MECA_NBPHASE_MAXI), fmix, coefHard
    aster_logical :: l_visc, l_elas, l_plas_tran, l_hard_line
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbPhase .le. META_MECA_NBPHASE_MAXI)
    nbPhaseCold = nbPhase-1
    trans = 0.d0
    coef = 1.d0
    fpt = 0.d0
    kpt = 0.d0
    eta = 0.d0
    n = 0.d0
    unsurn = 0.d0
    c = 0.d0
    m = 0.d0
    h0 = 0.d0
    phaseIncr = 0.d0
    phaseIncr = (phasCurr-phasPrev)

! - Cumulated cumulative plastic strain
    call jevech('PVARIPR', 'L', jvVari)
    epseq = 0.d0
    do iPhase = 1, nbPhase
        epseq(iPhase) = zr(jvVari+lgpg*(kpg-1)-1+iPhase)
    end do

! - Is elastic ?
    if (metaType .eq. META_STEEL .or. metaType .eq. META_STEEL_TEMPER) then
        l_elas = zr(jvVari+lgpg*(kpg-1)+(nbPhase+1)) .lt. 0.5d0
    elseif (metaType .eq. META_ZIRC) then
        l_elas = zr(jvVari+lgpg*(kpg-1)+nbPhase) .lt. 0.5d0
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get active mechanisms in behaviour
    call metaGetMechanism(metaRela, metaGlob, &
                          lVisc_=l_visc, &
                          lHardIsotLine_=l_hard_line, lPlasTran_=l_plas_tran)

! - Transformation plasticity parameters
    if (l_plas_tran) then
        call metaGetParaPlasTransf('+', fami, kpg, ksp, jvMaterCode, &
                                   metaType, nbPhase, nbPhaseCold, &
                                   phaseIncr, zColdCurr, &
                                   kpt, fpt)
    end if

! - Visco-plasticity parameters
    if (l_visc) then
        call metaGetParaVisc('+', fami, kpg, ksp, jvMaterCode, &
                             metaType, nbPhase, &
                             eta, n, unsurn, &
                             c, m)
    end if

! - Compute Sum(iphase) [kpt * fpt] on cold phases
    trans = 0.d0
    do iPhasesCold = 1, nbPhaseCold
        if (phaseIncr(iPhasesCold) .gt. 0) then
            trans = trans+kpt(iPhasesCold)*fpt(iPhasesCold)*phaseIncr(iPhasesCold)
        end if
    end do

! - Compute coefficient
    if (l_elas .or. l_visc) then
        coef = 1.d0
    else
! ----- Mixing law: yield
        call metaGetParaMixture('+', fami, kpg, ksp, jvMaterCode, &
                                l_visc, metaType, nbPhase, zColdCurr, &
                                fmix)

! ----- Get point on hardening curve
        if (l_hard_line) then
            coefHard = 1.d0
            call metaGetParaHardLine('+', fami, kpg, ksp, jvMaterCode, &
                                     metaType, nbPhase, &
                                     young, coefHard, h0)
        else
            call metaGetParaHardTrac(jvMaterCode, metaType, nbPhase, &
                                     l_temp, temp, &
                                     epseq, h0)
        end if

! ----- Compute coefficient
        rprim = 0.d0
        if (zColdCurr .gt. 0.d0) then
            do iPhasesCold = 1, nbPhaseCold
                rprim = rprim+phasCurr(iPhasesCold)*h0(iPhasesCold)
            end do
            rprim = rprim/zColdCurr
        else
            rprim = 0.d0
        end if
        rprim = (1.d0-fmix)*h0(nbPhase)+fmix*rprim
        coef = 1.d0-(1.5d0*deuxmu)/(1.5d0*deuxmu+rprim)
    end if
!
end subroutine
