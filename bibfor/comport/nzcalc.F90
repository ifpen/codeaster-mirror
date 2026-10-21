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
! aslint: disable=W0413
!
subroutine nzcalc(carcri, nbPhase, phase, zSumCold, &
                  fmix, seuil, timeIncr, trans, &
                  rprim, deuxmu, eta, unsurn, &
                  dp, codret)
!
    implicit none
!
#include "asterfort/Behaviour_type.h"
#include "asterfort/nzfpri.h"
#include "asterfort/zeroco.h"
!
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    integer(kind=8), intent(in) :: nbPhase
    real(kind=8), intent(in) :: phase(nbPhase), zSumCold
    real(kind=8), intent(in) :: fmix, seuil, timeIncr, trans
    real(kind=8), intent(in) :: rprim, deuxmu, eta(nbPhase), unsurn(nbPhase)
    real(kind=8), intent(out) :: dp
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Comportment - Metallurgy
!
! CALCUL DE DP PAR RESOLUTION DE L'EQUATION SCALAIRE FPRIM=0
!
! --------------------------------------------------------------------------------------------------
!
! FPRIM=F1-F(I) TEL QUE F1 : DROITE DECROISSANTE
!                       F(I):A*(DP**UNSURN(I))
!
! In  nbPhase          : total number of phase (cold and hot)
! In  phase            : phases
! In  zSumCold           : sum of cold phases
! In  fmix             : mixing function
! IN CRIT : CRITERES DE CONVERGENCE LOCAUX
! IN FPRIM   : FONCTION SEUIL
! IN DT    :TP-TM
!
! OUT DP
!     IRET   CODE RETOUR DE LA RESOLUTION DE L'EQUATION SCALAIRE
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ECHEC
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: fprim, fplas, fdevi, dpu(5), fp(5), fd(5), eps, dpnew, fpnew
    real(kind=8) :: spetra, r0, coefa(5), dpplas, dpmin, dpmax, x(4), y(4)
    integer(kind=8) :: iter, i_phase
!
! --------------------------------------------------------------------------------------------------
!
    eps = 1.d-6
    spetra = 1.5d0*deuxmu*trans+1.d0
    r0 = 1.5d0*deuxmu+rprim*spetra

! - Plastic try
    dpplas = seuil/r0
    dp = dpplas
    call nzfpri(deuxmu, trans, rprim, seuil, &
                nbPhase, phase, zSumCold, &
                fmix, eta, unsurn, &
                timeIncr, dp, &
                fplas, fp, fd, &
                fprim, fdevi)
    if (abs(fprim)/(1.d0+seuil) .lt. carcri(3)) then
        goto 999
    end if
    dpmax = dpplas
    dpmin = 0.d0
!
! - Get bounds
!
    do i_phase = 1, nbPhase
        if (phase(i_phase) .gt. eps) then
            dp = dpplas
            call nzfpri(deuxmu, trans, rprim, seuil, &
                        nbPhase, phase, zSumCold, &
                        fmix, eta, unsurn, &
                        timeIncr, dp, &
                        fplas, fp, fd, &
                        fprim, fdevi)
            if (abs(fprim)/(1.d0+seuil) .lt. carcri(3)) then
                goto 999
            end if
            if (abs(fp(i_phase))/(1.d0+seuil) .lt. carcri(3)) then
                goto 99
            end if
            dp = 0.d0
            coefa(i_phase) = (eta(i_phase)*spetra)/timeIncr**unsurn(i_phase)
            do iter = 1, nint(carcri(1))
                call nzfpri(deuxmu, trans, rprim, seuil, &
                            nbPhase, phase, zSumCold, &
                            fmix, eta, unsurn, &
                            timeIncr, dp, &
                            fplas, fp, fd, &
                            fprim, fdevi)
                if (abs(fprim)/(1+seuil) .lt. carcri(3)) then
                    goto 999
                end if
                if (abs(fp(i_phase))/(1+seuil) .lt. carcri(3)) then
                    goto 99
                end if
                dp = (fplas/coefa(i_phase))**(1.d0/unsurn(i_phase))
                if (dp .gt. dpplas) then
                    dp = dpplas
                end if
                call nzfpri(deuxmu, trans, rprim, seuil, &
                            nbPhase, phase, zSumCold, &
                            fmix, eta, unsurn, &
                            timeIncr, dp, &
                            fplas, fp, fd, &
                            fprim, fdevi)
                if (abs(fprim)/(1+seuil) .lt. carcri(3)) then
                    goto 999
                end if
                if (abs(fp(i_phase))/(1+seuil) .lt. carcri(3)) then
                    goto 99
                end if
                dp = dp-fp(i_phase)/fd(i_phase)
            end do
            codret = 1
            goto 999
99          continue
            dpu(i_phase) = dp
            if (dpmin .eq. 0.d0) then
                dpmin = dpu(i_phase)
                dpmax = dpu(i_phase)
            else
                dpmin = min(dpmin, dpu(i_phase))
                dpmax = max(dpmax, dpu(i_phase))
            end if
        end if
    end do
!
!
!  RESOLUTION DE FPRIM=0 - METHODE SECANTE AVEC NEWTON SI BESOIN
!
!------EXAMEN DE LA SOLUTION DP=DPMIN
    dp = dpmin
    call nzfpri(deuxmu, trans, rprim, seuil, &
                nbPhase, phase, zSumCold, &
                fmix, eta, unsurn, &
                timeIncr, dp, &
                fplas, fp, fd, &
                fprim, fdevi)

    if (fprim .lt. 0.d0) then
        codret = 1
        goto 999
    end if
    x(2) = dp
    y(2) = fprim
!
!------EXAMEN DE LA SOLUTION DP=DPMAX
!
    dp = dpmax
    call nzfpri(deuxmu, trans, rprim, seuil, &
                nbPhase, phase, zSumCold, &
                fmix, eta, unsurn, &
                timeIncr, dp, &
                fplas, fp, fd, &
                fprim, fdevi)
!
    x(1) = dp
    y(1) = fprim
!
!------CALCUL DE DP : EQUATION SCALAIRE FPRIM = 0 AVEC DPMIN< DP < DPMAX
    x(3) = x(1)
    y(3) = y(1)
    x(4) = x(2)
    y(4) = y(2)
    do iter = 1, int(carcri(1))
        call zeroco(x, y)
        dp = x(4)
        call nzfpri(deuxmu, trans, rprim, seuil, &
                    nbPhase, phase, zSumCold, &
                    fmix, eta, unsurn, &
                    timeIncr, dp, &
                    fplas, fp, fd, &
                    fprim, fdevi)
        if (abs(fprim)/(1.d0+seuil) .lt. carcri(3)) then
            goto 999
        end if
        dpnew = dp-fprim/fdevi
        if ((dpnew .ge. dpmin) .and. (dpnew .le. dpmax)) then
            call nzfpri(deuxmu, trans, rprim, seuil, &
                        nbPhase, phase, zSumCold, &
                        fmix, eta, unsurn, &
                        timeIncr, dpnew, &
                        fplas, fp, fd, &
                        fpnew, fdevi)
            if (abs(fpnew)/(1.d0+seuil) .lt. carcri(3)) then
                dp = dpnew
                goto 999
            end if
            if (abs(fpnew)/(1.d0+seuil) .lt. abs(fprim)/(1.d0+seuil)) then
                dp = dpnew
                fprim = fpnew
            end if
        end if
        y(4) = fprim
        x(4) = dp
    end do
    codret = 1
!
999 continue
end subroutine
