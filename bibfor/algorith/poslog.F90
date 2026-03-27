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
! aslint: disable=W1504
!
subroutine poslog(lCorr, lMatr, lSigm, lVari, tlogPrev, &
                  tlogCurr, fPrev, lgpg, vip, ndim, &
                  fCurr, kpg, dtde, sigm, cplan, &
                  fami, mate, instp, angmas, gn, &
                  lamb, logl, sigmCurr, dsidep, pk2Prev, &
                  pk2Curr, codret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/d1macp.h"
#include "asterfort/deflg2.h"
#include "asterfort/deflg3.h"
#include "asterfort/lcdetf.h"
#include "asterfort/pk2sig.h"
#include "asterfort/symt46.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
    aster_logical, intent(in) :: lCorr, lMatr, lSigm, lVari
    aster_logical, intent(in) :: cplan
    real(kind=8), intent(in) :: tlogPrev(6)
    real(kind=8), intent(in) :: tlogCurr(6)
    real(kind=8), intent(in) :: fPrev(3, 3)
    real(kind=8), intent(in) :: fCurr(3, 3)
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: lgpg
    real(kind=8), intent(out) :: vip(lgpg)
    integer(kind=8), intent(in) :: kpg
    real(kind=8), intent(in) :: dtde(6, 6)
    real(kind=8), intent(in) :: sigm(2*ndim)
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: mate
    real(kind=8), intent(in) :: instp
    real(kind=8), intent(in) :: angmas(*)
    real(kind=8), intent(in) :: gn(3, 3)
    real(kind=8), intent(in) :: lamb(3)
    real(kind=8), intent(in) :: logl(3)
    real(kind=8), intent(out) :: sigmCurr(2*ndim)
    real(kind=8), intent(out) :: dsidep(6, 6)
    real(kind=8), intent(out) :: pk2Prev(6)
    real(kind=8), intent(out) :: pk2Curr(6)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!     BUT:  POST TRAITEMENT GRANDES DEFORMATIONS 2D ET 3D LOG
!     SUIVANT ARTICLE MIEHE APEL LAMBRECHT CMAME 2002
!     CONFIGURATION LAGRANGIENNE (PK2)
!
! --------------------------------------------------------------------------------------------------
!
! in  resi    : .true. si full_meca/raph_meca .false. si rigi_meca_tang
! in  lMatr    : .true. si full_meca/rigi_meca_tang
! in  tlogPrev      : contraintes associees aux def. logarithmiques en t-
! in  tlogCurr      : contraintes associees aux def. logarithmiques en t+
! in  fPrev      : gradient transformation en t-
! in  lgpg    : dimension du vecteur des var. internes pour 1 pt gauss
! var vip     : variables internes en t+
! in  ndim    : dimension de l'espace
! in  fCurr      : gradient transformation en t+
! in  pes     : operateur de transformation tlogPrev (ou tlogCurr) en pk2
! in  kpg       : numero du points de gauss
! in  dtde    : operateur tangent issu de nmcomp (6,6)
! in  sigm    : contrainte de cauchy en t-
! in  gn      : termes utiles au calcul de tl dans poslog
! in  feta    : termes utiles au calcul de tl dans poslog
! in  xi      : termes utiles au calcul de tl dans poslog
! in  me      : termes utiles au calcul de tl dans poslog
! out sigmCurr    : contraintes de cauchy en t+
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, j, kl
    real(kind=8) :: sig(6)
    real(kind=8) :: pes(6, 6), tp2(6), fr(3, 3), detf
    real(kind=8) :: tl(3, 3, 3, 3), tls(6, 6), epse(4), d1(4, 4)
    real(kind=8) :: feta(4), xi(3, 3), me(3, 3, 3, 3)
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8), dimension(6), parameter :: vrac2 = (/1.d0, 1.d0, 1.d0, rac2, rac2, rac2/)
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    pk2Prev = 0.d0
    pk2Curr = 0.d0
    codret = 0
    sig = 0.d0
    sigmCurr = 0.d0
!
! - Get gradient
    if (lCorr) then
        fr = fCurr
    else
        fr = fPrev
    end if
    call lcdetf(ndim, fr, detf)
    if (detf .le. r8prem()) then
        codret = 1
        goto 999
    end if
!
! CORRECTION POUR LES CONTRAINTES PLANES
! NE FONCTIONNE QUE SI DET(F_PLAS)=1  SOIT DEF. PLAS. INCOMPRESSIBLES
! CF. COMP. METHODES FOR PLASTICITY - DE SOUZA-NIETO, PERIC, OWEN p.603
    if (cplan) then
        epse = 0.d0
        if (lCorr) then
            call d1macp(fami, mate, instp, '+', kpg, &
                        1, angmas(1), d1)
            do i = 1, 4
                do j = 1, 4
                    epse(i) = epse(i)+d1(i, j)*tlogCurr(j)
                end do
            end do
            epse(3) = d1(1, 2)*(tlogCurr(1)+tlogCurr(2))
        else
            call d1macp(fami, mate, instp, '-', kpg, &
                        1, angmas(1), d1)
            do i = 1, 4
                do j = 1, 4
                    epse(i) = epse(i)+d1(i, j)*tlogPrev(j)
                end do
            end do
            epse(3) = d1(1, 2)*(tlogPrev(1)+tlogPrev(2))
        end if
        detf = exp(epse(1)+epse(2)+epse(3))
    end if
!
! - Tensor to change stress(log) to stress(PK2)
    call deflg2(gn, lamb, logl, pes, feta, &
                xi, me)
!
! - Tangent matrix
    if (lMatr) then
        dsidep = 0.d0
!        POUR LA RIGIDITE GEOMETRIQUE : CALCUL AVEC LES PK2
        if (lCorr) then
            tp2 = tlogCurr
        else
            sig(1:2*ndim) = sigm(1:2*ndim)
            call pk2sig(ndim, fPrev, detf, pk2Prev, sig, &
                        -1)
            do kl = 4, 2*ndim
                pk2Prev(kl) = pk2Prev(kl)*rac2
            end do
            tp2 = tlogPrev
!
        end if
!
        call deflg3(gn, feta, xi, me, tp2, tl)
        call symt46(tl, tls)
!
        dsidep = matmul(matmul(transpose(pes), dtde), pes)
!
        b_n = to_blas_int(36)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, tls, b_incx, dsidep, b_incy)
!
    end if
!
! - Compute Cauchy stress
    if (lSigm) then
        sigmCurr = 0.d0
        do i = 1, 6
            do j = 1, 6
                pk2Curr(i) = pk2Curr(i)+tlogCurr(j)*pes(j, i)
            end do
        end do
        call pk2sig(ndim, fCurr, detf, pk2Curr, sigmCurr, 1)
    end if
!
! - On stocke TP comme variable interne
    if (lVari) then
        !vip(lgpg-1:lgpg) = 0.d0
        b_n = to_blas_int(2*ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, tlogCurr/vrac2, b_incx, vip(lgpg-6+1), b_incy)
    end if
!
999 continue
!
end subroutine
