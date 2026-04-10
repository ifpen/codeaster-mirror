! --------------------------------------------------------------------
! Copyright (C) LAPACK
! Copyright (C) 2007 - 2026 - EDF - www.code-aster.org
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
! ===============================================================
! THIS LAPACK 2.0 ROUTINE IS DEPRECATED
! DO NOT USE IT : YOU SHOULD PREFER UP-TO-DATE LAPACK ROUTINE
!
! BUT DO NOT REMOVE IT :
! THE PRESENT ROUTINE IS MANDATORY FOR ARPACK LIBRARY
! WHICH STICKS TO LAPACK 2.0 VERSION
! ==============================================================
subroutine as_ztrevc(side, howmny, select, n, t, &
                     ldt, vl, ldvl, vr, ldvr, &
                     mm, m, work, rwork, info)
!  -- LAPACK ROUTINE (VERSION 2.0) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     SEPTEMBER 30, 1994
!
!  PURPOSE
!  =======
!
!  ZTREVC COMPUTES SOME OR ALL OF THE RIGHT AND/OR LEFT EIGENVECTORS OF
!  A COMPLEX UPPER TRIANGULAR MATRIX T.
!
!  THE RIGHT EIGENVECTOR X AND THE LEFT EIGENVECTOR Y OF T CORRESPONDING
!  TO AN EIGENVALUE W ARE DEFINED BY:
!
!               T*X = W*X,     Y'*T = W*Y'
!
!  WHERE Y' DENOTES THE CONJUGATE TRANSPOSE OF THE VECTOR Y.
!
!  IF ALL EIGENVECTORS ARE REQUESTED, THE ROUTINE MAY EITHER RETURN THE
!  MATRICES X AND/OR Y OF RIGHT OR LEFT EIGENVECTORS OF T, OR THE
!  PRODUCTS Q*X AND/OR Q*Y, WHERE Q IS AN INPUT UNITARY
!  MATRIX. IF T WAS OBTAINED FROM THE SCHUR FACTORIZATION OF AN
!  ORIGINAL MATRIX A = Q*T*Q', THEN Q*X AND Q*Y ARE THE MATRICES OF
!  RIGHT OR LEFT EIGENVECTORS OF A.
!
!  ARGUMENTS
!  =========
!
!  SIDE    (INPUT) CHARACTER*1
!          = 'R':  COMPUTE RIGHT EIGENVECTORS ONLY;
!          = 'L':  COMPUTE LEFT EIGENVECTORS ONLY;
!          = 'B':  COMPUTE BOTH RIGHT AND LEFT EIGENVECTORS.
!
!  HOWMNY  (INPUT) CHARACTER*1
!          = 'A':  COMPUTE ALL RIGHT AND/OR LEFT EIGENVECTORS;
!          = 'B':  COMPUTE ALL RIGHT AND/OR LEFT EIGENVECTORS,
!                  AND BACKTRANSFORM THEM USING THE INPUT MATRICES
!                  SUPPLIED IN VR AND/OR VL;
!          = 'S':  COMPUTE SELECTED RIGHT AND/OR LEFT EIGENVECTORS,
!                  SPECIFIED BY THE LOGICAL ARRAY SELECT.
!
!  SELECT  (INPUT) LOGICAL ARRAY, DIMENSION (N)
!          IF HOWMNY = 'S', SELECT SPECIFIES THE EIGENVECTORS TO BE
!          COMPUTED.
!          IF HOWMNY = 'A' OR 'B', SELECT IS NOT REFERENCED.
!          TO SELECT THE EIGENVECTOR CORRESPONDING TO THE J-TH
!          EIGENVALUE, SELECT(J) MUST BE SET TO .TRUE..
!
!  N       (INPUT) INTEGER
!          THE ORDER OF THE MATRIX T. N >= 0.
!
!  T       (INPUT/OUTPUT) COMPLEX*16 ARRAY, DIMENSION (LDT,N)
!          THE UPPER TRIANGULAR MATRIX T.  T IS MODIFIED, BUT RESTORED
!          ON EXIT.
!
!  LDT     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY T. LDT >= MAX(1,N).
!
!  VL      (INPUT/OUTPUT) COMPLEX*16 ARRAY, DIMENSION (LDVL,MM)
!          ON ENTRY, IF SIDE = 'L' OR 'B' AND HOWMNY = 'B', VL MUST
!          CONTAIN AN N-BY-N MATRIX Q (USUALLY THE UNITARY MATRIX Q OF
!          SCHUR VECTORS RETURNED BY ZHSEQR).
!          ON EXIT, IF SIDE = 'L' OR 'B', VL CONTAINS:
!          IF HOWMNY = 'A', THE MATRIX Y OF LEFT EIGENVECTORS OF T;
!          IF HOWMNY = 'B', THE MATRIX Q*Y;
!          IF HOWMNY = 'S', THE LEFT EIGENVECTORS OF T SPECIFIED BY
!                           SELECT, STORED CONSECUTIVELY IN THE COLUMNS
!                           OF VL, IN THE SAME ORDER AS THEIR
!                           EIGENVALUES.
!          IF SIDE = 'R', VL IS NOT REFERENCED.
!
!  LDVL    (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY VL.  LDVL >= MAX(1,N) IF
!          SIDE = 'L' OR 'B'; LDVL >= 1 OTHERWISE.
!
!  VR      (INPUT/OUTPUT) COMPLEX*16 ARRAY, DIMENSION (LDVR,MM)
!          ON ENTRY, IF SIDE = 'R' OR 'B' AND HOWMNY = 'B', VR MUST
!          CONTAIN AN N-BY-N MATRIX Q (USUALLY THE UNITARY MATRIX Q OF
!          SCHUR VECTORS RETURNED BY ZHSEQR).
!          ON EXIT, IF SIDE = 'R' OR 'B', VR CONTAINS:
!          IF HOWMNY = 'A', THE MATRIX X OF RIGHT EIGENVECTORS OF T;
!          IF HOWMNY = 'B', THE MATRIX Q*X;
!          IF HOWMNY = 'S', THE RIGHT EIGENVECTORS OF T SPECIFIED BY
!                           SELECT, STORED CONSECUTIVELY IN THE COLUMNS
!                           OF VR, IN THE SAME ORDER AS THEIR
!                           EIGENVALUES.
!          IF SIDE = 'L', VR IS NOT REFERENCED.
!
!  LDVR    (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY VR.  LDVR >= MAX(1,N) IF
!           SIDE = 'R' OR 'B'; LDVR >= 1 OTHERWISE.
!
!  MM      (INPUT) INTEGER
!          THE NUMBER OF COLUMNS IN THE ARRAYS VL AND/OR VR. MM >= M.
!
!  M       (OUTPUT) INTEGER
!          THE NUMBER OF COLUMNS IN THE ARRAYS VL AND/OR VR ACTUALLY
!          USED TO STORE THE EIGENVECTORS.  IF HOWMNY = 'A' OR 'B', M
!          IS SET TO N.  EACH SELECTED EIGENVECTOR OCCUPIES ONE
!          COLUMN.
!
!  WORK    (WORKSPACE) COMPLEX*16 ARRAY, DIMENSION (2*N)
!
!  RWORK   (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (N)
!
!  INFO    (OUTPUT) INTEGER
!          = 0:  SUCCESSFUL EXIT
!          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
!
!  FURTHER DETAILS
!  ===============
!
!  THE ALGORITHM USED IN THIS PROGRAM IS BASICALLY BACKWARD (FORWARD)
!  SUBSTITUTION, WITH SCALING TO MAKE THE THE CODE ROBUST AGAINST
!  POSSIBLE OVERFLOW.
!
!  EACH EIGENVECTOR IS NORMALIZED SO THAT THE ELEMENT OF LARGEST
!  MAGNITUDE HAS MAGNITUDE 1; HERE THE MAGNITUDE OF A COMPLEX NUMBER
!  (X,Y) IS TAKEN TO BE |X| + |Y|.
!
!  =====================================================================
!-----------------------------------------------------------------------
! ASTER INFORMATION
! 14/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            REMPLACEMENT DE 1 RETURN PAR GOTO 1000,
!            IMPLICIT NONE.
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!
!     .. SCALAR ARGUMENTS ..
#include "asterf_types.h"
#include "asterc/isbaem.h"
#include "asterc/matfpe.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/xerbla.h"
#include "blas/dzasum.h"
#include "blas/izamax.h"
#include "blas/lsame.h"
#include "blas/zcopy.h"
#include "blas/zdscal.h"
#include "blas/zgemv.h"
#include "blas/zlatrs.h"
    character(len=1) :: howmny, side
    integer(kind=8) :: info, ldt, ldvl, ldvr, m, mm, n
!     ..
!     .. ARRAY ARGUMENTS ..
    aster_logical :: select(*)
    real(kind=8) :: rwork(*)
    complex(kind=8) :: t(ldt, *), vl(ldvl, *), vr(ldvr, *), work(*)
!     ..
!     .. PARAMETERS ..
    real(kind=8) :: zero, one
    parameter(zero=0.0d+0, one=1.0d+0)
    complex(kind=8) :: cmzero, cmone
    parameter(cmzero=(0.0d+0, 0.0d+0),&
     &                   cmone=(1.0d+0, 0.0d+0))
!     ..
!     .. LOCAL SCALARS ..
    aster_logical :: allv, bothv, leftv, over, rightv, somev
    integer(kind=8) :: i, ii, is, j, k, ki
    integer(kind=4) :: info4
    real(kind=8) :: remax, scale, smin, smlnum, ulp, unfl
    blas_int :: b_incx, b_incy, b_n
    blas_int :: b_lda, b_m
!     ..
!     .. EXTERNAL FUNCTIONS ..
!     ..
!     .. STATEMENT FUNCTIONS ..
!     ..
!     .. STATEMENT FUNCTION DEFINITIONS ..
#define cabs1( cdum ) abs( dble( cdum ) ) + abs( dimag( cdum ) )
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
    call matfpe(-1)
!
!     DECODE AND TEST THE INPUT PARAMETERS
!
    bothv = lsame(side, 'B')
    rightv = lsame(side, 'R') .or. bothv
    leftv = lsame(side, 'L') .or. bothv
!
    allv = lsame(howmny, 'A')
    over = lsame(howmny, 'B') .or. lsame(howmny, 'O')
    somev = lsame(howmny, 'S')
!
!     SET M TO THE NUMBER OF COLUMNS REQUIRED TO STORE THE SELECTED
!     EIGENVECTORS.
!
    if (somev) then
        m = 0
        do j = 1, n
            if (select(j)) m = m+1
        end do
    else
        m = n
    end if
!
    info4 = 0
    if (.not. rightv .and. .not. leftv) then
        info4 = -1
    else if (.not. allv .and. .not. over .and. .not. somev) then
        info4 = -2
    else if (n .lt. 0) then
        info4 = -4
    else if (ldt .lt. max(1, n)) then
        info4 = -6
    else if (ldvl .lt. 1 .or. (leftv .and. ldvl .lt. n)) then
        info4 = -8
    else if (ldvr .lt. 1 .or. (rightv .and. ldvr .lt. n)) then
        info4 = -10
    else if (mm .lt. m) then
        info4 = -11
    end if
    if (info4 .ne. 0) then
        info = info4
        call xerbla('ZTREVC', -info)
        goto 1000
    end if
!
!     QUICK RETURN IF POSSIBLE.
!
    if (n .eq. 0) goto 1000
!
!     SET THE CONSTANTS TO CONTROL OVERFLOW.
!
    unfl = r8miem()
! DUE TO CRS512      OVFL = ONE / UNFL
    ulp = r8prem()*0.5d0*isbaem()
    smlnum = unfl*(n/ulp)
!
!     STORE THE DIAGONAL ELEMENTS OF T IN WORKING ARRAY WORK.
!
    do i = 1, n
        work(i+n) = t(i, i)
    end do
!
!     COMPUTE 1-NORM OF EACH COLUMN OF STRICTLY UPPER TRIANGULAR
!     PART OF T TO CONTROL OVERFLOW IN TRIANGULAR SOLVER.
!
    rwork(1) = zero
    do j = 2, n
        b_n = to_blas_int(j-1)
        b_incx = to_blas_int(1)
        rwork(j) = dzasum(b_n, t(1, j), b_incx)
    end do
!
    if (rightv) then
!
!        COMPUTE RIGHT EIGENVECTORS.
!
        is = m
        do ki = n, 1, -1
!
            if (somev) then
                if (.not. select(ki)) goto 80
            end if
            smin = max(ulp*(cabs1(t(ki, ki))), smlnum)
!
            work(1) = cmone
!
!           FORM RIGHT-HAND SIDE.
!
            do k = 1, ki-1
                work(k) = -t(k, ki)
            end do
!
!           SOLVE THE TRIANGULAR SYSTEM:
!              (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.
!
            do k = 1, ki-1
                t(k, k) = t(k, k)-t(ki, ki)
                if (cabs1(t(k, k)) .lt. smin) t(k, k) = smin
            end do
!
            if (ki .gt. 1) then
                b_lda = to_blas_int(ldt)
                b_n = to_blas_int(ki-1)
                call zlatrs('U', 'N', 'N', 'Y', b_n, &
                            t, b_lda, work(1), scale, rwork, &
                            info4)
                work(ki) = scale
            end if
!
!           COPY THE VECTOR X OR Q*X TO VR AND NORMALIZE.
!
            if (.not. over) then
                b_n = to_blas_int(ki)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call zcopy(b_n, work(1), b_incx, vr(1, is), b_incy)
!
                b_n = to_blas_int(ki)
                b_incx = to_blas_int(1)
                ii = izamax(b_n, vr(1, is), b_incx)
                remax = one/cabs1(vr(ii, is))
                b_n = to_blas_int(ki)
                b_incx = to_blas_int(1)
                call zdscal(b_n, remax, vr(1, is), b_incx)
!
                do k = ki+1, n
                    vr(k, is) = cmzero
                end do
            else
                if (ki .gt. 1) then
                    b_lda = to_blas_int(ldvr)
                    b_m = to_blas_int(n)
                    b_n = to_blas_int(ki-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call zgemv('N', b_m, b_n, cmone, vr, &
                               b_lda, work(1), b_incx, dcmplx(scale), vr(1, ki), &
                               b_incy)
                end if
!
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                ii = izamax(b_n, vr(1, ki), b_incx)
                remax = one/cabs1(vr(ii, ki))
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                call zdscal(b_n, remax, vr(1, ki), b_incx)
            end if
!
!           SET BACK THE ORIGINAL DIAGONAL ELEMENTS OF T.
!
            do k = 1, ki-1
                t(k, k) = work(k+n)
            end do
!
            is = is-1
80          continue
        end do
    end if
!
    if (leftv) then
!
!        COMPUTE LEFT EIGENVECTORS.
!
        is = 1
        do ki = 1, n
!
            if (somev) then
                if (.not. select(ki)) goto 130
            end if
            smin = max(ulp*(cabs1(t(ki, ki))), smlnum)
!
            work(n) = cmone
!
!           FORM RIGHT-HAND SIDE.
!
            do k = ki+1, n
                work(k) = -dconjg(t(ki, k))
            end do
!
!           SOLVE THE TRIANGULAR SYSTEM:
!              (T(KI+1:N,KI+1:N) - T(KI,KI))'*X = SCALE*WORK.
!
            do k = ki+1, n
                t(k, k) = t(k, k)-t(ki, ki)
                if (cabs1(t(k, k)) .lt. smin) t(k, k) = smin
            end do
!
            if (ki .lt. n) then
                b_lda = to_blas_int(ldt)
                b_n = to_blas_int(n-ki)
                call zlatrs('U', 'C', 'N', 'Y', b_n, &
                            t(ki+1, ki+1), b_lda, work(ki+1), scale, rwork, &
                            info4)
                work(ki) = scale
            end if
!
!           COPY THE VECTOR X OR Q*X TO VL AND NORMALIZE.
!
            if (.not. over) then
                b_n = to_blas_int(n-ki+1)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call zcopy(b_n, work(ki), b_incx, vl(ki, is), b_incy)
!
                b_n = to_blas_int(n-ki+1)
                b_incx = to_blas_int(1)
                ii = izamax(b_n, vl(ki, is), b_incx)+ki-1
                remax = one/cabs1(vl(ii, is))
                b_n = to_blas_int(n-ki+1)
                b_incx = to_blas_int(1)
                call zdscal(b_n, remax, vl(ki, is), b_incx)
!
                do k = 1, ki-1
                    vl(k, is) = cmzero
                end do
            else
                if (ki .lt. n) then
                    b_lda = to_blas_int(ldvl)
                    b_m = to_blas_int(n)
                    b_n = to_blas_int(n-ki)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call zgemv('N', b_m, b_n, cmone, vl(1, ki+1), &
                               b_lda, work(ki+1), b_incx, dcmplx(scale), vl(1, ki), &
                               b_incy)
                end if
!
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                ii = izamax(b_n, vl(1, ki), b_incx)
                remax = one/cabs1(vl(ii, ki))
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                call zdscal(b_n, remax, vl(1, ki), b_incx)
            end if
!
!           SET BACK THE ORIGINAL DIAGONAL ELEMENTS OF T.
!
            do k = ki+1, n
                t(k, k) = work(k+n)
            end do
!
            is = is+1
130         continue
        end do
    end if
!
1000 continue
    info = info4
    call matfpe(1)
!
!     END OF ZTREVC
!
end subroutine
