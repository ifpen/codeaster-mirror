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
subroutine as_dlrsyl(trana, tranb, isgn, m, n, &
                     a, lda, b, ldb, c, &
                     ldc, scale, info)
!
!     SUBROUTINE LAPACK RESOLVANT L'EQUATION DE SYLVESTER.
!-----------------------------------------------------------------------
!  -- LAPACK ROUTINE (VERSION 2.0) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     MARCH 31, 1993
!
!  PURPOSE
!  =======
!
!  DLRSYL SOLVES THE REAL SYLVESTER MATRIX EQUATION:
!
!     OP(A)*X + X*OP(B) = SCALE*C OR
!     OP(A)*X - X*OP(B) = SCALE*C,
!
!  WHERE OP(A) = A OR A**T, AND  A AND B ARE BOTH UPPER QUASI-
!  TRIANGULAR. A IS M-BY-M AND B IS N-BY-N, THE RIGHT HAND SIDE C AND
!  THE SOLUTION X ARE M-BY-N, AND SCALE IS AN OUTPUT SCALE FACTOR, SET
!  <= 1 TO AVOID OVERFLOW IN X.
!
!  A AND B MUST BE IN SCHUR CANONICAL FORM (AS RETURNED BY DHSEQR), THAT
!  IS, BLOCK UPPER TRIANGULAR WITH 1-BY-1 AND 2-BY-2 DIAGONAL BLOCKS,
!  EACH 2-BY-2 DIAGONAL BLOCK HAS ITS DIAGONAL ELEMENTS EQUAL AND ITS
!  OFF-DIAGONAL ELEMENTS OF OPPOSITE SIGN.
!
!  ARGUMENTS
!  =========
!
!  TRANA   (INPUT) CHARACTER*1
!          SPECIFIES THE OPTION OP(A):
!          = 'N': OP(A) = A    (NO TRANSPOSE)
!          = 'T': OP(A) = A**T (TRANSPOSE)
!          = 'C': OP(A) = A**H (CONJUGATE TRANSPOSE = TRANSPOSE)
!
!  TRANB   (INPUT) CHARACTER*1
!          SPECIFIES THE OPTION OP(B):
!          = 'N': OP(B) = B    (NO TRANSPOSE)
!          = 'T': OP(B) = B**T (TRANSPOSE)
!          = 'C': OP(B) = B**H (CONJUGATE TRANSPOSE = TRANSPOSE)
!
!  ISGN    (INPUT) INTEGER
!          SPECIFIES THE SIGN IN THE EQUATION:
!          = +1: SOLVE OP(A)*X + X*OP(B) = SCALE*C
!          = -1: SOLVE OP(A)*X - X*OP(B) = SCALE*C
!
!  M       (INPUT) INTEGER
!          THE ORDER OF THE MATRIX A, AND THE NUMBER OF ROWS IN THE
!          MATRICES X AND C. M >= 0.
!
!  N       (INPUT) INTEGER
!          THE ORDER OF THE MATRIX B, AND THE NUMBER OF COLUMNS IN THE
!          MATRICES X AND C. N >= 0.
!
!  A       (INPUT) REAL*8 ARRAY, DIMENSION (LDA,M)
!          THE UPPER QUASI-TRIANGULAR MATRIX A, IN SCHUR CANONICAL FORM.
!
!  LDA     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY A. LDA >= MAX(1,M).
!
!  B       (INPUT) REAL*8 ARRAY, DIMENSION (LDB,N)
!          THE UPPER QUASI-TRIANGULAR MATRIX B, IN SCHUR CANONICAL FORM.
!
!  LDB     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY B. LDB >= MAX(1,N).
!
!  C       (INPUT/OUTPUT) REAL*8 ARRAY, DIMENSION (LDC,N)
!          ON ENTRY, THE M-BY-N RIGHT HAND SIDE MATRIX C.
!          ON EXIT, C IS OVERWRITTEN BY THE SOLUTION MATRIX X.
!
!  LDC     (INPUT) INTEGER
!          THE LEADING DIMENSION OF THE ARRAY C. LDC >= MAX(1,M)
!
!  SCALE   (OUTPUT) REAL*8
!          THE SCALE FACTOR, SCALE, SET <= 1 TO AVOID OVERFLOW IN X.
!
!  INFO    (OUTPUT) INTEGER
!          = 0: SUCCESSFUL EXIT
!          < 0: IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
!          = 1: A AND B HAVE COMMON OR VERY CLOSE EIGENVALUES, PERTURBED
!               VALUES WERE USED TO SOLVE THE EQUATION (BUT THE MATRICES
!               A AND B ARE UNCHANGED).
!
! ASTER INFORMATION
! 07/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!         REMPLACEMENT DE DLAMCH ET DLABAD PAR R8PREM,R8MIEM ET ISBAEM,
!            REMPLACEMENT DE 2 RETURN PAR UN GOTO 1000,
!            MODIFICATION DES APPELS BLAS (ROUTINE ASTER BL...),
!            IMPLICIT NONE.
! INTRINSIC FUNCTION
!   MAX, MIN, DBLE, ABS
! ENDLIB
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
! aslint: disable=W1501
    implicit none
!
!     .. SCALAR ARGUMENTS ..
#include "asterf_types.h"
#include "asterc/isbaem.h"
#include "asterc/matfpe.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/as_dlaln2.h"
#include "asterfort/as_dlasy2.h"
#include "asterfort/xerbla.h"
#include "blas/ddot.h"
#include "blas/dlange.h"
#include "blas/dscal.h"
#include "blas/lsame.h"
    character(len=1) :: trana, tranb
    integer(kind=8) :: info, isgn, lda, ldb, ldc, m, n
    real(kind=8) :: scale
!     ..
!     .. ARRAY ARGUMENTS ..
    real(kind=8) :: a(lda, *), b(ldb, *), c(ldc, *)
!     ..
!     .. PARAMETERS ..
    real(kind=8) :: zero, one
    parameter(zero=0.0d+0, one=1.0d+0)
!     ..
!     .. LOCAL SCALARS ..
    aster_logical :: notrna, notrnb
    integer(kind=8) :: ierr, j, k, k1, k2, knext, l, l1, l2, lnext
    real(kind=8) :: a11, bignum, da11, db, eps, scaloc, sgn, smin, smlnum, suml
    real(kind=8) :: sumr, xnorm
!     ..
!     .. LOCAL ARRAYS ..
    real(kind=8) :: dum(1), vec(2, 2), x(2, 2)
    blas_int :: b_incx, b_n
    blas_int :: b_incy
    blas_int :: b_lda, b_ldabis, b_m, b_mbis, b_nbis
!     ..
!     .. EXTERNAL FUNCTIONS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
    call matfpe(-1)
!
! DUE TO CRS513
    dum(1) = 0.d0
!
!     DECODE AND TEST INPUT PARAMETERS
!
    notrna = lsame(trana, 'N')
    notrnb = lsame(tranb, 'N')
!
    info = 0
    if (.not. notrna .and. .not. lsame(trana, 'T') .and. .not. lsame(trana, 'C')) then
        info = -1
    else if (.not. notrnb .and. .not. lsame(tranb, 'T') .and. .not. &
             lsame(tranb, 'C')) then
        info = -2
    else if (isgn .ne. 1 .and. isgn .ne. -1) then
        info = -3
    else if (m .lt. 0) then
        info = -4
    else if (n .lt. 0) then
        info = -5
    else if (lda .lt. max(1, m)) then
        info = -7
    else if (ldb .lt. max(1, n)) then
        info = -9
    else if (ldc .lt. max(1, m)) then
        info = -11
    end if
    if (info .ne. 0) then
        call xerbla('DLRSYL', -info)
        goto 1000
    end if
!
!     QUICK RETURN IF POSSIBLE
!
    if (m .eq. 0 .or. n .eq. 0) goto 1000
!
!     SET CONSTANTS TO CONTROL OVERFLOW
!
    eps = r8prem()*0.5d0*isbaem()
    smlnum = r8miem()
    bignum = one/smlnum
    smlnum = smlnum*dble(m*n)/eps
    bignum = one/smlnum
!
    b_lda = to_blas_int(lda)
    b_m = to_blas_int(m)
    b_n = to_blas_int(m)
    smin = max( &
           smlnum, eps*dlange('M', b_m, b_n, a, b_lda, dum), &
           eps*dlange('M', b_mbis, b_nbis, b, b_ldabis, dum) &
           )
!
    scale = one
    sgn = isgn
!
    if (notrna .and. notrnb) then
!
!        SOLVE    A*X + ISGN*X*B = SCALE*C.
!
!        THE (K,L)TH BLOCK OF X IS DETERMINED STARTING FROM
!        BOTTOM-LEFT CORNER COLUMN BY COLUMN BY
!
!         A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        WHERE
!                  M                         L-1
!        R(K,L) = SUM (A(K,I)*X(I,L)) + ISGN*SUM (X(K,J)*B(J,L)).
!                I=K+1                       J=1
!
!        START COLUMN LOOP (INDEX = L)
!        L1 (L2) : COLUMN INDEX OF THE FIRST (FIRST) ROW OF X(K,L).
!
        lnext = 1
        do l = 1, n
            if (l .lt. lnext) goto 60
            if (l .eq. n) then
                l1 = l
                l2 = l
            else
                if (b(l+1, l) .ne. zero) then
                    l1 = l
                    l2 = l+1
                    lnext = l+2
                else
                    l1 = l
                    l2 = l
                    lnext = l+1
                end if
            end if
!
!           START ROW LOOP (INDEX = K)
!           K1 (K2): ROW INDEX OF THE FIRST (LAST) ROW OF X(K,L).
!
            knext = m
            do k = m, 1, -1
                if (k .gt. knext) goto 50
                if (k .eq. 1) then
                    k1 = k
                    k2 = k
                else
                    knext = k-1
                    if (a(k, knext) .ne. zero) then
                        k1 = knext
                        k2 = k
                        knext = k-2
                    else
                        k1 = k
                        k2 = k
                    end if
                end if
!
                if (l1 .eq. l2 .and. k1 .eq. k2) then
                    b_n = to_blas_int(m-k1)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k1+1, m)), b_incx, c(min(k1+1, m), l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l1), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
                    scaloc = one
!
                    a11 = a(k1, k1)+sgn*b(l1, l1)
                    da11 = abs(a11)
                    if (da11 .le. smin) then
                        a11 = smin
                        da11 = smin
                        info = 1
                    end if
                    db = abs(vec(1, 1))
                    if (da11 .lt. one .and. db .gt. one) then
                        if (db .gt. bignum*da11) scaloc = one/db
                    end if
                    x(1, 1) = (vec(1, 1)*scaloc)/a11
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
!
                else if (l1 .eq. l2 .and. k1 .ne. k2) then
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k2+1, m)), b_incx, c(min(k2+1, m), l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l1), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k2, min(k2+1, m)), b_incx, c(min(k2+1, m), l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k2, 1), b_incx, b(1, l1), b_incy)
                    vec(2, 1) = c(k2, l1)-(suml+sgn*sumr)
!
                    call as_dlaln2(.false._1, 2, 1, smin, one, &
                                   a(k1, k1), lda, one, one, vec, &
                                   2, -sgn*b(l1, l1), zero, x, 2, &
                                   scaloc, xnorm, ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k2, l1) = x(2, 1)
!
                else if (l1 .ne. l2 .and. k1 .eq. k2) then
!
                    b_n = to_blas_int(m-k1)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k1+1, m)), b_incx, c(min(k1+1, m), l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l1), b_incy)
                    vec(1, 1) = sgn*(c(k1, l1)-(suml+sgn*sumr))
!
                    b_n = to_blas_int(m-k1)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k1+1, m)), b_incx, c(min(k1+1, m), l2), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l2), b_incy)
                    vec(2, 1) = sgn*(c(k1, l2)-(suml+sgn*sumr))
!
                    call as_dlaln2(.true._1, 2, 1, smin, one, &
                                   b(l1, l1), ldb, one, one, vec, &
                                   2, -sgn*a(k1, k1), zero, x, 2, &
                                   scaloc, xnorm, ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k1, l2) = x(2, 1)
!
                else if (l1 .ne. l2 .and. k1 .ne. k2) then
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k2+1, m)), b_incx, c(min(k2+1, m), l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l1), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k2+1, m)), b_incx, c(min(k2+1, m), l2), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l2), b_incy)
                    vec(1, 2) = c(k1, l2)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k2, min(k2+1, m)), b_incx, c(min(k2+1, m), l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k2, 1), b_incx, b(1, l1), b_incy)
                    vec(2, 1) = c(k2, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k2, min(k2+1, m)), b_incx, c(min(k2+1, m), l2), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k2, 1), b_incx, b(1, l2), b_incy)
                    vec(2, 2) = c(k2, l2)-(suml+sgn*sumr)
!
                    call as_dlasy2(.false._1, .false._1, isgn, 2, 2, &
                                   a(k1, k1), lda, b(l1, l1), ldb, vec, &
                                   2, scaloc, x, 2, xnorm, &
                                   ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k1, l2) = x(1, 2)
                    c(k2, l1) = x(2, 1)
                    c(k2, l2) = x(2, 2)
                end if
!
50              continue
            end do
!
60          continue
        end do
!
    else if (.not. notrna .and. notrnb) then
!
!        SOLVE    A' *X + ISGN*X*B = SCALE*C.
!
!        THE (K,L)TH BLOCK OF X IS DETERMINED STARTING FROM
!        UPPER-LEFT CORNER COLUMN BY COLUMN BY
!
!          A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        WHERE
!                   K-1                        L-1
!          R(K,L) = SUM (A(I,K)'*X(I,L)) +ISGN*SUM (X(K,J)*B(J,L))
!                   I=1                        J=1
!
!        START COLUMN LOOP (INDEX = L)
!        L1 (L2): COLUMN INDEX OF THE FIRST (LAST) ROW OF X(K,L)
!
        lnext = 1
        do l = 1, n
            if (l .lt. lnext) goto 120
            if (l .eq. n) then
                l1 = l
                l2 = l
            else
                if (b(l+1, l) .ne. zero) then
                    l1 = l
                    l2 = l+1
                    lnext = l+2
                else
                    l1 = l
                    l2 = l
                    lnext = l+1
                end if
            end if
!
!           START ROW LOOP (INDEX = K)
!           K1 (K2): ROW INDEX OF THE FIRST (LAST) ROW OF X(K,L)
!
            knext = 1
            do k = 1, m
                if (k .lt. knext) goto 110
                if (k .eq. m) then
                    k1 = k
                    k2 = k
                else
                    if (a(k+1, k) .ne. zero) then
                        k1 = k
                        k2 = k+1
                        knext = k+2
                    else
                        k1 = k
                        k2 = k
                        knext = k+1
                    end if
                end if
!
                if (l1 .eq. l2 .and. k1 .eq. k2) then
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l1), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
                    scaloc = one
!
                    a11 = a(k1, k1)+sgn*b(l1, l1)
                    da11 = abs(a11)
                    if (da11 .le. smin) then
                        a11 = smin
                        da11 = smin
                        info = 1
                    end if
                    db = abs(vec(1, 1))
                    if (da11 .lt. one .and. db .gt. one) then
                        if (db .gt. bignum*da11) scaloc = one/db
                    end if
                    x(1, 1) = (vec(1, 1)*scaloc)/a11
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
!
                else if (l1 .eq. l2 .and. k1 .ne. k2) then
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l1), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k2), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k2, 1), b_incx, b(1, l1), b_incy)
                    vec(2, 1) = c(k2, l1)-(suml+sgn*sumr)
!
                    call as_dlaln2(.true._1, 2, 1, smin, one, &
                                   a(k1, k1), lda, one, one, vec, &
                                   2, -sgn*b(l1, l1), zero, x, 2, &
                                   scaloc, xnorm, ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k2, l1) = x(2, 1)
!
                else if (l1 .ne. l2 .and. k1 .eq. k2) then
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l1), b_incy)
                    vec(1, 1) = sgn*(c(k1, l1)-(suml+sgn*sumr))
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l2), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l2), b_incy)
                    vec(2, 1) = sgn*(c(k1, l2)-(suml+sgn*sumr))
!
                    call as_dlaln2(.true._1, 2, 1, smin, one, &
                                   b(l1, l1), ldb, one, one, vec, &
                                   2, -sgn*a(k1, k1), zero, x, 2, &
                                   scaloc, xnorm, ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k1, l2) = x(2, 1)
!
                else if (l1 .ne. l2 .and. k1 .ne. k2) then
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l1), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l2), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k1, 1), b_incx, b(1, l2), b_incy)
                    vec(1, 2) = c(k1, l2)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k2), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k2, 1), b_incx, b(1, l1), b_incy)
                    vec(2, 1) = c(k2, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k2), b_incx, c(1, l2), b_incy)
                    b_n = to_blas_int(l1-1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(1)
                    sumr = ddot(b_n, c(k2, 1), b_incx, b(1, l2), b_incy)
                    vec(2, 2) = c(k2, l2)-(suml+sgn*sumr)
!
                    call as_dlasy2(.true._1, .false._1, isgn, 2, 2, &
                                   a(k1, k1), lda, b(l1, l1), ldb, vec, &
                                   2, scaloc, x, 2, xnorm, &
                                   ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k1, l2) = x(1, 2)
                    c(k2, l1) = x(2, 1)
                    c(k2, l2) = x(2, 2)
                end if
!
110             continue
            end do
120         continue
        end do
!
    else if (.not. notrna .and. .not. notrnb) then
!
!        SOLVE    A'*X + ISGN*X*B' = SCALE*C.
!
!        THE (K,L)TH BLOCK OF X IS DETERMINED STARTING FROM
!        TOP-RIGHT CORNER COLUMN BY COLUMN BY
!
!           A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
!
!        WHERE
!                     K-1                          N
!            R(K,L) = SUM (A(I,K)'*X(I,L)) + ISGN*SUM (X(K,J)*B(L,J)').
!                     I=1                        J=L+1
!
!        START COLUMN LOOP (INDEX = L)
!        L1 (L2): COLUMN INDEX OF THE FIRST (LAST) ROW OF X(K,L)
!
        lnext = n
        do l = n, 1, -1
            if (l .gt. lnext) goto 180
            if (l .eq. 1) then
                l1 = l
                l2 = l
            else
                lnext = l-1
                if (b(l, lnext) .ne. zero) then
                    l1 = lnext
                    l2 = l
                    lnext = l-2
                else
                    l1 = l
                    l2 = l
                end if
            end if
!
!           START ROW LOOP (INDEX = K)
!           K1 (K2): ROW INDEX OF THE FIRST (LAST) ROW OF X(K,L)
!
            knext = 1
            do k = 1, m
                if (k .lt. knext) goto 170
                if (k .eq. m) then
                    k1 = k
                    k2 = k
                else
                    if (a(k+1, k) .ne. zero) then
                        k1 = k
                        k2 = k+1
                        knext = k+2
                    else
                        k1 = k
                        k2 = k
                        knext = k+1
                    end if
                end if
!
                if (l1 .eq. l2 .and. k1 .eq. k2) then
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(n-l1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l1+1, n)), b_incx, b(l1, min(l1+1, n)), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
                    scaloc = one
!
                    a11 = a(k1, k1)+sgn*b(l1, l1)
                    da11 = abs(a11)
                    if (da11 .le. smin) then
                        a11 = smin
                        da11 = smin
                        info = 1
                    end if
                    db = abs(vec(1, 1))
                    if (da11 .lt. one .and. db .gt. one) then
                        if (db .gt. bignum*da11) scaloc = one/db
                    end if
                    x(1, 1) = (vec(1, 1)*scaloc)/a11
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
!
                else if (l1 .eq. l2 .and. k1 .ne. k2) then
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l2+1, n)), b_incx, b(l1, min(l2+1, n)), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k2), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k2, min(l2+1, n)), b_incx, b(l1, min(l2+1, n)), b_incy)
                    vec(2, 1) = c(k2, l1)-(suml+sgn*sumr)
!
                    call as_dlaln2(.true._1, 2, 1, smin, one, &
                                   a(k1, k1), lda, one, one, vec, &
                                   2, -sgn*b(l1, l1), zero, x, 2, &
                                   scaloc, xnorm, ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k2, l1) = x(2, 1)
!
                else if (l1 .ne. l2 .and. k1 .eq. k2) then
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l2+1, n)), b_incx, b(l1, min(l2+1, n)), b_incy)
                    vec(1, 1) = sgn*(c(k1, l1)-(suml+sgn*sumr))
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l2), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l2+1, n)), b_incx, b(l2, min(l2+1, n)), b_incy)
                    vec(2, 1) = sgn*(c(k1, l2)-(suml+sgn*sumr))
!
                    call as_dlaln2(.false._1, 2, 1, smin, one, &
                                   b(l1, l1), ldb, one, one, vec, &
                                   2, -sgn*a(k1, k1), zero, x, 2, &
                                   scaloc, xnorm, ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k1, l2) = x(2, 1)
!
                else if (l1 .ne. l2 .and. k1 .ne. k2) then
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l2+1, n)), b_incx, b(l1, min(l2+1, n)), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k1), b_incx, c(1, l2), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l2+1, n)), b_incx, b(l2, min(l2+1, n)), b_incy)
                    vec(1, 2) = c(k1, l2)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k2), b_incx, c(1, l1), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k2, min(l2+1, n)), b_incx, b(l1, min(l2+1, n)), b_incy)
                    vec(2, 1) = c(k2, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(k1-1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(1, k2), b_incx, c(1, l2), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k2, min(l2+1, n)), b_incx, b(l2, min(l2+1, n)), b_incy)
                    vec(2, 2) = c(k2, l2)-(suml+sgn*sumr)
!
                    call as_dlasy2(.true._1, .true._1, isgn, 2, 2, &
                                   a(k1, k1), lda, b(l1, l1), ldb, vec, &
                                   2, scaloc, x, 2, xnorm, &
                                   ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k1, l2) = x(1, 2)
                    c(k2, l1) = x(2, 1)
                    c(k2, l2) = x(2, 2)
                end if
!
170             continue
            end do
180         continue
        end do
!
    else if (notrna .and. .not. notrnb) then
!
!        SOLVE    A*X + ISGN*X*B' = SCALE*C.
!
!        THE (K,L)TH BLOCK OF X IS DETERMINED STARTING FROM
!        BOTTOM-RIGHT CORNER COLUMN BY COLUMN BY
!
!            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
!
!        WHERE
!                      M                          N
!            R(K,L) = SUM (A(K,I)*X(I,L)) + ISGN*SUM (X(K,J)*B(L,J)').
!                    I=K+1                      J=L+1
!
!        START COLUMN LOOP (INDEX = L)
!        L1 (L2): COLUMN INDEX OF THE FIRST (LAST) ROW OF X(K,L)
!
        lnext = n
        do l = n, 1, -1
            if (l .gt. lnext) goto 240
            if (l .eq. 1) then
                l1 = l
                l2 = l
            else
                lnext = l-1
                if (b(l, lnext) .ne. zero) then
                    l1 = lnext
                    l2 = l
                    lnext = l-2
                else
                    l1 = l
                    l2 = l
                end if
            end if
!
!           START ROW LOOP (INDEX = K)
!           K1 (K2): ROW INDEX OF THE FIRST (LAST) ROW OF X(K,L)
!
            knext = m
            do k = m, 1, -1
                if (k .gt. knext) goto 230
                if (k .eq. 1) then
                    k1 = k
                    k2 = k
                else
                    knext = k-1
                    if (a(k, knext) .ne. zero) then
                        k1 = knext
                        k2 = k
                        knext = k-2
                    else
                        k1 = k
                        k2 = k
                    end if
                end if
!
                if (l1 .eq. l2 .and. k1 .eq. k2) then
                    b_n = to_blas_int(m-k1)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k1+1, m)), b_incx, c(min(k1+1, m), l1), b_incy)
                    b_n = to_blas_int(n-l1)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l1+1, n)), b_incx, b(l1, min(l1+1, n)), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
                    scaloc = one
!
                    a11 = a(k1, k1)+sgn*b(l1, l1)
                    da11 = abs(a11)
                    if (da11 .le. smin) then
                        a11 = smin
                        da11 = smin
                        info = 1
                    end if
                    db = abs(vec(1, 1))
                    if (da11 .lt. one .and. db .gt. one) then
                        if (db .gt. bignum*da11) scaloc = one/db
                    end if
                    x(1, 1) = (vec(1, 1)*scaloc)/a11
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
!
                else if (l1 .eq. l2 .and. k1 .ne. k2) then
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k2+1, m)), b_incx, c(min(k2+1, m), l1), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l2+1, n)), b_incx, b(l1, min(l2+1, n)), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k2, min(k2+1, m)), b_incx, c(min(k2+1, m), l1), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k2, min(l2+1, n)), b_incx, b(l1, min(l2+1, n)), b_incy)
                    vec(2, 1) = c(k2, l1)-(suml+sgn*sumr)
!
                    call as_dlaln2(.false._1, 2, 1, smin, one, &
                                   a(k1, k1), lda, one, one, vec, &
                                   2, -sgn*b(l1, l1), zero, x, 2, &
                                   scaloc, xnorm, ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k2, l1) = x(2, 1)
!
                else if (l1 .ne. l2 .and. k1 .eq. k2) then
!
                    b_n = to_blas_int(m-k1)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k1+1, m)), b_incx, c(min(k1+1, m), l1), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l2+1, n)), b_incx, b(l1, min(l2+1, n)), b_incy)
                    vec(1, 1) = sgn*(c(k1, l1)-(suml+sgn*sumr))
!
                    b_n = to_blas_int(m-k1)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k1+1, m)), b_incx, c(min(k1+1, m), l2), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l2+1, n)), b_incx, b(l2, min(l2+1, n)), b_incy)
                    vec(2, 1) = sgn*(c(k1, l2)-(suml+sgn*sumr))
!
                    call as_dlaln2(.false._1, 2, 1, smin, one, &
                                   b(l1, l1), ldb, one, one, vec, &
                                   2, -sgn*a(k1, k1), zero, x, 2, &
                                   scaloc, xnorm, ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k1, l2) = x(2, 1)
!
                else if (l1 .ne. l2 .and. k1 .ne. k2) then
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k2+1, m)), b_incx, c(min(k2+1, m), l1), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l2+1, n)), b_incx, b(l1, min(l2+1, n)), b_incy)
                    vec(1, 1) = c(k1, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k1, min(k2+1, m)), b_incx, c(min(k2+1, m), l2), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k1, min(l2+1, n)), b_incx, b(l2, min(l2+1, n)), b_incy)
                    vec(1, 2) = c(k1, l2)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k2, min(k2+1, m)), b_incx, c(min(k2+1, m), l1), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k2, min(l2+1, n)), b_incx, b(l1, min(l2+1, n)), b_incy)
                    vec(2, 1) = c(k2, l1)-(suml+sgn*sumr)
!
                    b_n = to_blas_int(m-k2)
                    b_incx = to_blas_int(lda)
                    b_incy = to_blas_int(1)
                    suml = ddot(b_n, a(k2, min(k2+1, m)), b_incx, c(min(k2+1, m), l2), b_incy)
                    b_n = to_blas_int(n-l2)
                    b_incx = to_blas_int(ldc)
                    b_incy = to_blas_int(ldb)
                    sumr = ddot(b_n, c(k2, min(l2+1, n)), b_incx, b(l2, min(l2+1, n)), b_incy)
                    vec(2, 2) = c(k2, l2)-(suml+sgn*sumr)
!
                    call as_dlasy2(.false._1, .true._1, isgn, 2, 2, &
                                   a(k1, k1), lda, b(l1, l1), ldb, vec, &
                                   2, scaloc, x, 2, xnorm, &
                                   ierr)
                    if (ierr .ne. 0) info = 1
!
                    if (scaloc .ne. one) then
                        do j = 1, n
                            b_n = to_blas_int(m)
                            b_incx = to_blas_int(1)
                            call dscal(b_n, scaloc, c(1, j), b_incx)
                        end do
                        scale = scale*scaloc
                    end if
                    c(k1, l1) = x(1, 1)
                    c(k1, l2) = x(1, 2)
                    c(k2, l1) = x(2, 1)
                    c(k2, l2) = x(2, 2)
                end if
!
230             continue
            end do
240         continue
        end do
!
    end if
!
1000 continue
    call matfpe(1)
!
!     END OF DLRSYL
!
end subroutine
