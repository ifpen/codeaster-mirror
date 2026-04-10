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
subroutine as_zneigh(rnorm, n, h, ldh, ritz, &
                     bounds, q, ldq, workl, rwork, &
                     ierr)
!
!     SUBROUTINE ARPACK CALCULANT LES MODES PROPRES DE LA MATRICE DE
!     HESSENBERG AINSI QUE LEURS ERREURS.
!-----------------------------------------------------------------------
!\BEGINDOC
!
!\NAME: ZNEIGH
!
!\DESCRIPTION:
!  COMPUTE THE EIGENVALUES OF THE CURRENT UPPER HESSENBERG MATRIX
!  AND THE CORRESPONDING RITZ ESTIMATES GIVEN THE CURRENT RESIDUAL NORM.
!
!\USAGE:
!  CALL ZNEIGH
!     ( RNORM, N, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, RWORK, IERR )
!
!\ARGUMENTS
!  RNORM   DOUBLE PRECISION SCALAR.  (INPUT)
!          RESIDUAL NORM CORRESPONDING TO THE CURRENT UPPER HESSENBERG
!          MATRIX H.
!
!  N       INTEGER.  (INPUT)
!          SIZE OF THE MATRIX H.
!
!  H       COMPLEX*16 N BY N ARRAY.  (INPUT)
!          H CONTAINS THE CURRENT UPPER HESSENBERG MATRIX.
!
!  LDH     INTEGER.  (INPUT)
!          LEADING DIMENSION OF H EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  RITZ    COMPLEX*16 ARRAY OF LENGTH N.  (OUTPUT)
!          ON OUTPUT, RITZ(1:N) CONTAINS THE EIGENVALUES OF H.
!
!  BOUNDS  COMPLEX*16 ARRAY OF LENGTH N.  (OUTPUT)
!          ON OUTPUT, BOUNDS CONTAINS THE RITZ ESTIMATES ASSOCIATED WITH
!          THE EIGENVALUES HELD IN RITZ.  THIS IS EQUAL TO RNORM
!          TIMES THE LAST COMPONENTS OF THE EIGENVECTORS CORRESPONDING
!          TO THE EIGENVALUES IN RITZ.
!
!  Q       COMPLEX*16 N BY N ARRAY.  (WORKSPACE)
!          WORKSPACE NEEDED TO STORE THE EIGENVECTORS OF H.
!
!  LDQ     INTEGER.  (INPUT)
!          LEADING DIMENSION OF Q EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  WORKL   COMPLEX*16 WORK ARRAY OF LENGTH N**2 + 3*N.  (WORKSPACE)
!          PRIVATE (REPLICATED) ARRAY ON EACH PE OR ARRAY ALLOCATED ON
!          THE FRONT END.  THIS IS NEEDED TO KEEP THE FULL SCHUR FORM
!          OF H AND ALSO IN THE CALCULATION OF THE EIGENVECTORS OF H.
!
!  RWORK   DOUBLE PRECISION  WORK ARRAY OF LENGTH N (WORKSPACE)
!          PRIVATE (REPLICATED) ARRAY ON EACH PE OR ARRAY ALLOCATED ON
!          THE FRONT END.
!
!  IERR    INTEGER.  (OUTPUT)
!          ERROR EXIT FLAG FROM ZLAHQR OR ZTREVC.
!
!\ENDDOC
!
!-----------------------------------------------------------------------
!
!\BEGINLIB
!
!\LOCAL VARIABLES:
!     XXXXXX  COMPLEX*16
!
!\ROUTINES CALLED:
!     IVOUT   ARPACK UTILITY ROUTINE THAT PRINTS INTEGERS.
!     ZMOUT   ARPACK UTILITY ROUTINE THAT PRINTS MATRICES
!     ZVOUT   ARPACK UTILITY ROUTINE THAT PRINTS VECTORS.
!     DVOUT   ARPACK UTILITY ROUTINE THAT PRINTS VECTORS.
!     ZLACPY  LAPACK MATRIX COPY ROUTINE.
!     ZLAHQR  LAPACK ROUTINE TO COMPUTE THE SCHUR FORM OF AN
!             UPPER HESSENBERG MATRIX.
!     ZLASET  LAPACK MATRIX INITIALIZATION ROUTINE.
!     ZTREVC  LAPACK ROUTINE TO COMPUTE THE EIGENVECTORS OF A MATRIX
!             IN UPPER TRIANGULAR FORM
!     ZCOPY   LEVEL 1 BLAS THAT COPIES ONE VECTOR TO ANOTHER.
!     ZDSCAL  LEVEL 1 BLAS THAT SCALES A COMPLEX VECTOR BY A REAL
!             NUMBER.
!     DZNRM2  LEVEL 1 BLAS THAT COMPUTES THE NORM OF A VECTOR.
!
!
!\AUTHOR
!     DANNY SORENSEN               PHUONG VU
!     RICHARD LEHOUCQ              CRPC / RICE UNIVERSITY
!     DEPT. OF COMPUTATIONAL &     HOUSTON, TEXAS
!     APPLIED MATHEMATICS
!     RICE UNIVERSITY
!     HOUSTON, TEXAS
!
!\SCCS INFORMATION: @(#)
! FILE: NEIGH.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
!
!\REMARKS
!     NONE
!
!\ENDLIB
!
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     %-----------------------------%
!     | INCLUDE FILES FOR DEBUGGING |
!     %-----------------------------%
!
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterfort/as_ztrevc.h"
#include "asterfort/as_zmout.h"
#include "asterfort/as_zvout.h"
#include "blas/dznrm2.h"
#include "blas/zcopy.h"
#include "blas/zdscal.h"
#include "blas/zlacpy.h"
#include "blas/zlahqr.h"
#include "blas/zlaset.h"
    integer(kind=8) :: logfil, ndigit, mgetv0, mnaupd, mnaup2, mnaitr, mneigh, mnapps
    integer(kind=8) :: mngets, mneupd
    common/debug/&
     &  logfil, ndigit, mgetv0,&
     &  mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd
!
!     %------------------%
!     | SCALAR ARGUMENTS |
!     %------------------%
!
    integer(kind=8) :: ierr, n, ldh, ldq
    integer(kind=4) :: ierr4
    real(kind=8) :: rnorm
!
!     %-----------------%
!     | ARRAY ARGUMENTS |
!     %-----------------%
!
    complex(kind=8) :: bounds(n), h(ldh, n), q(ldq, n), ritz(n), workl(n*(n+3))
    real(kind=8) :: rwork(n)
!
!     %------------%
!     | PARAMETERS |
!     %------------%
!
    complex(kind=8) :: one, zero
    real(kind=8) :: rone
    parameter(one=(1.0d+0, 0.0d+0), zero=(0.0d+0, 0.0d+0),&
     &           rone=1.0d+0)
!
!     %------------------------%
!     | LOCAL SCALARS & ARRAYS |
!     %------------------------%
!
    aster_logical :: select(1)
    integer(kind=8) :: j, msglvl
    complex(kind=8) :: vl(1)
    real(kind=8) :: temp
    blas_int :: b_incx, b_incy, b_n
    blas_int :: b_lda, b_ldb, b_m
    blas_int :: b_ihi, b_ihiz, b_ilo, b_iloz, b_ldh, b_ldz
!
!     %--------------------%
!     | EXTERNAL FUNCTIONS |
!     %--------------------%
!
!
!     %-----------------------%
!     | EXECUTABLE STATEMENTS |
!     %-----------------------%
!
    call matfpe(-1)
!
!     %-------------------------------%
!     | INITIALIZE TIMING STATISTICS  |
!     | & MESSAGE LEVEL FOR DEBUGGING |
!     %-------------------------------%
!
    msglvl = mneigh
!
    if (msglvl .gt. 2) then
        call as_zmout(logfil, n, n, h, ldh, &
                      ndigit, '_NEIGH: ENTERING UPPER HESSENBERG MATRIX H ')
    end if
!
!     %----------------------------------------------------------%
!     | 1. COMPUTE THE EIGENVALUES, THE LAST COMPONENTS OF THE   |
!     |    CORRESPONDING SCHUR VECTORS AND THE FULL SCHUR FORM T |
!     |    OF THE CURRENT UPPER HESSENBERG MATRIX H.             |
!     |    ZLAHQR RETURNS THE FULL SCHUR FORM OF H               |
!     |    IN WORKL(1:N**2), AND THE SCHUR VECTORS IN Q.         |
!     %----------------------------------------------------------%
!
    b_ldb = to_blas_int(n)
    b_lda = to_blas_int(ldh)
    b_m = to_blas_int(n)
    b_n = to_blas_int(n)
    call zlacpy('A', b_m, b_n, h, b_lda, &
                workl, b_ldb)
    b_lda = to_blas_int(ldq)
    b_m = to_blas_int(n)
    b_n = to_blas_int(n)
    call zlaset('A', b_m, b_n, zero, one, &
                q, b_lda)
!
    b_ldz = to_blas_int(ldq)
    b_ldh = to_blas_int(ldh)
    b_n = to_blas_int(n)
    b_ilo = to_blas_int(1)
    b_ihi = to_blas_int(n)
    b_iloz = to_blas_int(1)
    b_ihiz = to_blas_int(n)
    call zlahqr(.true._1, .true._1, b_n, b_ilo, b_ihi, &
                workl, b_ldh, ritz, b_iloz, b_ihiz, &
                q, b_ldz, ierr4)
    if (ierr4 .ne. 0) goto 9000
!
    b_n = to_blas_int(n)
    b_incx = to_blas_int(ldq)
    b_incy = to_blas_int(1)
    call zcopy(b_n, q(n-1, 1), b_incx, bounds, b_incy)
    if (msglvl .gt. 1) then
        call as_zvout(logfil, n, bounds, ndigit, '_NEIGH: LAST ROW OF THE SCHUR MATRIX FOR H')
    end if
!
!     %----------------------------------------------------------%
!     | 2. COMPUTE THE EIGENVECTORS OF THE FULL SCHUR FORM T AND |
!     |    APPLY THE SCHUR VECTORS TO GET THE CORRESPONDING      |
!     |    EIGENVECTORS.                                         |
!     %----------------------------------------------------------%
!
    call as_ztrevc('R', 'B', select, n, workl, &
                   n, vl, n, q, ldq, &
                   n, n, workl(n*n+1), rwork, ierr)
!
    if (ierr .ne. 0) goto 9000
!
!     %------------------------------------------------%
!     | SCALE THE RETURNING EIGENVECTORS SO THAT THEIR |
!     | EUCLIDEAN NORMS ARE ALL ONE. LAPACK SUBROUTINE |
!     | ZTREVC RETURNS EACH EIGENVECTOR NORMALIZED SO  |
!     | THAT THE ELEMENT OF LARGEST MAGNITUDE HAS      |
!     | MAGNITUDE 1; HERE THE MAGNITUDE OF A COMPLEX   |
!     | NUMBER (X,Y) IS TAKEN TO BE |X| + |Y|.         |
!     %------------------------------------------------%
!
    do j = 1, n
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        temp = dznrm2(b_n, q(1, j), b_incx)
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        call zdscal(b_n, rone/temp, q(1, j), b_incx)
    end do
!
    if (msglvl .gt. 1) then
        b_n = to_blas_int(n)
        b_incx = to_blas_int(ldq)
        b_incy = to_blas_int(1)
        call zcopy(b_n, q(n, 1), b_incx, workl, b_incy)
        call as_zvout(logfil, n, workl, ndigit, '_NEIGH: LAST ROW OF THE EIGENVECTOR MATRIX FOR H')
    end if
!
!     %----------------------------%
!     | COMPUTE THE RITZ ESTIMATES |
!     %----------------------------%
!
    b_n = to_blas_int(n)
    b_incx = to_blas_int(n)
    b_incy = to_blas_int(1)
    call zcopy(b_n, q(n, 1), b_incx, bounds, b_incy)
    b_n = to_blas_int(n)
    b_incx = to_blas_int(1)
    call zdscal(b_n, rnorm, bounds, b_incx)
!
    if (msglvl .gt. 2) then
        call as_zvout(logfil, n, ritz, ndigit, '_NEIGH: THE EIGENVALUES OF H')
        call as_zvout(logfil, n, bounds, ndigit, '_NEIGH: RITZ ESTIMATES FOR THE EIGENVALUES OF H')
    end if
!
!
9000 continue
    call matfpe(1)
!
!     %---------------%
!     | END OF ZNEIGH |
!     %---------------%
!
end subroutine
