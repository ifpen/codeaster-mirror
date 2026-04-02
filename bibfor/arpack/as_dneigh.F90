! --------------------------------------------------------------------
! Copyright (C) LAPACK
! Copyright (C) 2007 - 2025 - EDF R&D - www.code-aster.org
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
subroutine as_dneigh(rnorm, n, h, ldh, ritzr, &
                  ritzi, bounds, q, ldq, workl, &
                  ierr)
!
!     SUBROUTINE ARPACK CALCULANT LES MODES PROPRES DE LA MATRICE DE
!     HESSENBERG AINSI QUE LEURS ERREURS.
!-----------------------------------------------------------------------
! BEGINDOC
!
! DESCRIPTION:
!  COMPUTE THE EIGENVALUES OF THE CURRENT UPPER HESSENBERG MATRIX
!  AND THE CORRESPONDING RITZ ESTIMATES GIVEN THE CURRENT RESIDUAL NORM.
!
! ARGUMENTS
!  RNORM   REAL*8 SCALAR.  (INPUT)
!          RESIDUAL NORM CORRESPONDING TO THE CURRENT UPPER HESSENBERG
!          MATRIX H.
!
!  N       INTEGER.  (INPUT)
!          SIZE OF THE MATRIX H.
!
!  H       REAL*8 N BY N ARRAY.  (INPUT)
!          H CONTAINS THE CURRENT UPPER HESSENBERG MATRIX.
!
!  LDH     INTEGER.  (INPUT)
!          LEADING DIMENSION OF H EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  RITZR,  REAL*8 ARRAYS OF LENGTH N.  (OUTPUT)
!  RITZI   ON OUTPUT, RITZR(1:N) (RESP. RITZI(1:N)) CONTAINS THE REAL
!          (RESPECTIVELY IMAGINARY) PARTS OF THE EIGENVALUES OF H.
!
!  BOUNDS  REAL*8 ARRAY OF LENGTH N.  (OUTPUT)
!         ON OUTPUT, BOUNDS CONTAINS THE RITZ ESTIMATES ASSOCIATED WITH
!         THE EIGENVALUES RITZR AND RITZI.  THIS IS EQUAL TO RNORM
!         TIMES THE LAST COMPONENTS OF THE EIGENVECTORS CORRESPONDING
!         TO THE EIGENVALUES IN RITZR AND RITZI.
!
!  Q       REAL*8 N BY N ARRAY.  (WORKSPACE)
!          WORKSPACE NEEDED TO STORE THE EIGENVECTORS OF H.
!
!  LDQ     INTEGER.  (INPUT)
!          LEADING DIMENSION OF Q EXACTLY AS DECLARED IN THE CALLING
!          PROGRAM.
!
!  WORKL   REAL*8 WORK ARRAY OF LENGTH N**2 + 3*N.  (WORKSPACE)
!          PRIVATE (REPLICATED) ARRAY ON EACH PE OR ARRAY ALLOCATED ON
!          THE FRONT END.  THIS IS NEEDED TO KEEP THE FULL SCHUR FORM
!          OF H AND ALSO IN THE CALCULATION OF THE EIGENVECTORS OF H.
!
!  IERR    INTEGER.  (OUTPUT)
!          ERROR EXIT FLAG FROM DLAQRB OR DTREVC.
!
! ENDDOC
!-----------------------------------------------------------------------
!
! BEGINLIB
!
! ROUTINES CALLED:
!     DLAQRB  ARPACK ROUTINE TO COMPUTE THE REAL SCHUR FORM OF AN
!             UPPER HESSENBERG MATRIX AND LAST ROW OF THE SCHUR VECTORS.
!     DMOUT   ARPACK UTILITY ROUTINE THAT PRINTS MATRICES
!     DVOUT   ARPACK UTILITY ROUTINE THAT PRINTS VECTORS.
!     DLACPY  LAPACK MATRIX COPY ROUTINE.
!     DLAPY2  LAPACK ROUTINE TO COMPUTE SQRT(X**2+Y**2) CAREFULLY.
!     DTREVC  LAPACK ROUTINE TO COMPUTE THE EIGENVECTORS OF A MATRIX
!             IN UPPER QUASI-TRIANGULAR FORM
!     DGEMV   LEVEL 2 BLAS ROUTINE FOR MATRIX VECTOR MULTIPLICATION.
!     DCOPY   LEVEL 1 BLAS THAT COPIES ONE VECTOR TO ANOTHER .
!     DNRM2   LEVEL 1 BLAS THAT COMPUTES THE NORM OF A VECTOR.
!     DSCAL   LEVEL 1 BLAS THAT SCALES A VECTOR.
!
! INTRINSIC FUNCTIONS
!
!     ABS
!
! AUTHOR
!     DANNY SORENSEN               PHUONG VU
!     RICHARD LEHOUCQ              CRPC / RICE UNIVERSITY
!     DEPT. OF COMPUTATIONAL &     HOUSTON, TEXAS
!     APPLIED MATHEMATICS
!     RICE UNIVERSITY
!     HOUSTON, TEXAS
!
! REVISION HISTORY:
!     XX/XX/92: VERSION ' 2.1'
!
! FILE: NEIGH.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
! ASTER INFORMATION
! 07/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            DISPARITION DE SECOND,
!            COMMON TIMING REMPLACE PAR COMMON INFOR,
!            MODIFICATION DES APPELS BLAS (ROUTINE ASTER BL...),
!            IMPLICIT NONE.
!
! ENDLIB
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
#include "asterfort/as_dtrevc.h"
#include "asterfort/as_dlaqrb.h"
#include "asterfort/as_dmout.h"
#include "asterfort/as_dvout.h"
#include "blas/dgemv.h"
#include "blas/dlacpy.h"
#include "blas/dlapy2.h"
#include "blas/dnrm2.h"
#include "blas/dscal.h"
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
    real(kind=8) :: rnorm
!
!     %-----------------%
!     | ARRAY ARGUMENTS |
!     %-----------------%
!
    real(kind=8) :: bounds(n), h(ldh, n), q(ldq, n), ritzi(n), ritzr(n)
    real(kind=8) :: workl(n*(n+3))
!
!     %------------%
!     | PARAMETERS |
!     %------------%
!
    real(kind=8) :: one, zero
    parameter(one=1.0d+0, zero=0.0d+0)
!
!     %------------------------%
!     | LOCAL SCALARS & ARRAYS |
!     %------------------------%
!
    aster_logical :: select(1)
    integer(kind=8) :: i, iconj, msglvl
    real(kind=8) :: temp, vl(1)
    blas_int :: b_incx, b_incy, b_lda, b_m, b_n
    blas_int :: b_ldb
!
!     %-----------%
!     | FUNCTIONS |
!     %-----------%
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
        call as_dmout(logfil, n, n, h, ldh, &
                   ndigit, '_NEIGH: ENTERING UPPER HESSENBERG MATRIX H ')
    end if
!
!     %-----------------------------------------------------------%
!     | 1. COMPUTE THE EIGENVALUES, THE LAST COMPONENTS OF THE    |
!     |    CORRESPONDING SCHUR VECTORS AND THE FULL SCHUR FORM T  |
!     |    OF THE CURRENT UPPER HESSENBERG MATRIX H.              |
!     | DLAQRB RETURNS THE FULL SCHUR FORM OF H IN WORKL(1:N**2)  |
!     | AND THE LAST COMPONENTS OF THE SCHUR VECTORS IN BOUNDS.   |
!     %-----------------------------------------------------------%
! DUE TO CRP_102 CALL DLACPY ('ALL', N, N, H, LDH, WORKL, N)
    b_ldb = to_blas_int(n)
    b_lda = to_blas_int(ldh)
    b_m = to_blas_int(n)
    b_n = to_blas_int(n)
    call dlacpy('A', b_m, b_n, h, b_lda, &
                workl, b_ldb)
    call as_dlaqrb(.true._1, n, 1, n, workl, &
                n, ritzr, ritzi, bounds, ierr)
    if (ierr .ne. 0) goto 9000
!
    if (msglvl .gt. 1) then
        call as_dvout(logfil, n, bounds, ndigit, '_NEIGH: LAST ROW OF THE SCHUR MATRIX FOR H')
    end if
!
!     %-----------------------------------------------------------%
!     | 2. COMPUTE THE EIGENVECTORS OF THE FULL SCHUR FORM T AND  |
!     |    APPLY THE LAST COMPONENTS OF THE SCHUR VECTORS TO GET  |
!     |    THE LAST COMPONENTS OF THE CORRESPONDING EIGENVECTORS. |
!     | REMEMBER THAT IF THE I-TH AND (I+1)-ST EIGENVALUES ARE    |
!     | COMPLEX CONJUGATE PAIRS, THEN THE REAL & IMAGINARY PART   |
!     | OF THE EIGENVECTOR COMPONENTS ARE SPLIT ACROSS ADJACENT   |
!     | COLUMNS OF Q.                                             |
!     %-----------------------------------------------------------%
!
    call as_dtrevc('R', 'A', select, n, workl, &
                   n, vl, n, q, ldq, &
                   n, n, workl(n*n+1), ierr)
!
    if (ierr .ne. 0) goto 9000
!
!     %------------------------------------------------%
!     | SCALE THE RETURNING EIGENVECTORS SO THAT THEIR |
!     | EUCLIDEAN NORMS ARE ALL ONE. LAPACK SUBROUTINE |
!     | DTREVC RETURNS EACH EIGENVECTOR NORMALIZED SO  |
!     | THAT THE ELEMENT OF LARGEST MAGNITUDE HAS      |
!     | MAGNITUDE 1, HERE THE MAGNITUDE OF A COMPLEX   |
!     | NUMBER (X,Y) IS TAKEN TO BE |X| + |Y|.         |
!     %------------------------------------------------%
!
    iconj = 0
    do i = 1, n
        if (abs(ritzi(i)) .le. zero) then
!
!           %----------------------%
!           | REAL EIGENVALUE CASE |
!           %----------------------%
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            temp = dnrm2(b_n, q(1, i), b_incx)
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            call dscal(b_n, one/temp, q(1, i), b_incx)
        else
!
!           %-------------------------------------------%
!           | COMPLEX CONJUGATE PAIR CASE. NOTE THAT    |
!           | SINCE THE REAL AND IMAGINARY PART OF      |
!           | THE EIGENVECTOR ARE STORED IN CONSECUTIVE |
!           | COLUMNS, WE FURTHER NORMALIZE BY THE      |
!           | SQUARE ROOT OF TWO.                       |
!           %-------------------------------------------%
!
            if (iconj .eq. 0) then
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                temp = dlapy2(dnrm2(b_n, q(1, i), b_incx), dnrm2(b_n, q(1, i+1), b_incx))
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                call dscal(b_n, one/temp, q(1, i), b_incx)
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                call dscal(b_n, one/temp, q(1, i+1), b_incx)
                iconj = 1
            else
                iconj = 0
            end if
        end if
    end do
!
    b_lda = to_blas_int(ldq)
    b_m = to_blas_int(n)
    b_n = to_blas_int(n)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dgemv('T', b_m, b_n, one, q, &
               b_lda, bounds, b_incx, zero, workl, &
               b_incy)
!
    if (msglvl .gt. 1) then
        call as_dvout(logfil, n, workl, ndigit, '_NEIGH: LAST ROW OF THE EIGENVECTOR MATRIX FOR H')
    end if
!
!     %----------------------------%
!     | COMPUTE THE RITZ ESTIMATES |
!     %----------------------------%
!
    iconj = 0
    do i = 1, n
        if (abs(ritzi(i)) .le. zero) then
!
!           %----------------------%
!           | REAL EIGENVALUE CASE |
!           %----------------------%
!
            bounds(i) = rnorm*abs(workl(i))
        else
!
!           %-------------------------------------------%
!           | COMPLEX CONJUGATE PAIR CASE. NOTE THAT    |
!           | SINCE THE REAL AND IMAGINARY PART OF      |
!           | THE EIGENVECTOR ARE STORED IN CONSECUTIVE |
!           | COLUMNS, WE NEED TO TAKE THE MAGNITUDE    |
!           | OF THE LAST COMPONENTS OF THE TWO VECTORS |
!           %-------------------------------------------%
!
            if (iconj .eq. 0) then
                bounds(i) = rnorm*dlapy2(workl(i), workl(i+1))
                bounds(i+1) = bounds(i)
                iconj = 1
            else
                iconj = 0
            end if
        end if
    end do
!
    if (msglvl .gt. 2) then
        call as_dvout(logfil, n, ritzr, ndigit, '_NEIGH: REAL PART OF THE EIGENVALUES OF H')
        call as_dvout(logfil, n, ritzi, ndigit, '_NEIGH: IMAGINARY PART OF THE EIGENVALUES OF H')
        call as_dvout(logfil, n, bounds, ndigit, '_NEIGH: RITZ ESTIMATES FOR THE EIGENVALUES OF H')
    end if
!
9000 continue
!
    call matfpe(1)
!
!     %---------------%
!     | END OF DNEIGH |
!     %---------------%
!
end subroutine
