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
subroutine as_dmout(lout, m, n, a, lda, &
                    idigit, ifmt)
!
!     SUBROUTINE ARPACK ECRIVANT DES MATRICES.
!-----------------------------------------------------------------------
!  PURPOSE:    REAL MATRIX OUTPUT ROUTINE.
!
!  USAGE:      CALL DMOUT (LOUT, M, N, A, LDA, IDIGIT, IFMT)
!
!  ARGUMENTS
!     M      - NUMBER OF ROWS OF A.  (INPUT)
!     N      - NUMBER OF COLUMNS OF A.  (INPUT)
!     A      - REAL M BY N MATRIX TO BE PRINTED.  (INPUT)
!     LDA    - LEADING DIMENSION OF A EXACTLY AS SPECIFIED IN THE
!              DIMENSION STATEMENT OF THE CALLING PROGRAM.  (INPUT)
!     IFMT   - FORMAT TO BE USED IN PRINTING MATRIX A.  (INPUT)
!     IDIGIT - PRINT UP TO IABS(IDIGIT) DECIMAL DIGITS PER NUMBER.  (IN)
!              IF IDIGIT .LT. 0, PRINTING IS DONE WITH 72 COLUMNS.
!              IF IDIGIT .GT. 0, PRINTING IS DONE WITH 132 COLUMNS.
!  INTRINSIC FUNCTIONS
!     MIN, LEN.
!-----------------------------------------------------------------------
    implicit none
!
!     .. SCALAR ARGUMENTS ..
    character(len=*) :: ifmt
    integer(kind=8) :: idigit, lda, lout, m, n
!
!     .. ARRAY ARGUMENTS ..
    real(kind=8) :: a(lda, *)
!
!     .. LOCAL SCALARS ..
    character(len=80) :: line
    integer(kind=8) :: i, j, k1, k2, lll, ndigit
!
!     .. LOCAL ARRAYS ..
    character(len=1) :: icol(3)
!
!     .. DATA STATEMENTS ..
    data icol(1), icol(2), icol(3)/'C', 'O',&
     &                   'L'/
!
!     .. EXECUTABLE STATEMENTS ..
!
!     ... FIRST EXECUTABLE STATEMENT
!
    lll = min(len(ifmt), 80)
    do i = 1, lll
        line(i:i) = '-'
    end do
!
    do i = lll+1, 80
        line(i:i) = ' '
    end do
!
    write (lout, fmt=999) ifmt, line(1:lll)
999 format(/1x, a, /1x, a)
!
    if (m .le. 0 .or. n .le. 0 .or. lda .le. 0) goto 100
    ndigit = idigit
    if (idigit .eq. 0) ndigit = 4
!
!=======================================================================
!             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
!=======================================================================
!
    if (idigit .lt. 0) then
        ndigit = -idigit
        if (ndigit .le. 4) then
            do k1 = 1, n, 5
                k2 = min(n, k1+4)
                write (lout, fmt=998) (icol, i, i=k1, k2)
                do i = 1, m
                    write (lout, fmt=994) i, (a(i, j), j=k1, &
                                              k2)
                end do
            end do
!
        else if (ndigit .le. 6) then
            do k1 = 1, n, 4
                k2 = min(n, k1+3)
                write (lout, fmt=997) (icol, i, i=k1, k2)
                do i = 1, m
                    write (lout, fmt=993) i, (a(i, j), j=k1, &
                                              k2)
                end do
            end do
!
        else if (ndigit .le. 10) then
            do k1 = 1, n, 3
                k2 = min(n, k1+2)
                write (lout, fmt=996) (icol, i, i=k1, k2)
                do i = 1, m
                    write (lout, fmt=992) i, (a(i, j), j=k1, &
                                              k2)
                end do
            end do
!
        else
            do k1 = 1, n, 2
                k2 = min(n, k1+1)
                write (lout, fmt=995) (icol, i, i=k1, k2)
                do i = 1, m
                    write (lout, fmt=991) i, (a(i, j), j=k1, &
                                              k2)
                end do
            end do
        end if
!
!=======================================================================
!             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
!=======================================================================
!
    else
        if (ndigit .le. 4) then
            do k1 = 1, n, 10
                k2 = min(n, k1+9)
                write (lout, fmt=998) (icol, i, i=k1, k2)
                do i = 1, m
                    write (lout, fmt=994) i, (a(i, j), j=k1, &
                                              k2)
                end do
            end do
!
        else if (ndigit .le. 6) then
            do k1 = 1, n, 8
                k2 = min(n, k1+7)
                write (lout, fmt=997) (icol, i, i=k1, k2)
                do i = 1, m
                    write (lout, fmt=993) i, (a(i, j), j=k1, &
                                              k2)
                end do
            end do
!
        else if (ndigit .le. 10) then
            do k1 = 1, n, 6
                k2 = min(n, k1+5)
                write (lout, fmt=996) (icol, i, i=k1, k2)
                do i = 1, m
                    write (lout, fmt=992) i, (a(i, j), j=k1, &
                                              k2)
                end do
            end do
!
        else
            do k1 = 1, n, 5
                k2 = min(n, k1+4)
                write (lout, fmt=995) (icol, i, i=k1, k2)
                do i = 1, m
                    write (lout, fmt=991) i, (a(i, j), j=k1, &
                                              k2)
                end do
            end do
        end if
    end if
    write (lout, fmt=990)
!
998 format(10x, 10(4x, 3a1, i4, 1x))
997 format(10x, 8(5x, 3a1, i4, 2x))
996 format(10x, 6(7x, 3a1, i4, 4x))
995 format(10x, 5(9x, 3a1, i4, 6x))
994 format(1x, ' ROW', i4, ':', 1x, 1p, 10d12.3)
993 format(1x, ' ROW', i4, ':', 1x, 1p, 8d14.5)
992 format(1x, ' ROW', i4, ':', 1x, 1p, 6d18.9)
991 format(1x, ' ROW', i4, ':', 1x, 1p, 5d22.13)
990 format(1x, ' ')
!
100 continue
end subroutine
