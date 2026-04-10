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
subroutine as_dvout(lout, n, sx, idigit, ifmt)
!
!     SUBROUTINE ARPACK ECRIVANT DES VECTEURS DE REELS.
!-----------------------------------------------------------------------
!  PURPOSE:    REAL VECTOR OUTPUT ROUTINE.
!
!  USAGE:      CALL DVOUT (LOUT, N, SX, IDIGIT, IFMT)
!
!  ARGUMENTS
!     N      - LENGTH OF ARRAY SX.  (INPUT)
!     SX     - REAL ARRAY TO BE PRINTED.  (INPUT)
!     IFMT   - FORMAT TO BE USED IN PRINTING ARRAY SX.  (INPUT)
!     IDIGIT - PRINT UP TO IABS(IDIGIT) DECIMAL DIGITS PER NUMBER.  (IN)
!              IF IDIGIT .LT. 0, PRINTING IS DONE WITH 72 COLUMNS.
!              IF IDIGIT .GT. 0, PRINTING IS DONE WITH 132 COLUMNS.
!
!  INTRINSIC FUNCTIONS
!     LEN, MIN.
!-----------------------------------------------------------------------
    implicit none
!
!     .. SCALAR ARGUMENTS ..
    character(len=*) :: ifmt
    integer(kind=8) :: idigit, lout, n
!
!     .. ARRAY ARGUMENTS ..
    real(kind=8) :: sx(*)
!
!     .. LOCAL SCALARS ..
    character(len=80) :: line
    integer(kind=8) :: i, k1, k2, lll, ndigit
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
    if (n .le. 0) goto 100
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
                write (lout, fmt=998) k1, k2, (sx(i), i=k1, &
                                               k2)
            end do
        else if (ndigit .le. 6) then
            do k1 = 1, n, 4
                k2 = min(n, k1+3)
                write (lout, fmt=997) k1, k2, (sx(i), i=k1, &
                                               k2)
            end do
        else if (ndigit .le. 10) then
            do k1 = 1, n, 3
                k2 = min(n, k1+2)
                write (lout, fmt=996) k1, k2, (sx(i), i=k1, &
                                               k2)
            end do
        else
            do k1 = 1, n, 2
                k2 = min(n, k1+1)
                write (lout, fmt=995) k1, k2, (sx(i), i=k1, &
                                               k2)
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
                write (lout, fmt=998) k1, k2, (sx(i), i=k1, &
                                               k2)
            end do
        else if (ndigit .le. 6) then
            do k1 = 1, n, 8
                k2 = min(n, k1+7)
                write (lout, fmt=997) k1, k2, (sx(i), i=k1, &
                                               k2)
            end do
        else if (ndigit .le. 10) then
            do k1 = 1, n, 6
                k2 = min(n, k1+5)
                write (lout, fmt=996) k1, k2, (sx(i), i=k1, &
                                               k2)
            end do
        else
            do k1 = 1, n, 5
                k2 = min(n, k1+4)
                write (lout, fmt=995) k1, k2, (sx(i), i=k1, &
                                               k2)
            end do
        end if
    end if
    write (lout, fmt=994)
998 format(1x, i4, ' - ', i4, ':', 1p, 10d12.3)
997 format(1x, i4, ' - ', i4, ':', 1x, 1p, 8d14.5)
996 format(1x, i4, ' - ', i4, ':', 1x, 1p, 6d18.9)
995 format(1x, i4, ' - ', i4, ':', 1x, 1p, 5d24.13)
994 format(1x, ' ')
!
100 continue
end subroutine
