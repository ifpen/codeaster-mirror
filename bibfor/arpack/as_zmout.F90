! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
subroutine as_zmout(lout, m, n, a, lda, &
                 idigit, ifmt)
!
!     SUBROUTINE ARPACK ECRIVANT DES MATRICES COMPLEXES.
!-----------------------------------------------------------------------
!
!  ROUTINE:    ZMOUT
!
!  PURPOSE:    COMPLEX*16 MATRIX OUTPUT ROUTINE.
!
!  USAGE:      CALL ZMOUT (LOUT, M, N, A, LDA, IDIGIT, IFMT)
!
!  ARGUMENTS
!     M      - NUMBER OF ROWS OF A.  (INPUT)
!     N      - NUMBER OF COLUMNS OF A.  (INPUT)
!     A      - COMPLEX*16 M BY N MATRIX TO BE PRINTED.  (INPUT)
!     LDA    - LEADING DIMENSION OF A EXACTLY AS SPECIFIED IN THE
!              DIMENSION STATEMENT OF THE CALLING PROGRAM.  (INPUT)
!     IFMT   - FORMAT TO BE USED IN PRINTING MATRIX A.  (INPUT)
!     IDIGIT - PRINT UP TO IABS(IDIGIT) DECIMAL DIGITS PER NUMBER.  (IN)
!              IF IDIGIT .LT. 0, PRINTING IS DONE WITH 72 COLUMNS.
!              IF IDIGIT .GT. 0, PRINTING IS DONE WITH 132 COLUMNS.
!
!\SCCS INFORMATION: @(#)
! FILE: ZMOUT.F   SID: 2.1   DATE OF SID: 11/16/95   RELEASE: 2
!
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
!     ... SPECIFICATIONS FOR ARGUMENTS
    integer(kind=8) :: m, n, idigit, lda, lout
    complex(kind=8) :: a(lda, *)
    character(len=*) :: ifmt
!
!     ... SPECIFICATIONS FOR LOCAL VARIABLES
    integer(kind=8) :: i, j, ndigit, k1, k2, lll
    character(len=1) :: icol(3)
    character(len=80) :: line
!     ...
!
    data icol(1), icol(2), icol(3)/'C', 'O',&
     &                   'L'/
!     ...
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
    write (lout, 9999) ifmt, line(1:lll)
9999 format(/1x, a/1x, a)
!
    if (m .le. 0 .or. n .le. 0 .or. lda .le. 0) goto 1000
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
            do k1 = 1, n, 2
                k2 = min(n, k1+1)
                write (lout, 9998) (icol, i, i=k1, k2)
                do i = 1, m
                    if (k1 .ne. n) then
                        write (lout, 9994) i, (a(i, j), j=k1, k2 &
                                               )
                    else
                        write (lout, 9984) i, (a(i, j), j=k1, k2 &
                                               )
                    end if
                end do
            end do
!
        else if (ndigit .le. 6) then
            do k1 = 1, n, 2
                k2 = min(n, k1+1)
                write (lout, 9997) (icol, i, i=k1, k2)
                do i = 1, m
                    if (k1 .ne. n) then
                        write (lout, 9993) i, (a(i, j), j=k1, k2 &
                                               )
                    else
                        write (lout, 9983) i, (a(i, j), j=k1, k2 &
                                               )
                    end if
                end do
            end do
!
        else if (ndigit .le. 8) then
            do k1 = 1, n, 2
                k2 = min(n, k1+1)
                write (lout, 9996) (icol, i, i=k1, k2)
                do i = 1, m
                    if (k1 .ne. n) then
                        write (lout, 9992) i, (a(i, j), j=k1, k2 &
                                               )
                    else
                        write (lout, 9982) i, (a(i, j), j=k1, k2 &
                                               )
                    end if
                end do
            end do
!
        else
            do k1 = 1, n
                write (lout, 9995) icol, k1
                do i = 1, m
                    write (lout, 9991) i, a(i, k1)
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
            do k1 = 1, n, 4
                k2 = min(n, k1+3)
                write (lout, 9998) (icol, i, i=k1, k2)
                do i = 1, m
                    if ((k1+3) .le. n) then
                        write (lout, 9974) i, (a(i, j), j=k1, k2 &
                                               )
                    else if ((k1+3-n) .eq. 1) then
                        write (lout, 9964) i, (a(i, j), j=k1, k2 &
                                               )
                    else if ((k1+3-n) .eq. 2) then
                        write (lout, 9954) i, (a(i, j), j=k1, k2 &
                                               )
                    else if ((k1+3-n) .eq. 3) then
                        write (lout, 9944) i, (a(i, j), j=k1, k2 &
                                               )
                    end if
                end do
            end do
!
        else if (ndigit .le. 6) then
            do k1 = 1, n, 3
                k2 = min(n, k1+2)
                write (lout, 9997) (icol, i, i=k1, k2)
                do i = 1, m
                    if ((k1+2) .le. n) then
                        write (lout, 9973) i, (a(i, j), j=k1, k2 &
                                               )
                    else if ((k1+2-n) .eq. 1) then
                        write (lout, 9963) i, (a(i, j), j=k1, k2 &
                                               )
                    else if ((k1+2-n) .eq. 2) then
                        write (lout, 9953) i, (a(i, j), j=k1, k2 &
                                               )
                    end if
                end do
            end do
!
        else if (ndigit .le. 8) then
            do k1 = 1, n, 3
                k2 = min(n, k1+2)
                write (lout, 9996) (icol, i, i=k1, k2)
                do i = 1, m
                    if ((k1+2) .le. n) then
                        write (lout, 9972) i, (a(i, j), j=k1, k2 &
                                               )
                    else if ((k1+2-n) .eq. 1) then
                        write (lout, 9962) i, (a(i, j), j=k1, k2 &
                                               )
                    else if ((k1+2-n) .eq. 2) then
                        write (lout, 9952) i, (a(i, j), j=k1, k2 &
                                               )
                    end if
                end do
            end do
!
        else
            do k1 = 1, n, 2
                k2 = min(n, k1+1)
                write (lout, 9995) (icol, i, i=k1, k2)
                do i = 1, m
                    if ((k1+1) .le. n) then
                        write (lout, 9971) i, (a(i, j), j=k1, k2 &
                                               )
                    else
                        write (lout, 9961) i, (a(i, j), j=k1, k2 &
                                               )
                    end if
                end do
            end do
        end if
    end if
    write (lout, 9990)
!
9998 format(11x, 4(9x, 3a1, i4, 9x))
9997 format(10x, 4(11x, 3a1, i4, 11x))
9996 format(10x, 3(13x, 3a1, i4, 13x))
9995 format(12x, 2(18x, 3a1, i4, 18x))
!
!========================================================
!              FORMAT FOR 72 COLUMN
!========================================================
!
!            DISPLAY 4 SIGNIFICANT DIGITS
!
9994 format(1x, ' ROW', i4, ':', 1x, 1p, 2('(', d10.3, ',', d10.3, ')  '))
9984 format(1x, ' ROW', i4, ':', 1x, 1p, 1('(', d10.3, ',', d10.3, ')  '))
!
!            DISPLAY 6 SIGNIFICANT DIGITS
!
9993 format(1x, ' ROW', i4, ':', 1x, 1p, 2('(', d12.5, ',', d12.5, ')  '))
9983 format(1x, ' ROW', i4, ':', 1x, 1p, 1('(', d12.5, ',', d12.5, ')  '))
!
!            DISPLAY 8 SIGNIFICANT DIGITS
!
9992 format(1x, ' ROW', i4, ':', 1x, 1p, 2('(', d14.7, ',', d14.7, ')  '))
9982 format(1x, ' ROW', i4, ':', 1x, 1p, 1('(', d14.7, ',', d14.7, ')  '))
!
!            DISPLAY 13 SIGNIFICANT DIGITS
!
9991 format(1x, ' ROW', i4, ':', 1x, 1p, 1('(', d20.13, ',', d20.13, ')'))
9990 format(1x, ' ')
!
!
!========================================================
!              FORMAT FOR 132 COLUMN
!========================================================
!
!            DISPLAY 4 SIGNIFICANT DIGIT
!
9974 format(1x, ' ROW', i4, ':', 1x, 1p, 4('(', d10.3, ',', d10.3, ')  '))
9964 format(1x, ' ROW', i4, ':', 1x, 1p, 3('(', d10.3, ',', d10.3, ')  '))
9954 format(1x, ' ROW', i4, ':', 1x, 1p, 2('(', d10.3, ',', d10.3, ')  '))
9944 format(1x, ' ROW', i4, ':', 1x, 1p, 1('(', d10.3, ',', d10.3, ')  '))
!
!            DISPLAY 6 SIGNIFICANT DIGIT
!
9973 format(1x, ' ROW', i4, ':', 1x, 1p, 3('(', d12.5, ',', d12.5, ')  '))
9963 format(1x, ' ROW', i4, ':', 1x, 1p, 2('(', d12.5, ',', d12.5, ')  '))
9953 format(1x, ' ROW', i4, ':', 1x, 1p, 1('(', d12.5, ',', d12.5, ')  '))
!
!            DISPLAY 8 SIGNIFICANT DIGIT
!
9972 format(1x, ' ROW', i4, ':', 1x, 1p, 3('(', d14.7, ',', d14.7, ')  '))
9962 format(1x, ' ROW', i4, ':', 1x, 1p, 2('(', d14.7, ',', d14.7, ')  '))
9952 format(1x, ' ROW', i4, ':', 1x, 1p, 1('(', d14.7, ',', d14.7, ')  '))
!
!            DISPLAY 13 SIGNIFICANT DIGIT
!
9971 format(1x, ' ROW', i4, ':', 1x, 1p, 2('(', d20.13, ',', d20.13,&
&        ')  '))
9961 format(1x, ' ROW', i4, ':', 1x, 1p, 1('(', d20.13, ',', d20.13,&
&        ')  '))
!
!
!
1000 continue
end subroutine
