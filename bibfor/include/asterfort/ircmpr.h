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
interface
    subroutine ircmpr(nofimd, typech, nbimpr, ncaimi, ncaimk,&
                      ncmprf, ncmpve, ntlcmp, nbvato, nbenec,&
                      lienec, adsd, adsl, nomaas, ligrel,&
                      typgeo, nomtyp, ntproa, chanom, sdcarm,&
                      field_type, nosdfu, nosdpr)
        character(len=*) :: nofimd
        character(len=8) :: typech
        integer(kind=8) :: nbimpr
        character(len=24) :: ncaimi
        character(len=24) :: ncaimk
        integer(kind=8) :: ncmprf
        integer(kind=8) :: ncmpve
        character(len=*) :: ntlcmp
        integer(kind=8) :: nbvato
        integer(kind=8) :: nbenec
        integer(kind=8) :: lienec(*)
        integer(kind=8) :: adsd
        integer(kind=8) :: adsl
        character(len=8) :: nomaas
        character(len=19) :: ligrel
        integer(kind=8) :: typgeo(*)
        character(len=8) :: nomtyp(*)
        character(len=*) :: ntproa
        character(len=19) :: chanom
        character(len=8) :: sdcarm
        character(len=16), intent(in) :: field_type
        character(len=8) :: nosdfu
        character(len=*) :: nosdpr
    end subroutine ircmpr
end interface
