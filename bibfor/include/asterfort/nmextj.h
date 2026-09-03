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
!
interface
    subroutine nmextj(field_type, nb_elem  , nb_cmp  , list_cmp, type_extr_cmp, type_sele_cmp,&
                      poin_nume , spoi_nume, nb_vale , i_elem       , elem_nume    ,&
                      jcesd     , jcesv    , jcesl   , jcesc        , vale_resu)
        character(len=24), intent(in) :: field_type
        integer(kind=8), intent(in) :: nb_elem
        integer(kind=8), intent(in) :: nb_cmp
        character(len=24), intent(in) :: list_cmp
        character(len=8), intent(in) :: type_extr_cmp
        character(len=8), intent(in) :: type_sele_cmp
        integer(kind=8), intent(in) :: poin_nume
        integer(kind=8), intent(in):: i_elem
        integer(kind=8), intent(in):: elem_nume
        integer(kind=8), intent(in) :: spoi_nume
        integer(kind=8), intent(in) :: jcesd
        integer(kind=8), intent(in) :: jcesv
        integer(kind=8), intent(in) :: jcesl
        integer(kind=8), intent(in) :: jcesc
        integer(kind=8), intent(out) :: nb_vale
        real(kind=8), intent(out) :: vale_resu(*)
    end subroutine nmextj
end interface
