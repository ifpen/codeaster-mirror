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
    subroutine mfront_get_mater_value(extern_addr, BEHInteg, relaComp, &
                                      fami, kpg, ksp, jvMaterCode, &
                                      propVale, nbProp)
        use Behaviour_type
        character(len=16), intent(in) :: extern_addr
        type(Behaviour_Integ), intent(in) :: BEHInteg
        character(len=16), intent(in) :: relaComp
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode, nbProp
        real(kind=8), intent(out) :: propVale(nbProp)
    end subroutine mfront_get_mater_value
end interface
