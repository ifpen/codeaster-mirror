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
    subroutine mat_proto(BEHinteg, &
                         fami, kpg, ksp, poum, jvMaterCode, relaComp, &
                         nbProp, propVale)
        use Behaviour_type
        type(Behaviour_Integ), intent(in) :: BEHinteg
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp
        character(len=1), intent(in) :: poum
        integer(kind=8), intent(in) :: jvMaterCode
        character(len=*), intent(in) :: relaComp
        integer(kind=8), intent(inout) :: nbProp
        real(kind=8), intent(out) :: propVale(MGIS_MAX_PROPS)
    end subroutine mat_proto
end interface
