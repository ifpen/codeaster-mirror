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

subroutine digdel_wrap(modelo, nb_cmp)
    implicit none
!
!     ARGUMENTS:
!     ----------
!
#include "asterfort/digdel.h"
    integer(kind=8), intent(in) :: modelo
    integer(kind=8), intent(out) :: nb_cmp
! ----------------------------------------------------------------------
!     ENTREES:
!        MODELO : MODE_LOCAL (SON INDICE DANS &CATA.TE.MODELOC )
!
!     SORTIES:
!        nb_cmp : NOMBRE DE SCALAIRES REPRESENTANT LA GRANDEUR POUR LE
!                 MODE_LOCAL
!
! ---------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    nb_cmp = digdel(modelo)
end subroutine
