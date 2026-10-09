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
#include "asterf_types.h"
#include "MeshTypes_type.h"
!
interface
    subroutine xcinem(axi, igeom, nnop, nnos, idepl, ndim, he,&
                      nfiss, nfh, nfe, ddls, ddlm,&
                      fk, dkdgl, ff, dfdi, f, eps, grad, heavn)
        aster_logical, intent(in) :: axi
        integer(kind=8), intent(in) :: igeom
        integer(kind=8), intent(in) :: nnop
        integer(kind=8), intent(in) :: nnos
        integer(kind=8), intent(in) :: idepl
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(in) :: he(nfiss)
        integer(kind=8), intent(in) :: nfiss
        integer(kind=8), intent(in) :: nfh
        integer(kind=8), intent(in) :: nfe
        integer(kind=8), intent(in) :: ddls
        integer(kind=8), intent(in) :: ddlm
        real(kind=8), intent(in)::  fk(MT_NNOMAX,3,3)
        real(kind=8), intent(in)::  dkdgl(MT_NNOMAX,3,3,3)
        real(kind=8), intent(in) :: ff(nnop)
        real(kind=8), intent(in) :: dfdi(nnop, ndim)
        real(kind=8), intent(out) :: f(3, 3)
        real(kind=8), intent(out) :: eps(6)
        real(kind=8), intent(out) :: grad(ndim, ndim)
        integer(kind=8), intent(in) :: heavn(nnop, 5)
    end subroutine xcinem
end interface
