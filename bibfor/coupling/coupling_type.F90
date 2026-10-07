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
module coupling_type
!
    implicit none
!
    private
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/coupling_type.h"
#include "asterfort/HHO_size_module.h"
#include "FE_module.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! Coupling - Generic
!
! --------------------------------------------------------------------------------------------------
!
!
    type CouplingMap
!
        integer(kind=8) :: nbDoFs = 0, nbDoFsFECellSl = 0, nbDoFsFEFaceSl = 0
        integer(kind=8) :: nbDoFsFEFaceLagSl = 0, nbDoFsFace = 0
        integer(kind=8) :: nbDoFshhoFaceMa = 0, nbDoFsFEFaceMa = 0
        integer(kind=8) :: mapDoFsFECellSl(MAX_BV_CG) = 0
        integer(kind=8) :: mapDoFsFEFaceSl(MAX_BV_CG) = 0
        integer(kind=8) :: mapDoFsFEFaceLagSl(MAX_BV_CG) = 0
        integer(kind=8) :: mapDoFsFEFaceMa(MAX_BV_CG) = 0
        integer(kind=8) :: mapDoFsFace(2*MAX_BV_CG) = 0
        integer(kind=8) :: mapDoFshhoFaceMa(MSIZE_FACE_VEC) = 0
        integer(kind=8) :: face2cellSl(MT_NNOMAX2D) = 0
    end type
!
    type CouplingData
!
        integer(kind=8) :: nbDoFs = 0
        real(kind=8) :: disp_prev(MSIZE_FACE_VEC+MAX_BV_CG) = 0.d0
        real(kind=8) :: dispFECellSl_prev(MAX_BV_CG) = 0.d0
        real(kind=8) :: dispFEFaceSl_prev(MAX_BV_CG) = 0.d0
        real(kind=8) :: disphhoFaceMa_prev(MSIZE_FACE_VEC) = 0.d0
!
        real(kind=8) :: disp_curr(MSIZE_FACE_VEC+MAX_BV_CG) = 0.d0
        real(kind=8) :: dispFECellSl_curr(MAX_BV_CG) = 0.d0
        real(kind=8) :: dispFEFaceSl_curr(MAX_BV_CG) = 0.d0
        real(kind=8) :: disphhoFaceMa_curr(MSIZE_FACE_VEC) = 0.d0
!
        real(kind=8) :: coef_pena = 0.d0
!
        real(kind=8), dimension(MT_NNOMAX3D) :: E = 0.d0, nu = 0.d0
    end type
!
    public :: CouplingMap, CouplingData
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
end module
