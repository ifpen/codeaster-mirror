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
subroutine te0519(option, nomte)
!
    use coupling_nitsche_module
    use coupling_quadrature_module
    use coupling_type
    use FE_quadrature_module
    use FE_topo_module
    use HHO_matrix_module
    use HHO_quadrature_module
    use HHO_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/coupling_type.h"
#include "asterfort/HHO_size_module.h"
#include "FE_module.h"
!
    character(len=16), intent(in) :: option, nomte
!
!===================================================================================================
!
!
! Mechanics - Coupling FEM/HHO with Nitsche
!
!===================================================================================================
!
    type(FE_Cell):: FECellSl
    type(FE_Skin):: FEFaceSl
    type(HHO_Face):: hhoFaceMa
    type(HHO_Data):: hhoDataMa
    type(FE_Quadrature) :: FEQuadSl
    type(HHO_Quadrature) :: hhoQuadMa
    type(CouplingMap) :: cplNitsMap
    type(CouplingData) :: cplNitsData
    type(HHO_matrix) :: lhs
!
! --- Initialize topologie
!
    call cplNitsInitTopoFEMHHO(FEFaceSl, FECellSl, hhoFaceMa, hhoDataMa)
!
! --- Initialize mapping
!
    call cplNitsInitMapFEMHHOVec(FECellSl, hhoFaceMa, hhoDataMa, cplNitsMap)
!
! --- Initialize quadrature
!
    call cplGetQuadFEMHHO(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa)
!
! --- Initialize data
!
    call cplNitsInitData(option, FECellSl, cplNitsMap, cplNitsData)
!
! - Computation
!
    if (option == "RIGI_ELAS_CPL") then
!
! --- Compute coupling matrix
!
        call cplNitsLhsFEMHHO(FEFaceSl, FECellSl, hhoFaceMa, hhoDataMa, &
                              FEQuadSl, hhoQuadMa, cplNitsMap, cplNitsData, lhs)
!
! - Write matrix
!
        call lhs%write("PMATUUR", ASTER_TRUE)
        call lhs%free()
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
