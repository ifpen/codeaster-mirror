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

subroutine te0527(nomopt, nomte)
!
    use coupling_computation_module
    use coupling_penalisation_module
    use coupling_quadrature_module
    use coupling_nitsche_module
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
    character(len=16) :: nomte, nomopt
!
!===================================================================================================
!
!  Option: HHO_CALC_PROJ
!  Coupling FEM/HHO with Nitsche or Penalization
!  Compute RHS of L2-projector from FEM to HHO
!
!===================================================================================================
!
    type(FE_Cell):: FECellSl
    type(FE_Skin):: FEFaceSl
    type(HHO_Face):: hhoFaceMa
    type(HHO_Data):: hhoDataMa
    type(FE_Quadrature) :: FEQuadSl
    type(HHO_Quadrature) :: hhoQuadMa
    type(CouplingMap) :: cplMap
    type(HHO_matrix) :: lhs, projFEM2HHO
    integer(kind=8) :: iFe, irow, jFE, jcol
!
    if (nomte(1:2) == 'CP') then
!
! --- Initialize topologie
!
        call cplPenaInitTopoFEMHHO(FEFaceSl, hhoFaceMa, hhoDataMa)
!
! --- Initialize mapping
!
        call cplPenaInitMapFEMHHOVec(FEFaceSl, hhoFaceMa, hhoDataMa, cplMap)
!
    elseif (nomte(1:2) == 'CN') then
!
! --- Initialize topologie
!
        call cplNitsInitTopoFEMHHO(FEFaceSl, FECellSl, hhoFaceMa, hhoDataMa)
!
! --- Initialize mapping
!
        call cplNitsInitMapFEMHHOVec(FECellSl, hhoFaceMa, hhoDataMa, cplMap)
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
! --- Initialize quadrature
!
    call cplGetQuadFEMHHO(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa)
!
! ----- Compute RHS of Pi_F^k(phi_i^CG) on the whole HHO face F.
!       This is only the local part of the projector since a cG node
!       can be shared by two face. To compute the non-local operator,
!       it is necessary to assemble them together.

    call cplCmpProjFEM2HHOFaceLocVec(FEFaceSl, hhoFaceMa, hhoDataMa, &
                                     FEQuadSl, hhoQuadMa, projFEM2HHO)
!
!  Assemble matrix
!
    call lhs%initialize(cplMap%nbDoFs, cplMap%nbDoFs, 0.d0)
!
    do jFE = 1, cplMap%nbDoFsFEFaceSl
        jcol = cplMap%mapDoFsFEFaceSl(jFE)
        do iFE = 1, cplMap%nbDoFshhoFaceMa
            irow = cplMap%mapDoFshhoFaceMa(iFE)
            lhs%m(irow, jcol) = projFEM2HHO%m(iFE, jFE)
        end do
    end do
    call projFEM2HHO%free()
!
! - Write matrix
!
    call lhs%write("PMATUNS", ASTER_FALSE)
    call lhs%free()
!
end subroutine
