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

subroutine te0547(nomopt, nomte)
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
!  Option: HHO_EXTR_PROJ
!  Coupling FEM/HHO with Nitsche or Penalization
!  Extract scalar projection matrix.
!
!===================================================================================================
!
    type(FE_Cell):: FECellSl
    type(FE_Skin):: FEFaceSl
    type(HHO_Face):: hhoFaceMa
    type(HHO_Data):: hhoDataMa
    type(CouplingMap) :: cplMap
    type(HHO_matrix) :: lhs, projFEM2HHO
    integer(kind=8) :: iFe, irow, jFE, jcol, ndim, nrows, ncols
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
! - Read matrix
!
    call lhs%initialize(cplMap%nbDoFs, cplMap%nbDoFs, 0.d0)
    call lhs%read("PMAELNS1", ASTER_FALSE)
!
! - Extract matrix - vector valued to scalar valued
!
    ndim = hhoFaceMa%ndim+1
    nrows = cplMap%nbDoFshhoFaceMa/ndim
    ncols = cplMap%nbDoFsFEFaceSl/ndim
    call projFEM2HHO%initialize(nrows, ncols)
!
    do jFE = 1, ncols
        jcol = cplMap%mapDoFsFEFaceSl((jFE-1)*ndim+1)
        do iFE = 1, nrows
            irow = cplMap%mapDoFshhoFaceMa(iFE)
            projFEM2HHO%m(iFE, jFE) = lhs%m(irow, jcol)
        end do
    end do
!
    call projFEM2HHO%write("PMPHHOR", ASTER_FALSE)
!
    call projFEM2HHO%free()
    call lhs%free()
!
end subroutine
