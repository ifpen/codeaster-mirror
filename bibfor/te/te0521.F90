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
subroutine te0521(option, nomte)
!
    use coupling_penalisation_module
    use coupling_quadrature_module
    use coupling_type
    use FE_quadrature_module
    use FE_topo_module
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/coupling_type.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/writeVector.h"
#include "FE_module.h"
!
    character(len=16), intent(in) :: option, nomte
!
!===================================================================================================
!
!
! Mechanics - Coupling FEM/FEM with PENALISATION
!
!===================================================================================================
!
    type(FE_Skin):: FEFaceSl, FEFaceMa
    type(FE_Quadrature) :: FEQuadSl, FEQuadMa

    type(CouplingMap) :: cplPenaMap
    type(CouplingData) :: cplPenaData
    real(kind=8) :: rhs(MSIZE_CPL_PENA_FEM)
    real(kind=8) :: lhs(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
!
! --- Initialize topologie
!
    call cplPenaInitTopoFEMFEM(FEFaceSl, FEFaceMa)
!
! --- Initialize mapping
!
    call cplPenaInitMapFEMFEMVec(FEFaceSl, FEFaceMa, cplPenaMap)
!
! --- Initialize quadrature
!
    call cplGetQuadFEMFEM(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa)
!
! --- Initialize data
!
    call cplPenaInitData(option, cplPenaMap, cplPenaData)
!
! - Computation
!
    if (option == "CHAR_MECA_CPL") then
!
! --- Compute coupling residual
!
        call cplPenaRhsFEMFEM(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, &
                              cplPenaData, rhs)
!
! --- Write vector
!
        call writeVector('PVECTUR', cplPenaMap%nbDoFs, rhs)
!
    elseif (option == "RIGI_CPL" .or. option == "RIGI_ELAS_CPL") then
!
! --- Compute coupling matrix
!
        call cplPenaLhsFEMFEM(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, &
                              cplPenaData, lhs)
!
! - Write matrix
!
        call writeMatrix("PMATUUR", cplPenaMap%nbDoFs, cplPenaMap%nbDoFs, &
                         ASTER_TRUE, lhs)
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
