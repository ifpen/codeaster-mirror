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
module coupling_nitsche_module
!
    use coupling_computation_module
    use coupling_penalisation_module
    use coupling_quadrature_module
    use coupling_type
    use FE_algebra_module
    use FE_Basis_module
    use FE_quadrature_module
    use FE_topo_module
    use HHO_init_module
    use HHO_matrix_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
!
    implicit none
!
    private
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/coupling_type.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "FE_basis_module.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! Coupling - Nitsche methods
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: cplNitsInitTopoFEMHHO, cplNitsInitMapFEMHHOVec, cplNitsInitData
    public :: cplNitsLhsFEMHHO, cplNitsInitMapFEMHHOScal
    private :: readNitSlavMap, cplNitsAssLhsFEMHHO
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplNitsInitTopoFEMHHO(FEFaceSl, FECellSl, hhoFaceMa, hhoDataMa)
!
        implicit none
!
        type(FE_Cell), intent(out) :: FECellSl
        type(FE_Skin), intent(out) :: FEFaceSl
        type(HHO_Face), intent(out) :: hhoFaceMa
        type(HHO_Data), intent(out) :: hhoDataMa
!
!===================================================================================================
!    Initialize FEM and HHO topology
!===================================================================================================
!
        integer(kind=8) :: nbNodes, nodes(MT_NNOMAX2D)
!
        call FECellSl%init()
!
        call readNitSlavMap(nbNodes, nodes)
        FEFaceSl = FECellSl%getSkin(nbNodes, nodes)
!
        call hhoInfoInitFace(hhoFaceMa, hhoDataMa)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine readNitSlavMap(nbNodes, nodes)
!
        implicit none
!
        integer(kind=8), intent(out) :: nbNodes, nodes(MT_NNOMAX2D)
!
!===================================================================================================
!    Initialize FEM and HHO topology
!===================================================================================================
!
        integer(kind=8) :: jpair, iNode
!
! - Map Volu to Surf
!
        call jevech('PPAIRR', 'L', jpair)
        nbNodes = nint(zr(jpair-1+OFFSET_NB_NODES_NITSCHE))
!
        do iNode = 1, nbNodes
            nodes(iNode) = nint(zr(jpair-1+OFFSET_NODE_IDX_FACE+iNode-1))
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplNitsInitMapFEMHHOVec(FECellSl, hhoFaceMa, hhoDataMa, cplMap)
!
        implicit none
!
        type(FE_Cell), intent(in) :: FECellSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(CouplingMap), intent(inout) :: cplMap
!
!===================================================================================================
!    Initialize FEM and HHO topology
!===================================================================================================
!
        type(FE_basis) :: FEBasisSl
        integer(kind=8) :: nbNodes, nodes(MT_NNOMAX2D), iBase, idim, fbs, node, iNode, iDoF
!
        call readNitSlavMap(nbNodes, nodes)
        cplMap%face2cellSl(1:nbNodes) = nodes(1:nbNodes)
        call FEBasisSl%initCell(FECellSl)
        ASSERT(FEBasisSl%typeEF == EF_LAGRANGE)
!
        cplMap%nbDoFs = 0
        cplMap%nbDoFsFECellSl = 0
!
        ASSERT(FEBasisSl%size == FECellSl%nbnodes)
        do iBase = 1, FEBasisSl%size
            do idim = 1, FEBasisSl%ndim
                cplMap%nbDoFs = cplMap%nbDoFs+1
                cplMap%nbDoFsFECellSl = cplMap%nbDoFsFECellSl+1
                cplMap%mapDoFsFECellSl(cplMap%nbDoFsFECellSl) = cplMap%nbDoFs
            end do
        end do
!
        cplMap%nbDoFsFEFaceSl = 0
        cplMap%nbDoFsFace = 0
!
        do iNode = 1, nbNodes
            node = nodes(iNode)
            do idim = 1, FEBasisSl%ndim
                cplMap%nbDoFsFace = cplMap%nbDoFsFace+1
                cplMap%nbDoFsFEFaceSl = cplMap%nbDoFsFEFaceSl+1
                iDoF = cplMap%mapDoFsFECellSl(FEBasisSl%ndim*(node-1)+idim)
                cplMap%mapDoFsFEFaceSl(cplMap%nbDoFsFEFaceSl) = iDoF
                cplMap%mapDoFsFace(cplMap%nbDoFsFace) = iDoF
            end do
        end do
!
        call hhoMecaFaceDofs(hhoFaceMa, hhoDataMa, fbs)
        cplMap%nbDoFshhoFaceMa = 0
!
        do iBase = 1, fbs
            cplMap%nbDoFs = cplMap%nbDoFs+1
            cplMap%nbDoFsFace = cplMap%nbDoFsFace+1
            cplMap%nbDoFshhoFaceMa = cplMap%nbDoFshhoFaceMa+1
            cplMap%mapDoFshhoFaceMa(cplMap%nbDoFshhoFaceMa) = cplMap%nbDoFs
            cplMap%mapDoFsFace(cplMap%nbDoFsFace) = cplMap%nbDoFs
        end do
!
        ASSERT(cplMap%nbDoFs <= MSIZE_CPL_NITS)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplNitsInitMapFEMHHOScal(FECellSl, hhoFaceMa, hhoDataMa, cplMap)
!
        implicit none
!
        type(FE_Cell), intent(in) :: FECellSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(CouplingMap), intent(inout) :: cplMap
!
!===================================================================================================
!    Initialize FEM and HHO topology
!===================================================================================================
!
        type(FE_basis) :: FEBasisSl
        integer(kind=8) :: nbNodes, nodes(MT_NNOMAX2D), iBase, fbs, node, iNode, iDoF
!
        call readNitSlavMap(nbNodes, nodes)
        cplMap%face2cellSl(1:nbNodes) = nodes(1:nbNodes)
        call FEBasisSl%initCell(FECellSl)
        ASSERT(FEBasisSl%typeEF == EF_LAGRANGE)
!
        cplMap%nbDoFs = 0
        cplMap%nbDoFsFECellSl = 0
!
        ASSERT(FEBasisSl%size == FECellSl%nbnodes)
        do iBase = 1, FEBasisSl%size
            cplMap%nbDoFs = cplMap%nbDoFs+1
            cplMap%nbDoFsFECellSl = cplMap%nbDoFsFECellSl+1
            cplMap%mapDoFsFECellSl(cplMap%nbDoFsFECellSl) = cplMap%nbDoFs
        end do
!
        cplMap%nbDoFsFEFaceSl = 0
        cplMap%nbDoFsFace = 0
!
        do iNode = 1, nbNodes
            node = nodes(iNode)
            cplMap%nbDoFsFace = cplMap%nbDoFsFace+1
            cplMap%nbDoFsFEFaceSl = cplMap%nbDoFsFEFaceSl+1
            iDoF = cplMap%mapDoFsFECellSl(node)
            cplMap%mapDoFsFEFaceSl(cplMap%nbDoFsFEFaceSl) = iDoF
            cplMap%mapDoFsFace(cplMap%nbDoFsFace) = iDoF
        end do
!
        call hhoTherFaceDofs(hhoFaceMa, hhoDataMa, fbs)
        cplMap%nbDoFshhoFaceMa = 0
!
        do iBase = 1, fbs
            cplMap%nbDoFs = cplMap%nbDoFs+1
            cplMap%nbDoFsFace = cplMap%nbDoFsFace+1
            cplMap%nbDoFshhoFaceMa = cplMap%nbDoFshhoFaceMa+1
            cplMap%mapDoFshhoFaceMa(cplMap%nbDoFshhoFaceMa) = cplMap%nbDoFs
            cplMap%mapDoFsFace(cplMap%nbDoFsFace) = cplMap%nbDoFs
        end do
!
        ASSERT(cplMap%nbDoFs <= MSIZE_CPL_NITS)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplNitsInitData(option, FECellSl, cplMap, cplData)
!
        implicit none
!
        character(len=16), intent(in) :: option
        type(FE_Cell), intent(in) :: FECellSl
        type(CouplingMap), intent(in) :: cplMap
        type(CouplingData), intent(inout) :: cplData
!
!===================================================================================================
!    Initialize FEM and HHO data
!===================================================================================================
!
        integer(kind=8), parameter :: nbPara = 9
        integer(kind=8) :: jmater, iNode
!
        call cplPenaInitData(option, cplMap, cplData)
!
        call jevech('PMATERR', 'L', jmater)
!
        do iNode = 1, FECellSl%nbnodes
            cplData%E(iNode) = zr(jmater-1+nbPara*(iNode-1)+4)
            cplData%nu(iNode) = zr(jmater-1+nbPara*(iNode-1)+5)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplNitsLhsFEMHHO(FEFaceSl, FECellSl, hhoFaceMa, hhoDataMa, &
                                FEQuadSl, hhoQuadMa, cplMap, cplData, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(FE_Cell), intent(in) :: FECellSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuadMa
        type(CouplingMap), intent(in) :: cplMap
        type(CouplingData), intent(in) :: cplData
        type(HHO_matrix), intent(inout) :: lhs
!
!===================================================================================================
!    NITSCHE - FEM/HHO - Compute matrix
!===================================================================================================
!
        type(HHO_matrix) :: lhs_pena, lhs_stress, lhs_gradlift
!
        call cplCmpPenaFEMHHOMatVec(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                    cplData%coef_pena, lhs_pena)
!
        call cplCmpStressFEMHHOMat(FEFaceSl, FECellSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                   cplMap, cplData, lhs_stress)
!
        call cplCmpGradLiftFEMHHOMatVec(FEFaceSl, FECellSl, hhoFaceMa, hhoDataMa, &
                                        FEQuadSl, hhoQuadMa, cplData, lhs_gradlift)
!
        call cplNitsAssLhsFEMHHO(cplMap, lhs_pena, lhs_stress, lhs_gradlift, lhs)
!
        call lhs_pena%free()
        call lhs_stress%free()
        call lhs_gradlift%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplNitsAssLhsFEMHHO(cplMap, lhs_pena, lhs_stress, lhs_gradlift, lhs)
!
        implicit none
!
        type(CouplingMap), intent(in) :: cplMap
        type(HHO_matrix), intent(in) :: lhs_pena, lhs_stress, lhs_gradlift
        type(HHO_matrix), intent(inout) :: lhs
!
!===================================================================================================
!    NITSCHE - FEM/HHO - Assembly matrix
!===================================================================================================
!
        integer(kind=8) :: jFE, iFE, irow, jcol
!
        call lhs%initialize(cplMap%nbDoFs, cplMap%nbDoFs, 0.d0)
!
!      Terms acting on the projected interface jump:
!        gamma/h*(v^S-v^M, u^S-u^M)
!        +(C:sym(R_T(u)), sym(R_T(v)))_T
        ASSERT(lhs_pena%nrows == cplMap%nbDoFsFace)
        ASSERT(lhs_pena%ncols == cplMap%nbDoFsFace)
        ASSERT(lhs_gradlift%nrows == cplMap%nbDoFsFace)
        ASSERT(lhs_gradlift%ncols == cplMap%nbDoFsFace)
        do jFE = 1, cplMap%nbDoFsFace
            jcol = cplMap%mapDoFsFace(jFE)
            do iFE = 1, cplMap%nbDoFsFace
                irow = cplMap%mapDoFsFace(iFE)
                lhs%m(irow, jcol) = lhs_pena%m(iFE, jFE)+lhs_gradlift%m(iFE, jFE)
            end do
        end do
!
!      Term: - [(stress(u^S).n^S, v^S-v^M) + (stress(v^S).n^S, u^S-u^M)]
!      Rmq: This is a minus since in the paper they used n_I = n_hho = -n_cG
!          and by convention we used the outward normal for all intermediate computation
        ASSERT(lhs_stress%nrows == cplMap%nbDoFsFECellSl)
        ASSERT(lhs_stress%ncols == cplMap%nbDoFsFace)
        do jFE = 1, cplMap%nbDoFsFace
            jcol = cplMap%mapDoFsFace(jFE)
            do iFE = 1, cplMap%nbDoFsFECellSl
                irow = cplMap%mapDoFsFECellSl(iFE)
                lhs%m(irow, jcol) = lhs%m(irow, jcol)-lhs_stress%m(iFE, jFE)
                lhs%m(jcol, irow) = lhs%m(jcol, irow)-lhs_stress%m(iFE, jFE)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
end module
