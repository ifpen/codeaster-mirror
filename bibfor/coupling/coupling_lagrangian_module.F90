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
module coupling_lagrangian_module
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
#include "asterfort/coupling_type.h"
#include "asterfort/HHO_size_module.h"
#include "FE_module.h"
!
! --------------------------------------------------------------------------------------------------
!
! Coupling - Augmented Lagrangian methods
!
! --------------------------------------------------------------------------------------------------
!
!
!
    public :: cplLagrInitData
    public :: cplLagrInitTopoFEMFEM, cplLagrInitMapFEMFEMVec
    public :: cplLagrRhsFEMFEM, cplLagrLhsFEMFEM
    private :: cplLagrAssLhsFEMFEM
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplLagrInitTopoFEMFEM(FEFaceSl, FEFaceMa, FEFaceLagSl)
!
        implicit none
!
        type(FE_Skin), intent(out) :: FEFaceSl, FEFaceMa, FEFaceLagSl
!
!===================================================================================================
!    Initialize FEM and FEM topology - Lagrangian
!===================================================================================================
!
        call cplPenaInitTopoFEMFEM(FEFaceSl, FEFaceMa)
        ! Same degree used for DEPL and LAGR on Slave side
        FEFaceLagSl = FEFaceSl
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplLagrInitMapFEMFEMVec(FEFaceSl, FEFaceMa, FEFaceLagSl, cplMap)
!
        implicit none
!
        type(CouplingMap), intent(inout) :: cplMap
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa, FEFaceLagSl
!
!===================================================================================================
!    Initialize FEM and HHO topology
!===================================================================================================
!
        type(FE_basis) :: FEBasisSl, FEBasisMa, FEBasisLagSl
        integer(kind=8) :: iBase, idim
!
        call FEBasisSl%initFace(FEFaceSl)
        call FEBasisMa%initFace(FEFaceMa)
        call FEBasisLagSl%initFace(FEFaceLagSl)
!
        cplMap%nbDoFs = 0
        cplMap%nbDoFsFEFaceSl = 0
        cplMap%nbDoFsFEFaceLagSl = 0
        cplMap%nbDoFsFace = 0
!
        ASSERT(FEBasisSl%size == FEFaceSl%nbnodes)
        ASSERT(FEBasisSl%size >= FEBasisLagSl%size)
        do iBase = 1, FEBasisSl%size
            do idim = 1, FEBasisSl%ndim
                cplMap%nbDoFs = cplMap%nbDoFs+1
                cplMap%nbDoFsFEFaceSl = cplMap%nbDoFsFEFaceSl+1
                cplMap%mapDoFsFEFaceSl(cplMap%nbDoFsFEFaceSl) = cplMap%nbDoFs
                cplMap%nbDoFsFace = cplMap%nbDoFsFace+1
                cplMap%mapDoFsFace(cplMap%nbDoFsFace) = cplMap%nbDoFs
            end do
            if (iBase <= FEBasisLagSl%size) then
                do idim = 1, FEBasisLagSl%ndim
                    cplMap%nbDoFs = cplMap%nbDoFs+1
                    cplMap%nbDoFsFEFaceLagSl = cplMap%nbDoFsFEFaceLagSl+1
                    cplMap%mapDoFsFEFaceLagSl(cplMap%nbDoFsFEFaceLagSl) = cplMap%nbDoFs
                end do
            end if
        end do
!
        cplMap%nbDoFsFEFaceMa = 0
!
        ASSERT(FEBasisMa%size == FEFaceMa%nbnodes)
        do iBase = 1, FEBasisma%size
            do idim = 1, FEBasisMa%ndim
                cplMap%nbDoFs = cplMap%nbDoFs+1
                cplMap%nbDoFsFEFaceMa = cplMap%nbDoFsFEFaceMa+1
                cplMap%mapDoFsFEFaceMa(cplMap%nbDoFsFEFaceMa) = cplMap%nbDoFs
                cplMap%nbDoFsFace = cplMap%nbDoFsFace+1
                cplMap%mapDoFsFace(cplMap%nbDoFsFace) = cplMap%nbDoFs
            end do
        end do
!
        ASSERT(cplMap%nbDoFs <= MSIZE_CPL_LAGR_FEM)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplLagrInitData(option, cplMap, cplData)
!
        implicit none
!
        character(len=16), intent(in) :: option
        type(CouplingMap), intent(in) :: cplMap
        type(CouplingData), intent(inout) :: cplData
!
!===================================================================================================
!    Initialize FEM and HHO data
!===================================================================================================
!
        call cplPenaInitData(option, cplMap, cplData)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplLagrRhsFEMFEM(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, &
                                cplMap, cplData, rhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa, FEFaceLagSl
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        type(CouplingMap), intent(in) :: cplMap
        type(CouplingData), intent(in) :: cplData
        real(kind=8), intent(inout) :: rhs(MSIZE_CPL_LAGR_FEM)
!
!===================================================================================================
!    Lagrangian - FEM/FEM - Compute residual
!===================================================================================================
!
        type(HHO_matrix) :: mat
!
        rhs = 0.d0
!
        call cplLagrLhsFEMFEM(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, &
                              cplMap, cplData, mat)
!
        call mat%dot(cplData%disp_curr, rhs)
        call mat%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplLagrLhsFEMFEM(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, &
                                cplMap, cplData, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa, FEFaceLagSl
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        type(CouplingMap), intent(in) :: cplMap
        type(CouplingData), intent(in) :: cplData
        type(HHO_matrix), intent(inout) :: lhs
!
!===================================================================================================
!    Lagrangian - FEM/FEM - Compute matrix
!===================================================================================================
!
        real(kind=8) :: matPena(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
        real(kind=8) :: matLagr(MSIZE_LAGR_FEM, MSIZE_CPL_PENA_FEM)
!
        call cplCmpPenaFEMFEMMatVec(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, &
                                    cplData%coef_pena, matPena)
!
        call cplCmpLagrFEMFEMMatVec(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, matLagr)
!
        call cplLagrAssLhsFEMFEM(cplMap, matPena, matLagr, lhs)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplLagrAssLhsFEMFEM(cplMap, matPena, matLagr, lhs)
!
        implicit none
!
        type(CouplingMap), intent(in) :: cplMap
        real(kind=8), intent(in) :: matPena(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
        real(kind=8), intent(in) :: matLagr(MSIZE_LAGR_FEM, MSIZE_CPL_PENA_FEM)
        type(HHO_matrix), intent(inout) :: lhs
!
!===================================================================================================
!    Lagrangian - FEM/FEM - Assembly matrix
!===================================================================================================
!
        integer(kind=8) :: jFESlLg, jFE, iFE, irow, jcol
!
        call lhs%initialize(cplMap%nbDoFs, cplMap%nbDoFs, 0.d0)
!
!      Term: gamma/h*(v^S - v^M, u^S-u^M)
        do jFE = 1, cplMap%nbDoFsFace
            jcol = cplMap%mapDoFsFace(jFE)
            do iFE = 1, cplMap%nbDoFsFace
                irow = cplMap%mapDoFsFace(iFE)
                lhs%m(irow, jcol) = matPena(iFE, jFE)
            end do
        end do
!
!     Term: (v^S-v^M, lag^S) and (mu^S, u^S-u^M) by symmetry
        do jFESlLg = 1, cplMap%nbDoFsFEFaceLagSl
            jcol = cplMap%mapDoFsFEFaceLagSl(jFESlLg)
            do iFE = 1, cplMap%nbDoFsFace
                irow = cplMap%mapDoFsFace(iFE)
                lhs%m(jcol, irow) = matLagr(jFESlLg, iFE)
                lhs%m(irow, jcol) = lhs%m(jcol, irow)
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
