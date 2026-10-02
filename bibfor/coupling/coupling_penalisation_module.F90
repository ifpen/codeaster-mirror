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
module coupling_penalisation_module
!
    use coupling_computation_module
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
#include "asterfort/readVector.h"
#include "blas/dsymv.h"
#include "FE_module.h"
!
! --------------------------------------------------------------------------------------------------
!
! Coupling - PENALISATION methods
!
! --------------------------------------------------------------------------------------------------
!
!
!
    public :: cplPenaInitTopoFEMHHO, cplPenaInitMapFEMHHOVec, cplPenaInitData
    public :: cplPenaInitTopoFEMFEM, cplPenaInitMapFEMFEMVec
    public :: cplPenaRhsFEMHHO, cplPenaLhsFEMHHO
    public :: cplPenaRhsFEMFEM, cplPenaLhsFEMFEM
    public :: cplPenaInitMapFEMHHOScal
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaInitTopoFEMHHO(FEFaceSl, hhoFaceMa, hhoDataMa)
!
        implicit none
!
        type(FE_Skin), intent(out) :: FEFaceSl
        type(HHO_Face), intent(out) :: hhoFaceMa
        type(HHO_Data), intent(out) :: hhoDataMa
!
!===================================================================================================
!    Initialize FEM and HHO topology - Penalisation
!===================================================================================================
!
        call FEFaceSl%init("SLAVE")
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
    subroutine cplPenaInitTopoFEMFEM(FEFaceSl, FEFaceMa)
!
        implicit none
!
        type(FE_Skin), intent(out) :: FEFaceSl, FEFaceMa
!
!===================================================================================================
!    Initialize FEM and FEM topology - Penalisation
!===================================================================================================
!
        call FEFaceSl%init("SLAVE")
        call FEFaceMa%init("MASTER")
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaInitMapFEMHHOVec(FEFaceSl, hhoFaceMa, hhoDataMa, cplMap)
!
        implicit none
!
        type(CouplingMap), intent(inout) :: cplMap
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
!
!===================================================================================================
!    Initialize FEM and HHO topology
!===================================================================================================
!
        type(FE_basis) :: FEBasisSl
        integer(kind=8) :: iBase, idim, fbs
!
        call FEBasisSl%initFace(FEFaceSl)
        ASSERT(FEBasisSl%typeEF == EF_LAGRANGE)
!
        cplMap%nbDoFs = 0
        cplMap%nbDoFsFEFaceSl = 0
!
        ASSERT(FEBasisSl%size == FEFaceSl%nbnodes)
        do iBase = 1, FEBasisSl%size
            do idim = 1, FEBasisSl%ndim
                cplMap%nbDoFs = cplMap%nbDoFs+1
                cplMap%nbDoFsFEFaceSl = cplMap%nbDoFsFEFaceSl+1
                cplMap%mapDoFsFEFaceSl(cplMap%nbDoFsFEFaceSl) = cplMap%nbDoFs
            end do
        end do
!
        call hhoMecaFaceDofs(hhoFaceMa, hhoDataMa, fbs)
        cplMap%nbDoFshhoFaceMa = 0
!
        do iBase = 1, fbs
            cplMap%nbDoFs = cplMap%nbDoFs+1
            cplMap%nbDoFshhoFaceMa = cplMap%nbDoFshhoFaceMa+1
            cplMap%mapDoFshhoFaceMa(cplMap%nbDoFshhoFaceMa) = cplMap%nbDoFs
        end do
!
        ASSERT(cplMap%nbDoFs <= MSIZE_CPL_PENA)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaInitMapFEMHHOScal(FEFaceSl, hhoFaceMa, hhoDataMa, cplMap)
!
        implicit none
!
        type(CouplingMap), intent(inout) :: cplMap
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
!
!===================================================================================================
!    Initialize FEM and HHO topology
!===================================================================================================
!
        type(FE_basis) :: FEBasisSl
        integer(kind=8) :: iBase, fbs
!
        call FEBasisSl%initFace(FEFaceSl)
        ASSERT(FEBasisSl%typeEF == EF_LAGRANGE)
!
        cplMap%nbDoFs = 0
        cplMap%nbDoFsFEFaceSl = 0
!
        ASSERT(FEBasisSl%size == FEFaceSl%nbnodes)
        do iBase = 1, FEBasisSl%size
            cplMap%nbDoFs = cplMap%nbDoFs+1
            cplMap%nbDoFsFEFaceSl = cplMap%nbDoFsFEFaceSl+1
            cplMap%mapDoFsFEFaceSl(cplMap%nbDoFsFEFaceSl) = cplMap%nbDoFs
        end do
!
        call hhoTherFaceDofs(hhoFaceMa, hhoDataMa, fbs)
        cplMap%nbDoFshhoFaceMa = 0
!
        do iBase = 1, fbs
            cplMap%nbDoFs = cplMap%nbDoFs+1
            cplMap%nbDoFshhoFaceMa = cplMap%nbDoFshhoFaceMa+1
            cplMap%mapDoFshhoFaceMa(cplMap%nbDoFshhoFaceMa) = cplMap%nbDoFs
        end do
!
        ASSERT(cplMap%nbDoFs <= MSIZE_CPL_PENA)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaInitMapFEMFEMVec(FEFaceSl, FEFaceMa, cplMap)
!
        implicit none
!
        type(CouplingMap), intent(inout) :: cplMap
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa
!
!===================================================================================================
!    Initialize FEM and HHO topology
!===================================================================================================
!
        type(FE_basis) :: FEBasisSl, FEBasisMa
        integer(kind=8) :: iBase, idim
!
        call FEBasisSl%initFace(FEFaceSl)
        call FEBasisMa%initFace(FEFaceMa)
!
        cplMap%nbDoFs = 0
        cplMap%nbDoFsFEFaceSl = 0
!
        ASSERT(FEBasisSl%size == FEFaceSl%nbnodes)
        do iBase = 1, FEBasisSl%size
            do idim = 1, FEBasisSl%ndim
                cplMap%nbDoFs = cplMap%nbDoFs+1
                cplMap%nbDoFsFEFaceSl = cplMap%nbDoFsFEFaceSl+1
                cplMap%mapDoFsFEFaceSl(cplMap%nbDoFsFEFaceSl) = cplMap%nbDoFs
            end do
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
            end do
        end do
!
        ASSERT(cplMap%nbDoFs <= MSIZE_CPL_PENA_FEM)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaInitData(option, cplMap, cplData)
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
        integer(kind=8) :: jpair
!
        cplData%nbDoFs = cplMap%nbDoFs
        if (option == "CHAR_MECA_CPL") then
            call readVector('PDEPLMR', cplMap%nbDoFs, cplData%disp_prev)
            call readVector('PDEPLPR', cplMap%nbDoFs, cplData%disp_curr)
            call daxpy_1(cplMap%nbDoFs, 1.d0, cplData%disp_prev, cplData%disp_curr)
        end if
!
        call jevech('PPAIRR', 'L', jpair)
        cplData%coef_pena = zr(jpair-1+OFFSET_COEF_PENA)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaRhsFEMHHO(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuad, cplData, rhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuad
        type(CouplingData), intent(in) :: cplData
        real(kind=8), intent(inout) :: rhs(MSIZE_CPL_PENA)
!
!===================================================================================================
!    PENALISATION - FEM/HHO - Compute residual
!===================================================================================================
!
        type(HHO_matrix) :: mat
!
        rhs = 0.d0
!
        call cplCmpPenaFEMHHOMatVec(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuad, &
                                    cplData%coef_pena, mat)
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
    subroutine cplPenaLhsFEMHHO(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuad, cplData, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuad
        type(CouplingData), intent(in) :: cplData
        type(HHO_matrix), intent(inout) :: lhs
!
!===================================================================================================
!    PENALISATION - FEM/HHO - Compute matrix
!===================================================================================================
!
        call cplCmpPenaFEMHHOMatVec(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuad, &
                                    cplData%coef_pena, lhs)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaRhsFEMFEM(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, cplData, rhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        type(CouplingData), intent(in) :: cplData
        real(kind=8), intent(inout) :: rhs(MSIZE_CPL_PENA_FEM)
!
!===================================================================================================
!    PENALISATION - FEM/FEM - Compute residual
!===================================================================================================
!
        real(kind=8) :: mat(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
        blas_int :: b_n, b_lda
        blas_int, parameter :: b_one = to_blas_int(1)
!
        rhs = 0.d0
!
        call cplCmpPenaFEMFEMMatVec(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, &
                                    cplData%coef_pena, mat)
!
        b_n = to_blas_int(cplData%nbDoFs)
        b_lda = to_blas_int(MSIZE_CPL_PENA_FEM)
        call dsymv("U", b_n, 1.d0, mat, b_lda, cplData%disp_curr, b_one, 0.d0, rhs, b_one)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaLhsFEMFEM(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, cplData, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        type(CouplingData), intent(in) :: cplData
        real(kind=8), intent(inout) :: lhs(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
!
!===================================================================================================
!    PENALISATION - FEM/FEM - Compute matrix
!===================================================================================================
!
        call cplCmpPenaFEMFEMMatVec(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, &
                                    cplData%coef_pena, lhs)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
end module
