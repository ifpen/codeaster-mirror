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
module HHO_massmat_module
!
    use HHO_basis_module
    use HHO_quadrature_module
    use HHO_type
    use HHO_utils_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterf_debug.h"
#include "asterfort/assert.h"
#include "asterfort/HHO_size_module.h"
#include "blas/dsyr.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Module to compute mass matrix
!
! --------------------------------------------------------------------------------------------------
    type HHO_massmat_cell
        integer(kind=8) :: max_nrows = MSIZE_CELL_SCAL, max_ncols = MSIZE_CELL_SCAL
        integer(kind=8) :: nrows = 0, ncols = 0
        aster_logical :: isIdentity = ASTER_FALSE
        real(kind=8) :: m(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL)
!
! ----- member function
    contains
        procedure, pass :: compute => hhoMassMatCellScal
!
    end type
!
    type HHO_massmat_face
        integer(kind=8) :: max_nrows = MSIZE_FACE_SCAL, max_ncols = MSIZE_FACE_SCAL
        integer(kind=8) :: nrows = 0, ncols = 0
        aster_logical :: isIdentity = ASTER_FALSE
        real(kind=8) :: m(MSIZE_FACE_SCAL, MSIZE_FACE_SCAL)
!
! ----- member function
    contains
        procedure, pass :: compute => hhoMassMatFaceScal
!
    end type
!
    public :: HHO_massmat_cell, HHO_massmat_face
    private :: hhoMassMatCellScal, hhoMassMatFaceScal
!    private  ::
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMassMatCellScal(this, hhoCell, min_order, max_order, hhoQuad_, type_basis_)
!
        implicit none
!
        class(HHO_massmat_cell), intent(inout) :: this
        type(HHO_Cell), intent(in) :: hhoCell
        integer(kind=8), intent(in) :: min_order
        integer(kind=8), intent(in) :: max_order
        type(HHO_quadrature), intent(in), optional :: hhoQuad_
        integer(kind=8), optional, intent(in) :: type_basis_

!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the scalar mass matrix of a Cell form order "min_order" to order "max_order"
!   In hhoCell     : the current HHO Cell
!   In min_order   : minimum order to evaluate
!   In max_order   : maximum order to evaluate
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_quadrature) :: hhoQuad
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: basisScalEval
        integer(kind=8) :: dimMat, ipg, i
        aster_logical :: dbg
        blas_int :: b_incx, b_lda, b_n
! --------------------------------------------------------------------------------------------------
!
! ----- init basis
        if (present(type_basis_)) then
            call hhoBasisCell%initialize(hhoCell, type_basis_)
        else
            call hhoBasisCell%initialize(hhoCell)
        end if
! ----- dimension of massMat
        dimMat = hhoBasisCell%BSSize(min_order, max_order)
        this%m = 0.d0
        this%nrows = dimMat
        this%ncols = dimMat
        ASSERT(this%nrows <= this%max_nrows)
        ASSERT(this%ncols <= this%max_ncols)
!
#ifdef ASTER_DEBUG_FC
        dbg = ASTER_TRUE
#else
        dbg = ASTER_FALSE
#endif
!
        if (hhoBasisCell%isOrthonormal() .and. .not. dbg) then
            do i = 1, dimMat
                this%m(i, i) = 1.d0
            end do
            this%isIdentity = ASTER_TRUE
        else
!
! ----- get quadrature
            if (present(hhoQuad_)) then
                hhoQuad = hhoQuad_
                ASSERT(2*max_order <= hhoQuad%order)
            else
                call hhoQuad%GetQuadCell(hhoCell, 2*max_order)
            end if
!
! ----- Loop on quadrature point
            do ipg = 1, hhoQuad%nbQuadPoints
! --------- Eval bais function at the quadrature point
                call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), min_order, max_order, &
                                         basisScalEval)
! --------  Eval massMat
                b_n = to_blas_int(dimMat)
                b_incx = to_blas_int(1)
                b_lda = to_blas_int(this%max_nrows)
                call dsyr('U', b_n, hhoQuad%weights(ipg), basisScalEval, b_incx, &
                          this%m, b_lda)
            end do
!
! ----- Copy the lower part
!
            call hhoCopySymPartMat('U', this%m(1:dimMat, 1:dimMat))
!
            if (hhoBasisCell%isOrthonormal() .and. dbg) then
                if (hhoCell%measure .ge. 1.d-8) then
                    if (.not. hhoIsIdentityMat(this%m, dimMat)) then
                        call hhoPrintMat(this%m(1:dimMat, 1:dimMat))
                        ASSERT(ASTER_FALSE)
                    end if
                end if
            end if
!
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMassMatFaceScal(this, hhoFace, min_order, max_order, hhoQuad_)
!
        implicit none
!
        class(HHO_massmat_face), intent(inout) :: this
        type(HHO_Face), intent(in) :: hhoFace
        integer(kind=8), intent(in) :: min_order
        integer(kind=8), intent(in) :: max_order
        type(HHO_quadrature), intent(in), optional :: hhoQuad_
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the scalar mass matrix of a Face form order "min_order" to order "max_order"
!   In hhoFace     : the current HHO Face
!   In min_order   : minimum order to evaluate
!   In max_order   : maximum order to evaluate
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_face) :: hhoBasisFace
        type(HHO_quadrature) :: hhoQuad
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: basisScalEval
        integer(kind=8) :: dimMat, ipg, i
        aster_logical :: dbg
        blas_int :: b_incx, b_lda, b_n
! --------------------------------------------------------------------------------------------------
!
! ----- init basis
        call hhoBasisFace%initialize(hhoFace)
! ----- dimension of massMat
        dimMat = hhoBasisFace%BSSize(min_order, max_order)
!
        this%m = 0.d0
        this%nrows = dimMat
        this%ncols = dimMat
        ASSERT(this%nrows <= this%max_nrows)
        ASSERT(this%ncols <= this%max_ncols)
!
#ifdef ASTER_DEBUG_FC
        dbg = ASTER_TRUE
#else
        dbg = ASTER_FALSE
#endif
!
        if (hhoBasisFace%isOrthonormal() .and. .not. dbg) then
            do i = 1, dimMat
                this%m(i, i) = 1.d0
            end do
            this%isIdentity = ASTER_TRUE
        else
!
! ----- get quadrature
            if (present(hhoQuad_)) then
                hhoQuad = hhoQuad_
                ASSERT(2*max_order <= hhoQuad%order)
            else
                call hhoQuad%GetQuadFace(hhoFace, 2*max_order)
            end if
!
! ----- Loop on quadrature point
            do ipg = 1, hhoQuad%nbQuadPoints
! --------- Eval bais function at the quadrature point
                call hhoBasisFace%BSEval(hhoQuad%points(1:3, ipg), min_order, max_order, &
                                         basisScalEval)
! --------  Eval massMat
                b_n = to_blas_int(dimMat)
                b_incx = to_blas_int(1)
                b_lda = to_blas_int(this%max_nrows)
                call dsyr('U', b_n, hhoQuad%weights(ipg), basisScalEval, b_incx, &
                          this%m, b_lda)
            end do
!
! ----- Copy the lower part
!
            call hhoCopySymPartMat('U', this%m(1:dimMat, 1:dimMat))
!
            if (hhoBasisFace%isOrthonormal() .and. dbg) then
                if (hhoFace%measure .ge. 1.d-7) then
                    ASSERT(hhoIsIdentityMat(this%m, dimMat))
                end if
            end if
        end if
!
! call hhoPrintMat(this%m)
!
    end subroutine
!
end module
