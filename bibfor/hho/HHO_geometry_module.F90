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
module HHO_geometry_module
!
    use HHO_type
    use HHO_utils_module, only: CellNameL2S
    use compensated_ops_module, only: dot_product, norm2, addFMA, det

!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/apnorm.h"
#include "asterfort/assert.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfvf.h"
#include "asterfort/provec.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - Geometry module
!
! Compute outward normal of a surface
! Compute barycenter of an element
!
! --------------------------------------------------------------------------------------------------
!
    public :: barycenter, hhoFaceInitCoor, hhoGeomBasis, hhoGeomDerivBasis
    public :: hhoLocalTangentBasisFace, hhoNormalFace2, hhoNormalFace3, hhoNormalFaceQP
    public :: hhoSplitSimplex, hho_transfo_3d, hho_transfo_quad, hhoIsJacobCst, prod_vec
    public :: hhoTangentFaceQP, hhoLocalBasisFaceQP
    private :: hho_jaco_cst_quad, hho_jaco_cst_3d
    private :: hhoNormalFace2d, well_oriented, hhoNormalFace1d, find_lowest_vertex
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    function find_lowest_vertex(v1, v2) result(ind)
!
        implicit none
!
        real(kind=8), dimension(3), intent(in) :: v1, v2
        integer(kind=8) :: ind
!
! --------------------------------------------------------------------------------------------------
! Find the nodes at bottom left between two nodes
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), parameter :: prec = 1.d-10
        real(kind=8) :: norm, v(3)
! --------------------------------------------------------------------------------------------------
!
!
        norm = max(max(norm2(v1), norm2(v2)), 1.d0)
        v = (v1-v2)/norm
!
        if (abs(v(1)) > prec) then
            if (v(1) > 0.d0) then
                ind = 1
            else
                ind = 2
            end if
        else if (abs(v(2)) > prec) then
            if (v(2) > 0.d0) then
                ind = 1
            else
                ind = 2
            end if
        else if (abs(v(3)) > prec) then
            if (v(3) > 0.d0) then
                ind = 1
            else
                ind = 2
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function barycenter(nodes, nbnodes) result(bar)
!
        implicit none
!
        integer(kind=8), intent(in) :: nbnodes
        real(kind=8), dimension(3, nbnodes), intent(in) :: nodes
        real(kind=8), dimension(3) :: bar
!
! --------------------------------------------------------------------------------------------------
!  In nodes        :: list of nodes
!  In nbnodes      :: number of nodes
!  In ndim         :: topological dimension
!  Out bar         :: barycenter of nodes
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: inode
! --------------------------------------------------------------------------------------------------
!
        bar = 0.d0
!
        do inode = 1, nbnodes
            bar(1:3) = bar(1:3)+nodes(1:3, inode)
        end do
!
        bar(1:3) = bar(1:3)/real(nbnodes, kind=8)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function prod_vec(v0, v1) result(v2)
!
        implicit none
!
        real(kind=8), dimension(3), intent(in) :: v0
        real(kind=8), dimension(3), intent(in) :: v1
        real(kind=8), dimension(3) :: v2
!
! --------------------------------------------------------------------------------------------------
!  In v0        :: vector 0
!  In v1        :: vector 1
!  Out v2       :: vector 2 = v0 x v1
! --------------------------------------------------------------------------------------------------
!
        v2 = 0.d0
!
        v2(1) = v0(2)*v1(3)-v0(3)*v1(2)
        v2(2) = -(v0(1)*v1(3)-v0(3)*v1(1))
        v2(3) = v0(1)*v1(2)-v0(2)*v1(1)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function well_oriented(v0, normal)
!
        implicit none
!
        real(kind=8), dimension(3), intent(in) :: v0
        real(kind=8), dimension(3), intent(in) :: normal
        aster_logical :: well_oriented
!
! --------------------------------------------------------------------------------------------------
!  In v0        :: vector 0
!  In normal    :: normal
!  Out logical  :: is well oriented
! --------------------------------------------------------------------------------------------------
!
        if (dot_product(v0, normal) .le. r8prem()) then
            well_oriented = ASTER_FALSE
        else
            well_oriented = ASTER_TRUE
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoNormalFaceQP(hhoFace, qp_param) result(normal)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        real(kind=8), dimension(2), intent(in) :: qp_param
        real(kind=8), dimension(3) :: normal
!
! --------------------------------------------------------------------------------------------------
!  In coorno             :: coordinates of the nodes
!  In nbnodes            :: number of nodes
!  In barycenter_face    :: barycenter of the face
!  In barycenter_cell    :: barycenter of the cell
!  Out normal            :: outward normal
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: tau1(3), tau2(3)
!  -------------------------------------------------------------------------------------------------

        call hhoLocalBasisFaceQP(hhoFace, qp_param, tau1, tau2, normal)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoNormalFace2d(coorno, nbnodes, barycenter_face, barycenter_cell) result(normal)
!
        implicit none
!
        real(kind=8), dimension(3, 4), intent(in) :: coorno
        integer(kind=8), intent(in) :: nbnodes
        real(kind=8), dimension(3), optional, intent(in) :: barycenter_face
        real(kind=8), dimension(3), optional, intent(in) :: barycenter_cell
        real(kind=8), dimension(3) :: normal
!
! --------------------------------------------------------------------------------------------------
!  In coorno             :: coordinates of the nodes
!  In nbnodes            :: number of nodes
!  In barycenter_face    :: barycenter of the face
!  In barycenter_cell    :: barycenter of the cell
!  Out normal            :: outward normal
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) :: v0, v1, vbar
!  -------------------------------------------------------------------------------------------------
        normal = 0.d0
!
        v0(1:3) = coorno(1:3, 2)-coorno(1:3, 1)
        v1(1:3) = coorno(1:3, nbnodes)-coorno(1:3, 1)
!
        normal = prod_vec(v0, v1)
        normal = normal/norm2(normal)
!
! ---- Test normal
        if (present(barycenter_cell)) then
            vbar(1:3) = barycenter_face(1:3)-barycenter_cell(1:3)
            if (.not. well_oriented(vbar, normal)) then
                normal = -normal
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoNormalFace1d(coorno, barycenter_face, barycenter_cell) result(normal)
!
        implicit none
!
        real(kind=8), dimension(3, 4), intent(in) :: coorno
        real(kind=8), dimension(3), optional, intent(in) :: barycenter_face
        real(kind=8), dimension(3), optional, intent(in) :: barycenter_cell
        real(kind=8), dimension(3) :: normal
!
! --------------------------------------------------------------------------------------------------
!  In coorno             :: coordinates of the node
!  In barycenter_face    :: barycenter of the face
!  In barycenter_cell    :: barycenter of the cell
!  Out normal            :: outward normal
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) :: tangente, vbar
! --------------------------------------------------------------------------------------------------
! -- Normal to the face
        tangente = coorno(1:3, 2)-coorno(1:3, 1)
!
        normal = (/-tangente(2), tangente(1), 0.d0/)
        normal = -normal/norm2(normal)
!
! ---- Test normal
        if (present(barycenter_cell)) then
            vbar(1:3) = barycenter_face(1:3)-barycenter_cell(1:3)
            if (.not. well_oriented(vbar, normal)) then
                normal = -normal
            end if
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoNormalFace3(hhoFace, barycenter_cell) result(normal)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        real(kind=8), dimension(3), intent(in) :: barycenter_cell
        real(kind=8), dimension(3) :: normal
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Face           :: face HHO
!  In barycenter_cell    :: barycenter of the cell
!  In normal cell        :: normal to the cell (use for 2D cell)
!  Out normal            :: outward normal of the face
! --------------------------------------------------------------------------------------------------
!
! --------------------------------------------------------------------------------------------------
        normal = 0.d0
!
        if (hhoFace%typema == MT_QUAD4) then
            normal = hhoNormalFace2d(hhoFace%coorno, 4, hhoFace%barycenter, barycenter_cell)
        else if (hhoFace%typema == MT_TRIA3) then
            normal = hhoNormalFace2d(hhoFace%coorno, 3, hhoFace%barycenter, barycenter_cell)
        else if (hhoFace%typema == MT_SEG2) then
            normal = hhoNormalFace1d(hhoFace%coorno, hhoFace%barycenter, barycenter_cell)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoNormalFace2(typma, nodes_coor) result(normal)
!
        implicit none
!
        integer(kind=8), intent(in) :: typma
        real(kind=8), dimension(3, 4), intent(in) :: nodes_coor
        real(kind=8), dimension(3) :: normal
!
! --------------------------------------------------------------------------------------------------
!  In typma              :: type of face
!  In nodes_coor         :: coordinates of the face
!  In normal cell        :: normal to the cell (use for 2D cell)
!  Out normal            :: outward normal of the face
! --------------------------------------------------------------------------------------------------
!
! --------------------------------------------------------------------------------------------------
        normal = 0.d0
!
        if (typma == MT_QUAD4) then
            normal = hhoNormalFace2d(nodes_coor, 4)
        else if (typma == MT_TRIA3) then
            normal = hhoNormalFace2d(nodes_coor, 3)
        else if (typma == MT_SEG2) then
            normal = hhoNormalFace1d(nodes_coor)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoLocalBasisFaceQP(hhoFace, qp_param, tau1, tau2, normal, jaco)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        real(kind=8), dimension(2), intent(in) :: qp_param
        real(kind=8), dimension(3), intent(out) :: tau1, tau2, normal
        real(kind=8), optional, intent(out) :: jaco
!
! --------------------------------------------------------------------------------------------------
!  In coorno             :: coordinates of the nodes
!  In nbnodes            :: number of nodes
!  In barycenter_face    :: barycenter of the face
!  In barycenter_cell    :: barycenter of the cell
!  Out normal            :: outward normal
! --------------------------------------------------------------------------------------------------
!
        character(len=8) :: ts
        real(kind=8) :: coor(3, MT_NNOMAX2D), tau(3), noor
!  -------------------------------------------------------------------------------------------------
        normal = 0.d0
        call CellNameL2S(hhoFace%typema, ts)
        coor(1:3, 1:4) = hhoFace%coorno(1:3, 1:4)
!
        call apnorm(hhoFace%nbnodes, ts, hhoFace%ndim+1, coor, qp_param(1), &
                    qp_param(2), normal, tau1, tau2, noor)
!
        if (norm2(hhoFace%normal) > 0.5d0) then
            if (.not. well_oriented(hhoFace%normal, normal)) then
                normal = -normal
                if (hhoFace%ndim == 1) then
                    tau1 = -tau1
                    tau2 = 0.d0
                else
                    tau = tau1
                    tau1 = tau2
                    tau2 = tau
                end if
            end if
        end if
!
        if (present(jaco)) then
            jaco = noor
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoTangentFaceQP(hhoFace, qp_param, tau1, tau2)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        real(kind=8), dimension(2), intent(in) :: qp_param
        real(kind=8), dimension(3), intent(out) :: tau1, tau2
!
! --------------------------------------------------------------------------------------------------
!  In coorno             :: coordinates of the nodes
!  In nbnodes            :: number of nodes
!  In barycenter_face    :: barycenter of the face
!  In barycenter_cell    :: barycenter of the cell
!  Out normal            :: outward normal
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) ::  normal(3)
!  -------------------------------------------------------------------------------------------------
!
        call hhoLocalBasisFaceQP(hhoFace, qp_param, tau1, tau2, normal)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoLocalTangentBasisFace(hhoFace) result(axes)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        real(kind=8), dimension(3, 2) :: axes
!
! --------------------------------------------------------------------------------------------------
!   HHO - geometry
!
!   compute orthonormal local basis of the face
!   In hhoFace              : the current HHO Face
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) :: v0, v1
!
        axes = 0.d0
!
        if (hhoFace%ndim == 2) then
            v0 = hhoFace%coorno(1:3, 2)-hhoFace%coorno(1:3, 1)
            v1 = hhoFace%coorno(1:3, hhoFace%nbnodes)-hhoFace%coorno(1:3, 1)
!
            v1 = v1-(dot_product(v0, v1)*v0)/(dot_product(v0, v0))
!
            axes(1:3, 1) = v0/norm2(v0)
            axes(1:3, 2) = v1/norm2(v1)
!
        else if (hhoFace%ndim == 1) then
            v0 = hhoFace%coorno(1:3, 2)-hhoFace%coorno(1:3, 1)
!
            axes(1:3, 1) = v0/norm2(v0)
!
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoFaceInitCoor(coorno, nnos, nbnodes, ndimF, numsorted_) result(nodes_face)
!
        implicit none
!
        integer(kind=8), intent(in) :: ndimF
        real(kind=8), dimension(3, 9), intent(in) :: coorno
        integer(kind=8), intent(in) :: nbnodes, nnos
        real(kind=8), dimension(3, 9) :: nodes_face
        integer(kind=8), intent(out), optional :: numsorted_(9)
!
! --------------------------------------------------------------------------------------------------
!   We have to reorder the nodes of the face to use the same basis functions for a face
!     which is shared by two cells
!  In HHO_Face           :: face HHO
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ino, minnum1, minnum2, numsorted(9), ind, candidate(2)
!
        numsorted = 0
!
        if (ndimF == 1) then
            ASSERT(nnos == 2)
            ind = find_lowest_vertex(coorno(:, 1), coorno(:, 2))
!
            if (ind == 1) then
                numsorted(1:3) = (/1, 2, 3/)
            else if (ind == 2) then
                numsorted(1:3) = (/2, 1, 3/)
            else
                ASSERT(ASTER_FALSE)
            end if
        else if (ndimF == 2) then
            minnum1 = 1
            candidate = (/2, nnos/)
            do ino = 2, nnos
                ind = find_lowest_vertex(coorno(:, minnum1), coorno(:, ino))
                if (ind == 2) then
                    minnum1 = ino
                    if (ino == nnos) then
                        candidate = (/ino-1, 1/)
                    else
                        candidate = (/ino-1, ino+1/)
                    end if
                end if
            end do
!
            minnum2 = candidate( &
                      find_lowest_vertex(coorno(:, candidate(1)), coorno(:, candidate(2))))
!
            if (nnos == 3) then
                if (minnum1 == 1) then
                    if (minnum2 == 2) then
                        numsorted(1:7) = (/1, 2, 3, 4, 5, 6, 7/)
                    else if (minnum2 == 3) then
                        numsorted(1:7) = (/1, 3, 2, 6, 5, 4, 7/)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                else if (minnum1 == 2) then
                    if (minnum2 == 1) then
                        numsorted(1:7) = (/2, 1, 3, 4, 6, 5, 7/)
                    else if (minnum2 == 3) then
                        numsorted(1:7) = (/2, 3, 1, 5, 6, 4, 7/)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                else if (minnum1 == 3) then
                    if (minnum2 == 1) then
                        numsorted(1:7) = (/3, 1, 2, 6, 4, 5, 7/)
                    else if (minnum2 == 2) then
                        numsorted(1:7) = (/3, 2, 1, 5, 4, 6, 7/)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                else
                    ASSERT(ASTER_FALSE)
                end if
            else if (nnos == 4) then
                if (minnum1 == 1) then
                    if (minnum2 == 2) then
                        numsorted(1:9) = (/1, 2, 3, 4, 5, 6, 7, 8, 9/)
                    else if (minnum2 == 4) then
                        numsorted(1:9) = (/1, 4, 3, 2, 8, 7, 6, 5, 9/)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                else if (minnum1 == 2) then
                    if (minnum2 == 1) then
                        numsorted(1:9) = (/2, 1, 4, 3, 5, 8, 7, 6, 9/)
                    else if (minnum2 == 3) then
                        numsorted(1:9) = (/2, 3, 4, 1, 6, 7, 8, 5, 9/)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                else if (minnum1 == 3) then
                    if (minnum2 == 2) then
                        numsorted(1:9) = (/3, 2, 1, 4, 6, 5, 8, 7, 9/)
                    else if (minnum2 == 4) then
                        numsorted(1:9) = (/3, 4, 1, 2, 7, 8, 5, 6, 9/)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                else if (minnum1 == 4) then
                    if (minnum2 == 1) then
                        numsorted(1:9) = (/4, 1, 2, 3, 8, 5, 6, 7, 9/)
                    else if (minnum2 == 3) then
                        numsorted(1:9) = (/4, 3, 2, 1, 7, 6, 5, 8, 9/)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                else
                    ASSERT(ASTER_FALSE)
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
!
        else
            ASSERT(ASTER_FALSE)
        end if
!
! --- Copy the coordinates
        nodes_face = 0.d0
        do ino = 1, nbnodes
            nodes_face(1:3, ino) = coorno(1:3, numsorted(ino))
        end do
!
        if (present(numsorted_)) then
            numsorted_ = numsorted
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGeomBasis(typema, pt, basis)
!
        implicit none
!
        integer(kind=8), intent(in) :: typema
        real(kind=8), intent(in) :: pt(3)
        real(kind=8), intent(out) :: basis(8)
!
! ---------------------------------------------------------------------------------
!  HHO - geometrie
!  Compute the functions which descript the geometry
!
! In typema : type of element
! In pt     : coordinate of the point
! In basis  : evaluation of the basis at the point pt
! ---------------------------------------------------------------------------------
!
        basis = 0.d0
!
        select case (typema)
        case (MT_SEG2)
            call elrfvf('SE2', pt, basis)
        case (MT_TRIA3)
            call elrfvf('TR3', pt, basis)
        case (MT_QUAD4)
            call elrfvf('QU4', pt, basis)
        case (MT_TETRA4)
            call elrfvf('TE4', pt, basis)
        case (MT_PYRAM5)
            call elrfvf('PY5', pt, basis)
        case (MT_HEXA8)
            call elrfvf('HE8', pt, basis)
        case (MT_PENTA6)
            call elrfvf('PE6', pt, basis)
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGeomDerivBasis(typema, pt, dbasis)
!
        implicit none
!
        integer(kind=8), intent(in) :: typema
        real(kind=8), intent(in) :: pt(3)
        real(kind=8), intent(out) :: dbasis(3, 8)
!
! ---------------------------------------------------------------------------------
!  HHO - geometrie
!  Compute the derivative of functions which descript the geometry
!
! In typema : type of element
! In pt     : coordinate of the point
! In dbasis  : evaluation of the derivative of basis at the point pt
! ---------------------------------------------------------------------------------
!
!
        dbasis = 0.d0
!
        select case (typema)
        case (MT_SEG2)
            call elrfdf('SE2', pt, dbasis)
        case (MT_TRIA3)
            call elrfdf('TR3', pt, dbasis)
        case (MT_QUAD4)
            call elrfdf('QU4', pt, dbasis)
        case (MT_TETRA4)
            call elrfdf('TE4', pt, dbasis)
        case (MT_PYRAM5)
            call elrfdf('PY5', pt, dbasis)
        case (MT_HEXA8)
            call elrfdf('HE8', pt, dbasis)
        case (MT_PENTA6)
            call elrfdf('PE6', pt, dbasis)
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoSplitSimplex(typema, n_simpl, indice_simpl)
!
        implicit none
!
        integer(kind=8), intent(in) :: typema
        integer(kind=8), intent(out) :: n_simpl
        integer(kind=8), intent(out) :: indice_simpl(6, 4)
!
! ---------------------------------------------------------------------------------
!  HHO - geometrie
!  Split an element in simplexe
!
! In typema : type of element
! ---------------------------------------------------------------------------------
!
!
        n_simpl = 0
        indice_simpl = 0
!
        select case (typema)
        case (MT_TRIA3)
            n_simpl = 1
            indice_simpl(1, 1:3) = [1, 2, 3]
        case (MT_QUAD4)
            n_simpl = 2
            indice_simpl(1, 1:3) = [1, 2, 3]
            indice_simpl(2, 1:3) = [1, 3, 4]
        case (MT_TETRA4)
            n_simpl = 1
            indice_simpl(1, 1:4) = [1, 2, 3, 4]
        case (MT_PYRAM5)
            n_simpl = 2
            indice_simpl(1, 1:4) = [1, 2, 3, 5]
            indice_simpl(2, 1:4) = [1, 3, 4, 5]
        case (MT_HEXA8)
            n_simpl = 5
            indice_simpl(1, 1:4) = [1, 2, 4, 5]
            indice_simpl(2, 1:4) = [2, 3, 4, 7]
            indice_simpl(3, 1:4) = [2, 4, 5, 7]
            indice_simpl(4, 1:4) = [2, 5, 6, 7]
            indice_simpl(5, 1:4) = [4, 5, 7, 8]
        case (MT_PENTA6)
            n_simpl = 3
            indice_simpl(1, 1:4) = [1, 2, 3, 4]
            indice_simpl(2, 1:4) = [2, 3, 4, 5]
            indice_simpl(3, 1:4) = [3, 4, 5, 6]
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_transfo_3d(coorno, nbnodes, typema, coorref, coorac, &
                              jacob)
!
        implicit none
!
        integer(kind=8), intent(in) :: nbnodes
        real(kind=8), dimension(3, nbnodes), intent(in) :: coorno
        integer(kind=8), intent(in) :: typema
        real(kind=8), dimension(3), intent(in) :: coorref
        real(kind=8), dimension(3), optional, intent(out) :: coorac
        real(kind=8), optional, intent(out) :: jacob
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   From reference element to current element
!   In coorno       : coordinates of the nodes
!   In coorref      : coordinates in the reference conf
!   Out coorac      : coordinates in the current conf
!   Out jacob       : determiant of the jacobienne of the transformation
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(8) :: basis
        real(kind=8), dimension(3, 8) :: dbasis
        real(kind=8), dimension(3, 3) :: jaco
        integer(kind=8) :: i, idim
!
        if (present(coorac)) then
!
! ----- shape function
!
            call hhoGeomBasis(typema, coorref, basis)
!
            coorac = 0.d0
!
            do i = 1, nbnodes
                do idim = 1, 3
                    call addFMA(coorno(idim, i), basis(i), coorac(idim))
                end do
            end do
        end if
!
        if (present(jacob)) then
!
! ----- derivative of shape function
!
            call hhoGeomDerivBasis(typema, coorref, dbasis)
!
! ---  Compute the jacobienne
            jaco = 0.d0
            do i = 1, nbnodes
                do idim = 1, 3
                    call addFMA(coorno(1, i), dbasis(idim, i), jaco(idim, 1))
                    call addFMA(coorno(2, i), dbasis(idim, i), jaco(idim, 2))
                    call addFMA(coorno(3, i), dbasis(idim, i), jaco(idim, 3))
                end do
            end do
!
            jacob = det(jaco)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_transfo_quad(coorno, coorref, ndim, coorac, jacob)
!
        implicit none
!
        real(kind=8), dimension(3, 4), intent(in) :: coorno
        real(kind=8), dimension(2), intent(in) :: coorref
        integer(kind=8), intent(in) :: ndim
        real(kind=8), optional, intent(out) :: jacob
        real(kind=8), dimension(3), intent(out), optional :: coorac
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   From reference element to current element
!   In coorno       : coordinates of the nodes
!   In coorref      : coordinates in the reference conf
!   In ndim         : topological dimenison (space where live the quad)
!   Out coorac      : coordinates in the current conf
!   Out jacob       : determiant of the jacobienne of the transformation
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), parameter :: typema = MT_QUAD4
        real(kind=8), dimension(8) :: basis
        real(kind=8), dimension(3, 8) :: dbasis
        real(kind=8), dimension(2, 2) :: jaco
        real(kind=8), dimension(3) :: da, db, normal
        integer(kind=8) :: i, idim
!
        if (present(coorac)) then
!
! ----      shape function
!
            call hhoGeomBasis(typema, (/coorref(1), coorref(2), 0.d0/), basis)
!
            coorac = 0.d0
!
            do i = 1, 4
                do idim = 1, 3
                    call addFMA(coorno(idim, i), basis(i), coorac(idim))
                end do
            end do
        end if
!
        if (present(jacob)) then
!
!       derivative of shape function
!
            call hhoGeomDerivBasis(typema, (/coorref(1), coorref(2), 0.d0/), dbasis)
!
! ---  Compute the jacobienne
            jacob = 0.d0
            select case (ndim)
            case (2)
                jaco = 0.d0
                do i = 1, 4
                    do idim = 1, 2
                        call addFMA(coorno(1, i), dbasis(idim, i), jaco(idim, 1))
                        call addFMA(coorno(2, i), dbasis(idim, i), jaco(idim, 2))
                    end do
                end do
!
                jacob = det(jaco)
            case (3)
                da = 0.d0
                db = 0.d0
                do i = 1, 4
                    do idim = 1, 3
                        call addFMA(coorno(idim, i), dbasis(1, i), da(idim))
                        call addFMA(coorno(idim, i), dbasis(2, i), db(idim))
                    end do
                end do
                call provec(da, db, normal)
                jacob = norm2(normal)
            case default
                ASSERT(ASTER_FALSE)
            end select
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function hho_jaco_cst_quad(coorno, ndim) result(l_cst)
!
        implicit none
!
        real(kind=8), dimension(3, 4), intent(in) :: coorno
        integer(kind=8), intent(in) :: ndim
        aster_logical :: l_cst
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Return True if the jacobienne is constant on cell
!   In coorno       : coordinates of the nodes
!   In ndim         : topological dimenison (space where live the quad)
!
! --------------------------------------------------------------------------------------------------
!
!
        real(kind=8) :: coor_nno(3, MT_NNOMAX), jac1, jac, tole
        integer(kind=8) :: i_node
!
        l_cst = ASTER_TRUE
        call elrfno('QU4', nodeCoor=coor_nno)
!
        call hho_transfo_quad(coorno, coor_nno(1:2, 1), ndim, jacob=jac1)
        tole = 1.d-8*abs(jac1)
        do i_node = 2, 4
            call hho_transfo_quad(coorno, coor_nno(1:2, i_node), ndim, jacob=jac)
            if (abs(jac-jac1) > tole) then
                l_cst = ASTER_FALSE
                exit
            end if
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hho_jaco_cst_3d(typema, coorno) result(l_cst)
!
        implicit none
!
        real(kind=8), dimension(3, 8), intent(in) :: coorno
        integer(kind=8), intent(in) :: typema
        aster_logical :: l_cst
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Return True if the jacobienne is constant on cell
!   In coorno       : coordinates of the nodes
!
! --------------------------------------------------------------------------------------------------
!
!
        real(kind=8) :: coor_nno(3, MT_NNOMAX), jac1, jac, tole
        character(len=8) :: sname
        integer(kind=8) :: i_node, nno
!
        l_cst = ASTER_TRUE
        call CellNameL2S(typema, sname)
        call elrfno(sname, nno=nno, nodeCoor=coor_nno)
!
        call hho_transfo_3d(coorno, nno, typema, coor_nno(1:3, 1), jacob=jac1)
        tole = 1.d-8*abs(jac1)
        do i_node = 2, nno
            call hho_transfo_3d(coorno, nno, typema, coor_nno(1:3, i_node), jacob=jac)
            if (abs(jac-jac1) > tole) then
                l_cst = ASTER_FALSE
                exit
            end if
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoIsJacobCst(typema, coorno, ndim) result(l_cst)
!
        implicit none
!
        real(kind=8), dimension(3, *), intent(in) :: coorno
        integer(kind=8), intent(in) :: typema, ndim
        aster_logical :: l_cst
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Return True if the jacobienne is constant on cell
!   In coorno       : coordinates of the nodes
!
! --------------------------------------------------------------------------------------------------
!
!
        select case (typema)
        case (MT_SEG2)
            l_cst = ASTER_TRUE
        case (MT_TRIA3)
            l_cst = ASTER_TRUE
        case (MT_QUAD4)
            l_cst = hho_jaco_cst_quad(coorno, ndim)
        case (MT_TETRA4)
            l_cst = ASTER_TRUE
        case (MT_HEXA8)
            l_cst = hho_jaco_cst_3d(typema, coorno)
        case (MT_PYRAM5)
            l_cst = hho_jaco_cst_3d(typema, coorno)
        case (MT_PENTA6)
            l_cst = hho_jaco_cst_3d(typema, coorno)
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end function
!
end module
