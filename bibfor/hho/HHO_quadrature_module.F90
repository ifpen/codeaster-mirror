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
module HHO_quadrature_module
!
    use HHO_type
    use HHO_measure_module
    use HHO_geometry_module
    use HHO_utils_module, only: CellNameL2S
    use compensated_ops_module, only: sum, dot_product
!
    implicit none
!
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterf_debug.h"
#include "asterfort/assert.h"
#include "asterfort/elraga.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/reereg.h"
#include "asterfort/utmess.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic
!
! Module to generate quadratures used for HHO
!
! --------------------------------------------------------------------------------------------------
!
    type HHO_Quadrature
        integer(kind=8)                     :: order = 0
        integer(kind=8)                     :: nbQuadPoints = 0
        real(kind=8), dimension(3, MSIZE_QP)  :: points = 0.d0
        real(kind=8), dimension(MSIZE_QP)     :: weights = 0.d0
        real(kind=8), dimension(3, MSIZE_QP)  :: points_param = 0.d0
        real(kind=8), dimension(MSIZE_QP)     :: weights_param = 0.d0
        aster_logical                       :: l_point_param = ASTER_TRUE
! ----- member functions
    contains
        procedure, private, pass :: hho_edge_rules
        procedure, private, pass :: hho_tri_rules
        procedure, private, pass :: hho_quad_rules
        procedure, private, pass :: hho_tetra_rules
        procedure, private, pass :: hho_hexa_rules
        procedure, private, pass :: hho_pyram_rules
        procedure, private, pass :: hho_prism_rules
        procedure, private, pass :: hho_subtetra_rules
        procedure, private, pass :: hho_subtri_rules
        procedure, public, pass :: getQuadCell => hhoGetQuadCell
        procedure, public, pass :: getQuadFace => hhoGetQuadFace
        procedure, public, pass :: print => hhoQuadPrint
        procedure, public, pass :: initCell => hhoinitCellFamiQ
        procedure, public, pass :: initFace => hhoinitFaceFamiQ
    end type
!
!===================================================================================================
!
!===================================================================================================
!
    public   :: HHO_Quadrature
    private  :: hho_edge_rules, hho_hexa_rules, hho_tri_rules, hho_tetra_rules, &
                hho_prism_rules, hho_pyram_rules, hho_subtetra_rules, hho_subtri_rules, &
                hho_quad_rules, hhoGetQuadCell, hhoGetQuadFace, hhoQuadPrint, &
                hhoinitCellFamiQ, hhoinitFaceFamiQ, check_order, hhoSelectOrder
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine check_order(order, maxAutorized)
!
        implicit none
!
        integer(kind=8), intent(in) :: order
        integer(kind=8), intent(in) :: maxAutorized
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   control the order of integration of order has to be included in [0, maxAutorized]
!   In order  : order of integration
!   In maxAutorized : maximum autorized order
!
! --------------------------------------------------------------------------------------------------
!
        if (order > maxAutorized) then
            call utmess('F', 'HHO1_2', ni=2, vali=(/order, maxAutorized/))
        end if
!
        if (order < 0) then
            call utmess('F', 'HHO1_3', si=order)
        end if
!
    end subroutine

!
!===================================================================================================
!
!===================================================================================================
!
! --------------------------------------------------------------------------------------------------
!       !!!! BE CAREFULL ALL THE TRANSFORMATIONS ARE LINEAR (PLANAR ELEMENTS) !!!!
! --------------------------------------------------------------------------------------------------
    subroutine hho_edge_rules(this, coorno, measure, barycenter)
!
        implicit none
!
        real(kind=8), dimension(3, 2), intent(in)        :: coorno
        real(kind=8), intent(in)                        :: measure
        real(kind=8), dimension(3), intent(in)          :: barycenter
        class(HHO_quadrature), intent(inout)            :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for an edge
!   In coorno       : coordinates of the nodes
!   In measure      : length of the edge
!   In barycenter   : barycenter
!   Out this        : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: meas_2
        real(kind=8), dimension(3) :: v1
        integer(kind=8), parameter :: max_order = 15
        integer(kind=8), parameter :: max_pg = 8
        character(len=8), dimension(0:max_order) ::rules
        integer(kind=8) :: dimp, nbpg, ipg
        real(kind=8), dimension(max_pg) :: xpg, poidpg
!
! ----- check order of integration
        call check_order(this%order, max_order)
!
        rules = (/'FPG1', 'FPG1', 'FPG2', 'FPG2', 'FPG3', 'FPG3', 'FPG4', 'FPG4', &
                  'FPG5', 'FPG5', 'FPG6', 'FPG6', 'FPG7', 'FPG7', 'FPG8', 'FPG8'/)
!
        meas_2 = measure/2.d0
        v1 = (coorno(1:3, 2)-coorno(1:3, 1))/2.d0
!
!------ get quadrature points
        xpg = 0.d0
        poidpg = 0.d0
        call elraga('SE2', rules(this%order), dimp, nbpg, xpg, poidpg)
!
! ----- fill hhoQuad
        ASSERT(nbpg <= MSIZE_QP)
        this%nbQuadPoints = nbpg
!
        do ipg = 1, nbpg
            this%points_param(1, ipg) = xpg(ipg)
            this%weights_param(ipg) = poidpg(ipg)
            this%points(1:3, ipg) = barycenter+xpg(ipg)*v1
            this%weights(ipg) = meas_2*poidpg(ipg)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_quad_rules(this, coorno, ndim)
!
        implicit none
!
        real(kind=8), dimension(3, 4), intent(in)       :: coorno
        integer(kind=8), intent(in)                     :: ndim
        class(HHO_quadrature), intent(inout)            :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for an quad
!   In coorno       : coordinates of the nodes
!   In ndim         : topological dimenison (space where live the quad)
!   Out this        : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) ::  coorpgglo
        integer(kind=8), parameter :: max_order = 15
        integer(kind=8), parameter :: max_pg = 64
        character(len=8), dimension(0:max_order) ::rules
        integer(kind=8) :: dimp, nbpg, ipg
        real(kind=8) :: coorpg(max_pg*2), poidpg(max_pg), x, y, jaco
!
! ----- check order of integration
        call check_order(this%order, max_order)
!
        rules = (/'FPG1 ', 'FPG1 ', 'FPG4 ', 'FPG4 ', 'FPG9 ', 'FPG9 ', 'FPG16', 'FPG16', &
                  'FPG25', 'FPG25', 'FPG36', 'FPG36', 'FPG49', 'FPG49', 'FPG64', 'FPG64'/)
!
!------ get quadrature points
        coorpg = 0.d0
        poidpg = 0.d0
        call elraga('QU4', rules(this%order), dimp, nbpg, coorpg, poidpg)
!
! ----- fill hhoQuad
        ASSERT(nbpg <= MSIZE_QP)
        this%nbQuadPoints = nbpg
!
        do ipg = 1, nbpg
            x = coorpg(dimp*(ipg-1)+1)
            y = coorpg(dimp*(ipg-1)+2)
            call hho_transfo_quad(coorno, (/x, y/), ndim, coorpgglo, jaco)
            this%points_param(1:2, ipg) = (/x, y/)
            this%weights_param(ipg) = poidpg(ipg)
            this%points(1:3, ipg) = coorpgglo
            this%weights(ipg) = abs(jaco)*poidpg(ipg)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_hexa_rules(this, coorno)
!
        implicit none
!
        real(kind=8), dimension(3, 8), intent(in)        :: coorno
        class(HHO_quadrature), intent(inout)            :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for an hexhedra
!   In coorno       : coordinates of the nodes
!   Out this        : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: jaco
        real(kind=8), dimension(3) :: coorac
        integer(kind=8), parameter :: max_order = 15
        integer(kind=8), parameter :: max_pg = 512
        character(len=8), dimension(0:max_order) :: rules
        integer(kind=8) :: dimp, nbpg, ipg
        real(kind=8) :: coorpg(max_pg*3), poidpg(max_pg), x, y, z
!
! ----- check order of integration
        call check_order(this%order, max_order)
!
        rules = (/'FPG1  ', 'FPG1  ', 'FPG8  ', 'FPG8  ', 'FPG27 ', 'FPG27 ', 'FPG64 ', 'FPG64 ', &
                  'FPG125', 'FPG125', 'FPG216', 'FPG216', 'FPG343', 'FPG343', 'FPG512', 'FPG512'/)

!
!------ get quadrature points
        coorpg = 0.d0
        poidpg = 0.d0
        call elraga('HE8', rules(this%order), dimp, nbpg, coorpg, poidpg)
!
! ----- fill hhoQuad
        ASSERT(nbpg <= MSIZE_QP)
        this%nbQuadPoints = nbpg
!
        do ipg = 1, nbpg
            x = coorpg(dimp*(ipg-1)+1)
            y = coorpg(dimp*(ipg-1)+2)
            z = coorpg(dimp*(ipg-1)+3)
            call hho_transfo_3d(coorno, 8, MT_HEXA8, (/x, y, z/), coorac, jaco)
            this%points_param(1:3, ipg) = (/x, y, z/)
            this%weights_param(ipg) = poidpg(ipg)
            this%points(1:3, ipg) = coorac(1:3)
            this%weights(ipg) = abs(jaco)*poidpg(ipg)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_tri_rules(this, coorno, measure)
!
        implicit none
!
        real(kind=8), dimension(3, 3), intent(in)        :: coorno
        real(kind=8), intent(in)                        :: measure
        class(HHO_quadrature), intent(inout)            :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for a triangle
!   In coorno       : coordinates of the nodes
!   In measure      : surface of the triangle
!   Out this        : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) ::  coorac
        integer(kind=8), parameter :: max_order = 14
        integer(kind=8), parameter :: max_pg = 42
        character(len=8), dimension(0:max_order) ::rules
        integer(kind=8) :: dimp, nbpg, ipg, idim
        real(kind=8) :: coorpg(max_pg*2), poidpg(max_pg), x, y, basis(8), jaco
!
! ----- check order of integration
        call check_order(this%order, max_order)
!
        rules = (/'FPG1 ', 'FPG1 ', 'FPG3 ', 'FPG4 ', 'FPG6 ', 'FPG7 ', 'FPG12', 'FPG13', 'FPG16', &
                  'FPG19', 'FPG25', 'FPG28', 'FPG33', 'FPG37', 'FPG42'/)
!
!------ get quadrature points
        coorpg = 0.d0
        poidpg = 0.d0
        call elraga('TR3', rules(this%order), dimp, nbpg, coorpg, poidpg)
!
! ----- fill hhoQuad
        ASSERT(nbpg <= MSIZE_QP)
        this%nbQuadPoints = nbpg
        jaco = 2.d0*measure
!
        do ipg = 1, nbpg
            x = coorpg(dimp*(ipg-1)+1)
            y = coorpg(dimp*(ipg-1)+2)
            coorac = 0.d0
            call hhoGeomBasis(MT_TRIA3, (/x, y, 0.d0/), basis)
            do idim = 1, 3
                coorac(idim) = dot_product(coorno(idim, 1:3), basis(1:3))
            end do
            this%points_param(1:2, ipg) = (/x, y/)
            this%weights_param(ipg) = poidpg(ipg)
            this%points(1:3, ipg) = coorac
            this%weights(ipg) = jaco*poidpg(ipg)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_subtri_rules(this, typema, ndim, nbnode, coorno, param)
!
        implicit none
!
        integer(kind=8), intent(in)                     :: typema, nbnode, ndim
        real(kind=8), dimension(3, 4), intent(in)       :: coorno
        aster_logical, intent(in)                       :: param
        class(HHO_quadrature), intent(inout)            :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for a atriangle
!   In coorno       : coordinates of the nodes
!   In measure      : surface of the triangle
!   Out this        : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) ::  coorac
        integer(kind=8), parameter :: max_order = 14
        integer(kind=8), parameter :: max_pg = 42
        character(len=8), dimension(0:max_order) :: rules
        character(len=8) :: typma_s
        integer(kind=8) :: dimp, nbpg, ipg, ino, iret, idim
        integer(kind=8) :: itri, n_simp, ind_simp(6, 4)
        real(kind=8) :: coor_tri(3, 3), xe(3)
        real(kind=8) :: coorpg(max_pg*2), poidpg(max_pg), x, y, basis(8), jaco
!
! ----- split into tria
        call cellNameL2S(typema, typma_s)
        call hhoSplitSimplex(typema, n_simp, ind_simp)
        ASSERT(n_simp <= 2)
!
! ----- check order of integration
        call check_order(this%order, max_order)
!
        rules = (/'FPG1 ', 'FPG1 ', 'FPG3 ', 'FPG4 ', 'FPG6 ', 'FPG7 ', 'FPG12', 'FPG13', 'FPG16', &
                  'FPG19', 'FPG25', 'FPG28', 'FPG33', 'FPG37', 'FPG42'/)
!
!------ get quadrature points
        coorpg = 0.d0
        poidpg = 0.d0
        call elraga('TR3', rules(this%order), dimp, nbpg, coorpg, poidpg)
!
! ----- fill hhoQuad
        ASSERT(n_simp*nbpg <= MSIZE_QP)
        this%nbQuadPoints = 0
        this%l_point_param = param
!
        do itri = 1, n_simp
            do ino = 1, 3
                coor_tri(1:3, ino) = coorno(1:3, ind_simp(itri, ino))
            end do
            jaco = 2.d0*hho_surface_tri(coor_tri)
            do ipg = 1, nbpg
                x = coorpg(dimp*(ipg-1)+1)
                y = coorpg(dimp*(ipg-1)+2)
                coorac = 0.d0
                call hhoGeomBasis(MT_TRIA3, (/x, y, 0.d0/), basis)
                do idim = 1, 3
                    coorac(idim) = dot_product(coor_tri(idim, 1:3), basis(1:3))
                end do
                this%points(1:3, this%nbQuadPoints+ipg) = coorac
                this%weights(this%nbQuadPoints+ipg) = jaco*poidpg(ipg)
                if (param) then
                    call reereg("S", typma_s, nbnode, coorno, coorac, ndim, &
                                xe, iret, 1.d-8, 3)
                    this%points_param(1:2, this%nbQuadPoints+ipg) = xe(1:2)
                    this%weights_param(this%nbQuadPoints+ipg) = poidpg(ipg)
                end if
            end do
            this%nbQuadPoints = this%nbQuadPoints+nbpg
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_tetra_rules(this, coorno, measure)
!
        implicit none
!
        real(kind=8), dimension(3, 4), intent(in)       :: coorno
        real(kind=8), intent(in)                        :: measure
        class(HHO_quadrature), intent(inout)            :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for a tetrahedra
!   In coorno       : coordinates of the nodes
!   In measure      : measure of the cell
!   Out this        : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) ::  coorac
        integer(kind=8), parameter :: max_order = 14
        integer(kind=8), parameter :: max_pg = 204
        character(len=8), dimension(0:max_order) :: rules
        integer(kind=8) :: dimp, nbpg, ipg
        real(kind=8) :: coorpg(max_pg*3), poidpg(max_pg), x, y, z, jaco
!
! ----- check order of integration
        call check_order(this%order, max_order)
!
        rules = (/'FPG1  ', 'FPG1  ', 'FPG4  ', 'FPG5  ', 'FPG11 ', 'FPG15 ', 'FPG24 ', "FPG35 ", &
                  'FPG46 ', 'FPG59 ', 'FPG74 ', 'FPG94 ', 'FPG117', 'FPG144', 'FPG204'/)
!
!------ get quadrature points
        coorpg = 0.d0
        poidpg = 0.d0
        call elraga('TE4', rules(this%order), dimp, nbpg, coorpg, poidpg)
!
! ----- fill hhoQuad
        ASSERT(nbpg <= MSIZE_QP)
        this%nbQuadPoints = nbpg
        jaco = 6.d0*measure
!
        do ipg = 1, nbpg
            x = coorpg(dimp*(ipg-1)+1)
            y = coorpg(dimp*(ipg-1)+2)
            z = coorpg(dimp*(ipg-1)+3)
            call hho_transfo_3d(coorno, 4, MT_TETRA4, (/x, y, z/), coorac)
            this%points_param(1:3, ipg) = (/x, y, z/)
            this%weights_param(ipg) = poidpg(ipg)
            this%points(1:3, ipg) = coorac
            this%weights(ipg) = jaco*poidpg(ipg)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_subtetra_rules(this, typema, nbnode, coorno, param)
!
        implicit none
!
        integer(kind=8), intent(in)                             :: typema, nbnode
        real(kind=8), dimension(3, 8), intent(in)       :: coorno
        aster_logical, intent(in)                       :: param
        class(HHO_quadrature), intent(inout)            :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for a tetrahedra
!   In coorno       : coordinates of the nodes
!   Out this        : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) ::  coorac
        integer(kind=8), parameter :: max_order = 14
        integer(kind=8), parameter :: max_pg = 204
        character(len=8), dimension(0:max_order) :: rules
        character(len=8) :: typma_s
        integer(kind=8) :: dimp, nbpg, ipg, ino, iret
        integer(kind=8) :: itet, n_simp, ind_simp(6, 4)
        real(kind=8) :: coorpg(max_pg*3), poidpg(max_pg), x, y, z, jaco
        real(kind=8) :: coor_tet(3, 4), xe(3)
!
! ----- split into tetra
        call cellNameL2S(typema, typma_s)
        call hhoSplitSimplex(typema, n_simp, ind_simp)
        ASSERT(n_simp <= 6)
!
! ----- check order of integration
        call check_order(this%order, max_order)
!
        rules = (/'FPG1  ', 'FPG1  ', 'FPG4  ', 'FPG5  ', 'FPG11 ', 'FPG15 ', 'FPG24 ', "FPG35 ", &
                  'FPG46 ', 'FPG59 ', 'FPG74 ', 'FPG94 ', 'FPG117', 'FPG144', 'FPG204'/)
!
!------ get quadrature points
        coorpg = 0.d0
        poidpg = 0.d0
        call elraga('TE4', rules(this%order), dimp, nbpg, coorpg, poidpg)
!
! ----- fill hhoQuad
        ASSERT(n_simp*nbpg <= MSIZE_QP)
        this%nbQuadPoints = 0
        this%l_point_param = param
!
        do itet = 1, n_simp
            do ino = 1, 4
                coor_tet(1:3, ino) = coorno(1:3, ind_simp(itet, ino))
            end do
            jaco = 6.d0*hho_vol_tetra(coor_tet)
            do ipg = 1, nbpg
                x = coorpg(dimp*(ipg-1)+1)
                y = coorpg(dimp*(ipg-1)+2)
                z = coorpg(dimp*(ipg-1)+3)
                call hho_transfo_3d(coor_tet, 4, MT_TETRA4, (/x, y, z/), coorac)
                this%points(1:3, this%nbQuadPoints+ipg) = coorac
                this%weights(this%nbQuadPoints+ipg) = jaco*poidpg(ipg)
                if (param) then
                    call reereg("S", typma_s, nbnode, coorno, coorac, 3, &
                                xe, iret, 1.d-8, 3)
                    this%points_param(1:3, this%nbQuadPoints+ipg) = xe(1:3)
                    this%weights_param(this%nbQuadPoints+ipg) = poidpg(ipg)
                end if
            end do
            this%nbQuadPoints = this%nbQuadPoints+nbpg
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_pyram_rules(this, coorno)
!
        implicit none
!
        real(kind=8), dimension(3, 5), intent(in)        :: coorno
        class(HHO_quadrature), intent(inout)            :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for a pyramid
!   In coorno       : coordinates of the nodes
!   Out this        : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: jaco, x, y, z
        real(kind=8), dimension(3) :: coorac
        integer(kind=8), parameter :: max_order = 10
        integer(kind=8), parameter :: max_pg = 83
        character(len=8), dimension(0:max_order) :: rules
        integer(kind=8) :: dimp, nbpg, ipg
        real(kind=8) :: coorpg(max_pg*3), poidpg(max_pg)
!
! ----- check order of integration
        call check_order(this%order, max_order)
!
        rules = (/'FPG1B', 'FPG1B', 'FPG5 ', 'FPG6 ', 'FPG10', 'FPG15', 'FPG24', 'FPG31', &
                  'FPG47', 'FPG62', 'FPG83'/)
!
!------ get quadrature points
        coorpg = 0.d0
        poidpg = 0.d0
        call elraga('PY5', rules(this%order), dimp, nbpg, coorpg, poidpg)
!
! ----- fill hhoQuad
        ASSERT(nbpg <= MSIZE_QP)
        this%nbQuadPoints = nbpg
!
        do ipg = 1, nbpg
            x = coorpg(dimp*(ipg-1)+1)
            y = coorpg(dimp*(ipg-1)+2)
            z = coorpg(dimp*(ipg-1)+3)
            call hho_transfo_3d(coorno, 5, MT_PYRAM5, (/x, y, z/), coorac, jaco)
            this%points_param(1:3, ipg) = (/x, y, z/)
            this%weights_param(ipg) = poidpg(ipg)
            this%points(1:3, ipg) = coorac(1:3)
            this%weights(ipg) = abs(jaco)*poidpg(ipg)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_prism_rules(this, coorno)
!
        implicit none
!
        real(kind=8), dimension(3, 6), intent(in)        :: coorno
        class(HHO_quadrature), intent(inout)            :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for a prism
!   In coorno       : coordinates of the nodes
!   Out this        : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: jaco
        real(kind=8), dimension(3) :: coorac
        integer(kind=8), parameter :: max_order = 11
        integer(kind=8), parameter :: max_pg = 168
        character(len=8), dimension(0:max_order) :: rules
        integer(kind=8) :: dimp, nbpg, ipg
        real(kind=8) :: coorpg(max_pg*3), poidpg(max_pg), x, y, z
!
! ----- check order of integration
        call check_order(this%order, max_order)
!
        rules = (/'FPG1  ', 'FPG6B ', 'FPG6B ', 'FPG8  ', 'FPG21 ', 'FPG21 ', 'FPG29 ', &
                  'FPG52 ', 'FPG95 ', 'FPG95 ', 'FPG168', 'FPG168'/)
!
!------ get quadrature points
        coorpg = 0.d0
        poidpg = 0.d0
        call elraga('PE6', rules(this%order), dimp, nbpg, coorpg, poidpg)
!
! ----- fill hhoQuad
        ASSERT(nbpg <= MSIZE_QP)
        this%nbQuadPoints = nbpg
!
        do ipg = 1, nbpg
            x = coorpg(dimp*(ipg-1)+1)
            y = coorpg(dimp*(ipg-1)+2)
            z = coorpg(dimp*(ipg-1)+3)
            call hho_transfo_3d(coorno, 6, MT_PENTA6, (/x, y, z/), coorac, jaco)
            this%points_param(1:3, ipg) = (/x, y, z/)
            this%weights_param(ipg) = poidpg(ipg)
            this%points(1:3, ipg) = coorac(1:3)
            this%weights(ipg) = abs(jaco)*poidpg(ipg)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGetQuadCell(this, hhoCell, order, split, param, adapt)
!
        implicit none
!
        type(HHO_cell), intent(in)            :: hhoCell
        integer(kind=8), intent(in)           :: order
        aster_logical, intent(in), optional   :: split, param, adapt
        class(HHO_quadrature), intent(out)    :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for the current cell
!   In hhoCell      : a HHO cell
!   In order        : quadrature order
!   Out this        : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ipg
        aster_logical :: split_simpl, param_, adapt_, axis_
!
        this%order = order
!
        adapt_ = ASTER_TRUE
        if (present(adapt)) then
            adapt_ = adapt
        end if
!
        axis_ = hhoCell%l_axis
        if (axis_ .and. adapt_) then
            this%order = this%order+1
        end if
!
        if (present(split)) then
            split_simpl = split
        else
            split_simpl = ASTER_TRUE
        end if
!
        if (hhoCell%l_jaco_cst) then
            split_simpl = ASTER_FALSE
        end if
!
        if (split_simpl) then
!
            if (present(param)) then
                param_ = param
            else
                param_ = ASTER_FALSE
            end if
!
            if (hhoCell%ndim == 3) then
                call this%hho_subtetra_rules(hhoCell%typema, hhoCell%nbnodes, &
                                             hhoCell%coorno(1:3, 1:8), param_)
            elseif (hhoCell%ndim == 2) then
                call this%hho_subtri_rules(hhoCell%typema, 2, hhoCell%nbnodes, &
                                           hhoCell%coorno(1:3, 1:4), param_)
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            select case (hhoCell%typema)
            case (MT_HEXA8)
                call this%hho_hexa_rules(hhoCell%coorno(1:3, 1:8))
            case (MT_TETRA4)
                call this%hho_tetra_rules(hhoCell%coorno(1:3, 1:4), hhoCell%measure)
            case (MT_PYRAM5)
                call this%hho_pyram_rules(hhoCell%coorno(1:3, 1:5))
            case (MT_PENTA6)
                call this%hho_prism_rules(hhoCell%coorno(1:3, 1:6))
            case (MT_QUAD4)
                call this%hho_quad_rules(hhoCell%coorno(1:3, 1:4), 2)
            case (MT_TRIA3)
                call this%hho_tri_rules(hhoCell%coorno(1:3, 1:3), hhoCell%measure)
            case default
                ASSERT(ASTER_FALSE)
            end select
        end if
!
        if (axis_) then
            do ipg = 1, this%nbQuadPoints
                this%weights(ipg) = this%weights(ipg)*this%points(1, ipg)
            end do
        end if
!
        ASSERT(this%nbQuadPoints <= MSIZE_QP)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGetQuadFace(this, hhoFace, order, split, param, adapt)
!
        implicit none
!
        type(HHO_face), intent(in)          :: hhoFace
        integer(kind=8), intent(in)         :: order
        aster_logical, intent(in), optional :: adapt, split, param
        class(HHO_quadrature), intent(out)  :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules for the current face
!   In hhoFace      : a HHO face
!   In order        : quadrature order
!   Out this     : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ipg
        aster_logical :: split_simpl, param_, adapt_, axis_
!
        this%order = order
!
        adapt_ = ASTER_TRUE
        if (present(adapt)) then
            adapt_ = adapt
        end if
!
        axis_ = hhoFace%l_axis
        if (axis_) then
            ! We modifiy quadrature for edges on axe
            ! Do no multiply by r since it remains an
            ! edge after revolution
            axis_ = .not. hhoFace%l_axis_on_axe
        end if
!
        if (axis_ .and. adapt_) then
            this%order = this%order+1
        end if
!
        if (present(split)) then
            split_simpl = split
        else
            split_simpl = ASTER_TRUE
        end if
!
        if (hhoFace%l_jaco_cst) then
            split_simpl = ASTER_FALSE
        end if
!
        if (split_simpl) then
!
            if (present(param)) then
                param_ = param
            else
                param_ = ASTER_FALSE
            end if
!
            if (hhoFace%ndim == 2) then
                call this%hho_subtri_rules(hhoFace%typema, 3, hhoFace%nbnodes, &
                                           hhoFace%coorno(1:3, 1:4), param_)
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            select case (hhoFace%typema)
            case (MT_QUAD4)
                call this%hho_quad_rules(hhoFace%coorno(1:3, 1:4), 3)
            case (MT_TRIA3)
                call this%hho_tri_rules(hhoFace%coorno(1:3, 1:3), hhoFace%measure)
            case (MT_SEG2)
                call this%hho_edge_rules(hhoFace%coorno(1:3, 1:2), hhoFace%measure, &
                                         hhoFace%barycenter)
            case default
                ASSERT(ASTER_FALSE)
            end select
        end if
!
        if (axis_) then
            do ipg = 1, this%nbQuadPoints
                this%weights(ipg) = this%weights(ipg)*this%points(1, ipg)
            end do
        end if
!
        ASSERT(this%nbQuadPoints <= MSIZE_QP)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoSelectOrder(typema, npg, order)
!
        implicit none
!
        integer(kind=8), intent(in)  :: typema
        integer(kind=8), intent(in)           :: npg
        integer(kind=8), intent(out)          :: order
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the order of the quadrature rules from a familly definied in the catalogue of code_aster
!   In typema   : type of Cell or Face
!   In npg      : number of quadrature points
!   Out order   : order of the quadrature
!
! --------------------------------------------------------------------------------------------------
!
        order = 0
!
        select case (typema)
        case (MT_HEXA8)
            select case (npg)
            case (1)
                order = 1
            case (8)
                order = 3
            case (27)
                order = 5
            case (64)
                order = 7
            case (125)
                order = 9
            case (216)
                order = 11
            case (343)
                order = 13
            case default
                ASSERT(ASTER_FALSE)
            end select
        case (MT_PENTA6)
            select case (npg)
            case (1)
                order = 0
            case (6)
                order = 2
            case (8)
                order = 3
            case (21)
                order = 5
            case (29)
                order = 6
            case (52)
                order = 6
            case (95)
                order = 9
            case (168)
                order = 11
            case default
                ASSERT(ASTER_FALSE)
            end select
        case (MT_PYRAM5)
            select case (npg)
            case (1)
                order = 1
            case (5)
                order = 2
            case (6)
                order = 3
            case (10)
                order = 4
            case (15)
                order = 5
            case (24)
                order = 6
            case (31)
                order = 7
            case (47)
                order = 8
            case (62)
                order = 9
            case (83)
                order = 10
            case default
                ASSERT(ASTER_FALSE)
            end select
        case (MT_TETRA4)
            select case (npg)
            case (1)
                order = 1
            case (4)
                order = 2
            case (5)
                order = 3
            case (11)
                order = 4
            case (15)
                order = 5
            case (24)
                order = 6
            case (35)
                order = 7
            case (46)
                order = 8
            case (59)
                order = 9
            case (74)
                order = 10
            case (94)
                order = 11
            case (117)
                order = 12
            case (144)
                order = 13
            case (204)
                order = 14
            case default
                ASSERT(ASTER_FALSE)
            end select
        case (MT_QUAD4)
            select case (npg)
            case (1)
                order = 1
            case (4)
                order = 3
            case (9)
                order = 5
            case (16)
                order = 7
            case (25)
                order = 9
            case (36)
                order = 11
            case (49)
                order = 13
            case default
                ASSERT(ASTER_FALSE)
            end select
        case (MT_TRIA3)
            select case (npg)
            case (1)
                order = 1
            case (3)
                order = 2
            case (4)
                order = 3
            case (6)
                order = 4
            case (7)
                order = 5
            case (12)
                order = 6
            case (13)
                order = 7
            case (16)
                order = 8
            case (19)
                order = 9
            case (25)
                order = 10
            case (28)
                order = 11
            case default
                ASSERT(ASTER_FALSE)
            end select
        case (MT_SEG2)
            select case (npg)
            case (1)
                order = 1
            case (2)
                order = 3
            case (3)
                order = 5
            case (4)
                order = 7
            case (5)
                order = 9
            case (6)
                order = 11
            case default
                ASSERT(ASTER_FALSE)
            end select
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
    subroutine hhoinitCellFamiQ(this, hhoCell, npg)
!
        implicit none
!
        type(HHO_cell), intent(in)          :: hhoCell
        integer(kind=8), intent(in)         :: npg
        class(HHO_quadrature), intent(out)  :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules from a familly definied in the catalogue of code_aster
!   In hhoCell  : hhoCell
!   In npg      : number of quadrature points
!   Out this    : hho quadrature
!
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: order
!
        ASSERT(npg .le. MSIZE_QP)
        this%nbQuadPoints = npg
!
        call hhoSelectOrder(hhoCell%typema, npg, order)
!
        call hhoGetQuadCell(this, hhoCell, order, split=ASTER_FALSE, adapt=ASTER_FALSE)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoinitFaceFamiQ(this, hhoFace, npg)
!
        implicit none
!
        type(HHO_Face), intent(in)          :: hhoFace
        integer(kind=8), intent(in)         :: npg
        class(HHO_quadrature), intent(out)  :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Get the quadrature rules from a familly definied in the catalogue of code_aster
!   In hhoFace  : hhoFace
!   In npg      : number of quadrature points
!   Out this    : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: order
!
        ASSERT(npg .le. MSIZE_QP)
        this%nbQuadPoints = npg
!
        call hhoSelectOrder(hhoFace%typema, npg, order)
!
        call hhoGetQuadFace(this, hhoFace, order, &
                            split=ASTER_FALSE, adapt=ASTER_FALSE)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoQuadPrint(this)
!
        implicit none
!
        class(HHO_quadrature), intent(in)  :: this
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Print quadrature informations
!   In this    : hho quadrature
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ipg
!
        write (6, *) "QUADRATURE INFORMATIONS"
        write (6, *) "number of qp: ", this%nbQuadPoints
        write (6, *) "order: ", this%order
!
        do ipg = 1, this%nbQuadPoints
            write (6, *) "coordo qp ", ipg, ":", this%points(1:3, ipg)
            write (6, *) "weight qp ", ipg, ":", this%weights(ipg)
        end do
        write (6, *) "END QUADRATURE INFORMATIONS"
!
    end subroutine
!
end module
