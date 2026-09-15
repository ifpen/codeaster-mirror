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
module HHO_type
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/tecach.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Define types for datastructures
!
! --------------------------------------------------------------------------------------------------
!
    public :: HHO_Field, HHO_Data, HHO_Cell, HHO_Face
!
    type HHO_Field
!
        aster_logical      :: l_debug = ASTER_FALSE
! ----- Fields for Dirichlet loads
        aster_logical      :: l_cine_f = ASTER_FALSE
        character(len=19)  :: fieldCineFunc = ''
        character(len=19)  :: fieldCineVale = ''
        integer(kind=8), pointer :: v_info_cine(:) => null()
!
    end type HHO_Field
!
!===================================================================================================
!
!===================================================================================================
!
    type HHO_Data
        aster_logical      :: l_debug = ASTER_FALSE
! ----- Ordre inconnues face
        integer(kind=8)            :: m_face_degree
! ----- Ordre inconnues cellule
        integer(kind=8)            :: m_cell_degree
! ----- Ordre inconnues grad
        integer(kind=8)            :: m_grad_degree
! ----- Adapt stabilization coefficient ?
        aster_logical      :: l_adapt_coef
! ----- Coefficient stabilisation
        real(kind=8)       :: m_coeff_stab
! ----- member function
    contains
        procedure, pass    :: initialize => initialize_data
        procedure, pass    :: face_degree
        procedure, pass    :: cell_degree
        procedure, pass    :: grad_degree
        procedure, pass    :: debug => debug_data
        procedure, pass    :: adapt
        procedure, pass    :: coeff_stab
        procedure, pass    :: precompute
        procedure, pass    :: setCoeffStab
        procedure, pass    :: print => print_data
    end type HHO_Data
!
!===================================================================================================
!
!===================================================================================================
!
    type HHO_Face
! ----- Dimension topologique
        integer(kind=8)              :: ndim = 0
! ----- Type maille
        integer(kind=8)              :: typema = 0
! ----- Nombre de noeuds
        integer(kind=8)              :: nbnodes = 0
        integer(kind=8)              :: nbnodes_post = 0
! ----- Coordonnees des noeuds   (max 9 noeuds pour quad)
        real(kind=8), dimension(3, 9):: coorno = 0.d0
! ----- Coordonnees du barycentre de la face
        real(kind=8), dimension(3)  :: barycenter = 0.d0
! ----- Diametre de la face
        real(kind=8)                :: diameter = 0.d0
! ----- Surface ou longueur de la face
        real(kind=8)                :: measure = 0.d0
! ----- Normale sortante au barycentre
        real(kind=8), dimension(3)  :: normal = 0.d0
! ----- Axes locaux de la cellule
        real(kind=8), dimension(3, 2) :: axes = 0.d0
! ----- Index locale des noeuds de la face
        integer(kind=8), dimension(9) :: nodes_loc = 0
        integer(kind=8)               :: node_bar_loc = 0
! ----- Index locale de la face pour une cellule
        integer(kind=8)               :: face_loc = 0
! ----- Jacobien is constant
        aster_logical               :: l_jaco_cst = ASTER_FALSE
! ----- For AXIS, the edge is on the axis r=0
        aster_logical               :: l_axis_on_axe = ASTER_FALSE
! ----- Axis problem
        aster_logical                :: l_axis = ASTER_FALSE
! ----- ReOrientation (True if local renumbering change orientation of the normal )
        aster_logical                :: reorientation = ASTER_FALSE
! ----- member function
    contains
        procedure, public, pass :: print => print_face
    end type HHO_Face
!
!===================================================================================================
!
!===================================================================================================
!
    type HHO_Cell
! ----- Dimension topologique
        integer(kind=8)                     :: ndim = 0
! ----- Type maille
        integer(kind=8)                     :: typema = 0
! ----- Nombre de noeuds
        integer(kind=8)                     :: nbnodes = 0
! ----- Coordonnees des noeuds   (max 27 noeuds pour hexa)
        real(kind=8), dimension(3, 27):: coorno = 0.d0
! ----- Coordonnees du barycentre de la cellule
        real(kind=8), dimension(3)  :: barycenter = 0.d0
! ----- Diametre de la cellule
        real(kind=8)                :: diameter = 0.d0
! ----- Axes locaux de la cellule
        real(kind=8), dimension(3, 3) :: axes = 0.d0
! ----- Volume ou Surface de la cellule
        real(kind=8)                :: measure = 0.d0
! ----- Nombre de faces de la cellule
        integer(kind=8)             :: nbfaces = 0
! ----- Donnees sur les faces (max 6 faces pour hexa)
        type(HHO_Face), dimension(6) :: faces
! ----- Index locale des noeuds de la cellule
        integer(kind=8)              :: node_bar_loc = 0
! ----- Jacobien is constant
        aster_logical                :: l_jaco_cst = ASTER_FALSE
! ----- Initialize face
        aster_logical                :: l_face_init = ASTER_FALSE
! ----- Axis problem
        aster_logical                :: l_axis = ASTER_FALSE
! ----- member function
    contains
        procedure, public, pass :: print => print_cell
    end type HHO_Cell
!
contains
!
!---------------------------------------------------------------------------------------------------
! -- member function for HHO_Data type
!---------------------------------------------------------------------------------------------------
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_data(this, face_degree, cell_degree, grad_degree, coeff, &
                               l_debug, l_adapt_coef)
!
        implicit none
!
        class(HHO_Data), intent(inout)  :: this
        integer(kind=8), intent(in)             :: face_degree
        integer(kind=8), intent(in)             :: cell_degree
        integer(kind=8), intent(in)             :: grad_degree
        aster_logical, intent(in)       :: l_adapt_coef
        real(kind=8), intent(in)        :: coeff
        aster_logical, intent(in)       :: l_debug
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a HHO_DATA type
!   In this     : a HHo Data
!   In face_degree
!   In cell_degree
!   In grad_degree
!   In l_adapt_coef: adapt automatically the coefficient yes or no ?
!   In coeff    : coefficient of stabilization
!   In l_debug  : debug yes or no ?
! --------------------------------------------------------------------------------------------------
!
        this%m_face_degree = face_degree
        this%m_cell_degree = cell_degree
        this%m_grad_degree = grad_degree
        this%l_adapt_coef = l_adapt_coef
        this%m_coeff_stab = coeff
        this%l_debug = l_debug
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function face_degree(this) result(degree)
!
        implicit none
!
        class(HHO_Data), intent(in) :: this
        integer(kind=8)             :: degree
!
! --------------------------------------------------------------------------------------------------
!
!   return the face degree of a HHO_DATA type
!   In this     : a HHo Data
!   Out degree  : face_degree
! --------------------------------------------------------------------------------------------------
!
        degree = this%m_face_degree
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function cell_degree(this) result(degree)
!
        implicit none
!
        class(HHO_Data), intent(in) :: this
        integer(kind=8)                     :: degree
!
! --------------------------------------------------------------------------------------------------
!
!   return the cell degree of a HHO_DATA type
!   In this     : a HHo Data
!   Out degree  : cell_degree
! --------------------------------------------------------------------------------------------------
!
        degree = this%m_cell_degree
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function grad_degree(this) result(degree)
!
        implicit none
!
        class(HHO_Data), intent(in) :: this
        integer(kind=8)                     :: degree
!
! --------------------------------------------------------------------------------------------------
!
!   return the grad degree of a HHO_DATA type
!   In this     : a HHo Data
!   Out degree  : grad_degree
! --------------------------------------------------------------------------------------------------
!
        degree = this%m_grad_degree
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function debug_data(this) result(logic)
!
        implicit none
!
        class(HHO_Data), intent(in) :: this
        logical                     :: logic
!
! --------------------------------------------------------------------------------------------------
!
!   return the debug logical of a HHO_DATA type
!   In this     : a HHo Data
!   Out logic   : debug yes or no ?
! --------------------------------------------------------------------------------------------------
!
        logic = this%l_debug
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function adapt(this) result(logic)
!
        implicit none
!
        class(HHO_Data), intent(in) :: this
        logical                     :: logic
!
! --------------------------------------------------------------------------------------------------
!
!   adapt automatically the stabilization parameter
!   In this     : a HHo Data
!   Out logic   : adapt yes or no ?
! --------------------------------------------------------------------------------------------------
!
        logic = this%l_adapt_coef
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function precompute(this) result(precomp)
!
        implicit none
!
        class(HHO_Data), intent(in) :: this
        aster_logical                :: precomp
!
! --------------------------------------------------------------------------------------------------
!
!   Operators are precomputed ?
!   In this     : a HHo Data
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: iret, itab(8)
!
        call tecach('NNO', 'PCHHOGT', 'L', iret, nval=1, itab=itab)
        precomp = iret == 0
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function coeff_stab(this) result(coeff)
!
        implicit none
!
        class(HHO_Data), intent(in) :: this
        real(kind=8)                :: coeff
!
! --------------------------------------------------------------------------------------------------
!
!   return the stabilization coefficient of a HHO_DATA type
!   In this     : a HHo Data
!   Out coeff   : coefficient of stabilization
! --------------------------------------------------------------------------------------------------
!
        coeff = this%m_coeff_stab
!
    end function
!===================================================================================================
!
!===================================================================================================
!
    subroutine setCoeffStab(this, coeff)
!
        implicit none
!
        class(HHO_Data), intent(inout) :: this
        real(kind=8), intent(in)       :: coeff
!
! --------------------------------------------------------------------------------------------------
!
!   give the stabilization coefficient of a HHO_DATA type
!   In this     : a HHo Data
!   In coeff   : coefficient of stabilization
! --------------------------------------------------------------------------------------------------
!
        this%m_coeff_stab = coeff
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine print_data(this)
!
        implicit none
!
        class(HHO_Data), intent(in) :: this
!
! --------------------------------------------------------------------------------------------------
!
!   print informations of a HHO_DATA type
!   In this     : a HHo Data
! --------------------------------------------------------------------------------------------------
!
        write (6, *) "hhoData debug"
        write (6, *) "face degree: ", this%face_degree()
        write (6, *) "cell degree: ", this%cell_degree()
        write (6, *) "grad degree: ", this%grad_degree()
        write (6, *) "coeff stab: ", this%coeff_stab()
        write (6, *) "end hhoData debug"
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!---------------------------------------------------------------------------------------------------
! -- member function for HHO_Face type
!---------------------------------------------------------------------------------------------------
!
    subroutine print_face(this)
!
        implicit none
!
        class(HHO_Face), intent(in) :: this
!
! --------------------------------------------------------------------------------------------------
! Print information about the face
! In this              : a HHO Face
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: inode
! --------------------------------------------------------------------------------------------------
        write (6, *) "Informations on HHO Face"
        write (6, *) "Type maille: ", this%typema
        write (6, *) "Dimension topo: ", this%ndim
        write (6, *) "Number of nodes: ", this%nbnodes
        do inode = 1, this%nbnodes
            write (6, *) "    node", inode, ": ", this%coorno(1:3, inode)
        end do
        write (6, *) "Barycenter: ", this%barycenter
        write (6, *) "Normal: ", this%normal
        write (6, *) "Measure: ", this%measure
        write (6, *) "Diameter: ", this%diameter
        write (6, *) "Cst Jacobien: ", this%l_jaco_cst
        write (6, *) "Local axis: "
        write (6, *) "    a1: ", this%axes(1:3, 1)
        if (this%ndim > 1) then
            write (6, *) "    a2: ", this%axes(1:3, 2)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!---------------------------------------------------------------------------------------------------
! -- member function for HHO_Cell type
!---------------------------------------------------------------------------------------------------
!
    subroutine print_cell(this)
!
        implicit none
!
        class(HHO_Cell), intent(in) :: this
!
! --------------------------------------------------------------------------------------------------
!   Print information about the cell
! In this              : a HHO Cell
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: inode, iface
!
        write (6, *) "Informations on HHO Cell"
        write (6, *) "Type maille: ", this%typema
        write (6, *) "Dimension topo: ", this%ndim
        write (6, *) "Number of nodes: ", this%nbnodes
        do inode = 1, this%nbnodes
            write (6, *) "    node", inode, ": ", this%coorno(1:3, inode)
        end do
        write (6, *) "Barycenter: ", this%barycenter
        write (6, *) "Measure: ", this%measure
        write (6, *) "Diameter: ", this%diameter
        write (6, *) "Cst Jacobien: ", this%l_jaco_cst
        write (6, *) "Local axis: "
        write (6, *) "    a1: ", this%axes(1:3, 1)
        if (this%ndim > 1) then
            write (6, *) "    a2: ", this%axes(1:3, 2)
            if (this%ndim > 2) then
                write (6, *) "    a3: ", this%axes(1:3, 3)
            end if
        end if
        write (6, *) "Number of face: ", this%nbfaces
        do iface = 1, this%nbfaces
            write (6, *) "    face", iface, ": "
            call this%faces(iface)%print()
        end do
!
    end subroutine
!
end module
