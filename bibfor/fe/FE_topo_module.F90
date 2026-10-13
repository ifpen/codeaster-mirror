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
module fe_topo_module
!
    use HHO_utils_module, only: CellNameS2L
!
    implicit none
!
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/apnorm.h"
#include "asterfort/assert.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfvf.h"
#include "asterfort/getFESkinSubType.h"
#include "asterfort/jevech.h"
#include "asterfort/teattr.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! FE
!
! Define types for datastructures
!
! --------------------------------------------------------------------------------------------------
!
    public :: FE_Cell, FE_Skin
    private :: normal_face, init_face, func_face, print_face, diameter_face
!
!===================================================================================================
!
!===================================================================================================
!
    type FE_Skin
! ----- Dimension topologique
        integer(kind=8)             :: ndim = 0
! ----- Type maille
        character(len=8)            :: typema = ''
! ----- Type maille short
        character(len=8)            :: typemas = ''
! ----- Nombre de noeuds
        integer(kind=8)             :: nbnodes = 0
! ----- Coordonnees des noeuds
        real(kind=8), dimension(3, MT_NNOMAX2D):: coorno = 0.d0
! ----- member function
    contains
        procedure, public, pass :: print => print_face
        procedure, public, pass :: func => func_face
        procedure, public, pass :: init => init_face
        procedure, public, pass :: normal => normal_face
        procedure, public, pass :: updateCoordinates => updateCoordinatesFace
        procedure, public, pass :: diameter => diameter_face

    end type FE_Skin
!
!===================================================================================================
!
!===================================================================================================
!
    type FE_Cell
! ----- Dimension topologique
        integer(kind=8)             :: ndim = 0
! ----- Type maille
        character(len=8)            :: typema = ''
! ----- Type maille
        character(len=8)            :: typemas = ''
! ----- Nombre de noeuds
        integer(kind=8)             :: nbnodes = 0
! ----- Coordonnees des noeuds
        real(kind=8), dimension(3, MT_NNOMAX3D):: coorno = 0.d0
! ----- member function
    contains
        procedure, public, pass :: print => print_cell
        procedure, public, pass :: init => init_cell
        procedure, public, pass :: func => func_cell
        procedure, public, pass :: evalCoor => coor_cell
        procedure, public, pass :: barycenter => bary_cell
        procedure, public, pass :: updateCoordinates
        procedure, public, pass :: getSkin
    end type FE_Cell
!
contains
!===================================================================================================
!
!===================================================================================================
!
!---------------------------------------------------------------------------------------------------
! -- member function for FE_Skin type
!---------------------------------------------------------------------------------------------------
!
    subroutine print_face(this)
!
        implicit none
!
        class(FE_Skin), intent(in) :: this
!
! --------------------------------------------------------------------------------------------------
! Print information about the face
! In this              : a FE Face
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: inode
! --------------------------------------------------------------------------------------------------
        write (6, *) "Informations on FE Face"
        write (6, *) "Type maille: ", this%typema
        write (6, *) "Dimension topo: ", this%ndim
        write (6, *) "Number of nodes: ", this%nbnodes
        do inode = 1, this%nbnodes
            write (6, *) "    node", inode, ": ", this%coorno(1:3, inode)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!---------------------------------------------------------------------------------------------------
! -- member function for FE_Cell type
!---------------------------------------------------------------------------------------------------
!
    subroutine print_cell(this)
!
        implicit none
!
        class(FE_Cell), intent(in) :: this
!
! --------------------------------------------------------------------------------------------------
!   Print information about the cell
! In this              : a FE Cell
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: inode
!
        write (6, *) "Informations on FE Cell"
        write (6, *) "Type maille: ", this%typema
        write (6, *) "Dimension topo: ", this%ndim
        write (6, *) "Number of nodes: ", this%nbnodes
        do inode = 1, this%nbnodes
            write (6, *) "    node", inode, ": ", this%coorno(1:3, inode)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine init_cell(this)
!
        implicit none
!
        class(FE_Cell), intent(out) :: this
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Initialize a FE cell
!
! --------------------------------------------------------------------------------------------------
!
! Out FECell           : a FE cell
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
        integer(kind=8) :: inode, idim, iret, jv_geom
        character(len=8) :: typ
! --------------------------------------------------------------------------------------------------
!
        call teattr('S', 'TYPMA', typ, iret)
        ASSERT(iret == 0)
!
        select case (typ)
        case ("TS1", "TS2")
            this%typemas = "TR3"
        case ("TS3", "TS4")
            this%typemas = "TR6"
        case ("TS7")
            this%typemas = "TR1"
        case ("QS1", "QS2")
            this%typemas = "QU4"
        case ("QS3", "QS4")
            this%typemas = "QU8"
        case ("QS5", "QS6")
            this%typemas = "QU9"
        case ("QS7")
            this%typemas = "Q12"
        case ("T11")
            this%typemas = "TE4"
        case ("T17")
            this%typemas = "T10"
        case ("T27")
            this%typemas = "T20"
        case ("H18")
            this%typemas = "HE8"
        case ("H26")
            this%typemas = "H20"
        case ("H34")
            this%typemas = "H27"
        case default
            this%typemas = typ
        end select
!
        call elrfno(this%typemas, nno=this%nbnodes, ndim=this%ndim)
        call CellNameS2L(this%typemas, this%typema)
!
! - Get coordinates
!
        call jevech('PGEOMER', 'L', jv_geom)
        do inode = 1, this%nbnodes
            do idim = 1, this%ndim
                this%coorno(idim, inode) = zr(jv_geom+(inode-1)*this%ndim+idim-1)
            end do
        end do
!
        if (l_debug) then
            call this%print()
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function getSkin(this, nbNodes, nodes) result(skin)
!
        implicit none
!
        class(FE_Cell), intent(in) :: this
        integer(kind=8), intent(in) :: nbNodes, nodes(MT_NNOMAX2D)
        type(FE_Skin) :: skin
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Initialize a FE Skin from a FE_Cell
!
! --------------------------------------------------------------------------------------------------
!
! Out FESkin           : a FE face
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
        integer(kind=8) :: inode, idim, node
! --------------------------------------------------------------------------------------------------
!
        skin%nbnodes = nbNodes
        skin%ndim = this%ndim-1
!
        if (this%ndim == 2) then
            select case (nbNodes)
            case (2)
                skin%typemas = "SE2"
            case (3)
                skin%typemas = "SE3"
            case (4)
                skin%typemas = "SE4"
            case default
                ASSERT(ASTER_FALSE)
            end select
        else
            select case (nbNodes)
            case (3)
                skin%typemas = "TR3"
            case (6)
                skin%typemas = "TR6"
            case (7)
                skin%typemas = "TR7"
            case (10)
                skin%typemas = "TR1"
            case (4)
                skin%typemas = "QU4"
            case (8)
                skin%typemas = "QU8"
            case (9)
                skin%typemas = "QU9"
            case (12)
                skin%typemas = "Q12"
            case default
                ASSERT(ASTER_FALSE)
            end select
        end if
!
        call CellNameS2L(skin%typemas, skin%typema)
!
! - Get coordinates
!
        do inode = 1, skin%nbnodes
            node = nodes(inode)
            do idim = 1, this%ndim
                skin%coorno(idim, inode) = this%coorno(idim, node)
            end do
        end do
!
        if (l_debug) then
            call skin%print()
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine init_face(this, side)
!
        implicit none
!
        class(FE_Skin), intent(out) :: this
        character(len=*), optional, intent(in) :: side
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Initialize a FE face
!
! --------------------------------------------------------------------------------------------------
!
! Out FESkin           : a FE face
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
        integer(kind=8) :: inode, idim, iret, jv_geom, node_init, node_idx
        character(len=8) :: side_, sub_type
! --------------------------------------------------------------------------------------------------
!
        call teattr('S', 'TYPMA', this%typemas, iret)
        ASSERT(iret == 0)
!
        node_init = 1
        if (present(side)) then
            side_ = side
            call getFESkinSubType(this%typemas, side_, sub_type, node_init)
            this%typemas = sub_type
        end if
!
        call elrfno(this%typemas, nno=this%nbnodes, ndim=this%ndim)
        call CellNameS2L(this%typemas, this%typema)
!
! - Get coordinates
!
        call jevech('PGEOMER', 'L', jv_geom)
        do inode = 1, this%nbnodes
            node_idx = node_init-1+inode-1
            do idim = 1, this%ndim+1
                this%coorno(idim, inode) = zr(jv_geom+node_idx*(this%ndim+1)+idim-1)
            end do
        end do
!
        if (l_debug) then
            call this%print()
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function coor_cell(this, pt) result(coor)
!
        implicit none
!
        class(FE_Cell), intent(in) :: this
        real(kind=8), intent(in) :: pt(3)
        real(kind=8) :: coor(3)
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Initialize a FE cell
!
! --------------------------------------------------------------------------------------------------
!
! Out FECell           : a FE cell
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: funcGeom(MT_NNOMAX)
        integer(kind=8) :: i
!
        funcGeom = this%func(pt)
        coor = 0.0
        do i = 1, this%nbnodes
            coor(1) = coor(1)+this%coorno(1, i)*funcGeom(i)
            coor(2) = coor(2)+this%coorno(2, i)*funcGeom(i)
            coor(3) = coor(3)+this%coorno(3, i)*funcGeom(i)
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function func_cell(this, pt) result(func)
!
        implicit none
!
        class(FE_Cell), intent(in) :: this
        real(kind=8), intent(in) :: pt(3)
        real(kind=8) :: func(MT_NNOMAX)
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Initialize a FE cell
!
! --------------------------------------------------------------------------------------------------
!
! Out FECell           : a FE cell
! --------------------------------------------------------------------------------------------------
!
        func = 0.d0
        call elrfvf(this%typemas, pt, func)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function func_face(this, pt) result(func)
!
        implicit none
!
        class(FE_Skin), intent(in) :: this
        real(kind=8), intent(in) :: pt(3)
        real(kind=8) :: func(MT_NNOMAX2D)
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Initialize a FE face
!
! --------------------------------------------------------------------------------------------------
!
! Out FESkin           : a FE cell
! --------------------------------------------------------------------------------------------------
!
        func = 0.d0
        call elrfvf(this%typemas, pt, func)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function bary_cell(this) result(bary)
!
        implicit none
!
        class(FE_Cell), intent(in) :: this
        real(kind=8), dimension(3) :: bary
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Compute barycenter
!
! --------------------------------------------------------------------------------------------------
!
! Out FESkin           : a FE cell
! --------------------------------------------------------------------------------------------------
!
!
        integer(kind=8) :: i_node
!
        bary = 0.d0
        do i_node = 1, this%nbnodes
            bary(1:3) = bary(1:3)+this%coorno(1:3, i_node)
        end do
!
        bary = bary/real(this%nbnodes, kind=8)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function normal_face(this, qp_param) result(normal)
!
        implicit none
!
        class(FE_Skin), intent(in) :: this
        real(kind=8), dimension(2), intent(in)    :: qp_param
        real(kind=8), dimension(3)                :: normal
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Compute normal at quadrature point
!
! --------------------------------------------------------------------------------------------------
!
! Out FESkin           : a FE cell
! --------------------------------------------------------------------------------------------------
!
!
        call apnorm(this%nbnodes, this%typemas, this%ndim+1, this%coorno, &
                    qp_param(1), qp_param(2), normal)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine updateCoordinates(this, disp)
!
        implicit none
!
        class(FE_Cell), intent(inout) :: this
        real(kind=8), intent(in) :: disp(*)
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Update coordinates with displacement
!
! --------------------------------------------------------------------------------------------------
!
! Out FECell           : a FE cell
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: inode, idim
! --------------------------------------------------------------------------------------------------
!
!
! - Update coordinates
!
        do inode = 1, this%nbnodes
            do idim = 1, this%ndim
                this%coorno(idim, inode) = this%coorno(idim, inode)+disp((inode-1)*this%ndim+idim)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine updateCoordinatesFace(this, disp)
!
        implicit none
!
        class(FE_Skin), intent(inout) :: this
        real(kind=8), intent(in) :: disp(*)
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Update coordinates with displacement
!
! --------------------------------------------------------------------------------------------------
!
! Out FECell           : a FE cell
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: inode, idim
! --------------------------------------------------------------------------------------------------
!
!
! - Update coordinates
!
        do inode = 1, this%nbnodes
            do idim = 1, this%ndim+1
                this%coorno(idim, inode) = this%coorno(idim, inode)+ &
                                           disp((inode-1)*(this%ndim+1)+idim)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function diameter_face(this)
!
        implicit none
!
        class(FE_Skin), intent(in) :: this
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic tools
!
! Initialize a FE Skin from a FE_Cell
!
! --------------------------------------------------------------------------------------------------
!
! In FESkin           : a FE face
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: inode, jnode, nnos
        real(kind=8), dimension(3) :: vector
! --------------------------------------------------------------------------------------------------
!
        diameter_face = 0.d0
!
        call elrfno(this%typemas, nnos=nnos)
!
        do inode = 1, nnos
            do jnode = inode+1, nnos
                vector(1:3) = this%coorno(1:3, inode)-this%coorno(1:3, jnode)
                diameter_face = max(norm2(vector), diameter_face)
            end do
        end do
!
    end function
!
end module
