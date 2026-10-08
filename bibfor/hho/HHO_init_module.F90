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
module HHO_init_module
!
    use HHO_type
    use HHO_geometry_module
    use HHO_measure_module
    use HHO_quadrature_module
    use HHO_basis_module
    use HHO_utils_module
    use HHO_inertia_module
!
    implicit none
!
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/rcvala.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic
!
! generic routines to initialize data for HHO model
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: hhoInfoInitCell, hhoInfoInitFace, hhoInitFacesOfCell, hhoInfoInitCellAndFace
    public :: hhoUpdateCoordinatesFace, hhoCellInit
    private :: hhoGeomData, hhoGeomFace, hhoDataInit, hhoFaceInit
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGeomData(nodes_coor, nbnodes, typma, elem_dim)
!
        implicit none
!
        real(kind=8), dimension(3, 27), intent(out)  :: nodes_coor
        integer(kind=8), intent(out)                 :: nbnodes
        character(len=8), intent(out)                :: typma
        integer(kind=8), intent(out)                 :: elem_dim
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tools
!
! Get geometric data of the current element
!
! --------------------------------------------------------------------------------------------------
!
! Out nodes_coor        : coordinates of the nodes
! Out nbnodes           : number of nodes
! Out typma             : type of the element
! Out elem_dim          : dimension of the element
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
        integer(kind=8) :: jv_geom
        integer(kind=8) :: inode, idim, iret
! --------------------------------------------------------------------------------------------------
!
! --- Init
        nodes_coor = 0.d0
        nbnodes = 0
        typma = ''
        elem_dim = 0
!
        call teattr('S', 'TYPMA', typma, iret)
        ASSERT(iret == 0)
!
        if (typma == 'H27') then
            typma = 'HEXA27'
            nbnodes = 27
            elem_dim = 3
        elseif (typma == 'QU8') then
            typma = 'QUAD8'
            nbnodes = 8
            elem_dim = 2
        elseif (typma == 'QU9') then
            typma = 'QUAD9'
            nbnodes = 9
            elem_dim = 2
        elseif (typma == 'TR6') then
            typma = 'TRIA6'
            nbnodes = 6
            elem_dim = 2
        elseif (typma == 'TR7') then
            typma = 'TRIA7'
            nbnodes = 7
            elem_dim = 2
        elseif (typma == 'T15') then
            typma = 'TETRA15'
            nbnodes = 15
            elem_dim = 3
        elseif (typma == 'P19') then
            typma = 'PYRAM19'
            nbnodes = 19
            elem_dim = 3
        elseif (typma == 'P21') then
            typma = 'PENTA21'
            nbnodes = 21
            elem_dim = 3
        else
            ASSERT(ASTER_FALSE)
        end if
!
        ASSERT((nbnodes .ge. 4) .and. (nbnodes .le. 27))
        ASSERT((elem_dim .eq. 2) .or. (elem_dim .eq. 3))
!
        call jevech('PGEOMER', 'L', jv_geom)
!
! --- Get coordinates
!
        do inode = 1, nbnodes
            do idim = 1, elem_dim
                nodes_coor(idim, inode) = zr(jv_geom+(inode-1)*elem_dim+idim-1)
            end do
        end do
!
        if (l_debug) then
            write (6, *) "hhoGeomData debug"
            write (6, *) "typma: ", typma
            write (6, *) "nbnodes: ", nbnodes
            do inode = 1, nbnodes
                write (6, *) "node ", inode, ": ", nodes_coor(1:3, inode)
            end do
            write (6, *) "end hhoGeomData debug"
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGeomFace(nodes_coor, nbnodes, nbnodes_post, typma, elem_dim)
!
        implicit none
!
        real(kind=8), dimension(3, 9), intent(out)  :: nodes_coor
        integer(kind=8), intent(out)                        :: nbnodes, nbnodes_post
        character(len=8), intent(out)               :: typma
        integer(kind=8), intent(out)                        :: elem_dim
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tools
!
! Get geometric data of the current element
!
! --------------------------------------------------------------------------------------------------
!
! Out nodes_coor        : coordinates of the nodes
! Out nbnodes           : number of nodes
! Out typma             : type of the element
! Out elem_dim          : dimension of the element
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
        integer(kind=8) :: jv_geom, node_i
        integer(kind=8) :: inode, idim, iret, node_init
!
! --- Init
        nodes_coor = 0.d0
        nbnodes = 0
        nbnodes_post = 0
        typma = ''
        elem_dim = 0
        node_init = 1
!
        call teattr('S', 'TYPMA', typma, iret)
        ASSERT(iret == 0)
!
        if (typma == 'QU9') then
            typma = 'QUAD4'
            nbnodes = 4
            nbnodes_post = 9
            elem_dim = 2
        elseif (typma == 'TR7') then
            typma = 'TRIA3'
            nbnodes = 3
            nbnodes_post = 7
            elem_dim = 2
        elseif (typma == 'SE3') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            elem_dim = 1
        elseif (typma == 'TS2') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            node_init = 4
            elem_dim = 1
        elseif (typma == 'TS4') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            node_init = 7
            elem_dim = 1
        elseif (typma == 'TS7') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            node_init = 11
            elem_dim = 1
        elseif (typma == 'QS2') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            node_init = 5
            elem_dim = 1
        elseif (typma == 'QS4') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            node_init = 9
            elem_dim = 1
        elseif (typma == 'QS6') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            node_init = 10
            elem_dim = 1
        elseif (typma == 'QS7') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            node_init = 13
            elem_dim = 1
        elseif (typma == 'S23') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            node_init = 3
            elem_dim = 1
        elseif (typma == 'S33') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            node_init = 4
            elem_dim = 1
        elseif (typma == 'S43') then
            typma = 'SEG2'
            nbnodes = 2
            nbnodes_post = 3
            node_init = 5
            elem_dim = 1
        elseif (typma == 'QQ6') then
            typma = 'QUAD4'
            nbnodes = 4
            nbnodes_post = 9
            node_init = 5
            elem_dim = 2
        elseif (typma == 'QQ3') then
            typma = 'QUAD4'
            nbnodes = 4
            nbnodes_post = 9
            node_init = 9
            elem_dim = 2
        elseif (typma == 'Q99') then
            typma = 'QUAD4'
            nbnodes = 4
            nbnodes_post = 9
            node_init = 10
            elem_dim = 2
        elseif (typma == 'TQ5') then
            typma = 'QUAD4'
            nbnodes = 4
            nbnodes_post = 9
            node_init = 4
            elem_dim = 2
        elseif (typma == 'TQ3') then
            typma = 'QUAD4'
            nbnodes = 4
            nbnodes_post = 9
            node_init = 7
            elem_dim = 2
        elseif (typma == 'QT8') then
            typma = 'TRIA3'
            nbnodes = 3
            nbnodes_post = 7
            node_init = 5
            elem_dim = 2
        elseif (typma == 'QT9') then
            typma = 'TRIA3'
            nbnodes = 3
            nbnodes_post = 7
            node_init = 9
            elem_dim = 2
        elseif (typma == 'QT0') then
            typma = 'TRIA3'
            nbnodes = 3
            nbnodes_post = 7
            node_init = 10
            elem_dim = 2
        elseif (typma == 'TT6') then
            typma = 'TRIA3'
            nbnodes = 3
            nbnodes_post = 7
            node_init = 4
            elem_dim = 2
        elseif (typma == 'TT7') then
            typma = 'TRIA3'
            nbnodes = 3
            nbnodes_post = 7
            node_init = 7
            elem_dim = 2
        elseif (typma == 'T11') then
            typma = 'TRIA3'
            nbnodes = 3
            nbnodes_post = 7
            node_init = 5
            elem_dim = 2
        elseif (typma == 'T17') then
            typma = 'TRIA3'
            nbnodes = 3
            nbnodes_post = 7
            node_init = 11
            elem_dim = 2
        elseif (typma == 'T27') then
            typma = 'TRIA3'
            nbnodes = 3
            nbnodes_post = 7
            node_init = 21
            elem_dim = 2
        else
            print *, typma
            ASSERT(ASTER_FALSE)
        end if
!
        ASSERT((nbnodes .ge. 2) .and. (nbnodes .le. 4))
        ASSERT((nbnodes_post .ge. 2) .and. (nbnodes_post .le. 9))
        ASSERT((elem_dim .eq. 1) .or. (elem_dim .eq. 2))
!
        call jevech('PGEOMER', 'L', jv_geom)
!
! - Get coordinates
!
        node_i = 0
        do inode = node_init, node_init-1+nbnodes_post
            node_i = node_i+1
            do idim = 1, elem_dim+1
                nodes_coor(idim, node_i) = zr(jv_geom+(inode-1)*(elem_dim+1)+idim-1)
            end do
        end do
!
        if (l_debug) then
            write (6, *) "hhoGeomFace debug"
            write (6, *) "typma: ", typma
            write (6, *) "nbnodes: ", nbnodes
            write (6, *) "nbnodes_post: ", nbnodes_post
            do inode = 1, nbnodes
                write (6, *) "node ", inode, ": ", nodes_coor(1:3, inode)
            end do
            write (6, *) "end hhoGeomFace debug"
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoFaceInit(hhoFace, typma, ndim, nbnodes, nbnodes_post, nodes_coor, &
                           num_nodes_loc, num_face_loc, barycenter_cell)
!
        implicit none
!
        integer(kind=8), intent(in)                     :: typma
        integer(kind=8), intent(in)                     :: ndim
        real(kind=8), dimension(3, 9), intent(in)       :: nodes_coor
        integer(kind=8), intent(in)                     :: nbnodes, nbnodes_post
        integer(kind=8), dimension(9), intent(in)       :: num_nodes_loc
        integer(kind=8), intent(in)                     :: num_face_loc
        real(kind=8), dimension(3), optional, intent(in):: barycenter_cell
        type(HHO_Face), intent(out)                     :: hhoFace
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tools
!
! Initialize a HHO Face
!
! --------------------------------------------------------------------------------------------------
!
! In typma              : type of element
! In ndim               : dimension of the element
! In nodes_coor         : coordinates of nodes
! In nbnodes            : number of nodes
! In barycenter_cell    : barycenter of the cell
! Out hhoFace           : a HHO Face
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: numsorted(9), ino
        real(kind=8) :: normal_curr(3)
        ASSERT(nbnodes .le. 4)
!
! --- Init: be carefull the order to call the functions is important
!
        hhoFace%typema = typma
        hhoFace%nbnodes = nbnodes
        hhoFace%nbnodes_post = nbnodes_post
        hhoFace%ndim = ndim
        hhoFace%coorno = hhoFaceInitCoor(nodes_coor, nbnodes, nbnodes_post, ndim, numsorted)
        hhoFace%barycenter = barycenter(hhoFace%coorno, hhoFace%nbnodes)
        if (present(barycenter_cell)) then
            hhoFace%normal = hhoNormalFace3(hhoFace, barycenter_cell)
        else
            hhoFace%normal = hhoNormalFace2(typma, nodes_coor)
        end if
        normal_curr = hhoNormalFace2(hhoFace%typema, hhoFace%coorno(1:3, 1:4))
        if (dot_product(hhoFace%normal, normal_curr) >= 0.d0) then
            hhoFace%reorientation = ASTER_FALSE
        else
            hhoFace%reorientation = ASTER_TRUE
        end if
        hhoFace%measure = hhoMeasureFace(hhoFace)
        hhoFace%diameter = hhoDiameterFace(hhoFace)
        hhoFace%axes = hhoLocalAxesFace(hhoFace)
        hhoFace%face_loc = num_face_loc
        do ino = 1, hhoFace%nbnodes_post
            hhoFace%nodes_loc(ino) = num_nodes_loc(numsorted(ino))
        end do
        ! Last node is the index of barycenter
        hhoFace%node_bar_loc = num_nodes_loc(hhoFace%nbnodes_post)
        hhoFace%l_jaco_cst = hhoIsJacobCst(hhoFace%typema, hhoFace%coorno, hhoFace%ndim+1)
!
        if (ndim == 1) then
            if (lteatt("TYPMOD", "AXIS")) then
                ! Special treatment is applied to these edges (quadrature, stabilization, gradient)
                hhoFace%l_axis = ASTER_TRUE
                hhoFace%l_axis_on_axe = (hhoFace%barycenter(1) < 1.d-12)
            end if
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCellInit(typma, nodes_coor, hhoCell)
!
        implicit none
!
        character(len=8), intent(in)                :: typma
        real(kind=8), dimension(3, 27), intent(in)  :: nodes_coor
        type(HHO_Cell), intent(out)                 :: hhoCell
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tools
!
! Initialize a HHO cell
!
! --------------------------------------------------------------------------------------------------
!
! In typma              : type of element
! In nodes_coor         : coordinates of nodes
! Out hhoCell           : a HHO cell
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
! --------------------------------------------------------------------------------------------------
!
!
        if (typma == 'HEXA27') then
            hhoCell%typema = MT_HEXA8
            hhoCell%nbnodes = 8
            hhoCell%ndim = 3
            hhoCell%nbfaces = 6
            hhoCell%node_bar_loc = 27
        else if (typma == 'TETRA15') then
            hhoCell%typema = MT_TETRA4
            hhoCell%nbnodes = 4
            hhoCell%ndim = 3
            hhoCell%nbfaces = 4
            hhoCell%node_bar_loc = 15
        else if (typma == 'PYRAM19') then
            hhoCell%typema = MT_PYRAM5
            hhoCell%nbnodes = 5
            hhoCell%ndim = 3
            hhoCell%nbfaces = 5
            hhoCell%node_bar_loc = 19
        else if (typma == 'PENTA21') then
            hhoCell%typema = MT_PENTA6
            hhoCell%nbnodes = 6
            hhoCell%ndim = 3
            hhoCell%nbfaces = 5
            hhoCell%node_bar_loc = 21
        else if (typma == 'QUAD9') then
            hhoCell%typema = MT_QUAD4
            hhoCell%nbnodes = 4
            hhoCell%ndim = 2
            hhoCell%nbfaces = 4
            hhoCell%node_bar_loc = 9
        else if (typma == 'TRIA7') then
            hhoCell%typema = MT_TRIA3
            hhoCell%nbnodes = 3
            hhoCell%ndim = 2
            hhoCell%nbfaces = 3
            hhoCell%node_bar_loc = 7
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ----- Copy coordinates
        hhoCell%coorno = nodes_coor
        hhoCell%barycenter = barycenter(hhoCell%coorno, hhoCell%nbnodes)
        hhoCell%measure = hhoMeasureCell(hhoCell)
        hhoCell%diameter = hhoDiameterCell(hhoCell)
!
        hhoCell%axes = hhoLocalAxesCell(hhoCell)
!
! ----- Jacobienne
!
        if (hhoCell%typema == MT_TETRA4 .or. hhoCell%typema == MT_TRIA3) then
            hhoCell%l_jaco_cst = ASTER_TRUE
        else
            hhoCell%l_jaco_cst = hhoIsJacobCst(hhoCell%typema, hhoCell%coorno, hhoCell%ndim)
        end if
!
        if (hhoCell%ndim == 2) then
            hhoCell%l_axis = lteatt("TYPMOD", "AXIS")
        end if
!
        if (l_debug) then
            call hhoCell%print()
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoDataInit(hhoData)
!
        implicit none
!
        type(HHO_Data), intent(out)            :: hhoData
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Get hho data on the modelisation
!
! --------------------------------------------------------------------------------------------------
!
! Out hhoData  :: information on the modelisation
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
        aster_logical :: l_adapt_coeff
        integer(kind=8) :: grad_deg, face_deg, cell_deg, jv_mater
        integer(kind=8) ::  jtab(1), iret, icodre(1)
        real(kind=8) :: coef_stab
        real(kind=8) :: resu_vale(1)
        character(len=16), parameter :: resu_name(1) = (/"COEF_STAB"/)
!
! --- Read Parameters
!
        ASSERT(lteatt('TYPMOD2', 'HHO') .or. lteatt('TYPMOD2', 'HHO_GRAD'))
        if (lteatt('FORMULATION', 'HHO_CSTE')) then
            face_deg = 0
            cell_deg = 0
            grad_deg = 0
        elseif (lteatt('FORMULATION', 'HHO_LINE') .or. lteatt('FORMULATION', 'HHO_VLINE')) then
            face_deg = 1
            cell_deg = 1
            grad_deg = 1
        elseif (lteatt('FORMULATION', 'HHO_QUAD') .or. lteatt('FORMULATION', 'HHO_VQUAD')) then
            face_deg = 2
            cell_deg = 2
            grad_deg = 2
        elseif (lteatt('FORMULATION', 'HHO_CUBI')) then
            face_deg = 3
            cell_deg = 3
            grad_deg = 3
        elseif (lteatt('FORMULATION', 'HHO_QUAR')) then
            face_deg = 4
            cell_deg = 4
            grad_deg = 4
        elseif (lteatt('FORMULATION', 'HHO_MCSTE')) then
            face_deg = 0
            cell_deg = 1
            grad_deg = 0
        elseif (lteatt('FORMULATION', 'HHO_MLINE')) then
            face_deg = 1
            cell_deg = 2
            grad_deg = 1
        elseif (lteatt('FORMULATION', 'HHO_MQUAD')) then
            face_deg = 2
            cell_deg = 3
            grad_deg = 2
        elseif (lteatt('FORMULATION', 'HHO_MCUBI')) then
            face_deg = 3
            cell_deg = 4
            grad_deg = 3
        elseif (lteatt('FORMULATION', 'HHO_MQUAR')) then
            face_deg = 4
            cell_deg = 5
            grad_deg = 4
        else
            ASSERT(ASTER_FALSE)
        end if
!
! --- Get coefficient (if possible)
!
        coef_stab = 0.d0
        iret = -1
        l_adapt_coeff = ASTER_TRUE
        call tecach('NNN', 'PMATERC', 'L', iret, nval=1, itab=jtab)
        if (iret .eq. 0) then
            call jevech('PMATERC', 'L', jv_mater)
            call rcvala(zi(jv_mater), ' ', 'HHO', 0, ' ', [0.0], &
                        1, resu_name, resu_vale, icodre, 0)
            l_adapt_coeff = icodre(1) == 1
            if (.not. l_adapt_coeff) then
                coef_stab = resu_vale(1)
            end if
        end if
!
! --- Init
!
        call hhoData%initialize(face_deg, cell_deg, grad_deg, coef_stab, l_debug, l_adapt_coeff)
!
        if (hhoData%debug()) then
            call hhoData%print()
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuad)
!
        implicit none
!
        type(HHO_Cell), intent(out)                 :: hhoCell
        type(HHO_Data), intent(out)                 :: hhoData
        integer(kind=8), intent(in), optional       :: npg
        type(HHO_Quadrature), optional, intent(out) :: hhoQuad

!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Get hho data on the modelisation
!
! --------------------------------------------------------------------------------------------------
!
!   Initialize information for a HHO cell
!   Out hhoCell  : the current HHO Cell
!   Out hhoData  : information on HHO methods
!   Out hhoQuad  : Quadrature for the cell
!   In npg (optional)      : number of quadrature point for the face
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: nbnodes, elem_dim
        character(len=8) :: typma
        real(kind=8) :: coor(3, 27)
!
        coor = 0.d0
!
! --- Get HHO informations
!
        call hhoGeomData(coor, nbnodes, typma, elem_dim)
        call hhoDataInit(hhoData)
!
! --- Initialize HHO Cell
        call hhoCellInit(typma, coor, hhoCell)
!
! --- Get quadrature (optional)
!
        if (present(hhoQuad)) then
            if (present(npg)) then
                call hhoQuad%initCell(hhoCell, npg)
            else
                call hhoQuad%getQuadCell(hhoCell, 2*hhoData%cell_degree(), &
                                         param=ASTER_TRUE)
            end if
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoInfoInitCellAndFace(hhoCell, hhoData, npg, hhoQuad)
!
        implicit none
!
        type(HHO_Cell), intent(out)                 :: hhoCell
        type(HHO_Data), intent(out)                 :: hhoData
        integer(kind=8), intent(in), optional       :: npg
        type(HHO_Quadrature), optional, intent(out) :: hhoQuad

!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Get hho data on the modelisation
!
! --------------------------------------------------------------------------------------------------
!
!   Initialize information for a HHO cell
!   Out hhoCell  : the current HHO Cell
!   Out hhoData  : information on HHO methods
!   Out hhoQuad  : Quadrature for the cell
!   In npg (optional)      : number of quadrature point for the face
! --------------------------------------------------------------------------------------------------
!
        if (present(npg)) then
            call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuad)
        else
            call hhoInfoInitCell(hhoCell, hhoData)
        end if
        call hhoInitFacesOfCell(hhoCell)
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoInfoInitFace(hhoFace, hhoData, npg, hhoQuadFace)
!
        implicit none
!
        type(HHO_Face), intent(out)   :: hhoFace
        type(HHO_Data), intent(out)   :: hhoData
        integer(kind=8), intent(in), optional :: npg
        type(HHO_Quadrature), intent(out), optional :: hhoQuadFace
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tool
!
! Get hho info on the modelisation and initialize
!
! --------------------------------------------------------------------------------------------------
!
!   Initialize information for surfacique load
!   Out hhoFace  : the current HHO Face
!   Out hhoData  : information on HHO methods
!   Out hhoQuad  : Quadrature for the face
!   In npg (optional)      : number of quadrature point for the face
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: nbnodes, elem_dim, numnodes(9), enumf, nbnodes_post
        real(kind=8) :: nodes_coor(3, 9)
        character(len=8) :: typma
!
! --- Get HHO informations
!
        call hhoGeomFace(nodes_coor, nbnodes, nbnodes_post, typma, elem_dim)
        ASSERT(elem_dim == 1 .or. elem_dim == 2)
        call hhoDataInit(hhoData)
!
        if (typma == 'QUAD4') then
            numnodes(1:9) = (/1, 2, 3, 4, 5, 6, 7, 8, 9/)
            enumf = MT_QUAD4
        else if (typma == 'TRIA3') then
            numnodes(1:7) = (/1, 2, 3, 4, 5, 6, 7/)
            enumf = MT_TRIA3
        else if (typma == 'SEG2 ') then
            numnodes(1:3) = (/1, 2, 3/)
            enumf = MT_SEG2
        else
            ASSERT(ASTER_FALSE)
        end if
!
! --- Initialize HHO Face
!
        call hhoFaceInit(hhoFace, enumf, elem_dim, nbnodes, nbnodes_post, &
                         nodes_coor, numnodes, 1)

! --- Get quadrature (optional)
!
        if (present(hhoQuadFace)) then
            if (present(npg)) then
                call hhoQuadFace%initFace(hhoFace, npg)
            else
                call hhoQuadFace%GetQuadFace(hhoFace, 2*hhoData%face_degree(), &
                                             param=ASTER_TRUE)
            end if
        end if
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoInitFacesOfCell(hhoCell)
!
        implicit none
!
        type(HHO_Cell), intent(inout)                 :: hhoCell
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tools
!
! Initialize a HHO cell - only face
!
! --------------------------------------------------------------------------------------------------
!
! Out hhoCell           : a HHO cell
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), parameter                      :: max_faces = 6
        integer(kind=8), parameter                      :: max_nodes = 9
        integer(kind=8), dimension(max_nodes, max_faces) :: nodes_faces
        integer(kind=8), dimension(max_faces)           :: nbnodes_faces
        integer(kind=8), dimension(max_faces)           :: nbnodes_post_faces
        integer(kind=8), dimension(max_faces)           :: type_faces
        real(kind=8), dimension(3, max_nodes) :: coor_face
        integer(kind=8) :: i_face, i_node

!
        if (.not. hhoCell%l_face_init) then
            ! --- Init
            nodes_faces = 0
            nbnodes_faces = 0
            nbnodes_post_faces = 0
            type_faces = 0
            coor_face = 0.d0

            if (hhoCell%typema == MT_HEXA8) then
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N21
                nodes_faces(1:9, 1) = (/1, 4, 3, 2, 12, 11, 10, 9, 21/)
                nbnodes_faces(1) = 4
                nbnodes_post_faces(1) = 9
                type_faces(1) = MT_QUAD4
! ----- Face 2 -> N22
                nodes_faces(1:9, 2) = (/1, 2, 6, 5, 9, 14, 17, 13, 22/)
                nbnodes_faces(2) = 4
                nbnodes_post_faces(2) = 9
                type_faces(2) = MT_QUAD4
! ----- Face 3 -> N23
                nodes_faces(1:9, 3) = (/2, 3, 7, 6, 10, 15, 18, 14, 23/)
                nbnodes_faces(3) = 4
                nbnodes_post_faces(3) = 9
                type_faces(3) = MT_QUAD4
! ----- Face 4 -> N24
                nodes_faces(1:9, 4) = (/3, 4, 8, 7, 11, 16, 19, 15, 24/)
                nbnodes_faces(4) = 4
                nbnodes_post_faces(4) = 9
                type_faces(4) = MT_QUAD4
! ----- Face 5 -> N25
                nodes_faces(1:9, 5) = (/1, 5, 8, 4, 13, 20, 16, 12, 25/)
                nbnodes_faces(5) = 4
                nbnodes_post_faces(5) = 9
                type_faces(5) = MT_QUAD4
! ----- Face 6 -> N26
                nodes_faces(1:9, 6) = (/5, 6, 7, 8, 17, 18, 19, 20, 26/)
                nbnodes_faces(6) = 4
                nbnodes_post_faces(6) = 9
                type_faces(6) = MT_QUAD4
!
            else if (hhoCell%typema == MT_TETRA4) then
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N11
                nodes_faces(1:7, 1) = (/1, 3, 2, 7, 6, 5, 11/)
                nbnodes_faces(1) = 3
                nbnodes_post_faces(1) = 7
                type_faces(1) = MT_TRIA3
! ----- Face 2 -> N12
                nodes_faces(1:7, 2) = (/1, 2, 4, 5, 9, 8, 12/)
                nbnodes_faces(2) = 3
                nbnodes_post_faces(2) = 7
                type_faces(2) = MT_TRIA3
! ----- Face 3 -> N13
                nodes_faces(1:7, 3) = (/1, 4, 3, 8, 10, 7, 13/)
                nbnodes_faces(3) = 3
                nbnodes_post_faces(3) = 7
                type_faces(3) = MT_TRIA3
! ----- Face 4 -> N14
                nodes_faces(1:7, 4) = (/2, 3, 4, 6, 10, 9, 14/)
                nbnodes_faces(4) = 3
                nbnodes_post_faces(4) = 7
                type_faces(4) = MT_TRIA3
!
            else if (hhoCell%typema == MT_PYRAM5) then
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N14
                nodes_faces(1:9, 1) = (/1, 4, 3, 2, 9, 8, 7, 6, 14/)
                nbnodes_faces(1) = 4
                nbnodes_post_faces(1) = 9
                type_faces(1) = MT_QUAD4
! ----- Face 2 -> N15
                nodes_faces(1:7, 2) = (/1, 2, 5, 6, 11, 10, 15/)
                nbnodes_faces(2) = 3
                nbnodes_post_faces(2) = 7
                type_faces(2) = MT_TRIA3
! ----- Face 3 -> N16
                nodes_faces(1:7, 3) = (/2, 3, 5, 7, 12, 11, 16/)
                nbnodes_faces(3) = 3
                nbnodes_post_faces(3) = 7
                type_faces(3) = MT_TRIA3
! ----- Face 4 -> N17
                nodes_faces(1:7, 4) = (/3, 4, 5, 8, 13, 12, 17/)
                nbnodes_faces(4) = 3
                nbnodes_post_faces(4) = 7
                type_faces(4) = MT_TRIA3
! ----- Face 5 -> N18
                nodes_faces(1:7, 5) = (/4, 1, 5, 9, 10, 13, 18/)
                nbnodes_faces(5) = 3
                nbnodes_post_faces(5) = 7
                type_faces(5) = MT_TRIA3
!
            else if (hhoCell%typema == MT_PENTA6) then
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N16
                nodes_faces(1:9, 1) = (/1, 2, 5, 4, 7, 11, 13, 10, 16/)
                nbnodes_faces(1) = 4
                nbnodes_post_faces(1) = 9
                type_faces(1) = MT_QUAD4
! ----- Face 2 -> N17
                nodes_faces(1:9, 2) = (/2, 3, 6, 5, 8, 12, 14, 11, 17/)
                nbnodes_faces(2) = 4
                nbnodes_post_faces(2) = 9
                type_faces(2) = MT_QUAD4
! ----- Face 3 -> N18
                nodes_faces(1:9, 3) = (/1, 4, 6, 3, 10, 15, 12, 9, 18/)
                nbnodes_faces(3) = 4
                nbnodes_post_faces(3) = 9
                type_faces(3) = MT_QUAD4
! ----- Face 4 -> N19
                nodes_faces(1:7, 4) = (/1, 3, 2, 9, 8, 7, 19/)
                nbnodes_faces(4) = 3
                nbnodes_post_faces(4) = 7
                type_faces(4) = MT_TRIA3
! ----- Face 5 -> N20
                nodes_faces(1:7, 5) = (/4, 5, 6, 15, 14, 13, 20/)
                nbnodes_faces(5) = 3
                nbnodes_post_faces(5) = 7
                type_faces(5) = MT_TRIA3
!
            else if (hhoCell%typema == MT_QUAD4) then
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N5
                nodes_faces(1:3, 1) = (/1, 2, 5/)
                nbnodes_faces(1) = 2
                nbnodes_post_faces(1) = 3
                type_faces(1) = MT_SEG2
! ----- Face 2 -> N6
                nodes_faces(1:3, 2) = (/2, 3, 6/)
                nbnodes_faces(2) = 2
                nbnodes_post_faces(2) = 3
                type_faces(2) = MT_SEG2
! ----- Face 3 -> N7
                nodes_faces(1:3, 3) = (/3, 4, 7/)
                nbnodes_faces(3) = 2
                nbnodes_post_faces(3) = 3
                type_faces(3) = MT_SEG2
! ----- Face 4 -> N8
                nodes_faces(1:3, 4) = (/4, 1, 8/)
                nbnodes_faces(4) = 2
                nbnodes_post_faces(4) = 3
                type_faces(4) = MT_SEG2
!
            else if (hhoCell%typema == MT_TRIA3) then
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N4
                nodes_faces(1:3, 1) = (/1, 2, 4/)
                nbnodes_faces(1) = 2
                nbnodes_post_faces(1) = 3
                type_faces(1) = MT_SEG2
! ----- Face 2 -> N5
                nodes_faces(1:3, 2) = (/2, 3, 5/)
                nbnodes_faces(2) = 2
                nbnodes_post_faces(2) = 3
                type_faces(2) = MT_SEG2
! ----- Face 3 -> N6
                nodes_faces(1:3, 3) = (/3, 1, 6/)
                nbnodes_faces(3) = 2
                nbnodes_post_faces(3) = 3
                type_faces(3) = MT_SEG2
!
            else
                ASSERT(ASTER_FALSE)
            end if

            do i_face = 1, hhoCell%nbfaces
                coor_face = 0.d0
                do i_node = 1, nbnodes_post_faces(i_face)
                    coor_face(1:3, i_node) = hhoCell%coorno(1:3, nodes_faces(i_node, i_face))
                end do
                call hhoFaceInit(hhoCell%faces(i_face), type_faces(i_face), hhoCell%ndim-1, &
                                 nbnodes_faces(i_face), nbnodes_post_faces(i_face), &
                                 coor_face, nodes_faces(:, i_face), &
                                 i_face, hhoCell%barycenter)
            end do
            hhoCell%l_face_init = ASTER_TRUE
        end if
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoUpdateCoordinatesFace(disp_nodes, hhoFace)
!
        implicit none
!
        real(kind=8), intent(in) :: disp_nodes(3, 9)
        type(HHO_Face), intent(inout) :: hhoFace
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tools
!
! Update coordinate of a face and recompute quantities
!
! --------------------------------------------------------------------------------------------------
!
! Out hhoFace           : a HHO face
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ino
        real(kind=8), dimension(3, 9) :: nodes_coor, normal_curr(3)

        do ino = 1, hhoFace%nbnodes
            hhoFace%coorno(1:3, ino) = hhoFace%coorno(1:3, ino)+ &
                                       disp_nodes(1:3, ino)
            nodes_coor(1:3, hhoFace%nodes_loc(ino)) = hhoFace%coorno(1:3, ino)
        end do
!
        hhoFace%barycenter = barycenter(hhoFace%coorno, hhoFace%nbnodes)
        hhoFace%normal = hhoNormalFace2(hhoFace%typema, nodes_coor)
        normal_curr = hhoNormalFace2(hhoFace%typema, hhoFace%coorno(1:3, 1:4))
        if (dot_product(hhoFace%normal, normal_curr) >= 0.d0) then
            hhoFace%reorientation = ASTER_FALSE
        else
            hhoFace%reorientation = ASTER_TRUE
        end if
        hhoFace%measure = hhoMeasureFace(hhoFace)
        hhoFace%diameter = hhoDiameterFace(hhoFace)
        hhoFace%axes = hhoLocalAxesFace(hhoFace)
        hhoFace%l_jaco_cst = hhoIsJacobCst(hhoFace%typema, hhoFace%coorno, hhoFace%ndim+1)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
end module
