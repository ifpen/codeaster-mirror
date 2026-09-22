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
! aslint: disable=W1306
!
! ==================================================================================================
!
! Module for pairing meshes
!
! ==================================================================================================
!
module MeshPairing_module
! ==================================================================================================
    use mesh_type
    use mesh_cell_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: fastPair, robustPair, pairDeallocate, pairAllocate, pairResize
    public :: inteCellArea, intePoinCoor, quadPoinCoor
    public :: getPairJV, getInteJV, pairAdd
    private :: pairGetStartCells, getClosestNodesFromCell, getCellsFromNode
    private :: cellInteProj, cellProjOnCell, nodeProjOnCell
    private :: poinProjOnCell, poinProjByVect, poinIsInsideOtherCell
    private :: inteCellChck, addPoinInte2D, addPoinInte3D, addPoinOnEdge
    private :: inteCellSegm, intePoinSort, intePoinInCell
    private :: isFatalError
    private :: getCritParamArea
! ==================================================================================================
    private
#include "asterc/r8gaem.h"
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/mesh_pairing_type.h"
#include "asterfort/ordr8.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
! ==================================================================================================
! ==================================================================================================
! Global variables
! ==================================================================================================
! - Index of next nodes for segment, triangle and quadrangle
    integer(kind=8), parameter :: nodeNextSEG(2) = (/2, 1/)
    integer(kind=8), parameter :: nodeNextTRIA(3) = (/2, 3, 1/)
    integer(kind=8), parameter :: nodeNextQUAD(4) = (/2, 3, 4, 1/)
! - Index of previous nodes for triangle and quadrangle
    integer(kind=8), parameter :: nodePrevTRIA(3) = (/3, 1, 2/)
    integer(kind=8), parameter :: nodePrevQUAD(4) = (/4, 1, 2, 3/)
! ==================================================================================================
! Type for pairing
! ==================================================================================================
    type MESH_PAIRING
! ----- Dimension of space (2 or 3)
        integer(kind=8) :: spaceDime = 0
! ----- Flag for debug (text)
        aster_logical :: debug = ASTER_FALSE
! ----- Main tolerance for pairing
        real(kind=8) :: pairTole = 0.d0
! ----- Tolerance to remove some intersection cells
        real(kind=8) :: areaTole = 0.d0
! ----- Tolerance for extension of cell
        real(kind=8) :: distRatio = 0.d0
! ----- Number of pairs
        integer(kind=8) :: nbPair = 0
        integer(kind=8) :: nbPairMax = 0
! ----- Pointer to pairs
        integer(kind=8), pointer :: pair(:) => null()
! ----- Pointer to list of intersection points
        integer(kind=8), pointer :: nbPoinInte(:) => null()
        real(kind=8), pointer :: poinInte(:) => null()
    end type MESH_PAIRING
! ==================================================================================================
! Type for parameters of projection algorithm
! ==================================================================================================
    type MESH_PROJ_PARA
! ----- Parameters for algorithm
        aster_logical :: debug = ASTER_FALSE
        integer(kind=8) :: newtIterMaxi = 0
        real(kind=8) :: newtTole = 0.d0
! ----- Coordinate of point to project (in global frame)
        real(kind=8) :: coorPoinGlob(3) = 0.d0
! ----- Dimension of space (2 or 3)
        integer(kind=8) :: spaceDime = 0
! ----- Vector to project
        real(kind=8):: projVect(3) = 0.d0
! ----- Coordinate of point after projection (in cell parametric space)
        real(kind=8) :: coorProjPara(2) = 0.d0
! ----- Tangents at projection of point
        real(kind=8) :: tau1(3) = 0.d0, tau2(3) = 0.d0
! ----- Return code
        integer(kind=8):: errorCode = 0
    end type MESH_PROJ_PARA
! ==================================================================================================
    public :: MESH_PAIRING, nodeNextTRIA, nodeNextQUAD, nodePrevTRIA, nodePrevQUAD
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! fastPair
!
! Fast pairing on zone
!
! In  mesh             : mesh
! In  newgeo           : updated coordinates of nodes
! In  mastConxInvName  : name of object for inverse connectivity of master cells on current zone
! In  mastNeighName    : name of object for neighbours of master cells
! In  slavNeighName    : name of object for neighbours of slave cells
! In  nbCellSlav       : number of slave cells
! In  nbCellMast       : number of master cells
! In  nbNodeMast       : number of master nodes
! In  listCellSlav     : list of slave cells
! In  listCellMast     : list of master cells
! In  listNodeMast     : list of master nodes
! In  meshPairing      : main datastructure for pairing
!
! --------------------------------------------------------------------------------------------------
    subroutine fastPair(mesh, newgeo, mastConxInvName, &
                        mastNeighName, slavNeighName, &
                        nbCellSlav, nbCellMast, nbNodeMast, &
                        listCellSlav, listCellMast, listNodeMast, &
                        meshPairing)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: mesh
        character(len=24), intent(in) :: newgeo, mastConxInvName
        character(len=24), intent(in) :: mastNeighName, slavNeighName
        integer(kind=8), intent(in) :: nbCellSlav, nbCellMast
        integer(kind=8), intent(in) :: listCellMast(nbCellMast), listCellSlav(nbCellSlav)
        integer(kind=8), intent(in) :: nbNodeMast
        integer(kind=8), intent(in) :: listNodeMast(nbNodeMast)
        type(MESH_PAIRING), intent(inout) :: meshPairing
! ----- Local
        aster_logical :: pair_exist, isFatal, l_recup
        integer(kind=8) :: nbSlavStart, nbMastStart
        integer(kind=8), pointer :: cellSlavStart(:) => null()
        integer(kind=8), pointer :: cellMastStart(:) => null()
        integer(kind=8), pointer :: cellMastFlag(:) => null()
        integer(kind=8), pointer :: cellSlavFlag(:) => null(), cellMastPaired(:) => null()
        integer(kind=8) :: slavIndxMini, slavIndxMaxi
        integer(kind=8) :: mastIndxMini, mastIndxMaxi
        integer(kind=8) :: iCell, iMastNeigh, iSlavNeigh
        real(kind=8), pointer :: nodeCoor(:) => null()
        integer(kind=8), pointer :: meshTypeGeom(:) => null()
        integer(kind=8), pointer :: meshConx(:) => null(), meshConxCumu(:) => null()
        integer(kind=8), pointer :: mastConxInv(:) => null(), mastConxInvCumu(:) => null()
        type(CELL_GEOM) :: cellSlav, cellMast
        type(CELL_GEOM) :: cellSlavLine, cellMastLine
        integer(kind=8) :: nbSlavNeigh, nbMastPaired, inteNeigh(MAX_NB_NEIGH), nbMastNeigh
        integer(kind=8) :: cellNeigh(MAX_NB_NEIGH), iret
        real(kind=8) :: inteArea
        integer(kind=8) :: mastFindIndx, cellMastNume, cellMastIndx, cellSlavIndx, cellSlavNume
        integer(kind=8) :: nbPoinInte
        real(kind=8) :: poinInte(2, MAX_NB_INTE), poinInteReal(3, MAX_NB_INTE)
        integer(kind=8) :: cellNeighNume, cellNeighIndx
        integer(kind=8), pointer :: meshMastNeigh(:) => null(), meshSlavNeigh(:) => null()
        real(kind=8) :: critCellArea
!   ------------------------------------------------------------------------------------------------
!
        pair_exist = ASTER_TRUE
        inteNeigh = 0
        cellNeigh = 0
        mastIndxMini = minval(listCellMast)
        mastIndxMaxi = maxval(listCellMast)
        slavIndxMaxi = maxval(listCellSlav)
        slavIndxMini = minval(listCellSlav)

        critCellArea = 0d0
! ----- Management of debug
        if (meshPairing%debug) then
            write (6, *) "================"
            write (6, *) "= Fast pairing ="
            write (6, *) "================"
            write (6, *) " "
        end if

! ----- Access to updated geometry
        call jeveuo(newgeo(1:19)//'.VALE', 'L', vr=nodeCoor)

! ----- Access to mesh
        call jeveuo(mesh(1:8)//'.TYPMAIL', 'L', vi=meshTypeGeom)
        call jeveuo(mesh(1:8)//'.CONNEX', 'L', vi=meshConx)
        call jeveuo(jexatr(mesh(1:8)//'.CONNEX', 'LONCUM'), 'L', vi=meshConxCumu)
        call jeveuo(mastConxInvName, 'L', vi=mastConxInv)
        call jeveuo(jexatr(mastConxInvName, 'LONCUM'), 'L', vi=mastConxInvCumu)
        call jeveuo(mastNeighName, 'L', vi=meshMastNeigh)
        call jeveuo(slavNeighName, 'L', vi=meshSlavNeigh)

! ----- Protection
        if (nbCellSlav .eq. 0 .or. nbCellMast .eq. 0) then
            call utmess('F', 'MESH4_1')
        end if

! ----- List of starting cells
        AS_ALLOCATE(vi=cellSlavStart, size=nbCellSlav)
        AS_ALLOCATE(vi=cellMastStart, size=nbCellSlav)

! ----- Flag for slave cell usage status
!                        0 - Never used
!                        1 - Used as starting point
!                        2 - Used
        AS_ALLOCATE(vi=cellSlavFlag, size=slavIndxMaxi+1-slavIndxMini)

! ----- Flag for master cell usage status
!                        0 - Never used
!                        1 - Used as starting point
        AS_ALLOCATE(vi=cellMastFlag, size=mastIndxMaxi+1-mastIndxMini)

! ----- Master elements paired with the current slave cell (number of master elements: nbMastPaired)
        AS_ALLOCATE(vi=cellMastPaired, size=nbCellMast)

! ----- While loop on the existence of a pair slave-master
        do while (pair_exist)
            if (meshPairing%debug) then
                WRITE (6, *) "Get cells for starting search"
                WRITE (6, *) "============================="
            end if
! --------- Get initial start cells
            call pairGetStartCells(meshPairing, nodeCoor, &
                                   meshTypeGeom, meshConx, meshConxCumu, &
                                   mastConxInv, mastConxInvCumu, &
                                   nbCellSlav, listCellSlav, cellSlavFlag, &
                                   nbCellMast, listCellMast, &
                                   nbNodeMast, listNodeMast, &
                                   nbMastStart, cellMastStart, &
                                   nbSlavStart, cellSlavStart)
            ASSERT(nbMastStart .le. 1)
            ASSERT(nbSlavStart .le. 1)
            if (meshPairing%debug) then
                if (nbSlavStart .eq. 1 .and. nbMastStart .eq. 1) then
                    WRITE (6, *) "  New starting point (M/S): ", cellMastStart(1), cellSlavStart(1)
                else
                    WRITE (6, *) "  No new starting point (M/S)"
                end if
            end if

! --------- No more slave cell
            if (nbSlavStart == 0) then
                pair_exist = ASTER_FALSE
            end if

            do while (nbSlavStart > 0)
! ------------- Get slave element
                cellSlavNume = cellSlavStart(1)
                cellSlavIndx = cellSlavNume+1-slavIndxMini

! ------------- Shift list of starting slave cells
                do iCell = 1, nbSlavStart-1
                    cellSlavStart(iCell) = cellSlavStart(iCell+1)
                end do
                nbSlavStart = nbSlavStart-1

! ------------- Create slave cell
                call cellCreate(cellSlavNume, nodeCoor, &
                                meshTypeGeom, meshConx, meshConxCumu, &
                                cellSlav, cellSlavLine)

                if (meshPairing%debug) then
                    write (6, *) "Current slave element      : ", cellSlavNume
                    write (6, *) " Coordinates (global frame): ", &
                        cellSlav%coorNodeGlob(1:meshPairing%spaceDime, 1:cellSlav%nbNode)
                end if

! ------------- Get number of neighbours
                nbSlavNeigh = cellSlav%nbNeigh
                ASSERT(nbSlavNeigh .le. MAX_NB_NEIGH)
                if (meshPairing%debug) then
                    write (6, *) "Potential number of neighbours: ", nbSlavNeigh
                    do iSlavNeigh = 1, nbSlavNeigh
                        cellNeighIndx = MAX_NB_NEIGH*(cellSlavIndx-1)+iSlavNeigh
                        if (cellNeighIndx .ne. 0) then
                            cellNeighNume = meshSlavNeigh(cellNeighIndx)
                            if (cellNeighNume .eq. 0) then
                                write (6, *) "Neighbour (", iSlavNeigh, "): "
                            else

                                write (6, *) "Neighbour (", iSlavNeigh, "): ", cellNeighNume
                            end if
                        end if
                    end do
                end if
                cellNeigh = 0

! ------------- Get master element to start (don't use this starting cell for pairing)
                cellMastNume = cellMastStart(1)
                mastFindIndx = cellMastNume+1-mastIndxMini
                cellMastFlag(mastFindIndx) = 1
                if (meshPairing%debug) then
                    WRITE (6, *) "Avant décalage"
                    WRITE (6, *) " => ", nbMastStart, cellMastStart(1:nbMastStart)
                end if

! ------------- Delete cell in the list of master element start
                do iCell = 1, nbMastStart-1
                    cellMastStart(iCell) = cellMastStart(iCell+1)
                end do
                nbMastStart = nbMastStart-1
                if (meshPairing%debug) then
                    WRITE (6, *) "Après décalage"
                    WRITE (6, *) " => ", nbMastStart, cellMastStart(1:nbMastStart)
                end if

! ------------- Management of list of master elements: first element to seek
                cellMastPaired(1) = cellMastNume
                nbMastPaired = 1
                if (meshPairing%debug) then
                    write (6, *) "Master cell to start: ", cellMastNume
                end if

! ------------- Initialization list of pairs
                l_recup = ASTER_TRUE

! ------------- Loop on master elements => Look for the master elements
                do while (nbMastPaired > 0)
                    inteArea = 0.d0
                    nbPoinInte = 0
                    poinInte = 0.d0
                    poinInteReal = 0.d0

! ----------------- Get master element
                    cellMastNume = cellMastPaired(1)
                    cellMastIndx = cellMastNume+1-mastIndxMini

! ----------------- Shift list of master element to pair
                    do iCell = 1, nbMastPaired-1
                        cellMastPaired(iCell) = cellMastPaired(iCell+1)
                    end do
                    nbMastPaired = nbMastPaired-1

! ----------------- Create master cell
                    call cellCreate(cellMastNume, nodeCoor, &
                                    meshTypeGeom, meshConx, meshConxCumu, &
                                    cellMast, cellMastLine)
                    if (meshPairing%debug) then
                        write (6, *) "Current master element: ", cellMastNume
                        write (6, *) " Coordinates (global frame): ", &
                            cellMast%coorNodeGlob(1:meshPairing%spaceDime, 1:cellMast%nbNode)
                    end if

! ----------------- Get number of neighbours
                    nbMastNeigh = cellMast%nbNeigh
                    ASSERT(nbMastNeigh .le. MAX_NB_NEIGH)
                    if (meshPairing%debug) then
                        write (6, *) "Potential number of neighbours: ", nbMastNeigh
                    end if

! ----------------- Compute intersection of the two cells
                    if (meshPairing%debug) then
                        WRITE (6, *) "Compute intersection and projection in master space"
                    end if
                    call cellInteProj(meshPairing, &
                                      cellSlav, cellSlavLine, cellMastLine, &
                                      iret, &
                                      nbPoinInte, poinInte, &
                                      inteArea, inteNeigh)
                    isFatal = isFatalError(iret)
                    if (.not. isFatal .and. iret .ne. ERR_PAIR_NONE) then
                        call utmess('A', 'MESH4_3')
                        inteArea = 0.d0
                        nbPoinInte = 0
                    end if
                    ASSERT(nbPoinInte .le. MAX_NB_INTE)
                    if (meshPairing%debug) then
                        WRITE (6, *) "Intersection area: ", inteArea
                    end if
                    call getCritParamArea(cellSlavLine, critCellArea)
! ----------------- Add pair
                    if (inteArea > meshPairing%areaTole*critCellArea .and. &
                        iret == ERR_PAIR_NONE) then

! --------------------- Debug
                        if (meshPairing%debug) then
! ------------------------- Compute coordinates in global space for intersection points
                            call intePoinCoor(cellSlav, nbPoinInte, poinInte, poinInteReal)
                            WRITE (6, *) "Add pair: ", meshPairing%nbPair+1, &
                                "(", cellSlavNume, "-", cellMastNume, ")"
                            WRITE (6, *) "Nb points integrations                : ", &
                                nbPoinInte
                            WRITE (6, *) "Coor. points integrations (parametric): ", &
                                poinInte(:, 1:nbPoinInte)
                            WRITE (6, *) "Coef. points integrations (global)    : ", &
                                poinInteReal(:, 1:nbPoinInte)
                            WRITE (6, *) "Area of intersection                  : ", &
                                inteArea
                        end if
! --------------------- Add pair
                        call pairAdd(cellSlavNume, cellMastNume, &
                                     nbPoinInte, poinInte, &
                                     meshPairing)
                    end if

! ----------------- Find neighbour of current master element
                    if (inteArea > meshPairing%pairTole .or. l_recup) then
! --------------------- Prepare next master element
                        if (meshPairing%debug) then
                            WRITE (6, *) "Prepare next master element"
                        end if
                        do iMastNeigh = 1, nbMastNeigh
                            cellNeighIndx = MAX_NB_NEIGH*(cellMastIndx-1)+iMastNeigh
                            cellNeighNume = meshMastNeigh(cellNeighIndx)
                            if (cellNeighNume .ne. 0) then
                                if (cellMastFlag(cellNeighNume+1-mastIndxMini) .ne. 1) then
                                    nbMastPaired = nbMastPaired+1
                                    cellMastPaired(nbMastPaired) = cellNeighNume
                                    cellMastFlag(cellNeighNume+1-mastIndxMini) = 1
                                    if (meshPairing%debug) then
                                        WRITE (6, *) " => added: ", nbMastPaired, &
                                            cellMastPaired(nbMastPaired)
                                    end if
                                end if
                            end if
                        end do

! --------------------- Prepare next slave element
                        if (meshPairing%debug) then
                            WRITE (6, *) "Prepare next slave element"
                        end if
                        do iSlavNeigh = 1, nbSlavNeigh
                            cellNeighIndx = MAX_NB_NEIGH*(cellSlavIndx-1)+iSlavNeigh
                            cellNeighNume = meshSlavNeigh(cellNeighIndx)
                            if (meshPairing%debug) then
                                WRITE (6, *) " < Index:", iSlavNeigh
                                WRITE (6, *) " < cellNeighNume:", cellNeighNume
                                WRITE (6, *) " < inteNeigh: ", inteNeigh(iSlavNeigh)
                            end if
                            if (cellNeighNume .ne. 0) then
                                if (cellSlavFlag(cellNeighNume+1-slavIndxMini) .ne. 1 .and. &
                                    inteNeigh(iSlavNeigh) == 1) then
                                    cellNeigh(iSlavNeigh) = cellMastNume
                                    if (meshPairing%debug) then
                                        WRITE (6, *) " => added: ", cellMastNume
                                    end if
                                end if
                            end if
                        end do
                        l_recup = ASTER_FALSE
                    end if
                end do

! ------------- Next elements
                if (meshPairing%debug) then
                    WRITE (6, *) "Prepare next elements - Nb: ", nbSlavNeigh
                end if
                do iSlavNeigh = 1, nbSlavNeigh
                    cellNeighIndx = MAX_NB_NEIGH*(cellSlavIndx-1)+iSlavNeigh
                    cellNeighNume = meshSlavNeigh(cellNeighIndx)
                    if (cellNeighNume .ne. 0) then
                        if (cellNeigh(iSlavNeigh) .ne. 0 .and. &
                            cellSlavFlag(cellNeighNume+1-slavIndxMini) .ne. 1) then
                            nbSlavStart = nbSlavStart+1
                            cellSlavStart(nbSlavStart) = cellNeighNume
                            cellSlavFlag(cellNeighNume+1-slavIndxMini) = 1
                            nbMastStart = nbMastStart+1
                            cellMastStart(nbMastStart) = cellNeigh(iSlavNeigh)
                        end if
                    end if
                end do

! ------------- All master cells are candidates for pairing
                cellMastFlag(1:mastIndxMaxi+1-mastIndxMini) = 0
            end do
            !pair_exist = ASTER_FALSE
        end do

! ----- Clean memory
        AS_DEALLOCATE(vi=cellSlavStart)
        AS_DEALLOCATE(vi=cellMastStart)
        AS_DEALLOCATE(vi=cellSlavFlag)
        AS_DEALLOCATE(vi=cellMastFlag)
        AS_DEALLOCATE(vi=cellMastPaired)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! nodeProjOnCell
!
! Project nodes of cell (origin)  on another cell (target) in target parametric space using
! raytracing
!
! In  meshPairing      : main datastructure for pairing
! In  cellOrig         : geometric properties of cell to project
! In  cellOrigLine     : geometric properties of cell to project (linearized)
! In  cellTargLine     : geometric properties of cell where to project (linearized)
! IO  cellProj         : geometric properties of projected cell
! Out nbNodeProj       : number of projected nodes
! Out iret             : return code error
!
! --------------------------------------------------------------------------------------------------
    subroutine nodeProjOnCell(meshPairing, &
                              cellOrig, cellOrigLine, cellTargLine, &
                              cellProj, nbNodeProj, iret)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        type(CELL_GEOM), intent(in) :: cellOrig, cellOrigLine, cellTargLine
        type(CELL_GEOM), intent(inout) :: cellProj
        integer(kind=8), intent(out) :: nbNodeProj, iret
! ----- Local
        integer(kind=8) :: iNode, nbNode
        real(kind=8) :: coorNodeGlob(3), coorNodePara(3), projNodePara(3), projNodeGlob(3)
        type(CELL_SKIN_BASE) :: base
        aster_logical :: errorProj
!   ------------------------------------------------------------------------------------------------
!
        iret = ERR_PAIR_NONE
        nbNodeProj = 0
        ASSERT(cellOrigLine%isLinear)
        ASSERT(cellTargLine%isLinear)

        if (meshPairing%debug) then
            WRITE (6, *) "        Project nodes of cell on another cell using raytracing"
            WRITE (6, *) "        ======================================================"
        end if

! ----- Projected cell: same type as original one
        call cellCopyType(cellOrigLine, cellProj)
        cellProj%coorNodePara = 0.d0
        cellProj%coorNodeGlob = 0.d0

! ----- Project only vertices of original cell (linearized cell)
        nbNode = cellOrigLine%nbNode
        nbNodeProj = nbNode
        if (meshPairing%debug) then
            WRITE (6, *) "        Number of nodes on original cell: ", nbNode
        end if
        do iNode = 1, nbNode
            if (meshPairing%debug) then
                WRITE (6, *) "        Current node: ", iNode
            end if

! --------- Get coordinates of current node (real space and parametric space)
            coorNodeGlob = cellOrigLine%coorNodeGlob(:, iNode)
            coorNodePara = cellOrigLine%coorNodePara(:, iNode)
            if (meshPairing%debug) then
                WRITE (6, *) "        Coordinate of node to project (local frame): ", &
                    coorNodePara
                WRITE (6, *) "        Coordinate of node to project (global frame): ", &
                    coorNodeGlob
            end if

! --------- Compute base on this node (outward normal)
            call cellCompBaseAtPoint(cellOrig, meshPairing%spaceDime, coorNodePara, base)
            if (meshPairing%debug) then
                WRITE (6, *) "        Normal at node to project: ", base%norm
            end if

! --------- Project node on target cell
            call poinProjOnCell(meshPairing%spaceDime, meshPairing%pairTole, &
                                coorNodeGlob, base%norm, &
                                cellTargLine, &
                                projNodePara, errorProj)

! --------- Careful ! projNodePara can be outside of reference frame !

! --------- Failure
            if (errorProj) then
                if (meshPairing%debug) then
                    WRITE (6, *) "        Failure of projection of node"
                end if
                nbNodeProj = 0
                cellProj%coorNodePara = 0.d0
                cellProj%coorNodeGlob = 0.d0
                iret = ERR_PAIR_PROJ
                exit
            end if

! --------- Transform coordinates of the projection in global space
            call cellPoinParaToGlob(cellTargLine, projNodePara, projNodeGlob)
            if (meshPairing%debug) then
                WRITE (6, *) "        Success of projection of node"
                WRITE (6, *) "        Coordinates of the projection of the node (local frame): ", &
                    projNodePara(1:2)
                WRITE (6, *) "        Coordinates of the projection of the node(global frame): ", &
                    projNodeGlob
            end if

! --------- Save projection of node for this cell
            cellProj%coorNodePara(1:2, iNode) = projNodePara(1:2)
            cellProj%coorNodeGlob(:, iNode) = projNodeGlob

        end do

! ----- Update parameters for projected cell
        call cellCompBary(cellProj)
        call cellCompDiam(cellProj)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! poinProjOnCell
!
! Project a point on cell - Cell is linear
!
! In  coorPoinGlob     : coordinates of point in global frame
! In  projTole         : tolerance for projection
! In  normPoint        : normal at point
! In  cellGeomTarget   : cell where to project
! Out coorProjPara     : coordinates of point in reference frame of target cell
! Out errorProj        : error code
!
! --------------------------------------------------------------------------------------------------
    subroutine poinProjOnCell(spaceDime, projTole, &
                              coorPoinGlob, normPoint, &
                              cellGeomTarget, &
                              coorProjPara, errorProj)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: spaceDime
        real(kind=8), intent(in) :: projTole
        real(kind=8), intent(in) :: coorPoinGlob(3), normPoint(3)
        type(CELL_GEOM), intent(in) :: cellGeomTarget
        real(kind=8), intent(out) :: coorProjPara(3)
        aster_logical, intent(out) :: errorProj
! ----- Local
        aster_logical, parameter :: debug = ASTER_FALSE
        type(MESH_PROJ_PARA) :: projPara
!   ------------------------------------------------------------------------------------------------
!
        coorProjPara = 0.d0
        errorProj = ASTER_FALSE
        ASSERT(cellGeomTarget%isLinear)
        ASSERT(cellGeomTarget%isSkin)

! ----- Set parameters of Newton algorithm for projection
        projPara%debug = debug
        projPara%newtIterMaxi = 75
        projPara%newtTole = projTole

! ----- Prepare object for projection on target cell
        projPara%spaceDime = spaceDime
        projPara%coorPoinGlob = coorPoinGlob
        projPara%projVect = normPoint

! ----- Projection by given vector
        call poinProjByVect(projPara, cellGeomTarget)

! ----- Copy to output results
        coorProjPara(1:2) = projPara%coorProjPara
        errorProj = projPara%errorCode .ne. 0
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! poinProjByVect
!
! Projection of a point on _skin_ cell by given vector
!
! IO  projPara         : parameters of projection
! In  cellTarget       : target cell
!
! --------------------------------------------------------------------------------------------------
    subroutine poinProjByVect(projPara, cellTarget)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PROJ_PARA), intent(inout) :: projPara
        type(CELL_GEOM), intent(in) :: cellTarget
! ----- Local
        real(kind=8), parameter :: zero = 0.d0, one = 1.d0
        integer(kind=8) :: iNode, iDime, spaceDime
        real(kind=8) :: shapeFunc(MT_NNOMAX3D), dShapeFunc(3, MT_NNOMAX3D)
        real(kind=8) :: vect_posi(3), dist
        real(kind=8) :: residu(3), matrix(3, 3), det
        real(kind=8) :: dksi(2), dbeta, beta
        integer(kind=8) :: iterNewt
        real(kind=8) :: toleAbso, toleRela, toleNewt
        real(kind=8) :: distMini, ksiMini(2), betaMini
        real(kind=8) :: refe, test
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(cellTarget%isSkin)
        spaceDime = projPara%spaceDime
        projPara%errorCode = 0
        projPara%coorProjPara = zero
        beta = one
        iterNewt = 0
        toleAbso = projPara%newtTole/100.d0
        toleRela = projPara%newtTole
        distMini = r8gaem()

! ----- Newton loop
20      continue
        vect_posi = zero
        projPara%tau1 = zero
        projPara%tau2 = zero
        matrix = zero
        residu = zero
        dksi = zero
        dbeta = zero

! ----- Get shape functions and derivated shape functions
        call elrfvf(cellTarget%cellCode, projPara%coorProjPara, shapeFunc)
        call elrfdf(cellTarget%cellCode, projPara%coorProjPara, dShapeFunc)

! ----- Position vector of current point
        do iDime = 1, 3
            do iNode = 1, cellTarget%nbNode
                vect_posi(iDime) = &
                    cellTarget%coorNodeGlob(iDime, iNode)*shapeFunc(iNode)+vect_posi(iDime)
            end do
        end do

! ----- Compute local base
        call cellCompTang(cellTarget, spaceDime, dShapeFunc, projPara%tau1, projPara%tau2)

! ----- Quantity to minimize
        do iDime = 1, 3
            vect_posi(iDime) = projPara%coorPoinGlob(iDime)-vect_posi(iDime)
        end do
        dist = sqrt(vect_posi(1)*vect_posi(1)+vect_posi(2)*vect_posi(2)+vect_posi(3)*vect_posi(3))

! ----- Newton residual
        do iDime = 1, 3
            residu(iDime) = vect_posi(iDime)-beta*projPara%projVect(iDime)
        end do

! ----- Tangent matrix (Newton)
        do iDime = 1, 3
            matrix(iDime, 1) = projPara%tau1(iDime)
            if (spaceDime .eq. 2) then
                matrix(iDime, 2) = projPara%projVect(iDime)
            elseif (spaceDime .eq. 3) then
                matrix(iDime, 2) = projPara%tau2(iDime)
                matrix(iDime, 3) = projPara%projVect(iDime)
            else
                ASSERT(ASTER_FALSE)
            end if
        end do

! ----- System determinant
        if (spaceDime .eq. 2) then
            det = matrix(1, 1)*matrix(2, 2)-matrix(1, 2)*matrix(2, 1)
        else if (spaceDime .eq. 3) then
            det = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3)-matrix(3, 2)*matrix(2, 3))- &
                  matrix(2, 1)*(matrix(1, 2)*matrix(3, 3)-matrix(3, 2)*matrix(1, 3))+ &
                  matrix(3, 1)*(matrix(1, 2)*matrix(2, 3)-matrix(2, 2)*matrix(1, 3))
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (abs(det) .le. r8prem()) then
            projPara%errorCode = 1
            goto 99
        end if

! ----- Solve system
        if (spaceDime .eq. 2) then
            dksi(1) = (residu(1)*matrix(2, 2)-residu(2)*matrix(1, 2))/det
            dksi(2) = 0.d0
            dbeta = (residu(2)*matrix(1, 1)-residu(1)*matrix(2, 1))/det
        else if (spaceDime .eq. 3) then
            dksi(1) = (residu(1)*(matrix(2, 2)*matrix(3, 3)-matrix(3, 2)*matrix(2, 3))+ &
                       residu(2)*(matrix(3, 2)*matrix(1, 3)-matrix(1, 2)*matrix(3, 3))+ &
                       residu(3)*(matrix(1, 2)*matrix(2, 3)-matrix(2, 2)*matrix(1, 3)))/det
            dksi(2) = (residu(1)*(matrix(3, 1)*matrix(2, 3)-matrix(2, 1)*matrix(3, 3))+ &
                       residu(2)*(matrix(1, 1)*matrix(3, 3)-matrix(3, 1)*matrix(1, 3))+ &
                       residu(3)*(matrix(2, 1)*matrix(1, 3)-matrix(2, 3)*matrix(1, 1)))/det
            dbeta = (residu(1)*(matrix(2, 1)*matrix(3, 2)-matrix(3, 1)*matrix(2, 2))+ &
                     residu(2)*(matrix(3, 1)*matrix(1, 2)-matrix(1, 1)*matrix(3, 2))+ &
                     residu(3)*(matrix(1, 1)*matrix(2, 2)-matrix(2, 1)*matrix(1, 2)))/det
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Update
        projPara%coorProjPara = projPara%coorProjPara+dksi
        beta = beta+dbeta

! ----- Save values if Newton avoids
        if (dist .le. distMini) then
            distMini = dist
            ksiMini = projPara%coorProjPara
            betaMini = beta
        end if

! ----- Convergence
        refe = (projPara%coorProjPara(1)*projPara%coorProjPara(1)+ &
                projPara%coorProjPara(2)*projPara%coorProjPara(2)+beta*beta)
        if (refe .le. toleRela) then
            toleNewt = toleAbso
            test = sqrt(dksi(1)*dksi(1)+dksi(2)*dksi(2)+dbeta*dbeta)
        else
            toleNewt = toleRela
            test = sqrt(dksi(1)*dksi(1)+dksi(2)*dksi(2)+dbeta*dbeta)/sqrt(refe)
        end if

! ----- Continue or not ?
        if ((test .gt. toleNewt) .and. (iterNewt .lt. projPara%newtIterMaxi)) then
            iterNewt = iterNewt+1
            goto 20
        else if ((iterNewt .ge. projPara%newtIterMaxi) .and. (test .gt. toleNewt)) then
            projPara%coorProjPara = ksiMini
            call elrfvf(cellTarget%cellCode, projPara%coorProjPara, shapeFunc)
            call elrfdf(cellTarget%cellCode, projPara%coorProjPara, dShapeFunc)
            call cellCompTang(cellTarget, spaceDime, dShapeFunc, projPara%tau1, projPara%tau2)
            projPara%errorCode = 1
        end if

! ----- End of loop
99      continue
!
        if (projPara%debug .and. projPara%errorCode .eq. 1) then
            WRITE (6, *) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            WRITE (6, *) "Newton algorithm for projection"
            WRITE (6, *) " => failure "
            WRITE (6, *) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            write (6, *) "Dimension de l'espace: ", spaceDime
            write (6, *) "Coordonnées du point à projeter: ", &
                projPara%coorPoinGlob(1), &
                projPara%coorPoinGlob(2), &
                projPara%coorPoinGlob(3)
            write (6, *) "Direction de projection: ", &
                projPara%projVect(1), &
                projPara%projVect(2), &
                projPara%projVect(3)
            write (6, *) 'Maille cible de la projection: ', &
                cellTarget%cellCode, &
                cellTarget%nbNode
            do iNode = 1, cellTarget%nbNode
                write (6, *) '  Noeud ', iNode
                write (6, *) '   (X,Y,Z)', &
                    cellTarget%coorNodeGlob(1, iNode), &
                    cellTarget%coorNodeGlob(2, iNode), &
                    cellTarget%coorNodeGlob(3, iNode)
            end do
            write (6, *) 'KSI   : ', projPara%coorProjPara(1), projPara%coorProjPara(2)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! inteCellChck
!
! Check intersection of cells (linearized)
!
! In  meshPairing      : main datastructure for pairing
! In  cellProj         : geometric properties of projected cell
! In  nbNodeProj       : number of projected nodes
! In  cellOrigLine     : geometric properties of projected cell in origin reference frame
! In  cellTargLine     : geometric properties of cell where cell has been projected
! In  normOrig         : normal to origin cell
! In  normTarg         : normal to target cell
! Out inteIsEmpty      : flag for empty intersection
!
! --------------------------------------------------------------------------------------------------
    subroutine inteCellChck(meshPairing, &
                            cellProj, nbNodeProj, &
                            cellOrigLine, cellTargLine, &
                            normOrig, normTarg, &
                            inteIsEmpty)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        type(CELL_GEOM), intent(in) :: cellProj
        integer(kind=8), intent(in) :: nbNodeProj
        type(CELL_GEOM), intent(in) :: cellOrigLine, cellTargLine
        real(kind=8), intent(in) :: normOrig(3), normTarg(3)
        aster_logical, intent(out) :: inteIsEmpty
! ----- Local
        integer(kind=8) :: iNodeProj
        integer(kind=8) :: listNodeNext(4)
        real(kind=8) :: edgeOrig(3), edgeProj(3), vectProjOrig(3)
        real(kind=8) :: sig, sp_ns_vsp, sp_np_vsp
!   ------------------------------------------------------------------------------------------------
!
        inteIsEmpty = ASTER_FALSE
        ASSERT(cellProj%isLinear)
        ASSERT(cellOrigLine%isLinear)
        ASSERT(cellTargLine%isLinear)
        ASSERT(cellProj%cellCode .eq. cellOrigLine%cellCode)

        if (meshPairing%debug) then
            WRITE (6, *) "    Check intersection of cells"
            WRITE (6, *) "    ==========================="
        end if

! ----- Set index of next nodes
        ASSERT(nbNodeProj .le. 4)
        listNodeNext = 0
        if (cellOrigLine%cellCode .eq. "SE2") then
            listNodeNext(1:2) = nodeNextSEG
        elseif (cellOrigLine%cellCode .eq. "TR3") then
            listNodeNext(1:3) = nodeNextTRIA
        elseif (cellOrigLine%cellCode .eq. "QU4") then
            listNodeNext(1:4) = nodeNextQUAD
        else
            ASSERT(ASTER_FALSE)
        end if

        do iNodeProj = 1, nbNodeProj

! --------- Original edge
            edgeOrig = 0.d0
            edgeOrig = &
                cellOrigLine%coorNodeGlob(:, listNodeNext(iNodeProj))- &
                cellOrigLine%coorNodeGlob(:, iNodeProj)

! --------- Projected edge
            edgeProj = 0.d0
            edgeProj = &
                cellProj%coorNodeGlob(:, listNodeNext(iNodeProj))- &
                cellProj%coorNodeGlob(:, iNodeProj)

! --------- Compute vector start side to projected side
            vectProjOrig = 0.d0
            vectProjOrig = &
                cellProj%coorNodeGlob(:, iNodeProj)- &
                cellOrigLine%coorNodeGlob(:, iNodeProj)

! --------- Sign of colinear product edgeProj . edgeOrig
            sig = 0.d0
            if (meshPairing%spaceDime .eq. 3) then
                sig = edgeProj(1)*edgeOrig(1)+ &
                      edgeProj(2)*edgeOrig(2)+ &
                      edgeProj(3)*edgeOrig(3)
            elseif (meshPairing%spaceDime .eq. 2) then
                sig = edgeProj(1)*edgeOrig(1)+ &
                      edgeProj(2)*edgeOrig(2)
            else
                ASSERT(ASTER_FALSE)
            end if

! --------- Sign of colinear product: direction of normal to edges
            sp_ns_vsp = 0.d0
            sp_np_vsp = 0.d0
            if (meshPairing%spaceDime .eq. 3) then
                sp_ns_vsp = normOrig(1)*vectProjOrig(1)+ &
                            normOrig(2)*vectProjOrig(2)+ &
                            normOrig(3)*vectProjOrig(3)
                sp_np_vsp = normTarg(1)*vectProjOrig(1)+ &
                            normTarg(2)*vectProjOrig(2)+ &
                            normTarg(3)*vectProjOrig(3)
            elseif (meshPairing%spaceDime .eq. 2) then
                sp_ns_vsp = normOrig(1)*vectProjOrig(1)+ &
                            normOrig(2)*vectProjOrig(2)
                sp_np_vsp = normTarg(1)*vectProjOrig(1)+ &
                            normTarg(2)*vectProjOrig(2)
            else
                ASSERT(ASTER_FALSE)
            end if

! --------- Check signs
            if (sig .lt. (0.d0-meshPairing%pairTole)) then
! ------------- The two edges are opposite
                inteIsEmpty = ASTER_TRUE
                if (meshPairing%debug) then
                    WRITE (6, *) "      Edge on original cell   : ", edgeOrig
                    WRITE (6, *) "      Edge on projected cell: ", edgeProj
                    WRITE (6, *) "      Vector between edges  : ", vectProjOrig
                    WRITE (6, *) "      Normal to edge on origin cell  : ", normOrig
                    WRITE (6, *) "      Normal to edge on target cell  : ", normTarg
                    WRITE (6, *) "      Error: the two edges are opposite"
                end if
                exit
            elseif (sp_ns_vsp .lt. (0.d0-meshPairing%pairTole) .and. &
                    sp_np_vsp .lt. (0.d0-meshPairing%pairTole)) then
                inteIsEmpty = ASTER_TRUE
                if (meshPairing%debug) then
                    WRITE (6, *) "      Edge on original cell   : ", edgeOrig
                    WRITE (6, *) "      Edge on projected cell: ", edgeProj
                    WRITE (6, *) "      Vector between edges  : ", vectProjOrig
                    WRITE (6, *) "      Normal to edge on original cell  : ", normOrig
                    WRITE (6, *) "      Normal to edge on target cell  : ", normTarg
                    WRITE (6, *) "      Error: les normales ne sont pas opposées"
                end if
                exit
            elseif (sp_ns_vsp .gt. (0.d0+meshPairing%pairTole) .and. &
                    sp_np_vsp .gt. (0.d0+meshPairing%pairTole)) then
                inteIsEmpty = ASTER_TRUE
                if (meshPairing%debug) then
                    WRITE (6, *) "      Edge on original cell   : ", edgeOrig
                    WRITE (6, *) "      Edge on projected cell: ", edgeProj
                    WRITE (6, *) "      Vector between edges  : ", vectProjOrig
                    WRITE (6, *) "      Normal to edge on original cell  : ", normOrig
                    WRITE (6, *) "      Normal to edge on target cell  : ", normTarg
                    WRITE (6, *) "      Error: les normales ne sont pas opposées"
                end if
                exit
            end if
        end do
        if (meshPairing%debug) then
        if (inteIsEmpty) then
            WRITE (6, *) "      Empty intersection"
        else
            WRITE (6, *) "      Not empty intersection"
        end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellProjOnCell
!
! Project a cell to another one: cellOrig ===> cellProj
!
! In  meshPairing      : main datastructure for pairing
! In  cellOrig         : geometric properties of cell to project
! In  cellOrigLine     : geometric properties of cell to project (linearized)
! In  cellTargLine     : geometric properties of cell where to project (linearized)
! Out cellProj         : geometric properties of projected cell
! Out nbNodeProj       : number of projected nodes
! Out iret             : return code error
!
! --------------------------------------------------------------------------------------------------
    subroutine cellProjOnCell(meshPairing, &
                              cellOrig, cellOrigLine, cellTargLine, &
                              cellProj, nbNodeProj, iret)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        type(CELL_GEOM), intent(in) :: cellOrig, cellOrigLine, cellTargLine
        type(CELL_GEOM), intent(out) :: cellProj
        integer(kind=8), intent(out) :: nbNodeProj
        integer(kind=8), intent(out) :: iret
! ----- Local
        real(kind=8) :: normOrig(3), normTarg(3)
        real(kind=8) :: ps, distCell
        aster_logical :: inteIsEmpty
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(cellTargLine%isSkin)
        ASSERT(cellTargLine%isSkin)
        ASSERT(cellTargLine%isLinear)
        ASSERT(cellOrigLine%isLinear)

        iret = ERR_PAIR_NONE
        nbNodeProj = 0
        if (meshPairing%debug) then
            WRITE (6, *) "      Compute projection of cell in another one"
            WRITE (6, *) "      ========================================="
        end if

! ----- Compute norms at barycenter
        call cellCompNormAtBary(meshPairing%spaceDime, cellTargLine, normTarg)
        call cellCompNormAtBary(meshPairing%spaceDime, cellOrigLine, normOrig)
        if (meshPairing%debug) then
            WRITE (6, *) "      Barycenter normal (target): ", normTarg
            WRITE (6, *) "      Barycenter normal (origin): ", normOrig
        end if

! ----- If cells are orthogonal -> exit
        ps = dot_product(normTarg(1:3), normOrig(1:3))
        if (abs(ps) <= meshPairing%pairTole) then
            if (meshPairing%debug) then
                WRITE (6, *) "      Cells are normal  (", abs(ps), ") => STOP"
            end if
            iret = ERR_CELL_ORTH
            goto 99
        end if

! ----- If distance between barycenter is too high -> exit
        if (meshPairing%distRatio > 0) then
            ASSERT(cellTargLine%diameter .ge. 0.d0)
            ASSERT(cellOrig%diameter .ge. 0.d0)
            distCell = norm2(cellTargLine%baryGlob-cellOrig%baryGlob)
            if (distCell >= &
                2*meshPairing%distRatio*max(cellTargLine%diameter, cellOrig%diameter)) then
                iret = ERR_CELL_OOR
                if (meshPairing%debug) then
                    WRITE (6, *) "      Cells are too far away (", distCell, ") => STOP"
                end if
                go to 99
            end if
        end if

! ----- Project nodes in target cell parametric space using raytracing
        call nodeProjOnCell(meshPairing, &
                            cellOrig, cellOrigLine, cellTargLine, &
                            cellProj, nbNodeProj, iret)
        ASSERT(cellProj%isLinear)
        if (iret .eq. ERR_PAIR_PROJ) then
            goto 99
        else
            if (meshPairing%debug) then
                WRITE (6, *) "      Success of the projection of cell"
            end if
        end if

! ----- Check if intersection is empty or not
        call inteCellChck(meshPairing, &
                          cellProj, nbNodeProj, &
                          cellOrigLine, cellTargLine, &
                          normOrig, normTarg, &
                          inteIsEmpty)
        if (inteIsEmpty) then
            iret = ERR_INTE_VOID
        end if

99      continue
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellInteProj - prjint_ray
!
! Compute intersection and projection from cell to target cell
!
! In  meshPairing      : main datastructure for pairing
! In  cellOrig         : general geometric properties of origin cell
! In  cellOrigLine     : general geometric properties of linearized origin cell
! In  cellTargLine     : general geometric properties of linearized target cell
! Out iret             : return code error
! Out nbPoinInte       : number of intersection points
! Out poinInteOrig     : coordinates of intersection points (in origin cell parametric space)
!
! --------------------------------------------------------------------------------------------------
    subroutine cellInteProj(meshPairing, &
                            cellOrig, cellOrigLine, cellTargLine, &
                            iret_, &
                            nbPoinInte_, poinInteOrig_, &
                            inteArea_, inteNeigh_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        type(CELL_GEOM), intent(in) :: cellOrig, cellOrigLine, cellTargLine
        integer(kind=8), optional, intent(out) :: iret_
        integer(kind=8), optional, intent(out) :: nbPoinInte_
        real(kind=8), optional, intent(out) :: poinInteOrig_(2, MAX_NB_INTE)
        real(kind=8), optional, intent(out) :: inteArea_
        integer(kind=8), optional, intent(out) :: inteNeigh_(MAX_NB_NEIGH)
! ----- Local
        aster_logical, parameter :: debug = ASTER_FALSE
        integer(kind=8) :: iret, inteNeigh(MAX_NB_NEIGH)
        integer(kind=8) :: nbPoinInteI
        real(kind=8) :: poinInteTargI(2, 2*MAX_NB_INTE), poinInteOrigI(2, 2*MAX_NB_INTE)
        integer(kind=8) :: nbPoinInteS
        real(kind=8) :: poinInteTargS(2, MAX_NB_INTE), poinInteOrigS(2, MAX_NB_INTE)
        integer(kind=8) :: nbPoinInte
        real(kind=8) :: poinInteOrig(2, MAX_NB_INTE)
        real(kind=8) :: inteArea
        aster_logical :: errorProj, errorInte
        type(CELL_GEOM) :: cellProj
        integer(kind=8) :: nbCmpPara, nbNodeProj
!   ------------------------------------------------------------------------------------------------
!
        inteNeigh = 0
        nbPoinInte = 0
        poinInteOrig = 0.d0
        inteArea = 0.d0
        iret = ERR_PAIR_NONE
        nbCmpPara = meshPairing%spaceDime-1

! ----- Compute projection of cell on target cell
        errorProj = ASTER_FALSE
        call cellProjOnCell(meshPairing, &
                            cellOrig, cellOrigLine, cellTargLine, &
                            cellProj, nbNodeProj, iret)
        if (iret .ne. ERR_PAIR_NONE) then
            errorProj = ASTER_TRUE
            goto 100
        end if
        if (meshPairing%debug) then
            WRITE (6, *) "     Intersection is OK"
        end if

! ----- Compute intersection in parametric space of target cell
        nbPoinInteI = 0
        poinInteTargI = 0.d0
        poinInteOrigI = 0.d0
        if (meshPairing%spaceDime .eq. 2) then
            ASSERT(nbNodeProj .eq. 2)
            call addPoinInte2D(meshPairing, cellProj, &
                               nbPoinInteI, poinInteTargI, poinInteOrigI, &
                               inteNeigh)
        elseif (meshPairing%spaceDime .eq. 3) then
            call addPoinInte3D(meshPairing, cellProj, nbNodeProj, &
                               cellOrigLine, cellTargLine, &
                               nbPoinInteI, poinInteTargI, poinInteOrigI, &
                               inteNeigh, iret)
            call addPoinOnEdge(meshPairing, cellProj, nbNodeProj, &
                               cellOrigLine, cellTargLine, &
                               nbPoinInteI, poinInteTargI, poinInteOrigI, &
                               inteNeigh)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Trop de points d'intersection (avant tri) => on ne sait pas gérer
        ASSERT(nbPoinInteI .le. 2*MAX_NB_INTE)

! ----- Error management
        errorInte = ASTER_FALSE
        if (nbPoinInteI == 0 .or. iret == ERR_PAIR_MAST) then
            errorInte = ASTER_TRUE
            goto 100
        end if

! ----- Debug print
        if (meshPairing%debug) then
            WRITE (6, *) "     Intersection (before re-ordering): ", nbPoinInteI, &
                " points of intersection"
            WRITE (6, *) "     Target side : ", poinInteTargI(1:nbCmpPara, 1:nbPoinInteI)
            WRITE (6, *) "     Origin side : ", poinInteOrigI(1:nbCmpPara, 1:nbPoinInteI)
        end if

! ----- Sort list of intersection points
        if ((nbPoinInteI .gt. 2 .and. meshPairing%spaceDime == 3) .or. &
            (nbPoinInteI .ge. 2 .and. meshPairing%spaceDime == 2)) then
            call intePoinSort(meshPairing, &
                              nbPoinInteI, poinInteTargI, poinInteOrigI, &
                              nbPoinInteS, poinInteTargS, poinInteOrigS)
        else
            ASSERT(nbPoinInteI .le. MAX_NB_INTE)
            nbPoinInteS = nbPoinInteI
            poinInteTargS(:, 1:nbPoinInteI) = poinInteTargI(:, 1:nbPoinInteI)
            poinInteOrigS(:, 1:nbPoinInteI) = poinInteOrigI(:, 1:nbPoinInteI)
        end if

! ----- Error management
        if (nbPoinInteS > MAX_NB_INTE) then
            if (meshPairing%debug) then
                WRITE (6, *) "     Intersection (after re-ordering): ", nbPoinInteS, &
                    " points of intersection"
                WRITE (6, *) "     Too many points !"
            end if
            iret = ERR_PAIR_SLAV
            go to 99
        end if

! ----- Debug print
        if (meshPairing%debug) then
            WRITE (6, *) "     Intersection (after re-ordering): ", nbPoinInteS, &
                " points of intersection"
            WRITE (6, *) "     Target side : ", poinInteTargS(1:nbCmpPara, 1:nbPoinInteS)
            WRITE (6, *) "     Origin side : ", poinInteOrigS(1:nbCmpPara, 1:nbPoinInteS)
        end if

! ----- All nodes have to be inside original cell
        call intePoinInCell(meshPairing, cellOrigLine, &
                            nbPoinInteS, poinInteOrigS, poinInteOrig)
        nbPoinInte = nbPoinInteS

! ----- Compute weight of intersection
        inteArea = 0.d0
        if (nbPoinInte .ne. 0) then
            call inteCellArea(meshPairing%spaceDime, nbPoinInte, poinInteOrig, inteArea)
            if (meshPairing%debug) then
                WRITE (6, *) "     Area of intersection : ", inteArea
            end if
        end if

! ----- Error
100     continue
        if (errorProj .or. errorInte) then
            if (meshPairing%debug) then
                WRITE (6, *) "     Failure of intersection"
            end if
            nbPoinInte = 0
            poinInteOrig = 0.d0
            inteNeigh = 0
            inteArea = 0.d0
        else
            if (meshPairing%debug) then
                WRITE (6, *) "     Success of intersection"
            end if
        end if

! ----- Fatal error only in debug mode
99      continue
        if (debug) then
            ASSERT(iret .ne. ERR_PAIR_NONE)
        end if

! ----- Outputs
        if (present(inteNeigh_)) then
            inteNeigh_ = inteNeigh
        end if
        if (present(nbPoinInte_)) then
            nbPoinInte_ = nbPoinInte
        end if
        if (present(inteArea_)) then
            inteArea_ = inteArea
        end if
        if (present(poinInteOrig_)) then
            poinInteOrig_ = poinInteOrig
        end if
        if (present(iret_)) then
            iret_ = iret
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! addPoinInte2D
!
! Compute intersection in parametric space of target cell - For 2D case (only segments)
!
! In  meshPairing      : main datastructure for pairing
! In  cellProj         : geometric properties of original cell
! Out nbPoinInte       : number of intersection points
! Out poinInteTarg     : coordinates of intersection points in target cell
! Out poinInteOrig     : coordinates of intersection points in original cell
! Out inteNeigh        : active cell neighbours for each node of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine addPoinInte2D(meshPairing, cellProj, &
                             nbPoinInte, poinInteTarg, poinInteOrig, &
                             inteNeigh)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        type(CELL_GEOM), intent(in) :: cellProj
        integer(kind=8), intent(out) :: nbPoinInte
        real(kind=8), intent(out) :: poinInteTarg(2, 2*MAX_NB_INTE), poinInteOrig(2, 2*MAX_NB_INTE)
        integer(kind=8), intent(out) :: inteNeigh(MAX_NB_NEIGH)
! ----- Local
        integer(kind=8) :: iNode, nbNode
        real(kind=8) :: ksi
        real(kind=8), parameter :: coorSegPara(2) = (/-1.d0, +1.d0/)
        real(kind=8) :: ksiMini, ksiMaxi, t1
!   ------------------------------------------------------------------------------------------------
!
        nbPoinInte = 0
        poinInteTarg = 0.d0
        poinInteOrig = 0.d0
        inteNeigh = 0
        nbNode = cellProj%nbNode
        ASSERT(nbNode .eq. 2)

! ----- Add projected nodes inside target cell
        do iNode = 1, nbNode
            ksi = cellProj%coorNodePara(1, iNode)
            if (ksi > (-1.d0-meshPairing%pairTole) .and. &
                ksi < (1.d0+meshPairing%pairTole)) then
                nbPoinInte = nbPoinInte+1
                poinInteTarg(1, nbPoinInte) = ksi
                poinInteOrig(1, nbPoinInte) = coorSegPara(iNode)
                inteNeigh(iNode) = 1
            end if
        end do

! ----- Add target nodes if they are inside projected cell
        ksiMini = min(cellProj%coorNodePara(1, 1), cellProj%coorNodePara(1, 2))
        ksiMaxi = max(cellProj%coorNodePara(1, 1), cellProj%coorNodePara(1, 2))
        if ((ksiMini-meshPairing%pairTole) <= -1.d0 .and. &
            -1.d0 <= (ksiMaxi+meshPairing%pairTole)) then
            nbPoinInte = nbPoinInte+1
            poinInteTarg(1, nbPoinInte) = coorSegPara(1)
            if (abs(-1.d0-cellProj%coorNodePara(1, 1)) < meshPairing%pairTole) then
                poinInteOrig(1, nbPoinInte) = coorSegPara(1)
            else
                t1 = (-1.d0-cellProj%coorNodePara(1, 1))/ &
                     (cellProj%coorNodePara(1, 2)-cellProj%coorNodePara(1, 1))
                poinInteOrig(1, nbPoinInte) = 2*t1-1
            end if
        end if
        if ((ksiMini-meshPairing%pairTole) <= 1.d0 .and. &
            1.d0 <= (ksiMaxi+meshPairing%pairTole)) then
            nbPoinInte = nbPoinInte+1
            poinInteTarg(1, nbPoinInte) = coorSegPara(2)
            if (abs(1.d0-cellProj%coorNodePara(1, 1)) .lt. meshPairing%pairTole) then
                poinInteOrig(1, nbPoinInte) = coorSegPara(1)
            else
                t1 = (1.d0-cellProj%coorNodePara(1, 1))/ &
                     (cellProj%coorNodePara(1, 2)-cellProj%coorNodePara(1, 1))
                poinInteOrig(1, nbPoinInte) = 2*t1-1
            end if
        end if
        if (meshPairing%debug) then
            WRITE (6, *) "     Neighbours  : ", inteNeigh
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! addPoinInte3D
!
! Compute intersection in parametric master space - For 3D case
!
! In  meshPairing      : main datastructure for pairing
! In  cellProj         : geometric properties of projected cell from slave cell
! In  nbNodeProj       : number of projected nodes
! In  cellOrigLine     : geometric properties of cell to project (linearized)
! In  cellTargLine     : geometric properties of cell where to project (linearized)
! Out nbPoinInte       : number of intersection points
! Out poinInteTarg     : coordinates of intersection points in target cell
! Out poinInteOrig     : coordinates of intersection points in original cell
! Out inteNeigh        : active cell neighbours for each node of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine addPoinInte3D(meshPairing, cellProj, nbNodeProj, &
                             cellOrigLine, cellTargLine, &
                             nbPoinInte, poinInteTarg, poinInteOrig, &
                             inteNeigh, iret)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        type(CELL_GEOM), intent(in) :: cellProj
        integer(kind=8), intent(in) :: nbNodeProj
        type(CELL_GEOM), intent(in) :: cellOrigLine, cellTargLine
        integer(kind=8), intent(out) :: nbPoinInte
        real(kind=8), intent(out) :: poinInteTarg(2, 2*MAX_NB_INTE), poinInteOrig(2, 2*MAX_NB_INTE)
        integer(kind=8), intent(out) :: inteNeigh(MAX_NB_NEIGH)
        integer(kind=8), intent(out) :: iret
! ----- Local
        aster_logical :: pointIsInside
        real(kind=8) :: nodeProjPara(2), origCoorPara(2)
        integer(kind=8) :: listNodePrev(4)
        integer(kind=8) :: iNodeTarg, iNodeProj, errorInside
!   ------------------------------------------------------------------------------------------------
!
        inteNeigh = 0
        nbPoinInte = 0
        iret = ERR_PAIR_NONE
        if (meshPairing%debug) then
            WRITE (6, *) "Inter. - Détection de la position des noeuds projetés."
        end if

! ----- Liste orientée des noeuds projetés de la cellule
        ASSERT(nbNodeProj .le. 4)
        listNodePrev = 0
        if (cellOrigLine%cellCode .eq. "TR3") then
            listNodePrev(1:3) = nodePrevTRIA
        elseif (cellOrigLine%cellCode .eq. "QU4") then
            listNodePrev(1:4) = nodePrevQUAD
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Add projection of original nodes inside target cell (parametric space)
        do iNodeProj = 1, nbNodeProj
            nodeProjPara(1:2) = cellProj%coorNodePara(1:2, iNodeProj)
            if (meshPairing%debug) then
                WRITE (6, *) "Project original nodes inside target cell: ", &
                    iNodeProj, nodeProjPara(1:2)
                WRITE (6, *) " Target cell type:", cellTargLine%cellCode
            end if
            pointIsInside = cellPoinInside(cellTargLine, meshPairing%pairTole, nodeProjPara)
            if (meshPairing%debug) then
                if (pointIsInside) then
                    WRITE (6, *) " => this node is inside target cell"
                else
                    WRITE (6, *) " => this node is NOT inside target cell"
                end if
            end if
            if (pointIsInside) then
                nbPoinInte = nbPoinInte+1
                ASSERT(nbPoinInte .le. MAX_NB_INTE)
                poinInteTarg(:, nbPoinInte) = nodeProjPara
                poinInteOrig(:, nbPoinInte) = cellOrigLine%coorNodePara(1:2, iNodeProj)
                if (meshPairing%debug) then
                    WRITE (6, *) " nbPoinInte :", nbPoinInte
                    WRITE (6, *) " Coor. targ.:", nodeProjPara
                    WRITE (6, *) " Coor. orig.:", cellOrigLine%coorNodePara(1:2, iNodeProj)
                end if
                inteNeigh(iNodeProj) = 1
                inteNeigh(listNodePrev(iNodeProj)) = 1
            end if
        end do

! ----- Add targeted nodes if they are inside projected original cell
        do iNodeTarg = 1, cellTargLine%nbNode
! --------- Current parametric coordinates of targeted node
            nodeProjPara = cellTargLine%coorNodePara(1:2, iNodeTarg)
            if (meshPairing%debug) then
                WRITE (6, *) "Project targeted nodes inside original cell: ", &
                    iNodeTarg, nodeProjPara(1:2)
                WRITE (6, *) " on cell ", cellProj%cellCode
            end if

! --------- Test
            call poinIsInsideOtherCell(meshPairing%pairTole, meshPairing%debug, &
                                       nodeProjPara, cellProj, &
                                       origCoorPara, pointIsInside, &
                                       errorInside)
            if (meshPairing%debug) then
                if (pointIsInside) then
                    WRITE (6, *) " => this node is inside original cell"
                else
                    WRITE (6, *) " => this node is NOT inside original cell"
                end if
            end if

            if (pointIsInside) then
                nbPoinInte = nbPoinInte+1
                poinInteTarg(:, nbPoinInte) = nodeProjPara
                poinInteOrig(:, nbPoinInte) = origCoorPara
                if (meshPairing%debug) then
                    WRITE (6, *) " nbPoinInte :", nbPoinInte
                    WRITE (6, *) " Coor. targ.:", nodeProjPara
                    WRITE (6, *) " Coor. orig.:", origCoorPara
                end if
            end if
            if (errorInside == ERR_CELL_DEGE) then
                iret = ERR_PAIR_MAST
                if (meshPairing%debug) then
                    WRITE (6, *) "Error during projection"
                end if
            end if
        end do
        if (meshPairing%debug) then
            WRITE (6, *) "Inter. - Détection de la position des noeuds projetés."
            WRITE (6, *) "     Neighbours after points intersection : ", inteNeigh
            WRITE (6, *) "     Nb pts intersection after points intersection : ", nbPoinInte
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! poinIsInsideOtherCell
!
! Detect if point is inside other cell
!
! In  projTole         : tolerance for projection
! In  poinCoorPara     : coordinates of point in parametric space
! In  cellGeom         : geometric properties of cell
! Out origCoorPara
! Out poinIsInside
! Out error
!
! --------------------------------------------------------------------------------------------------
    subroutine poinIsInsideOtherCell(projTole, debug, &
                                     poinCoorPara, cellGeom, &
                                     origCoorPara, poinIsInside, &
                                     error)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: projTole
        aster_logical, intent(in) :: debug
        type(CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(in) :: poinCoorPara(2)
        real(kind=8), intent(out) :: origCoorPara(2)
        aster_logical, intent(out) :: poinIsInside
        integer(kind=8), intent(out) :: error
! ----- Local
        character(len=8) :: cellCode
        real(kind=8) :: cellCoorPara(2, MT_NNOMAX3D)
        real(kind=8) :: v0(2), v1(2), v2(2), d00, d10, d11, m, u, v
        real(kind=8) :: d02, d12
!   ------------------------------------------------------------------------------------------------
!
        poinIsInside = ASTER_FALSE
        error = ERR_PAIR_NONE
        cellCode = cellGeom%cellCode
        cellCoorPara = cellGeom%coorNodePara(1:2, :)
        ASSERT(celLGeom%isSkin)
        ASSERT(cellGeom%isLinear)
        ASSERT(cellGeom%nbNode .le. MT_NNOMAX2D)
        origCoorPara = 0.d0

! ----- First detect with triangle

! ----- Vectorial basis for element
        v0 = cellCoorPara(:, 2)-cellCoorPara(:, 1)
        v1 = cellCoorPara(:, 3)-cellCoorPara(:, 1)
        d00 = v0(1)*v0(1)+v0(2)*v0(2)
        d10 = v0(1)*v1(1)+v0(2)*v1(2)
        d11 = v1(1)*v1(1)+v1(2)*v1(2)
        m = (d00*d11-d10*d10)
        if (debug) then
            WRITE (6, *) "Base vectorielle: ", v0, v1, m
        end if

! ----- Degenerated vectorial basis for element (colinear vectors) => exit
        if (abs(m) .le. projTole) then
            poinIsInside = ASTER_FALSE
            error = ERR_CELL_DEGE
            goto 99
        end if

! ----- Coordinates of point in element's basis
        v2 = poinCoorPara(:)-cellCoorPara(:, 1)
        d02 = v0(1)*v2(1)+v0(2)*v2(2)
        d12 = v1(1)*v2(1)+v1(2)*v2(2)
        if (debug) then
            WRITE (6, *) "Base vectorielle 2: ", v2
        end if

! ----- Point is in element => exit
        if (sqrt(v2(1)**2+v2(2)**2) .le. 0.d0+projTole) then
            poinIsInside = ASTER_TRUE
            if (cellCode .eq. 'TR3') then
                origCoorPara(1) = 0.d0
                origCoorPara(2) = 0.d0
            elseif (cellCode .eq. 'QU4') then
                origCoorPara(1) = -1.d0
                origCoorPara(2) = -1.d0
            else
                ASSERT(ASTER_FALSE)
            end if
            if (debug) then
                WRITE (6, *) "Dedans sans extension"
            end if
            goto 99
        end if

! ----- Extension with projTole
        u = 1/m*(d11*d02-d10*d12)
        v = 1/m*(d00*d12-d10*d02)
        if (debug) then
            WRITE (6, *) "Extension T3: ", u, v, projTole
        end if
        if (u .ge. (0.d0-projTole) .and. &
            v .ge. (0.d0-projTole) .and. &
            (u+v) .le. (1.d0+projTole)) then
            poinIsInside = ASTER_TRUE
            if (cellCode .eq. 'TR3') then
                origCoorPara(1) = 0.d0+u
                origCoorPara(2) = 0.d0+v
            elseif (cellCode .eq. 'QU4') then
                origCoorPara(1) = -1.d0+u*2+v*2
                origCoorPara(2) = -1.d0+v*2
            else
                ASSERT(ASTER_FALSE)
            end if
            if (debug) then
                WRITE (6, *) "Dedans avec extension"
            end if
            goto 99
        else
            poinIsInside = ASTER_FALSE
        end if

        if (cellCode .eq. 'QU4') then
! --------- Vectorial basis for element
            v0(:) = cellCoorPara(:, 3)-cellCoorPara(:, 1)
            v1(:) = cellCoorPara(:, 4)-cellCoorPara(:, 1)
            d00 = v0(1)*v0(1)+v0(2)*v0(2)
            d10 = v0(1)*v1(1)+v0(2)*v1(2)
            d11 = v1(1)*v1(1)+v1(2)*v1(2)
            m = (d00*d11-d10*d10)
            if (debug) then
                WRITE (6, *) "Base vectorielle Q4-1: ", v0, v1, m
            end if

! --------- Degenerated vectorial basis for element (colinear vectors) => exit
            if (abs(m) .le. projTole) then
                poinIsInside = ASTER_FALSE
                error = ERR_CELL_DEGE
                goto 99
            end if

! --------- Coordinates of point in element's basis
            v2(:) = poinCoorPara(:)-cellCoorPara(:, 1)
            d02 = v0(1)*v2(1)+v0(2)*v2(2)
            d12 = v1(1)*v2(1)+v1(2)*v2(2)
            if (debug) then
                WRITE (6, *) "Base vectorielle Q4-2: ", v2
            end if

! --------- Point is in element => exit
            if (sqrt(v2(1)**2+v2(2)**2) .le. 0.d0+projTole) then
                poinIsInside = ASTER_TRUE
                if (cellCode .eq. 'QU4') then
                    origCoorPara(1) = -1.d0
                    origCoorPara(2) = -1.d0
                end if
                if (debug) then
                    WRITE (6, *) "Dedans sans extension"
                end if
                goto 99
            end if

! --------- Extension with projTole
            u = 1/m*(d11*d02-d10*d12)
            v = 1/m*(d00*d12-d10*d02)
            if (debug) then
                WRITE (6, *) "Extension Q4: ", u, v, projTole
            end if
            if (u .ge. (0.d0-projTole) .and. &
                v .ge. (0.d0-projTole) .and. &
                (u+v) .le. (1.d0+projTole)) then
                poinIsInside = ASTER_TRUE
                if (cellCode .eq. 'QU4') then
                    origCoorPara(1) = -1.d0+u*2
                    origCoorPara(2) = -1.d0+v*2+u*2
                end if
                if (debug) then
                    WRITE (6, *) "Dedans avec extension"
                end if
                goto 99
            else
                poinIsInside = ASTER_FALSE
            end if
        end if
!
99      continue
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! inteCellSegm
!
! Compute intersection between segment and (linearized) cell
!
! In  meshPairing      : main datastructure for pairing
! In  coorSegm         : coordinates of segment
! In  cellLine         : general geometric properties of linearized cell
! IO  nbPoinInte       : number of intersection points
! IO  poinInte         : coordinates of intersection points
!
! --------------------------------------------------------------------------------------------------
    subroutine inteCellSegm(meshPairing, &
                            coorSegm, cellLine, &
                            nbPoinInte, poinInte)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        real(kind=8), intent(in) :: coorSegm(2, 2)
        type(CELL_GEOM), intent(in) :: cellLine
        integer(kind=8), intent(inout) :: nbPoinInte
        real(kind=8), intent(inout) :: poinInte(2, 2*MAX_NB_INTE)
! ----- Local
        integer(kind=8) :: nbNode, iNode
        real(kind=8) :: a, b, c, d
        real(kind=8) :: t1, t2, det, aux(2), norm
        real(kind=8) :: x1, y1, x2, y2
        integer(kind=8) :: listNodeNext(4)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(cellLine%isLinear)
        nbNode = cellLine%nbNode

! ----- Set index of next nodes
        ASSERT(nbNode .le. 4)
        listNodeNext = 0
        if (cellLine%cellCode .eq. "SE2") then
            listNodeNext(1:2) = nodeNextSEG
        elseif (cellLine%cellCode .eq. "TR3") then
            listNodeNext(1:3) = nodeNextTRIA
        elseif (cellLine%cellCode .eq. "QU4") then
            listNodeNext(1:4) = nodeNextQUAD
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Coefficients for parametric equation of segment
        x1 = coorSegm(1, 1)
        y1 = coorSegm(2, 1)
        a = coorSegm(1, 2)-x1
        b = coorSegm(2, 2)-y1

        if (meshPairing%debug) then
            WRITE (6, *) " Intersection"
            WRITE (6, *) "   Segment 1: ", x1, y1
            WRITE (6, *) "    Parametric equation: ", a, b
        end if

! ----- Loop on edges of element
        do iNode = 1, nbNode
! --------- Current segment in cell
            x2 = cellLine%coorNodePara(1, iNode)
            y2 = cellLine%coorNodePara(2, iNode)

! --------- Coefficients for parametric equation of segment in cell
            c = cellLine%coorNodePara(1, listNodeNext(iNode))-x2
            d = cellLine%coorNodePara(2, listNodeNext(iNode))-y2

! --------- Compute intersection
            det = b*c-a*d
            if (meshPairing%debug) then
                WRITE (6, *) "   Segment : ", iNode, x2, y2
                WRITE (6, *) "    Parametric equation: ", c, d
            end if

            if (sqrt(det**2) .gt. meshPairing%pairTole) then
                t1 = 1/det*(d*(x1-x2)-c*(y1-y2))
                t2 = 1/det*(b*(x1-x2)-a*(y1-y2))
                aux(1) = (-t1*a-t2*c)-(x1-x2)
                aux(2) = (t1*b+t2*d)-(y1-y2)
                norm = sqrt(aux(1)**2+aux(2)**2)
            else
                t1 = -1.d0
                t2 = -1.d0
            end if

! --------- Test intersection
            if (t1 .lt. 1.d0+meshPairing%pairTole .and. &
                t1 .gt. 0.d0-meshPairing%pairTole .and. &
                t2 .lt. 1.d0+meshPairing%pairTole .and. &
                t2 .gt. 0.d0-meshPairing%pairTole) then
                nbPoinInte = nbPoinInte+1
                ASSERT(nbPoinInte .gt. 0)
                ASSERT(nbPoinInte .le. 2*MAX_NB_INTE)
                poinInte(1, nbPoinInte) = (t2*c+x2+t1*a+x1)/2.d0
                poinInte(2, nbPoinInte) = (t2*d+y2+t1*b+y1)/2.d0
                if (meshPairing%debug) then
                    WRITE (6, *) " -> intersection détectée: ", &
                        nbPoinInte, poinInte(:, nbPoinInte)
                end if
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! addPoinOnEdge
!
! Compute intersection of edges for 3D case
!
! In  meshPairing      : main datastructure for pairing
! In  cellProj         : geometric properties of projected cell from slave cell
! In  nbNodeProj       : number of projected nodes
! In  cellOrigLine     : general geometric properties of linearized origin cell
! In  cellTargLine     : general geometric properties of linearized target cell
! Out nbPoinInte       : number of intersection points
! Out poinInteTarg     : coordinates of intersection points on target cell
! Out poinInteOrig     : coordinates of intersection points on origin cell
! IO  inteNeigh        : active cell neighbours for each node of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine addPoinOnEdge(meshPairing, cellProj, nbNodeProj, &
                             cellOrigLine, cellTargLine, &
                             nbPoinInte, poinInteTarg, poinInteOrig, &
                             inteNeigh)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        type(CELL_GEOM), intent(in) :: cellProj
        integer(kind=8), intent(in) :: nbNodeProj
        type(CELL_GEOM), intent(in) :: cellOrigLine, cellTargLine
        integer(kind=8), intent(inout) :: nbPoinInte
        real(kind=8), intent(inout) :: poinInteTarg(2, 2*MAX_NB_INTE)
        real(kind=8), intent(inout) :: poinInteOrig(2, 2*MAX_NB_INTE)
        integer(kind=8), intent(inout) :: inteNeigh(MAX_NB_NEIGH)
! ----- Local
        real(kind=8) :: coorSegmProj(2, 2), t1, t2
        real(kind=8) :: coorSegmOrig(2, 2)
        integer(kind=8) :: nbPoinAdd, nbPoinInteAv
        integer(kind=8) :: iPoinAdd, listNodeNext(4), iNodeProj
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(nbNodeProj .le. 4)

        if (meshPairing%debug) then
            WRITE (6, *) "Intersection des mailles de peau - Compute intersection of edges."
        end if

! ----- Set index of next nodes
        listNodeNext = 0
        if (cellOrigLine%cellCode .eq. "SE2") then
            listNodeNext(1:2) = nodeNextSEG
        elseif (cellOrigLine%cellCode .eq. "TR3") then
            listNodeNext(1:3) = nodeNextTRIA
        elseif (cellOrigLine%cellCode .eq. "QU4") then
            listNodeNext(1:4) = nodeNextQUAD
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Intersection of edges
        do iNodeProj = 1, nbNodeProj

! --------- Coordinates of segment from projected cell (in target cell reference frame)
            coorSegmProj(1:2, 1) = cellProj%coorNodePara(1:2, iNodeProj)
            coorSegmProj(1:2, 2) = cellProj%coorNodePara(1:2, listNodeNext(iNodeProj))

! --------- Coordinates of segment from original cell (in original cell reference frame)
            coorSegmOrig(1:2, 1) = cellOrigLine%coorNodePara(1:2, iNodeProj)
            coorSegmOrig(1:2, 2) = cellOrigLine%coorNodePara(1:2, listNodeNext(iNodeProj))

            if (meshPairing%debug) then
                WRITE (6, *) "Intersection of edges: ", &
                    iNodeProj, coorSegmProj
                WRITE (6, *) " Intersection points before: ", &
                    nbPoinInte, poinInteTarg
            end if

! --------- Compute intersection between edge of master and projected slave cells
            nbPoinInteAv = nbPoinInte
            call inteCellSegm(meshPairing, &
                              coorSegmProj, cellTargLine, &
                              nbPoinInte, poinInteTarg)
            if (meshPairing%debug) then
                WRITE (6, *) " => intersection points: ", &
                    nbPoinInte, poinInteTarg
            end if

! --------- Number of intersection points to add
            nbPoinAdd = nbPoinInte-nbPoinInteAv
            if (nbPoinAdd .gt. 0) then
                inteNeigh(iNodeProj) = 1
                do iPoinAdd = 1, nbPoinAdd
                    t1 = 0.d0
                    t2 = 0.d0
                    if (abs(poinInteTarg(1, nbPoinInteAv+iPoinAdd)-coorSegmProj(1, 1)) .gt. &
                        meshPairing%pairTole) then
                        t1 = (coorSegmProj(1, 2)-coorSegmProj(1, 1))/ &
                             (poinInteTarg(1, nbPoinInteAv+iPoinAdd)-coorSegmProj(1, 1))
                        poinInteOrig(1, nbPoinInteAv+iPoinAdd) = &
                            (1.0/t1)*(coorSegmOrig(1, 2)-coorSegmOrig(1, 1))+coorSegmOrig(1, 1)
                    else
                        poinInteOrig(1, nbPoinInteAv+iPoinAdd) = &
                            coorSegmOrig(1, 1)
                    end if
                    if (abs(poinInteTarg(2, nbPoinInteAv+iPoinAdd)-coorSegmProj(2, 1)) .gt. &
                        meshPairing%pairTole) then
                        t2 = (coorSegmProj(2, 2)-coorSegmProj(2, 1))/ &
                             (poinInteTarg(2, nbPoinInteAv+iPoinAdd)-coorSegmProj(2, 1))
                        poinInteOrig(2, nbPoinInteAv+iPoinAdd) = &
                            (1.0/t2)*(coorSegmOrig(2, 2)-coorSegmOrig(2, 1))+coorSegmOrig(2, 1)
                        poinInteOrig(1, nbPoinInteAv+iPoinAdd) = &
                            (1.0/t2)*(coorSegmOrig(1, 2)-coorSegmOrig(1, 1))+coorSegmOrig(1, 1)
                    else
                        if (t1 .lt. meshPairing%pairTole) then
                            poinInteOrig(2, nbPoinInteAv+iPoinAdd) = &
                                coorSegmOrig(2, 1)
                        else
                            poinInteOrig(2, nbPoinInteAv+iPoinAdd) = &
                                (1.0/t1)*(coorSegmOrig(2, 2)-coorSegmOrig(2, 1))+coorSegmOrig(2, 1)
                        end if
                    end if
                end do
            end if
        end do
        if (meshPairing%debug) then
            WRITE (6, *) "Intersection des mailles de peau - Compute intersection of edges."
            WRITE (6, *) "     Neighbours after segments intersection : ", inteNeigh
            WRITE (6, *) "     Nb pts intersection after points intersection : ", nbPoinInte
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pairGetStartCells - ap_infast_n
!
! Get cell for starting search
!
! In  meshPairing      : main datastructure for pairing
! Ptr nodeCoor         : updated coordinates of nodes of mesh
! Ptr meshTypeGeom     : pointer to type of cells in mesh
! Ptr meshConx         : pointers to connectivity of mesh
!     meshConxCumu
! Ptr mastConxInv      : pointers to inverse conenctivity of master cells
!     mastConxInvCumu
! In  nbCellSlav       : number of slave cells
! In  listCellSlav     : list of slave cells
! In  cellSlavFlag     : flag for slave cell has been used in a pair
!                        0 - Never used
!                        1 - Used as starting point
!                        2 - Used
! In  nbCellMast       : number of master cells
! In  listCellMast     : list of master cells
! In  nbNodeMast       : number of master nodes
! In  listNodeMast     : list of master nodes
! In  nbMastStart      : number of master cells used as a start cell
! In  cellMastStart    : list of master cells used as a start cell
! In  nbSlavStart      : number of slave cells used as a start cell
! In  cellSlavStart    : list of slave cells used as a start cell
!
! --------------------------------------------------------------------------------------------------
    subroutine pairGetStartCells(meshPairing, nodeCoor, &
                                 meshTypeGeom, meshConx, meshConxCumu, &
                                 mastConxInv, mastConxInvCumu, &
                                 nbCellSlav, listCellSlav, cellSlavFlag, &
                                 nbCellMast, listCellMast, &
                                 nbNodeMast, listNodeMast, &
                                 nbMastStart, cellMastStart, &
                                 nbSlavStart, cellSlavStart)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        real(kind=8), pointer :: nodeCoor(:)
        integer(kind=8), pointer :: meshTypeGeom(:)
        integer(kind=8), pointer :: meshConx(:), meshConxCumu(:)
        integer(kind=8), pointer :: mastConxInv(:), mastConxInvCumu(:)
        integer(kind=8), intent(in) :: nbCellSlav, listCellSlav(nbCellSlav)
        integer(kind=8), pointer :: cellSlavFlag(:)
        integer(kind=8), intent(in) :: nbCellMast, listCellMast(nbCellMast)
        integer(kind=8), intent(in) :: nbNodeMast, listNodeMast(nbNodeMast)
        integer(kind=8), intent(out) :: nbMastStart, cellMastStart(nbCellSlav)
        integer(kind=8), intent(out) :: nbSlavStart, cellSlavStart(nbCellSlav)
! ----- Local
        integer(kind=8) :: iCellSlav, iCellMast
        integer(kind=8) :: cellSlavNume, cellSlavIndx
        integer(kind=8) :: cellMastNume
        type(CELL_GEOM) :: cellSlav, cellSlavLine
        type(CELL_GEOM) :: cellMast, cellMastLine
        integer(kind=8) :: nodeNumeClosest
        integer(kind=8) :: slavIndxMini
        integer(kind=8) :: nbCellToNode, cellToNode(nbCellMast)
        real(kind=8) :: inteArea
!   ------------------------------------------------------------------------------------------------
!
        slavIndxMini = minval(listCellSlav)
        nbMastStart = 0
        cellMastStart = 0
        nbSlavStart = 0
        cellSlavStart = 0

! ----- Loop on slave elements
        do iCellSlav = 1, nbCellSlav
! --------- Get current slave element
            cellSlavNume = listCellSlav(iCellSlav)
            cellSlavIndx = cellSlavNume+1-slavIndxMini
            if (meshPairing%debug) then
                WRITE (6, *) "  Current slave cell : ", cellSlavNume
            end if

! --------- Already tracked ?
            if (cellSlavFlag(cellSlavIndx) .eq. 0) then
                if (meshPairing%debug) then
                    WRITE (6, *) "  Current slave cell not yet tracked"
                end if

! ------------- Create slave cell
                call cellCreate(cellSlavNume, nodeCoor, &
                                meshTypeGeom, meshConx, meshConxCumu, &
                                cellSlav, cellSlavLine)
                if (meshPairing%debug) then
                    WRITE (6, *) "Properties of slave cell"
                    WRITE (6, *) "========================"
                    call cellDebug(cellSlav)
                end if

! ------------- Find the closest master node from center of slave cell
                if (meshPairing%debug) then
                    WRITE (6, *) "  Seek for closest master node from center of slave cell"
                end if
                call getClosestNodesFromCell(cellSlav, nodeCoor, &
                                             nbNodeMast, listNodeMast, &
                                             nodeNumeClosest)
                if (meshPairing%debug) then
                    WRITE (6, *) "  => ", nodeNumeClosest
                end if

! ------------- Construct list of master cells attached to this node
                if (meshPairing%debug) then
                    WRITE (6, *) "  Get list of cells attached to this node"
                end if
                call getCellsFromNode(nodeNumeClosest, &
                                      mastConxInv, mastConxInvCumu, &
                                      nbCellMast, listCellMast, &
                                      nbCellToNode, cellToNode)
                if (meshPairing%debug) then
                    WRITE (6, *) "  => Number of cells: ", nbCellToNode
                    WRITE (6, *) "  => Cells connected: ", cellToNode(1:nbCellToNode)
                end if

! ------------- Loop on master elements linked to the closest master node
                do iCellMast = 1, nbCellToNode
                    if (meshPairing%debug) then
                        WRITE (6, *) "  Seek for closest master cell"
                    end if
! ----------------- Get current master cell
                    cellMastNume = cellToNode(iCellMast)
                    if (meshPairing%debug) then
                        WRITE (6, *) "  Current master cell : ", cellMastNume
                    end if

! ----------------- Create master cell
                    call cellCreate(cellMastNume, nodeCoor, &
                                    meshTypeGeom, meshConx, meshConxCumu, &
                                    cellMast, cellMastLine)
                    if (meshPairing%debug) then
                        WRITE (6, *) "Properties of master cell"
                        WRITE (6, *) "========================"
                        call cellDebug(cellMast)
                    end if

! ----------------- Projection/intersection of the two cells
                    if (meshPairing%debug) then
                        WRITE (6, *) "  Compute intersection and projection in slave space"
                    end if
                    call cellInteProj(meshPairing, &
                                      cellMast, cellMastLine, cellSlavLine, &
                                      inteArea_=inteArea)
                    if (meshPairing%debug) then
                        WRITE (6, *) "  Intersection area: ", inteArea
                    end if

! ----------------- Set start elements
                    if (inteArea .gt. 100*meshPairing%pairTole) then
                        cellMastStart(1) = cellMastNume
                        nbMastStart = 1
                        cellSlavStart(1) = cellSlavNume
                        nbSlavStart = 1
                        cellSlavFlag(cellSlavIndx) = 1
                        goto 100
                    end if
                end do
            else
                if (meshPairing%debug) then
                    WRITE (6, *) "  Current slave cell is already tracked"
                end if
            end if
            cellSlavFlag(cellSlavIndx) = 2
        end do
100     continue
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! intePoinSort
!
! Sort points of intersection
!
! In  meshPairing      : main datastructure for pairing
! In  nbPoinInteIn     : number of intersection points (input)
! In  poinInteTargIn   : intersection points on target cell (input)
! In  poinInteOrigIn   : intersection points on original cell (input)
! Out nbPoinInteOut    : number of intersection points (output)
! Out poinInteTargOut  : intersection points on target cell (output)
! Out poinInteOrigOut  : intersection points on original cell (output)
!
! --------------------------------------------------------------------------------------------------
    subroutine intePoinSort(meshPairing, &
                            nbPoinInteIn, poinInteTargIn, poinInteOrigIn, &
                            nbPoinInteOut, poinInteTargOut, poinInteOrigOut)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        integer(kind=8), intent(in) :: nbPoinInteIn
        real(kind=8), intent(in) :: poinInteTargIn(2, 2*MAX_NB_INTE)
        real(kind=8), intent(in) :: poinInteOrigIn(2, 2*MAX_NB_INTE)
        integer(kind=8), intent(out) :: nbPoinInteOut
        real(kind=8), intent(out) :: poinInteTargOut(2, MAX_NB_INTE)
        real(kind=8), intent(out) :: poinInteOrigOut(2, MAX_NB_INTE)
! ----- Local
        real(kind=8) :: poinInteSort1(2, 2*MAX_NB_INTE), poinInteSort2(2, 2*MAX_NB_INTE)
        real(kind=8) :: angle(2*MAX_NB_INTE)
        real(kind=8) :: v(2), norm, bary(2)
        integer(kind=8) :: iPoinInte, angle_sorted(2*MAX_NB_INTE), listPoinNext(2*MAX_NB_INTE)
!   ------------------------------------------------------------------------------------------------
!
        bary = 0.d0
        v = 0.d0
        nbPoinInteOut = 0
        poinInteSort1 = 0.d0
        poinInteSort2 = 0.d0
        poinInteTargOut = 0.d0
        poinInteOrigOut = 0.d0

! ----- Set index of next points
        do iPoinInte = 2, nbPoinInteIn
            listPoinNext(iPoinInte-1) = iPoinInte
        end do
        listPoinNext(nbPoinInteIn) = 1

        if (meshPairing%spaceDime .eq. 3) then
! --------- Coordinates of barycenter
            do iPoinInte = 1, nbPoinInteIn
                bary(:) = bary(:)+poinInteTargIn(:, iPoinInte)/real(nbPoinInteIn)
            end do

! --------- Compute angles
            do iPoinInte = 1, nbPoinInteIn
                v(:) = poinInteTargIn(:, iPoinInte)-bary(:)
                angle(iPoinInte) = atan2(v(1), v(2))
            end do

! --------- Sort angles
            call ordr8(angle, nbPoinInteIn, angle_sorted)

! --------- Sort
            nbPoinInteOut = 0
            do iPoinInte = 1, nbPoinInteIn
                norm = sqrt((angle(angle_sorted(iPoinInte))- &
                             angle(angle_sorted(listPoinNext(iPoinInte))))**2)
                if (norm .gt. 10*meshPairing%pairTole) then
                    nbPoinInteOut = nbPoinInteOut+1
                    poinInteSort1(1:2, nbPoinInteOut) = poinInteTargIn(1:2, angle_sorted(iPoinInte))
                    poinInteSort2(1:2, nbPoinInteOut) = poinInteOrigIn(1:2, angle_sorted(iPoinInte))
                end if
            end do

        elseif (meshPairing%spaceDime .eq. 2) then
            poinInteSort1(1, 1) = poinInteTargIn(1, 1)
            poinInteSort1(1, 2) = poinInteTargIn(1, 1)
            poinInteSort2(1, 1) = poinInteOrigIn(1, 1)
            poinInteSort2(1, 2) = poinInteOrigIn(1, 1)

            do iPoinInte = 2, nbPoinInteIn
                if (poinInteTargIn(1, iPoinInte) .le. poinInteSort1(1, 1) .and. &
                    poinInteTargIn(1, iPoinInte) .ge. (-1.d0)) then
                    poinInteSort1(1, 1) = poinInteTargIn(1, iPoinInte)
                    poinInteSort2(1, 1) = poinInteOrigIn(1, iPoinInte)

                elseif (poinInteTargIn(1, iPoinInte) .ge. poinInteSort1(1, 2) .and. &
                        poinInteTargIn(1, iPoinInte) .le. (1.d0)) then
                    poinInteSort1(1, 2) = poinInteTargIn(1, iPoinInte)
                    poinInteSort2(1, 2) = poinInteOrigIn(1, iPoinInte)

                end if
            end do
            nbPoinInteOut = 2

        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Copy
        if (nbPoinInteOut .le. MAX_NB_INTE) then
            do iPoinInte = 1, nbPoinInteOut
                poinInteTargOut(1, iPoinInte) = poinInteSort1(1, iPoinInte)
                poinInteOrigOut(1, iPoinInte) = poinInteSort2(1, iPoinInte)
                if (meshPairing%spaceDime == 3) then
                    poinInteTargOut(2, iPoinInte) = poinInteSort1(2, iPoinInte)
                    poinInteOrigOut(2, iPoinInte) = poinInteSort2(2, iPoinInte)
                end if
            end do
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! intePoinInCell
!
! Check if intersection points are in cell and adjust them with tolerance
!
! In  meshPairing   : main datastructure for pairing
! In  cellGeom      : general geometric properties of cell
! In  nbPoinInte    : number of intersection points
! In  poinInteIn    : list of intersection points
! Out poinInteOut   : list of intersection points after adjustements
!
! --------------------------------------------------------------------------------------------------
    subroutine intePoinInCell(meshPairing, cellGeom, &
                              nbPoinInte, poinInteIn, poinInteOut)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        type(CELL_GEOM), intent(in) :: cellGeom
        integer(kind=8), intent(in) :: nbPoinInte
        real(kind=8), intent(in) :: poinInteIn(2, MAX_NB_INTE)
        real(kind=8), intent(out) :: poinInteOut(2, MAX_NB_INTE)
! ----- Local
        integer(kind=8) :: iPoinInte
        aster_logical :: pointIsInside
!   ------------------------------------------------------------------------------------------------
!
        poinInteOut = 0.d0
        do iPoinInte = 1, nbPoinInte
            poinInteOut(:, iPoinInte) = poinInteIn(:, iPoinInte)

!---------- Adjust point inside element
            call cellPoinAdjust(cellGeom, meshPairing%pairTole, poinInteOut(:, iPoinInte))

! --------- Test if point is inside element
            pointIsInside = &
                cellPoinInside(cellGeom, meshPairing%pairTole, poinInteOut(:, iPoinInte))
            ASSERT(pointIsInside)

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! inteCellArea
!
! Compute area of intersection
!
! In  spaceDime        : dimension of space (2 or 3)
! In  nbPoinInte       : number of intersection points
! In  poinInte         : list of intersection points
! Out inteArea         : area of intersection
!
! --------------------------------------------------------------------------------------------------
    subroutine inteCellArea(spaceDime, nbPoinInte, poinInte, &
                            inteArea)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: spaceDime, nbPoinInte
        real(kind=8), intent(in) :: poinInte(2, MAX_NB_INTE)
        real(kind=8), intent(out) :: inteArea
! ----- Local
        integer(kind=8) :: listPoinNext(MAX_NB_INTE), iPoinInte
!   ------------------------------------------------------------------------------------------------
!
        inteArea = 0.d0
        ASSERT(nbPoinInte .gt. 0)

! ----- Set index of next points
        if (spaceDime .eq. 3) then
            do iPoinInte = 2, nbPoinInte
                listPoinNext(iPoinInte-1) = iPoinInte
            end do
            listPoinNext(nbPoinInte) = 1
        end if

! ----- Compute area
        if ((nbPoinInte .gt. 2 .and. spaceDime .eq. 3) .or. &
            (nbPoinInte .ge. 2 .and. spaceDime .eq. 2)) then
            if (spaceDime .eq. 3) then
                do iPoinInte = 1, nbPoinInte
                    inteArea = inteArea+ &
                               poinInte(1, iPoinInte)* &
                               poinInte(2, listPoinNext(iPoinInte))- &
                               poinInte(1, listPoinNext(iPoinInte))* &
                               poinInte(2, iPoinInte)
                end do
                inteArea = 1.d0/2.d0*inteArea
                inteArea = sqrt(inteArea**2)
            elseif (spaceDime .eq. 2) then
                inteArea = sqrt((poinInte(1, 2)-poinInte(1, 1))**2)
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isFatalError
!
! Is the error is fatal ?
!
! --------------------------------------------------------------------------------------------------
    function isFatalError(errorPair)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: isFatalError
        integer(kind=8), intent(in) :: errorPair
!   ------------------------------------------------------------------------------------------------
!
        isFatalError = ASTER_TRUE
        if (errorPair .eq. ERR_CELL_ORTH) then
            isFatalError = ASTER_TRUE
        elseif (errorPair .eq. ERR_CELL_OOR) then
            isFatalError = ASTER_TRUE
        elseif (errorPair .eq. ERR_PAIR_PROJ) then
            isFatalError = ASTER_FALSE
        elseif (errorPair .eq. ERR_INTE_VOID) then
            isFatalError = ASTER_TRUE
        elseif (errorPair .eq. ERR_PAIR_SLAV) then
            isFatalError = ASTER_TRUE
        elseif (errorPair .eq. ERR_PAIR_MAST) then
            isFatalError = ASTER_TRUE
        end if
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! pairAllocate
!
! Allocate objects for pairing
!
! --------------------------------------------------------------------------------------------------
    subroutine pairAllocate(pairDime, meshPairing)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: pairDime
        type(MESH_PAIRING), intent(inout) :: meshPairing
!   ------------------------------------------------------------------------------------------------
!
        meshPairing%nbPairMax = pairDime
        AS_ALLOCATE(vi=meshPairing%pair, size=2*pairDime)
        AS_ALLOCATE(vi=meshPairing%nbPoinInte, size=pairDime)
        AS_ALLOCATE(vr=meshPairing%poinInte, size=2*pairDime*MAX_NB_INTE)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pairDeallocate
!
! Deallocate objects for pairing
!
! --------------------------------------------------------------------------------------------------
    subroutine pairDeallocate(meshPairing)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(inout) :: meshPairing
!   ------------------------------------------------------------------------------------------------
!
        meshPairing%nbPairMax = 0
        meshPairing%nbPair = 0
        AS_DEALLOCATE(vi=meshPairing%pair)
        AS_DEALLOCATE(vi=meshPairing%nbPoinInte)
        AS_DEALLOCATE(vr=meshPairing%poinInte)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pairResize
!
! Resize objects for pairing
!
! --------------------------------------------------------------------------------------------------
    subroutine pairResize(meshPairing, new_size)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(inout) :: meshPairing
        integer(kind=8), intent(in) :: new_size
!   ------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: nbPairOld
        integer(kind=8), pointer :: pair(:) => null()
        integer(kind=8), pointer :: nbPoinInte(:) => null()
        real(kind=8), pointer :: poinInte(:) => null()
!
        if (new_size > meshPairing%nbPairMax) then
            nbPairOld = meshPairing%nbPair
!
            AS_ALLOCATE(vi=pair, size=2*nbPairOld)
            AS_ALLOCATE(vi=nbPoinInte, size=nbPairOld)
            AS_ALLOCATE(vr=poinInte, size=2*nbPairOld*MAX_NB_INTE)
!
            pair(1:2*nbPairOld) = meshPairing%pair(1:2*nbPairOld)
            nbPoinInte(1:nbPairOld) = meshPairing%nbPoinInte(1:nbPairOld)
            poinInte(1:2*nbPairOld*MAX_NB_INTE) = meshPairing%poinInte(1:2*nbPairOld*MAX_NB_INTE)
!
            call pairDeallocate(meshPairing)
            call pairAllocate(new_size, meshPairing)
!
            meshPairing%nbPair = nbPairOld
            meshPairing%pair(1:2*nbPairOld) = pair(1:2*nbPairOld)
            meshPairing%nbPoinInte(1:nbPairOld) = nbPoinInte(1:nbPairOld)
            meshPairing%poinInte = poinInte(1:2*nbPairOld*MAX_NB_INTE)
!
            AS_DEALLOCATE(vi=pair)
            AS_DEALLOCATE(vi=nbPoinInte)
            AS_DEALLOCATE(vr=poinInte)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pairAdd
!
! Add a pair
!
! --------------------------------------------------------------------------------------------------
    subroutine pairAdd(cellSlavNume, cellMastNume, &
                       nbPoinInte, poinInte, &
                       meshPairing)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: cellSlavNume, cellMastNume
        integer(kind=8), intent(in) :: nbPoinInte
        real(kind=8), intent(in) :: poinInte(2, MAX_NB_INTE)
        type(MESH_PAIRING), intent(inout) :: meshPairing
! ----- Locals
        integer(kind=8) :: iPoinInte, iPair
!   ------------------------------------------------------------------------------------------------
!
        iPair = meshPairing%nbPair+1
        if (meshPairing%nbPairMax < iPair) then
            call pairResize(meshPairing, 2*meshPairing%nbPairMax)
        end if
        ASSERT(nbPoinInte .le. MAX_NB_INTE)
        meshPairing%nbPair = meshPairing%nbPair+1
        meshPairing%pair(2*(meshPairing%nbPair-1)+1) = cellSlavNume
        meshPairing%pair(2*(meshPairing%nbPair-1)+2) = cellMastNume
        meshPairing%nbPoinInte(meshPairing%nbPair) = nbPoinInte
        do iPoinInte = 1, MAX_NB_INTE
            meshPairing%poinInte(2*MAX_NB_INTE*(iPair-1)+iPoinInte) = &
                poinInte(1, iPoinInte)
            meshPairing%poinInte(2*MAX_NB_INTE*(iPair-1)+MAX_NB_INTE+iPoinInte) = &
                poinInte(2, iPoinInte)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! intePoinCoor
!
! Compute coordinates of points in real space
!
! In  cellGeom         : general geometric properties of cell
! In  nbPoinInte       : number of intersection points
! In  poinInte         : list of intersection points (parametric space of cellGeom)
! Out poinInteReal     : list of intersection points after transformation (global space)
!
! --------------------------------------------------------------------------------------------------
    subroutine intePoinCoor(cellGeom, nbPoinInte, poinInte, poinInteReal)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellGeom
        integer(kind=8), intent(in) :: nbPoinInte
        real(kind=8), intent(in) :: poinInte(2, MAX_NB_INTE)
        real(kind=8), intent(out) :: poinInteReal(3, MAX_NB_INTE)
! ----- Local
        integer(kind=8) :: iPoinInte
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(nbPoinInte .le. MAX_NB_INTE)
        poinInteReal = 0.d0
        do iPoinInte = 1, nbPoinInte
            call cellPoinParaToGlob(cellGeom, &
                                    poinInte(:, iPoinInte), poinInteReal(:, iPoinInte))
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! quadPoinCoor
!
! Compute coordinates of quadrature points in real space
!
! In  cellGeom         : general geometric properties of cell
! In  nbPoinQuad       : number of quadrature points
! In  quadPoinSlav     : list of quadrature points (parametric space of cellGeom)
! Out quadPoinReal     : list of quadrature points after transformation (global space)
!
! --------------------------------------------------------------------------------------------------
    subroutine quadPoinCoor(cellSlav, nbPoinQuad, quadPoinSlav, quadPoinReal)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellSlav
        integer(kind=8), intent(in) :: nbPoinQuad
        real(kind=8), intent(in) :: quadPoinSlav(2, MAX_NB_QUAD)
        real(kind=8), intent(out) :: quadPoinReal(3, MAX_NB_QUAD)
! ----- Local
        integer(kind=8) :: iPoinQuad
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(nbPoinQuad .le. MAX_NB_QUAD)
        quadPoinReal = 0.d0
        do iPoinQuad = 1, nbPoinQuad
            call cellPoinParaToGlob(cellSlav, &
                                    quadPoinSlav(:, iPoinQuad), quadPoinReal(:, iPoinQuad))
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getClosestNodesFromCell
!
! Get the closest node in a list from given cell
!
! In  cellGeom         : geometric properties of cell
! Ptr nodeCoor         : pointer to coordinates of nodes
! In  nbNode           : number of nodes
! In  listNode         : list of nodes
! Out nodeNumeClosest  : index of the closest node (index in mesh)
!
! --------------------------------------------------------------------------------------------------
    subroutine getClosestNodesFromCell(cellGeom, nodeCoor, &
                                       nbNode, listNode, &
                                       nodeNumeClosest)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), pointer :: nodeCoor(:)
        integer(kind=8), intent(in) :: nbNode, listNode(nbNode)
        integer(kind=8), intent(out) :: nodeNumeClosest
! ----- Local
        real(kind=8) :: cellCentGlob(3)
        integer(kind=8) :: iDime, iNode, nodeNume
        real(kind=8) :: vect_pm(3), distMini, dist
!   ------------------------------------------------------------------------------------------------
!
        nodeNumeClosest = 0

! ----- Get center of cell in global space
        call cellCompCenterGlob(cellGeom, cellCentGlob)

! ----- Find closest node
        distMini = 0.d0
        do iNode = 1, nbNode
            nodeNume = listNode(iNode)
            do iDime = 1, 3
                vect_pm(iDime) = nodeCoor(3*(nodeNume-1)+iDime)-cellCentGlob(iDime)
            end do
            dist = sqrt(vect_pm(1)**2+vect_pm(2)**2+vect_pm(3)**2)
            if (dist .lt. distMini .or. iNode .eq. 1) then
                distMini = dist
                nodeNumeClosest = nodeNume
            end if
        end do
        ASSERT(nodeNumeClosest .ne. 0)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getCellsFromNode
!
! Get list of cells attached to a node
!
! In  nodeNume         : index in mesh of reference node
! Ptr meshConxInve     : pointers to inverse connectivity
!     meshConxInveCumu
! In  nbCell           : length of list of cells
! In  listCell         : list of cells
! Out nbCellToNode     : length of list of cells attached to node
! Out cellToNode       : list of cells attached to node
!
! --------------------------------------------------------------------------------------------------
    subroutine getCellsFromNode(nodeNume, &
                                meshConxInve, meshConxInveCumu, &
                                nbCell, listCell, &
                                nbCellToNode, cellToNode)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nodeNume
        integer(kind=8), pointer :: meshConxInve(:), meshConxInveCumu(:)
        integer(kind=8), intent(in) :: nbCell
        integer(kind=8), intent(in) :: listCell(nbCell)
        integer(kind=8), intent(out) :: nbCellToNode
        integer(kind=8), intent(out) :: cellToNode(nbCell)
! ----- Local
        integer(kind=8) :: iCell
!   ------------------------------------------------------------------------------------------------
!
        cellToNode = 0
        nbCellToNode = meshConxInveCumu(nodeNume+1)-meshConxInveCumu(nodeNume)
        do iCell = 1, nbCellToNode
            cellToNode(iCell) = listCell(meshConxInve(meshConxInveCumu(nodeNume)-1+iCell))
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getPairJV
!
! Get pair from JEVEUX object
!
! In  mesh             : mesh
! In  baseName         : JEVEUX base name for output objects
! Ptr nodeCoor         : pointer to coordinates of nodes
! In  iPair            : index of pair
! Out cellSlav         : slave cell
! Out cellMast         : master cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getPairJV(mesh, baseName, nodeCoor, iPair, cellSlav_, cellMast_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: mesh
        character(len=24), intent(in) :: baseName
        integer(kind=8), intent(in) :: iPair
        real(kind=8), pointer :: nodeCoor(:)
        type(CELL_GEOM), optional, intent(out) :: cellSlav_, cellMast_
! ----- Local
        integer(kind=8) :: jvData
        type(CELL_GEOM) :: cellSlav, cellMast
        character(len=24) :: zonePair
        integer(kind=8), pointer :: meshTypeGeom(:) => null()
        integer(kind=8), pointer :: meshConx(:) => null(), meshConxCumu(:) => null()
        integer(kind=8) :: nbPair, cellSlavNume, cellMastNume
!   ------------------------------------------------------------------------------------------------
!
! ----- Access to mesh
        call jeveuo(mesh(1:8)//'.TYPMAIL', 'L', vi=meshTypeGeom)
        call jeveuo(mesh(1:8)//'.CONNEX', 'L', vi=meshConx)
        call jeveuo(jexatr(mesh(1:8)//'.CONNEX', 'LONCUM'), 'L', vi=meshConxCumu)

! ----- Access to pair objects
        zonePair = baseName(1:8)//".LISTPAIRS"
        call jeexin(zonePair, jvData)
        if (jvData == 0) then
            call utmess("F", "MESH4_4")
        end if
        call jeveuo(zonePair, 'L', jvData)
        call jelira(zonePair, 'LONMAX', nbPair)
        if (iPair .le. 0 .or. iPair .gt. nbPair) then
            call utmess("F", "MESH4_5")
        end if

! ----- Get current cells in pair
        cellSlavNume = zi(jvData+2*(iPair-1)-1+1)
        call cellCreate(cellSlavNume, nodeCoor, &
                        meshTypeGeom, meshConx, meshConxCumu, &
                        cellSlav)

        !call cellSetType(meshTypeGeom, cellSlavNume, cellSlav)
        cellMastNume = zi(jvData+2*(iPair-1)-1+2)
        !call cellSetType(meshTypeGeom, cellMastNume, cellMast)
        call cellCreate(cellMastNume, nodeCoor, &
                        meshTypeGeom, meshConx, meshConxCumu, &
                        cellMast)

! ----- Outputs
        if (present(cellSlav_)) then
            cellSlav_ = cellSlav
        end if
        if (present(cellMast_)) then
            cellMast_ = cellMast
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getInteJV
!
! Get integration points from JEVEUX object
!
! In  baseName         : JEVEUX base name for output objects
! In  iPair            : index of pair
! Out nbPoinInte       : number of intersection points
! Out poinInte         : coordinates of intersection points (parametric slave coordinates)
!
! --------------------------------------------------------------------------------------------------
    subroutine getInteJV(baseName, iPair, nbPoinInte, poinInte)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=24), intent(in) :: baseName
        integer(kind=8), intent(in) :: iPair
        integer(kind=8), intent(out) :: nbPoinInte
        real(kind=8), intent(out)  :: poinInte(2, MAX_NB_INTE)
! ----- Local
        integer(kind=8) :: jvData, iPoinInte
        character(len=24) :: zonePair, zoneNbPoinInte, zonePoinInte
        integer(kind=8) :: nbPair
!   ------------------------------------------------------------------------------------------------
!
        nbPoinInte = 0
        poinInte = 0.d0

! ----- Access to pair objects
        zonePair = baseName(1:8)//".LISTPAIRS"
        call jeexin(zonePair, jvData)
        if (jvData == 0) then
            call utmess("F", "MESH4_4")
        end if
        call jeveuo(zonePair, 'L', jvData)
        call jelira(zonePair, 'LONMAX', nbPair)
        if (iPair .le. 0 .or. iPair .gt. nbPair) then
            call utmess("F", "MESH4_5")
        end if

! ----- Get intersection points
        zoneNbPoinInte = baseName(1:8)//".NBPOIN"
        zonePoinInte = baseName(1:8)//".INTERSLPTS"
        call jeveuo(zoneNbPoinInte, 'L', jvData)
        nbPoinInte = zi(jvData-1+iPair)
        call jeveuo(zonePoinInte, 'L', jvData)
        do iPoinInte = 1, MAX_NB_INTE
            poinInte(1, iPoinInte) = zr(jvData-1+2*MAX_NB_INTE*(iPair-1)+iPoinInte)
            poinInte(2, iPoinInte) = zr(jvData-1+2*MAX_NB_INTE*(iPair-1)+MAX_NB_INTE+iPoinInte)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! robustPair
!
! Robust pairing on zone
!
! In  mesh             : mesh
! In  newgeo           : updated coordinates of nodes
! In  nbCellSlav       : number of slave cells
! In  nbCellMast       : number of master cells
! In  listCellSlav     : list of slave cells
! In  listCellMast     : list of master cells
! IO  meshPairing      : main datastructure for pairing
!
! --------------------------------------------------------------------------------------------------
    subroutine robustPair(mesh, newgeo, &
                          nbCellSlav, nbCellMast, &
                          listCellSlav, listCellMast, &
                          meshPairing)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: mesh
        character(len=24), intent(in) :: newgeo
        integer(kind=8), intent(in) :: nbCellSlav, nbCellMast
        integer(kind=8), intent(in) :: listCellMast(nbCellMast), listCellSlav(nbCellSlav)
        type(MESH_PAIRING), intent(inout) :: meshPairing
! ----- Local
        aster_logical :: isFatal
        integer(kind=8) :: iCellSlav, cellSlavNume, iCellMast, cellMastNume
        real(kind=8), pointer :: nodeCoor(:) => null()
        integer(kind=8), pointer :: meshTypeGeom(:) => null()
        integer(kind=8), pointer :: meshConx(:) => null(), meshConxCumu(:) => null()
        type(CELL_GEOM) :: cellSlav, cellMast
        type(CELL_GEOM) :: cellSlavLine, cellMastLine
        integer(kind=8) :: iret
        real(kind=8) :: inteArea
        integer(kind=8) :: nbPoinInte
        real(kind=8) :: poinInte(2, MAX_NB_INTE), poinInteReal(3, MAX_NB_INTE)
        real(kind=8) :: critCellArea
!   ------------------------------------------------------------------------------------------------
!
        critCellArea = 0d0
        if (meshPairing%debug) then
            write (6, *) "=================="
            write (6, *) "= Robust pairing ="
            write (6, *) "=================="
            write (6, *) " "
        end if

! ----- Access to updated geometry
        call jeveuo(newgeo(1:19)//'.VALE', 'L', vr=nodeCoor)

! ----- Access to mesh
        call jeveuo(mesh(1:8)//'.TYPMAIL', 'L', vi=meshTypeGeom)
        call jeveuo(mesh(1:8)//'.CONNEX', 'L', vi=meshConx)
        call jeveuo(jexatr(mesh(1:8)//'.CONNEX', 'LONCUM'), 'L', vi=meshConxCumu)

! ----- Protection
        if (nbCellSlav .eq. 0 .or. nbCellMast .eq. 0) then
            call utmess('F', 'MESH4_1')
        end if

! ----- Loop on slave cells
        do iCellSlav = 1, nbCellSlav
! --------- Get slave element
            cellSlavNume = listCellSlav(iCellSlav)

! --------- Create slave cell
            call cellCreate(cellSlavNume, nodeCoor, &
                            meshTypeGeom, meshConx, meshConxCumu, &
                            cellSlav, cellSlavLine)

            if (meshPairing%debug) then
                write (6, *) "Current slave element      : ", cellSlavNume
                write (6, *) " Coordinates (global frame): ", &
                    cellSlav%coorNodeGlob(1:meshPairing%spaceDime, 1:cellSlav%nbNode)
            end if

! --------- Loop on master cells
            do iCellMast = 1, nbCellMast
! ------------- Get master element
                cellMastNume = listCellMast(iCellMast)

! ------------- Create master cell
                call cellCreate(cellMastNume, nodeCoor, &
                                meshTypeGeom, meshConx, meshConxCumu, &
                                cellMast, cellMastLine)
                if (meshPairing%debug) then
                    write (6, *) "Current master element: ", cellMastNume
                    write (6, *) " Coordinates (global frame): ", &
                        cellMast%coorNodeGlob(1:meshPairing%spaceDime, 1:cellMast%nbNode)
                end if

! ------------- Compute intersection of the two cells
                inteArea = 0.d0
                nbPoinInte = 0
                poinInte = 0.d0
                poinInteReal = 0.d0
                if (meshPairing%debug) then
                    WRITE (6, *) "Compute intersection and projection in master space"
                end if
                call cellInteProjSH(meshPairing, &
                                    cellSlav, cellSlavLine, cellMastLine, &
                                    iret, &
                                    nbPoinInte, poinInte, &
                                    inteArea)
                isFatal = isFatalError(iret)
                if (.not. isFatal .and. iret .ne. ERR_PAIR_NONE) then
                    call utmess('A', 'MESH4_3')
                    inteArea = 0.d0
                    nbPoinInte = 0
                end if
                ASSERT(nbPoinInte .le. MAX_NB_INTE)
                if (meshPairing%debug) then
                    WRITE (6, *) "Intersection area: ", inteArea
                end if
                call getCritParamArea(cellSlavLine, critCellArea)
! ------------- Add pair
                if (inteArea > meshPairing%areaTole*critCellArea .and. iret == ERR_PAIR_NONE) then
                    if (meshPairing%debug) then
                        call intePoinCoor(cellSlav, nbPoinInte, poinInte, poinInteReal)
                        WRITE (6, *) "Add pair: ", meshPairing%nbPair+1, &
                            "(", cellSlavNume, "-", cellMastNume, ")"
                        WRITE (6, *) "Nb points integrations                : ", &
                            nbPoinInte
                        WRITE (6, *) "Coor. points integrations (parametric): ", &
                            poinInte(:, 1:nbPoinInte)
                        WRITE (6, *) "Coef. points integrations (global)    : ", &
                            poinInteReal(:, 1:nbPoinInte)
                        WRITE (6, *) "Area of intersection                  : ", &
                            inteArea
                    end if
                    call pairAdd(cellSlavNume, cellMastNume, &
                                 nbPoinInte, poinInte, &
                                 meshPairing)
                end if
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine

! --------------------------------------------------------------------------------------------------
!
! getCritParamArea
!
! Return the reference area of the slave cell
!
! In  meshPairing      : main datastructure for pairing
! In  cellProj         : geometric properties of projected cell
! In  nbNodeProj       : number of projected nodes
! In  cellOrigLine     : geometric properties of projected cell in origin reference frame
! In  cellTargLine     : geometric properties of cell where cell has been projected
! In  normOrig         : normal to origin cell
! In  normTarg         : normal to target cell
! Out inteIsEmpty      : flag for empty intersection
!
! --------------------------------------------------------------------------------------------------
    subroutine getCritParamArea(cellOrigLine, refCritArea)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellOrigLine
        real(kind=8), intent(out) :: refCritArea
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(cellOrigLine%isLinear)
        if (cellOrigLine%cellCode .eq. "SE2") then
            refCritArea = 2.0d0
        elseif (cellOrigLine%cellCode .eq. "TR3") then
            refCritArea = 0.5d0
        elseif (cellOrigLine%cellCode .eq. "QU4") then
            refCritArea = 4.0d0
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------------
!
! addPoinInte3DSH
!
! Compute intersection in parametric master space and return points coordinates in slave space
!    - For 3D case with Sutherland-Hodgman algorithm
!
! In  meshPairing      : main datastructure for pairing
! In  cellProj         : geometric properties of projected cell from slave cell
! In  nbNodeProj       : number of projected nodes
! In  cellOrigLine     : geometric properties of cell to project (linearized)
! In  cellTargLine     : geometric properties of cell where to project (linearized)
! Out nbPoinInte       : number of intersection points
! Out poinInteTarg     : coordinates of intersection points in target cell
! Out poinInteOrig     : coordinates of intersection points in original cell
!
! --------------------------------------------------------------------------------------------------
    subroutine addPoinInte3DSH(meshPairing, cellProj, nbNodeProj, &
                               cellOrigLine, cellTargLine, &
                               nbPoinInte, poinInteTarg, poinInteOrig)
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        type(CELL_GEOM), intent(in) :: cellProj
        integer(kind=8), intent(in) :: nbNodeProj
        type(CELL_GEOM), intent(in) :: cellOrigLine, cellTargLine
        integer(kind=8), intent(out) :: nbPoinInte
        real(kind=8), intent(out) :: poinInteTarg(2, 2*MAX_NB_INTE), poinInteOrig(2, 2*MAX_NB_INTE)
        ! Compute polygon intersection in the master plane
        call polygonIntersectionSH(cellProj%coorNodePara(1:2, :), nbNodeProj, &
                                   cellTargLine%coorNodePara(1:2, :), cellTargLine%nbNode, &
                                   poinInteTarg, nbPoinInte)
        ! Compute the intersection in the slave plane
        call computeIntePtsOrigCell(cellOrigLine, nbPoinInte, poinInteOrig, poinInteTarg, &
                                    cellProj, meshPairing, cellTargLine)
    end subroutine addPoinInte3DSH

! --------------------------------------------------------------------------------------------------
!
! polygonIntersectionSH
!
! Compute intersection between two convex polygons using Sutherland-Hodgman algorithm
!
! In  polygon1          : nodes of the first polygon (subject)
! In  n1                : number of nodes of the first polygon
! In  polygon2          : nodes of the second polygon (clip)
! In  n2                : number of nodes of the second polygon
! Out  intersection     : nodes of the intersection (output)
! Out  ninter           : number of nodes of the intersection
! --------------------------------------------------------------------------------------------------
    subroutine polygonIntersectionSH(polygon1, n1, polygon2, n2, intersection, ninter)
! ----- Parameters
        integer(kind=8), intent(in) :: n1, n2
        real(kind=8), intent(in) :: polygon1(2, n1), polygon2(2, n2)
        integer(kind=8), intent(out) :: ninter
        real(kind=8), intent(out) :: intersection(2, 2*MAX_NB_INTE)
! ----- Local
        real(kind=8) :: clipped1(2, 2*MAX_NB_INTE), clipped2(2, 2*MAX_NB_INTE)
        integer(kind=8) :: nclip1, nclip2
        !
        call sutherlandHodgmanClip(polygon1, n1, polygon2, n2, clipped1, nclip1)
        call sutherlandHodgmanClip(polygon1, n1, polygon2, n2, clipped2, nclip2)
        if (nclip1 > 0 .and. nclip2 == 0) then
            intersection(:, 1:nclip1) = clipped1(:, 1:nclip1)
            ninter = nclip1
        else if (nclip2 > 0 .and. nclip1 == 0) then
            intersection(:, 1:nclip2) = clipped2(:, 1:nclip2)
            ninter = nclip2
        else if (nclip1 > 0 .and. nclip2 > 0) then
            intersection(:, 1:nclip1) = clipped1(:, 1:nclip1)
            ninter = nclip1
        else
            ninter = 0
        end if
    end subroutine polygonIntersectionSH

! --------------------------------------------------------------------------------------------------
!
! sutherlandHodgmanClip
!
! Compute intersection between two convex polygons using Sutherland-Hodgman algorithm
!
! In  cellA             : nodes of the first polygon (subject)
! In  nA                : number of nodes of the first polygon
! In  cellB             : nodes of the second polygon (clip)
! In  nB                : number of nodes of the second polygon
! Out  intePts          : nodes of the intersection (output)
! Out  nintePts         : number of nodes of the intersection
! --------------------------------------------------------------------------------------------------
    subroutine sutherlandHodgmanClip(cellA, nA, cellB, nB, intePts, nintePts)
! ----- Parameters
        integer(kind=8), intent(in) :: nA, nB
        real(kind=8), intent(in) :: cellA(2, nA), cellB(2, nA)
        real(kind=8), intent(out) :: intePts(2, 2*MAX_NB_INTE)
        integer(kind=8), intent(out) :: nintePts
! ----- Local
        real(kind=8) :: inputPts(2, 2*MAX_NB_INTE), edgePt1(2), edgePt2(2)
        real(kind=8) :: S(2), E(2), intersection(2)
        integer(kind=8) :: i, j, ninput
        aster_logical :: inside_S, inside_E
        !
        intePts(:, 1:nA) = cellA(:, 1:nA)
        nintePts = nA

        do i = 1, nB
            edgePt2(:) = cellB(:, i)
            if (i .eq. 1) then
                edgePt1(:) = cellB(:, nB)
            else
                edgePt1(:) = cellB(:, i-1)
            end if

            ninput = nintePts

            inputPts(:, 1:ninput) = intePts(:, 1:ninput)

            nintePts = 0
            if (ninput == 0) exit

            S = inputPts(:, ninput)

            do j = 1, ninput
                E = inputPts(:, j)
                call isInside(E, edgePt1, edgePt2, inside_E)
                call isInside(S, edgePt1, edgePt2, inside_S)
                if (inside_E) then
                    if (.not. inside_S) then
                        call computeIntersectionPoint3D(S, E, edgePt1, edgePt2, intersection)
                        nintePts = nintePts+1
                        intePts(:, nintePts) = intersection
                    end if
                    nintePts = nintePts+1
                    intePts(:, nintePts) = E
                else if (inside_S) then
                    call computeIntersectionPoint3D(S, E, edgePt1, edgePt2, intersection)
                    nintePts = nintePts+1
                    intePts(:, nintePts) = intersection
                end if
                S = E
            end do
        end do
    end subroutine sutherlandHodgmanClip

! --------------------------------------------------------------------------------------------------
!
! isInside
!
!   Check if a point is inside a polygone, given one oriented edge
!
! In  point                    : point to test
! In  edgePt1                  : first point of the edge
! In  edgePt2                  : second point of the edge
! Out  inside                  : boolean, if True then the point is inside
! --------------------------------------------------------------------------------------------------
    subroutine isInside(point, edgePt1, edgePt2, inside)
! ----- Parameters
        real(kind=8), intent(in) :: point(2), edgePt1(2), edgePt2(2)
        aster_logical, intent(out) :: inside
! ----- Local
        real(kind=8) cross_product
        !
        cross_product = (edgePt2(1)-edgePt1(1))*(point(2)-edgePt1(2))- &
                        (edgePt2(2)-edgePt1(2))*(point(1)-edgePt1(1))
        inside = cross_product > 0.0
    end subroutine isInside

! --------------------------------------------------------------------------------------------------
!
! computeIntersectionPoint3D
!
!   Compute the intersection point between two lines
!
! In  point1                    : init point of the first line
! In  point2                    : end point of the first line
! In  edgePt1                   : init poinr of the second line
! In  edgePt2                   : end point of the second line
! Out  intersection             : intersection point of the two lines
!
! --------------------------------------------------------------------------------------------------
    subroutine computeIntersectionPoint3D(point1, point2, edgePt1, edgePt2, intersection)
! ----- Parameters
        real(kind=8), intent(in) :: point1(2), point2(2), edgePt1(2), edgePt2(2)
        real(kind=8), intent(out) :: intersection(2)
! ----- Local
        real(kind=8) :: dcx, dcy, dpx, dpy, n1, n2, n3
        !
        dcx = edgePt1(1)-edgePt2(1)
        dcy = edgePt1(2)-edgePt2(2)
        dpx = point1(1)-point2(1)
        dpy = point1(2)-point2(2)
        n1 = edgePt1(1)*edgePt2(2)-edgePt1(2)*edgePt2(1)
        n2 = point1(1)*point2(2)-point1(2)*point2(1)
        n3 = 1.0/(dcx*dpy-dcy*dpx)
        intersection(1) = (n1*dpx-n2*dcx)*n3
        intersection(2) = (n1*dpy-n2*dcy)*n3
    end subroutine computeIntersectionPoint3D

! --------------------------------------------------------------------------------------------------
!
! computeIntePtsOrigCell
!
!   Compute coordinates of the intersection points in the slave cell (origin)
!
! In  cellOrigLine              : geometric properties of cell to project (linearized)
! In  nbPoinInte                : number of intersection points
! In  cellProj                  : geometric properties of projected cell from slave cell
! In  meshPairing               : main datastructure for pairing
! In  cellTargLine              : geometric properties of cell where to project (linearized)
! Out  poinInteOrig             : coordinates of intersection points in original cell
! Out  poinInteTarg             : coordinates of intersection points in target cell
!
! --------------------------------------------------------------------------------------------------
    subroutine computeIntePtsOrigCell(cellOrigLine, nbPoinInte, poinInteOrig, &
                                      poinInteTarg, cellProj, meshPairing, cellTargLine)
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellOrigLine, cellProj, cellTargLine
        integer(kind=8), intent(in) :: nbPoinInte
        real(kind=8), intent(out) :: poinInteOrig(2, 2*MAX_NB_INTE), poinInteTarg(2, 2*MAX_NB_INTE)
        type(MESH_PAIRING), intent(in) :: meshPairing
! ----- Local
        real(kind=8) :: lambda_tria(3)
        real(kind=8) :: lambda_quad(4)
        integer(kind=8) :: iPtInte
        real(kind=8) :: normOrig(3), normOrigOppo(3)
        type(MESH_PROJ_PARA) :: projPara
        real(kind=8) :: poinInteGlob(3)
        aster_logical :: isOnEdge
        !
        if (cellOrigLine%cellCode .eq. "TR3") then
            do iPtInte = 1, nbPoinInte
                ! Compute barycentric coordinates
                call computeBarycentricCoordsTRIA(cellProj%coorNodePara(1:2, 1), &
                                                  cellProj%coorNodePara(1:2, 2), &
                                                  cellProj%coorNodePara(1:2, 3), &
                                                  poinInteTarg(1:2, iPtInte), &
                                                  lambda_tria)
                ! Compute the coordinates of the point
                poinInteOrig(1:2, iPtInte) = lambda_tria(1)*cellOrigLine%coorNodePara(1:2, 1)+ &
                                             lambda_tria(2)*cellOrigLine%coorNodePara(1:2, 2)+ &
                                             lambda_tria(3)*cellOrigLine%coorNodePara(1:2, 3)
            end do
        elseif (cellOrigLine%cellCode .eq. "QU4") then
            do iPtInte = 1, nbPoinInte
                call computeBaryCoordsIfOnEdgeQUAD(cellProj%coorNodePara(1:2, 1), &
                                                   cellProj%coorNodePara(1:2, 2), &
                                                   cellProj%coorNodePara(1:2, 3), &
                                                   cellProj%coorNodePara(1:2, 4), &
                                                   poinInteTarg(1:2, iPtInte), &
                                                   lambda_quad, isOnEdge)
                if (isOnEdge) then
                    poinInteOrig(1:2, iPtInte) = lambda_quad(1)*cellOrigLine%coorNodePara(1:2, 1)+ &
                                                 lambda_quad(2)*cellOrigLine%coorNodePara(1:2, 2)+ &
                                                 lambda_quad(3)*cellOrigLine%coorNodePara(1:2, 3)+ &
                                                 lambda_quad(4)*cellOrigLine%coorNodePara(1:2, 4)
                else
                    call cellPoinParaToGlob(cellTargLine, poinInteTarg(1:2, iPtInte), poinInteGlob)
                    call cellCompNormAtBary(meshPairing%spaceDime, cellOrigLine, normOrig)
                    normOrigOppo = -normOrig
                    ! ----- Set parameters of Newton algorithm for projection
                    projPara%debug = ASTER_FALSE
                    projPara%newtIterMaxi = 75
                    projPara%newtTole = meshPairing%pairTole
                    ! ----- Prepare object for projection on target cell
                    projPara%spaceDime = meshPairing%spaceDime
                    projPara%coorPoinGlob = poinInteGlob
                    projPara%projVect = normOrig
                    ! ----- Projection by given vector
                    call poinProjByVect(projPara, cellOrigLine)
                    poinInteOrig(1:2, iPtInte) = projPara%coorProjPara
                end if
            end do
        else
            ASSERT(ASTER_FALSE)
        end if

    end subroutine computeIntePtsOrigCell

! --------------------------------------------------------------------------------------------------
!
! computeBarycentricCoordsTRIA
!
!   Compute barycentric coordinates of a point in a triangle
!
! In  pointA            : first vertex of the triangle
! In  pointB            : second vertex of the triangle
! In  pointC            : third vertex of the triangle
! In  pointX            : point of interest
! Out lambda            : barycentric coordinates of pointX
!
! --------------------------------------------------------------------------------------------------
    subroutine computeBarycentricCoordsTRIA(pointA, pointB, pointC, pointX, lambda)
! ----- Parameters
        real(kind=8), intent(in) :: pointA(2), pointB(2), pointC(2), pointX(2)
        real(kind=8), intent(out) :: lambda(3)
! ----- Local
        real(kind=8) :: Area, Area_a, Area_b, Area_c
        !
        ! Compute the area of the main triangle
        Area = 0.5*abs((pointB(1)-pointA(1))*(pointC(2)-pointA(2))- &
                       (pointB(2)-pointA(2))*(pointC(1)-pointA(1)))
        ! Compute the area of subtriangles
        Area_a = 0.5*abs((pointB(1)-pointX(1))*(pointC(2)-pointX(2))- &
                         (pointB(2)-pointX(2))*(pointC(1)-pointX(1)))
        Area_b = 0.5*abs((pointC(1)-pointX(1))*(pointA(2)-pointX(2))- &
                         (pointC(2)-pointX(2))*(pointA(1)-pointX(1)))
        Area_c = 0.5*abs((pointA(1)-pointX(1))*(pointB(2)-pointX(2))- &
                         (pointA(2)-pointX(2))*(pointB(1)-pointX(1)))
        lambda(1) = Area_a/Area
        lambda(2) = Area_b/Area
        lambda(3) = Area_c/Area
    end subroutine computeBarycentricCoordsTRIA

! --------------------------------------------------------------------------------------------------
!
! computeBaryCoordsIfOnEdgeQUAD
!
!   Compute the barycentric coordinates in quadrangle if point is on one edge
!
! In  pointA            : first vertex of the quadrangle
! In  pointB            : second vertex of the quadrangle
! In  pointC            : third vertex of the quadrangle
! In  pointD            : third vertex of the quadrangle
! In  pointX            : point of interest
! Out lambda_quad       : barycentric coordinates
! Out isOnEdge          : boolean, if True, point is on one edge
!
! --------------------------------------------------------------------------------------------------
    subroutine computeBaryCoordsIfOnEdgeQUAD(pointA, pointB, pointC, pointD, pointX, &
                                             lambda_quad, isOnEdge)
! ----- Parameters
        real(kind=8), intent(in) :: pointA(2), pointB(2), pointC(2), pointD(2), pointX(2)
        real(kind=8), intent(out) :: lambda_quad(4)
        aster_logical :: isOnEdge
! ----- Local
        real(kind=8) :: A1, A2, A3, A4
        real(kind=8) :: tolArea
        real(kind=8) :: lambda_tri(3)
        aster_logical :: inside_tria_1, inside_tria_2
        !
        ! Compute signed areas
        call computeSignedAreaTRIA(pointX, pointA, pointB, A1)
        call computeSignedAreaTRIA(pointX, pointB, pointC, A2)
        call computeSignedAreaTRIA(pointX, pointC, pointD, A3)
        call computeSignedAreaTRIA(pointX, pointD, pointA, A4)
        ! Compute areas
        A1 = abs(A1)
        A2 = abs(A2)
        A3 = abs(A3)
        A4 = abs(A4)
        ! Test
        tolArea = 1e-12
        call isInsideTRIA(pointX, pointA, pointB, pointC, inside_tria_1)
        call isInsideTRIA(pointX, pointA, pointC, pointD, inside_tria_2)
        if ((A1 <= tolArea) .or. (A2 <= tolArea)) then
            call isInsideTRIA(pointX, pointA, pointB, pointC, inside_tria_1)
            if (inside_tria_1) then
                call computeBarycentricCoordsTRIA(pointA, pointB, pointC, pointX, lambda_tri)
                lambda_quad(1) = lambda_tri(1)
                lambda_quad(2) = lambda_tri(2)
                lambda_quad(3) = lambda_tri(3)
                lambda_quad(4) = 0.0
                isOnEdge = ASTER_TRUE
            else
                lambda_quad = -1.0
                isOnEdge = ASTER_FALSE
            end if
        else if ((A3 <= tolArea) .or. (A4 <= tolArea)) then
            call isInsideTRIA(pointX, pointA, pointC, pointD, inside_tria_2)
            if (inside_tria_2) then
                call computeBarycentricCoordsTRIA(pointA, pointC, pointD, pointX, lambda_tri)
                lambda_quad(1) = lambda_tri(1)
                lambda_quad(2) = 0.0
                lambda_quad(3) = lambda_tri(2)
                lambda_quad(4) = lambda_tri(3)
                isOnEdge = ASTER_TRUE
            else
                lambda_quad = -1.0
                isOnEdge = ASTER_FALSE
            end if
        else
            lambda_quad = -1.0
            isOnEdge = ASTER_FALSE
        end if

    end subroutine computeBaryCoordsIfOnEdgeQUAD

! --------------------------------------------------------------------------------------------------
!
! computeSignedAreaTRIA
!
!   Compute signed area of a triangle
!
! In  pointA            : first vertex of the triangle
! In  pointB            : second vertex of the triangle
! In  pointC            : third vertex of the triangle
! Out signArea          : signed area of the triangle
!
! --------------------------------------------------------------------------------------------------
    subroutine computeSignedAreaTRIA(pointA, pointB, pointC, signArea)
! ----- Parameters
        real(kind=8), intent(in) :: pointA(2), pointB(2), pointC(2)
        real(kind=8), intent(out) :: signArea
        !
        signArea = 0.5*((pointB(1)-pointA(1))*(pointC(2)-pointA(2)) &
                        -(pointB(2)-pointA(2))*(pointC(1)-pointA(1)))
    end subroutine computeSignedAreaTRIA

! --------------------------------------------------------------------------------------------------
!
! isInsideTRIA
!
!   Check if a point is inside a triangle
!
! In  pointA            : first vertex of the triangle
! In  pointB            : second vertex of the triangle
! In  pointC            : third vertex of the triangle
! In  pointX            : point of interest
! Out inside            : boolean, if True the point is inside
!
! --------------------------------------------------------------------------------------------------
    subroutine isInsideTRIA(pointX, pointA, pointB, pointC, inside)
! ----- Parameters
        real(kind=8), intent(in) :: pointX(2), pointA(2), pointB(2), pointC(2)
        aster_logical, intent(out) :: inside
! ----- Local
        real(kind=8) A1, A2, A3, prod, tole
        !
        tole = 1e-12
        call computeSignedAreaTRIA(pointX, pointA, pointB, A1)
        call computeSignedAreaTRIA(pointX, pointB, pointC, A2)
        call computeSignedAreaTRIA(pointX, pointC, pointA, A3)
        prod = abs(A1*A2*A3)
        inside = (A1 >= 0.0 .and. A2 >= 0.0 .and. A3 >= 0.0) .or. &
                 (A1 <= 0.0 .and. A2 <= 0.0 .and. A3 <= 0.0) .or. (prod <= tole)
    end subroutine isInsideTRIA
! --------------------------------------------------------------------------------------------------
!
! cellInteProjSH - prjint_ray
!
! Compute intersection and projection from cell to target cell (only used for robust for now)
!
! In  meshPairing      : main datastructure for pairing
! In  cellOrig         : general geometric properties of origin cell
! In  cellOrigLine     : general geometric properties of linearized origin cell
! In  cellTargLine     : general geometric properties of linearized target cell
! Out iret             : return code error
! Out nbPoinInte       : number of intersection points
! Out poinInteOrig     : coordinates of intersection points (in origin cell parametric space)
!
! --------------------------------------------------------------------------------------------------
    subroutine cellInteProjSH(meshPairing, &
                              cellOrig, cellOrigLine, cellTargLine, &
                              iret_, &
                              nbPoinInte_, poinInteOrig_, &
                              inteArea_, inteNeigh_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_PAIRING), intent(in) :: meshPairing
        type(CELL_GEOM), intent(in) :: cellOrig, cellOrigLine, cellTargLine
        integer(kind=8), optional, intent(out) :: iret_
        integer(kind=8), optional, intent(out) :: nbPoinInte_
        real(kind=8), optional, intent(out) :: poinInteOrig_(2, MAX_NB_INTE)
        real(kind=8), optional, intent(out) :: inteArea_
        integer(kind=8), optional, intent(out) :: inteNeigh_(MAX_NB_NEIGH)
! ----- Local
        aster_logical, parameter :: debug = ASTER_FALSE
        integer(kind=8) :: iret, inteNeigh(MAX_NB_NEIGH)
        integer(kind=8) :: nbPoinInteI
        real(kind=8) :: poinInteTargI(2, 2*MAX_NB_INTE), poinInteOrigI(2, 2*MAX_NB_INTE)
        integer(kind=8) :: nbPoinInteS
        real(kind=8) :: poinInteTargS(2, MAX_NB_INTE), poinInteOrigS(2, MAX_NB_INTE)
        integer(kind=8) :: nbPoinInte
        real(kind=8) :: poinInteOrig(2, MAX_NB_INTE)
        real(kind=8) :: inteArea
        aster_logical :: errorProj, errorInte
        type(CELL_GEOM) :: cellProj
        integer(kind=8) :: nbCmpPara, nbNodeProj
!   ------------------------------------------------------------------------------------------------
!
        inteNeigh = 0
        nbPoinInte = 0
        poinInteOrig = 0.d0
        inteArea = 0.d0
        iret = ERR_PAIR_NONE
        nbCmpPara = meshPairing%spaceDime-1

! ----- Compute projection of cell on target cell
        errorProj = ASTER_FALSE
        call cellProjOnCell(meshPairing, &
                            cellOrig, cellOrigLine, cellTargLine, &
                            cellProj, nbNodeProj, iret)
        if (iret .ne. ERR_PAIR_NONE) then
            errorProj = ASTER_TRUE
            goto 100
        end if
        if (meshPairing%debug) then
            WRITE (6, *) "     Intersection is OK"
        end if

! ----- Compute intersection in parametric space of target cell
        nbPoinInteI = 0
        poinInteTargI = 0.d0
        poinInteOrigI = 0.d0
        if (meshPairing%spaceDime .eq. 2) then
            ASSERT(nbNodeProj .eq. 2)
            call addPoinInte2D(meshPairing, cellProj, &
                               nbPoinInteI, poinInteTargI, poinInteOrigI, &
                               inteNeigh)
        elseif (meshPairing%spaceDime .eq. 3) then
            call addPoinInte3DSH(meshPairing, cellProj, nbNodeProj, &
                                 cellOrigLine, cellTargLine, &
                                 nbPoinInteI, poinInteTargI, poinInteOrigI)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Trop de points d'intersection (avant tri) => on ne sait pas gérer
        ASSERT(nbPoinInteI .le. 2*MAX_NB_INTE)

! ----- Error management
        errorInte = ASTER_FALSE
        if (nbPoinInteI == 0 .or. iret == ERR_PAIR_MAST) then
            errorInte = ASTER_TRUE
            goto 100
        end if

! ----- Debug print
        if (meshPairing%debug) then
            WRITE (6, *) "     Intersection (before re-ordering): ", nbPoinInteI, &
                " points of intersection"
            WRITE (6, *) "     Target side : ", poinInteTargI(1:nbCmpPara, 1:nbPoinInteI)
            WRITE (6, *) "     Origin side : ", poinInteOrigI(1:nbCmpPara, 1:nbPoinInteI)
        end if

! ----- Sort list of intersection points
        if ((nbPoinInteI .gt. 2 .and. meshPairing%spaceDime == 3) .or. &
            (nbPoinInteI .ge. 2 .and. meshPairing%spaceDime == 2)) then
            call intePoinSort(meshPairing, &
                              nbPoinInteI, poinInteTargI, poinInteOrigI, &
                              nbPoinInteS, poinInteTargS, poinInteOrigS)
        else
            ASSERT(nbPoinInteI .le. MAX_NB_INTE)
            nbPoinInteS = nbPoinInteI
            poinInteTargS(:, 1:nbPoinInteI) = poinInteTargI(:, 1:nbPoinInteI)
            poinInteOrigS(:, 1:nbPoinInteI) = poinInteOrigI(:, 1:nbPoinInteI)
        end if

! ----- Error management
        if (nbPoinInteS > MAX_NB_INTE) then
            if (meshPairing%debug) then
                WRITE (6, *) "     Intersection (after re-ordering): ", nbPoinInteS, &
                    " points of intersection"
                WRITE (6, *) "     Too many points !"
            end if
            iret = ERR_PAIR_SLAV
            go to 99
        end if

! ----- Debug print
        if (meshPairing%debug) then
            WRITE (6, *) "     Intersection (after re-ordering): ", nbPoinInteS, &
                " points of intersection"
            WRITE (6, *) "     Target side : ", poinInteTargS(1:nbCmpPara, 1:nbPoinInteS)
            WRITE (6, *) "     Origin side : ", poinInteOrigS(1:nbCmpPara, 1:nbPoinInteS)
        end if

! ----- All nodes have to be inside original cell
        call intePoinInCell(meshPairing, cellOrigLine, &
                            nbPoinInteS, poinInteOrigS, poinInteOrig)
        nbPoinInte = nbPoinInteS

! ----- Compute weight of intersection
        inteArea = 0.d0
        if (nbPoinInte .ne. 0) then
            call inteCellArea(meshPairing%spaceDime, nbPoinInte, poinInteOrig, inteArea)
            if (meshPairing%debug) then
                WRITE (6, *) "     Area of intersection : ", inteArea
            end if
        end if

! ----- Error
100     continue
        if (errorProj .or. errorInte) then
            if (meshPairing%debug) then
                WRITE (6, *) "     Failure of intersection"
            end if
            nbPoinInte = 0
            poinInteOrig = 0.d0
            inteNeigh = 0
            inteArea = 0.d0
        else
            if (meshPairing%debug) then
                WRITE (6, *) "     Success of intersection"
            end if
        end if

! ----- Fatal error only in debug mode
99      continue
        if (debug) then
            ASSERT(iret .ne. ERR_PAIR_NONE)
        end if

! ----- Outputs
        if (present(inteNeigh_)) then
            inteNeigh_ = inteNeigh
        end if
        if (present(nbPoinInte_)) then
            nbPoinInte_ = nbPoinInte
        end if
        if (present(inteArea_)) then
            inteArea_ = inteArea
        end if
        if (present(poinInteOrig_)) then
            poinInteOrig_ = poinInteOrig
        end if
        if (present(iret_)) then
            iret_ = iret
        end if
    end subroutine cellInteProjSH
!
!   ------------------------------------------------------------------------------------------------
end module
