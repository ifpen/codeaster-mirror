/**
 * @file BaseMeshInterface.cxx
 * @brief Interface python de BaseMesh
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "PythonBindings/BaseMeshInterface.h"

#include "aster_pybind.h"

#include <Meshes/BaseMesh.h>

void exportBaseMeshToPython( py::module_ &mod ) {

    py::class_< BaseMesh, BaseMesh::BaseMeshPtr, DataStructure >( mod, "BaseMesh" )
        // fake initFactoryPtr: created by subclass
        // fake initFactoryPtr: created by subclass
        .def( "build", &BaseMesh::build, R"(
Build list of Tables based on the mesh

Returns:
    bool: true if building is ok
        )" )
        .def( "getNumberOfNodes", &BaseMesh::getNumberOfNodes, R"(
Return the number of nodes of the mesh.

Returns:
    int: Number of nodes.
        )" )
        .def( "getNumberOfCells", &BaseMesh::getNumberOfCells, R"(
Return the number of cells of the mesh.

Returns:
    int: Number of cells.
        )" )
        .def( "getCoordinates", &BaseMesh::getCoordinates, R"(
Return the coordinates of the mesh.

Returns:
    MeshCoordinatesField: Field of the coordinates.
        )" )
        .def( "isConnection", &BaseMesh::isConnection, R"(
Function to know if a mesh is a ConnectionMesh)" )
        .def( "isIncomplete", &BaseMesh::isIncomplete, R"(
Tell if the mesh is complete on parallel instances.

Returns:
    bool: *False* for a centralized or parallel mesh, *True* for an incomplete mesh.
        )" )
        .def( "isParallel", &BaseMesh::isParallel, R"(
Tell if the mesh is distributed on parallel instances.

Returns:
    bool: *False* for a centralized mesh, *True* for a parallel mesh.
        )" )
        .def( "getDimension", &BaseMesh::getDimension, R"(
Return the dimension of the mesh.

Returns:
    int: 2 or 3
        )" )
        .def( "getConnectivity", &BaseMesh::getConnectivityZeroBased, R"(
Return the connectivity of the mesh as Python lists.

Returns:
    list[list[int]]: List of, for each cell, a list of the nodes indexes.
        )" )
        .def( "getCellTypeName", &BaseMesh::getCellTypeName, R"(
Return the type name of the given cell

Arguments:
    index (int) : index of the cell (0-based)

Returns:
    str : name of the cell type (stripped)
        )",
              py::arg( "index" ) )
        .def( "getCellType", &BaseMesh::getCellType, R"(
Return the type of the given cell

Arguments:
    index (int) : index of the cell (0-based)

Returns:
    int : the cell type
        )",
              py::arg( "index" ) )
        .def( "hasCellsOfType", &BaseMesh::hasCellsOfType, R"(
Return True if mesh contains at least one cell of given type

Arguments:
    type (str) : cell type

Returns:
    bool : *True* if mesh contains at least one cell of given type, else *False*
        )",
              py::arg( "type" ) )
        .def( "getMedConnectivity", &BaseMesh::getMedConnectivityZeroBased, R"(
Return the connectivity of the mesh as Python lists following the Med IDs.

Returns:
    list[list[int]]: List of, for each cell, a list of the nodes indexes.
        )" )
        .def( "getMedCellsTypes", &BaseMesh::getMedCellsTypes, R"(
Return the Med type of each cell.

Returns:
    list[int]: List of Med types.
        )" )
        .def( "getCellsType", &BaseMesh::getCellsType, R"(
Return the Aster type of each cell.

Returns:
    list[int]: List of Aster types.
        )" )
        .def( "updateInternalState", &BaseMesh::update_tables, R"(
Update the internal state of the datastructure.

Returns:
    bool: *True* in case of success, *False* otherwise.
        )" )
        .def( "getTable", &BaseMesh::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              py::arg( "identifier" ) )
        .def( "printMedFile", &BaseMesh::printMedFile, R"(
Print the mesh in the MED format

Arguments:
    filename (Path|str): Name of the file
    local (bool=True) : print local values only (relevant for a ParallelMesh only)
    version (list): list of size 3 ([major, minor, release])
    mode (str): 'a' or 'w' ('w' to overwrite file)

Returns:
    Bool: True if of
            )",
              py::arg( "fileName" ), py::arg( "local" ) = true,
              py::arg( "version" ) = std::array< int, 3 >( { 0, 0, 0 } ), py::arg( "mode" ) = 'a' )

        /* Mesh builder functions */
        .def( "_initDefinition", &BaseMesh::initDefinition, R"(
Initialize the mesh creation (*for internal use only*).

Arguments:
    dimension (int): Dimension of the mesh.
    coordinates (list[float]): Nodes coordinates (size: 3 x number of nodes).
    connectivity (list[list[int]]): Connectivity of each cells (size: number of cells).
    types (list[int]): Types of each cells (size: number of cells).
    nbGrpCells (int): Number of groups of cells to be defined.
    nbGrpNodes (int): Number of groups of nodes to be defined.
        )",
              py::arg( "dimension" ), py::arg( "coordinates" ), py::arg( "connectivity" ),
              py::arg( "types" ), py::arg( "nbGrpCells" ), py::arg( "nbGrpNodes" ) )

        .def( "_endDefinition", &BaseMesh::endDefinition, R"(
Terminate the mesh creation (*for internal use only*).
        )" )

        .def( "_addGroupsOfNodes", &BaseMesh::addGroupsOfNodes, R"(
Add groups of nodes (*for internal use only*).

Arguments:
    names (str): Groups names.
    groupsOfNodes (list[list[int]]): List of nodes indexes for each group.
        )",
              py::arg( "names" ), py::arg( "groupsOfNodes" ) )

        .def( "_addGroupsOfCells", &BaseMesh::addGroupsOfCells, R"(
Add groups of cells (*for internal use only*).

Arguments:
    name (str): Groups names.
    groupsOfCells (list[list[int]]): List of cells indexes for each group.
        )",
              py::arg( "names" ), py::arg( "groupsOfCells" ) )

        .def( "show", &BaseMesh::show, R"(
Show mesh informations.

Arguments:
    verbosity (int): Verbosity level (default: 1)
        )",
              py::arg( "verbosity" ) = 1 )

        .def( "check", &BaseMesh::check, R"(
Check some properties of the mesh.

Arguments:
    tolerance (float): Tolerance used to detect flat cells.
        )",
              py::arg( "tolerance" ) )
        .def( "getLocalToGlobalNodeIds", &BaseMesh::getLocalToGlobalNodeIds,
              R"(
Returns local to global node Ids mapping

Returns:
    list[int]: local to global IDs mapping.
        )" )
        .def( "getLocalToGlobalCellIds", &BaseMesh::getLocalToGlobalCellIds,
              R"(
Returns local to global IDs mapping for cells

Returns:
    list[int]: local to global IDs mapping.
        )" )
        .def( "getRestrictedToOriginalNodesIds", &BaseMesh::getRestrictedToOriginalNodesIds,
              R"(
If the mesh is created as a restriction of an initial mesh,
It returns for each nodes, the node id of the initial mesh.

Returns:
    list[int]: for each nodes, the node id of the initial mesh.
        )" )
        .def( "getRestrictedToOriginalCellsIds", &BaseMesh::getRestrictedToOriginalCellsIds,
              R"(
If the mesh is created as restriction of an initial mesh,
It returns for each cells, the cell id of the initial mesh.

Returns:
    list[int]: for each cells, the cell id of the initial mesh.
        )" )
        .def( "getOriginalToRestrictedNodesIds", &BaseMesh::getOriginalToRestrictedNodesIds,
              R"(
If the mesh is created as a restriction of an initial mesh,
It returns a dict betweenn the node id of the initial mesh and the current node id.

Returns:
    dict[int]: a dict betweenn the node id of the initial mesh and the current node id.
        )" )
        .def( "getOriginalToRestrictedCellsIds", &BaseMesh::getOriginalToRestrictedCellsIds,
              R"(
If the mesh is created as restriction of an initial mesh,
It returns a dict between the cell id of the initial mesh and the current cell id.

Returns:
    dict[int]: a dict between the cell id of the initial mesh and the current cell id.
        )" )
        .def( "getMinMaxEdgeSizes", &BaseMesh::getMinMaxEdgeSizes, R"(
Get minimum and maximum length of edges in group of cells

Returns:
    tuple(real): values of min and max edges
        )" );
};
