/**
 * @file AssemblyMatrixInterface.cxx
 * @brief Interface python de AssemblyMatrix
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

#include "PythonBindings/BaseAssemblyMatrixInterface.h"

#include "aster_pybind.h"

void exportBaseAssemblyMatrixToPython( py::module_ &mod ) {

    py::class_< BaseAssemblyMatrix, BaseAssemblyMatrixPtr, DataStructure >( mod,
                                                                            "BaseAssemblyMatrix" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< BaseAssemblyMatrix, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< BaseAssemblyMatrix, std::string, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< BaseAssemblyMatrix, PhysicalProblemPtr, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "getDOFNumbering", &BaseAssemblyMatrix::getDOFNumbering )
        // -----------------------------------------------------------------------------------------
        .def( "getMesh", &BaseAssemblyMatrix::getMesh, R"(
Return the mesh.

Returns:
    Mesh: a pointer to the mesh
        )" )

        // -----------------------------------------------------------------------------------------
        .def( "isBuilt", &BaseAssemblyMatrix::isBuilt, R"(
Tell if the matrix has already been built.

Returns:
    bool: *True* if the matrix has been built.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "isFactorized", &BaseAssemblyMatrix::isFactorized,
              R"(
Tell if the matrix is factorized.

Returns:
    bool: *True* if the matrix is factorized, *False* otherwise.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "setDOFNumbering", &BaseAssemblyMatrix::setDOFNumbering )
        // -----------------------------------------------------------------------------------------
        .def( "updateDOFNumbering", &BaseAssemblyMatrix::updateDOFNumbering )
        // -----------------------------------------------------------------------------------------
        .def( "setSolverName", &BaseAssemblyMatrix::setSolverName )
        // -----------------------------------------------------------------------------------------
        .def( "hasDirichletEliminationDOFs", &BaseAssemblyMatrix::hasDirichletEliminationDOFs, R"(
Tell if matrix has some DOFs eliminated by Dirichlet boundaries conditions.

Returns:
    bool: *True* if matrix has some DOFs eliminated by Dirichlet boundaries conditions else *False*
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "getDirichletBCDOFs", &BaseAssemblyMatrix::getDirichletBCDOFs, R"(
Return a vector with DOFs eliminated by Dirichlet boundaries conditions (if it exists)

Returns:
    tuple(int): a vector with DOFs eliminated by Dirichlet boundaries conditions of
        size = neq + 1,
        tuple(ieq = 0, neq - 1) = 1 then DOF eliminated else 0,
        tuple(neq) = number of DOFs eliminated.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "zeroDirichletBCDOFs", &BaseAssemblyMatrix::zeroDirichletBCDOFs, R"(
Zero the Dirichlet BC DOFs in the given field

Arguments:
    field (FieldOnNodesReal): the field to modify
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "getLagrangeScaling", &BaseAssemblyMatrix::getLagrangeScaling, R"(
Return the scaling used for Lagrange multipliers. It returns 1 if no Lagrange.

Returns:
    float: scaling used for Lagrange multipliers. It returns 1 if no Lagrange
    are present.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "print", &BaseAssemblyMatrix::print,
              R"(
            Print the matrix in code_aster or matlab format (with information on the DOF).

            Arguments:
                format (str): 'ASTER' (default) or 'MATLAB'
            )",
              py::arg( "format" ) = "ASTER", py::arg( "unit" ) = 6 )
        // -----------------------------------------------------------------------------------------
        .def( "symmetrize", &BaseAssemblyMatrix::symmetrize,
              R"(
            Make the assembly matrix symmetric in place
            )" )
        // -----------------------------------------------------------------------------------------
        .def( "isSymmetric", &BaseAssemblyMatrix::isSymmetric,
              R"(
            Return True if matrix is symmetric
            )" )
        // -----------------------------------------------------------------------------------------
        .def( "getCalculOption", &BaseAssemblyMatrix::getCalculOption,
              R"(
            Return the option of CALCUL.

            Returns:
                str: Name of the option.
            )" )
        // -----------------------------------------------------------------------------------------
        .def( "transpose", &BaseAssemblyMatrix::transpose,
              R"(
            Return the transpose of the matrix in place.
            )" );
};
