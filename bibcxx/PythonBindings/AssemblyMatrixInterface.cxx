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

#include "PythonBindings/AssemblyMatrixInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/FieldOnNodesInterface.h"

void exportAssemblyMatrixToPython( py::module_ &mod ) {

    py::class_< AssemblyMatrixDisplacementReal, AssemblyMatrixDisplacementRealPtr,
                BaseAssemblyMatrix >( mod, "AssemblyMatrixDisplacementReal" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementReal > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementReal, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementReal, PhysicalProblemPtr > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementReal,
                                         const AssemblyMatrixDisplacementReal & > ) )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast< const AssemblyMatrixDisplacementReal::ElementaryMatrixPtr &,
                                 const ListOfLoadsPtr & >(
                  &AssemblyMatrixDisplacementReal::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixDisplacementReal::ElementaryMatrixPtr > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixDisplacementReal::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< std::pair< AssemblyMatrixDisplacementReal::ElementaryMatrixPtr,
                                                ASTERDOUBLE > > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixDisplacementReal::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast< const AssemblyMatrixDisplacementReal::ElementaryMatrixPtr &,
                                 const DirichletBCPtr & >(
                  &AssemblyMatrixDisplacementReal::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixDisplacementReal::ElementaryMatrixPtr > &,
                  const DirichletBCPtr & >( &AssemblyMatrixDisplacementReal::assemble ),
              R"(
                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixDisplacementReal::ElementaryMatrixPtr > &,
                  const std::vector< DirichletBCPtr > & >(
                  &AssemblyMatrixDisplacementReal::assemble ),
              R"(
                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    dirichlet (list[DirichletBC]) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< std::pair< AssemblyMatrixDisplacementReal::ElementaryMatrixPtr,
                                                ASTERDOUBLE > > &,
                  const DirichletBCPtr & >( &AssemblyMatrixDisplacementReal::assemble ),
              R"(
                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "applyDirichletBC", &AssemblyMatrixDisplacementReal::applyDirichletBC, R"(
Apply the DirichletBC into the Rhs (aka kinematic aka no Lagrange multipliers).

Arguments:
    DirichletBC [FieldOnNodes] the values on the DirichletBC.
    Rhs [FieldOnNodes] The residual to be modified.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixDisplacementReal::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "scale", &AssemblyMatrixDisplacementReal::scale, R"(
Scale the matrix in place using right and left multiplication by diagonal matrices stored as vectors

Arguments:
    lvect (list[float]): List of the values.
    rvect (list[float]): List of the values.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "defineSolver", &AssemblyMatrixDisplacementReal::defineSolver )
        // -----------------------------------------------------------------------------------------
        .def( "size", &AssemblyMatrixDisplacementReal::size, R"(
Get the size of the matrix

Arguments:
    local (bool) local or global size
        )",
              py::arg( "local" ) = true )
        // -----------------------------------------------------------------------------------------
        .def( "copy", &AssemblyMatrixDisplacementReal::copy )
        // -----------------------------------------------------------------------------------------
        .def( "getUpperValues", &AssemblyMatrixDisplacementReal::getUpperValues )
        .def( "getLowerValues", &AssemblyMatrixDisplacementReal::getLowerValues )
        // -----------------------------------------------------------------------------------------
        .def(
            "__call__",
            []( const AssemblyMatrixDisplacementReal &self, const ASTERINTEGER i_row,
                const ASTERINTEGER i_col ) { return self( i_row, i_col ); },
            R"(
Return one matrix coefficient.
Only supports Morse storage. A coefficient absent from the profile is zero.
UpdateValuePointers() has to be called first.

Arguments:
    i_row [int] Local row DOF index, zero-based.
    j_col [int] Local column DOF index, zero-based.

Returns:
    [float]: value of the matritx at (i_row, j_col)
        )",
            py::arg( "i_row" ), py::arg( "i_col" ) )
        // -----------------------------------------------------------------------------------------
        .def( "updateValuePointers", &AssemblyMatrixDisplacementReal::updateValuePointers, R"(
Update pointers.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "extractSubMatrix", &AssemblyMatrixDisplacementReal::extractSubMatrix, R"(
Extract a sub-matrix.

Arguments:
    fed [FiniteElementDescriptor] extract sub-matrix for element of fed.
    option [str] Option to known which dofs to extract.
    parameter [str] parameter of the option to known which dofs to extract.

Returns:
    [ElementaryMatrixReal]: elementary matrix extracted.
        )",
              py::arg( "fed" ), py::arg( "option" ), py::arg( "parameter" ) )
        // -----------------------------------------------------------------------------------------
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( py::self -= py::self )
        .def( py::self += py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( -py::self )
        .def(
            "__mul__", +[]( const AssemblyMatrixDisplacementReal &M, const FieldOnNodesReal &v ) {
                return M * v;
            } );

    py::class_< AssemblyMatrixEliminatedReal, AssemblyMatrixEliminatedRealPtr,
                AssemblyMatrixDisplacementReal >( mod, "AssemblyMatrixEliminatedReal" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixEliminatedReal > ) )
        .def( py::init( &initFactoryPtr< AssemblyMatrixEliminatedReal, std::string > ) );

    py::class_< AssemblyMatrixDisplacementComplex, AssemblyMatrixDisplacementComplexPtr,
                BaseAssemblyMatrix >( mod, "AssemblyMatrixDisplacementComplex" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementComplex > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixDisplacementComplex, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast< const AssemblyMatrixDisplacementComplex::ElementaryMatrixPtr &,
                                 const ListOfLoadsPtr & >(
                  &AssemblyMatrixDisplacementComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixDisplacementComplex::ElementaryMatrixPtr > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixDisplacementComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< std::pair<
                      AssemblyMatrixDisplacementComplex::ElementaryMatrixPtr, ASTERDOUBLE > > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixDisplacementComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast< const AssemblyMatrixDisplacementComplex::ElementaryMatrixPtr &,
                                 const DirichletBCPtr & >(
                  &AssemblyMatrixDisplacementComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixDisplacementComplex::ElementaryMatrixPtr > &,
                  const DirichletBCPtr & >( &AssemblyMatrixDisplacementComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixDisplacementComplex::ElementaryMatrixPtr > &,
                  const std::vector< DirichletBCPtr > & >(
                  &AssemblyMatrixDisplacementComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    dirichlet (list[DirichletBC]) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< std::pair<
                      AssemblyMatrixDisplacementComplex::ElementaryMatrixPtr, ASTERDOUBLE > > &,
                  const DirichletBCPtr & >( &AssemblyMatrixDisplacementComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "transposeConjugate", &AssemblyMatrixDisplacementComplex::transposeConjugate )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixDisplacementComplex::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "defineSolver", &AssemblyMatrixDisplacementComplex::defineSolver )
        // -----------------------------------------------------------------------------------------
        .def( "copy", &AssemblyMatrixDisplacementComplex::copy )
        // -----------------------------------------------------------------------------------------
        .def( "getUpperValues", &AssemblyMatrixDisplacementComplex::getUpperValues )
        .def( "getLowerValues", &AssemblyMatrixDisplacementComplex::getLowerValues )
        // -----------------------------------------------------------------------------------------
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( -py::self )
        .def(
            "__mul__", +[]( const AssemblyMatrixDisplacementComplex &M,
                            const FieldOnNodesComplex &v ) { return M * v; } )
        // -----------------------------------------------------------------------------------------
        .def(
            "__call__",
            []( const AssemblyMatrixDisplacementComplex &self, const ASTERINTEGER i_row,
                const ASTERINTEGER i_col ) { return self( i_row, i_col ); },
            R"(
Return one matrix coefficient.
Only supports Morse storage. A coefficient absent from the profile is zero.
UpdateValuePointers() has to be called first.

Arguments:
    i_row [int] Local row DOF index, zero-based.
    j_col [int] Local column DOF index, zero-based.

Returns:
    [cmplx]: value of the matritx at (i_row, j_col)
        )",
            py::arg( "i_row" ), py::arg( "i_col" ) )
        // -----------------------------------------------------------------------------------------
        .def( "updateValuePointers", &AssemblyMatrixDisplacementComplex::updateValuePointers, R"(
Update pointers.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "extractSubMatrix", &AssemblyMatrixDisplacementComplex::extractSubMatrix, R"(
Extract a sub-matrix.

Arguments:
    fed [FiniteElementDescriptor] extract sub-matrix for element of fed.
    option [str] Option to known which dofs to extract.
    parameter [str] parameter of the option to known which dofs to extract.

Returns:
    [ElementaryMatrixComplex]: elementary matrix extracted.
        )",
              py::arg( "fed" ), py::arg( "option" ), py::arg( "parameter" ) );
    // -----------------------------------------------------------------------------------------

    py::class_< AssemblyMatrixTemperatureReal, AssemblyMatrixTemperatureRealPtr,
                BaseAssemblyMatrix >( mod, "AssemblyMatrixTemperatureReal" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixTemperatureReal > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixTemperatureReal, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixTemperatureReal, PhysicalProblemPtr > ) )
        // -----------------------------------------------------------------------------------------
        .def(
            "assemble",
            py::overload_cast< const AssemblyMatrixTemperatureReal::ElementaryMatrixPtr &,
                               const ListOfLoadsPtr & >( &AssemblyMatrixTemperatureReal::assemble ),
            R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
            py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixTemperatureReal::ElementaryMatrixPtr > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixTemperatureReal::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< std::pair< AssemblyMatrixTemperatureReal::ElementaryMatrixPtr,
                                                ASTERDOUBLE > > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixTemperatureReal::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        //-----------------------------------------------------------------------------------------
        .def(
            "assemble",
            py::overload_cast< const AssemblyMatrixTemperatureReal::ElementaryMatrixPtr &,
                               const DirichletBCPtr & >( &AssemblyMatrixTemperatureReal::assemble ),
            R"(
               Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
            py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixTemperatureReal::ElementaryMatrixPtr > &,
                  const DirichletBCPtr & >( &AssemblyMatrixTemperatureReal::assemble ),
              R"(
                 Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def(
            "assemble",
            py::overload_cast<
                const std::vector< AssemblyMatrixTemperatureReal::ElementaryMatrixPtr > &,
                const std::vector< DirichletBCPtr > & >( &AssemblyMatrixTemperatureReal::assemble ),
            R"(
                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    dirichlet (list[DirichletBC]) : dirichlet BC to impose.
            )",
            py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< std::pair< AssemblyMatrixTemperatureReal::ElementaryMatrixPtr,
                                                ASTERDOUBLE > > &,
                  const DirichletBCPtr & >( &AssemblyMatrixTemperatureReal::assemble ),
              R"(
                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "applyDirichletBC", &AssemblyMatrixTemperatureReal::applyDirichletBC, R"(
Apply the DirichletBC into the Rhs (aka kinematic aka no Lagrange multipliers).

Arguments:
    DirichletBC [FieldOnNodes] the values on the DirichletBC.
    Rhs [FieldOnNodes] The residual to be modified.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixTemperatureReal::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "scale", &AssemblyMatrixTemperatureReal::scale, R"(
Scale the matrix in place using right and left multiplication by diagonal matrices stored as vectors

Arguments:
    lvect (list[float]): List of the values.
    rvect (list[float]): List of the values.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "defineSolver", &AssemblyMatrixTemperatureReal::defineSolver )
        // -----------------------------------------------------------------------------------------
        .def( "size", &AssemblyMatrixTemperatureReal::size, R"(
Get the size of the matrix

Arguments:
    local (bool) local or global size
        )",
              py::arg( "local" ) = true )
        // -----------------------------------------------------------------------------------------
        .def( "copy", &AssemblyMatrixTemperatureReal::copy )
        // -----------------------------------------------------------------------------------------
        .def( "getUpperValues", &AssemblyMatrixTemperatureReal::getUpperValues )
        .def( "getLowerValues", &AssemblyMatrixTemperatureReal::getLowerValues )
        // -----------------------------------------------------------------------------------------
        .def(
            "__call__",
            []( const AssemblyMatrixTemperatureReal &self, const ASTERINTEGER i_row,
                const ASTERINTEGER i_col ) { return self( i_row, i_col ); },
            R"(
Return one matrix coefficient.
Only supports Morse storage. A coefficient absent from the profile is zero.
UpdateValuePointers() has to be called first.

Arguments:
    i_row [int] Local row DOF index, zero-based.
    j_col [int] Local column DOF index, zero-based.

Returns:
    [float]: value of the matritx at (i_row, j_col)
        )",
            py::arg( "i_row" ), py::arg( "i_col" ) )
        // -----------------------------------------------------------------------------------------
        .def( "updateValuePointers", &AssemblyMatrixTemperatureReal::updateValuePointers, R"(
Update pointers.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "extractSubMatrix", &AssemblyMatrixTemperatureReal::extractSubMatrix, R"(
Extract a sub-matrix.

Arguments:
    fed [FiniteElementDescriptor] extract sub-matrix for element of fed.
    option [str] Option to known which dofs to extract.
    parameter [str] parameter of the option to known which dofs to extract.

Returns:
    [ElementaryMatrixReal]: elementary matrix extracted.
        )",
              py::arg( "fed" ), py::arg( "option" ), py::arg( "parameter" ) )
        // -----------------------------------------------------------------------------------------
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( py::self -= py::self )
        .def( py::self += py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( -py::self )
        .def(
            "__mul__", +[]( const AssemblyMatrixTemperatureReal &M, const FieldOnNodesReal &v ) {
                return M * v;
            } );
    // -----------------------------------------------------------------------------------------

    py::class_< AssemblyMatrixTemperatureComplex, AssemblyMatrixTemperatureComplexPtr,
                BaseAssemblyMatrix >( mod, "AssemblyMatrixTemperatureComplex" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixTemperatureComplex > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixTemperatureComplex, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast< const AssemblyMatrixTemperatureComplex::ElementaryMatrixPtr &,
                                 const ListOfLoadsPtr & >(
                  &AssemblyMatrixTemperatureComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixTemperatureComplex::ElementaryMatrixPtr > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixTemperatureComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices added.

                Arguments:
                    clean (bool) : Clean elementary matrices after building (default = true)
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        //-----------------------------------------------------------------------------------------
        .def(
            "assemble",
            py::overload_cast<
                const std::vector< std::pair< AssemblyMatrixTemperatureComplex::ElementaryMatrixPtr,
                                              ASTERDOUBLE > > &,
                const ListOfLoadsPtr & >( &AssemblyMatrixTemperatureComplex::assemble ),
            R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
            py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast< const AssemblyMatrixTemperatureComplex::ElementaryMatrixPtr &,
                                 const DirichletBCPtr & >(
                  &AssemblyMatrixTemperatureComplex::assemble ),
              R"(
               Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixTemperatureComplex::ElementaryMatrixPtr > &,
                  const DirichletBCPtr & >( &AssemblyMatrixTemperatureComplex::assemble ),
              R"(
                 Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        //-----------------------------------------------------------------------------------------
        .def(
            "assemble",
            py::overload_cast<
                const std::vector< std::pair< AssemblyMatrixTemperatureComplex::ElementaryMatrixPtr,
                                              ASTERDOUBLE > > &,
                const DirichletBCPtr & >( &AssemblyMatrixTemperatureComplex::assemble ),
            R"(
                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
            py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "transposeConjugate", &AssemblyMatrixTemperatureComplex::transposeConjugate );

    py::class_< AssemblyMatrixPressureReal, AssemblyMatrixPressureRealPtr, BaseAssemblyMatrix >(
        mod, "AssemblyMatrixPressureReal" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixPressureReal > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixPressureReal, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast< const AssemblyMatrixPressureReal::ElementaryMatrixPtr &,
                                 const ListOfLoadsPtr & >( &AssemblyMatrixPressureReal::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixPressureReal::ElementaryMatrixPtr > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixPressureReal::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector<
                      std::pair< AssemblyMatrixPressureReal::ElementaryMatrixPtr, ASTERDOUBLE > > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixPressureReal::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast< const AssemblyMatrixPressureReal::ElementaryMatrixPtr &,
                                 const DirichletBCPtr & >( &AssemblyMatrixPressureReal::assemble ),
              R"(
               Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixPressureReal::ElementaryMatrixPtr > &,
                  const DirichletBCPtr & >( &AssemblyMatrixPressureReal::assemble ),
              R"(
                 Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector<
                      std::pair< AssemblyMatrixPressureReal::ElementaryMatrixPtr, ASTERDOUBLE > > &,
                  const DirichletBCPtr & >( &AssemblyMatrixPressureReal::assemble ),
              R"(
                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "applyDirichletBC", &AssemblyMatrixPressureReal::applyDirichletBC, R"(
Apply the DirichletBC into the Rhs (aka kinematic aka no Lagrange multipliers).

Arguments:
    DirichletBC [FieldOnNodes] the values on the DirichletBC.
    Rhs [FieldOnNodes] The residual to be modified.
        )" )
        // -----------------------------------------------------------------------------------------
        .def( "setValues", &AssemblyMatrixPressureReal::setValues, R"(
Erase the assembly matrix and set new values in it.

The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
indices are sumed according to an assembly process.

Arguments:
    idx (list[int]): List of the row indices.
    jdx (list[int]): List of the column indices.
    values (list[float]): List of the values.
        )" )
        .def( "copy", &AssemblyMatrixPressureReal::copy )
        // -----------------------------------------------------------------------------------------
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( py::self -= py::self )
        .def( py::self += py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( -py::self )
        .def(
            "__mul__", +[]( const AssemblyMatrixPressureReal &M, const FieldOnNodesReal &v ) {
                return M * v;
            } );
    // -----------------------------------------------------------------------------------------

    py::class_< AssemblyMatrixPressureComplex, AssemblyMatrixPressureComplexPtr,
                BaseAssemblyMatrix >( mod, "AssemblyMatrixPressureComplex" )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixPressureComplex > ) )
        // -----------------------------------------------------------------------------------------
        .def( py::init( &initFactoryPtr< AssemblyMatrixPressureComplex, std::string > ) )
        // -----------------------------------------------------------------------------------------
        .def(
            "assemble",
            py::overload_cast< const AssemblyMatrixPressureComplex::ElementaryMatrixPtr &,
                               const ListOfLoadsPtr & >( &AssemblyMatrixPressureComplex::assemble ),
            R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
            py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        // -----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixPressureComplex::ElementaryMatrixPtr > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixPressureComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< std::pair< AssemblyMatrixPressureComplex::ElementaryMatrixPtr,
                                                ASTERDOUBLE > > &,
                  const ListOfLoadsPtr & >( &AssemblyMatrixPressureComplex::assemble ),
              R"(
                Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    listOfLoads (ListOfLoads) : list of loads to assemble
            )",
              py::arg( "elemMatrix" ), py::arg( "listOfLoads" ) = nullptr )
        //-----------------------------------------------------------------------------------------
        .def(
            "assemble",
            py::overload_cast< const AssemblyMatrixPressureComplex::ElementaryMatrixPtr &,
                               const DirichletBCPtr & >( &AssemblyMatrixPressureComplex::assemble ),
            R"(
               Assembly matrix from elementar matrices and list of loads.

                Arguments:
                    elemMatrix (ElementaryMatrixReal) : elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
            py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< AssemblyMatrixPressureComplex::ElementaryMatrixPtr > &,
                  const DirichletBCPtr & >( &AssemblyMatrixPressureComplex::assemble ),
              R"(
                 Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        //-----------------------------------------------------------------------------------------
        .def( "assemble",
              py::overload_cast<
                  const std::vector< std::pair< AssemblyMatrixPressureComplex::ElementaryMatrixPtr,
                                                ASTERDOUBLE > > &,
                  const DirichletBCPtr & >( &AssemblyMatrixPressureComplex::assemble ),
              R"(
                Arguments:
                    elemMatrix (list[ElementaryMatrixReal, float]) : list of pair composed of an
                    elementary matrix and the multiplicatif coefficent to assemble.
                    dirichlet (DirichletBC) : dirichlet BC to impose.
            )",
              py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        .def(
            "assemble",
            py::overload_cast<
                const std::vector< AssemblyMatrixPressureComplex::ElementaryMatrixPtr > &,
                const std::vector< DirichletBCPtr > & >( &AssemblyMatrixPressureComplex::assemble ),
            R"(
                Arguments:
                    elemMatrix (list[ElementaryMatrixReal]) : list elementary matrix to assemble.
                    dirichlet (list[DirichletBC]) : dirichlet BC to impose.
            )",
            py::arg( "elemMatrix" ), py::arg( "dirichlet" ) )
        // -----------------------------------------------------------------------------------------
        .def( "transposeConjugate", &AssemblyMatrixPressureComplex::transposeConjugate )
        // -----------------------------------------------------------------------------------------
        .def( "defineSolver", &AssemblyMatrixPressureComplex::defineSolver )
        // -----------------------------------------------------------------------------------------
        .def( "copy", &AssemblyMatrixPressureComplex::copy )
        // -----------------------------------------------------------------------------------------
        .def( "getUpperValues", &AssemblyMatrixPressureComplex::getUpperValues )
        .def( "getLowerValues", &AssemblyMatrixPressureComplex::getLowerValues )
        // -----------------------------------------------------------------------------------------
        .def( float() * py::self )
        .def( py::self *= float() )
        .def( py::self -= py::self )
        .def( py::self += py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( -py::self )
        .def(
            "__mul__", +[]( const AssemblyMatrixPressureComplex &M, const FieldOnNodesComplex &v ) {
                return M * v;
            } );
};
