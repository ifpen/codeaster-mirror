/**
 * @file MechanicalLoadInterface.cxx
 * @brief Interface python de MechanicalLoad
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

#include "PythonBindings/MechanicalLoadInterface.h"

#include "aster_pybind.h"

void exportMechanicalLoadToPython( py::module_ &mod ) {

    py::class_< MechanicalLoadReal, MechanicalLoadReal::MechanicalLoadPtr, DSWithCppPickling,
                DataStructure >( mod, "MechanicalLoadReal" )
        .def( py::init( &initFactoryPtr< MechanicalLoadReal, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< MechanicalLoadReal, std::string, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< MechanicalLoadReal, const py::tuple & > ) )
        .def( define_pickling< MechanicalLoadReal >() )
        .def( "getFiniteElementDescriptor", &MechanicalLoadReal::getFiniteElementDescriptor )
        .def( "hasLoadField", &MechanicalLoadReal::hasLoadField )
        .def( "updateValuePointers", &MechanicalLoadReal::updateValuePointers )
        .def( "getMechanicalLoadDescription", &MechanicalLoadReal::getMechanicalLoadDescription )
        .def( "getModel", &MechanicalLoadReal::getModel )
        .def( "getMesh", &MechanicalLoadReal::getMesh )
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              py::arg( "identifier" ) );

    py::class_< MechanicalLoadFunction, MechanicalLoadFunction::MechanicalLoadPtr,
                DSWithCppPickling, DataStructure >( mod, "MechanicalLoadFunction" )
        .def( py::init( &initFactoryPtr< MechanicalLoadFunction, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< MechanicalLoadFunction, std::string, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< MechanicalLoadFunction, const py::tuple & > ) )
        .def( define_pickling< MechanicalLoadFunction >() )
        .def( "getFiniteElementDescriptor", &MechanicalLoadFunction::getFiniteElementDescriptor )
        .def( "hasLoadField", &MechanicalLoadFunction::hasLoadField )
        .def( "updateValuePointers", &MechanicalLoadFunction::updateValuePointers )
        .def( "getModel", &MechanicalLoadFunction::getModel )
        .def( "getMesh", &MechanicalLoadFunction::getMesh )
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              py::arg( "identifier" ) );

    py::class_< MechanicalLoadComplex, MechanicalLoadComplex::MechanicalLoadPtr, DSWithCppPickling,
                DataStructure >( mod, "MechanicalLoadComplex" )
        .def( py::init( &initFactoryPtr< MechanicalLoadComplex, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< MechanicalLoadComplex, std::string, ModelPtr > ) )
        .def( py::init( &initFactoryPtr< MechanicalLoadComplex, const py::tuple & > ) )
        .def( define_pickling< MechanicalLoadComplex >() )
        .def( "getFiniteElementDescriptor", &MechanicalLoadComplex::getFiniteElementDescriptor )
        .def( "hasLoadField", &MechanicalLoadComplex::hasLoadField )
        .def( "updateValuePointers", &MechanicalLoadComplex::updateValuePointers )
        .def( "getModel", &MechanicalLoadComplex::getModel )
        .def( "getMesh", &MechanicalLoadComplex::getMesh )
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              py::arg( "identifier" ) );
};
