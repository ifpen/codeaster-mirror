/**
 * @file ElementaryMatrixInterface.cxx
 * @brief Interface python de ElementaryMatrix
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

#include "PythonBindings/ElementaryTermInterface.h"

#include "aster_pybind.h"

void exportElementaryTermToPython( py::module_ &mod ) {

    py::class_< ElementaryTermReal, ElementaryTermRealPtr, DataField >( mod, "ElementaryTermReal" )
        .def( py::init( &initFactoryPtr< ElementaryTermReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryTermReal, std::string > ) )
        .def( "getFiniteElementDescriptor", &ElementaryTermReal::getFiniteElementDescriptor, R"(
            Return the finite element descriptor

            Returns:
                FiniteElementDescriptor: finite element descriptor
            )" )
        .def( "getOption", &ElementaryTermReal::getOption, R"(
            Return the optior used to compute it

            Returns:
                str: name of the option
            )" )
        .def( "getMesh", &ElementaryTermReal::getMesh, R"(
            Return the mesh

            Returns:
                BaseMesh: a pointer to the mesh
        )" )
        .def( "getLocatedComponent", &ElementaryTermReal::getLocatedComponent, R"(
            Return the local mode.

            Returns:
                str: the local mode
        )" )
        .def( "getValues", py::overload_cast<>( &ElementaryTermReal::getValues, py::const_ ), R"(
            Return the values of the field.

            Returns:
                list[list[float]]: values
        )" )
        .def( "getPhysicalQuantity", &ElementaryTermReal::getPhysicalQuantity, R"(
            Return the physical quantity

            Returns:
                str: name of the physical quantity
            )" )
        .def( "exists", &ElementaryTermReal::exists, R"(
            Return True if it exists

            Returns:
                bool: True if exist
            )" );

    py::class_< ElementaryTermComplex, ElementaryTermComplexPtr, DataField >(
        mod, "ElementaryTermComplex" )
        .def( py::init( &initFactoryPtr< ElementaryTermComplex > ) )
        .def( py::init( &initFactoryPtr< ElementaryTermComplex, std::string > ) )
        .def( "getFiniteElementDescriptor", &ElementaryTermComplex::getFiniteElementDescriptor, R"(
            Return the finite element descriptor

            Returns:
                FiniteElementDescriptor: finite element descriptor
            )" )
        .def( "getOption", &ElementaryTermComplex::getOption, R"(
            Return the optior used to compute it

            Returns:
                str: name of the option
            )" )
        .def( "getMesh", &ElementaryTermComplex::getMesh, R"(
            Return the mesh

            Returns:
                BaseMesh: a pointer to the mesh
            )" )
        .def( "getLocatedComponent", &ElementaryTermComplex::getLocatedComponent, R"(
            Return the local mode.

            Returns:
                str: the local mode
        )" )
        .def( "getPhysicalQuantity", &ElementaryTermComplex::getPhysicalQuantity, R"(
            Return the physical quantity

            Returns:
                str: name of the physical quantity
            )" )
        .def( "exists", &ElementaryTermComplex::exists, R"(
            Return True if it exists

            Returns:
                bool: True if exist
            )" );
};
