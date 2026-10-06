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

#include "PythonBindings/ElementaryMatrixInterface.h"

#include "aster_pybind.h"

void exportElementaryMatrixToPython( py::module_ &mod ) {

    py::class_< BaseElementaryMatrix, BaseElementaryMatrixPtr, DataStructure >(
        mod, "BaseElementaryMatrix" )
        // fake initFactoryPtr: not buildable
        // fake initFactoryPtr: not buildable
        .def( "getModel", &BaseElementaryMatrix::getModel )
        .def( "setModel", &BaseElementaryMatrix::setModel );

    py::class_< ElementaryMatrixDisplacementReal,
                ElementaryMatrixDisplacementReal::ElementaryMatrixPtr, BaseElementaryMatrix >(
        mod, "ElementaryMatrixDisplacementReal" )
        .def( py::init( &initFactoryPtr< ElementaryMatrixDisplacementReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryMatrixDisplacementReal, std::string > ) )
        .def(
            py::init( &initFactoryPtr< ElementaryMatrixDisplacementReal, ModelPtr, std::string > ) )
        .def( py::self *= float() )
        .def( "build", &ElementaryMatrixDisplacementReal::build )
        .def( "getMesh", &ElementaryMatrixDisplacementReal::getMesh )
        .def( "getFiniteElementDescriptors",
              &ElementaryMatrixDisplacementReal::getFiniteElementDescriptors )
        .def( "addElementaryTerm", py::overload_cast< const ElementaryTermRealPtr & >(
                                       &ElementaryMatrixDisplacementReal::addElementaryTerm ) )
        .def( "addElementaryTerm",
              py::overload_cast< const std::vector< ElementaryTermRealPtr > & >(
                  &ElementaryMatrixDisplacementReal::addElementaryTerm ) )
        .def( "getElementaryTerms", &ElementaryMatrixDisplacementReal::getElementaryTerms )
        .def( "getNumberOfElementaryTerms",
              &ElementaryMatrixDisplacementReal::getNumberOfElementaryTerms )
        .def( "hasElementaryTerms", &ElementaryMatrixDisplacementReal::hasElementaryTerms );

    py::class_< ElementaryMatrixDisplacementComplex,
                ElementaryMatrixDisplacementComplex::ElementaryMatrixPtr, BaseElementaryMatrix >(
        mod, "ElementaryMatrixDisplacementComplex" )
        .def( py::init( &initFactoryPtr< ElementaryMatrixDisplacementComplex > ) )
        .def( py::init( &initFactoryPtr< ElementaryMatrixDisplacementComplex, std::string > ) )
        .def( py::init(
            &initFactoryPtr< ElementaryMatrixDisplacementComplex, ModelPtr, std::string > ) )
        .def( "getMesh", &ElementaryMatrixDisplacementComplex::getMesh )
        .def( "build", &ElementaryMatrixDisplacementComplex::build )
        .def( "getFiniteElementDescriptors",
              &ElementaryMatrixDisplacementComplex::getFiniteElementDescriptors )
        .def( "addElementaryTerm", py::overload_cast< const ElementaryTermComplexPtr & >(
                                       &ElementaryMatrixDisplacementComplex::addElementaryTerm ) )
        .def( "addElementaryTerm",
              py::overload_cast< const std::vector< ElementaryTermComplexPtr > & >(
                  &ElementaryMatrixDisplacementComplex::addElementaryTerm ) )
        .def( "getElementaryTerms", &ElementaryMatrixDisplacementComplex::getElementaryTerms )
        .def( "getNumberOfElementaryTerms",
              &ElementaryMatrixDisplacementComplex::getNumberOfElementaryTerms )
        .def( "hasElementaryTerms", &ElementaryMatrixDisplacementComplex::hasElementaryTerms );

    py::class_< ElementaryMatrixTemperatureReal,
                ElementaryMatrixTemperatureReal::ElementaryMatrixPtr, BaseElementaryMatrix >(
        mod, "ElementaryMatrixTemperatureReal" )
        .def( py::init( &initFactoryPtr< ElementaryMatrixTemperatureReal > ) )
        .def( py::init( &initFactoryPtr< ElementaryMatrixTemperatureReal, std::string > ) )
        .def(
            py::init( &initFactoryPtr< ElementaryMatrixTemperatureReal, ModelPtr, std::string > ) )
        .def( py::self *= float() )
        .def( "getMesh", &ElementaryMatrixTemperatureReal::getMesh )
        .def( "build", &ElementaryMatrixTemperatureReal::build )
        .def( "getFiniteElementDescriptors",
              &ElementaryMatrixTemperatureReal::getFiniteElementDescriptors )
        .def( "addElementaryTerm", py::overload_cast< const ElementaryTermRealPtr & >(
                                       &ElementaryMatrixTemperatureReal::addElementaryTerm ) )
        .def( "addElementaryTerm",
              py::overload_cast< const std::vector< ElementaryTermRealPtr > & >(
                  &ElementaryMatrixTemperatureReal::addElementaryTerm ) )
        .def( "getElementaryTerms", &ElementaryMatrixTemperatureReal::getElementaryTerms )
        .def( "getNumberOfElementaryTerms",
              &ElementaryMatrixTemperatureReal::getNumberOfElementaryTerms )
        .def( "hasElementaryTerms", &ElementaryMatrixTemperatureReal::hasElementaryTerms );

    py::class_< ElementaryMatrixPressureComplex,
                ElementaryMatrixPressureComplex::ElementaryMatrixPtr, BaseElementaryMatrix >(
        mod, "ElementaryMatrixPressureComplex" )
        .def( py::init( &initFactoryPtr< ElementaryMatrixPressureComplex > ) )
        .def( py::init( &initFactoryPtr< ElementaryMatrixPressureComplex, std::string > ) )
        .def(
            py::init( &initFactoryPtr< ElementaryMatrixPressureComplex, ModelPtr, std::string > ) )
        .def( "getMesh", &ElementaryMatrixPressureComplex::getMesh )
        .def( "build", &ElementaryMatrixPressureComplex::build )
        .def( "getFiniteElementDescriptors",
              &ElementaryMatrixPressureComplex::getFiniteElementDescriptors )
        .def( "addElementaryTerm", py::overload_cast< const ElementaryTermComplexPtr & >(
                                       &ElementaryMatrixPressureComplex::addElementaryTerm ) )
        .def( "addElementaryTerm",
              py::overload_cast< const std::vector< ElementaryTermComplexPtr > & >(
                  &ElementaryMatrixPressureComplex::addElementaryTerm ) )
        .def( "getElementaryTerms", &ElementaryMatrixPressureComplex::getElementaryTerms )
        .def( "getNumberOfElementaryTerms",
              &ElementaryMatrixPressureComplex::getNumberOfElementaryTerms )
        .def( "hasElementaryTerms", &ElementaryMatrixPressureComplex::hasElementaryTerms );
};
