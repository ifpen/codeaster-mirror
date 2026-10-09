/**
 * @file ContactComputationInterface.cxx
 * @brief Interface python de ContactComputation
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

#include "PythonBindings/ContactComputationInterface.h"

#include "aster_pybind.h"

// aslint: disable=C3006

void exportContactComputationToPython( py::module_ &mod ) {

    py::class_< ContactComputation, ContactComputationPtr >( mod, "ContactComputation" )
        .def( py::init( &initFactoryPtr< ContactComputation, ContactNewPtr > ) )
        .def( define_pickling< ContactComputation >() )
        .def( "contactData", &ContactComputation::contactData, R"(
Compute contact data as input to compute contact forces and matrices.

Arguments:
    pairing (ContactPairing): pairing object
    material (MaterialField): material field
    initial_contact (bool): True to use value in contact definition (CONTACT_INIT).

Returns:
    FieldOnCellsReal: contact data
        )",
              py::arg( "pairing" ), py::arg( "material" ), py::arg( "initial_contact" ) )
        .def( "hasVirtualCells", &ContactComputation::hasVirtualCells, R"(
Detect if virtual cells has been created

Arguments:
    pairing (ContactPairing): pairing object

Returns:
    bool: surfaces have been paired
        )",
              py::arg( "pairing" ) )
        .def( "setVerbosity", &ContactComputation::setVerbosity, R"(
Set level of verbosity
0- without
1- normal (default)
2- detailled

Arguments:
    level (int) : level of verbosity
)",
              py::arg( "level" ) )
        .def( "getVerbosity", &ContactComputation::getVerbosity, R"(
Get level of verbosity
0- without
1- normal
2- detailled

Returns:
    int: level of verbosity
)" )
        .def( "contactCoefficient", &ContactComputation::contactCoefficient, R"(
Compute contact coefficients at the nodes of the slave surface based on values of COEF_CONT
and COEF_FROT

Returns:
    list[FieldOnNodesReal]: coefficients (COEF_CONT and COEF_FROT)
        )" );
};
