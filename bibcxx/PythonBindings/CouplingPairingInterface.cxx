/**
 * @file CouplingPairingInterface.cxx
 * @brief Interface python de CouplingPairing
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

#include "PythonBindings/CouplingPairingInterface.h"

#include "aster_pybind.h"
// aslint: disable=C3006

void exportCouplingPairingToPython( py::module_ &mod ) {

    py::class_< CouplingPairing, CouplingPairingPtr, DataStructure > class_( mod, "CouplingPairing",
                                                                             R"(
Object to create contact pairing.)" );
    class_.def( py::init( &initFactoryPtr< CouplingPairing, ModelPtr, ASTERINTEGER > ) );
    class_.def( "getMesh", &CouplingPairing::getMesh, R"(
Mesh

Returns:
    BaseMesh: the mesh
)" );
    class_.def( "compute", &CouplingPairing::compute,
                R"(
Compute the pairing quantitie

Returns:
    bool: True if the pairing quantities are updated appropriately
)" );
    class_.def( "getNumberOfPairs", &CouplingPairing::getNumberOfPairs, R"(
Return number of pairs

Returns:
    int: number of pairs
)" );
    class_.def( "getFiniteElementDescriptor", &CouplingPairing::getFiniteElementDescriptor, R"(
Return Finite Element Descriptor for virtual cells from pairing.

Returns:
    FiniteElementDescriptor: finite element for virtual cells
)" );

    class_.def( "getPairingField", &CouplingPairing::getPairingField, R"(
Get intersection points

Returns:
    FieldOnCellsReal: intersection points.
)" );
    class_.def( "getHHOProjector", &CouplingPairing::getHHOProjector, R"(
Get HHO L2-projector

Argument:
    pairing [FieldOnCellsReal]: pairing.

Returns:
    FieldOnCellsReal: L2-projector.
)",
                py::arg( "pairing" ) );
    class_.def( "addZone", &CouplingPairing::addZone, R"(
Add a new zone of coupling;

Argument:
    zone [CouplingZonePairing]: zone.
)",
                py::arg( "zone" ) );
};
