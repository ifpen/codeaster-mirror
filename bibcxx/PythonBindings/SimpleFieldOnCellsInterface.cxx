/**
 * @file SimpleFieldOnCellsInterface.cxx
 * @brief Interface python de SimpleFieldOnCells
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

#include "PythonBindings/SimpleFieldOnCellsInterface.h"

#include "aster_pybind.h"

#include "DataFields/FieldConverter.h"
#include "PythonBindings/DataStructureInterface.h"

void exportSimpleFieldOnCellsToPython( py::module_ &mod ) {
    py::class_< SimpleFieldOnCellsReal, SimpleFieldOnCellsRealPtr, DataField >(
        mod, "SimpleFieldOnCellsReal" )
        .def( py::init( &initFactoryPtr< SimpleFieldOnCellsReal, std::string > ) )
        .def( py::init( &initFactoryPtr< SimpleFieldOnCellsReal, BaseMeshPtr > ),
              py::arg( "mesh" ) )
        .def( py::init( &initFactoryPtr< SimpleFieldOnCellsReal, BaseMeshPtr, std::string,
                                         std::string, VectorString > ),
              py::arg( "mesh" ), py::arg( "loc" ), py::arg( "quantity" ), py::arg( "cmps" ) )
        .def( py::init( &initFactoryPtr< SimpleFieldOnCellsReal, BaseMeshPtr, std::string,
                                         std::string, VectorString, bool > ),
              py::arg( "mesh" ), py::arg( "loc" ), py::arg( "quantity" ), py::arg( "cmps" ),
              py::arg( "prol_zero" ) )
        .def( py::init( &initFactoryPtr< SimpleFieldOnCellsReal, BaseMeshPtr, std::string,
                                         std::string, VectorString, ASTERINTEGER, ASTERINTEGER > ),
              py::arg( "mesh" ), py::arg( "loc" ), py::arg( "quantity" ), py::arg( "cmps" ),
              py::arg( "nbPoints" ), py::arg( "nbSubPoints" ) )
        .def( py::init(
                  &initFactoryPtr< SimpleFieldOnCellsReal, BaseMeshPtr, std::string, std::string,
                                   VectorString, ASTERINTEGER, ASTERINTEGER, bool > ),
              py::arg( "mesh" ), py::arg( "loc" ), py::arg( "quantity" ), py::arg( "cmps" ),
              py::arg( "nbPoints" ), py::arg( "nbSubPoints" ), py::arg( "prol_zero" ) )
        .def( py::init(
                  &initFactoryPtr< SimpleFieldOnCellsReal, BaseMeshPtr, std::string, std::string,
                                   VectorString, const VectorInt &, ASTERINTEGER, bool > ),
              py::arg( "mesh" ), py::arg( "loc" ), py::arg( "quantity" ), py::arg( "cmps" ),
              py::arg( "nbPoints" ), py::arg( "nbSubPoints" ), py::arg( "prol_zero" ) )
        .def(
            "__getitem__",
            +[]( const SimpleFieldOnCellsReal &v, const VectorLong &i ) {
                if ( i.size() == 3 ) {
                    return v( i[0], i[1], i[2] );
                }

                return v( i[0], i[1], i[2], i[3] );
            } )
        .def(
            "__setitem__",
            +[]( SimpleFieldOnCellsReal &v, const VectorLong &i, ASTERDOUBLE f ) {
                if ( i.size() == 3 ) {
                    return v( i[0], i[1], i[2] ) = f;
                }

                return v( i[0], i[1], i[2], i[3] ) = f;
            } )
        .def( "setMesh", &SimpleFieldOnCellsReal::setMesh,
              R"(
            Set mesh.

            Arguments:
                mesh [BaseMesh]: mesh to set.
            )",
              py::arg( "mesh" ) )
        .def( "getValue", &SimpleFieldOnCellsReal::getValue, py::return_value_policy::copy, R"(
Returns the value of the `icmp` component of the field on the `ima` cell,
at the `ipt` point, at the `ispt` sub-point.

Args:
    ima  (int): Index of cells.
    icmp (int): Index of component.
    ipt  (int): Index of point.
    ispt (int): Index of sub-point (default = 0).

Returns:
    float: Value of field at *ima*, of *icmp*, at *ipt*, at *ispt*;
             NaN if the position is not allocated.
        )",
              py::arg( "ima" ), py::arg( "icmp" ), py::arg( "ipt" ), py::arg( "ispt" ) = 0 )
        .def( "hasValue",
              py::overload_cast< const ASTERINTEGER &, const ASTERINTEGER &, const ASTERINTEGER &,
                                 const ASTERINTEGER & >( &SimpleFieldOnCellsReal::hasValue,
                                                         py::const_ ),
              R"(
Returns True  if the value of the `icmp` component of the field on the `ima` cell,
at the `ipt` point, at the `ispt` sub-point is affected.

Args:
    ima  (int): Index of cells.
    icmp (int): Index of component.
    ipt  (int): Index of point.
    ispt (int): Index of sub-point (default = 0).

Returns:
    bool: True  if the value is affected
        )",
              py::arg( "ima" ), py::arg( "icmp" ), py::arg( "ipt" ), py::arg( "ispt" ) = 0 )
        .def( "hasValue",
              py::overload_cast< const ASTERINTEGER &, const std::string &, const ASTERINTEGER &,
                                 const ASTERINTEGER & >( &SimpleFieldOnCellsReal::hasValue,
                                                         py::const_ ),
              R"(
Returns True  if the value of the `icmp` component of the field on the `ima` cell,
at the `ipt` point, at the `ispt` sub-point is affected.

Args:
    ima  (int): Index of cells.
    cmp (str): name of component.
    ipt  (int): Index of point.
    ispt (int): Index of sub-point (default = 0).

Returns:
    bool: True  if the value is affected
        )",
              py::arg( "ima" ), py::arg( "cmp" ), py::arg( "ipt" ), py::arg( "ispt" ) = 0 )
        .def( "setValue",
              py::overload_cast< const ASTERINTEGER &, const ASTERINTEGER &, const ASTERINTEGER &,
                                 const ASTERINTEGER &, const ASTERDOUBLE & >(
                  &SimpleFieldOnCellsReal::setValue ),
              R"(
Set the value of the `icmp` component of the field on the `ima` cell,
at the `ipt` point, at the `ispt` sub-point.

Args:
    ima  (int): Index of cells.
    icmp (int): Index of component.
    ipt  (int): Index of point.
    ispt (int): Index of sub-point.
    val (float) : value to set
        )",
              py::arg( "ima" ), py::arg( "icmp" ), py::arg( "ipt" ), py::arg( "ispt" ),
              py::arg( "val" ) )

        .def( "setValue",
              py::overload_cast< const ASTERINTEGER &, const ASTERINTEGER &, const ASTERINTEGER &,
                                 const ASTERDOUBLE & >( &SimpleFieldOnCellsReal::setValue ),
              R"(
Set the value of the `icmp` component of the field on the `ima` cell,
at the `ipt` point, at the `ispt=0` sub-point.

Args:
    ima  (int): Index of cells.
    icmp (int): Index of component.
    ipt  (int): Index of point.
    val (float) : value to set
        )",
              py::arg( "ima" ), py::arg( "icmp" ), py::arg( "ipt" ), py::arg( "val" ) )

        .def( "toNumpy", &SimpleFieldOnCellsReal::toNumpy,
              R"(
Returns two numpy arrays with shape ( number_of_cells_with_components, number_of_components )
The first array contains the field values while the second one is a mask
which is `True` if the corresponding value exists, `False` otherwise.

Where the mask is `False` the corresponding value is set to zero.

Returns:
    ndarray (float): Field values.
    ndarray (bool): Mask for the field values.
        )" )
        .def( "getValuesWithDescription",
              py::overload_cast< const VectorString &, const VectorLong & >(
                  &SimpleFieldOnCellsReal::getValuesWithDescription, py::const_ ),
              R"(
Returns values and description corresponding to given cmp and given cells

Arguments:
    cmps (list[str]): components to extract.
    cells (list[int]): list of cells.

Returns:
    values[list[double],
    tuple[cells[list[int]], cmps[list[int]]  points[list[int]], subpoints[list[int]]]
        )",
              py::arg( "cmps" ), py::arg( "cells" ) )
        .def( "getValuesWithDescription",
              py::overload_cast< const VectorString &, const VectorString & >(
                  &SimpleFieldOnCellsReal::getValuesWithDescription, py::const_ ),
              R"(
Returns values and description corresponding to given cmp and given cells

Arguments:
    cmps (list[str]): components to extract.
    groupsOfCells (list[str]): list of groups of cells to use.

Returns:
    values[list[double],
    tuple[cells[list[int]], cmps[list[int]]  points[list[int]], subpoints[list[int]]]
        )",
              py::arg( "cmps" ) = VectorString(), py::arg( "cells" ) = VectorString() )
        .def( "getCellsWithValues", &SimpleFieldOnCellsReal::getCellsWithValues, R"(
Returns the list of cells where the field is defined.

Returns:
    tuple (int): Indexes of cells where the field is defined.
        )" )
        .def( "allocate", &SimpleFieldOnCellsReal::allocate,
              R"(
            Allocate the field.

            Arguments:
                loc [str]: localization like 'ELEM'
                quantity [str]: physical quantity like 'DEPL_R'
                cmps [list[str]]: list of components.
                nbPG [int]: number of Gauss Point by cell
                nbSP [int]: number of sub-point by point.
            )",
              py::arg( "loc" ), py::arg( "quantity" ), py::arg( "cmps" ), py::arg( "nbPG" ),
              py::arg( "nbSP" ) = 1, py::arg( "zero" ) = false )
        .def( "getNumberOfComponents", &SimpleFieldOnCellsReal::getNumberOfComponents )
        .def( "getComponent", &SimpleFieldOnCellsReal::getComponent )
        .def( "getComponents", &SimpleFieldOnCellsReal::getComponents )
        .def( "getNumberOfCells", &SimpleFieldOnCellsReal::getNumberOfCells )
        .def( "getMaxNumberOfPoints", &SimpleFieldOnCellsReal::getMaxNumberOfPoints )
        .def( "getMaxNumberOfSubPoints", &SimpleFieldOnCellsReal::getMaxNumberOfSubPoints )
        .def( "getMesh", &SimpleFieldOnCellsReal::getMesh, R"(Returns base mesh)" )
        .def( "getNumberOfPointsOfCell", &SimpleFieldOnCellsReal::getNumberOfPointsOfCell )
        .def( "getNumberOfSubPointsOfCell", &SimpleFieldOnCellsReal::getNumberOfSubPointsOfCell )
        .def( "getNumberOfComponentsForSubpointsOfCell",
              &SimpleFieldOnCellsReal::getNumberOfComponentsForSubpointsOfCell )
        .def( "setValues",
              py::overload_cast< const VectorLong &, const VectorString &, const VectorLong &,
                                 const VectorLong &, const VectorReal & >(
                  &SimpleFieldOnCellsReal::setValues ),
              R"(
            Set values for a given list of tuple (cell, cmp, ipg, isp, value).
            Each value of the tuple is given as a separated list.

            Arguments:
                cells (list[int]): list of cells.
                cmps (list[str)]: list of components
                npg (list[int]): list of point
                spt (list[int]): list of sub-point
                values (list[float]): list of values to set.
            )",
              py::arg( "cells" ), py::arg( "cmps" ), py::arg( "npg" ), py::arg( "spt" ),
              py::arg( "values" ) )
        .def( "setValues",
              py::overload_cast< const VectorReal & >( &SimpleFieldOnCellsReal::setValues ), R"(
             Set values for each cells and components as (cell_0_val_0, cell_0_val_1, ...)

            Arguments:
                values (list[float]): list of values to set.

            )",
              py::arg( "values" ) )
        .def( "getPhysicalQuantity", &SimpleFieldOnCellsReal::getPhysicalQuantity )
        .def( "getComponentsName2Index", &SimpleFieldOnCellsReal::getComponentsName2Index )
        .def( "getLocalization", &SimpleFieldOnCellsReal::getLocalization )
        .def( "updateValuePointers", &SimpleFieldOnCellsReal::updateValuePointers )
        .def(
            "toFieldOnCells",
            []( const SimpleFieldOnCellsRealPtr &f, const FiniteElementDescriptorPtr fed,
                const std::string option, const std::string nompar, bool extensionToZero ) {
                return toFieldOnCells( f, fed, option, nompar, extensionToZero );
            },
            R"(
            Converts to FieldOnCells

            Arguments:
                fed [FiniteElementDescriptor]: finite element descriptor
                option [str]: name of option like TOUT_INI_ELGA (default: " ")
                nompar [str]: name of parameter like DEPL_R (default: " ")
                zeroExtension [bool]: true if field could be extended to zero (default: True)

            Returns:
                FieldOnCellsReal: field converted.
            )",
            py::arg( "fed" ), py::arg( "option" ) = std::string(),
            py::arg( "nompar" ) = std::string(), py::arg( "zeroExtension" ) = true )
        .def(
            "toFieldOnNodes", []( const SimpleFieldOnCellsReal &f ) { return toFieldOnNodes( f ); },
            R"(
Convert to FieldOnNodes

Returns:
    FieldOnNodesReal: field converted
        )" )
        .def(
            "toSimpleFieldOnNodes",
            []( const SimpleFieldOnCellsReal &f ) { return toSimpleFieldOnNodes( f ); },
            R"(
Convert to SimpleFieldOnNodes

Returns:
    SimpleFieldOnNodesReal: field converted
        )" )
        .def( "restrict", &SimpleFieldOnCellsReal::restrict,
              R"(
            Return a new field restricted to the list of components and groups of cells given

            Arguments:
                cmps[list[str]]: filter on list of components
                If empty, all components are used
                groupsOfCells[list[str]]: filter on list of groups of cells (default=" ").
                If empty, the full mesh is used

            Returns:
                SimpleFieldOnCellsReal: field restricted.
            )",
              py::arg( "cmps" ) = VectorString(), py::arg( "groupsOfCells" ) = VectorString() )
        .def( "asPhysicalQuantity", &SimpleFieldOnCellsReal::asPhysicalQuantity,
              R"(
            Return a new field with a new physical quantity and renamed components.

            Arguments:
                physQuantity [str]: name of the new physical quantity
                map_cmps [dict[str, str]]: dict to rename components
                (only renamed component will be keeped)

            Returns:
                SimpleFieldOnCellsReal: field with name physical quantity.
            )",
              py::arg( "physQuantity" ), py::arg( "map_cmps" ) );
};
