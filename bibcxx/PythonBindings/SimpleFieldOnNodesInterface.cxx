/**
 * @file SimpleFieldOnNodesInterface.cxx
 * @brief Interface python de SimpleFieldOnNodes
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

#include "PythonBindings/SimpleFieldOnNodesInterface.h"

#include "aster_pybind.h"

#include "DataFields/FieldConverter.h"
#include "PythonBindings/DataStructureInterface.h"

void exportSimpleFieldOnNodesToPython( py::module_ &mod ) {
    py::class_< SimpleFieldOnNodesReal, SimpleFieldOnNodesRealPtr, DataField >(
        mod, "SimpleFieldOnNodesReal" )
        .def( py::init( &initFactoryPtr< SimpleFieldOnNodesReal, std::string > ) )
        .def( py::init( &initFactoryPtr< SimpleFieldOnNodesReal, BaseMeshPtr > ),
              py::arg( "mesh" ) )
        .def(
            py::init(
                &initFactoryPtr< SimpleFieldOnNodesReal, BaseMeshPtr, std::string, VectorString > ),
            py::arg( "mesh" ), py::arg( "quantity" ), py::arg( "cmps" ) )
        .def( py::init( &initFactoryPtr< SimpleFieldOnNodesReal, BaseMeshPtr, std::string,
                                         VectorString, bool > ),
              py::arg( "mesh" ), py::arg( "quantity" ), py::arg( "cmps" ), py::arg( "prol_zero" ) )
        .def(
            "__getitem__", +[]( const SimpleFieldOnNodesReal &v,
                                const PairLong &i ) { return v.operator()( i.first, i.second ); } )
        .def(
            "__getitem__",
            +[]( const SimpleFieldOnNodesReal &v,
                 const std::pair< ASTERINTEGER, std::string > &i ) {
                return v.operator()( i.first, i.second );
            } )
        .def(
            "__setitem__", +[]( SimpleFieldOnNodesReal &v, const PairLong &i,
                                ASTERDOUBLE f ) { return v.operator()( i.first, i.second ) = f; } )
        .def(
            "__setitem__",
            +[]( SimpleFieldOnNodesReal &v, const std::pair< ASTERINTEGER, std::string > &i ) {
                return v.operator()( i.first, i.second );
            } )
        .def( "allocate", &SimpleFieldOnNodesReal::allocate,
              R"(
            Allocate the field.

            Arguments:
                quantity [str]: physical quantity like 'DEPL_R'
                cmps [list[str]]: list of components.
            )",
              py::arg( "quantity" ), py::arg( "cmps" ), py::arg( "zero" ) = false )
        .def( "setMesh", &SimpleFieldOnNodesReal::setMesh,
              R"(
            Set mesh.

            Arguments:
                mesh [BaseMesh]: mesh to set.
            )",
              py::arg( "mesh" ) )
        .def(
            "toFieldOnNodes",
            []( const SimpleFieldOnNodesRealPtr &f ) { return toFieldOnNodes( f ); },
            R"(
Convert to FieldOnNodes

Returns:
    FieldOnNodesReal: field converted
        )" )
        .def(
            "toFieldOnNodes",
            []( const SimpleFieldOnNodesRealPtr &f, const BaseDOFNumberingPtr &dofNum,
                bool zeroExtension ) { return toFieldOnNodes( f, dofNum ); },
            R"(
Convert to FieldOnNodes

Arguments:
    dofNum (BaseDOFNumbering): DOF numbering used to build FieldOnNodes
    zeroExtension (bool): true if field can be extended to zero (when missing values)

Returns:
    FieldOnNodesReal: field converted
        )",
            py::arg( "dofNum" ), py::arg( "zeroExtension" ) = PythonBool::False )
        .def( "_restrict", &SimpleFieldOnNodesReal::restrict,
              R"(
            Return a new field restricted to the list of components and groups of nodes given

            Arguments:
                cmps[list[str]]: filter on list of components
                If empty, all components are used
                groupsOfNodes[list[str]]: filter on list of groups of nodes (default=" ").
                If empty, the full mesh is used
                same_rank : - None: keep all nodes (default: None)
                            - True: keep the nodes which are owned by the current MPI-rank
                            - False: keep the nodes which are not owned by the current MPI-rank

            Returns:
                SimpleFieldOnNodesReal: field restricted.
            )",
              py::arg( "cmps" ) = VectorString(), py::arg( "groupsOfNodes" ) = VectorString(),
              py::arg( "same_rank" ) = PythonBool::None )
        .def( "asPhysicalQuantity", &SimpleFieldOnNodesReal::asPhysicalQuantity,
              R"(
            Return a new field with a new physical quantity and renamed components.

            Arguments:
                physQuantity [str]: name of the new physical quantity
                map_cmps [dict[str, str]]: dict to rename components
                (only renamed component will be keeped)

            Returns:
                SimpleFieldOnNodesReal: field with name physical quantity.
            )",
              py::arg( "physQuantity" ), py::arg( "map_cmps" ) )
        .def( "toNumpy", &SimpleFieldOnNodesReal::toNumpy, R"(
Returns two numpy arrays with shape ( number_of_components, space_dimension )
The first array contains the field values while the second one is a mask
which is `True` if the corresponding value exists, `False` otherwise.

Where the mask is `False` the corresponding value is set to zero.

Returns:
    ndarray (float): Field values.
    ndarray (bool): Mask for the field values.
        )" )
        .def( "getNumberOfComponents", &SimpleFieldOnNodesReal::getNumberOfComponents )
        .def( "getNumberOfNodes", &SimpleFieldOnNodesReal::getNumberOfNodes )
        .def( "getComponents", &SimpleFieldOnNodesReal::getComponents )
        .def( "getComponent", &SimpleFieldOnNodesReal::getComponent )
        .def( "hasComponent", &SimpleFieldOnNodesReal::hasComponent )
        .def( "getMesh", &SimpleFieldOnNodesReal::getMesh, R"(Returns base mesh)" )
        .def( "getPhysicalQuantity", &SimpleFieldOnNodesReal::getPhysicalQuantity )
        .def( "getLocalization", &SimpleFieldOnNodesReal::getLocalization )
        .def( "setValues",
              py::overload_cast< const VectorLong &, const VectorString &, const VectorReal & >(
                  &SimpleFieldOnNodesReal::setValues ),
              R"(
            Set values for a given list of triplet (node, cmp, value).
            Each value of the triplet is given as a separated list.

            Arguments:
                nodes (list[int]): list of nodes.
                cmps (list[str]): list of comp components
                values (list[float]): list of values to set.
            )",
              py::arg( "nodes" ), py::arg( "cmps" ), py::arg( "values" ) )
        .def( "setValues",
              py::overload_cast< const VectorReal & >( &SimpleFieldOnNodesReal::setValues ),
              R"(
             Set values for each nodes and components as (node_0_val_0, node_0_val_1, ...)

            Arguments:
                values (list[float]): list of values to set.

            )",
              py::arg( "values" ) )
        .def( "setValues",
              py::overload_cast< const std::vector< VectorReal > & >(
                  &SimpleFieldOnNodesReal::setValues ),
              R"(
            Set values for each nodes and components.

            Arguments:
                values (list[list[float]]): list of values to set.
                For each node, give the values for all component is a list.
            )",
              py::arg( "values" ) )
        .def( "setValues",
              py::overload_cast< const std::map< std::string, ASTERDOUBLE > &, const VectorLong & >(
                  &SimpleFieldOnNodesReal::setValues ),
              R"(
            Set values of the field where components and values are given as a dict.
            If the component is not present in the field then it is discarded
            Example: { "X1" : 0.0, "X3" : 0.0 }

            Arguments:
                value (dict[str, float]): dict of values to set (key: str, value: float)
                nodes (list[int]): list of nodes.
            )",
              py::arg( "value" ), py::arg( "nodes" ) )
        .def(
            "setValues",
            py::overload_cast< const std::map< std::string, ASTERDOUBLE > &, const VectorString & >(
                &SimpleFieldOnNodesReal::setValues ),
            R"(
            Set values of the field where components and values are given as a dict.
            If the component is not present in the field then it is discarded
            Example: { "X1" : 0.0, "X3" : 0.0 }

            Arguments:
                value (dict[str, float]): dict of values to set (key: str, value: float)
                groupsOfNodes (list[str]): list of groups. If empty, the full mesh is considered
            )",
            py::arg( "value" ), py::arg( "groupsOfNodes" ) = VectorString() )
        .def( "setValues",
              py::overload_cast< const ASTERDOUBLE >( &SimpleFieldOnNodesReal::setValues ),
              R"(
            Set the value everywhere.

            Arguments:
                value [float]: value to set everywhere.
            )",
              py::arg( "value" ) )
        .def( "getValuesWithDescription",
              py::overload_cast< const VectorString &, const VectorString & >(
                  &SimpleFieldOnNodesReal::getValuesWithDescription, py::const_ ),
              R"(
            Return the values of components of the field.

            Arguments:
               cmps (list[str]) : Extracted components or all components if it is empty.
               groups (list[str]): The extraction is limited to the given groups of nodes.

            Returns:
               tuple( values, description ): List of values and description.
                The description provides a tuple with( nodes ids, components ).
                 )",
              py::arg( "cmps" ) = VectorString(), py::arg( "groupsOfNodes" ) = VectorString() )

        .def( "getValuesWithDescription",
              py::overload_cast< const VectorString &, const VectorLong & >(
                  &SimpleFieldOnNodesReal::getValuesWithDescription, py::const_ ),
              R"(
            Return the values of components of the field.

            Arguments:
               cmps (list[str]) : Extracted components or all components if it is empty.
               nodes (list[int]): The extraction is limited to the given nodes.

            Returns:
               tuple( values, description ): List of values and description.
                The description provides a tuple with( nodes ids, components ).
                 )",
              py::arg( "cmps" ), py::arg( "nodes" ) )

        .def( "updateValuePointers", &SimpleFieldOnNodesReal::updateValuePointers );

    py::class_< SimpleFieldOnNodesComplex, SimpleFieldOnNodesComplexPtr, DataField >(
        mod, "SimpleFieldOnNodesComplex" )
        .def( py::init( &initFactoryPtr< SimpleFieldOnNodesComplex, std::string > ) )
        .def( py::init( &initFactoryPtr< SimpleFieldOnNodesComplex, BaseMeshPtr, std::string,
                                         VectorString, bool > ) )

        .def(
            "__getitem__", +[]( const SimpleFieldOnNodesComplex &v,
                                const PairLong &i ) { return v.operator()( i.first, i.second ); } )
        .def(
            "__setitem__", +[]( SimpleFieldOnNodesComplex &v, const PairLong &i,
                                ASTERCOMPLEX f ) { return v.operator()( i.first, i.second ) = f; } )
        .def( "setMesh", &SimpleFieldOnNodesComplex::setMesh,
              R"(
            Set mesh.

            Arguments:
                mesh [BaseMesh]: mesh to set.
            )",
              py::arg( "mesh" ) )
        .def( "toNumpy", &SimpleFieldOnNodesComplex::toNumpy,
              R"(
Returns two numpy arrays with shape ( number_of_components, space_dimension )
The first array contains the field values while the second one is a mask
which is `True` if the corresponding value exists, `False` otherwise.

Where the mask is `False` the corresponding value is set to zero.

Returns:
    ndarray (complex): Field values.
    ndarray (bool): Mask for the field values.
        )" )
        .def( "getNumberOfComponents", &SimpleFieldOnNodesComplex::getNumberOfComponents )
        .def( "getNumberOfNodes", &SimpleFieldOnNodesComplex::getNumberOfNodes )
        .def( "getComponents", &SimpleFieldOnNodesComplex::getComponents )
        .def( "getComponent", &SimpleFieldOnNodesComplex::getComponent )
        .def( "hasComponent", &SimpleFieldOnNodesComplex::hasComponent )
        .def( "getMesh", &SimpleFieldOnNodesComplex::getMesh, R"(Returns base mesh)" )
        .def( "getPhysicalQuantity", &SimpleFieldOnNodesComplex::getPhysicalQuantity )
        .def( "getLocalization", &SimpleFieldOnNodesComplex::getLocalization )
        .def( "updateValuePointers", &SimpleFieldOnNodesComplex::updateValuePointers );
};
