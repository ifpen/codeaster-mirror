/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2026 - EDF - www.code-aster.org                 */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#include "PythonBindings/NumpyInterface.h"

#include "DataFields/SimpleFieldOnCells.h"
#include "DataFields/SimpleFieldOnNodes.h"

/**
 * @brief Get _valuesList as NumPy array
 */
py::array_t< ASTERDOUBLE > coordinatesToNumpy( MeshCoordinatesFieldPtr coords ) {
    py::object parent = py::cast( coords );

    coords->updateValuePointers();
    auto values = coords->getValues();
    int ncol = 3;
    auto nrow = values->size() / ncol;
    std::vector< ssize_t > shape = { nrow, ncol };
    std::vector< ssize_t > strides = { 3 * sizeof( ASTERDOUBLE ), sizeof( ASTERDOUBLE ) };
    py::array_t< ASTERDOUBLE > array( shape, strides, values->getDataPtr(), parent );
    return array;
}

template < typename ValueType >
py::tuple simpleFieldOnNodesToNumpy( std::shared_ptr< SimpleFieldOnNodes< ValueType > > field ) {
    py::object parent = py::cast( field );

    field->updateValuePointers();

    auto values = field->getValues();
    auto bools = field->getLogicalValues();
    auto size = values->size();
    auto ncol = (ssize_t)field->getNumberOfComponents();
    auto nrow = size / ncol;
    std::vector< ssize_t > shape = { nrow, ncol };
    std::vector< ssize_t > strides = { ncol * (ssize_t)sizeof( ValueType ), sizeof( ValueType ) };
    std::vector< ssize_t > strides_mask = { ncol * (ssize_t)sizeof( bool ), sizeof( bool ) };

    py::array_t< ValueType > array_values( shape, strides, values->getDataPtr(), parent );
    py::array_t< bool > array_mask( shape, strides_mask, bools->getDataPtr(), parent );
    py::tuple ret = py::make_tuple( array_values, array_mask );
    return ret;
};

template < typename ValueType >
py::tuple simpleFieldOnCellsToNumpy( std::shared_ptr< SimpleFieldOnCells< ValueType > > field ) {
    py::object parent = py::cast( field );

    field->updateValuePointers();

    auto values = field->getValues();
    auto sizes = field->getSizes();
    auto bools = field->getLogicalValues();
    auto size = values->size();
    auto ncol = (ssize_t)field->getNumberOfComponents();
    auto nrow = size / ncol;
    std::vector< ssize_t > shape = { nrow, ncol };
    std::vector< ssize_t > strides = { ncol * (ssize_t)sizeof( ValueType ), sizeof( ValueType ) };
    std::vector< ssize_t > strides_mask = { ncol * (ssize_t)sizeof( bool ), sizeof( bool ) };

    std::vector< ssize_t > shape_sizes = { sizes->size() };
    std::vector< ssize_t > strides_sizes = { sizeof( ASTERINTEGER ) };

    py::array_t< ValueType > array_values( shape, strides, values->getDataPtr(), parent );
    py::array_t< bool > array_mask( shape, strides_mask, bools->getDataPtr(), parent );
    py::array_t< ASTERINTEGER > array_sizes( shape_sizes, strides_sizes, sizes->getDataPtr(),
                                             parent );

    py::tuple ret = py::make_tuple( array_values, array_mask, array_sizes );
    return ret;
}

template py::tuple simpleFieldOnNodesToNumpy< ASTERDOUBLE >( SimpleFieldOnNodesRealPtr );
template py::tuple simpleFieldOnNodesToNumpy< ASTERCOMPLEX >( SimpleFieldOnNodesComplexPtr );

template py::tuple simpleFieldOnCellsToNumpy< ASTERDOUBLE >( SimpleFieldOnCellsRealPtr );

void exportNumpyToPython( py::module_ &mod ) {
    mod.def( "_coordinatesToNumpy", &coordinatesToNumpy );
    mod.def( "_simpleFieldOnNodesToNumpyReal", &simpleFieldOnNodesToNumpy< ASTERDOUBLE > );
    mod.def( "_simpleFieldOnNodesToNumpyComplex", &simpleFieldOnNodesToNumpy< ASTERCOMPLEX > );
    mod.def( "_simpleFieldOnCellsToNumpyReal", &simpleFieldOnCellsToNumpy< ASTERDOUBLE > );
}
