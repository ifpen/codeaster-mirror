
/**
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

#pragma once

#include "astercxx.h"

#include "Contact/ContactNew.h"
#include "Contact/ContactPairing.h"
#include "DataFields/FieldOnCells.h"
#include "LinearAlgebra/ElementaryMatrix.h"

class ContactComputation {
  private:
    // Contact Definition
    ContactNewPtr _contact;

    /** @brief Convert ELNO -> NOEU for virtual nodes */
    FieldOnNodesRealPtr convertVirtualField( const FieldOnCellsRealPtr field ) const;

    /** @brief Level of verbosity */
    ASTERINTEGER _verbosity;

  public:
    /** @brief Main constructor */
    ContactComputation( const ContactNewPtr contact ) : _contact( contact ) {};

    /** @brief Restricted constructor (Set) to support pickling */
    ContactComputation( const py::tuple &tup )
        : ContactComputation( tup[0].cast< ContactNewPtr >() ) {};

    /** @brief Method (Get) to support pickling */
    py::tuple _getState() const { return py::make_tuple( _contact ); };

    /** @brief Compute contact mortar matrix  */
    ElementaryMatrixDisplacementRealPtr contactMortarMatrix() const;

    /** @brief Compute field for input data in elementary computations of contact  */
    FieldOnCellsRealPtr contactData( const ContactPairingPtr contPairing,
                                     const MaterialFieldPtr mater,
                                     const bool &initial_contact ) const;

    /** @brief Compute contact coefficient field (COEF_CONT) */
    std::pair< FieldOnNodesRealPtr, FieldOnNodesRealPtr > contactCoefficient() const;

    /** @brief Set verbosity */
    void setVerbosity( const ASTERINTEGER &level ) { _verbosity = level; }

    /** @brief Get verbosity */
    ASTERINTEGER getVerbosity() const { return _verbosity; }

    /** @brief Get virtual cells */
    ASTERBOOL hasVirtualCells( const ContactPairingPtr contPairing ) const {
        return ( contPairing->hasVirtualCells() );
    }
};

using ContactComputationPtr = std::shared_ptr< ContactComputation >;
