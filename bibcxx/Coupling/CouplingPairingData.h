/**
 * @file CouplingPairingData.h
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

#include "DataFields/FieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Loads/PhysicalQuantity.h"

class CouplingPairingData : public DataStructure {
    /** Datastructure for pairing */
  protected:
    FieldOnCellsRealPtr _pairing;

    FieldOnCellsRealPtr _hhoProj;

  public:
    /** @brief No default constructor */
    CouplingPairingData() = delete;

    /** @brief Constructor with given name */
    CouplingPairingData( const std::string name, const FieldOnCellsRealPtr pairing,
                         const FieldOnCellsRealPtr hho_proj );

    /** @brief Constructor with automatic name */
    CouplingPairingData( const FieldOnCellsRealPtr pairing, const FieldOnCellsRealPtr hho_proj )
        : CouplingPairingData( ResultNaming::getNewResultName(), pairing, hho_proj ) {};

    /** @brief Get Finite Element Descriptor from pairing */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const {
        return _pairing->getDescription();
    };

    bool exists() const { return _pairing->exists(); };

    FieldOnCellsRealPtr getPairingField() const { return _pairing; };

    FieldOnCellsRealPtr getHHOProjector() const { return _hhoProj; };
};

using CouplingPairingDataPtr = std::shared_ptr< CouplingPairingData >;
