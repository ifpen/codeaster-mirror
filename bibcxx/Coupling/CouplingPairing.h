/**
 * @file CouplingPairing.h
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

#include "Coupling/CouplingZonePairing.h"
#include "DataFields/FieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Meshes/MeshEnum.h"
#include "Meshes/MeshPairing.h"
#include "Modeling/Model.h"

class CouplingPairing : public DataStructure {
    /** Datastructure for pairing */
  protected:
    /* verboisity */
    ASTERINTEGER _verbosity;

    /** @brief Model */
    ModelPtr _model;

    /** @brief Finite element descriptor for virtual elements of coupling */
    FiniteElementDescriptorPtr _fed;

    std::vector< CouplingZonePairingPtr > _zones;

    MapLong _cell2Zone;

  protected:
    /** @brief Create virtual elements for coupling */
    void createVirtualElemForCoupling( MapLong &cplElemType,
                                       const JeveuxContiguousCollectionLong meshConnectivity,
                                       std::vector< VectorLong > &listCplElem,
                                       std::vector< VectorPairLong > &listCplType,
                                       SetLong &slaveNodePaired, SetLong &slaveCellPaired );

    void createVirtualElemForOrphelanNodesForCoupling(
        MapLong &cplElemType, const JeveuxContiguousCollectionLong meshConnectivity,
        std::vector< VectorLong > &listCplElem, std::vector< VectorPairLong > &listCplType,
        SetLong &slaveNodePaired, SetLong &slaveCellPaired );

    /** @brief Get index of coupling cell */
    ASTERINTEGER getCplCellType( const CouplingMethod algo, const std::string &slavCellTypeName,
                                 const std::string &mastCellTypeName ) const;

    /* Check mesh */
    void check() const;

  public:
    /** @brief No default constructor */
    CouplingPairing() = delete;

    /** @brief Constructor with given name */
    CouplingPairing( const std::string name, const ModelPtr model,
                     const ASTERINTEGER verbosity = 1 );

    /** @brief Constructor with automatic name */
    CouplingPairing( const ModelPtr model, const ASTERINTEGER verbosity = 1 )
        : CouplingPairing( ResultNaming::getNewResultName(), model, verbosity ) {};

    /** @brief Get model */
    ModelPtr getModel() const { return _model; };

    /** @brief Get mesh */
    BaseMeshPtr getMesh() const { return _model->getMesh(); };

    /** @brief Compute pairing quantities of all zones */
    ASTERBOOL compute();

    /** @brief Set group of slave cells */
    void addZone( const CouplingZonePairingPtr zone ) { _zones.push_back( zone ); };

    /** @brief Build Finite Element Descriptor from pairing */
    virtual void buildFiniteElementDescriptor();

    /** @brief Get Finite Element Descriptor from pairing */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _fed; };

    ASTERINTEGER getNumberOfPairs() const;

    VectorPairLong getListOfPairs() const;

    FieldOnCellsRealPtr getPairingField() const;

    FieldOnCellsRealPtr getHHOProjector() const;

    void setVerbosity( const ASTERINTEGER verbosity );

    ASTERINTEGER getVerbosity() const { return _verbosity; };
};

using CouplingPairingPtr = std::shared_ptr< CouplingPairing >;
