/**
 * @file CouplingPairing.cxx
 * @brief Implementation de Coupling
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

#include "Coupling/CouplingPairing.h"

#include "aster_fort_ds.h"

#include "DataFields/FieldOnCellsBuilder.h"
#include "Messages/Messages.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

namespace {
// Description of a virtual coupling cell

const std::map< std::pair< std::string, std::string >, std::string > cplCellNits {
    { { "MEDPTR3", "MECA_2D_HHO1_F" }, "CN_T3S3_HHO1" },
    { { "MEDPTR6", "MECA_2D_HHO1_F" }, "CN_T6S3_HHO1" },
    { { "MEDPTR10", "MECA_2D_HHO1_F" }, "CN_T10S3_HHO1" },
    { { "MEDPQU4", "MECA_2D_HHO1_F" }, "CN_Q4S3_HHO1" },
    { { "MEDPQU8", "MECA_2D_HHO1_F" }, "CN_Q8S3_HHO1" },
    { { "MEDPQU9", "MECA_2D_HHO1_F" }, "CN_Q9S3_HHO1" },
    { { "MEDPQU12", "MECA_2D_HHO1_F" }, "CN_Q12S3_HHO1" },
    { { "MEDPTR3", "MECA_2D_HHO2_F" }, "CN_T3S3_HHO2" },
    { { "MEDPTR6", "MECA_2D_HHO2_F" }, "CN_T6S3_HHO2" },
    { { "MEDPTR10", "MECA_2D_HHO2_F" }, "CN_T10S3_HHO2" },
    { { "MEDPQU4", "MECA_2D_HHO2_F" }, "CN_Q4S3_HHO2" },
    { { "MEDPQU8", "MECA_2D_HHO2_F" }, "CN_Q8S3_HHO2" },
    { { "MEDPQU9", "MECA_2D_HHO2_F" }, "CN_Q9S3_HHO2" },
    { { "MEDPQU12", "MECA_2D_HHO2_F" }, "CN_Q12S3_HHO2" },
};

const std::map< std::pair< std::string, std::string >, std::string > cplCellPena {
    { { "MEPLSE2", "MECA_2D_HHO1_F" }, "CP_S2S3_HHO1" },
    { { "MEPLSE3", "MECA_2D_HHO1_F" }, "CP_S3S3_HHO1" },
    { { "MEPLSE4", "MECA_2D_HHO1_F" }, "CP_S4S3_HHO1" },
    { { "MEPLSE2", "MECA_2D_HHO2_F" }, "CP_S2S3_HHO2" },
    { { "MEPLSE3", "MECA_2D_HHO2_F" }, "CP_S3S3_HHO2" },
    { { "MEPLSE4", "MECA_2D_HHO2_F" }, "CP_S4S3_HHO2" },
    { { "MEPLSE2", "MECA_2D_HHO3_F" }, "CP_S2S3_HHO3" },
    { { "MEPLSE3", "MECA_2D_HHO3_F" }, "CP_S3S3_HHO3" },
    { { "MECA_FACE4", "MECA3DQU9_HHO1_F" }, "CP_Q4Q9_HHO1" },
    { { "MECA_FACE8", "MECA3DQU9_HHO1_F" }, "CP_Q8Q9_HHO1" },
    { { "MECA_FACE9", "MECA3DQU9_HHO1_F" }, "CP_Q9Q9_HHO1" },
    { { "MECA_FACE3", "MECA3DQU9_HHO1_F" }, "CP_T3Q9_HHO1" },
    { { "MECA_FACE6", "MECA3DQU9_HHO1_F" }, "CP_T6Q9_HHO1" },
    { { "MECA_FACE4", "MECA3DTR7_HHO1_F" }, "CP_Q4T7_HHO1" },
    { { "MECA_FACE8", "MECA3DTR7_HHO1_F" }, "CP_Q8T7_HHO1" },
    { { "MECA_FACE9", "MECA3DTR7_HHO1_F" }, "CP_Q9T7_HHO1" },
    { { "MECA_FACE3", "MECA3DTR7_HHO1_F" }, "CP_T3T7_HHO1" },
    { { "MECA_FACE6", "MECA3DTR7_HHO1_F" }, "CP_T6T7_HHO1" },
    { { "MECA_FACE4", "MECA3DQU9_HHO2_F" }, "CP_Q4Q9_HHO2" },
    { { "MECA_FACE8", "MECA3DQU9_HHO2_F" }, "CP_Q8Q9_HHO2" },
    { { "MECA_FACE9", "MECA3DQU9_HHO2_F" }, "CP_Q9Q9_HHO2" },
    { { "MECA_FACE3", "MECA3DQU9_HHO2_F" }, "CP_T3Q9_HHO2" },
    { { "MECA_FACE6", "MECA3DQU9_HHO2_F" }, "CP_T6Q9_HHO2" },
    { { "MECA_FACE4", "MECA3DTR7_HHO2_F" }, "CP_Q4T7_HHO2" },
    { { "MECA_FACE8", "MECA3DTR7_HHO2_F" }, "CP_Q8T7_HHO2" },
    { { "MECA_FACE9", "MECA3DTR7_HHO2_F" }, "CP_Q9T7_HHO2" },
    { { "MECA_FACE3", "MECA3DTR7_HHO2_F" }, "CP_T3T7_HHO2" },
    { { "MECA_FACE6", "MECA3DTR7_HHO2_F" }, "CP_T6T7_HHO2" },
    { { "MEPLSE2", "MEPLSE2" }, "CP_S2S2" },
    { { "MEPLSE2", "MEPLSE3" }, "CP_S2S3" },
    { { "MEPLSE3", "MEPLSE2" }, "CP_S3S2" },
    { { "MEPLSE3", "MEPLSE3" }, "CP_S3S3" },
    { { "MEPLSE4", "MEPLSE4" }, "CP_S4S4" },
    { { "MEPLSE4", "MEPLSE4" }, "CP_S4S3" },
    { { "MECA_FACE4", "MECA_FACE8" }, "CP_Q4Q8" },
    { { "MECA_FACE4", "MECA_FACE9" }, "CP_Q4Q9" },
    { { "MECA_FACE4", "MECA_FACE3" }, "CP_Q4T3" },
    { { "MECA_FACE4", "MECA_FACE6" }, "CP_Q4T6" },
    { { "MECA_FACE8", "MECA_FACE4" }, "CP_Q8Q4" },
    { { "MECA_FACE8", "MECA_FACE9" }, "CP_Q8Q9" },
    { { "MECA_FACE8", "MECA_FACE3" }, "CP_Q8T3" },
    { { "MECA_FACE8", "MECA_FACE6" }, "CP_Q8T6" },
    { { "MECA_FACE9", "MECA_FACE4" }, "CP_Q9Q4" },
    { { "MECA_FACE9", "MECA_FACE8" }, "CP_Q9Q8" },
    { { "MECA_FACE9", "MECA_FACE3" }, "CP_Q9T3" },
    { { "MECA_FACE9", "MECA_FACE6" }, "CP_Q9T6" },
    { { "MECA_FACE4", "MECA_FACE4" }, "CP_Q4Q4" },
    { { "MECA_FACE8", "MECA_FACE8" }, "CP_Q8Q8" },
    { { "MECA_FACE9", "MECA_FACE9" }, "CP_Q9Q9" },
    { { "MECA_FACE3", "MECA_FACE4" }, "CP_T3Q4" },
    { { "MECA_FACE3", "MECA_FACE8" }, "CP_T3Q8" },
    { { "MECA_FACE3", "MECA_FACE9" }, "CP_T3Q9" },
    { { "MECA_FACE3", "MECA_FACE6" }, "CP_T3T6" },
    { { "MECA_FACE6", "MECA_FACE4" }, "CP_T6Q4" },
    { { "MECA_FACE6", "MECA_FACE8" }, "CP_T6Q8" },
    { { "MECA_FACE6", "MECA_FACE9" }, "CP_T6Q9" },
    { { "MECA_FACE6", "MECA_FACE3" }, "CP_T6T3" },
    { { "MECA_FACE3", "MECA_FACE3" }, "CP_T3T3" },
    { { "MECA_FACE6", "MECA_FACE6" }, "CP_T6T6" },
};

const std::map< std::pair< std::string, std::string >, std::string > cplCellLagr {
    { { "MEPLSE2", "MEPLSE2" }, "CL_S2S2" },       { { "MEPLSE2", "MEPLSE3" }, "CL_S2S3" },
    { { "MEPLSE3", "MEPLSE2" }, "CL_S3S2" },       { { "MEPLSE3", "MEPLSE3" }, "CL_S3S3" },
    { { "MEPLSE4", "MEPLSE4" }, "CL_S4S4" },       { { "MEPLSE4", "MEPLSE3" }, "CL_S4S3" },
    { { "MECA_FACE4", "MECA_FACE8" }, "CL_Q4Q8" }, { { "MECA_FACE4", "MECA_FACE9" }, "CL_Q4Q9" },
    { { "MECA_FACE4", "MECA_FACE3" }, "CL_Q4T3" }, { { "MECA_FACE4", "MECA_FACE6" }, "CL_Q4T6" },
    { { "MECA_FACE8", "MECA_FACE4" }, "CL_Q8Q4" }, { { "MECA_FACE8", "MECA_FACE9" }, "CL_Q8Q9" },
    { { "MECA_FACE8", "MECA_FACE3" }, "CL_Q8T3" }, { { "MECA_FACE8", "MECA_FACE6" }, "CL_Q8T6" },
    { { "MECA_FACE9", "MECA_FACE4" }, "CL_Q9Q4" }, { { "MECA_FACE9", "MECA_FACE8" }, "CL_Q9Q8" },
    { { "MECA_FACE9", "MECA_FACE3" }, "CL_Q9T3" }, { { "MECA_FACE9", "MECA_FACE6" }, "CL_Q9T6" },
    { { "MECA_FACE4", "MECA_FACE4" }, "CL_Q4Q4" }, { { "MECA_FACE8", "MECA_FACE8" }, "CL_Q8Q8" },
    { { "MECA_FACE9", "MECA_FACE9" }, "CL_Q9Q9" }, { { "MECA_FACE3", "MECA_FACE4" }, "CL_T3Q4" },
    { { "MECA_FACE3", "MECA_FACE8" }, "CL_T3Q8" }, { { "MECA_FACE3", "MECA_FACE9" }, "CL_T3Q9" },
    { { "MECA_FACE3", "MECA_FACE6" }, "CL_T3T6" }, { { "MECA_FACE6", "MECA_FACE4" }, "CL_T6Q4" },
    { { "MECA_FACE6", "MECA_FACE8" }, "CL_T6Q8" }, { { "MECA_FACE6", "MECA_FACE9" }, "CL_T6Q9" },
    { { "MECA_FACE6", "MECA_FACE3" }, "CL_T6T3" }, { { "MECA_FACE3", "MECA_FACE3" }, "CL_T3T3" },
    { { "MECA_FACE6", "MECA_FACE6" }, "CL_T6T6" }, { { "CL_POI2D", "NO" }, "CL_POI2D" },
    { { "CL_POI3D", "NO" }, "CL_POI3D" },
};

// Update also coupling_type.h
enum Offset {
    COEF_PENA = 0,
    NB_PTS_INTER = 1,
    COOR_PTS_X = 2,
    COOR_PTS_Y = 10,
    NB_NODES_NITSCHE = 18,
    NODE_IDX_FACE = 19,
    SIZE = 27,
};
} // namespace

CouplingPairing::CouplingPairing( const std::string name, const ModelPtr model,
                                  const ASTERINTEGER verbosity )
    : DataStructure( name, 8, "PAIRING_SD" ), _model( model ), _verbosity( verbosity ) {
    if ( !_model )
        raiseAsterError( "Mesh is empty" );

    if ( !_model->isMechanical() )
        UTMESS( "F", "CONTACT1_2" );
};

void CouplingPairing::check() const {
    // Some checks
    for ( const auto &zone : _zones ) {
        zone->check( _model );
    }
}

ASTERBOOL CouplingPairing::compute() {

    if ( _zones.empty() ) {
        raiseAsterError( "CouplingZone vector is empty" );
    }
    this->check();

    // Pairing
    for ( auto &zone : _zones ) {
        AS_ASSERT( zone->compute() );
    }

    // Build FED
    this->buildFiniteElementDescriptor();

    return true;
}

VectorPairLong CouplingPairing::getListOfPairs() const {

    VectorPairLong returnValue;
    ASTERINTEGER nbPairs = this->getNumberOfPairs();

    if ( nbPairs == 0 ) {
        raiseAsterError( "No contact pairs: was the pairing performed correctly? " );
    }

    returnValue.reserve( nbPairs );

    for ( const auto &zone : _zones ) {
        VectorPairLong pairsOnZone = zone->getListOfPairs();
        auto nbPairsZone = pairsOnZone.size();

        for ( auto iPairZone = 0; iPairZone < nbPairsZone; iPairZone++ ) {
            returnValue.push_back( pairsOnZone[iPairZone] );
        };
    }

    return returnValue;
}

ASTERINTEGER CouplingPairing::getCplCellType( const CouplingMethod algo,
                                              const std::string &slavCellTypeName,
                                              const std::string &mastCellTypeName ) const {
    std::string cplTypeName;

    if ( slavCellTypeName.find( "HHO" ) != std::string::npos ) {
        UTMESS( "F", "COUPLAGE_2" );
    }

    const auto key = std::make_pair( slavCellTypeName, mastCellTypeName );

    if ( algo == CouplingMethod::Nitsche ) {
        if ( mastCellTypeName.find( "HHO" ) == std::string::npos ) {
            UTMESS( "F", "COUPLAGE_3" );
        }

        auto it = cplCellNits.find( key );
        if ( it != cplCellNits.end() ) {
            cplTypeName = it->second;
        }
    } else if ( algo == CouplingMethod::Penalization ) {
        auto it = cplCellPena.find( key );
        if ( it != cplCellPena.end() ) {
            cplTypeName = it->second;
        }
    } else if ( algo == CouplingMethod::Lagrangian ) {
        auto it = cplCellLagr.find( key );
        if ( it != cplCellLagr.end() ) {
            cplTypeName = it->second;
        }
    } else {
        AS_ABORT( "Not implemented" );
    };

    if ( cplTypeName.empty() ) {
        UTMESS( "F", "COUPLAGE_1" );
    }

    // Get index of type of coupling cell
    return _model->getFiniteElementDescriptor()->getElemTypeNume( cplTypeName );
}

void CouplingPairing::createVirtualElemForOrphelanNodesForCoupling(
    MapLong &cplElemType, const JeveuxContiguousCollectionLong meshConnectivity,
    std::vector< VectorLong > &listCplElem, std::vector< VectorPairLong > &listCplType,
    SetLong &slaveNodePaired, SetLong &slaveCellPaired ) {

    // cplElemType: the number of elements for a given type
    // listCplType: list of coupling cells attached to pair (cellType, iPair)
    // listCplElem: list of coupling cells

    // Get mesh
    auto mesh = getMesh();
    ASTERINTEGER modelDim = _model->getGeometricDimension();

    ASTERINTEGER iCplPair = this->getNumberOfPairs();

    const auto fed = _model->getFiniteElementDescriptor();

    const auto typeFE = fed->getFiniteElementType();
    typeFE->updateValuePointer();

    for ( ASTERINTEGER iZone = 0; iZone < _zones.size(); iZone++ ) {
        const auto zone = _zones[iZone];

        const auto algo = zone->getMethod();
        AS_ASSERT( algo != CouplingMethod::Undefined );

        if ( algo != CouplingMethod::Lagrangian ) {
            continue;
        }

        // Get pairing
        const auto listOfPairsZone = zone->getListOfPairs();
        const auto nbPairsZone = this->getNumberOfPairs();

        // Get slave cells on this zone
        auto slaveCells = zone->getSlaveCells();

        // Create vector of (virtual) coupling cells for this zone
        VectorPairLong listCplTypeZone;
        listCplTypeZone.reserve( slaveCells.size() );

        // Find slave cells that are not paired (because of LAGR_C)
        for ( auto &slavCellNume : slaveCells ) {
            if ( slaveCellPaired.count( slavCellNume ) == 0 ) {
                slaveCellPaired.insert( slavCellNume );

                // Get nodes of slave cell
                auto slav_cell_con = ( *meshConnectivity )[slavCellNume + 1]->toVector();

                // Loop on nodes on slave cell
                for ( auto &nodeNume : slav_cell_con ) {
                    // This node hasn't been paired
                    if ( slaveNodePaired.count( nodeNume ) == 0 ) {
                        slaveNodePaired.insert( nodeNume );

                        // Type of cell for slave: POI1
                        std::string slavCellTypeName;

                        if ( modelDim == 3 ) {
                            slavCellTypeName = "CL_POI3D";
                        } else {
                            slavCellTypeName = "CL_POI2D";
                        }

                        // Type of cell for master
                        std::string mastCellTypeName = "NO";

                        // Get index of type of contact cell
                        const ASTERINTEGER typeElemNume =
                            this->getCplCellType( algo, slavCellTypeName, mastCellTypeName );
                        AS_ASSERT( typeElemNume != -1 );

                        // Add type of contact element
                        if ( cplElemType.count( typeElemNume ) == 0 ) {
                            cplElemType[typeElemNume] = 0;
                        }
                        cplElemType[typeElemNume] += 1;

                        // New virtual element
                        iCplPair++;
                        listCplTypeZone.push_back( std::make_pair( typeElemNume, iCplPair ) );
                        _cell2Zone[iCplPair - 1] = iZone;

                        // Add the nodes of the new contact element
                        listCplElem.push_back( VectorLong( { nodeNume, typeElemNume } ) );
                    }
                }
            }
        }
        if ( !listCplTypeZone.empty() ) {
#ifdef ASTER_DEBUG_CXX
            std::cout << "Not paired nodes: " << listCplTypeZone.size() << std::endl;
#endif
            listCplType.push_back( listCplTypeZone );
        }
    }
}

void CouplingPairing::createVirtualElemForCoupling(
    MapLong &cplElemType, const JeveuxContiguousCollectionLong meshConnectivity,
    std::vector< VectorLong > &listCplElem, std::vector< VectorPairLong > &listCplType,
    SetLong &slaveNodePaired, SetLong &slaveCellPaired ) {

    // cplElemType: the number of elements for a given type
    // listCplType: list of coupling cells attached to pair (cellType, iPair)
    // listCplElem: list of coupling cells

    ASTERINTEGER iCplPair = 0;

    const auto fed = _model->getFiniteElementDescriptor();

    const auto typeFE = fed->getFiniteElementType();
    typeFE->updateValuePointer();

    for ( ASTERINTEGER iZone = 0; iZone < _zones.size(); iZone++ ) {
        const auto zone = _zones[iZone];

        const auto algo = zone->getMethod();
        AS_ASSERT( algo != CouplingMethod::Undefined );

        // Get pairing
        const auto listOfPairsZone = zone->getListOfPairs();
        const auto nbPairsZone = this->getNumberOfPairs();

        // Link between surface and volume
        MapLong surf2Volu;
        if ( algo == CouplingMethod::Nitsche ) {
            if ( getMesh()->isParallel() ) {
                UTMESS( "F", "CONTACT1_5" );
            }
            surf2Volu = zone->getSlaveCellsSurfToVolu();
            if ( surf2Volu.size() == 0 ) {
                UTMESS( "F", "CONTACT1_3" );
            }
        }

        // Create vector of (virtual) coupling cells
        VectorPairLong listCplTypeZone;
        listCplTypeZone.reserve( nbPairsZone );

        // Loop on pairs in zone
        for ( int iPair = 0; iPair < nbPairsZone; iPair++ ) {
            // Get pairing
            auto [slavCellNume, mastCellNume] = listOfPairsZone[iPair];

            // Get cell slave to construct (is volumic cell for Nitsche)
            auto slavCellUsedNume = slavCellNume;
            if ( algo == CouplingMethod::Nitsche ) {
                slavCellUsedNume = surf2Volu[slavCellNume];
            }

            // Get slave and master cell type

            if ( ( *typeFE )[slavCellUsedNume] <= 0 ) {
                UTMESS( "F", "COUPLAGE_4", {}, { slavCellUsedNume }, {} );
            }
            if ( ( *typeFE )[mastCellNume] <= 0 ) {
                UTMESS( "F", "COUPLAGE_5", {}, { mastCellNume }, {} );
            }
            const auto slavCellTypeName = fed->getElemTypeName( ( *typeFE )[slavCellUsedNume] );
            const auto mastCellTypeName = fed->getElemTypeName( ( *typeFE )[mastCellNume] );

            // Get index of type of coupling cell
            const ASTERINTEGER typeElemNume =
                this->getCplCellType( algo, slavCellTypeName, mastCellTypeName );
            AS_ASSERT( typeElemNume != -1 );

            // Add coupling element to list
            if ( cplElemType.count( typeElemNume ) == 0 ) {
                cplElemType[typeElemNume] = 0;
            }
            cplElemType[typeElemNume] += 1;

            // New virtual element
            iCplPair++;
            listCplTypeZone.push_back( std::make_pair( typeElemNume, iCplPair ) );
            _cell2Zone[iCplPair - 1] = iZone;

            // Get nodes
            auto slav_cell_con = ( *meshConnectivity )[slavCellUsedNume + 1];
            auto toAdd1 = slav_cell_con->toVector();
            slav_cell_con = JeveuxCollectionObject< ASTERINTEGER >();
            auto mast_cell_con = ( *meshConnectivity )[mastCellNume + 1];
            auto toAdd2 = mast_cell_con->toVector();

            // Coupling element on zone
            VectorLong cplElemZone;
            cplElemZone.reserve( toAdd1.size() + toAdd2.size() + 1 );

            // Copy slave nodes to coupling element
            cplElemZone.insert( cplElemZone.end(), toAdd1.begin(), toAdd1.end() );

            // Add slave nodes to list of paired nodes
            slaveNodePaired.insert( toAdd1.begin(), toAdd1.end() );

            // Add slave cell to list of paired cells
            slaveCellPaired.insert( slavCellNume );

            // Copy master nodes to coupling element
            cplElemZone.insert( cplElemZone.end(), toAdd2.begin(), toAdd2.end() );

            // Add type of coupling element
            cplElemZone.push_back( typeElemNume );

            // Add coupling element to all coupling elements
            listCplElem.push_back( cplElemZone );
        }
        // Add coupling elements of zone
        listCplType.push_back( listCplTypeZone );
    }

    AS_ASSERT( iCplPair == this->getNumberOfPairs() );
}

void CouplingPairing::buildFiniteElementDescriptor() {

    CALL_JEMARQ();

    // Mesh
    auto mesh = getMesh();
    const auto meshConnectivity = mesh->getConnectivity();

    // Get pairing parameters
    const ASTERINTEGER nbCplPairTot = this->getNumberOfPairs();

    // Create objets for nodes and cells
    VectorOfVectorsLong listCplElem;
    listCplElem.reserve( nbCplPairTot );
    std::vector< VectorPairLong > listCplType;
    listCplType.reserve( 2 );

    // Objects
    SetLong slaveNodePaired, slaveCellPaired;

    // Object for number of cells for each type of coupling cell
    MapLong cplElemType;

    // Clear map between zone and coupling elements
    _cell2Zone.clear();

    // Create virtual elements for coupling
    createVirtualElemForCoupling( cplElemType, meshConnectivity, listCplElem, listCplType,
                                  slaveNodePaired, slaveCellPaired );

    createVirtualElemForOrphelanNodesForCoupling( cplElemType, meshConnectivity, listCplElem,
                                                  listCplType, slaveNodePaired, slaveCellPaired );

    // Create finite element descriptor for virtual coupling elements
    _fed = std::make_shared< FiniteElementDescriptor >( mesh );

    // Create list of virtual nodes (none !)
    _fed->setNumberOfVirtualNodes( 0 );

    // Create list of virtual elements (NEMA object)
    auto cplResFEDNema = _fed->getVirtualCellsDescriptor();
    cplResFEDNema->allocate( listCplElem );

    // Number of groups of elements and length of FED for coupling element
    ASTERINTEGER nbGrel = 0, lielLont = 0;
    for ( auto &[type, size] : cplElemType ) {
        lielLont += size;
        nbGrel += 1;
    }
    lielLont += nbGrel;

    // Create list of elements (LIEL object)
    auto cplResFEDLiel = _fed->getListOfGroupsOfElements();
    cplResFEDLiel->allocate( nbGrel, lielLont, Variable );

    // Add virtual elements for each GREL
    for ( auto &[type, size] : cplElemType ) {
        VectorLong virtualCells;
        virtualCells.reserve( size );

        for ( auto &listCplTypeZone : listCplType ) {
            for ( auto &[typeElemNume, iCplPair] : listCplTypeZone ) {
                if ( typeElemNume == type ) {
                    // Virtual cells = index of cells is negative
                    virtualCells.push_back( -iCplPair );
                }
            }
        }
        virtualCells.push_back( type );
        cplResFEDLiel->push_back( virtualCells );
    }

    // Get parameters from standard model
    auto paramToCopy = _model->getFiniteElementDescriptor()->getParameters();
    paramToCopy->updateValuePointer();
    auto docu = paramToCopy->getInformationParameter();

    // Create LGRF object
    auto parameters = _fed->getParameters();
    parameters->allocate( 4 );
    ( *parameters )[0] = mesh->getName();
    ( *parameters )[1] = _model->getName();
    ( *parameters )[2] = ( *paramToCopy )[2];
    parameters->setInformationParameter( docu );

    // Adapt FED
    CALLO_ADALIG_WRAP( _fed->getName() );
    bool l_calc_rigi = false;
    CALLO_INITEL( _fed->getName(), (ASTERLOGICAL *)&l_calc_rigi );

    // Final building
    _fed->build();

    CALL_JEDEMA();
};

FieldOnCellsRealPtr CouplingPairing::getPairingField() const {
    CALL_JEMARQ();

    // Field for intersection points and other thing ...
    auto data = FieldOnCellsPtrBuilder< ASTERDOUBLE >( _fed, "RIGI_ELAS_CPL", "PPAIRR" );

    // Get pairing
    const ASTERINTEGER nbPairs = this->getNumberOfPairs();

    VectorPairLong listPairs = this->getListOfPairs();

    // Acces to list of cells
    const auto meshConnex = this->getMesh()->getConnectivity();
    auto grel = _fed->getListOfGroupsOfElements();
    auto nbGrel = data->getNumberOfGroupOfElements();

    // Create local mapping (for Nitsche)
    auto mapping = []( const VectorLong &surf_nodes, const VectorLong &volu_nodes ) {
        VectorLong mapping;
        mapping.reserve( surf_nodes.size() );

        for ( auto &nodeId : surf_nodes ) {
            auto i = 1;
            for ( auto &nId : volu_nodes ) {
                if ( nId == nodeId ) {
                    break;
                }
                i++;
            }
            mapping.push_back( i );
        }

        return mapping;
    };

    // Loop on groups of elements
    ASTERINTEGER nbPair = 0;
    for ( ASTERINTEGER iGrel = 0; iGrel < nbGrel; iGrel++ ) {
        auto nbElem = data->getNumberOfElements( iGrel );
        auto liel = ( *grel )[iGrel + 1];
        liel->updateValuePointer();
        const auto typeFE = ( *liel )[liel->size() - 1];
        const auto nameFE = _fed->getElemTypeName( typeFE );
        if ( nameFE == "CL_POI2D" || nameFE == "CL_POI3D" ) {
            AS_ASSERT( data->getSizeOfFieldOfElement( iGrel ) == 0 );
            continue;
        }
        AS_ASSERT( data->getSizeOfFieldOfElement( iGrel ) == Offset::SIZE );

        // Loop on elements
        for ( ASTERINTEGER iElem = 0; iElem < nbElem; iElem++ ) {
            // Get mesh cell index
            auto iPair = -( *liel )[iElem];

            if ( iPair <= nbPairs ) {
                const auto iLocaPair = iPair - 1;

                // Current coupling zone
                auto iZone = _cell2Zone.at( iLocaPair );
                const auto zone = _zones[iZone];

                // Adress in field
                auto shift = data->getShifting( iGrel, iElem );

                // Set coordinates of slave intersection points
                const auto inter =
                    zone->getIntersectionPoints( iLocaPair, CoordinatesSpace::Slave );

                // Set number of intersection points
                ( *data )[shift + Offset::NB_PTS_INTER] = inter.size();
                AS_ASSERT( inter.size() <= 8 );

                for ( ASTERINTEGER iInter = 0; iInter < inter.size(); iInter++ ) {
                    ( *data )[shift + Offset::COOR_PTS_X + iInter] = inter[iInter][0];
                    ( *data )[shift + Offset::COOR_PTS_Y + iInter] = inter[iInter][1];
                }

                // Coefficient of penalization
                ( *data )[shift + Offset::COEF_PENA] = zone->getCoefficient();

                // For Nitsche
                if ( zone->getMethod() == CouplingMethod::Nitsche ) {
                    auto [slavCellNume, mastCellNume] = listPairs[iLocaPair];
                    auto slav_surf_con = ( *meshConnex )[slavCellNume + 1]->toVector();
                    auto slavVoluNume = zone->getSlaveCellSurfToVolu( slavCellNume );
                    auto slav_volu_con = ( *meshConnex )[slavVoluNume + 1]->toVector();

                    auto mapLoc = mapping( slav_surf_con, slav_volu_con );
                    // Number of nodes
                    ( *data )[shift + Offset::NB_NODES_NITSCHE] = ASTERDOUBLE( mapLoc.size() );

                    // Mapping
                    auto i = 0;
                    for ( auto &nodeId : mapLoc ) {
                        ( *data )[shift + Offset::NODE_IDX_FACE + i++] = ASTERDOUBLE( nodeId );
                    }
                }
            }

            nbPair++;
        }
    }

    CALL_JEDEMA();

    return data;
};

ASTERINTEGER CouplingPairing::getNumberOfPairs() const {
    ASTERINTEGER nbPair = 0;

    for ( const auto &zone : _zones ) {
        nbPair += zone->getNumberOfPairs();
    }

    return nbPair;
};

void CouplingPairing::setVerbosity( const ASTERINTEGER verbosity ) {
    _verbosity = verbosity;

    for ( auto &zone : _zones ) {
        zone->setVerbosity( _verbosity );
    }
};
