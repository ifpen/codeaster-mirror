/**
 * @file ContactPairing.cxx
 * @brief Implementation de Contact
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

#include "Contact/ContactPairing.h"

#include "Messages/Messages.h"

ContactPairing::ContactPairing( const std::string name, const ContactNewPtr contDefi )
    : DataStructure( name, 8, "PAIRING_SD" ),
      _contDefi( contDefi ),
      _mesh( contDefi->getMesh() ),
      _fed( nullptr ),
      _fedCreated( false ),
      _verbosity( 1 ) {
    if ( !_mesh )
        raiseAsterError( "Mesh is empty" );

    // Set verbosity
    setVerbosity( contDefi->getVerbosity() );

    // Create object for mesh coordinates
    _currentCoordinates = std::make_shared< MeshCoordinatesField >( *( _mesh->getCoordinates() ) );
    _contDefi->setCoordinates( _currentCoordinates );

    // Be sure that zones is not empty and get size of zones
    int nbZoneCont = _contDefi->getNumberOfContactZones();
    if ( nbZoneCont == 0 )
        raiseAsterError( "ContactZone vector is empty" );
};

void ContactPairing::setVerbosity( const ASTERINTEGER &level ) {
    _verbosity = level;
    _contDefi->setVerbosity( getVerbosity() );
}

ASTERBOOL ContactPairing::compute( ASTERINTEGER &indexZone ) {

    if ( indexZone < 0 || indexZone >= _contDefi->getNumberOfContactZones() ) {
        throw std::out_of_range( "The zone index should be between 0  and " +
                                 std::to_string( _contDefi->getNumberOfContactZones() - 1 ) );
    }

    ASTERBOOL returnValue;

    // Get current zone
    auto zone = _contDefi->getContactZone( indexZone );
    AS_ASSERT( !zone->hasSmoothing() );

    // Get pairing method
    auto variant = zone->getPairingParameter()->getAlgorithm();
    if ( variant != PairingAlgo::Mortar ) {
        AS_ABORT( "Not expected" );
    }

    // Get distance ratio
    auto dist_pairing = zone->getPairingParameter()->getDistanceRatio();

    // Tolerance for pairing (zero value for geometrical operations)
    // ASTERDOUBLE pair_tole = 1e-8;
    auto pair_tole = zone->getPairingParameter()->getPairingTolerance();

    // Tolerance for removing intersection cells
    // ASTERDOUBLE area_tole = 1e-2;
    auto area_tole = zone->getPairingParameter()->getAreaIntersectionTolerance();

    // Pairing
    zone->setVerbosity( getVerbosity() );
    returnValue = zone->pairing( dist_pairing, pair_tole, area_tole );

    return returnValue;
}

ASTERBOOL ContactPairing::compute() {

    // Pairing
    for ( ASTERINTEGER i = 0; i < _contDefi->getNumberOfContactZones(); i++ ) {
        AS_ASSERT( compute( i ) );
    }

    // Build FED
    this->buildFiniteElementDescriptor();

    return true;
}

void ContactPairing::clearPairing( const ASTERINTEGER &indexZone ) {
    if ( indexZone < 0 || indexZone >= _contDefi->getNumberOfContactZones() ) {
        throw std::out_of_range( "The zone index should be between 0  and " +
                                 std::to_string( _contDefi->getNumberOfContactZones() - 1 ) );
    }

    _contDefi->clearPairing( indexZone );
}

VectorPairLong ContactPairing::getListOfPairs( const ASTERINTEGER &indexZone ) const {

    if ( indexZone < 0 || indexZone >= _contDefi->getNumberOfContactZones() ) {
        throw std::out_of_range( "The zone index should be between 0  and " +
                                 std::to_string( _contDefi->getNumberOfContactZones() - 1 ) );
    }

    auto zone = _contDefi->getContactZone( indexZone );

    return ( zone->getMeshPairing()->getListOfPairs() );
}

ASTERINTEGER ContactPairing::getNumberOfPairs( const ASTERINTEGER &indexZone ) const {
    if ( indexZone < 0 || indexZone >= _contDefi->getNumberOfContactZones() ) {
        throw std::out_of_range( "The zone index should be between 0  and " +
                                 std::to_string( _contDefi->getNumberOfContactZones() - 1 ) );
    }
    return ( _contDefi->getContactZone( indexZone )->getMeshPairing()->getNumberOfPairs() );
}

ASTERINTEGER ContactPairing::getNumberOfPairs() const {
    ASTERINTEGER returnValue;
    returnValue = 0;
    for ( auto indexZone = 0; indexZone < _contDefi->getNumberOfContactZones(); indexZone++ ) {
        returnValue += _contDefi->getMeshPairing( indexZone )->getNumberOfPairs();
    }

    return returnValue;
};

ASTERINTEGER ContactPairing::getNumberOfZones() const {
    return _contDefi->getNumberOfContactZones();
}

VectorPairLong ContactPairing::getListOfPairs() const {

    VectorPairLong returnValue;
    ASTERINTEGER nbPairs = getNumberOfPairs();

    if ( nbPairs == 0 ) {
        raiseAsterError( "No contact pairs: was the pairing performed correctly? " );
    }

    returnValue.reserve( nbPairs );

    for ( ASTERINTEGER indexZone = 0; indexZone < _contDefi->getNumberOfContactZones();
          indexZone++ ) {
        VectorPairLong pairsOnZone = getListOfPairs( indexZone );
        auto nbPairsZone = pairsOnZone.size();

        for ( auto iPairZone = 0; iPairZone < nbPairsZone; iPairZone++ ) {
            returnValue.push_back( pairsOnZone[iPairZone] );
        };
    }

    return returnValue;
}

VectorOfVectorsReal
ContactPairing::getIntersectionPoints( ASTERINTEGER &indexZone,
                                       const CoordinatesSpace coorSpace ) const {

    VectorOfVectorsReal returnValue;
    ASTERINTEGER nbPairs = getNumberOfPairs( indexZone );
    if ( nbPairs == 0 ) {
        raiseAsterError( "No contact pairs: was the pairing performed correctly? " );
    }
    returnValue.reserve( nbPairs );
    const auto meshPair = _contDefi->getContactZone( indexZone )->getMeshPairing();
    for ( auto iPair = 0; iPair < nbPairs; iPair++ ) {
        VectorOfVectorsReal interOnZone = meshPair->getIntersectionPoints( iPair, coorSpace );
        ASTERINTEGER nbInter = meshPair->getNumberOfIntersectionPoints( iPair );
        if ( nbInter != 0 ) {
            VectorReal vectVale;
            for ( auto iInter = 0; iInter < nbInter; iInter++ ) {
                vectVale.insert( vectVale.end(), interOnZone[iInter].begin(),
                                 interOnZone[iInter].end() );
            };
            returnValue.push_back( vectVale );
        }
    }
    return returnValue;
}

std::vector< VectorOfVectorsReal >
ContactPairing::getIntersectionPoints( const CoordinatesSpace coorSpace ) const {

    std::vector< VectorOfVectorsReal > returnValue;
    const ASTERINTEGER nbZoneCont = _contDefi->getNumberOfContactZones();
    if ( nbZoneCont == 0 )
        raiseAsterError( "ContactZone vector is empty " );
    returnValue.reserve( nbZoneCont );

    for ( ASTERINTEGER iZone = 0; iZone < nbZoneCont; iZone++ ) {
        VectorOfVectorsReal interSlave;
        interSlave = ContactPairing::getIntersectionPoints( iZone, coorSpace );
        returnValue.push_back( interSlave );
    }
    return returnValue;
}

ASTERINTEGER ContactPairing::getContCellIndx( const ContactAlgo contAlgo,
                                              std::string slavCellTypeName,
                                              std::string mastCellTypeName ) {
    ASTERINTEGER cellIndx = -1;

    if ( contAlgo == ContactAlgo::Lagrangian ) {
        for ( int iContType = 0; iContType < contLagrType; iContType++ ) {
            if ( slavCellTypeName == contCellLagr[iContType].slavCellType ) {
                if ( mastCellTypeName == contCellLagr[iContType].mastCellType ) {
                    AS_ASSERT( cellIndx == -1 )
                    cellIndx = iContType;
                }
            }
        }
    } else if ( contAlgo == ContactAlgo::Nitsche ) {
        for ( int iContType = 0; iContType < contNitsType; iContType++ ) {
            if ( slavCellTypeName == contCellNits[iContType].slavCellType ) {
                if ( mastCellTypeName == contCellNits[iContType].mastCellType ) {
                    AS_ASSERT( cellIndx == -1 )
                    cellIndx = iContType;
                }
            }
        }
    } else if ( contAlgo == ContactAlgo::Penalization ) {
        for ( int iContType = 0; iContType < contPenaType; iContType++ ) {
            if ( slavCellTypeName == contCellPena[iContType].slavCellType ) {
                if ( mastCellTypeName == contCellPena[iContType].mastCellType ) {
                    AS_ASSERT( cellIndx == -1 )
                    cellIndx = iContType;
                }
            }
        }
    } else {
        AS_ABORT( "Not implemented" );
    };

    return cellIndx;
}

ASTERINTEGER ContactPairing::getContCellType( const ContactAlgo contAlgo,
                                              const ASTERINTEGER cellIndx, const bool lAxis,
                                              const bool lFric ) {

    AS_ASSERT( cellIndx != -1 );

    ASTERINTEGER contTypeNume = -1;

    // Get finite element descriptor of model
    auto modelFEDesc = _contDefi->getModel()->getFiniteElementDescriptor();

    // Get name of type of contact cell
    std::string contTypeName;
    if ( contAlgo == ContactAlgo::Lagrangian ) {
        if ( lFric ) {
            contTypeName = contCellLagr[cellIndx].fricElemType;
        } else {
            contTypeName = contCellLagr[cellIndx].contElemType;
        }
    } else if ( contAlgo == ContactAlgo::Nitsche ) {
        if ( lFric ) {
            contTypeName = contCellNits[cellIndx].fricElemType;
        } else {
            contTypeName = contCellNits[cellIndx].contElemType;
        }
    } else if ( contAlgo == ContactAlgo::Penalization ) {
        contTypeName = contCellPena[cellIndx].contElemType;

    } else {
        AS_ABORT( "Not implemented" );
    };

    if ( lAxis ) {
        contTypeName.append( "A" );
    }

    // Get index of type of contact cell
    contTypeNume = modelFEDesc->getElemTypeNume( contTypeName );
    return contTypeNume;
}

void ContactPairing::createVirtualElemForContact(
    const ASTERLOGICAL lAxis, const int nbZoneCont, MapLong &contactElemType,
    const JeveuxContiguousCollectionLong meshConnectivity, std::vector< VectorLong > &listContElem,
    std::vector< VectorPairLong > &listContType, SetLong &slaveNodePaired,
    SetLong &slaveCellPaired ) {

    // contactElemType: the number of elements for a given type
    // listContType: list of contact cells attached to pair (cellType, iPair)
    // listContElem: list of contact cells

    auto mesh = getMesh();
    ASTERINTEGER iContPair = 0;

    // Loop on contact zones
    for ( int iZone = 0; iZone < nbZoneCont; iZone++ ) {
        // Get current zone
        auto zone = _contDefi->getContactZone( iZone );

        // Get parameters for this zone
        auto contAlgo = zone->getContactParameter()->getAlgorithm();
        auto lFric = zone->getFrictionParameter()->hasFriction();

        // Get pairing for this zone
        auto listOfPairsZone = this->getListOfPairs( iZone );
        auto nbPairsZone = this->getNumberOfPairs( iZone );

        // Link between surface and volume
        MapLong surf2Volu;
        if ( contAlgo == ContactAlgo::Nitsche ) {
            if ( mesh->isParallel() ) {
                UTMESS( "F", "CONTACT1_5" );
            }
            surf2Volu = zone->getSlaveCellsSurfToVolu();
        }

        // Create vector of (virtual) contact cells for this zone
        VectorPairLong listContTypeZone;
        listContTypeZone.reserve( nbPairsZone );

        // Loop on pairs in zone
        for ( int iPair = 0; iPair < nbPairsZone; iPair++ ) {
            // Get pairing of current zone
            auto [slavCellNume, mastCellNume] = listOfPairsZone[iPair];

            // Get cell slave to construct (is volumic cell for Nitsche)
            auto slavCellUsedNume = slavCellNume;
            if ( contAlgo == ContactAlgo::Nitsche ) {
                if ( surf2Volu.size() == 0 ) {
                    UTMESS( "F", "CONTACT1_3" );
                }
                slavCellUsedNume = surf2Volu[slavCellNume];
            }

            // Get slave and master cell type
            auto slavCellTypeName = mesh->getCellTypeName( slavCellUsedNume );
            auto mastCellTypeName = mesh->getCellTypeName( mastCellNume );

            // Get index of contact cell (in fixed lists)
            ASTERINTEGER cellIndx = getContCellIndx( contAlgo, slavCellTypeName, mastCellTypeName );
            AS_ASSERT( cellIndx != -1 );

            // Get index of type of contact cell
            ASTERINTEGER typeElemNume = getContCellType( contAlgo, cellIndx, lAxis, lFric );
            AS_ASSERT( typeElemNume != -1 );

            // Add contact element to list
            if ( contactElemType.count( typeElemNume ) == 0 ) {
                contactElemType[typeElemNume] = 0;
            }
            contactElemType[typeElemNume] += 1;

            // New virtual element
            iContPair++;
            listContTypeZone.push_back( std::make_pair( typeElemNume, iContPair ) );
            _cell2Zone[iContPair - 1] = iZone;

            // Get nodes
            auto slav_cell_con = ( *meshConnectivity )[slavCellUsedNume + 1];
            auto toAdd1 = slav_cell_con->toVector();
            slav_cell_con = JeveuxCollectionObject< ASTERINTEGER >();
            auto mast_cell_con = ( *meshConnectivity )[mastCellNume + 1];
            auto toAdd2 = mast_cell_con->toVector();

            // Contact element on zone
            VectorLong contactElemZone;
            contactElemZone.reserve( toAdd1.size() + toAdd2.size() + 1 );

            // Copy slave nodes to contact element
            contactElemZone.insert( contactElemZone.end(), toAdd1.begin(), toAdd1.end() );

            // Add slave nodes to list of paired nodes
            slaveNodePaired.insert( toAdd1.begin(), toAdd1.end() );

            // Add slave cell to list of paired cells
            slaveCellPaired.insert( slavCellNume );

            // Copy master nodes to contact element
            contactElemZone.insert( contactElemZone.end(), toAdd2.begin(), toAdd2.end() );

            // Add type of contact element
            contactElemZone.push_back( typeElemNume );

            // Add contact element to all contact elements
            listContElem.push_back( contactElemZone );
        }
        // Add contact elements of zone
        listContType.push_back( listContTypeZone );
    }
    AS_ASSERT( iContPair == this->getNumberOfPairs() )
}

void ContactPairing::createVirtualElemForOrphelanNodes(
    const ASTERLOGICAL lAxis, const int nbZoneCont, MapLong &contactElemType,
    const JeveuxContiguousCollectionLong meshConnectivity, std::vector< VectorLong > &listContElem,
    std::vector< VectorPairLong > &listContType, SetLong &slaveNodePaired,
    SetLong &slaveCellPaired ) {

    // Get mesh
    auto mesh = getMesh();

    ASTERINTEGER iContPair = 0;
    iContPair = this->getNumberOfPairs();

    // Get model
    auto model = _contDefi->getModel();
    ASTERINTEGER modelDim = model->getGeometricDimension();

    // Loop on contact zones
    for ( int iZone = 0; iZone < nbZoneCont; iZone++ ) {
        // Get current zone
        auto zone = _contDefi->getContactZone( iZone );

        // Get contact parameters for this zone
        auto contAlgo = zone->getContactParameter()->getAlgorithm();
        auto lFric = zone->getFrictionParameter()->hasFriction();

        // Get slave cells on this zone
        auto slaveCells = zone->getSlaveCells();

        // Create vector of (virtual) contact cells for this zone
        VectorPairLong listContTypeZone;
        listContTypeZone.reserve( slaveCells.size() );

        // Find slave cells that are not paired (because of LAGR_C)
        for ( auto &slavCellNume : slaveCells ) {
            if ( slaveCellPaired.count( slavCellNume ) == 0 ) {
                slaveCellPaired.insert( slavCellNume );

                // Get nodes of slave cell
                auto slav_cell_con = ( *meshConnectivity )[slavCellNume + 1]->toVector();

                // Get slave cell type
                auto slavCellTypeName = mesh->getCellTypeName( slavCellNume );

                // Select number of nodes with LAGR_C: P2/P1 for 3D and P2/P2 for 2D
                ASTERINTEGER nno_lgar = 0;

                if ( slavCellTypeName == "SEG2" ) {
                    nno_lgar = 2;
                } else if ( slavCellTypeName == "SEG3" ) {
                    nno_lgar = 3;
                } else if ( slavCellTypeName == "TRIA3" || slavCellTypeName == "TRIA6" ||
                            slavCellTypeName == "TRIA7" ) {
                    nno_lgar = 3;
                } else if ( slavCellTypeName == "QUAD4" || slavCellTypeName == "QUAD8" ||
                            slavCellTypeName == "QUAD9" ) {
                    nno_lgar = 4;
                } else {
                    AS_ABORT( slavCellTypeName + " not supported" );
                }

                // Loop on nodes on slave cell
                ASTERINTEGER nno = 0;
                for ( auto &nodeNume : slav_cell_con ) {
                    nno++;
                    // This node hasn't been paired
                    if ( slaveNodePaired.count( nodeNume ) == 0 ) {
                        slaveNodePaired.insert( nodeNume );

                        // Type of cell for slave: POI1
                        auto slavCellTypeName = "POI1";

                        // Type of cell for master
                        std::string mastCellTypeName;
                        if ( nno <= nno_lgar ) {
                            mastCellTypeName = "LAG" + std::to_string( modelDim );
                        } else {
                            mastCellTypeName = "NOLAG" + std::to_string( modelDim );
                        }

                        if ( contAlgo == ContactAlgo::Lagrangian ) {

                        } else if ( contAlgo == ContactAlgo::Nitsche ) {
                            continue;
                        } else if ( contAlgo == ContactAlgo::Penalization ) {
                            continue;
                        } else {
                            AS_ABORT( "Not implemented" );
                        }

                        // Get index of contact cell (in fixed lists)
                        ASTERINTEGER cellIndx =
                            getContCellIndx( contAlgo, slavCellTypeName, mastCellTypeName );
                        AS_ASSERT( cellIndx != -1 );

                        // Number of nodes
                        ASTERINTEGER nbNodesCell = contCellLagr[cellIndx].nbNode;
                        AS_ASSERT( nbNodesCell == 1 );

                        // Get index of type of contact cell
                        ASTERINTEGER typeElemNume =
                            getContCellType( contAlgo, cellIndx, lAxis, lFric );
                        AS_ASSERT( typeElemNume != -1 )

                        // Add type of contact element
                        if ( contactElemType.count( typeElemNume ) == 0 ) {
                            contactElemType[typeElemNume] = 0;
                        }
                        contactElemType[typeElemNume] += 1;

                        // New virtual element
                        iContPair++;
                        listContTypeZone.push_back( std::make_pair( typeElemNume, iContPair ) );
                        _cell2Zone[iContPair - 1] = iZone;

                        // Add the nodes of the new contact element
                        listContElem.push_back( VectorLong( { nodeNume, typeElemNume } ) );
                    }
                }
            }
        }
        if ( !listContTypeZone.empty() ) {
#ifdef ASTER_DEBUG_CXX
            std::cout << "Not paired nodes: " << listContTypeZone.size() << std::endl;
#endif
            listContType.push_back( listContTypeZone );
        }
    }
};

void ContactPairing::buildFiniteElementDescriptor() {

    CALL_JEMARQ();

    // Model
    auto model = _contDefi->getModel();
    const ASTERLOGICAL lAxis = model->existsAxis();

    // Mesh
    auto mesh = getMesh();
    const auto meshConnectivity = mesh->getConnectivity();

    // Get pairing parameters
    const ASTERINTEGER nbZoneCont = _contDefi->getNumberOfContactZones();
    const ASTERINTEGER nbContPairTot = this->getNumberOfPairs();

    // Create objets for nodes and cells
    VectorOfVectorsLong listContElem;
    listContElem.reserve( nbContPairTot );
    std::vector< VectorPairLong > listContType;
    listContType.reserve( 2 * nbZoneCont );

    // Objects
    SetLong slaveNodePaired, slaveCellPaired;

    // Object for number of cells for each type of contact cell
    MapLong contactElemType;

    // Clear map between zone and contact elements
    _cell2Zone.clear();

    // Create virtual elements for contact
    createVirtualElemForContact( lAxis, nbZoneCont, contactElemType, meshConnectivity, listContElem,
                                 listContType, slaveNodePaired, slaveCellPaired );

    // Create virtual elements for orphelan nodes
    createVirtualElemForOrphelanNodes( lAxis, nbZoneCont, contactElemType, meshConnectivity,
                                       listContElem, listContType, slaveNodePaired,
                                       slaveCellPaired );

    if ( listContElem.size() != 0 ) {
        // Create finite element descriptor for virtual contact elements
        _fed = std::make_shared< FiniteElementDescriptor >( mesh );

        // Create list of virtual nodes (none !)
        _fed->setNumberOfVirtualNodes( 0 );

        // Create list of virtual elements (NEMA object)
        auto ContactResFEDNema = _fed->getVirtualCellsDescriptor();
        ContactResFEDNema->allocate( listContElem );

        // Number of groups of elements and length of FED for contact element
        ASTERINTEGER nbGrel = 0, lielLont = 0;
        for ( auto &[type, size] : contactElemType ) {
            lielLont += size;
            nbGrel += 1;
        }
        lielLont += nbGrel;

        // Create list of elements (LIEL object)
        auto ContactResFEDLiel = _fed->getListOfGroupsOfElements();
        ContactResFEDLiel->allocate( nbGrel, lielLont, Variable );

        // Add virtual elements for each GREL
        for ( auto &[type, size] : contactElemType ) {
            VectorLong virtualCells;
            virtualCells.reserve( size );

            // ASTERINTEGER iZone = 0;
            for ( auto &listContTypeZone : listContType ) {
                for ( auto &[typeElemNume, iContPair] : listContTypeZone ) {
                    if ( typeElemNume == type ) {
                        // Virtual cells = index of cells is negative
                        virtualCells.push_back( -iContPair );
                    }
                }
            }
            virtualCells.push_back( type );
            ContactResFEDLiel->push_back( virtualCells );
        }

        // Map between index of global pair and index of local pair in zone
        _globPairToLocaPair.clear();
        for ( ASTERINTEGER indexZone = 0; indexZone < nbZoneCont; indexZone++ ) {
            if ( indexZone == 0 ) {
                _globPairToLocaPair[indexZone] = 0;
            } else {
                auto zone = _contDefi->getContactZone( indexZone - 1 );
                ASTERINTEGER nbPairs = zone->getMeshPairing()->getNumberOfPairs();
                _globPairToLocaPair[indexZone] = _globPairToLocaPair[indexZone - 1] + nbPairs;
            }
        }

        // Get parameters from standard model
        auto paramToCopy = model->getFiniteElementDescriptor()->getParameters();
        paramToCopy->updateValuePointer();
        auto docu = paramToCopy->getInformationParameter();

        // Create LGRF object
        auto parameters = _fed->getParameters();
        parameters->allocate( 3 );
        ( *parameters )[0] = mesh->getName();
        ( *parameters )[1] = model->getName();
        ( *parameters )[2] = ( *paramToCopy )[2];
        parameters->setInformationParameter( docu );

        // Adapt FED
        CALLO_ADALIG_WRAP( _fed->getName() );
        bool l_calc_rigi = false;
        CALLO_INITEL( _fed->getName(), (ASTERLOGICAL *)&l_calc_rigi );

        // Final building
        _fed->build();
        _fedCreated = true;

    } else {
        _fedCreated = false;
    }

    // Clean
    listContElem.clear();
    slaveNodePaired.clear();
    slaveCellPaired.clear();

    CALL_JEDEMA();
};

VectorLong ContactPairing::getNumberOfIntersectionPoints() const {
    VectorLong returnValue;
    const ASTERINTEGER nbZoneCont = _contDefi->getNumberOfContactZones();
    for ( ASTERINTEGER iZone = 0; iZone < nbZoneCont; iZone++ ) {
        VectorLong nbInter;
        nbInter = ContactPairing::getNumberOfIntersectionPoints( iZone );
        returnValue.insert( returnValue.end(), nbInter.begin(), nbInter.end() );
    }
    return returnValue;
}

VectorLong ContactPairing::getNumberOfIntersectionPoints( ASTERINTEGER &indexZone ) const {

    VectorLong returnValue;

    returnValue = _contDefi->getMeshPairing( indexZone )->getNumberOfIntersectionPoints();

    return returnValue;
}

void ContactPairing::updateCoordinates( const FieldOnNodesRealPtr &disp ) {
#ifdef ASTER_HAVE_MPI
    if ( _mesh->getName() != disp->getMesh()->getName() && _mesh->isConnection() ) {
        const auto cMesh = std::dynamic_pointer_cast< ConnectionMesh >( _mesh );
        const auto &cCoords = _mesh->getCoordinates();
        const auto &localNum = cMesh->getNodesLocalNumbering();
        const auto &owners = cMesh->getNodesOwner();
        const auto &sFONRed = toSimpleFieldOnNodes( disp )->restrict( { "DX", "DY", "DZ" } );
        const auto &cmpNb = sFONRed->getComponents().size();

        const int rank = getMPIRank();
        const auto &nbNodes = localNum->size();
        cCoords->updateValuePointers();
        _currentCoordinates->updateValuePointers();
        VectorReal tmp( 3 * cMesh->getNumberOfNodes(), 0. ),
            res( 3 * cMesh->getNumberOfNodes(), 0. );
        for ( int i = 0; i < nbNodes; ++i ) {
            const auto &curOwner = ( *owners )[i];
            if ( curOwner == rank ) {
                const auto &curPos = ( *localNum )[i] - 1;
                for ( int j = 0; j < cmpNb; ++j ) {
                    tmp[i * 3 + j] = ( *cCoords )[i][j] + ( *sFONRed )[curPos * cmpNb + j];
                }
            }
        }
        AsterMPI::all_reduce( tmp, res, MPI_SUM );
        for ( int i = 0; i < nbNodes; ++i ) {
            for ( int j = 0; j < 3; ++j ) {
                ( *_currentCoordinates )[i][j] = res[i * 3 + j];
            }
        }
    } else {
#endif /* ASTER_HAVE_MPI */
        *_currentCoordinates = *( _mesh->getCoordinates() ) + *disp;
#ifdef ASTER_HAVE_MPI
    }
#endif /* ASTER_HAVE_MPI */
};
