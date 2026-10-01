/**
 * @file ParallelContactFEDescriptor.cxx
 * @brief Implementation de ParallelContactFEDescriptor
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

#include "Modeling/ParallelContactFEDescriptor.h"

#include "aster_fort_utils.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/Tools.h"

#include <algorithm>
#include <limits>

#ifdef ASTER_HAVE_MPI

int componentNumber( ASTERINTEGER codedInt ) {
    int resu = 0;
    unsigned long one = 1;
    const unsigned long int uCodedInt = codedInt;
    constexpr int size = sizeof( unsigned long int );
    for ( int i = 0; i < size; ++i ) {
        if ( ( uCodedInt & one ) != 0 )
            ++resu;
        one = one << 1;
    }
    return resu;
};

ParallelContactFEDescriptor::ParallelContactFEDescriptor( const FiniteElementDescriptorPtr &FEDesc,
                                                          const ConnectionMeshPtr &mesh,
                                                          const ModelPtr &connectionModel,
                                                          const ModelPtr &model,
                                                          const VectorString &masters,
                                                          const VectorString &slaves )
    : ParallelContactFEDescriptor( ResultNaming::getNewResultName(), FEDesc, mesh, connectionModel,
                                   model, masters, slaves ) {};

ParallelContactFEDescriptor::ParallelContactFEDescriptor( const FiniteElementDescriptorPtr &FEDesc,
                                                          const ConnectionMeshPtr &mesh,
                                                          const ModelPtr &model )
    : ParallelContactFEDescriptor( ResultNaming::getNewResultName(), FEDesc, mesh, model ) {};

ParallelContactFEDescriptor::ParallelContactFEDescriptor( const std::string &name,
                                                          const FiniteElementDescriptorPtr &FEDesc,
                                                          const ConnectionMeshPtr &mesh,
                                                          const ModelPtr &model )
    : FiniteElementDescriptor( name, "LIGREL_CP", model->getMesh() ),
      _joints( std::make_shared< Joints >() ),
      _owner( JeveuxVectorLong( getName() + ".PNOE" ) ),
      _multiplicity( JeveuxVectorLong( getName() + ".MULT" ) ),
      _outerMultiplicity( JeveuxVectorLong( getName() + ".MUL2" ) ),
      _globalNumberingVirtualNodes( JeveuxVectorLong( getName() + ".NULG" ) ),
      _FEDesc( FEDesc ),
      _localVirtualIdToGlobalNodeId( JeveuxVectorLong( getName() + ".LOGL" ) ),
      _globalNodeIdToLocalVirtual( JeveuxVectorLong( getName() + ".GLLO" ) ),
      _cMeshName( JeveuxVectorChar24( getName() + ".CMNO" ) ) {
    _cMeshName->allocate( 1 );
    ( *_cMeshName )[0] = mesh->getName();
    struct DNInformations {
        int globalNum = 0;
        int localNum = -1;
        std::set< int > ranks;
        bool slaveDN = false;
    };
    const int rank = getMPIRank();
    const int nbProcs = getMPISize();

    const auto &owner = *( mesh->getNodesOwner() );
    _globalNodeIdToLocalVirtual->allocate( owner.size(), -1 );
    owner.updateValuePointer();
    const auto &explorer = FEDesc->getVirtualCellsExplorer();

    VectorInt virtualCellToKeep;
    VectorInt meshNodesToKeep( owner.size(), -1 );
    _virtualCellToKeep = VectorLong( explorer.size(), 1 );
    // Dans cette map, on stocke en face du numero physique :
    //  - le numero global du faux Lagrange
    //  - le numéro des processeurs qui possedent ce Lagrange
    std::map< int, DNInformations > addedVirtualNodes;
    std::map< int, int > vNodesSet;
    ASTERINTEGER globalVNodeNumber = 0, nbElemToKeep = 0, totalSizeToKeep = 0;
    ASTERINTEGER localVNodeNumber = 0;
    // 400 - maximum number of nodes of a cell
    VectorInt tmpVec( 400, 0 );
    for ( int num2 = 0; num2 < owner.size(); ++num2 ) {
        const int &curOwner = owner[num2];
        if ( curOwner == rank ) {
            meshNodesToKeep[num2] = rank;
        }
    }

    const auto &cMeshExplorer = mesh->getConnectivityExplorer();

    const auto &cellOwner = mesh->getCellsOwner();
    const auto &cellLocalId = mesh->getCellsLocalNumbering();
    const auto &liel = FEDesc->getListOfGroupsOfElementsExplorer();
    std::vector< VectorLong > toLiel( liel.size(), VectorLong() );
    _cellMatching = VectorOfVectorsLong( liel.size(), VectorLong() );
    _lielMatching = VectorLong( liel.size(), std::numeric_limits< ASTERINTEGER >::min() );

    // On commence par regarder les noeuds et elements qui doivent etre
    // gardes dans le nouveau ligrel
    for ( const auto meshElem : explorer ) {
        const auto numElem = meshElem.getId();
        bool keepElem = false;
        int pos = 0, curOwner = -1, oldCurOwner;
        VectorInt curOwners;
        const auto nbNodesElem = meshElem.getNumberOfNodes();
        for ( auto numNode : meshElem ) {
            // Si on est sur un noeud tardif...
            if ( numNode < 0 ) {
                throw std::runtime_error(
                    "No virtual nodes allowed in contact finite element descriptor" );
            }
            // Si on a un noeud physique...
            const ASTERINTEGER num2 = numNode - 1;
            curOwner = owner[num2];
            // ... et que le processeur courant le possede, on conserve ce noeud
            // et on conserve l'element
            if ( curOwner == rank ) {
                keepElem = true;
                meshNodesToKeep[num2] = rank;
            }
            tmpVec[pos] = num2;
            curOwners.push_back( curOwner );
            ++pos;
            oldCurOwner = curOwner;
        }
        // Si l'element est a conserver, on le note
        if ( keepElem ) {
            virtualCellToKeep.push_back( numElem );
            _virtualCellToKeep[numElem] = nbElemToKeep - 1;
            --nbElemToKeep;
            totalSizeToKeep += nbNodesElem;
            for ( ASTERINTEGER ino = 0; ino < nbNodesElem; ++ino ) {
                const auto num2 = tmpVec[ino];
                if ( addedVirtualNodes.count( num2 ) == 0 ) {
                    addedVirtualNodes[num2] = DNInformations();
                    if ( vNodesSet.count( num2 ) == 0 ) {
                        vNodesSet[num2] = globalVNodeNumber;
                        addedVirtualNodes[num2].globalNum = globalVNodeNumber;
                        ++globalVNodeNumber;
                    } else {
                        addedVirtualNodes[num2].globalNum = vNodesSet[num2];
                    }
                    if ( owner[num2] != rank ) {
                        addedVirtualNodes[num2].localNum = localVNodeNumber;
                        ++localVNodeNumber;
                    }
                }
                auto &rankSet = addedVirtualNodes[num2].ranks;
                for ( const auto &curOwner2 : curOwners ) {
                    rankSet.insert( curOwner2 );
                }
            }
        } else {
            for ( ASTERINTEGER ino = 0; ino < nbNodesElem; ++ino ) {
                const auto num2 = tmpVec[ino];
                if ( vNodesSet.count( num2 ) == 0 ) {
                    vNodesSet[num2] = globalVNodeNumber;
                    ++globalVNodeNumber;
                }
            }
        }
    }

    auto nbPartialNodes = mesh->getNumberOfNodes();
    const auto &localNum = mesh->getNodesLocalNumbering();
    localNum->updateValuePointer();
    const auto &pNodesComp = FEDesc->getPhysicalNodesComponentDescriptor();
    // Calcul du nombre d'entier code
    int nec = pNodesComp->size() / nbPartialNodes;
    // Si des noeuds tardifs sont a conserver, on peut creer le ligrel
    if ( localVNodeNumber > 0 ) {
        _owner->allocate( localVNodeNumber );
        _multiplicity->allocate( localVNodeNumber );
        _outerMultiplicity->allocate( localVNodeNumber );
        _globalNumberingVirtualNodes->allocate( localVNodeNumber );
    }
    int i = 0, j = 0;
    std::vector< VectorLong > toSend( nbProcs );
    std::vector< VectorLong > toReceive( nbProcs );
    for ( const auto &[numNodeM1, curInfo] : addedVirtualNodes ) {
        const auto &vNodeGlobNum = curInfo.globalNum;
        const auto &vNodeLocNum = curInfo.localNum;
        const auto &vNodeOwners = curInfo.ranks;
        const auto &curOwner = owner[numNodeM1];
        if ( curOwner == rank ) {
            const auto &idToSend = ( *localNum )[numNodeM1];
            for ( const auto &proc : vNodeOwners )
                if ( proc != rank )
                    toSend[proc].push_back( idToSend );
        } else {
            toReceive[curOwner].push_back( -vNodeLocNum - 1 );
        }
        if ( curInfo.localNum != -1 ) {
            ( *_globalNumberingVirtualNodes )[i] = -vNodeGlobNum - 1;
            ( *_owner )[i] = curOwner;
            ( *_outerMultiplicity )[i] = 1;
            ( *_multiplicity )[i] = 1;
            ++i;
        }
        ++j;
    }
    if ( i != localVNodeNumber )
        throw std::runtime_error( "Out of bound error" );

    // Creation des raccords
    const std::string cadre( "G" );
    const std::string error( "F" );
    VectorOfVectorsLong recv, send;
    recv.reserve( nbProcs );
    send.reserve( nbProcs );

    VectorLong joints;
    for ( i = 0; i < nbProcs; ++i ) {
        const auto taille1 = toSend[i].size();
        const auto taille2 = toReceive[i].size();
        if ( taille1 != 0 or taille2 != 0 ) {
            send.push_back( toSend[i] );
            recv.push_back( toReceive[i] );
            joints.push_back( i );
        }
    }

    _joints->setOppositeDomains( joints );
    _joints->setSendedElements( send );
    _joints->setReceivedElements( recv );

    if ( nbElemToKeep < 0 ) {
        // Allocation du .NEMA
        _virtualCellsDescriptor->allocate( -nbElemToKeep, totalSizeToKeep - nbElemToKeep );

        const auto &endIter = addedVirtualNodes.end();
        // Remplissage du .NEMA avec les elements tardifs a conserver
        for ( auto &numElem : virtualCellToKeep ) {
            const auto curElem = explorer[numElem];
            VectorLong toCopy;
            for ( const auto &numNode : curElem ) {
                const auto &numNodeM1 = numNode - 1;
                if ( numNode > 0 ) {
                    if ( owner[numNodeM1] == rank ) {
                        toCopy.push_back( ( *localNum )[numNodeM1] );
                    } else {
                        const auto &curIter = addedVirtualNodes.find( numNodeM1 );
                        if ( curIter == endIter ) {
                            throw std::runtime_error( "Incoherent situation" );
                        } else {
                            toCopy.push_back( -curIter->second.localNum - 1 );
                        }
                    }
                }
            }
            toCopy.push_back( explorer[numElem].getType() );
            _virtualCellsDescriptor->push_back( toCopy );
        }
    }

    int lielPos = 0;
    int nbCollObj = 0, totalCollSize = 0;
    for ( const auto &colObj : liel ) {
        bool addedElem = false;
        int curPos = 0;
        auto &curCM = _cellMatching[lielPos];
        for ( const auto &val : colObj ) {
            // Si on est sur une maille physique...
            if ( val >= 0 ) {
                throw std::runtime_error(
                    "No physical cells allowed in contact finite element descriptor" );
            }
            if ( _virtualCellToKeep[-val - 1] != 1 ) {
                curCM.push_back( curPos );
                toLiel[nbCollObj].push_back( _virtualCellToKeep[-val - 1] );
                addedElem = true;
                ++totalCollSize;
            }
            ++curPos;
        }
        if ( addedElem ) {
            const ASTERINTEGER &type = colObj.getType();
            toLiel[nbCollObj].push_back( type );
            _lielMatching[lielPos] = nbCollObj;
            ++nbCollObj;
        }
        ++lielPos;
    }

    if ( nbCollObj > 0 ) {
        _listOfGroupsOfElements->allocate( nbCollObj, totalCollSize + nbCollObj );
        for ( const auto &vec : toLiel ) {
            if ( vec.size() != 0 ) {
                _listOfGroupsOfElements->push_back( vec );
            }
        }
    }

    // Remplissage du .NBNO avec le nouveau nombre de noeuds tardifs
    _numberOfDelayedNumberedConstraintNodes->allocate( 1 );
    ( *_numberOfDelayedNumberedConstraintNodes )[0] = localVNodeNumber;

    const auto param = FEDesc->getParameters();
    param->updateValuePointer();
    // Creation du .LGRF en y mettant les noms du maillage et modele d'origine
    _parameters->allocate( 4 );
    const auto &pMesh = mesh->getParallelMesh();
    ( *_parameters )[0] = pMesh->getName();
    ( *_parameters )[1] = model->getName();
    ( *_parameters )[2] = ( *param )[2];
    ( *_parameters )[3] = _joints->getName();

    auto docu = FEDesc->getParameters()->getInformationParameter();
    _parameters->setInformationParameter( docu );
    /** @todo ajouter un assert sur le maillage sous-jacent au modele */

    // Creation du .PRNM sur les noeuds du getParallelMesh
    // Recopie en position locale
    _dofDescriptor->allocate( ( pMesh->getNumberOfNodes() ) * nec );
    for ( int i = 0; i < nbPartialNodes; ++i ) {
        if ( meshNodesToKeep[i] == rank )
            for ( int j = 0; j < nec; ++j ) {
                const int newPos = nec * ( ( *localNum )[i] - 1 ) + j;
                ( *_dofDescriptor )[newPos] = ( *pNodesComp )[i * nec + j];
            }
    }
    if ( localVNodeNumber > 0 ) {
        // Creation des .LGNS et .PRNS
        // Recopie des valeurs sur les noeuds tardifs du nouveau ligrel
        _virtualNodesNumbering->allocate( localVNodeNumber + 2, 0 );
        _dofOfDelayedNumberedConstraintNodes->allocate( localVNodeNumber * nec );
        _localVirtualIdToGlobalNodeId->allocate( localVNodeNumber );

        int checkCount = 0;
        for ( const auto &[numNodeM1, curInfo] : addedVirtualNodes ) {
            const auto &vNodeLocNum = curInfo.localNum;
            const auto &bSDN = curInfo.slaveDN;
            if ( vNodeLocNum == -1 )
                continue;
            ( *_localVirtualIdToGlobalNodeId )[vNodeLocNum] = numNodeM1 + 1;
            ( *_globalNodeIdToLocalVirtual )[numNodeM1] = vNodeLocNum + 1;
            if ( !bSDN )
                continue;
            ++checkCount;
            const auto &posInPrnm = numNodeM1 * nec;
            ASTERINTEGER nbCmp = 0;
            for ( int j = 0; j < nec; ++j ) {
                const auto newPos = vNodeLocNum * nec + j;
                ( *_dofOfDelayedNumberedConstraintNodes )[newPos] = ( *pNodesComp )[posInPrnm + j];
                nbCmp += componentNumber( ( *pNodesComp )[posInPrnm + j] );
            }
        }
    }

    build();

    // debugPrint();
};

ParallelContactFEDescriptor::ParallelContactFEDescriptor(
    const std::string &name, const FiniteElementDescriptorPtr &FEDesc,
    const ConnectionMeshPtr &mesh, const ModelPtr &connectionModel, const ModelPtr &model,
    const VectorString &masters, const VectorString &slaves )
    : FiniteElementDescriptor( name, "LIGREL_CP", model->getMesh() ),
      _joints( std::make_shared< Joints >() ),
      _owner( JeveuxVectorLong( getName() + ".PNOE" ) ),
      _multiplicity( JeveuxVectorLong( getName() + ".MULT" ) ),
      _outerMultiplicity( JeveuxVectorLong( getName() + ".MUL2" ) ),
      _globalNumberingVirtualNodes( JeveuxVectorLong( getName() + ".NULG" ) ),
      _FEDesc( FEDesc ),
      _localVirtualIdToGlobalNodeId( JeveuxVectorLong( getName() + ".LOGL" ) ),
      _globalNodeIdToLocalVirtual( JeveuxVectorLong( getName() + ".GLLO" ) ),
      _cMeshName( JeveuxVectorChar24( getName() + ".CMNO" ) ) {
    _cMeshName->allocate( 1 );
    ( *_cMeshName )[0] = mesh->getName();
    struct DNInformations {
        int globalNum = 0;
        int localNum = -1;
        std::set< int > ranks;
        bool slaveDN = false;
    };
    const int rank = getMPIRank();
    const int nbProcs = getMPISize();
    const auto &fEType = connectionModel->getFiniteElementType();
    fEType->updateValuePointer();
    const auto &lmasterCells = mesh->getCells( masters );
    std::map< ASTERINTEGER, VectorLong > cellsByType;
    for ( const auto &cellId : lmasterCells ) {
        const auto &curCellType = ( *fEType )[cellId];
        cellsByType[curCellType].push_back( cellId );
    }

    const auto &owner = *( mesh->getNodesOwner() );
    _globalNodeIdToLocalVirtual->allocate( owner.size(), -1 );
    owner.updateValuePointer();
    // const auto &explorer = FEDesc->getVirtualCellsExplorer();
    const auto &explorer = FEDesc->getListOfGroupsOfElementsExplorer();

    VectorInt virtualCellToKeep;
    VectorInt meshNodesToKeep( owner.size(), -1 );
    _virtualCellToKeep = VectorLong( explorer.size(), 1 );
    // Dans cette map, on stocke en face du numero physique :
    //  - le numero global du faux Lagrange
    //  - le numéro des processeurs qui possedent ce Lagrange
    std::map< int, DNInformations > addedVirtualNodes;
    std::map< int, int > vNodesSet;
    ASTERINTEGER globalVNodeNumber = 0, nbElemToKeep = 0, totalSizeToKeep = 0;
    ASTERINTEGER localVNodeNumber = 0;
    for ( int num2 = 0; num2 < owner.size(); ++num2 ) {
        const int &curOwner = owner[num2];
        if ( curOwner == rank ) {
            meshNodesToKeep[num2] = rank;
        }
    }

    const auto &cMeshExplorer = mesh->getConnectivityExplorer();

    const auto &cellOwner = mesh->getCellsOwner();
    const auto &cellLocalId = mesh->getCellsLocalNumbering();
    const auto &liel = FEDesc->getListOfGroupsOfElementsExplorer();
    int nbCollObj = 1, totalCollSize = 0;
    std::vector< VectorLong > toLiel( liel.size() + lmasterCells.size(), VectorLong() );
    _cellMatching = VectorOfVectorsLong( liel.size() + lmasterCells.size(), VectorLong() );

    VectorOfVectorsLong vCellToAdd;
    std::set< int > cellAldreadyAdd;
    const bool addSlaveNodes = true;

    for ( const auto &cellTypeIter : cellsByType ) {

        const ASTERINTEGER &type = cellTypeIter.first;

        const auto &cellVector = cellTypeIter.second;
        for ( const auto &cellId : cellVector ) {
            const auto &curCell = cMeshExplorer[cellId];
            bool nemaAdd = false;
            VectorLong cellToAddToNema;
            int nodeCount = 0;
            for ( const auto &nodeId : curCell ) {
                const auto &curOwner = owner[nodeId];
                if ( addedVirtualNodes.count( nodeId ) == 0 ) {
                    addedVirtualNodes[nodeId] = DNInformations();
                    addedVirtualNodes[nodeId].globalNum = globalVNodeNumber;
                    ++globalVNodeNumber;
                    if ( curOwner != rank ) {
                        addedVirtualNodes[nodeId].localNum = localVNodeNumber;
                        ++localVNodeNumber;
                    }
                    if ( addSlaveNodes )
                        addedVirtualNodes[nodeId].slaveDN = true;
                    auto &rankSet = addedVirtualNodes[nodeId].ranks;
                    for ( int i = 0; i < nbProcs; ++i ) {
                        rankSet.insert( i );
                    }
                }
                if ( curOwner != rank ) {
                    nemaAdd = true;
                    cellToAddToNema.push_back( -addedVirtualNodes[nodeId].localNum - 1 );
                } else {
                    cellToAddToNema.push_back( nodeId + 1 );
                }
                ++nodeCount;
            }
            if ( addSlaveNodes ) {
                if ( nemaAdd ) {
                    --nbElemToKeep;
                    totalSizeToKeep += nodeCount;
                    totalCollSize += 1;
                    cellToAddToNema.push_back( type );
                    vCellToAdd.push_back( cellToAddToNema );
                    toLiel[nbCollObj - 1].push_back( nbElemToKeep );
                } else {
                    toLiel[nbCollObj - 1].push_back( ( *cellLocalId )[cellId] );
                    totalCollSize += 1;
                }
                cellAldreadyAdd.insert( cellId );
            }
        }
        if ( addSlaveNodes )
            toLiel[nbCollObj - 1].push_back( type );
        nbCollObj++;
    }

    // On commence par regarder les noeuds et elements qui doivent etre
    // gardes dans le nouveau ligrel
    for ( const auto meshElem : explorer ) {
        const auto numElem = meshElem.getId();

        const auto &curCell = cMeshExplorer[numElem];
        bool nemaAdd = false;
        VectorLong cellToAddToNema;
        int nodeCount = 0;
        for ( const auto &nodeId : curCell ) {
            const auto &curOwner = owner[nodeId];
            if ( addedVirtualNodes.count( nodeId ) == 0 ) {
                addedVirtualNodes[nodeId] = DNInformations();
                addedVirtualNodes[nodeId].globalNum = globalVNodeNumber;
                ++globalVNodeNumber;
                if ( curOwner != rank ) {
                    addedVirtualNodes[nodeId].localNum = localVNodeNumber;
                    ++localVNodeNumber;
                }
                if ( addSlaveNodes )
                    addedVirtualNodes[nodeId].slaveDN = true;
                auto &rankSet = addedVirtualNodes[nodeId].ranks;
                for ( int i = 0; i < nbProcs; ++i ) {
                    rankSet.insert( i );
                }
            }
            if ( curOwner != rank ) {
                nemaAdd = true;
                cellToAddToNema.push_back( -addedVirtualNodes[nodeId].localNum - 1 );
            } else {
                cellToAddToNema.push_back( nodeId + 1 );
            }
            ++nodeCount;
        }
        if ( addSlaveNodes ) {
            if ( nemaAdd ) {
                --nbElemToKeep;
                totalSizeToKeep += nodeCount;
                totalCollSize += 1;
                cellToAddToNema.push_back( meshElem.getType() );
                vCellToAdd.push_back( cellToAddToNema );
                toLiel[nbCollObj - 1].push_back( nbElemToKeep );
            } else {
                toLiel[nbCollObj - 1].push_back( numElem + 1 );
                totalCollSize += 1;
            }
            cellAldreadyAdd.insert( numElem );
        }
        if ( addSlaveNodes )
            toLiel[nbCollObj - 1].push_back( meshElem.getType() );
        nbCollObj++;
    }

    auto nbPartialNodes = mesh->getNumberOfNodes();
    const auto &localNum = mesh->getNodesLocalNumbering();
    localNum->updateValuePointer();
    const auto &pNodesComp = FEDesc->getPhysicalNodesComponentDescriptor();
    // Calcul du nombre d'entier code
    int nec = pNodesComp->size() / nbPartialNodes;
    // Si des noeuds tardifs sont a conserver, on peut creer le ligrel
    if ( localVNodeNumber > 0 ) {
        _owner->allocate( localVNodeNumber );
        _multiplicity->allocate( localVNodeNumber );
        _outerMultiplicity->allocate( localVNodeNumber );
        _globalNumberingVirtualNodes->allocate( localVNodeNumber );
    }
    int i = 0, j = 0;
    std::vector< VectorLong > toSend( nbProcs );
    std::vector< VectorLong > toReceive( nbProcs );
    for ( const auto &[numNodeM1, curInfo] : addedVirtualNodes ) {
        const auto &vNodeGlobNum = curInfo.globalNum;
        const auto &vNodeLocNum = curInfo.localNum;
        const auto &vNodeOwners = curInfo.ranks;
        const auto &curOwner = owner[numNodeM1];
        if ( curOwner == rank ) {
            const auto &idToSend = ( *localNum )[numNodeM1];
            for ( const auto &proc : vNodeOwners )
                if ( proc != rank )
                    toSend[proc].push_back( idToSend );
        } else {
            toReceive[curOwner].push_back( -vNodeLocNum - 1 );
        }
        if ( curInfo.localNum != -1 ) {
            ( *_globalNumberingVirtualNodes )[i] = -vNodeGlobNum - 1;
            ( *_owner )[i] = curOwner;
            ( *_outerMultiplicity )[i] = 1;
            ( *_multiplicity )[i] = 1;
            ++i;
        }
        ++j;
    }
    if ( i != localVNodeNumber )
        throw std::runtime_error( "Out of bound error" );

    // Creation des raccords
    const std::string cadre( "G" );
    const std::string error( "F" );
    VectorOfVectorsLong recv, send;
    recv.reserve( nbProcs );
    send.reserve( nbProcs );

    VectorLong joints;
    for ( i = 0; i < nbProcs; ++i ) {
        const auto taille1 = toSend[i].size();
        const auto taille2 = toReceive[i].size();
        if ( taille1 != 0 or taille2 != 0 ) {
            send.push_back( toSend[i] );
            recv.push_back( toReceive[i] );
            joints.push_back( i );
        }
    }

    _joints->setOppositeDomains( joints );
    _joints->setSendedElements( send );
    _joints->setReceivedElements( recv );

    if ( nbElemToKeep < 0 ) {
        // Allocation du .NEMA
        _virtualCellsDescriptor->allocate( -nbElemToKeep, totalSizeToKeep - nbElemToKeep );

        int cmpt = 0;
        for ( const auto &vCell : vCellToAdd ) {
            _virtualCellsDescriptor->push_back( vCell );
            ++cmpt;
        }

        const auto &endIter = addedVirtualNodes.end();
        // Remplissage du .NEMA avec les elements tardifs a conserver
        for ( auto &numElem : virtualCellToKeep ) {
            const auto curElem = explorer[numElem];
            VectorLong toCopy;
            for ( const auto &numNode : curElem ) {
                const auto &numNodeM1 = numNode - 1;
                if ( numNode > 0 ) {
                    if ( owner[numNodeM1] == rank ) {
                        toCopy.push_back( ( *localNum )[numNodeM1] );
                    } else {
                        const auto &curIter = addedVirtualNodes.find( numNodeM1 );
                        if ( curIter == endIter ) {
                            throw std::runtime_error( "Incoherent situation" );
                        } else {
                            toCopy.push_back( -curIter->second.localNum - 1 );
                        }
                    }
                }
            }
            toCopy.push_back( explorer[numElem].getType() );
            _virtualCellsDescriptor->push_back( toCopy );
            ++cmpt;
        }
    }

    int lielPos = 0;
    for ( const auto &colObj : liel ) {
        bool addedElem = false;
        int curPos = 0;
        auto &curCM = _cellMatching[lielPos];
        for ( const auto &val : colObj ) {
            if ( val < 0 ) {
                if ( _virtualCellToKeep[-val - 1] != 1 ) {
                    curCM.push_back( curPos );
                    toLiel[nbCollObj - 1].push_back( _virtualCellToKeep[-val - 1] );
                    addedElem = true;
                    ++totalCollSize;
                }
            } else {
                if ( cellAldreadyAdd.count( val - 1 ) == 0 ) {
                    if ( ( *cellOwner )[val - 1] == rank ) {
                        toLiel[nbCollObj - 1].push_back( ( *cellLocalId )[val - 1] );
                        addedElem = true;
                        ++totalCollSize;
                        curCM.push_back( curPos );
                    }
                }
            }
            ++curPos;
        }
        ++lielPos;
        if ( addedElem ) {
            const ASTERINTEGER &type = colObj.getType();
            toLiel[nbCollObj - 1].push_back( type );
            ++nbCollObj;
        }
    }

    _listOfGroupsOfElements->allocate( nbCollObj, totalCollSize + nbCollObj );
    for ( const auto &vec : toLiel ) {
        if ( vec.size() != 0 ) {
            _listOfGroupsOfElements->push_back( vec );
        }
    }

    // Remplissage du .NBNO avec le nouveau nombre de noeuds tardifs
    _numberOfDelayedNumberedConstraintNodes->allocate( 1 );
    ( *_numberOfDelayedNumberedConstraintNodes )[0] = localVNodeNumber;

    const auto param = FEDesc->getParameters();
    param->updateValuePointer();
    // Creation du .LGRF en y mettant les noms du maillage et modele d'origine
    _parameters->allocate( 4 );
    const auto &pMesh = mesh->getParallelMesh();
    ( *_parameters )[0] = pMesh->getName();
    ( *_parameters )[1] = model->getName();
    ( *_parameters )[2] = ( *param )[2];
    ( *_parameters )[3] = _joints->getName();

    auto docu = FEDesc->getParameters()->getInformationParameter();
    _parameters->setInformationParameter( docu );
    /** @todo ajouter un assert sur le maillage sous-jacent au modele */

    // Creation du .PRNM sur les noeuds du getParallelMesh
    // Recopie en position locale
    _dofDescriptor->allocate( ( pMesh->getNumberOfNodes() ) * nec );
    for ( int i = 0; i < nbPartialNodes; ++i ) {
        if ( meshNodesToKeep[i] == rank )
            for ( int j = 0; j < nec; ++j ) {
                const int newPos = nec * ( ( *localNum )[i] - 1 ) + j;
                ( *_dofDescriptor )[newPos] = ( *pNodesComp )[i * nec + j];
            }
    }
    if ( localVNodeNumber > 0 ) {
        // Creation des .LGNS et .PRNS
        // Recopie des valeurs sur les noeuds tardifs du nouveau ligrel
        _virtualNodesNumbering->allocate( localVNodeNumber + 2, 0 );
        _dofOfDelayedNumberedConstraintNodes->allocate( localVNodeNumber * nec );
        _localVirtualIdToGlobalNodeId->allocate( localVNodeNumber );

        int checkCount = 0;
        for ( const auto &[numNodeM1, curInfo] : addedVirtualNodes ) {
            const auto &vNodeLocNum = curInfo.localNum;
            const auto &bSDN = curInfo.slaveDN;
            if ( vNodeLocNum == -1 )
                continue;
            ( *_localVirtualIdToGlobalNodeId )[vNodeLocNum] = numNodeM1 + 1;
            ( *_globalNodeIdToLocalVirtual )[numNodeM1] = vNodeLocNum + 1;
            if ( !bSDN )
                continue;
            ++checkCount;
            const auto &posInPrnm = numNodeM1 * nec;
            ASTERINTEGER nbCmp = 0;
            for ( int j = 0; j < nec; ++j ) {
                const auto newPos = vNodeLocNum * nec + j;
                ( *_dofOfDelayedNumberedConstraintNodes )[newPos] = ( *pNodesComp )[posInPrnm + j];
                nbCmp += componentNumber( ( *pNodesComp )[posInPrnm + j] );
            }
        }
    }

    build();
};

#endif /* ASTER_HAVE_MPI */
