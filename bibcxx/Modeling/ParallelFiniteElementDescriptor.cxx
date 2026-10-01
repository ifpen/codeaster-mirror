/**
 * @file ParallelFiniteElementDescriptor.cxx
 * @brief Implementation de ParallelFiniteElementDescriptor
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

#include "Modeling/ParallelFiniteElementDescriptor.h"

#include "aster_fort_utils.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

#include <algorithm>

#ifdef ASTER_HAVE_MPI

ParallelFiniteElementDescriptor::ParallelFiniteElementDescriptor(
    const std::string &name, const FiniteElementDescriptorPtr &FEDesc,
    const ConnectionMeshPtr &mesh, const ModelPtr &model )
    : FiniteElementDescriptor( name, model->getMesh() ),
      _joints( std::make_shared< Joints >() ),
      _owner( JeveuxVectorLong( getName() + ".PNOE" ) ),
      _multiplicity( JeveuxVectorLong( getName() + ".MULT" ) ),
      _outerMultiplicity( JeveuxVectorLong( getName() + ".MUL2" ) ),
      _globalNumberingVirtualNodes( JeveuxVectorLong( getName() + ".NULG" ) ) {
    const int rank = getMPIRank();
    const int nbProcs = getMPISize();

    const auto &owner = *( mesh->getNodesOwner() );
    const auto &explorer = FEDesc->getVirtualCellsExplorer();

    VectorInt virtualCellToKeep;
    VectorInt meshNodesToKeep( owner.size(), -1 );
    ASTERINTEGER nbOldVirtualNodes = FEDesc->getNumberOfVirtualNodes();
    _virtualCellToKeep = VectorLong( explorer.size(), 1 );
    VectorInt virtualNodesToKeep( nbOldVirtualNodes, -1 );
    VectorInt virtualNodesNumbering( nbOldVirtualNodes, -1 );
    VectorInt virtualNodesMult( nbOldVirtualNodes, 0 );
    VectorInt virtualNodesOuterMult( nbOldVirtualNodes, 0 );
    VectorInt virtualNodesOwner( nbOldVirtualNodes, -1 );
    VectorInt nbOwnedVirtualNodes( nbProcs, 0 );
    std::vector< std::set< int > > sharedVirtualNodes( nbOldVirtualNodes );
    ASTERINTEGER nbVirtualNodes = 0, nbElemToKeep = 0, totalSizeToKeep = 0;
    // On commence par regarder les noeuds et elements qui doivent etre
    // gardes dans le nouveau ligrel
    for ( const auto meshElem : explorer ) {
        const auto numElem = meshElem.getId();
        bool keepElem = false;
        int pos = 0, curOwner = -1;
        for ( auto numNode : meshElem ) {
            if ( pos == 0 && numNode < 0 )
                throw std::runtime_error( "First node is assumed to be a physical node" );
            // Si on a un noeud physique...
            if ( numNode > 0 ) {
                const ASTERINTEGER num2 = numNode - 1;
                curOwner = owner[num2];
                // ... et que le processeur courant le possede, on conserve ce noeud
                // et on conserve l'element
                if ( curOwner == rank ) {
                    keepElem = true;
                    meshNodesToKeep[num2] = rank;
                    ++totalSizeToKeep;
                }
            }
            // Si on est sur un noeud tardif...
            if ( numNode < 0 ) {
                ++virtualNodesMult[-numNode - 1];
                if ( keepElem ) {
                    // ... et qu'il faut conserver l'element
                    // Alors on conserve les noeuds tardifs aussi
                    if ( virtualNodesToKeep[-numNode - 1] == -1 ) {
                        virtualNodesNumbering[-numNode - 1] = nbVirtualNodes;
                        ++nbVirtualNodes;
                    }
                    ++totalSizeToKeep;
                    virtualNodesToKeep[-numNode - 1] = rank;
                } else
                    ++virtualNodesOuterMult[-numNode - 1];
                // On cherche a equilibrer la charge des noeuds tardifs
                auto curOwner2 = virtualNodesOwner[-numNode - 1];
                if ( curOwner2 == -1 ) {
                    virtualNodesOwner[-numNode - 1] = curOwner;
                    ++nbOwnedVirtualNodes[curOwner];
                } else {
                    if ( nbOwnedVirtualNodes[curOwner2] > nbOwnedVirtualNodes[curOwner] ) {
                        virtualNodesOwner[-numNode - 1] = curOwner;
                        ++nbOwnedVirtualNodes[curOwner];
                        --nbOwnedVirtualNodes[curOwner2];
                    }
                }
                // On note tous les procs qui possèdent un noeud tardif
                sharedVirtualNodes[-numNode - 1].insert( curOwner );
            }
            ++pos;
        }
        // Si l'element est a conserver, on le note
        if ( keepElem ) {
            virtualCellToKeep.push_back( numElem );
            _virtualCellToKeep[numElem] = nbElemToKeep - 1;
            --nbElemToKeep;
        }
    }

    auto nbPartialNodes = mesh->getNumberOfNodes();
    const auto &localNum = mesh->getNodesLocalNumbering();
    const auto &pNodesComp = FEDesc->getPhysicalNodesComponentDescriptor();
    // Calcul du nombre d'entier code
    int nec = pNodesComp->size() / nbPartialNodes;
    // Si des noeuds tardifs sont a conserver, on peut creer le ligrel
    if ( nbVirtualNodes > 0 ) {
        _owner->allocate( nbVirtualNodes );
        _multiplicity->allocate( nbVirtualNodes );
        _outerMultiplicity->allocate( nbVirtualNodes );
        int i = 0, nbJoins = 0, j = 0;
        std::vector< VectorLong > toSend( nbProcs );
        std::vector< VectorLong > toReceive( nbProcs );
        for ( const auto &curSet : sharedVirtualNodes ) {
            if ( virtualNodesToKeep[i] == rank ) {
                if ( virtualNodesOwner[i] == rank ) {
                    for ( const auto &proc : curSet )
                        if ( proc != rank )
                            toSend[proc].push_back( -virtualNodesNumbering[i] - 1 );
                } else {
                    toReceive[virtualNodesOwner[i]].push_back( -virtualNodesNumbering[i] - 1 );
                }
                ( *_owner )[j] = virtualNodesOwner[i];
                ( *_outerMultiplicity )[j] = virtualNodesOuterMult[i];
                ( *_multiplicity )[j] = virtualNodesMult[i];
                ++j;
            }
            ++i;
        }
        if ( j != nbVirtualNodes )
            throw std::runtime_error( "Out of bound error" );

        // Creation numérotation globale pour les noeuds tardifs
        _globalNumberingVirtualNodes->allocate( nbVirtualNodes );
        int count = 0;
        for ( int iNode = 0; iNode < nbOldVirtualNodes; iNode++ ) {
            if ( virtualNodesNumbering[iNode] != -1 ) {
                ( *_globalNumberingVirtualNodes )[virtualNodesNumbering[iNode]] = -( iNode + 1 );
                count++;
            }
        }
        AS_ASSERT( nbVirtualNodes == count );

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

        // Allocation du .NEMA
        _virtualCellsDescriptor->allocate( -nbElemToKeep, totalSizeToKeep - nbElemToKeep );

        // Remplissage du .NEMA avec les elements tardifs a conserver
        for ( auto &numElem : virtualCellToKeep ) {
            const auto curElem = explorer[numElem];
            VectorLong toCopy;
            for ( const auto &numNode : curElem ) {
                if ( numNode > 0 ) {
                    toCopy.push_back( ( *localNum )[numNode - 1] );
                } else {
                    toCopy.push_back( -virtualNodesNumbering[-numNode - 1] - 1 );
                }
            }
            toCopy.push_back( curElem.getType() );
            _virtualCellsDescriptor->push_back( toCopy );
        }

        const auto &liel = FEDesc->getListOfGroupsOfElementsExplorer();
        int nbCollObj = 0, totalCollSize = 0;
        std::vector< VectorLong > toLiel( liel.size(), VectorLong() );
        ASTERINTEGER type = 0;
        nbCollObj = 1;
        for ( const auto &colObj : liel ) {
            bool addedElem = false;
            for ( const auto &val : colObj ) {
                if ( _virtualCellToKeep[-val - 1] != 1 ) {
                    toLiel[nbCollObj - 1].push_back( _virtualCellToKeep[-val - 1] );
                    addedElem = true;
                    ++totalCollSize;
                }
            }
            if ( addedElem ) {
                type = colObj.getType();
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
    } else {
        _listOfGroupsOfElements->allocate( 1, 1 );
        VectorLong joints;
        _joints->setOppositeDomains( joints );
    }

    // Remplissage du .NBNO avec le nouveau nombre de noeuds tardifs
    _numberOfDelayedNumberedConstraintNodes->allocate( 1 );
    ( *_numberOfDelayedNumberedConstraintNodes )[0] = nbVirtualNodes;

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
    if ( nbVirtualNodes > 0 ) {
        // Creation des .LGNS et .PRNM
        // Recopie des valeurs sur les noeuds tardifs du nouveau ligrel
        _virtualNodesNumbering->allocate( nbVirtualNodes + 2 );
        _dofOfDelayedNumberedConstraintNodes->allocate( nbVirtualNodes * nec );
        const auto &dNodesComp = FEDesc->getVirtualNodesComponentDescriptor();
        dNodesComp->updateValuePointer();
        const auto &numbering = FEDesc->getVirtualNodesNumbering();
        numbering->updateValuePointer();

        for ( ASTERINTEGER num = 0; num < nbOldVirtualNodes; ++num ) {
            if ( virtualNodesToKeep[num] == rank ) {
                const ASTERINTEGER newNum = virtualNodesNumbering[num];
                for ( int j = 0; j < nec; ++j ) {
                    const int newPos = nec * newNum + j;
                    ( *_dofOfDelayedNumberedConstraintNodes )[newPos] =
                        ( *dNodesComp )[num * nec + j];
                }
                ( *_virtualNodesNumbering )[newNum] = ( *numbering )[num];
            }
        }
    }

    build();
};

ParallelFiniteElementDescriptor::ParallelFiniteElementDescriptor( const std::string &name,
                                                                  const std::string &jName,
                                                                  const BaseMeshPtr &mesh )
    : FiniteElementDescriptor( name, mesh ),
      _joints( std::make_shared< Joints >( jName ) ),
      _owner( JeveuxVectorLong( getName() + ".PNOE" ) ),
      _multiplicity( JeveuxVectorLong( getName() + ".MULT" ) ),
      _outerMultiplicity( JeveuxVectorLong( getName() + ".MUL2" ) ),
      _globalNumberingVirtualNodes( JeveuxVectorLong( getName() + ".NULG" ) ) {};

#endif /* ASTER_HAVE_MPI */
