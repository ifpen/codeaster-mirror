/**
 * @file DOFNumbering.cxx
 * @brief Implementation de DOFNumbering
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

#include "astercxx.h"

#include "Numbering/ParallelEquationNumbering.h"

#include "aster_fort_calcul.h"
#include "aster_fort_ds.h"

#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/Exceptions.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/Tools.h"

#ifdef ASTER_HAVE_MPI

ParallelEquationNumbering::ParallelEquationNumbering()
    : ParallelEquationNumbering( DataStructureNaming::getNewName() ) {};

ParallelEquationNumbering::ParallelEquationNumbering( const std::string &baseName )
    : EquationNumbering( baseName ),
      _localToGlobal( JeveuxVectorLong( getName() + ".NULG" ) ),
      _localToRank( JeveuxVectorLong( getName() + ".PDDL" ) ),
      _joints( nullptr ) {};

void ParallelEquationNumbering::_buildGlobal2LocalMap() {
    getLocalToGlobalMapping()->updateValuePointer();
    ASTERINTEGER nloc = getLocalToGlobalMapping()->size();

    _global2localMap.reserve( nloc );
    for ( auto j = 0; j < nloc; j++ )
        _global2localMap[( *getLocalToGlobalMapping() )[j]] = j;
};

bool ParallelEquationNumbering::useLagrangeDOF() const {

    bool local_answer = EquationNumbering::useLagrangeDOF();
    bool global_answer;

    AsterMPI::all_reduce( local_answer, global_answer, MPI_LOR );

    return global_answer;
};

VectorLong ParallelEquationNumbering::getPhysicalDOFs( const bool local ) const {
    auto dofInformation = this->getLagrangianInformations();
    dofInformation->updateValuePointer();
    ASTERINTEGER size = dofInformation->size();
    VectorLong physicalRows;
    ASTERINTEGER physicalIndicator;
    auto loc2glo = getLocalToGlobalMapping();
    loc2glo->updateValuePointer();

    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *dofInformation )[i];
        if ( physicalIndicator == 0 ) {
            if ( local )
                physicalRows.push_back( i );
            else {
                physicalRows.push_back( ( *loc2glo )[i] );
            }
        }
    }
    return physicalRows;
};

VectorLong ParallelEquationNumbering::getGhostDOFs( const bool local,
                                                    const bool lastLayerOnly ) const {
    auto localToRank = this->getLocalToRank();
    localToRank->updateValuePointer();
    VectorLong ghostRows;
    const auto rank = getMPIRank();
    ASTERINTEGER dofOwner;
    auto loc2glo = getLocalToGlobalMapping();
    if ( !local )
        loc2glo->updateValuePointer();
    for ( int i = 0; i < getNumberOfDOFs( true ); i++ ) {
        dofOwner = ( *localToRank )[i];
        if ( dofOwner != rank ) {
            if ( local )
                ghostRows.push_back( i );
            else {
                ghostRows.push_back( ( *loc2glo )[i] );
            }
        }
    }
    if ( lastLayerOnly ) {
        JeveuxVectorShort nodesids = _mesh->getLastGhostsLayer();
        VectorLong lastLayerDofs;
        for ( auto &dof : ghostRows ) {
            // Node id
            ASTERINTEGER nodeid = getNodeAndComponentIdFromDOF( dof ).first;
            auto it = std::find( nodesids.begin(), nodesids.end(), nodeid );
            if ( it != nodesids.end() )
                lastLayerDofs.push_back( dof );
        }
        return lastLayerDofs;
    } else {
        return ghostRows;
    }
};

VectorLong ParallelEquationNumbering::getNoGhostDOFs( const bool local ) const {
    auto localToRank = this->getLocalToRank();
    localToRank->updateValuePointer();
    const auto rank = getMPIRank();
    ASTERINTEGER dofOwner;
    VectorLong noGhostRows;
    auto loc2glo = getLocalToGlobalMapping();
    if ( !local )
        loc2glo->updateValuePointer();

    for ( int i = 0; i < getNumberOfDOFs( true ); i++ ) {
        dofOwner = ( *localToRank )[i];
        if ( dofOwner == rank ) {
            if ( local )
                noGhostRows.push_back( i );
            else {
                noGhostRows.push_back( ( *loc2glo )[i] );
            }
        }
    }
    return noGhostRows;
};

VectorLong ParallelEquationNumbering::getLagrangeDOFs( const bool local ) const {
    auto dofInformation = this->getLagrangianInformations();
    dofInformation->updateValuePointer();
    ASTERINTEGER size = dofInformation->size();
    VectorLong lagrangeRows;
    ASTERINTEGER physicalIndicator;
    auto loc2glo = getLocalToGlobalMapping();
    loc2glo->updateValuePointer();

    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *dofInformation )[i];
        if ( physicalIndicator != 0 ) {
            if ( local )
                lagrangeRows.push_back( i );
            else {
                lagrangeRows.push_back( ( *loc2glo )[i] );
            }
        }
    }
    return lagrangeRows;
};

std::map< ASTERINTEGER, VectorLong >
ParallelEquationNumbering::getDictOfLagrangeDOFs( const bool local ) const {
    std::map< ASTERINTEGER, VectorLong > ret;
    ret[1] = VectorLong();
    ret[2] = VectorLong();
    VectorLong &lag1 = ret[1], &lag2 = ret[2];
    auto lagrInfo = this->getLagrangianInformations();
    lagrInfo->updateValuePointer();
    ASTERINTEGER size = lagrInfo->size();
    ASTERINTEGER physicalIndicator;
    auto loc2glo = getLocalToGlobalMapping();
    loc2glo->updateValuePointer();

    for ( int i = 0; i < size; i++ ) {
        physicalIndicator = ( *lagrInfo )[i];
        if ( physicalIndicator == -1 ) {
            if ( local )
                lag1.push_back( i );
            else
                lag1.push_back( ( *loc2glo )[i] );
        }
        if ( physicalIndicator == -2 ) {
            if ( local )
                lag2.push_back( i );
            else
                lag2.push_back( ( *loc2glo )[i] );
        }
    }
    return ret;
};

std::string ParallelEquationNumbering::getComponentFromDOF( const ASTERINTEGER dof,
                                                            const bool local ) const {
    auto localrow = dof;
    if ( !local )
        localrow = globalToLocalDOF( dof );

    auto [nodeId, cmpName] = EquationNumbering::getNodeAndComponentFromDOF( localrow );
    return cmpName;
};

ASTERINTEGER
ParallelEquationNumbering::getNodeFromDOF( const ASTERINTEGER dof, const bool local ) const {
    auto localrow = dof;
    if ( !local )
        localrow = globalToLocalDOF( dof );

    auto [nodeId, cmpId] = this->getNodeAndComponentIdFromDOF( localrow, local );

    return nodeId;
};

bool ParallelEquationNumbering::isPhysicalDOF( const ASTERINTEGER dof, const bool local ) const {
    auto localrow = dof;
    if ( !local )
        localrow = globalToLocalDOF( dof );
    auto [nodeId, cmpId] = this->getNodeAndComponentIdFromDOF( localrow );

    return cmpId > 0;
};

ASTERINTEGER
ParallelEquationNumbering::getNumberOfDOFs( const bool local ) const {
    this->getNumberOfEquations()->updateValuePointer();
    if ( local )
        return ( *this->getNumberOfEquations() )[0];
    else
        return ( *this->getNumberOfEquations() )[1];
};

bool ParallelEquationNumbering::useSingleLagrangeDOF() const {
    bool local_answer = EquationNumbering::useSingleLagrangeDOF(), global_answer;
    AsterMPI::all_reduce( local_answer, global_answer, MPI_LAND );

    return global_answer;
};

ASTERINTEGER ParallelEquationNumbering::localToGlobalDOF( const ASTERINTEGER loc ) const {
    if ( loc < 0 or loc >= getNumberOfDOFs( true ) )
        throw std::out_of_range( "Invalid dof index" );
    auto loc2glo = getLocalToGlobalMapping();
    loc2glo->updateValuePointer();
    return ( *loc2glo )[loc];
}

ASTERINTEGER ParallelEquationNumbering::globalToLocalDOF( const ASTERINTEGER glob ) const {
    if ( _global2localMap.empty() )
        const_cast< ParallelEquationNumbering * >( this )->_buildGlobal2LocalMap();

    auto search = _global2localMap.find( glob );
    if ( search != _global2localMap.end() ) {
        return search->second;
    } else {
        auto rank = getMPIRank();
        throw std::out_of_range( "Global Dof number " + std::to_string( glob ) +
                                 " not found on rank " + std::to_string( rank ) );
        return -1;
    }
};

VectorPairLong ParallelEquationNumbering::getNodeAndComponentIdFromDOF( const bool local ) const {
    auto ret = EquationNumbering::getNodeAndComponentIdFromDOF();

    AS_ASSERT( _mesh->isParallel() );
    if ( !local ) {
        auto mapLG = _mesh->getLocalToGlobalNodeIds();
        mapLG->updateValuePointer();
        ASTERINTEGER nb_eq = ret.size();
        for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
            const auto node_id = ret[i_eq].first;
            if ( node_id >= 0 ) {
                ret[i_eq].first = ( *mapLG )[node_id];
            }
        }
    }
    return ret;
};

PairLong ParallelEquationNumbering::getNodeAndComponentIdFromDOF( const ASTERINTEGER dof,
                                                                  const bool local ) const {
    auto ret = EquationNumbering::getNodeAndComponentIdFromDOF( dof, local );

    AS_ASSERT( _mesh->isParallel() );
    if ( !local ) {
        auto mapLG = _mesh->getLocalToGlobalNodeIds();
        mapLG->updateValuePointer();
        auto node_id = ret.first;
        if ( node_id >= 0 ) {
            ret.first = ( *mapLG )[node_id];
        }
    }
    return ret;
};

VectorString ParallelEquationNumbering::getComponents() const {
    auto cmp_set = EquationNumbering::getComponents();

    VectorString cmp_glb;
    AsterMPI::all_gather( cmp_set, cmp_glb );

    return unique( cmp_glb );
};

SetLong ParallelEquationNumbering::getComponentsId() const {

    auto cmp_set = EquationNumbering::getComponentsId();

    SetLong cmp_glb;
    AsterMPI::all_gather( cmp_set, cmp_glb );

    return cmp_glb;
};

/**
 * @brief Maps between the names of components and their Ids
 */
std::map< std::string, ASTERINTEGER > ParallelEquationNumbering::getComponentsNameToId() const {

    auto cmp_std = EquationNumbering::getComponentsNameToId();
    std::map< std::string, ASTERINTEGER > cmp_glb;
    AsterMPI::all_gather( cmp_std, cmp_glb );

    return cmp_glb;
};

std::vector< std::pair< ASTERINTEGER, std::string > >
ParallelEquationNumbering::getNodeAndComponentFromDOF( const bool local ) const {
    auto ret = EquationNumbering::getNodeAndComponentFromDOF();

    AS_ASSERT( _mesh->isParallel() );
    if ( !local ) {
        auto mapLG = _mesh->getLocalToGlobalNodeIds();
        mapLG->updateValuePointer();
        ASTERINTEGER nb_eq = ret.size();
        for ( ASTERINTEGER i_eq = 0; i_eq < nb_eq; i_eq++ ) {
            auto node_id = ret[i_eq].first;
            ret[i_eq].first = ( *mapLG )[node_id];
        }
    }
    return ret;
};

std::pair< std::pair< VectorLong, VectorString >, VectorLong >
ParallelEquationNumbering::getDOFsWithDescription( const VectorString &cmps,
                                                   const VectorLong &nodes, const bool local,
                                                   const ASTERINTEGER same_rank ) const {

    VectorLong v_nodes;
    VectorString v_cmps;
    VectorLong dofs;

    std::set< ASTERINTEGER > set_nodes;
    if ( nodes.size() == 0 ) {
        auto group = _mesh->getNodes( std::string(), local, same_rank );
        set_nodes = toSet( group );
    } else {
        set_nodes = toSet( nodes );
    }

    auto idToName = getComponentsIdToName();

    std::set< std::string > set_cmps;
    if ( cmps.size() == 0 ) {
        set_cmps = toSet( getComponents() );
    } else {
        set_cmps = toSet( cmps );
    }

    v_cmps.reserve( set_cmps.size() * set_nodes.size() );
    v_nodes.reserve( set_cmps.size() * set_nodes.size() );
    dofs.reserve( set_cmps.size() * set_nodes.size() );

    const auto descr = getNodeAndComponentIdFromDOF( local );

    auto mapLG = getLocalToGlobalMapping();
    mapLG->updateValuePointer();
    for ( auto dof = 0; dof < descr.size(); ++dof ) {
        if ( set_nodes.find( descr[dof].first ) != set_nodes.end() ) {
            if ( descr[dof].second > 0 ) {
                auto cmp = idToName[descr[dof].second];
                if ( set_cmps.find( cmp ) != set_cmps.end() ) {
                    v_nodes.push_back( descr[dof].first );
                    v_cmps.push_back( cmp );
                    if ( local ) {
                        dofs.push_back( dof );
                    } else {
                        dofs.push_back( ( *mapLG )[dof] );
                    }
                }
            } else if ( descr[dof].second == 0 ) {
                const std::string cmp = "LAGR:MPC";
                if ( set_cmps.find( cmp ) != set_cmps.end() ) {
                    v_nodes.push_back( descr[dof].first );
                    v_cmps.push_back( cmp );
                    dofs.push_back( dof );
                }
            } else if ( descr[dof].second < 0 ) {
                const std::string cmp = "LAGR:" + idToName[-descr[dof].second];
                if ( set_cmps.find( cmp ) != set_cmps.end() ) {
                    v_nodes.push_back( descr[dof].first );
                    v_cmps.push_back( cmp );
                    if ( local ) {
                        dofs.push_back( dof );
                    } else {
                        dofs.push_back( ( *mapLG )[dof] );
                    }
                }
            }
        }
    }

    return std::make_pair( std::make_pair( v_nodes, v_cmps ), dofs );
};

std::pair< std::pair< VectorLong, VectorString >, VectorLong >
ParallelEquationNumbering::getDOFsWithDescription( const VectorString &cmps,
                                                   const VectorString &groupNames, const bool local,
                                                   const ASTERINTEGER same_rank ) const {
    VectorLong nodes;
    if ( groupNames.size() == 0 ) {
        nodes = _mesh->getNodes( std::string(), local, same_rank );
        // pour les MPC
        nodes.push_back( -1 );
    } else {
        nodes = _mesh->getNodes( groupNames, local, same_rank );
    }

    return this->getDOFsWithDescription( cmps, nodes, local, same_rank );
};

std::pair< ASTERINTEGER, std::string >
ParallelEquationNumbering::getNodeAndComponentFromDOF( const ASTERINTEGER dof,
                                                       const bool local ) const {
    auto localrow = dof;
    if ( !local ) {
        localrow = globalToLocalDOF( dof );
    }

    auto ret = EquationNumbering::getNodeAndComponentFromDOF( localrow );

    AS_ASSERT( _mesh->isParallel() );
    if ( !local ) {
        auto mapLG = _mesh->getLocalToGlobalNodeIds();
        mapLG->updateValuePointer();
        auto node_id = ret.first;
        ret.first = ( *mapLG )[node_id];
    }

    return ret;
};

bool ParallelEquationNumbering::build() {
    if ( !_joints ) {
        _informations->updateValuePointer();
        auto name = strip( ( *_informations )[4] );
        _joints = std::make_shared< Joints >( name );
    }
    return _joints->build();
};

#endif
