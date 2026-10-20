#ifndef SIMPLEFIELDONNODES_H_
#define SIMPLEFIELDONNODES_H_

/**
 * @file SimpleFieldOnNodes.h
 * @brief Fichier entete de la classe SimpleFieldOnNodes
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

#include "aster_fort_ds.h"
#include "aster_fort_utils.h"

#include "DataFields/DataField.h"
#include "MemoryManager/JeveuxVector.h"
#include "MemoryManager/NumpyAccess.h"
#include "Meshes/BaseMesh.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

/**
 * @class SimpleFieldOnNodes
 * @brief Cette classe template permet de definir un champ aux noeuds Aster
 * @author Nicolas Sellenet
 */
template < class ValueType >
class SimpleFieldOnNodes : public DataField {
  private:
    /** @brief Vecteur Jeveux '.CNSK' */
    JeveuxVectorChar8 _descriptor;
    /** @brief Vecteur Jeveux '.CNSD' */
    JeveuxVectorLong _size;
    /** @brief Vecteur Jeveux '.CNSC' */
    JeveuxVectorChar8 _component;
    /** @brief Vecteur Jeveux '.CNSV' */
    JeveuxVector< ValueType > _values;
    /** @brief Vecteur Jeveux '.CNSL' */
    JeveuxVectorLogical _allocated;
    /** @brief Nombre de noeuds */
    ASTERINTEGER _nbNodes;
    /** @brief Nombre de composantes */
    ASTERINTEGER _nbComp;
    /** @brief Mesh */
    BaseMeshPtr _mesh;

    std::map< std::string, ASTERINTEGER > _name2Index;

    void _buildComponentsName2Index() {
        if ( _name2Index.empty() ) {

            auto nbCmp = this->getNumberOfComponents();
            for ( ASTERINTEGER i = 0; i < nbCmp; i++ ) {
                _name2Index[this->getComponent( i )] = i;
            }
        }
    }

    /**
     * Functions to check an out-of-range condition
     */
    void _checkNodeOOR( const ASTERINTEGER &ino ) const {
        ASTERINTEGER nbNodes = this->getNumberOfNodes();
        if ( ino < 0 || ino >= nbNodes ) {
            throw std::runtime_error( "Node index " + std::to_string( ino ) +
                                      " is out of range [0, " + std::to_string( nbNodes - 1 ) +
                                      "]" );
        };
    }

    void _checkCmpOOR( const ASTERINTEGER &icmp ) const {
        ASTERINTEGER ncmp = this->getNumberOfComponents();
        if ( icmp < 0 || icmp >= ncmp ) {
            throw std::runtime_error( "Component " + std::to_string( icmp ) +
                                      " is out of range [0, " + std::to_string( ncmp - 1 ) + "]" );
        }
    }

    void _checkSize( const ASTERINTEGER &ino, const ASTERINTEGER &icmp ) const {
        if ( this->getNumberOfNodes() == 0 || this->getNumberOfComponents() == 0 )
            throw std::runtime_error( "First call of updateValuePointers is mandatory" );

        this->_checkNodeOOR( ino );
        this->_checkCmpOOR( icmp );
    }

    void _existsValue() const {
        if ( this->getNumberOfNodes() == 0 || this->getNumberOfComponents() == 0 )
            throw std::runtime_error( "First call of updateValuePointers is mandatory" );

        auto size = _allocated->size();
        bool isAlloc = false;
        for ( ASTERINTEGER i = 0; i < size; i++ ) {
            if ( ( *_allocated )[i] ) {
                isAlloc = true;
                break;
            }
        }

        if ( !isAlloc ) {
            throw std::runtime_error( "No values affected in the field" );
        }
    }

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux noeuds
     */
    SimpleFieldOnNodes( const std::string name, const BaseMeshPtr mesh = nullptr )
        : DataField( name, "CHAM_NO_S" ),
          _descriptor( JeveuxVectorChar8( getName() + ".CNSK" ) ),
          _size( JeveuxVectorLong( getName() + ".CNSD" ) ),
          _component( JeveuxVectorChar8( getName() + ".CNSC" ) ),
          _values( JeveuxVector< ValueType >( getName() + ".CNSV" ) ),
          _allocated( JeveuxVectorLogical( getName() + ".CNSL" ) ),
          _nbNodes( 0 ),
          _nbComp( 0 ),
          _mesh( mesh ) {};

  public:
    /**
     * @typedef SimpleFieldOnNodesPtr
     * @brief Pointeur intelligent vers un SimpleFieldOnNodes
     */
    typedef std::shared_ptr< SimpleFieldOnNodes > SimpleFieldOnNodesPtr;

    /**
     * @brief Constructeur

     */
    SimpleFieldOnNodes() = delete;

    SimpleFieldOnNodes( const BaseMeshPtr mesh )
        : SimpleFieldOnNodes( DataStructureNaming::getNewName( 19 ), mesh ) {};

    SimpleFieldOnNodes( const BaseMeshPtr mesh, const std::string quantity,
                        const VectorString &comp, bool zero = false )
        : SimpleFieldOnNodes( mesh ) {
        this->allocate( quantity, comp, zero );
    }

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    SimpleFieldOnNodes( const py::tuple &tup )
        : SimpleFieldOnNodes( tup[0].cast< std::string >(), tup[1].cast< BaseMeshPtr >() ) {
        build();
    };

    py::tuple _getState() const { return py::make_tuple( getName(), getMesh() ); };

    BaseMeshPtr getMesh() const { return _mesh; };

    void allocate( const std::string quantity, const VectorString &comp, bool zero = false ) {
        ASTERINTEGER nbComp = comp.size();
        std::string base = "G";

        char *tabNames = vectorStringAsFStrArray( comp, 8 );

        CALL_CNSCRE_WRAP( _mesh->getName().c_str(), quantity.c_str(), &nbComp, tabNames,
                          base.c_str(), getName().c_str(), (ASTERLOGICAL *)&zero );

        FreeStr( tabNames );

        build();
    }

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    inline ValueType &operator[]( const ASTERINTEGER &i ) { return _values->operator[]( i ); };

    inline const ValueType &operator[]( const ASTERINTEGER &i ) const {
        return _values->operator[]( i );
    };

    ValueType &operator()( const ASTERINTEGER &ino, const ASTERINTEGER &icmp ) {
#ifdef ASTER_DEBUG_CXX_LOW_LEVEL
        _checkSize( ino, icmp );
#endif

        const ASTERINTEGER position = ino * this->getNumberOfComponents() + icmp;

        ( *_allocated )[position] = true;
        return this->operator[]( position );
    };

    const ValueType &operator()( const ASTERINTEGER &ino, const ASTERINTEGER &icmp ) const {
#ifdef ASTER_DEBUG_CXX_LOW_LEVEL

        if ( !this->hasValue( ino, icmp ) ) {
            AS_ABORT( "DEBUG: Position (" + std::to_string( ino ) + ", " + std::to_string( icmp ) +
                      ") is valid but not allocated!" )
        };
#endif

        const ASTERINTEGER position = ino * this->getNumberOfComponents() + icmp;

        return this->operator[]( position );
    };

    ValueType &operator()( const ASTERINTEGER &ino, const std::string &cmp ) {
        auto icmp = _name2Index.at( cmp );
        return this->operator()( ino, icmp );
    };

    const ValueType &operator()( const ASTERINTEGER &ino, const std::string &cmp ) const {
        auto icmp = _name2Index.at( cmp );
        return this->operator()( ino, icmp );
    };

    bool hasValue( const ASTERINTEGER &ino, const ASTERINTEGER &icmp ) const {

        if ( ino < 0 || this->getNumberOfNodes() == 0 || ino >= this->getNumberOfNodes() ) {
            return false;
        }

        if ( icmp < 0 || this->getNumberOfComponents() == 0 ||
             icmp >= this->getNumberOfComponents() ) {
            return false;
        }

#ifdef ASTER_DEBUG_CXX_LOW_LEVEL
        _checkSize( ino, icmp );
#endif

        const ASTERINTEGER position = ino * this->getNumberOfComponents() + icmp;

        return ( *_allocated )[position];
    };

    bool hasValue( const ASTERINTEGER &ino, const std::string &cmp ) const {
        auto icmp = _name2Index.at( cmp );

        return this->hasValue( ino, icmp );
    };

    /**
     * @brief Get number of components
     */
    inline ASTERINTEGER getNumberOfComponents() const { return _nbComp; }

    /**
     * @brief Get number of nodes
     */
    inline ASTERINTEGER getNumberOfNodes() const { return _nbNodes; }

    /**
     * @brief Return a pointer to the vector of data
     */
    const ValueType *getDataPtr() const { return _values->getDataPtr(); }

    /**
     * @brief Return value vector
     */
    JeveuxVector< ValueType > getValues() const { return _values; }

    /**
     * @brief Return logical value vector
     */
    JeveuxVectorLogical getLogicalValues() const { return _allocated; }

    /**
     * @brief Get the name of the i-th component
     */
    std::string getComponent( const ASTERINTEGER &i ) const {

        if ( i < 0 || i >= _nbComp ) {
            throw std::runtime_error( "Out of range" );
        };

        std::string name = strip( ( *_component )[i].toString() );
        return name;
    };

    /**
     * @Brief Get the names of all the components
     */
    VectorString getComponents() const {

        ASTERINTEGER size = this->getNumberOfComponents();
        VectorString names;
        names.reserve( size );
        for ( ASTERINTEGER i = 0; i < size; i++ ) {
            names.push_back( this->getComponent( i ) );
        }
        return names;
    }

    /**
     * @brief Maps between name of components and the nimber
     */
    std::map< std::string, ASTERINTEGER > getComponentsName2Index() const {

        if ( _name2Index.empty() ) {
            raiseAsterError( "getComponentsName2Index: build field before" );
        }

        return _name2Index;
    };

    /**
     * @brief Get physical quantity
     */
    std::string getPhysicalQuantity() const { return strip( ( *_descriptor )[1].toString() ); }

    /**
     * @brief Get field location
     */
    std::string getLocalization() const { return strip( ( *_descriptor )[2].toString() ); }

    bool hasComponent( const std::string &cmp ) const {
        if ( _name2Index.count( cmp ) > 0 ) {
            return true;
        }
        return false;
    }

    void setValues( const VectorLong &nodes, const VectorString &cmps,
                    const std::vector< ValueType > &values ) {

        AS_ASSERT( nodes.size() == cmps.size() );
        AS_ASSERT( nodes.size() == values.size() );

        this->updateValuePointers();

        const int size = values.size();
        for ( int i = 0; i < size; i++ ) {
            ( *this )( nodes[i], cmps[i] ) = values[i];
        }
    }

    void setValues( const std::vector< std::vector< ValueType > > &values ) {

        AS_ASSERT( values.size() == this->getNumberOfNodes() );
        this->updateValuePointers();

        for ( int ino = 0; ino < this->getNumberOfNodes(); ino++ ) {
            auto &row = values[ino];
#ifdef ASTER_DEBUG_CXX_LOW_LEVEL
            if ( this->getNumberOfComponents() != row.size() )
                throw std::runtime_error( "Incompatible size" );
#endif
            for ( int icmp = 0; icmp < this->getNumberOfComponents(); icmp++ ) {
                ( *this )( ino, icmp ) = row[icmp];
            }
        }
    }

    void setValues( const std::vector< ValueType > &values ) {

        AS_ASSERT( values.size() == this->getNumberOfNodes() * this->getNumberOfComponents() );
        _allocated->assign( true );
        *_values = values;
    }

    void setValues( const ValueType value ) {

        _allocated->assign( true );
        _values->assign( value );
    }

    void setValues( const std::map< std::string, ValueType > &values, const VectorLong &nodes ) {
        this->updateValuePointers();

        for ( auto &[cmp, val] : values ) {
            if ( this->hasComponent( cmp ) ) {
                const auto icmp = _name2Index[cmp];
                for ( auto &node : nodes ) {
                    ( *this )( node, cmp ) = val;
                }
            }
        }
    }

    void setValues( const std::map< std::string, ValueType > &values,
                    const VectorString &groupsOfNodes = {} ) {
        setValues( values, this->getMesh()->getNodes( groupsOfNodes, true ) );
    }

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() const {
        _descriptor->updateValuePointer();
        _size->updateValuePointer();
        _component->updateValuePointer();
        _values->updateValuePointer();
        _allocated->updateValuePointer();
    };

    bool build() {
        updateValuePointers();

        _nbNodes = ( *_size )[0];
        _nbComp = ( *_size )[1];

        _buildComponentsName2Index();

        AS_ASSERT( _values->size() == _nbNodes * _nbComp );
        if ( !this->getMesh()->isParallel() ) {
            AS_ASSERT( _values->size() > 0 );
        }

        return true;
    }

    SimpleFieldOnNodesPtr restrict( const VectorString &cmps = {},
                                    const VectorString &groupsOfNodes = {},
                                    const ASTERINTEGER same_rank = PythonBool::None ) const {

        this->updateValuePointers();

        VectorString list_cmp;
        auto list_cmp_in = this->getComponents();
        if ( cmps.empty() ) {
            list_cmp = list_cmp_in;
        } else {
            auto set_cmps = toSet( cmps );

            for ( auto &cmp : list_cmp_in ) {
                if ( set_cmps.count( cmp ) > 0 ) {
                    list_cmp.push_back( cmp );
                }
            }
        }

        if ( list_cmp.empty() ) {
            raiseAsterError( "Restriction on list of components is empty" );
        }

        auto ret = std::make_shared< SimpleFieldOnNodes< ValueType > >(
            this->getMesh(), this->getPhysicalQuantity(), list_cmp );

        VectorLong nodes = _mesh->getNodes( groupsOfNodes, true, same_rank );

        ASTERINTEGER nbNodesGl = nodes.size();
#ifdef ASTER_HAVE_MPI
        if ( _mesh->isParallel() ) {
            nbNodesGl = AsterMPI::max( nbNodesGl );
        }
#endif

        if ( nbNodesGl == 0 ) {
            raiseAsterError( "Restriction on list of nodes is empty" );
        }

        for ( auto &cmp : list_cmp ) {
            auto icmp_in = ( *this )._name2Index.at( cmp );
            auto icmp = ( *ret )._name2Index.at( cmp );
            for ( auto &node : nodes ) {
                if ( this->hasValue( node, icmp_in ) ) {
                    ( *ret )( node, icmp ) = ( *this )( node, icmp_in );
                }
            }
        }

        return ret;
    };

    SimpleFieldOnNodesPtr asPhysicalQuantity( const std::string physQuant,
                                              const MapString &map_cmps ) const {

        if ( map_cmps.empty() ) {
            raiseAsterError( "Map to rename components is empty" );
        }

        VectorString cmps;
        cmps.reserve( map_cmps.size() );

        for ( auto const &[key, val] : map_cmps ) {
            cmps.push_back( key );
        }

        auto new_field = this->restrict( cmps, {} );

        if ( this->getPhysicalQuantity().back() != physQuant.back() ) {
            raiseAsterError( "Scalar type are differents" );
        }

        // Change names of components
        ( *new_field->_descriptor )[1] = physQuant;

        ASTERINTEGER nbCmp = new_field->getNumberOfComponents();
        auto list_cmps = new_field->getComponents();

        VectorString list_new_cmps;
        list_new_cmps.reserve( nbCmp );

        for ( auto const &cmp : list_cmps ) {
            list_new_cmps.push_back( map_cmps.at( cmp ) );
        }

        char *tabNames = vectorStringAsFStrArray( list_new_cmps, 8 );
        ASTERINTEGER iret;

        CALL_VERIGD_WRAP( physQuant.c_str(), tabNames, &nbCmp, &iret );

        FreeStr( tabNames );

        if ( iret != 0 ) {
            raiseAsterError( "Some components are not in the new PhysicalQuantity" );
        }

        for ( int i = 0; i < nbCmp; i++ ) {
            ( *new_field->_component )[i] = list_new_cmps[i];
        }

        new_field->build();
        return new_field;
    };

    /**
     * @brief Compute norm
     * @param normType Type of norm ("NORM_1","NORM_2","NORM_INFINITY")
     */
    ASTERDOUBLE norm( const std::string normType ) const {

        if constexpr ( !std::is_same_v< ValueType, ASTERDOUBLE > &&
                       !std::is_same_v< ValueType, ASTERINTEGER > &&
                       !std::is_same_v< ValueType, ASTERCOMPLEX > ) {
            raiseAsterError( " norm method not defined for type " +
                             std::string( typeid( ValueType ).name() ) );
        }
        CALL_JEMARQ();
        ASTERDOUBLE norme = 0.0;
        _values->updateValuePointer();

        const int rank = getMPIRank();

        JeveuxVectorLong nodesRank = getMesh()->getNodesOwner();
        nodesRank->updateValuePointer();
        const ASTERINTEGER nbNodes = _mesh->getNumberOfNodes();
        const ASTERINTEGER nbCmp = this->getNumberOfComponents();

        if ( normType == "NORM_1" ) {
            for ( ASTERINTEGER node = 0; node < nbNodes; ++node ) {
                if ( ( *nodesRank )[node] == rank ) {
                    for ( ASTERINTEGER icmp = 0; icmp < nbCmp; icmp++ ) {
                        if ( this->hasValue( node, icmp ) ) {
                            const auto val = ( *this )( node, icmp );
                            norme += std::abs( val );
                        }
                    }
                }
            }
        } else if ( normType == "NORM_2" ) {
            for ( ASTERINTEGER node = 0; node < nbNodes; ++node ) {
                if ( ( *nodesRank )[node] == rank ) {
                    for ( ASTERINTEGER icmp = 0; icmp < nbCmp; icmp++ ) {
                        if ( this->hasValue( node, icmp ) ) {
                            const auto val = ( *this )( node, icmp );
                            norme += val * val;
                        }
                    }
                }
            }
        } else if ( normType == "NORM_INFINITY" ) {
            for ( ASTERINTEGER node = 0; node < nbNodes; ++node ) {
                if ( ( *nodesRank )[node] == rank ) {
                    for ( ASTERINTEGER icmp = 0; icmp < nbCmp; icmp++ ) {
                        if ( this->hasValue( node, icmp ) ) {
                            const auto val = ( *this )( node, icmp );
                            norme = std::max( norme, std::abs( val ) );
                        }
                    }
                }
            }
        } else {
            raiseAsterError( "Unknown norm: " + normType );
        }

#ifdef ASTER_HAVE_MPI
        if ( this->getMesh()->isParallel() ) {
            ASTERDOUBLE norm2 = norme;
            if ( normType == "NORM_1" || normType == "NORM_2" )
                norme = AsterMPI::sum( norm2 );
            else
                norme = AsterMPI::max( norm2 );
        }
#endif

        if ( normType == "NORM_2" )
            norme = std::sqrt( norme );
        CALL_JEDEMA();
        return norme;
    };

    std::pair< std::vector< ValueType >, std::pair< VectorLong, VectorString > >
    getValuesWithDescription( const VectorString &cmps, const VectorString &groupsOfNodes ) const {

        VectorLong nodes = _mesh->getNodes( groupsOfNodes, true, PythonBool::None );

        return this->getValuesWithDescription( cmps, nodes );
    }

    std::pair< std::vector< ValueType >, std::pair< VectorLong, VectorString > >
    getValuesWithDescription( const VectorString &cmps, const VectorLong &nodes ) const {

        this->updateValuePointers();

        VectorLong nodes_ = nodes;
        if ( nodes_.empty() ) {
            nodes_ = _mesh->getNodes( VectorString(), true );
        }

        VectorString list_cmp;
        auto list_cmp_in = this->getComponents();
        if ( cmps.empty() ) {
            list_cmp = list_cmp_in;
        } else {
            auto set_cmps = toSet( cmps );

            for ( auto &cmp : list_cmp_in ) {
                if ( set_cmps.count( cmp ) > 0 ) {
                    list_cmp.push_back( cmp );
                }
            }
        }

        if ( list_cmp.empty() ) {
            raiseAsterError( "Restriction on list of components is empty" );
        }

        std::vector< ValueType > values;
        VectorLong v_nodes;
        VectorString v_cmps;

        values.reserve( nodes_.size() * list_cmp.size() );
        v_nodes.reserve( nodes_.size() * list_cmp.size() );
        v_cmps.reserve( nodes_.size() * list_cmp.size() );

        for ( auto &node : nodes_ ) {
            for ( auto &cmp : list_cmp ) {
                auto icmp = ( *this )._name2Index.at( cmp );
                if ( this->hasValue( node, icmp ) ) {
                    values.push_back( ( *this )( node, icmp ) );
                    v_nodes.push_back( node );
                    v_cmps.push_back( cmp );
                }
            }
        }

        return std::make_pair( values, std::make_pair( v_nodes, v_cmps ) );
    }
};

/** @typedef SimpleFieldOnNodesReal Class d'une champ simple de doubles */
typedef SimpleFieldOnNodes< ASTERDOUBLE > SimpleFieldOnNodesReal;

/**
 * @typedef SimpleFieldOnNodesPtrReal
 * @brief Definition d'un champ simple de doubles
 */
typedef std::shared_ptr< SimpleFieldOnNodesReal > SimpleFieldOnNodesRealPtr;

/** @typedef SimpleFieldOnNodesLong Class d'un champ simple de long */
typedef SimpleFieldOnNodes< long > SimpleFieldOnNodesLong;

/**
 * @typedef SimpleFieldOnNodesPtrLong
 * @brief Definition d'un champ simple de long
 */
typedef std::shared_ptr< SimpleFieldOnNodesLong > SimpleFieldOnNodesLongPtr;

/** @typedef SimpleFieldOnNodesComplex
    @brief Class d'un champ simple de complexes */
typedef SimpleFieldOnNodes< ASTERCOMPLEX > SimpleFieldOnNodesComplex;

/**
 * @typedef SimpleFieldOnNodesComplexPtr
 * @brief Definition d'un champ simple aux noeuds de complexes
 */
typedef std::shared_ptr< SimpleFieldOnNodesComplex > SimpleFieldOnNodesComplexPtr;
#endif /* SIMPLEFIELDONNODES_H_ */
