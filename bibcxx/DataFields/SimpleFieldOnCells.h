#ifndef SIMPLEFIELDONCELLS_H_
#define SIMPLEFIELDONCELLS_H_

/**
 * @file SimpleFieldOnCells.h
 * @brief Fichier entete de la classe SimpleFieldOnCells
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

#include "DataFields/DataField.h"
#include "MemoryManager/JeveuxVector.h"
#include "MemoryManager/NumpyAccess.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

/**
 * @class SimpleFieldOnCells
 * @brief Cette classe template permet de definir un champ aux éléments Aster
 * @author Nicolas Sellenet
 */
template < class ValueType >
class SimpleFieldOnCells : public DataField {
  private:
    /** @brief Vecteur Jeveux '.CESK' */
    JeveuxVectorChar8 _descriptor;
    /** @brief Vecteur Jeveux '.CESD' */
    JeveuxVectorLong _size;
    /** @brief Vecteur Jeveux '.CESC' */
    JeveuxVectorChar8 _component;
    /** @brief Vecteur Jeveux '.CESV' */
    JeveuxVector< ValueType > _values;
    /** @brief Vecteur Jeveux '.CESL' */
    JeveuxVectorLogical _allocated;
    /** @brief Mesh */
    BaseMeshPtr _mesh;
    /** @brief Nombre de éléments */
    ASTERINTEGER _nbCells;
    /** @brief Nombre de composantes */
    ASTERINTEGER _nbComp;
    /** @brief Number of points */
    ASTERINTEGER _nbPt;
    /** @brief Number of subpoints */
    ASTERINTEGER _nbSpt;

    std::map< std::string, ASTERINTEGER > _name2Index;

    /**
     * Some unsafe functions to access values without checking dimension
     * Their public version add an if statement
     */

    void _buildComponentsName2Index() {
        auto nbCmp = this->getNumberOfComponents();
        for ( ASTERINTEGER i = 0; i < nbCmp; i++ ) {
            _name2Index[this->getComponent( i )] = i;
        }
    }

    ASTERINTEGER _ptCell( const ASTERINTEGER &ima ) const { return ( *_size )[4 + 4 * ima + 1]; }
    ASTERINTEGER _sptCell( const ASTERINTEGER &ima ) const { return ( *_size )[4 + 4 * ima + 2]; }
    ASTERINTEGER _cmpsSptCell( const ASTERINTEGER &ima ) const {
        return ( *_size )[4 + 4 * ima + 3];
    }
    ASTERINTEGER _shiftCell( const ASTERINTEGER &ima ) const { return ( *_size )[4 + 4 * ima + 4]; }

    std::string _nameCmp( const ASTERINTEGER &icmp ) const {
        return strip( ( *_component )[icmp].toString() );
    }

    /**
     * Calculate the position of value in CESV array
     */
    ASTERINTEGER _positionInArray( const ASTERINTEGER &ima, const ASTERINTEGER &icmp,
                                   const ASTERINTEGER &ipt, const ASTERINTEGER &ispt ) const {

#ifdef ASTER_DEBUG_CXX_LOW_LEVEL
        this->_checkAllOOR( ima, icmp, ipt, ispt );
#endif

        ASTERINTEGER npt = this->_ptCell( ima );
        ASTERINTEGER nspt = this->_sptCell( ima );
        ASTERINTEGER ncmp = this->_cmpsSptCell( ima );
        ASTERINTEGER decal = this->_shiftCell( ima );
        return decal + ipt * nspt * ncmp + ispt * ncmp + icmp;
    };

    /**
     * Calculate the size of CESV array
     */
    ASTERINTEGER _nbValArray() const {
        ASTERINTEGER nbVal = 0;
        ASTERINTEGER ncmp_max = this->getNumberOfComponents();
        for ( ASTERINTEGER ima = 0; ima < this->getNumberOfCells(); ima++ ) {
            ASTERINTEGER npt = this->_ptCell( ima );
            ASTERINTEGER nspt = this->_sptCell( ima );
            ASTERINTEGER ncmp = this->_cmpsSptCell( ima );
            if ( ncmp > 0 )
                AS_ASSERT( ncmp == ncmp_max );
            nbVal = nbVal + npt * nspt * ncmp;
        }
        AS_ASSERT( nbVal > 0 );
        return nbVal;
    }

    /**
     * Functions to check an out-of-range condition
     */
    void _checkCellOOR( const ASTERINTEGER &ima ) const {
        ASTERINTEGER nbCells = this->getNumberOfCells();
        if ( ima < 0 || ima >= nbCells ) {
            throw std::runtime_error( "Cell index '" + std::to_string( ima ) +
                                      "' is out of range  ( " + std::to_string( nbCells ) + "  )" );
        };
    }

    void _checkPtOOR( const ASTERINTEGER &ima, const ASTERINTEGER &ipt ) const {
        ASTERINTEGER npt = this->_ptCell( ima );
        if ( ipt < 0 || ipt >= npt ) {
            throw std::runtime_error( "Point     '" + std::to_string( ipt ) +
                                      "' is out of range  ( " + std::to_string( npt ) +
                                      "  ) for cell '" + std::to_string( ima ) + "'" );
        }
    }

    void _checkSptOOR( const ASTERINTEGER &ima, const ASTERINTEGER &ispt ) const {
        ASTERINTEGER nspt = this->_sptCell( ima );
        if ( ispt < 0 || ispt >= nspt ) {
            throw std::runtime_error( "SubPoint  '" + std::to_string( ispt ) +
                                      "' is out of range ( " + std::to_string( nspt ) +
                                      "  ) for cell '" + std::to_string( ima ) + "'" );
        }
    }

    void _checkCmpAtCellOOR( const ASTERINTEGER &ima, const ASTERINTEGER &icmp ) const {
        ASTERINTEGER ncmp = this->_cmpsSptCell( ima );
        if ( icmp < 0 || icmp >= ncmp ) {
            throw std::runtime_error( "Component '" + std::to_string( icmp ) +
                                      "' is out of range  ( " + std::to_string( ncmp ) +
                                      "  ) for cell '" + std::to_string( ima ) + "'" );
        }
    }

    void _checkCmpOOR( const ASTERINTEGER &icmp ) const {
        ASTERINTEGER ncmp = this->getNumberOfComponents();
        if ( icmp < 0 || icmp >= ncmp ) {
            throw std::runtime_error( "Component '" + std::to_string( icmp ) +
                                      "' is out of range  ( " + std::to_string( ncmp ) + "  )" );
        }
    }

    void _checkAllOOR( const ASTERINTEGER &ima, const ASTERINTEGER &icmp, const ASTERINTEGER &ipt,
                       const ASTERINTEGER &ispt ) const {
        this->_checkCellOOR( ima );
        this->_checkPtOOR( ima, ipt );
        this->_checkSptOOR( ima, ispt );
        this->_checkCmpAtCellOOR( ima, icmp );
    };

    bool _hasValueCells( const ASTERINTEGER &ima ) const {

        ASTERINTEGER npt = getNumberOfPointsOfCell( ima );
        ASTERINTEGER nspt = getNumberOfSubPointsOfCell( ima );
        ASTERINTEGER ncmp = getNumberOfComponents();

        for ( ASTERINTEGER icmp = 0; icmp < ncmp; icmp++ ) {
            for ( ASTERINTEGER ipt = 0; ipt < npt; ipt++ ) {
                for ( ASTERINTEGER ispt = 0; ispt < nspt; ispt++ ) {
                    if ( hasValue( ima, icmp, ipt, ispt ) ) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

  public:
    /**
     * @typedef SimpleFieldOnCellsPtr
     * @brief Pointeur intelligent vers un SimpleFieldOnCells
     */
    typedef std::shared_ptr< SimpleFieldOnCells > SimpleFieldOnCellsPtr;

    /**
     * @brief Constructeur

     */
    SimpleFieldOnCells() = delete;

    /**
     * @brief Constructeur
     * @param name Nom Jeveux du champ aux éléments
     */
    SimpleFieldOnCells( const std::string name, const BaseMeshPtr mesh = nullptr )
        : DataField( name, "CHAM_ELEM_S" ),
          _descriptor( JeveuxVectorChar8( getName() + ".CESK" ) ),
          _size( JeveuxVectorLong( getName() + ".CESD" ) ),
          _component( JeveuxVectorChar8( getName() + ".CESC" ) ),
          _values( JeveuxVector< ValueType >( getName() + ".CESV" ) ),
          _allocated( JeveuxVectorLogical( getName() + ".CESL" ) ),
          _nbCells( 0 ),
          _nbComp( 0 ),
          _nbPt( 0 ),
          _nbSpt( 0 ),
          _mesh( mesh ) {};

    SimpleFieldOnCells( const BaseMeshPtr mesh )
        : SimpleFieldOnCells( DataStructureNaming::getNewName( 19 ), mesh ) {};

    SimpleFieldOnCells( const BaseMeshPtr mesh, const std::string &loc, const std::string &quantity,
                        const VectorString &comp, bool zero = false )
        : SimpleFieldOnCells( mesh, loc, quantity, comp, 1, 1, zero ) {

        AS_ASSERT( loc == "ELEM" );
    }

    SimpleFieldOnCells( const BaseMeshPtr mesh, const std::string &loc, const std::string &quantity,
                        const VectorString &comp, const ASTERINTEGER &nbPG,
                        const ASTERINTEGER &nbSP, bool zero = false )
        : SimpleFieldOnCells( mesh ) {

        this->allocate( loc, quantity, comp, nbPG, nbSP, zero );
    }

    SimpleFieldOnCells( const BaseMeshPtr mesh, const std::string &loc, const std::string &quantity,
                        const VectorString &comp, const VectorInt &nbPG, const ASTERINTEGER &nbSP,
                        bool zero = false )
        : SimpleFieldOnCells( mesh ) {

        VectorLong nbPgTmp;
        nbPgTmp.reserve( nbPG.size() );
        for ( const auto &val : nbPG ) {
            nbPgTmp.push_back( val );
        }

        const std::string base = "G";
        ASTERINTEGER nbComp = comp.size();
        ASTERINTEGER _nsp = nbSP;
        ASTERLOGICAL _zero = zero;

        char *tabNames = vectorStringAsFStrArray( comp, 8 );

        CALL_CESCRE_WRAP( base.c_str(), getName().c_str(), loc.c_str(), _mesh->getName().c_str(),
                          quantity.c_str(), &nbComp, tabNames, nbPgTmp.data(), &_nsp, &nbComp,
                          &_zero );

        FreeStr( tabNames );

        build();
    }

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    SimpleFieldOnCells( const py::tuple &tup )
        : SimpleFieldOnCells( tup[0].cast< std::string >(), tup[1].cast< BaseMeshPtr >() ) {
        build();
    };

    py::tuple _getState() const { return py::make_tuple( getName(), getMesh() ); };

    BaseMeshPtr getMesh() const { return _mesh; };

    void allocate( const std::string &loc, const std::string &quantity, const VectorString &comp,
                   const ASTERINTEGER &nbPG, ASTERINTEGER nbSP = 1, bool zero = false ) {
        const std::string base = "G";
        ASTERINTEGER nbComp = comp.size();
        ASTERINTEGER _npg = -nbPG;
        ASTERINTEGER _nsp = nbSP;
        ASTERLOGICAL _zero = zero;

        char *tabNames = vectorStringAsFStrArray( comp, 8 );

        CALL_CESCRE_WRAP( base.c_str(), getName().c_str(), loc.c_str(), _mesh->getName().c_str(),
                          quantity.c_str(), &nbComp, tabNames, &_npg, &_nsp, &nbComp, &_zero );

        FreeStr( tabNames );

        build();
    }

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ValueType &operator[]( const ASTERINTEGER &i ) { return _values->operator[]( i ); };

    inline const ValueType &operator[]( const ASTERINTEGER &i ) const {
        return _values->operator[]( i );
    };

    ValueType &operator()( const ASTERINTEGER &ima, const ASTERINTEGER &icmp,
                           const ASTERINTEGER &ipt, const ASTERINTEGER &ispt ) {
#ifdef ASTER_DEBUG_CXX_LOW_LEVEL
        if ( this->getNumberOfCells() == 0 || this->getNumberOfComponents() == 0 )
            throw std::runtime_error( "First call of updateValuePointers is mandatory" );
#endif

        ASTERINTEGER position = this->_positionInArray( ima, icmp, ipt, ispt );

        ( *_allocated )[position] = true;
        return this->operator[]( position );
    };

    ValueType &operator()( const ASTERINTEGER &ima, const std::string &cmp, const ASTERINTEGER &ipt,
                           const ASTERINTEGER &ispt ) {
        auto icmp = _name2Index.at( cmp );

        return this->operator()( ima, icmp, ipt, ispt );
    }

    const ValueType &operator()( const ASTERINTEGER &ima, const ASTERINTEGER &icmp,
                                 const ASTERINTEGER &ipt, const ASTERINTEGER &ispt ) const {
#ifdef ASTER_DEBUG_CXX_LOW_LEVEL
        if ( this->getNumberOfCells() == 0 || this->getNumberOfComponents() == 0 )
            throw std::runtime_error( "First call of updateValuePointers is mandatory" );
#endif

        ASTERINTEGER position = this->_positionInArray( ima, icmp, ipt, ispt );

        if ( !( *_allocated )[position] ) {
            std::string mess = "DEBUG: Position (" + std::to_string( icmp ) + ", " +
                               std::to_string( ima ) + ", " + std::to_string( ipt ) + ", " +
                               std::to_string( ispt ) + ") is valid but not allocated!";
            raiseAsterError( mess );
        };

        return this->operator[]( position );
    };

    const ValueType &operator()( const ASTERINTEGER &ima, const std::string &cmp,
                                 const ASTERINTEGER &ipt, const ASTERINTEGER &ispt ) const {
        auto icmp = _name2Index.at( cmp );

        return this->operator()( ima, icmp, ipt, ispt );
    }

    ValueType &operator()( const ASTERINTEGER &ima, const ASTERINTEGER &icmp,
                           const ASTERINTEGER &ipt ) {

        return ( *this )( ima, icmp, ipt, 0 );
    };

    const ValueType &operator()( const ASTERINTEGER &ima, const ASTERINTEGER &icmp,
                                 const ASTERINTEGER &ipt ) const {
        return ( *this )( ima, icmp, ipt, 0 );
    };

    /**
     * Calculate the position of value in CESV array
     */
    ASTERINTEGER getPositionInArray( const ASTERINTEGER &ima, const ASTERINTEGER &icmp,
                                     const ASTERINTEGER &ipt, const ASTERINTEGER &ispt ) {
        return _positionInArray( ima, icmp, ipt, ispt );
    }

    /**
     * @brief Get the values of the field
     *
     */
    const JeveuxVector< ValueType > &getValues() const { return _values; }

    /**
     * @brief Return logical value vector
     */
    JeveuxVectorLogical getLogicalValues() const { return _allocated; }

    /**
     * @brief Return vector of sizes
     */
    JeveuxVectorLong getSizes() const { return _size; }

    /**
     * @brief Access to the (icmp) component of the (ima) cell
              at the (ipt) point, at the (ispt) sub-point.
    */
    ValueType const &getValue( const ASTERINTEGER &ima, const ASTERINTEGER &icmp,
                               const ASTERINTEGER &ipt, const ASTERINTEGER ispt = 0 ) const {

        return ( *this )( ima, icmp, ipt, ispt );
    }

    void setValue( const ASTERINTEGER &ima, const ASTERINTEGER &icmp, const ASTERINTEGER &ipt,
                   const ASTERINTEGER &ispt, const ValueType &val ) {

        ( *this )( ima, icmp, ipt, ispt ) = val;
    }

    void setValue( const ASTERINTEGER &ima, const ASTERINTEGER &icmp, const ASTERINTEGER &ipt,
                   const ValueType &val ) {

        this->setValue( ima, icmp, ipt, 0, val );
    }

    /**
     * @brief tell if value exists for (icmp) component of the (ima) cell
              at the (ipt) point, at the (ispt) sub-point.
    */
    bool hasValue( const ASTERINTEGER &ima, const ASTERINTEGER &icmp, const ASTERINTEGER &ipt,
                   const ASTERINTEGER &ispt = 0 ) const {

#ifdef ASTER_DEBUG_CXX_LOW_LEVEL
        if ( this->getNumberOfCells() == 0 || this->getNumberOfComponents() == 0 )
            throw std::runtime_error( "First call of updateValuePointers is mandatory" );
#endif

        const auto nbCells = this->getNumberOfCells();
        if ( ima < 0 || nbCells == 0 || ima >= nbCells ) {
            return false;
        }

        const auto ncmp = this->_cmpsSptCell( ima );
        if ( icmp < 0 || ncmp == 0 || icmp >= ncmp ) {
            return false;
        }

        const auto npt = this->_ptCell( ima );
        if ( ipt < 0 || npt == 0 || ipt >= npt ) {
            return false;
        }

        const auto nspt = this->_sptCell( ima );
        if ( ispt < 0 || nspt == 0 || ispt >= nspt ) {
            return false;
        }

        ASTERINTEGER position = this->_positionInArray( ima, icmp, ipt, ispt );

        return ( *_allocated )[position];
    }

    bool hasValue( const ASTERINTEGER &ima, const std::string &cmp, const ASTERINTEGER &ipt,
                   const ASTERINTEGER &ispt = 0 ) const {
        auto icmp = _name2Index.at( cmp );

        return this->hasValue( ima, icmp, ipt, ispt );
    }

    /**
     * @brief Get number of points of the i-th cell
     */
    ASTERINTEGER getNumberOfPointsOfCell( const ASTERINTEGER &ima ) const {
        this->_checkCellOOR( ima );
        return this->_ptCell( ima );
    }

    /**
     * @brief Get number of sub-points of the i-th cell
     */
    ASTERINTEGER getNumberOfSubPointsOfCell( const ASTERINTEGER &ima ) const {
        this->_checkCellOOR( ima );
        return this->_sptCell( ima );
    }

    /**
     * @brief Get number of components of the i-th cell
     */
    ASTERINTEGER
    getNumberOfComponentsForSubpointsOfCell( const ASTERINTEGER &ima ) const {
        this->_checkCellOOR( ima );
        return this->_cmpsSptCell( ima );
    }

    /**
     * @brief Get number of components
     */
    ASTERINTEGER getNumberOfComponents() const { return _nbComp; }

    /**
     * @brief Get number of cells
     */
    ASTERINTEGER getNumberOfCells() const { return _nbCells; }

    /**
     * @brief Get number of points
     */
    ASTERINTEGER getMaxNumberOfPoints() const { return _nbPt; }

    /**
     * @brief Get number of sub-points
     */
    ASTERINTEGER getMaxNumberOfSubPoints() const { return _nbSpt; }

    /**
     * @brief Get the name of the i-th component
     */
    std::string getComponent( const ASTERINTEGER &icmp ) const {

        if ( icmp < 0 || icmp >= this->getNumberOfComponents() ) {
            throw std::runtime_error( "Component '" + std::to_string( icmp ) +
                                      "' is out of range" );
        };
        return this->_nameCmp( icmp );
    };

    /**
     * @Brief Get the names of all the components
     */
    VectorString getComponents() const {

        ASTERINTEGER size = this->getNumberOfComponents();
        VectorString names;
        names.reserve( size );
        for ( ASTERINTEGER icmp = 0; icmp < size; icmp++ ) {
            names.push_back( this->_nameCmp( icmp ) );
        }
        return names;
    }

    /**
     * @brief Get physical quantity
     */
    std::string getPhysicalQuantity() const { return strip( ( *_descriptor )[1].toString() ); }

    /**
     * @brief Get field location
     */
    std::string getLocalization() const { return strip( ( *_descriptor )[2].toString() ); }

    /**
     * @brief Get cells holding components
     */
    VectorLong getCellsWithValues() const {
        auto nbCells = this->getNumberOfCells();

        VectorLong values;
        values.reserve( nbCells );

        for ( ASTERINTEGER ima = 0; ima < nbCells; ima++ ) {
            if ( this->_hasValueCells( ima ) ) {
                values.push_back( ima );
            }
        }
        return values;
    }

    std::pair< std::vector< ValueType >,
               std::tuple< VectorLong, VectorString, VectorLong, VectorLong > >
    getValuesWithDescription( const VectorString &cmps, const VectorString &groupsOfCells ) const {

        VectorLong cells = _mesh->getCells( groupsOfCells );

        return this->getValuesWithDescription( cmps, cells );
    }

    std::pair< std::vector< ValueType >,
               std::tuple< VectorLong, VectorString, VectorLong, VectorLong > >
    getValuesWithDescription( const VectorString &cmps, const VectorLong &cells ) const {

        std::vector< ValueType > values;
        VectorLong v_cells;
        VectorString v_cmps;
        VectorLong points;
        VectorLong subpoints;

        VectorLong cells_ = cells;
        if ( cells_.empty() ) {
            cells_ = _mesh->getCells();
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

        this->updateValuePointers();

        ASTERINTEGER size =
            cells_.size() * cmps.size() * getMaxNumberOfPoints() * getMaxNumberOfSubPoints();
        v_cells.reserve( size );
        v_cmps.reserve( size );
        values.reserve( size );
        points.reserve( size );
        subpoints.reserve( size );

        for ( auto &cell : cells_ ) {
            ASTERINTEGER npt = getNumberOfPointsOfCell( cell );
            ASTERINTEGER nspt = getNumberOfSubPointsOfCell( cell );
            for ( auto &cmp : list_cmp ) {
                auto icmp = _name2Index.at( cmp );

                for ( ASTERINTEGER ipt = 0; ipt < npt; ipt++ ) {
                    for ( ASTERINTEGER ispt = 0; ispt < nspt; ispt++ ) {
                        if ( hasValue( cell, icmp, ipt, ispt ) ) {
                            v_cells.push_back( cell );
                            v_cmps.push_back( cmp );
                            values.push_back( getValue( cell, icmp, ipt, ispt ) );
                            points.push_back( ipt );
                            subpoints.push_back( ispt );
                        }
                    }
                }
            }
        }

        return make_pair( values, make_tuple( v_cells, v_cmps, points, subpoints ) );
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
    }

    bool build() {
        updateValuePointers();

        _nbCells = ( *_size )[0];
        _nbComp = ( *_size )[1];
        _nbPt = ( *_size )[2];
        _nbSpt = ( *_size )[3];

        AS_ASSERT( _values->size() == this->_nbValArray() );

        _buildComponentsName2Index();

        return true;
    };

    void setValues( const VectorLong &cells, const VectorString &cmps, const VectorLong &npg,
                    const VectorLong &spt, const std::vector< ValueType > &values ) {

        AS_ASSERT( cells.size() == cmps.size() );
        AS_ASSERT( cells.size() == values.size() );
        AS_ASSERT( cells.size() == npg.size() );
        AS_ASSERT( cells.size() == spt.size() );

        this->updateValuePointers();

        const int size = values.size();
        for ( int i = 0; i < size; i++ ) {
            ( *this )( cells[i], cmps[i], npg[i], spt[i] ) = values[i];
        }
    }

    void setValues( const std::vector< ValueType > &values ) {
        if ( values.size() != _values->size() ) {
            raiseAsterError( "Incompatible size in setValues, expecting " +
                             std::to_string( _values->size() ) + ", got " +
                             std::to_string( values.size() ) );
        }

        _allocated->assign( true );
        *_values = values;
    }

    SimpleFieldOnCellsPtr restrict( const VectorString &cmps = {},
                                    const VectorString &groupsOfCells = {} ) const {

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

        auto ret = std::make_shared< SimpleFieldOnCells< ValueType > >( getMesh() );

        VectorLong cells = _mesh->getCells( groupsOfCells );
        for ( auto &cell : cells ) {
            cell += 1;
        }

        ASTERINTEGER nbCellsGl = cells.size();
#ifdef ASTER_HAVE_MPI
        if ( _mesh->isParallel() ) {
            nbCellsGl = AsterMPI::max( nbCellsGl );
        }
#endif

        if ( nbCellsGl == 0 ) {
            raiseAsterError( "Restriction on list of cells is empty" );
        }

        char *tabNames = vectorStringAsFStrArray( list_cmp, 8 );
        ASTERINTEGER nbCells = cells.size(), nbCmp = list_cmp.size();
        const std::string base = "G";

        CALL_CESRED_WRAP( getName().c_str(), &nbCells, cells.data(), &nbCmp, tabNames, base.c_str(),
                          ret->getName().c_str() );

        FreeStr( tabNames );

        ret->build();
        return ret;
    };

    SimpleFieldOnCellsPtr asPhysicalQuantity( const std::string physQuant,
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

    auto getComponentsName2Index() const { return _name2Index; }

    ASTERDOUBLE norm( const std::string normType ) const {
        AS_ASSERT( normType == "NORM_1" || normType == "NORM_2" || normType == "NORM_INFINITY" );

        ASTERDOUBLE norme = 0.0;

        const int rank = getMPIRank();

        JeveuxVectorLong cellsRank = getMesh()->getCellsOwner();
        cellsRank->updateValuePointer();

        const ASTERINTEGER nbCells = getMesh()->getNumberOfCells();
        const ASTERINTEGER nbCmp = this->getNumberOfComponents();

        for ( ASTERINTEGER cell = 0; cell < nbCells; ++cell ) {
            if ( ( *cellsRank )[cell] == rank ) {
                ASTERINTEGER npt = getNumberOfPointsOfCell( cell );
                ASTERINTEGER nspt = getNumberOfSubPointsOfCell( cell );
                for ( ASTERINTEGER icmp = 0; icmp < nbCmp; icmp++ ) {
                    for ( ASTERINTEGER ipt = 0; ipt < npt; ipt++ ) {
                        for ( ASTERINTEGER ispt = 0; ispt < nspt; ispt++ ) {
                            if ( hasValue( cell, icmp, ipt, ispt ) ) {
                                const auto val = this->getValue( cell, icmp, ipt, ispt );

                                if ( normType == "NORM_1" ) {
                                    norme += std::abs( val );
                                } else if ( normType == "NORM_2" ) {
                                    norme += val * val;
                                } else if ( normType == "NORM_INFINITY" ) {
                                    norme = std::max( norme, std::abs( val ) );
                                }
                            }
                        }
                    }
                }
            }
        }

#ifdef ASTER_HAVE_MPI
        if ( getMesh()->isParallel() ) {
            ASTERDOUBLE norm2 = norme;
            if ( normType == "NORM_1" || normType == "NORM_2" )
                AsterMPI::all_reduce( norm2, norme, MPI_SUM );
            else
                AsterMPI::all_reduce( norm2, norme, MPI_MAX );
        }
#endif

        if ( normType == "NORM_2" )
            norme = std::sqrt( norme );

        return norme;
    }
};

using SimpleFieldOnCellsReal = SimpleFieldOnCells< ASTERDOUBLE >;
using SimpleFieldOnCellsRealPtr = std::shared_ptr< SimpleFieldOnCellsReal >;
using SimpleFieldOnCellsLong = SimpleFieldOnCells< ASTERINTEGER >;
using SimpleFieldOnCellsLongPtr = std::shared_ptr< SimpleFieldOnCellsLong >;

#endif /* SIMPLEFIELDONCELLS_H_ */
