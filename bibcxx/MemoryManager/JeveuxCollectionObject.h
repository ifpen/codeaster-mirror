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

#include "astercxx.h"

#include "aster_fort_jeveux.h"
#include "aster_utils.h"

#include "MemoryManager/JeveuxAllowedTypes.h"
#include "MemoryManager/JeveuxObject.h"
#include "MemoryManager/JeveuxString.h"
#include "MemoryManager/NamesMap.h"
#include "Utilities/Blas.h"
#include "Utilities/Tools.h"

#pragma once

/**
 * @class JeveuxCollectionObject
 * @brief Cette classe template permet de definir un objet de collection Jeveux
 * @author Nicolas Sellenet
 */
template < class ValueType >
class JeveuxCollectionObjectClass : private AllowedJeveuxType< ValueType > {
  private:
    /** @brief Nom Jeveux de la collection */
    const JeveuxObjectClass *_collection;
    /** @brief Position dans la collection */
    ASTERINTEGER _numberInCollection;
    /** @brief Pointeur vers le vecteur Jeveux */
    ValueType *_valuePtr;
    /** @brief Longueur du vecteur Jeveux */
    ASTERINTEGER _size;
    /** @brief Adresse du vecteur Jeveux */
    ASTERINTEGER _jeveuxAdress;

    /**
     * @brief Allocation
     */
    std::string getJeveuxName() const {
        std::string charJeveuxName( 32, ' ' );
        ASTERINTEGER num = _numberInCollection;
        CALLO_JEXNUM( charJeveuxName, _collection->getName(), &num );

        return charJeveuxName;
    }

    /**
     * @brief Allocation
     */
    void allocate( const ASTERINTEGER &size, const std::string name = "" ) {
        ASTERINTEGER taille = size, ibid = 0;

        std::string nameOfObject( "" );
        if ( name != "" )
            nameOfObject = name;
        else
            ibid = _numberInCollection;

        CALLO_JUCROC_WRAP( _collection->getName(), nameOfObject, &ibid, &taille,
                           (void *)( &_valuePtr ) );

        ASTERINTEGER valTmp;
        JeveuxChar8 param( "IADM" );
        std::string charval = std::string( 32, ' ' );
        auto charJeveuxName = getJeveuxName();
        CALLO_JELIRA( charJeveuxName, param, &valTmp, charval );
        _jeveuxAdress = valTmp;
    };

  public:
    /**
     * @brief Constructeur
     * @param collectionName Nom de collection
     * @param number Numero de l'objet dans la collection
     * @param objectName Nom de l'objet de collection
     */
    JeveuxCollectionObjectClass( const JeveuxObjectClass &coll, const ASTERINTEGER &number,
                                 bool isNamed )
        : _collection( &coll ),
          _numberInCollection( number ),
          _valuePtr( nullptr ),
          _size( 0 ),
          _jeveuxAdress( 0 ) {

        if ( !exists() ) {
            AS_ABORT( "Error in collection object " + _collection->getName() );
        }

        std::string charJeveuxName = getJeveuxName();
        ASTERINTEGER valTmp;
        JeveuxChar8 param( "LONMAX" );
        std::string charval = std::string( 32, ' ' );
        CALLO_JELIRA( charJeveuxName, param, &valTmp, charval );
        _size = valTmp;
    };

    /**
     * @brief Constructeur
     * @param collectionName Nom de collection
     * @param number Numero de l'objet dans la collection
     */
    JeveuxCollectionObjectClass( const JeveuxObjectClass &coll, const ASTERINTEGER &number,
                                 const ASTERINTEGER &size )
        : _collection( &coll ),
          _numberInCollection( number ),
          _valuePtr( nullptr ),
          _size( size ),
          _jeveuxAdress( 0 ) {
        allocate( size );
    };

    /**
     * @brief Constructeur
     * @param collectionName Nom de collection
     * @param number Numero de l'objet dans la collection
     * @param objectName Nom de l'objet de collection
     */
    JeveuxCollectionObjectClass( const JeveuxObjectClass &coll, const ASTERINTEGER &number,
                                 const std::string &objectName, const ASTERINTEGER &size )
        : _collection( &coll ),
          _numberInCollection( number ),
          _valuePtr( nullptr ),
          _size( size ),
          _jeveuxAdress( 0 ) {
        allocate( size, objectName );
    };

    struct const_iterator {
        ASTERINTEGER position;
        const ValueType &valuePtr;

        inline const_iterator( ASTERINTEGER memoryPosition, const ValueType &val )
            : position( memoryPosition ), valuePtr( val ) {};

        inline const_iterator( const const_iterator &iter )
            : position( iter.position ), valuePtr( iter.valuePtr ) {};

        inline const_iterator &operator=( const const_iterator &testIter ) {
            position = testIter.position;
            valuePtr = testIter.valuePtr;
            return *this;
        };

        inline const_iterator &operator++() {
            ++position;
            return *this;
        };

        inline bool operator==( const const_iterator &testIter ) const {
            if ( testIter.position != position )
                return false;
            return true;
        };

        inline bool operator!=( const const_iterator &testIter ) const {
            if ( testIter.position != position )
                return true;
            return false;
        };

        inline const ValueType &operator->() const { return ( &valuePtr )[position]; };

        inline const ValueType &operator*() const { return ( &valuePtr )[position]; };
    };

    /**
     * @brief
     */
    const_iterator begin() const { return const_iterator( 0, *_valuePtr ); };

    /**
     * @brief
     * @todo revoir le fonctionnement du end car il peut provoquer de segfault
     */
    const_iterator end() const { return const_iterator( _size, *_valuePtr ); };

    /**
     * @brief test de l'existance
     * @return true si l'objet existe
     */
    bool exists() {
        std::string charval = getJeveuxName();
        ASTERINTEGER iret = 0;
        CALLO_JEEXIN( charval, &iret );

        if ( iret == 0 )
            return false;

        return true;
    };

    /**
     * @brief Get index of object in the collection
     * @return index
     */
    ASTERINTEGER getIndex() const { return _numberInCollection; };

    /**
     * @brief Mise a jour du pointeur Jeveux
     * @return true si la mise a jour s'est bien passee
     */
    void updateValuePointer() const {
        if ( hasMoved() ) {
            const_cast< JeveuxCollectionObjectClass< ValueType > * >( this )->_valuePtr = NULL;

            const std::string read( "L" );
            std::string charJeveuxName = getJeveuxName();
            CALLO_JEVEUOC( charJeveuxName, read, (void *)( &_valuePtr ) );
            AS_ASSERT( _valuePtr );

            ASTERINTEGER valTmp = 0;
            JeveuxChar8 param( "IADM" );
            std::string charval = std::string( 32, ' ' );
            CALLO_JELIRA( charJeveuxName, param, &valTmp, charval );
            const_cast< JeveuxCollectionObjectClass< ValueType > * >( this )->_jeveuxAdress =
                valTmp;
        }
    };

    inline const ValueType &operator[]( const ASTERINTEGER &i ) const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
        if ( i < 0 || i >= this->size() ) {
            std::string error = "Out of range of JeveuxCollectionObjectClass '" + getStringName() +
                                "', index = " + std::to_string( i ) +
                                " ( size = " + std::to_string( this->size() ) + " )";
            AS_ABORT( error );
        }
#endif

        return _valuePtr[i];
    };

    inline ValueType &operator[]( const ASTERINTEGER &i ) {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
        if ( i < 0 || i >= this->size() ) {
            std::string error = "Out of range of JeveuxCollectionObjectClass '" + getStringName() +
                                "', index = " + std::to_string( i ) +
                                " ( size = " + std::to_string( this->size() ) + " )";
            AS_ABORT( error );
        }
#endif

        return _valuePtr[i];
    };

    JeveuxChar32 getName() const { return JeveuxChar32( getStringName() ); };

    std::string getStringName() const {
        std::string collectionObjectName( 32, ' ' );
        std::string charJeveuxName = getJeveuxName();
        CALLO_JENUNO( charJeveuxName, collectionObjectName );
        return strip( std::string( collectionObjectName ) );
    };

    bool hasMoved() const {
        std::string charJeveuxName = getJeveuxName();

        // check if object is in memory
        ASTERINTEGER idummy;
        JeveuxChar8 param0( "USAGE" );
        JeveuxChar32 usage( " " );
        CALLO_JELIRA( charJeveuxName, param0, &idummy, usage );

        if ( usage[0] == 'X' ) {
            const_cast< JeveuxCollectionObjectClass< ValueType > * >( this )->_valuePtr = NULL;
            return true;
        }

        ASTERINTEGER valTmp;
        JeveuxChar8 param( "IADM" );
        JeveuxChar32 dummy( " " );
        CALLO_JELIRA( charJeveuxName, param, &valTmp, dummy );
        if ( valTmp != _jeveuxAdress ) {
            const_cast< JeveuxCollectionObjectClass< ValueType > * >( this )->_valuePtr = NULL;
            return true;
        }

        return false;
    };

    /** @brief Set values of collection object */
    void setValues( const std::vector< ValueType > &toCopy ) {
        if ( toCopy.size() != size() ) {
            AS_ABORT( "Sizes do not match: " + std::to_string( size() ) + " vs " +
                      std::to_string( toCopy.size() ) );
        }
        this->updateValuePointer();
        ASTERINTEGER pos = 0;
        for ( const auto &val : toCopy ) {
            _valuePtr[pos] = val;
            ++pos;
        }
    };

    /** @brief Set values of collection object */
    void setValues( const ValueType &toCopy ) {
        if ( size() != 1 ) {
            AS_ABORT( "Sizes do not match: " + std::to_string( size() ) + " vs " +
                      std::to_string( 1 ) );
        }
        _valuePtr[0] = toCopy;
    };

    /** @brief Set values of collection object */
    void setValues( const JeveuxCollectionObjectClass< ValueType > &toCopy ) {
        if ( size() != toCopy.size() ) {
            AS_ABORT( "Sizes do not match: " + std::to_string( size() ) + " vs " +
                      std::to_string( toCopy.size() ) );
        }
        for ( int i = 0; i < size(); ++i ) {
            _valuePtr[i] = toCopy[i];
        }
    };

    /** @brief Get size of collection object */
    ASTERINTEGER size() const { return _size; };

    /** @brief Get capacity of collection object */
    ASTERINTEGER capacity() const { return _size; };

    /** @brief Convert to std::vector */
    std::vector< ValueType > toVector() const {
        this->updateValuePointer();
        std::vector< ValueType > toReturn;
        toReturn.reserve( size() );
        for ( ASTERINTEGER i = 0; i < size(); ++i )
            toReturn.push_back( _valuePtr[i] );
        return toReturn;
    };

    bool operator==( const JeveuxCollectionObjectClass< ValueType > &toCompar ) const {
        if ( this->size() != toCompar.size() )
            return false;

        auto size = this->size();
        for ( ASTERINTEGER i = 0; i < size; i++ ) {
            if ( this->operator[]( i ) != toCompar[i] )
                return false;
        }

        return true;
    };

    /** @brief checks whether the container is empty */
    bool empty() const {
        if ( this->size() == 0 )
            return true;

        return false;
    }

    /**
     * @brief Return a pointer to the vector
     */
    ValueType *getDataPtr() { return _valuePtr; };

    /**
     * @brief Return a pointer to the vector
     */
    const ValueType *getDataPtr() const { return _valuePtr; };

    /**
     * @brief TimesEqual overloading
     * @return Updated JeveuxCollectionObject
     */
    JeveuxCollectionObjectClass< ValueType > &operator*=( const ValueType &scal ) {
        this->updateValuePointer();
        const auto size = this->size();

        AsterBLAS::scal( size, scal, getDataPtr(), ASTERINTEGER( 1 ) );

        return *this;
    };

    /** @brief overload << operator */
    friend std::ostream &operator<<( std::ostream &os,
                                     const JeveuxCollectionObjectClass< ValueType > &toPrint ) {
        os << "JeveuxCollectionObject: " << toPrint.getStringName() << "\n";
        os << "Size: " << std::to_string( toPrint.size() )
           << ", and capacity: " << std::to_string( toPrint.capacity() ) << ".\n";

        const auto size = toPrint.size();
        os << "List of values: \n";
        os << "( ";
        for ( auto i = 0; i < size - 1; i++ ) {
            os << toPrint[i] << ", ";
        }
        os << toPrint[size - 1] << " )"
           << "\n";

        return os;
    }
};

template < class ValueType >
class JeveuxCollectionObject {
  public:
    typedef std::shared_ptr< JeveuxCollectionObjectClass< ValueType > >
        JeveuxCollectionObjectTypePtr;

  private:
    JeveuxCollectionObjectTypePtr _jeveuxCOPtr;

  public:
    /* Default constructor to be initialized with a null pointer
     * and really created later.
     */
    JeveuxCollectionObject() : _jeveuxCOPtr( nullptr ) {};

    JeveuxCollectionObject( const JeveuxObjectClass &coll, const ASTERINTEGER &number,
                            bool isNamed )
        : _jeveuxCOPtr( std::make_shared< JeveuxCollectionObjectClass< ValueType > >(
              coll, number, isNamed ) ) {};

    JeveuxCollectionObject( const JeveuxObjectClass &coll, const ASTERINTEGER &number,
                            const ASTERINTEGER &size )
        : _jeveuxCOPtr( std::make_shared< JeveuxCollectionObjectClass< ValueType > >( coll, number,
                                                                                      size ) ) {};

    JeveuxCollectionObject( const JeveuxObjectClass &coll, const ASTERINTEGER &number,
                            const std::string &objectName, const ASTERINTEGER &size )
        : _jeveuxCOPtr( std::make_shared< JeveuxCollectionObjectClass< ValueType > >(
              coll, number, objectName, size ) ) {};

    JeveuxCollectionObject( const JeveuxCollectionObject< ValueType > &tmp )
        : _jeveuxCOPtr( tmp._jeveuxCOPtr ) {};

    ~JeveuxCollectionObject() {};

    JeveuxCollectionObject &operator=( const JeveuxCollectionObject< ValueType > &tmp ) {
        _jeveuxCOPtr = tmp._jeveuxCOPtr;
        return *this;
    };

    const JeveuxCollectionObjectTypePtr &operator->( void ) const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _jeveuxCOPtr );
#endif

        return _jeveuxCOPtr;
    };

    JeveuxCollectionObjectClass< ValueType > &operator*( void ) const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _jeveuxCOPtr );
#endif
        return *_jeveuxCOPtr;
    };

    bool isEmpty() const {
        if ( _jeveuxCOPtr == nullptr )
            return true;
        if ( _jeveuxCOPtr.use_count() == 0 )
            return true;
        return false;
    };

    auto begin() const { return _jeveuxCOPtr->begin(); };

    auto end() const { return _jeveuxCOPtr->end(); };

    auto cbegin() const { return _jeveuxCOPtr->begin(); };

    auto cend() const { return _jeveuxCOPtr->end(); };

    inline const ValueType &operator[]( const ASTERINTEGER &i ) const {
        return ( *_jeveuxCOPtr )[i];
    };

    inline ValueType &operator[]( const ASTERINTEGER &i ) { return ( *_jeveuxCOPtr )[i]; };
};

/** @typedef Definition d'un objet de collection de type entier long */
typedef JeveuxCollectionObject< ASTERINTEGER > JeveuxCollectionObjectLong;
/** @typedef Definition d'un objet de collection de type entier court */
typedef JeveuxCollectionObject< ASTERINTEGER4 > JeveuxCollectionObjectShort;
/** @typedef Definition d'un objet de collection de type double */
typedef JeveuxCollectionObject< ASTERDOUBLE > JeveuxCollectionObjectReal;
/** @typedef Definition d'un objet de collection de type double complex */
typedef JeveuxCollectionObject< ASTERCOMPLEX > JeveuxCollectionObjectComplex;
/** @typedef Definition d'un objet de collection de JeveuxChar8 */
typedef JeveuxCollectionObject< JeveuxChar8 > JeveuxCollectionObjectChar8;
/** @typedef Definition d'un objet de collection de JeveuxChar16 */
typedef JeveuxCollectionObject< JeveuxChar16 > JeveuxCollectionObjectChar16;
/** @typedef Definition d'un objet de collection de JeveuxChar24 */
typedef JeveuxCollectionObject< JeveuxChar24 > JeveuxCollectionObjectChar24;
/** @typedef Definition d'un objet de collection de JeveuxChar32 */
typedef JeveuxCollectionObject< JeveuxChar32 > JeveuxCollectionObjectChar32;
/** @typedef Definition d'un objet de collection de JeveuxChar80 */
typedef JeveuxCollectionObject< JeveuxChar80 > JeveuxCollectionObjectChar80;
