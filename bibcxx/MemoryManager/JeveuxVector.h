#ifndef JEVEUXVECTOR_H_
#define JEVEUXVECTOR_H_

/**
 * @file JeveuxVector.h
 * @brief Fichier entete de la classe JeveuxVector
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

#include "aster_fort_jeveux.h"
#include "aster_utils.h"
#include "shared_vars.h"

#include "MemoryManager/JeveuxAllowedTypes.h"
#include "MemoryManager/JeveuxObject.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/Blas.h"

/**
 * @class JeveuxVectorClass
 * @brief Cette classe template permet de definir un vecteur Jeveux
 * @author Nicolas Sellenet
 * @todo rajouter un constructeur de nom jeveux pour que .NSLV soit bien placé (cf DOFNumbering.cxx)
 */
template < typename ValueType >
class JeveuxVectorClass : public JeveuxObjectClass, private AllowedJeveuxType< ValueType > {
  private:
    /** @brief Pointeur vers la premiere position du vecteur Jeveux */
    ValueType *_valuePtr;
    /** @brief Adresse du vecteur Jeveux */
    ASTERINTEGER _jeveuxAdress;

  public:
    /**
     * @brief Constructeur
     * @param name Nom jeveux du vecteur
     *   Attention, le pointeur est mis a zero. Avant d'utiliser ce vecteur,
     *   il faut donc faire appel a JeveuxVectorClass::updateValuePointer
     */
    JeveuxVectorClass( const std::string &nom )
        : JeveuxObjectClass( nom ),
          _valuePtr( nullptr ),
          _jeveuxAdress( 0 ) {
              // #ifdef ASTER_DEBUG_CXX_OBJECTS
              //         std::cout << "DEBUG: JeveuxVector.create: " << _name << std::endl;
              // #endif
          };

    /**
     * @brief Destructeur
     */
    ~JeveuxVectorClass() {
        // #ifdef ASTER_DEBUG_CXX_OBJECTS
        //         std::cout << "DEBUG: JeveuxVector.destr: " << _name << std::endl;
        // #endif
        this->deallocate();
    };

    /**
     * @brief Surcharge de l'operateur =
     */
    JeveuxVectorClass &operator=( JeveuxVectorClass< ValueType > &toCopy ) {
        CALL_JEMARQ();
        if ( this->size() != 0 )
            this->deallocate();
        const auto size = toCopy.size();
        const auto capa = toCopy.capacity();

        if ( capa > 0 ) {
            this->allocate( capa );
            this->setSize( size );
            toCopy.updateValuePointer();
            for ( ASTERINTEGER i = 0; i < size; ++i )
                this->operator[]( i ) = toCopy[i];
            // Copy Jeveux attribute
            std::string docu = toCopy.getInformationParameter();
            if ( docu != "" )
                this->setInformationParameter( docu );
        }
        CALL_JEDEMA();

        return *this;
    };

    /**
     * @brief Surcharge de l'operateur =
     */
    JeveuxVectorClass &operator=( const std::vector< ValueType > &toCopy ) {
        CALL_JEMARQ();

        if ( this->size() != 0 )
            this->deallocate();
        const ASTERINTEGER size = toCopy.size();
        if ( size > 0 ) {
            this->allocate( size );
            for ( ASTERINTEGER i = 0; i < size; ++i )
                this->operator[]( i ) = toCopy[i];
        }
        CALL_JEDEMA();

        return *this;
    };

    /**
     * @brief Surcharge de l'operateur =
     */
    JeveuxVectorClass &operator=( const ValueType &val ) {
        CALL_JEMARQ();
        this->updateValuePointer();
        this->assign( val );
        CALL_JEDEMA();

        return *this;
    };

    /**
     * @brief Surcharge de l'operateur =
     */
    bool operator==( JeveuxVectorClass< ValueType > &toCompare ) {
        if ( this->size() != toCompare.size() )
            return false;
        CALL_JEMARQ();
        bool ret = true;

        toCompare.updateValuePointer();
        this->updateValuePointer();
        const auto size = toCompare.size();
        for ( ASTERINTEGER i = 0; i < size; ++i ) {
            if ( this->operator[]( i ) != toCompare[i] ) {
                ret = false;
                break;
            }
        }

        CALL_JEDEMA();

        return ret;
    };

    /**
     * @brief Surcharge de l'operateur ==
     */
    bool operator==( const std::vector< ValueType > &toCompare ) {
        if ( this->size() != toCompare.size() )
            return false;

        CALL_JEMARQ();
        bool ret = true;
        this->updateValuePointer();

        const ASTERINTEGER size = toCompare.size();
        for ( ASTERINTEGER i = 0; i < size; ++i ) {
            if ( this->operator[]( i ) != toCompare[i] ) {
                ret = false;
                break;
            }
        }

        CALL_JEDEMA();

        return ret;
    };

    /**
     * @brief Surcharge de l'operateur [] avec des const
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    inline const ValueType &operator[]( const ASTERINTEGER &i ) const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
        AS_ASSERT( !this->hasMoved() );
        if ( i < 0 || i >= this->size() ) {
            std::string error = "Out of range of JeveuxVector '" + _name +
                                "', index = " + std::to_string( i ) +
                                " ( size = " + std::to_string( this->size() ) + " )";
            AS_ABORT( error );
        }
#endif
        return _valuePtr[i];
    };

    /**
     * @brief Surcharge de l'operateur [] sans const (pour les lvalue)
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    inline ValueType &operator[]( const ASTERINTEGER &i ) {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
        AS_ASSERT( !this->hasMoved() );
        if ( i < 0 || i >= this->size() ) {
            std::string error = "Out of range of JeveuxVector '" + _name +
                                "', index = " + std::to_string( i ) +
                                " ( size = " + std::to_string( this->size() ) + " )";
            AS_ABORT( error );
        }
#endif
        return _valuePtr[i];
    };

    /**
     * @brief Fonction d'allocation d'un vecteur Jeveux
     * @param length Longueur du vecteur Jeveux a allouer
     * @return true si l'allocation s'est bien passee
     */
    void allocate( const ASTERINTEGER &length ) {
        bool ret = false;
        if ( _name != "" && length > 0 ) {
            ret = true;
            std::string strJeveuxBase = JeveuxMemoryTypesNames[_mem];
            const auto intType = AllowedJeveuxType< ValueType >::numTypeJeveux;
            std::string carac = strJeveuxBase + " V " + JeveuxTypesNames[intType];
            ASTERINTEGER taille = length;
            CALLO_WKVECTC( _name, carac, &taille, (void *)( &_valuePtr ) );
            if ( _valuePtr == NULL )
                ret = false;
        } else {
            AS_ABORT( "JeveuxVector " + _name + " allocate null size" );
        }

        AS_ASSERT( ret );

        updateValuePointer();
    };

    void allocate( const ASTERINTEGER &length, const ValueType &val ) {
        this->allocate( length );
        this->assign( val );
    };

    /**
     * @brief Desallocation d'un vecteur Jeveux
     */
    void deallocate() {
        // #ifdef ASTER_DEBUG_CXX
        //         std::cout << "DEBUG: JeveuxVector.dealloc: " << _name << std::endl;
        // #endif
        if ( _name != "" && get_sh_jeveux_status() == 1 )
            CALLO_JEDETR( _name );

        _jeveuxAdress = 0;
        _valuePtr = NULL;
        // #ifdef ASTER_DEBUG_CXX
        //         std::cout << "DEBUG: JeveuxVector.dealloc: " << _name << std::endl;
        // #endif
    };

    /**
     * @brief Return a pointer to the vector
     */
    ValueType *getDataPtr() { return _valuePtr; };

    /**
     * @brief Return a pointer to the vector
     */
    const ValueType *getDataPtr() const { return _valuePtr; };

    /**
     * @brief Get the value of DOCU parameter of jeveux object
     */
    std::string getInformationParameter() const {
        if ( !exists() )
            return "";
        const std::string param( "DOCU" );
        std::string charval( 4, ' ' );
        ASTERINTEGER valTmp;
        CALLO_JELIRA( _name, param, &valTmp, charval );
        std::string toReturn( charval );
        return toReturn;
    };

    /**
     * @brief Return the name
     */
    std::string getName() const { return _name; };

    /**
     * @brief Set the value of DOCU parameter of jeveux object
     */
    void setInformationParameter( const std::string value ) {
        AS_ASSERT( exists() );

        const std::string param( "DOCU" );
        CALLO_JEECRA_STRING_WRAP( _name, param, value );
    };

    /**
     * @brief Set the number of value used in the vector (=LONUTI)
     */
    void setSize( ASTERINTEGER value ) {
        AS_ASSERT( exists() );

        const std::string param( "LONUTI" );
        CALLO_JEECRA_WRAP( _name, param, &value );
        AS_ASSERT( this->size() <= this->capacity() );
    };

    /**
     * @brief Return the capacity of the vector (=LONMAX)
     */
    ASTERINTEGER capacity() const {
        if ( !exists() )
            return 0;

        ASTERINTEGER vectSize;
        JeveuxChar8 param( "LONMAX" );
        JeveuxChar32 dummy( " " );
        CALLO_JELIRA( _name, param, &vectSize, dummy );
        return vectSize;
    };

    /**
     * @brief Return the size of the vector (=LONUTI)
     */
    ASTERINTEGER size() const {
        if ( !exists() )
            return 0;

        ASTERINTEGER vectSize;
        JeveuxChar8 param( "LONUTI" );
        JeveuxChar32 dummy( " " );
        CALLO_JELIRA( _name, param, &vectSize, dummy );

        return vectSize;
    };

    /**
     * @brief Mise a jour du pointeur Jeveux
     * @return true si la mise a jour s'est bien passee
     */
    void updateValuePointer() const {
        if ( hasMoved() ) {
            const_cast< JeveuxVectorClass< ValueType > * >( this )->_valuePtr = NULL;
            const_cast< JeveuxVectorClass< ValueType > * >( this )->_jeveuxAdress = 0;

            const std::string read( "L" );
            CALLO_JEVEUOC( _name, read, (void *)( &_valuePtr ) );
            AS_ASSERT( _valuePtr );

            JeveuxChar8 param( "IADM" );
            JeveuxChar32 dummy( " " );
            CALLO_JELIRA( _name, param,
                          &const_cast< JeveuxVectorClass< ValueType > * >( this )->_jeveuxAdress,
                          dummy );
        }
    };

    /**
     * @brief Est-ce que le pointeur est toujours bon ?
     * @return true si pointeur toujours valide
     */
    bool hasMoved() const {
        if ( !_valuePtr ) {
            return true;
        }

        // check if object is in memory
        ASTERINTEGER idummy;
        JeveuxChar8 param0( "USAGE" );
        JeveuxChar32 usage( " " );
        CALLO_JELIRA( _name, param0, &idummy, usage );

        if ( usage[0] == 'X' ) {
            const_cast< JeveuxVectorClass< ValueType > * >( this )->_valuePtr = NULL;
            return true;
        }

        ASTERINTEGER valTmp;
        JeveuxChar8 param( "IADM" );
        JeveuxChar32 dummy( " " );
        CALLO_JELIRA( _name, param, &valTmp, dummy );
        if ( valTmp != _jeveuxAdress ) {
            const_cast< JeveuxVectorClass< ValueType > * >( this )->_valuePtr = NULL;
            return true;
        }

        return false;
    };

    /** @brief Convert to std::vector */
    std::vector< ValueType > toVector() const {
        if ( !exists() ) {
            return std::vector< ValueType >();
        }

        return slice( 0, size() );
    };

    /** @brief Vector containing the first n elements */
    std::vector< ValueType > head( const ASTERINTEGER &n ) const { return slice( 0, n ); };

    /** @brief Vector containing n elements, starting at position i */
    std::vector< ValueType > slice( const ASTERINTEGER &first, const ASTERINTEGER &n ) const {

        CALL_JEMARQ();

        updateValuePointer();

        std::vector< ValueType > toReturn;
        toReturn.reserve( n );

        const ASTERINTEGER total = first + n;

        for ( ASTERINTEGER i = first; i < total; ++i ) {
            toReturn.push_back( this->operator[]( i ) );
        }

        CALL_JEDEMA();

        return toReturn;
    };

    /** @brief checks whether the container is empty */
    bool empty() const {
        if ( this->size() == 0 )
            return true;

        return false;
    }

    /** @brief changes the number of elements stored
     * @param size new size of the vector
     */
    void resize( const ASTERINTEGER &size ) {
        if ( !this->exists() ) {
            this->allocate( size );
        } else if ( size > this->capacity() ) {
            ASTERINTEGER taille = size;
            CALLO_JUVECA( _name, &taille );
            updateValuePointer();
        }

        this->setSize( size );
    };

    /** @brief adds an element to the end
     * @param elem element to add
     */
    void push_back( const ValueType &elem ) {
        const auto size = this->size();
        const auto capa = this->capacity();

        if ( ( size + 1 ) > capa )
            this->resize( 2 * ( size + 1 ) );

        _valuePtr[size] = elem;
        this->setSize( size + 1 );
    };

    /** @brief replace in-place an element to the end
     * @param elem element to add
     */
    void emplace_back( const ValueType &elem ) {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
#endif

        _valuePtr[this->size() - 1] = elem;
    };

    /** @brief reserves storage. does not change the size of the vector.
     * @param capacity new capacity.
     */
    void reserve( const ASTERINTEGER &capacity ) {
        if ( this->exists() ) {
            const auto size = this->size();
            this->resize( capacity );
            this->setSize( size );
        } else {
            this->allocate( capacity );
            this->setSize( 0 );
        }
    };

    /** @brief clears the contents */
    void clear() { this->deallocate(); }

    /** @brief removes the last element  */
    void pop_back() { this->setSize( this->size() - 1 ); }

    /** @brief access the last element */
    ValueType back() const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
        AS_ASSERT( this->size() > 0 );
#endif
        return _valuePtr[this->size() - 1];
    }

    /** @brief access the first element */
    ValueType front() const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
        AS_ASSERT( this->size() > 0 );
#endif
        return _valuePtr[0];
    }

    /** @brief assign the vector with a given value */
    void assign( const ValueType &val ) {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
#endif
        const auto size = this->size();
        for ( auto i = 0; i < size; i++ ) {
            this->operator[]( i ) = val;
        }
    }

    /** @brief overload << operator */
    friend std::ostream &operator<<( std::ostream &os,
                                     const JeveuxVectorClass< ValueType > &toPrint ) {
        toPrint.updateValuePointer();
        os << "JeveuxVector: " << toPrint.getName() << "\n";
        os << "Size: " << std::to_string( toPrint.size() )
           << ", and capacity: " << std::to_string( toPrint.capacity() ) << ".\n";

        const auto size = toPrint.size();
        os << "List of values: \n";
        if ( size > 0 ) {
            os << "( ";
            for ( auto i = 0; i < size - 1; i++ ) {
                os << toPrint[i] << ", ";
            }
            os << toPrint[size - 1] << " )"
               << "\n";
        }

        return os;
    }

    struct iterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = ValueType;
        using pointer = ValueType *;
        using reference = ValueType &;

        iterator( pointer ptr ) : m_ptr( ptr ) {}

        reference operator*() const { return *m_ptr; }

        pointer operator->() { return m_ptr; }

        pointer operator&() const { return m_ptr; }

        iterator &operator++() {
            m_ptr++;
            return *this;
        }

        friend bool operator==( const iterator &a, const iterator &b ) {
            return a.m_ptr == b.m_ptr;
        };

        friend bool operator!=( const iterator &a, const iterator &b ) {
            return a.m_ptr != b.m_ptr;
        };

        iterator operator+( difference_type count ) const { return iterator( m_ptr + count ); }

        // Subtraction: it - count
        iterator operator-( difference_type count ) const { return iterator( m_ptr - count ); }

      private:
        pointer m_ptr;
    };

    iterator begin() {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
#endif

        return iterator( ( &_valuePtr )[0] );
    };

    iterator end() {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
#endif

        return iterator( (ValueType *)( _valuePtr + this->size() ) );
    };

    struct const_iterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = ValueType;
        using pointer = ValueType *;
        using reference = ValueType &;

        const_iterator( pointer ptr ) : m_ptr( ptr ) {}

        reference operator*() const { return *m_ptr; }

        pointer operator->() const { return m_ptr; }

        pointer operator&() const { return m_ptr; }

        const_iterator &operator++() {
            m_ptr++;
            return *this;
        }

        const_iterator operator+( difference_type count ) const {
            return const_iterator( m_ptr + count );
        }

        // Subtraction: it - count
        const_iterator operator-( difference_type count ) const {
            return const_iterator( m_ptr - count );
        }

        friend bool operator==( const const_iterator &a, const const_iterator &b ) {
            return a.m_ptr == b.m_ptr;
        };

        friend bool operator!=( const const_iterator &a, const const_iterator &b ) {
            return a.m_ptr != b.m_ptr;
        };

      private:
        pointer m_ptr;
    };

    const_iterator cbegin() {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
#endif

        return const_iterator( ( &_valuePtr )[0] );
    }

    const_iterator cend() {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        AS_ASSERT( _valuePtr != nullptr );
#endif

        return const_iterator( (ValueType *)( _valuePtr + this->size() ) );
    }

    /**
     * @brief TimesEqual overloading
     * @return Updated JeveuxVector
     */
    JeveuxVectorClass< ValueType > &operator*=( const ValueType &scal ) {
        CALL_JEMARQ();
        this->updateValuePointer();
        const auto size = this->size();

        AsterBLAS::scal( size, scal, getDataPtr(), ASTERINTEGER( 1 ) );

        CALL_JEDEMA();

        return *this;
    };

    JeveuxVectorClass< ValueType > &operator/=( const ValueType &scal ) {
        CALL_JEMARQ();
        this->updateValuePointer();
        const auto size = this->size();

        AsterBLAS::scal( size, 1.0 / scal, getDataPtr(), ASTERINTEGER( 1 ) );

        CALL_JEDEMA();

        return *this;
    };

    JeveuxVectorClass< ValueType > &operator/=( const JeveuxVectorClass< ValueType > &other ) {
        CALL_JEMARQ();
        this->updateValuePointer();
        const auto size = this->size();
        AS_ASSERT( size == other.size() );

        other.updateValuePointer();
        for ( ASTERINTEGER i = 0; i < size; ++i )
            this->operator[]( i ) = this->operator[]( i ) / other.operator[]( i );

        CALL_JEDEMA();

        return *this;
    };

    /**
     * @brief MinusEqual overloading
     * @return Updated JeveuxVector
     */
    JeveuxVectorClass< ValueType > &operator-=( const JeveuxVectorClass< ValueType > &vect ) {
        CALL_JEMARQ();
        vect.updateValuePointer();
        this->updateValuePointer();
        const auto size = this->size();

        if ( size != vect.size() ) {
            AS_ABORT( "Incompatible sizes: " + std::to_string( size ) + " vs " +
                      std::to_string( vect.size() ) );
        }

        for ( auto i = 0; i < size; i++ ) {
            this->operator[]( i ) -= vect[i];
        }

        CALL_JEDEMA();

        return *this;
    };

    /**
     * @brief MinusEqual overloading
     * @return Updated JeveuxVector
     */
    JeveuxVectorClass< ValueType > &operator-=( const ValueType &val ) {
        CALL_JEMARQ();
        this->updateValuePointer();
        const auto size = this->size();

        for ( auto i = 0; i < size; i++ ) {
            this->operator[]( i ) -= val;
        }

        CALL_JEDEMA();

        return *this;
    };

    /**
     * @brief PlusEqual overloading
     * @return Updated JeveuxVector
     */
    JeveuxVectorClass< ValueType > &operator+=( const JeveuxVectorClass< ValueType > &vect ) {
        CALL_JEMARQ();
        vect.updateValuePointer();
        this->updateValuePointer();
        const auto size = this->size();

        if ( size != vect.size() ) {
            AS_ABORT( "Incompatible sizes: " + std::to_string( size ) + " vs " +
                      std::to_string( vect.size() ) );
        }

        for ( auto i = 0; i < size; i++ ) {
            this->operator[]( i ) += vect[i];
        }

        CALL_JEDEMA();

        return *this;
    };

    /**
     * @brief PlusEqual overloading
     * @return Updated JeveuxVector
     */
    JeveuxVectorClass< ValueType > &operator+=( const ValueType &val ) {
        CALL_JEMARQ();
        this->updateValuePointer();
        const auto size = this->size();

        for ( auto i = 0; i < size; i++ ) {
            this->operator[]( i ) += val;
        }

        CALL_JEDEMA();

        return *this;
    };
};

/**
 * @class JeveuxVector
 * @brief Enveloppe d'un pointeur intelligent vers un JeveuxVectorClass
 * @author Nicolas Sellenet
 * @todo Supprimer la classe enveloppe
 */
template < class ValueType >
class JeveuxVector {
  public:
    typedef std::shared_ptr< JeveuxVectorClass< ValueType > > JeveuxVectorTypePtr;

  private:
    JeveuxVectorTypePtr _jeveuxVectorPtr;

  public:
    /* Default constructor to be initialized with a null pointer
     * and really created later.
     */
    JeveuxVector() : _jeveuxVectorPtr( nullptr ) {};

    JeveuxVector( std::string nom )
        : _jeveuxVectorPtr( std::make_shared< JeveuxVectorClass< ValueType > >( nom ) ) {};

    JeveuxVector( std::string nom, const std::vector< ValueType > &vect ) : JeveuxVector( nom ) {
        ( *_jeveuxVectorPtr ) = vect;
    };

    JeveuxVector( std::string nom, const ASTERINTEGER &size ) : JeveuxVector( nom ) {
        _jeveuxVectorPtr->allocate( size );
    };

    JeveuxVector( std::string nom, const ASTERINTEGER &size, const ValueType &val )
        : JeveuxVector( nom ) {
        _jeveuxVectorPtr->allocate( size, val );
    };

    ~JeveuxVector() {};

    JeveuxVector &operator=( const JeveuxVector< ValueType > &tmp ) {
        _jeveuxVectorPtr = tmp._jeveuxVectorPtr;
        return *this;
    };

    const JeveuxVectorTypePtr &operator->() const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        if ( !_jeveuxVectorPtr ) {
            AS_ABORT( "Null Jeveux vector" );
        }
#endif

        return _jeveuxVectorPtr;
    };

    JeveuxVectorClass< ValueType > &operator*() const {
#ifdef ASTER_DEBUG_CXX_OBJECTS
        if ( !_jeveuxVectorPtr ) {
            AS_ABORT( "Null Jeveux vector" );
        }
#endif
        return *_jeveuxVectorPtr;
    };

    bool exists() const {

        if ( _jeveuxVectorPtr == nullptr )
            return false;
        if ( _jeveuxVectorPtr.use_count() == 0 )
            return false;

        return _jeveuxVectorPtr->exists();
    }

    auto begin() { return _jeveuxVectorPtr->begin(); };

    auto end() { return _jeveuxVectorPtr->end(); };

    const auto cbegin() { return _jeveuxVectorPtr->cbegin(); };

    const auto cend() { return _jeveuxVectorPtr->cend(); };
};

/** @typedef Definition d'un vecteur Jeveux entier long */
typedef JeveuxVector< ASTERINTEGER > JeveuxVectorLong;
/** @typedef Definition d'un vecteur Jeveux entier court */
typedef JeveuxVector< ASTERINTEGER4 > JeveuxVectorShort;
/** @typedef Definition d'un vecteur Jeveux double */
typedef JeveuxVector< ASTERDOUBLE > JeveuxVectorReal;
/** @typedef Definition d'un vecteur Jeveux double complex */
typedef JeveuxVector< ASTERCOMPLEX > JeveuxVectorComplex;
/** @typedef Definition d'un vecteur de JeveuxChar8 */
typedef JeveuxVector< JeveuxChar8 > JeveuxVectorChar8;
/** @typedef Definition d'un vecteur JeveuxChar16 */
typedef JeveuxVector< JeveuxChar16 > JeveuxVectorChar16;
/** @typedef Definition d'un vecteur JeveuxChar24 */
typedef JeveuxVector< JeveuxChar24 > JeveuxVectorChar24;
/** @typedef Definition d'un vecteur JeveuxChar32 */
typedef JeveuxVector< JeveuxChar32 > JeveuxVectorChar32;
/** @typedef Definition d'un vecteur JeveuxChar80 */
typedef JeveuxVector< JeveuxChar80 > JeveuxVectorChar80;
/** @typedef Definition d'un vecteur JeveuxLogical */
typedef JeveuxVector< ASTERBOOL > JeveuxVectorLogical;

#endif /* JEVEUXVECTOR_H_ */
