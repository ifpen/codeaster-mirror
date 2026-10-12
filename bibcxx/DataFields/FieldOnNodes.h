#ifndef FIELDONNODES_H_
#define FIELDONNODES_H_

/**
 * @file FieldOnNodes.h
 * @brief Header of class for FieldOnNodes
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

#include "aster_fort_petsc.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"

#include "DataFields/DataField.h"
#include "DataFields/FieldConverter.h"
#include "MemoryManager/JeveuxAllowedTypes.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Messages/Messages.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/EquationNumbering.h"
#include "Numbering/ParallelDOFNumbering.h"
#include "ParallelUtilities/AsterMPI.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Blas.h"
#ifdef ASTER_HAVE_PETSC
#include <petscvec.h>
#endif
#include "ParallelUtilities/AsterMPI.h"

#include <filesystem>
#include <typeinfo>

/**
 * @struct AllowedFieldType
 * @brief Structure template permettant de limiter le type instanciable de
 * JeveuxVector
 * @tparam T Type autorise
 */
template < typename T >
struct AllowedFieldType; // undefined for bad types!

template <>
struct AllowedFieldType< ASTERINTEGER > {
    static const unsigned short numTypeJeveux = Integer;
};

template <>
struct AllowedFieldType< ASTERDOUBLE > {
    static const unsigned short numTypeJeveux = Real;
};

template <>
struct AllowedFieldType< ASTERCOMPLEX > {
    static const unsigned short numTypeJeveux = Complex;
};

template <>
struct AllowedFieldType< JeveuxChar8 > {
    static const unsigned short numTypeJeveux = Char8;
};

class FieldBuilder;

/**
 * @class FieldOnNodes
 * @brief Template class for FieldOnnodes
 */
template < class ValueType >
class FieldOnNodes : public DataField, private AllowedFieldType< ValueType > {
  private:
    /** @brief Vecteur Jeveux '.REFE' */
    JeveuxVectorChar24 _reference;
    /** @brief Vecteur Jeveux '.VALE' */
    JeveuxVector< ValueType > _values;
    /** @brief Dof description */
    EquationNumberingPtr _dofDescription;

  public:
    /** @typedef FieldOnNodesPtr */
    typedef std::shared_ptr< FieldOnNodes > FieldOnNodesPtr;

    /**
     * @brief Constructor
     * @param name Jeveux name of the field
     */
    FieldOnNodes( const std::string name, const EquationNumberingPtr nume = nullptr )
        : DataField( name, "CHAM_NO" ),
          _reference( JeveuxVectorChar24( getName() + ".REFE" ) ),
          _values( JeveuxVector< ValueType >( getName() + ".VALE" ) ),
          _dofDescription( nume ) {};

    /** @brief Constructor with automatic name */
    FieldOnNodes() : FieldOnNodes( DataStructureNaming::getNewName() ) {};

    /** @brief Copy constructor */
    FieldOnNodes( const std::string &name, const FieldOnNodes &toCopy )
        : FieldOnNodes( name, toCopy.getDescription() ) {
        // JeveuxVector to be duplicated
        *( _reference ) = *( toCopy._reference );
        *( _values ) = *( toCopy._values );
        *( _title ) = *( toCopy._title );

        this->updateValuePointers();
    }

    /** @brief Move constructor */
    FieldOnNodes( FieldOnNodes &&other ) : DataField( std::move( other ) ) {
        // Pointers to be moved$
        _reference = other._reference;
        _values = other._values;
        _title = other._title;
        _dofDescription = other._dofDescription;

        this->updateValuePointers();
    }

    /**
     * @brief Copy constructor
     */
    FieldOnNodes( const FieldOnNodes &toCopy )
        : FieldOnNodes( DataStructureNaming::getNewName(), toCopy ) {};

    /**
     * @brief Constructor with DOFNumbering
     */
    FieldOnNodes( const FiniteElementDescriptorPtr fed, const std::string &locatedComponent )
        : FieldOnNodes() {

        if ( !fed ) {
            raiseAsterError( " FEDescriptor is null" );
        }

        // Create numbering from Local mode (see Cata)
        std::string dofName;
        if ( fed->getMesh()->isParallel() ) {
#ifdef ASTER_HAVE_MPI
            auto dofNume = std::make_shared< ParallelDOFNumbering >();
            dofNume->computeNumbering( { fed }, locatedComponent, false );

            _dofDescription = dofNume->getEquationNumbering();
#endif
        } else {
            auto dofNume = std::make_shared< DOFNumbering >();

            dofNume->computeNumbering( { fed }, locatedComponent, false );

            _dofDescription = dofNume->getEquationNumbering();
        }

        const auto intType = AllowedFieldType< ValueType >::numTypeJeveux;
        CALLO_VTCREB_WRAP( getName(), JeveuxMemoryTypesNames[Permanent], JeveuxTypesNames[intType],
                           _dofDescription->getName() );

        this->updateValuePointers();
    };

    /**
     * @brief Constructor with DOFNumbering
     */
    FieldOnNodes( const ModelPtr &model ) : FieldOnNodes() {

        if ( !model ) {
            raiseAsterError( " Model is null" );
        }

        // Create numbering from Local mode (see Cata)
        auto dofNume = std::make_shared< DOFNumbering >();

        ListOfLoadsPtr lOL = std::make_shared< ListOfLoads >( model );

        dofNume->computeNumbering( model, lOL, false );

        _dofDescription = dofNume->getEquationNumbering();

        const auto intType = AllowedFieldType< ValueType >::numTypeJeveux;
        CALLO_VTCREB_WRAP( getName(), JeveuxMemoryTypesNames[Permanent], JeveuxTypesNames[intType],
                           _dofDescription->getName() );

        this->updateValuePointers();
    };

    /**
     * @brief Constructor with DOFNumbering
     */
    FieldOnNodes( const BaseDOFNumberingPtr &dofNum )
        : FieldOnNodes( dofNum->getEquationNumbering() ) {};

    /**
     * @brief Constructor with all attributes
     */
    FieldOnNodes( const EquationNumberingPtr &eqNum, const JeveuxVectorChar24 &reference,
                  const JeveuxVector< ValueType > &values )
        : FieldOnNodes( eqNum ) {
        *( _reference ) = *( reference );
        *( _values ) = *( values );

        this->updateValuePointers();
    };

    FieldOnNodes( const EquationNumberingPtr &equaNume ) : FieldOnNodes() {

        if ( !equaNume ) {
            raiseAsterError( "EquatioNumbering is null" );
        }

        _dofDescription = equaNume;

        const auto intType = AllowedFieldType< ValueType >::numTypeJeveux;
        CALLO_VTCREB_WRAP( getName(), JeveuxMemoryTypesNames[Permanent], JeveuxTypesNames[intType],
                           _dofDescription->getName() );

        _reference->updateValuePointer();
        const std::string name2 = strip( ( *_reference )[1].toString() );
        if ( strip( _dofDescription->getName() ) != name2 ) {
            _dofDescription = std::make_shared< EquationNumbering >( name2 );
            _dofDescription->setMesh( equaNume->getMesh() );
        }

        this->updateValuePointers();
    };

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    FieldOnNodes( const py::tuple &tup )
        : FieldOnNodes( tup[0].cast< std::string >(), tup[1].cast< EquationNumberingPtr >() ) {
        this->updateValuePointers();
    };

    py::tuple _getState() const { return py::make_tuple( getName(), getDescription() ); };

    /**
     * @brief Wrap of copy constructor
     * @return new field, copy of the calling field
     */
    FieldOnNodes copy() { return *this; }

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    ValueType &operator[]( const ASTERINTEGER &i ) { return _values->operator[]( i ); };

    /**
     * @brief Surcharge de l'operateur []
     * @param i Indice dans le tableau Jeveux
     * @return la valeur du tableau Jeveux a la position i
     */
    const ValueType &operator[]( const ASTERINTEGER &i ) const { return _values->operator[]( i ); };

    /**
     * @brief Check if fields are OK for +, +=, ...
     * @return true if compatible
     */
    bool isSimilarTo( const FieldOnNodes< ValueType > &tmp2 ) const {
        CALL_JEMARQ();
        bool similar = ( this->_reference->size() == tmp2._reference->size() );
        similar = ( similar && ( this->_values->size() == tmp2._values->size() ) );
        similar = ( similar && ( this->getMesh() == tmp2.getMesh() ) );

        if ( similar ) {
            similar = ( similar && ( ( *_dofDescription ) == ( *tmp2._dofDescription ) ) );
        }
        CALL_JEDEMA();
        return similar;
    }

    /**
     * @brief Scale a vector by a diagonal matrix stored as a vector
     * @return Updated field
     */
    FieldOnNodes< ValueType > &scale( const VectorReal &vect ) {
        _values->updateValuePointer();
        if ( _values->size() != vect.size() )
            raiseAsterError( "Field and vector have incompatible shapes" );
        for ( std::size_t i = 0; i < vect.size(); ++i ) {
            ( *_values )[i] = ( *_values )[i] * vect[i];
        }
        return *this;
    };

    /**
     * @brief PlusEqual overloading
     * @return Updated field
     */
    FieldOnNodes< ValueType > &operator+=( FieldOnNodes< ValueType > const &rhs ) {
        if ( !this->isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        ( *_values ) += ( *rhs._values );
        return *this;
    };

    /**
     * @brief MinusEqual overloading
     * @return Updated field
     * @todo ajouter une vérification sur la structure des champs
     */
    FieldOnNodes< ValueType > &operator-=( FieldOnNodes< ValueType > const &rhs ) {
        if ( !this->isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        ( *_values ) -= ( *rhs._values );
        return *this;
    };

    /**
     * @brief TimesEqual overloading
     * @return Updated field
     */
    FieldOnNodes< ValueType > &operator*=( const ASTERDOUBLE &scal ) {
        ( *_values ) *= scal;

        return *this;
    };

    /**
     * @brief DivideEqual overloading
     * @return Updated field
     */
    FieldOnNodes< ValueType > &operator/=( const ASTERDOUBLE &scal ) {
        ( *_values ) /= scal;

        return *this;
    };

    FieldOnNodes< ValueType > operator/=( const FieldOnNodes< ValueType > &rhs ) {
        if ( !this->isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        ( *_values ) /= ( *rhs._values );

        return *this;
    };

    /**
     * @brief Unary Minus overloading
     * @return Updated field
     */
    FieldOnNodes< ValueType > operator-() const {
        FieldOnNodes< ValueType > tmp( *this );
        ( *tmp._values ) *= ValueType( -1 );
        return tmp;
    };

    // NOTE ON FOLLOWING OPERATOR OVERLOADING
    // friends methods defined inside class body are inline and are hidden from
    // non-ADL lookup passing lhs by value helps optimize chained a+b+c otherwise,
    // both parameters may be const references return the result by value in order
    // to use copy constructor

    /**
     * @brief Multiply by a scalar on right overloading
     * @return New field
     */
    friend FieldOnNodes< ValueType > operator*( FieldOnNodes< ValueType > lhs,
                                                const ASTERDOUBLE &scal ) {

        lhs *= scal;
        return lhs;
    };

    /**
     * @brief Multiply by a scalar on left overloading
     * @return New field
     */
    friend FieldOnNodes< ValueType > operator*( const ASTERDOUBLE &scal,
                                                const FieldOnNodes< ValueType > &rhs ) {
        return rhs * scal;
    };

    /**
     * @brief Dividing by a scalar on right overloading
     * @return New field
     */
    friend FieldOnNodes< ValueType > operator/( FieldOnNodes< ValueType > lhs,
                                                const ASTERDOUBLE &scal ) {

        lhs /= scal;
        return lhs;
    };

    /**
     * @brief Plus overloading
     * @return New field
     */
    friend FieldOnNodes< ValueType > operator+( FieldOnNodes< ValueType > lhs,
                                                const FieldOnNodes< ValueType > &rhs ) {
        if ( !lhs.isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        lhs += rhs;
        return lhs;
    };

    /**
     * @brief Minus overloading
     * @return New field
     * @todo ajouter une vérification sur la structure des champs
     */
    friend FieldOnNodes< ValueType > operator-( FieldOnNodes< ValueType > lhs,
                                                const FieldOnNodes< ValueType > &rhs ) {
        if ( !lhs.isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        lhs -= rhs;
        return lhs;
    };

    std::string getLocalization() const { return "NOEU"; };

    VectorString getComponents() const { return _dofDescription->getComponents(); };

    ASTERINTEGER getNumberOfComponents() const { return getComponents().size(); }

    bool exists() const { return _reference.exists() && _values.exists(); };

    /**
     * @brief Get mesh
     */
    BaseMeshPtr getMesh() const {
        if ( _dofDescription ) {
            return _dofDescription->getMesh();
        }

        return nullptr;
    };

    bool printMedFile( const std::filesystem::path &fileName, bool local = true,
                       std::string version = std::string() ) const;

    /**
     * @brief Import a PETSc vector into a ParallelFieldOnNodes
     *
     * @param dofNmbrg The numbering of the DOFs
     * @param vec The PETSc vector
     * @param scaling The scaling of the Lagrange DOFs
     */
#ifdef ASTER_HAVE_PETSC
    FieldOnNodes *fromPetsc( const Vec &vec, const ASTERDOUBLE &scaling = 1.0,
                             bool local = false ) {
        const ASTERINTEGER local2 = local ? 1 : 0;
        CALLO_VECT_ASSE_FROM_PETSC( getName(), _dofDescription->getName(), &vec, &scaling,
                                    &local2 );
        _values->updateValuePointer();

        return this;
    };
#endif

    /**
     * @brief Communicates the values of the ghosts DOFs
     *
     */
    FieldOnNodes *updateGhostValues() {
        CALLO_VECT_ASSE_UPDATE_GHOST_VALUES( getName(), _dofDescription->getName() );
        _values->updateValuePointer();

        return this;
    };

    void applyLagrangeScaling( const ValueType scaling ) {
        _values->updateValuePointer();

        const auto descr = _dofDescription->getNodeAndComponentIdFromDOF();

        auto nbDofs = this->size();

        for ( ASTERINTEGER dof = 0; dof < nbDofs; dof++ ) {
            if ( descr[dof].second < 0 ) {
                ( *this )[dof] = ( *this )[dof] * scaling;
            }
        }
    };

    /**
     * @brief Set the Values object
     *
     * @param value Value to affect
     */
    void setValues( const ValueType &value ) {
        CALL_JEMARQ();
        _values->updateValuePointer();
        _values->assign( value );
        CALL_JEDEMA();
    };

    void setValues( const std::vector< ValueType > &values ) {
        AS_ASSERT( values.size() == size() );

        *_values = values;
    };

    void setValues( const std::map< std::string, ValueType > &values,
                    VectorString groupsOfNodes = {} ) {
        _values->updateValuePointer();

        auto num2name = _dofDescription->getComponentsIdToName();
        const auto descr = _dofDescription->getNodeAndComponentIdFromDOF();

        auto nbDofs = this->size();

        auto nodes = this->getMesh()->getNodes( groupsOfNodes, true );
        SetInt nodes_set;
        std::copy( nodes.begin(), nodes.end(), std::inserter( nodes_set, nodes_set.end() ) );

        for ( ASTERINTEGER dof = 0; dof < nbDofs; dof++ ) {
            auto search = values.find( num2name[descr[dof].second] );
            if ( search != values.end() ) {
                if ( nodes_set.count( descr[dof].first ) > 0 ) {
                    ( *this )[dof] = search->second;
                }
            }
        }
    };

    /**
     * @brief Get the values of the field
     *
     */
    const JeveuxVector< ValueType > &getValues() const { return _values; }

    /**
     * @brief Get the title of the field
     *
     */
    const JeveuxVectorChar80 &getTitle() const { return _title; }

    /**
     * @brief Get the reference of the field
     *
     */
    const JeveuxVectorChar24 &getReference() const { return _reference; }

    /**
     * @brief Get the equation numbering of the field
     *
     */
    const EquationNumberingPtr &getEquationNumbering() const { return _dofDescription; }

    /**
     * @brief Get the values of the field for some components and groups
     * @param cmps the components to extract
     * @param groupsOfNodes the nodes to extract
     *
     */
    std::vector< ValueType > getValues( const VectorString &cmps,
                                        const VectorString &groupsOfNodes = {} ) const {
        std::vector< ValueType > val;

        auto usedDof = _dofDescription->getDOFs( false, cmps, groupsOfNodes );
        val.reserve( usedDof.size() );
        _values->updateValuePointer();

        for ( auto &dof : usedDof ) {
            val.push_back( ( *_values )[dof] );
        }

        return val;
    }

    /**
     * @brief Get the values of the field for some dofs
     * @param cmps the dofs to extract
     *
     */
    std::vector< ValueType > getValues( const VectorLong &dofs ) const {
        std::vector< ValueType > result;
        result.reserve( dofs.size() );
        _values->updateValuePointer();
        for ( auto &dof : dofs )
            result.push_back( ( *_values )[dof] );
        return result;
    }

    /**
     * @brief Set FieldOnNodes description
     * @param desc object EquationNumberingPtr
     */
    void setDescription( const EquationNumberingPtr &desc ) {
        if ( !desc )
            raiseAsterError( "Empty EquationNumbering" );
        if ( _dofDescription && _dofDescription != desc )
            raiseAsterError( "EquationNumbering inconsistents" );
        if ( _reference.exists() ) {
            _reference->updateValuePointer();
            const auto descName = std::string( ( *_reference )[1].toString(), 0, 19 );
            if ( descName != desc->getName() ) {
                std::string mess;
                mess = "Description are incompatible: " + descName + " vs " + desc->getName();
                raiseAsterError( mess );
            }
        }
        _dofDescription = desc;
    };

    FieldOnNodesPtr restrict( const VectorString &cmps = {}, const VectorString &groupsOfNodes = {},
                              const ASTERINTEGER same_rank = PythonBool::None ) const {

        auto simpField = toSimpleFieldOnNodes( *this );
        auto simpFieldRest = simpField->restrict( cmps, groupsOfNodes, same_rank );

        return toFieldOnNodes( simpFieldRest );
    };

    FieldOnNodesPtr asPhysicalQuantity( const std::string physQuant,
                                        const MapString &map_cmps ) const {

        auto simpField = toSimpleFieldOnNodes( *this );
        auto simpFieldRest = simpField->asPhysicalQuantity( physQuant, map_cmps );

        return toFieldOnNodes( simpFieldRest );
    };

    /**
     * @brief Wrap of copy constructor
     * @param desc Description that is used for the copy
     * @param warn If set to true, raises an alarm if values are set to zero
     * @return new field, copy of the calling field
     */
    FieldOnNodesPtr copyUsingDescription( const EquationNumberingPtr desc,
                                          const bool warn = true ) {
        auto field = std::make_shared< FieldOnNodes< ValueType > >( desc );

        ASTERINTEGER iret = -1;
        CALLO_VTCOPY( getName(), field->getName(), &iret );

        if ( warn && iret != 0 ) {
            UTMESS( "A", "FIELD0_30" );
        }

        return field;
    }

    /**
     * @brief Set mesh
     * @param mesh object BaseMeshPtr
     */
    void setMesh( const BaseMeshPtr &mesh ) {
        AS_ASSERT( _dofDescription );
        _dofDescription->setMesh( mesh );
    };

    /**
     * @brief Compute norm
     * @param normType Type of norm ("NORM_1","NORM_2","NORM_INFINITY")
     */
    ASTERDOUBLE norm( const std::string normType, VectorString list_cmp = VectorString() ) const {

        if constexpr ( !std::is_same_v< ValueType, ASTERDOUBLE > &&
                       !std::is_same_v< ValueType, ASTERINTEGER > &&
                       !std::is_same_v< ValueType, ASTERCOMPLEX > ) {
            raiseAsterError( " norm method not defined for type " +
                             std::string( typeid( ValueType ).name() ) );
        }
        CALL_JEMARQ();
        ASTERDOUBLE norme = 0.0;
        _values->updateValuePointer();
        auto dofUsed = _dofDescription->getDOFs( true, list_cmp );

        if ( normType == "NORM_1" ) {
            for ( auto &dof : dofUsed ) {
                norme += std::abs( ( *this )[dof] );
            }
        } else if ( normType == "NORM_2" ) {
            for ( auto &dof : dofUsed ) {
                norme += std::pow( std::abs( ( *this )[dof] ), 2 );
            }
        } else if ( normType == "NORM_INFINITY" ) {
            for ( auto &dof : dofUsed ) {
                norme = std::max( norme, std::abs( ( *this )[dof] ) );
            }
        } else
            raiseAsterError( "Unknown norm: " + normType );

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

    /**
     * @brief Dot product
     * @param tmp object EquationNumberingPtr
     */
    ValueType dot( const FieldOnNodesPtr &tmp ) const {
        CALL_JEMARQ();
        tmp->updateValuePointers();
        _values->updateValuePointer();
        const auto taille = _values->size();

        if ( taille != tmp->size() )
            raiseAsterError( "Incompatible size" );

        auto dofUsed = _dofDescription->getDOFs( true );

        ValueType ret;
        if constexpr ( std::is_same_v< ValueType, ASTERDOUBLE > ||
                       std::is_same_v< ValueType, ASTERINTEGER > ) {
            ret = (ValueType)0.0;
            for ( auto &dof : dofUsed ) {
                ret += ( *this )[dof] * ( *tmp )[dof];
            }
        } else if constexpr ( std::is_same_v< ValueType, ASTERCOMPLEX > ) {
            ret = (ValueType)0.0;
            for ( auto &dof : dofUsed ) {
                ret += ( *this )[dof] * std::conj( ( *tmp )[dof] );
            }
        } else {
            raiseAsterError( " dot method not defined for type " +
                             std::string( typeid( ValueType ).name() ) );
        }

#ifdef ASTER_HAVE_MPI
        if ( this->getMesh()->isParallel() ) {
            ret = AsterMPI::sum( ret );
        }
#endif
        CALL_JEDEMA();
        return ret;
    }

    /**
     * @brief Get physical quantity
     */
    std::string getPhysicalQuantity() const {
        const std::string typeco( "CHAM_NO" );
        ASTERINTEGER repi = 0, ier = 0;
        JeveuxChar32 repk( " " );
        const std::string arret( "F" );
        const std::string questi( "NOM_GD" );
        CALLO_DISMOI( questi, this->getName(), typeco, &repi, repk, arret, &ier );
        auto retour = strip( repk.toString() );
        return retour;
    }

    /**
     * @brief Size of the FieldOnNodes
     */
    ASTERINTEGER size( void ) const { return _values->size(); }

    /**
     * @brief Get EquationNumbering
     */
    EquationNumberingPtr getDescription( void ) const { return _dofDescription; };

    /**
     * @brief Update field and build EquationNumbering if necessary
     */
    bool build( const BaseMeshPtr mesh = nullptr ) {
        if ( !_dofDescription ) {
            CALL_JEMARQ();

            _reference->updateValuePointer();
            const std::string name2 = strip( ( *_reference )[1].toString() );
            if ( !name2.empty() ) {
                if ( mesh ) {
                    if ( mesh->isParallel() ) {
#ifdef ASTER_HAVE_MPI
                        _dofDescription = std::make_shared< ParallelEquationNumbering >( name2 );
#endif
                    } else {
                        _dofDescription = std::make_shared< EquationNumbering >( name2 );
                    }
                    this->setMesh( mesh );
                } else {
                    _dofDescription = std::make_shared< EquationNumbering >( name2 );
                }
            } else {
                AS_ABORT( "NUME_EQUA is empty" );
            }
            CALL_JEDEMA();
        }
        if ( this->exists() ) {
            this->updateValuePointers();
        }
        return true;
    };

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() const {
        _reference->updateValuePointer();
        _values->updateValuePointer();
    };

    /**
     * @brief Transformer les valeurs de _values en appliquant
     *        la fonction "func" à chaque valeur
     * @return renvoie un nouvel objet de FieldOnNodes
     *         avec les valeurs transformées
     */
    FieldOnNodes< ValueType > transform( py::object &func ) const;

    friend class FieldBuilder;
};

template < class ValueType >
bool FieldOnNodes< ValueType >::printMedFile( const std::filesystem::path &fileName, bool local,
                                              std::string version ) const {
    const auto rank = getMPIRank();
    LogicalUnitFile a;
    ASTERINTEGER retour = -1;
    // In case that the print file (single and absolute path) is unique between processors,
    // it must only be created on proc 0.
    if ( getMesh()->isParallel() || ( !getMesh()->isParallel() && rank == 0 ) ) {
        if ( rank == 0 )
            a.openFile( fileName, Binary, New );
        if ( getMesh()->isParallel() ) {
#ifdef ASTER_HAVE_MPI
            AsterMPI::barrier();
#endif /* ASTER_HAVE_MPI */
            if ( rank != 0 )
                a.openFile( fileName, Binary, Old );
        }
        retour = a.getLogicalUnit();
    }

    CommandSyntax cmdSt( "IMPR_RESU" );

    SyntaxMapContainer dict;
    dict.container["FORMAT"] = "MED";
    dict.container["UNITE"] = (ASTERINTEGER)retour;
    if ( !version.empty() ) {
        dict.container["VERSION_MED"] = version;
    }

    if ( getMesh()->isParallel() ) {
        dict.container["PROC0"] = "NON";
        if ( !local )
            dict.container["FICHIER_UNIQUE"] = "OUI";
    } else
        dict.container["PROC0"] = "OUI";

    ListSyntaxMapContainer listeResu;
    SyntaxMapContainer dict2;
    dict2.container["CHAM_GD"] = getName();
    listeResu.push_back( dict2 );
    dict.container["RESU"] = listeResu;

    cmdSt.define( dict );

    ASTERINTEGER op = 39;
    CALL_EXECOP( &op );

    return true;
};

/** @typedef FieldOnNodesReal */
using FieldOnNodesReal = FieldOnNodes< ASTERDOUBLE >;
using FieldOnNodesRealPtr = std::shared_ptr< FieldOnNodesReal >;

/** @typedef FieldOnNodesLong */
using FieldOnNodesLong = FieldOnNodes< ASTERINTEGER >;
using FieldOnNodesLongPtr = std::shared_ptr< FieldOnNodesLong >;

/** @typedef FieldOnNodesComplex*/
using FieldOnNodesComplex = FieldOnNodes< ASTERCOMPLEX >;
using FieldOnNodesComplexPtr = std::shared_ptr< FieldOnNodesComplex >;

/** @typedef FieldOnNodesChar8 */
using FieldOnNodesChar8 = FieldOnNodes< JeveuxChar8 >;
using FieldOnNodesChar8Ptr = std::shared_ptr< FieldOnNodesChar8 >;

#endif /* FIELDONNODES_H_ */
