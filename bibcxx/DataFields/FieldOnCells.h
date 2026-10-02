#ifndef FIELDONCELLS_H_
#define FIELDONCELLS_H_

/**
 * @file FieldOnCells.h
 * @brief Header of class for FieldOnCells
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

#include "aster_fort_calcul.h"
#include "aster_fort_ds.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "DataFields/DataField.h"
#include "DataFields/SimpleFieldOnCells.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "ParallelUtilities/AsterMPI.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"

#include <filesystem>

/**
 * @class FieldOnCells
 * @brief Template class for FieldOnCells
 */
template < class ValueType >
class FieldOnCells : public DataField {
  private:
    using SimpleFieldOnCellsValueType = SimpleFieldOnCells< ValueType >;
    using SimpleFieldOnCellsValueTypePtr = std::shared_ptr< SimpleFieldOnCellsValueType >;

    /** @brief Vecteur Jeveux '.CELD' */
    JeveuxVectorLong _descriptor;
    /** @brief Vecteur Jeveux '.CELK' */
    JeveuxVectorChar24 _reference;
    /** @brief Vecteur Jeveux '.CELV' */
    JeveuxVector< ValueType > _values;
    /** @brief Finite element description */
    FiniteElementDescriptorPtr _dofDescription;
    /** @brief Object for dynamic fields (as VARI_ELGA) */
    std::shared_ptr< SimpleFieldOnCells< ASTERINTEGER > > _DCEL;

  public:
    using FieldOnCellsPtr = std::shared_ptr< FieldOnCells >;

    /**
     * @brief Constructor with given name
     * @param name Jeveux name of the field
     */
    FieldOnCells( const std::string name, const FiniteElementDescriptorPtr fed = nullptr )
        : DataField( name, "CHAM_ELEM" ),
          _descriptor( JeveuxVectorLong( getName() + ".CELD" ) ),
          _reference( JeveuxVectorChar24( getName() + ".CELK" ) ),
          _values( JeveuxVector< ValueType >( getName() + ".CELV" ) ),
          _dofDescription( fed ),
          _DCEL( nullptr ) {};

    /** @brief Constructor with automatic name */
    FieldOnCells() : FieldOnCells( ResultNaming::getNewResultName() ) {};

    /** @brief Constructor with automatic name and FE Descriptor */
    FieldOnCells( const FiniteElementDescriptorPtr FEDesc )
        : FieldOnCells( ResultNaming::getNewResultName(), FEDesc ) {};

    /** @brief Constructor with automatic name and model */
    FieldOnCells( const ModelPtr model ) : FieldOnCells( model->getFiniteElementDescriptor() ) {};

    /** @brief Copy constructor */
    FieldOnCells( const std::string &name, const FieldOnCells &toCopy ) : FieldOnCells( name ) {
        *this = toCopy;
    }

    /** @brief constructor */
    FieldOnCells( const FiniteElementDescriptorPtr FEDesc, const std::string &loc,
                  const std::string &quantity )
        : FieldOnCells( FEDesc ) {
        std::string option;
        std::string nompar;

        if ( loc == "ELGA" ) {
            option = "TOU_INI_ELGA";
        } else if ( loc == "ELNO" ) {
            option = "TOU_INI_ELNO";
        } else if ( loc == "ELEM" ) {
            option = "TOU_INI_ELEM";
        } else {
            option = loc;
        };

        if ( quantity[0] != 'P' ) {
            nompar = "P" + quantity;
        } else {
            nompar = quantity;
        };

        ASTERINTEGER iret = 0;

        std::string dcel = " ";

        CALLO_ALCHML( getDescription()->getName(), option, nompar,
                      JeveuxMemoryTypesNames[Permanent], getName(), &iret, dcel );

        if ( iret != 0 ) {
            throw std::runtime_error( "Error creating field " + quantity + " at " + loc );
        }

        updateValuePointers();
    }

    FieldOnCells( const ModelPtr model, const std::string &loc, const std::string &quantity )
        : FieldOnCells( model->getFiniteElementDescriptor(), loc, quantity ) {};

    /** @brief Move constructor */
    FieldOnCells( FieldOnCells &&other ) : DataField( std::move( other ) ) {
        _descriptor = other._descriptor;
        _reference = other._reference;
        _values = other._values;
        _title = other._title;
        setDescription( other._dofDescription );
        updateValuePointers();
    }

    /**
     * @brief Copy constructor
     */
    FieldOnCells( const FieldOnCells &toCopy )
        : FieldOnCells( DataStructureNaming::getNewName(), toCopy ) {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    FieldOnCells( const py::tuple &tup )
        : FieldOnCells( tup[0].cast< std::string >(),
                        tup[1].cast< FiniteElementDescriptorPtr >() ) {
        this->updateValuePointers();
    };

    py::tuple _getState() const { return py::make_tuple( getName(), getDescription() ); };

    /**
     * @brief Wrap of copy constructor
     */
    FieldOnCells copy() { return *this; }

    /**
     * @brief
     * @return
     */
    void deallocate() {
        _descriptor->deallocate();
        _reference->deallocate();
        _values->deallocate();
        _dofDescription = nullptr;
    };

    /**
     * @brief Surcharge de l'operateur =
     */
    FieldOnCells &operator=( const FieldOnCells &toCopy ) {
        *_descriptor = *toCopy._descriptor;
        *_reference = *toCopy._reference;
        *_values = *toCopy._values;
        *_title = *toCopy._title;
        setDescription( toCopy._dofDescription );
        updateValuePointers();
        return *this;
    };

    /**
     * @brief Check if fields are OK for +, +=, ...
     * @return true if compatible
     */
    ASTERBOOL isSimilarTo( const FieldOnCells< ValueType > &tmp2 ) const {
        bool similar = ( _descriptor->size() == tmp2._descriptor->size() );
        similar = ( similar && ( _reference->size() == tmp2._reference->size() ) );
        similar = ( similar && ( _values->size() == tmp2._values->size() ) );
        return similar;
    }

    /** @brief Get the mesh */
    BaseMeshPtr getMesh() const {
        if ( _dofDescription ) {
            return _dofDescription->getMesh();
        }
        return nullptr;
    };

    /** @brief Get datastructure for dynamic fields (as VARI_ELGA) */
    auto getExtentedInformations() const { return _DCEL; };

    /** @brief Set datastructure for dynamic fields (as VARI_ELGA) */
    void
    setExtentedInformations( const std::shared_ptr< SimpleFieldOnCells< ASTERINTEGER > > dcel ) {
        if ( dcel ) {
            _DCEL = dcel;
        }
    };

    /**
     * @brief Set the description of finite elements
     * @param curDesc object FiniteElementDescriptorPtr
     */
    void setDescription( const FiniteElementDescriptorPtr &curDesc ) {
        if ( !curDesc && _dofDescription ) {
            AS_ABORT( "FiniteElementDescriptor is empty" );
        }

        if ( _dofDescription && curDesc && _dofDescription != curDesc ) {
            std::string msg =
                "FiniteElementDescriptor inconsistents: " + _dofDescription->getName() + " vs " +
                curDesc->getName();
            AS_ABORT( msg );
        }

        _dofDescription = curDesc;
    };

    /** @brief Get the description of finite elements */
    FiniteElementDescriptorPtr getDescription() const { return _dofDescription; };

    bool exists() const { return _descriptor.exists() && _values.exists(); };

    /**  @brief Update field and build FiniteElementDescriptor if necessary */
    ASTERBOOL build( std::vector< FiniteElementDescriptorPtr > FEDs =
                         std::vector< FiniteElementDescriptorPtr >() ) {
        if ( !_dofDescription ) {
            CALL_JEMARQ();
            _reference->updateValuePointer();
            const std::string ligrel = strip( ( *_reference )[0].toString() );
            CALL_JEDEMA();

            if ( ligrel.substr( 0, 8 ) == getName().substr( 0, 8 ) ) {
                setDescription( std::make_shared< FiniteElementDescriptor >( ligrel, getMesh() ) );
            } else {
                for ( auto &fed : FEDs ) {
                    if ( fed && strip( fed->getName() ) == ligrel ) {
                        setDescription( fed );
                        break;
                    }
                }
            }
        }

        return true;
    };

    /**
     * @brief Mise a jour des pointeurs Jeveux
     * @return renvoie true si la mise a jour s'est bien deroulee, false sinon
     */
    void updateValuePointers() const {
        _descriptor->updateValuePointer();
        _reference->updateValuePointer();
        _values->updateValuePointer();
    };

    /**
     * @brief Transformer les valeurs de _values en appliquant
     *        la fonction "func" à chaque valeur
     * @return renvoie un nouvel objet de FieldOnCells
     *         avec les valeurs transformées
     */
    FieldOnCells< ValueType > transform( py::object &func ) const;

    // OVERLOADING C++ OPERATORS

    /**
     * @brief Unary Minus overloading
     * @return Updated field
     */
    FieldOnCells< ValueType > operator-() const {
        FieldOnCells< ValueType > tmp( *this );
        ( *tmp._values ) *= ValueType( -1 );
        return tmp;
    };

    /**
     * @brief Shorthand + operator assignement
     * @return Updated field
     */
    FieldOnCells< ValueType > &operator+=( const FieldOnCells< ValueType > &rhs ) {
        if ( !this->isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        ( *_values ) += ( *rhs._values );
        return *this;
    };

    /**
     * @brief Shorthand - operator assignement
     * @return Updated field
     */
    FieldOnCells< ValueType > &operator-=( const FieldOnCells< ValueType > &rhs ) {
        if ( !this->isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        ( *_values ) -= ( *rhs._values );
        return *this;
    };

    /**
     * @brief Subscripting overloading
     * @param i subscript
     * @return value at position i
     */
    ValueType &operator[]( int i ) { return _values->operator[]( i ); };

    const ValueType &operator[]( int i ) const { return _values->operator[]( i ); };

    /**
     * @brief Plus overloading
     * @return New field
     */
    friend FieldOnCells< ValueType > operator+( FieldOnCells< ValueType > lhs,
                                                const FieldOnCells< ValueType > &rhs ) {
        if ( !lhs.isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        lhs += rhs;
        return lhs;
    };

    /**
     * @brief Minus overloading
     * @return New field
     */
    friend FieldOnCells< ValueType > operator-( FieldOnCells< ValueType > lhs,
                                                const FieldOnCells< ValueType > &rhs ) {
        if ( !lhs.isSimilarTo( rhs ) )
            raiseAsterError( "Fields have incompatible shapes" );
        lhs -= rhs;
        return lhs;
    };

    /**
     * @brief Multiply by a scalar on right overloading
     * @return New field
     */

    friend FieldOnCells< ValueType > operator*( FieldOnCells< ValueType > lhs,
                                                const ASTERDOUBLE &scal ) {
        ( *lhs._values ) *= scal;
        return lhs;
    };

    /**
     * @brief Multiply by a scalar on left overloading
     * @return New field
     */

    friend FieldOnCells< ValueType > operator*( const ASTERDOUBLE &scal,
                                                const FieldOnCells< ValueType > &rhs ) {
        return rhs * scal;
    };

    // some getters

    /**
     * @brief Get values of the field
     */
    const JeveuxVector< ValueType > &getValues() const {
        _values->updateValuePointer();
        return _values;
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
     * @brief Get descriptor of the field
     */
    JeveuxVectorLong getDescriptor() const { return _descriptor; };

    /**
     * @brief Set the Values object
     *
     * @param value Value to affect
     */
    void setValues( const ValueType &value ) {
        _values->updateValuePointer();
        _values->assign( value );
    };

    void setValues( const std::vector< ValueType > &values ) {
        if ( values.size() != size() ) {
            raiseAsterError( "Incompatible size in setValues, expected: " +
                             std::to_string( size() ) );
        }

        *_values = values;
    };

    std::string getPhysicalQuantity() const {
        _descriptor->updateValuePointer();
        auto gd = ( *_descriptor )[0];
        return PhysicalQuantityManager::getPhysicalQuantityName( gd );
    }

    std::string getLocalization() const {
        _reference->updateValuePointer();

        return strip( ( *_reference )[2] );
    }

    VectorString getComponents() const {
        JeveuxVectorChar8 cmp( "&CMP" );
        ASTERINTEGER ncmp;
        VectorString cmps;

        if ( getPhysicalQuantity() == "VARI_R" ) {
            ncmp = ( *_descriptor )[3];
            cmps.reserve( ncmp );
            for ( auto icmp = 0; icmp < ncmp; icmp++ ) {
                cmps.push_back( "V" + std::to_string( icmp + 1 ) );
            }
        } else {
            CALL_UTNCMP( getName().c_str(), &ncmp, cmp->getName().c_str() );

            cmp->updateValuePointer();

            cmps.reserve( cmp->size() );

            for ( auto &cm : cmp ) {
                cmps.push_back( strip( cm.toString() ) );
            }
        }

        return cmps;
    }

    ASTERINTEGER getNumberOfComponents() const {
        if ( getPhysicalQuantity() == "VARI_R" ) {
            return ( *_descriptor )[3];
        } else {
            return getComponents().size();
        }
    }

    ASTERINTEGER getNumberOfGroupOfElements() const {
        _descriptor->updateValuePointer();

        return ( *_descriptor )[1];
    }

    ASTERINTEGER getNumberOfElements( const ASTERINTEGER &iGrel ) {
#ifdef ASTER_DEBUG_CXX
        if ( iGrel < 0 || iGrel >= getNumberOfGroupOfElements() ) {
            AS_ABORT( "Out of bounds" );
        }
#endif
        _descriptor->updateValuePointer();
        return ( *_descriptor )[( *_descriptor )[4 + iGrel] - 1 + 1];
    }

    ASTERINTEGER getSizeOfFieldOfElement( const ASTERINTEGER &iGrel ) {
#ifdef ASTER_DEBUG_CXX
        if ( iGrel < 0 || iGrel >= getNumberOfGroupOfElements() ) {
            AS_ABORT( "Out of bounds" );
        }
#endif
        _descriptor->updateValuePointer();
        return ( *_descriptor )[( *_descriptor )[4 + iGrel] - 1 + 3];
    }

    ASTERINTEGER getLocatedComponentIdOfGroupOfElements( const ASTERINTEGER &iGrel ) {
#ifdef ASTER_DEBUG_CXX
        if ( iGrel < 0 || iGrel >= getNumberOfGroupOfElements() ) {
            AS_ABORT( "Out of bounds" );
        }
#endif
        _descriptor->updateValuePointer();
        return ( *_descriptor )[( *_descriptor )[4 + iGrel] - 1 + 2];
    }

    ASTERINTEGER getShifting( const ASTERINTEGER &iGrel, const ASTERINTEGER &iElem ) {
#ifdef ASTER_DEBUG_CXX
        if ( iGrel < 0 || iGrel >= getNumberOfGroupOfElements() ) {
            AS_ABORT( "Out of bounds" );
        }
        if ( iElem < 0 || iElem >= getNumberOfElements( iGrel ) ) {
            AS_ABORT( "Out of bounds" );
        }
#endif
        _descriptor->updateValuePointer();
        return ( *_descriptor )[( *_descriptor )[4 + iGrel] - 1 + 4 + 4 * iElem + 4] - 1;
    }

    std::string getLocatedComponent() const {
        const auto nbGrel = getNumberOfGroupOfElements();
        _descriptor->updateValuePointer();
        const std::string cata = "&CATA.TE.NOMMOLOC";
        JeveuxChar24 objName, charName;

        std::string modeName;
        for ( auto igr = 0; igr < nbGrel; igr++ ) {
            auto mode = getLocatedComponentIdOfGroupOfElements( igr );
            if ( mode > 0 ) {
                CALLO_JEXNUM( objName, cata, &mode );
                CALLO_JENUNO( objName, charName );
                auto modeN = strip( charName.toString().substr( 14, 10 ) );
                if ( modeName.empty() ) {
                    modeName = modeN;
                } else {
                    if ( modeName != modeN ) {
                        AS_ABORT( "Multiple names." );
                    }
                }
            }
        }

        return modeName;
    }

    /**  @brief Size of the field */
    ASTERINTEGER size() const { return _values->size(); }

    // norm and dot methods

    /**
     * @brief Comput norm
     * @param normType Type of norm ("NORM_1","NORM_2","NORM_INFINITY")
     */
    ASTERDOUBLE norm( const std::string normType ) const;

    /**
     * @brief Dot product
     * @param tmp object FieldOnCellsPtr
     */
    ASTERDOUBLE dot( const FieldOnCellsPtr &tmp ) const;

    bool printMedFile( const std::filesystem::path &fileName, bool local = true,
                       std::string version = std::string() ) const;

    FieldOnCellsPtr asLocalization( const std::string &loc ) const {
        if ( loc == getLocalization() ) {
            return std::make_shared< FieldOnCells< ValueType > >( *this );
        }

        auto cham_model = std::make_shared< FieldOnCells< ValueType > >( _dofDescription, loc,
                                                                         getPhysicalQuantity() );

        auto cham_elem = std::make_shared< FieldOnCells< ValueType > >();

        std::string base = "G";
        std::string prol = "OUI", model = " ";

        CALLO_CHPCHD( getName(), loc, cham_model->getName(), prol, base, cham_elem->getName(),
                      model );

        cham_elem->build( { _dofDescription } );

        return cham_elem;
    }

    /**
     * @brief Check field for internal state variables
     * @param tmp object FieldOnCellsPtr
     */
    void checkInternalStateVariables( const ConstantFieldOnCellsChar16Ptr prevBehaviour,
                                      const ConstantFieldOnCellsChar16Ptr currBehaviour,
                                      const FiniteElementDescriptorPtr newFEDesc );

    /**
     * @brief Compare shape of field and project if it required
     * @param tmp object FieldOnCellsPtr
     */
    ASTERINTEGER compareShape( const FieldOnCellsRealPtr fieldModel, const bool projectOnLigrel,
                               const std::string paraName );

    friend class FieldBuilder;
};

template < class ValueType >
bool FieldOnCells< ValueType >::printMedFile( const std::filesystem::path &fileName, bool local,
                                              std::string version ) const {
    const auto rank = getMPIRank();
    LogicalUnitFile a;
    ASTERINTEGER retour = -1;
    // In case that the print file (single and absolute path) is unique between processors,
    // it must only be created on proc 0.
    if ( getMesh()->isParallel() || ( !getMesh()->isParallel() && rank == 0 ) ) {
        if ( rank == 0 )
            a.openFile( fileName, Binary, New );
#ifdef ASTER_HAVE_MPI
        AsterMPI::barrier();
#endif /* ASTER_HAVE_MPI */
        if ( rank != 0 )
            a.openFile( fileName, Binary, Old );
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

/** @typedef FieldOnCellsReal */
using FieldOnCellsReal = FieldOnCells< ASTERDOUBLE >;
using FieldOnCellsRealPtr = std::shared_ptr< FieldOnCellsReal >;

/** @typedef FieldOnCellsLong */
using FieldOnCellsLong = FieldOnCells< ASTERINTEGER >;
using FieldOnCellsLongPtr = std::shared_ptr< FieldOnCellsLong >;

/** @typedef FieldOnCellsComplex */
using FieldOnCellsComplex = FieldOnCells< ASTERCOMPLEX >;
using FieldOnCellsComplexPtr = std::shared_ptr< FieldOnCellsComplex >;

/** @typedef FieldOnCellsChar8 */
using FieldOnCellsChar8 = FieldOnCells< JeveuxChar8 >;
using FieldOnCellsChar8Ptr = std::shared_ptr< FieldOnCellsChar8 >;

#endif /* FIELDONCELLS_H_ */
