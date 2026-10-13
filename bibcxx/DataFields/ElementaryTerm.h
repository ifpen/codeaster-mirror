#ifndef ELEMENTARYTERM_H_
#define ELEMENTARYTERM_H_

/**
 * @file ElementaryTerm.h
 * @brief Fichier entete de la classe ElementaryTerm
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

#include "aster_fort_calcul.h"
#include "aster_fort_utils.h"

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/PhysicalQuantityManager.h"
#include "Supervis/Exceptions.h"

/**
 * @class ElementaryTerm
 * @brief Class which describe a RESUELEM (which is part of MATR_ELEM and VECT_ELEM )
 */
template < typename ValueType >
class ElementaryTerm : public DataField {
  private:
    /** @brief Objet Jeveux '.NOLI' */
    JeveuxVectorChar24 _noli;
    /** @brief Objet Jeveux '.DESC' */
    JeveuxVectorLong _descriptor;
    /** @brief Objet Jeveux '.RESL' */
    JeveuxCollection< ValueType > _resl;
    /** @brief Ligel */
    FiniteElementDescriptorPtr _FEDesc;

  public:
    /**
     * @brief Constructor only with predefined name
     * @param name predefined name
     */
    ElementaryTerm( const std::string name )
        : DataField( name, "RESUELEM" ),
          _noli( JeveuxVectorChar24( getName() + ".NOLI" ) ),
          _descriptor( JeveuxVectorLong( getName() + ".DESC" ) ),
          _resl( JeveuxCollection< ValueType >( getName() + ".RESL" ) ) {};

    ElementaryTerm() : ElementaryTerm( DataStructureNaming::getNewName() ) {};

    void allocate( const FiniteElementDescriptorPtr fEDesc, const std::string &option,
                   const ASTERINTEGER &physicalQuantityId,
                   const ASTERINTEGER &locatedComponentId ) {
        _FEDesc = fEDesc;
        _noli->allocate( 4 );
        ( *_noli )[0] = fEDesc->getName();
        ( *_noli )[1] = option;
        ( *_noli )[2] = "MPI_COMPLET";
        const auto &nbGrel = fEDesc->getListOfGroupsOfElements()->size();
        _descriptor->allocate( 2 + nbGrel );
        _descriptor->setInformationParameter( "RESL" );
        ( *_descriptor )[0] = physicalQuantityId;
        ( *_descriptor )[1] = nbGrel;
        for ( int i = 0; i < nbGrel; ++i )
            ( *_descriptor )[i + 2] = locatedComponentId;
        if ( nbGrel > 0 ) {
            _resl->allocate( nbGrel );
        }
    };

    void allocate( const FiniteElementDescriptorPtr fEDesc, const std::string &option,
                   const std::string &physicalQuantity, const ASTERINTEGER &locatedComponentId ) {

        const auto physicalQuantityId =
            PhysicalQuantityManager::getPhysicalQuantityId( physicalQuantity );

        allocate( fEDesc, option, physicalQuantityId, locatedComponentId );
    };

    void allocate( const FiniteElementDescriptorPtr fEDesc, const std::string &option,
                   const std::string &parameter ) {

        const auto physicalQuantityId =
            PhysicalQuantityManager::getPhysicalQuantityId( option, parameter );

        allocate( fEDesc, option, physicalQuantityId, -1 );

        const auto &nbGrel = fEDesc->getListOfGroupsOfElements()->size();

        ASTERINTEGER modloc = -1;
        std::string optionStr = ljust( strip( option ), 16, ' ' );
        std::string parameterStr = ljust( strip( parameter ), 8, ' ' );

        for ( int i = 0; i < nbGrel; ++i ) {
            const auto typeFE = fEDesc->getFiniteElementTypeIdOfGroupOfElements( i );

            CALL_GET_LOCAL_MODE( optionStr, parameterStr, &typeFE, &modloc );
            ( *_descriptor )[i + 2] = modloc;
        }
    };

    void setFiniteElementDescriptor( const FiniteElementDescriptorPtr FEDesc ) {
        if ( FEDesc ) {
            if ( _FEDesc && _FEDesc != FEDesc ) {
                std::string mess =
                    "Incompatible FED: " + _FEDesc->getName() + " vs " + FEDesc->getName();
                AS_ABORT( mess );
            } else {
                _noli->updateValuePointer();
                auto FEDname = strip( ( *_noli )[0].toString() );
                if ( FEDname != strip( FEDesc->getName() ) ) {
                    std::string mess = "Incompatible FED: " + FEDname + " vs " + FEDesc->getName();
                    raiseAsterError( mess );
                }

                _FEDesc = FEDesc;
            }
        }
    };

    const JeveuxCollection< ValueType > &getValues() const { return _resl; };

    JeveuxCollection< ValueType > &getValues() { return _resl; };

    FiniteElementDescriptorPtr getFiniteElementDescriptor() const { return _FEDesc; };

    bool exists() const { return _noli.exists() && _descriptor.exists() && _resl.exists(); };

    std::string getOption() {
        _noli->updateValuePointer();
        return strip( ( *_noli )[1].toString() );
    };

    BaseMeshPtr getMesh() const {
        if ( _FEDesc ) {
            return _FEDesc->getMesh();
        }
        return nullptr;
    }

    std::string getPhysicalQuantity() {
        const std::string typeco( "RESUELEM" );
        ASTERINTEGER repi = 0, ier = 0;
        JeveuxChar32 repk( " " );
        const std::string arret( "F" );
        const std::string questi( "NOM_GD" );
        CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
        return strip( repk.toString() );
    }

    ASTERINTEGER getPhysicalQuantityId() {
        _descriptor->updateValuePointer();
        return ( *_descriptor )[0];
    }

    ASTERINTEGER getNumberOfGroupOfCells() const {
        _descriptor->updateValuePointer();

        return ( *_descriptor )[1];
    }

    /**
     * @brief Return MODE_LOCAL
     */
    std::string getLocatedComponent() const {
        const auto nbGrel = getNumberOfGroupOfCells();
        _descriptor->updateValuePointer();
        const std::string cata = "&CATA.TE.NOMMOLOC";
        JeveuxChar24 objName, charName;

        std::string modeName;
        for ( auto igr = 0; igr < nbGrel; igr++ ) {
            auto mode = ( *_descriptor )[2 + igr];
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

    /**
     * @brief Return MODE_LOCAL
     */
    ASTERINTEGER getLocatedComponentId() const {
        const auto nbGrel = getNumberOfGroupOfCells();
        _descriptor->updateValuePointer();

        ASTERINTEGER modeId = -1;
        for ( auto igr = 0; igr < nbGrel; igr++ ) {
            auto mode = ( *_descriptor )[2 + igr];
            if ( mode > 0 ) {
                if ( modeId == -1 ) {
                    modeId = mode;
                } else {
                    if ( modeId != mode ) {
                        AS_ABORT( "Multiple names." );
                    }
                }
            }
        }

        return modeId;
    }

    ASTERINTEGER getLocatedComponentIdOfGroupOfElements( const ASTERINTEGER grel_idx ) const {
        _descriptor->updateValuePointer();
        return ( *_descriptor )[2 + grel_idx];
    }

    ASTERINTEGER getNumberOfComponentsOfGroupOfElements( const ASTERINTEGER grel_idx ) const {

        const auto mode_id = getLocatedComponentIdOfGroupOfElements( grel_idx );
        ASTERINTEGER nbCmp = -1;
        CALL_DIGDEL_WRAP( &mode_id, &nbCmp );
        return nbCmp;
    }

    bool isMPIFull() {
        AS_ASSERT( _noli.exists() );
        _noli->updateValuePointer();
        return strip( ( *_noli )[2].toString() ) == "MPI_COMPLET";
    }

    bool isEmpty() { return !( _noli.exists() && _descriptor.exists() ); }

    bool build() { return _resl->build( true ); };

    /**
     * @brief TimesEqual overloading
     */
    ElementaryTerm< ValueType > &operator*=( const ValueType &scal ) {

        ( *_resl ) *= scal;

        return *this;
    };
};

/** @typedef ElementaryTermRealPtr */
using ElementaryTermReal = ElementaryTerm< ASTERDOUBLE >;
using ElementaryTermRealPtr = std::shared_ptr< ElementaryTermReal >;

/** @typedef ElementaryTermComplexPtr */
using ElementaryTermComplex = ElementaryTerm< ASTERCOMPLEX >;
using ElementaryTermComplexPtr = std::shared_ptr< ElementaryTermComplex >;

#endif /* ELEMENTARYTERM_H_ */
