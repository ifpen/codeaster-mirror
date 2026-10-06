#ifndef ELEMENTARYMATRIX_H_
#define ELEMENTARYMATRIX_H_

/**
 * @file ElementaryMatrix.h
 * @brief Definition of elementary matrices
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

#include "aster_fort_calcul.h"

#include "DataFields/ElementaryTerm.h"
#include "LinearAlgebra/BaseElementaryMatrix.h"

/**
 * @class ElementaryMatrix
 * @brief Class for sd_matr_elem template
 */
template < typename ValueType, PhysicalQuantityEnum PhysicalQuantity >
class ElementaryMatrix : public BaseElementaryMatrix {
  private:
    /** @brief Vectors of RESUELEM */
    std::vector< std::shared_ptr< ElementaryTerm< ValueType > > > _elemTerm;

    void checkMesh( const BaseMeshPtr mesh ) const {
        const auto _mesh = getMesh();
        if ( _mesh && mesh ) {
            if ( _mesh != mesh ) {
                AsterErrorCpp( "Meshes are incompatible." );
            }
        }
    };

  public:
    /** @typedef ElementaryMatrixPtr */
    typedef std::shared_ptr< ElementaryMatrix< ValueType, PhysicalQuantity > > ElementaryMatrixPtr;

    /** @brief Constructor with a name */
    ElementaryMatrix( const std::string name )
        : BaseElementaryMatrix(
              name, "MATR_ELEM_" + std::string( PhysicalQuantityNames[PhysicalQuantity] ) +
                        ( typeid( ValueType ) == typeid( ASTERDOUBLE ) ? "_R" : "_C" ) ) {
        _elemTerm.clear();
    };

    /** @brief Constructor with automatic name */
    ElementaryMatrix() : ElementaryMatrix( ResultNaming::getNewResultName() ) {};

    ElementaryMatrix( const ModelPtr model ) : ElementaryMatrix() { this->setModel( model ); };

    ElementaryMatrix( const ModelPtr model, const std::string option ) : ElementaryMatrix() {
        this->setModel( model );
        this->prepareCompute( option );
    };

    /** @brief Get the mesh */
    BaseMeshPtr getMesh( void ) const {
        if ( _model )
            return _model->getMesh();

        for ( const auto &elemTerm : _elemTerm ) {
            if ( elemTerm->getMesh() ) {
                return elemTerm->getMesh();
            }
        }

        return nullptr;
    };

    /**
     * @brief Function to update ElementaryTerm
     */
    bool build() {
        if ( _elemComp->hasElementaryTerm() ) {
            CALL_REDETR( getName() );

            SetString elemSave;
            for ( auto &elemTerm : _elemTerm ) {
                elemSave.insert( strip( elemTerm->getName() ) );
            }

            auto elemTermNames = _elemComp->getNameOfElementaryTerms();

            SetString elemKeep;
            for ( auto &elemTerm : elemTermNames ) {
                const std::string name = strip( elemTerm.toString() );
                elemKeep.insert( name );
                if ( name != " " && elemSave.count( name ) == 0 ) {
                    _elemTerm.push_back( std::make_shared< ElementaryTerm< ValueType > >( name ) );
                }
            }

            // clean ElementaryTerm
            std::vector< std::shared_ptr< ElementaryTerm< ValueType > > > elemTermNew;
            elemTermNew.reserve( _elemTerm.size() );
            for ( auto &elemTerm : _elemTerm ) {
                auto name = strip( elemTerm->getName() );
                if ( elemKeep.count( name ) > 0 && elemTerm->exists() ) {
                    elemTerm->build();
                    elemTermNew.push_back( elemTerm );
                }
            }
            _elemTerm = std::move( elemTermNew );

            if ( _elemTerm.size() > 0 && !isMPIFull() && !getMesh()->isParallel() ) {
                std::string type = "MATR_ELEM";
                CALLO_SDMPIC( type, getName() );
            }
        }
        _isBuilt = true;
        return true;
    };

    /**
     * @brief is MPI_COMPLET ?
     */
    bool isMPIFull() {
        for ( auto &elemTerm : _elemTerm ) {
            if ( !elemTerm->isMPIFull() ) {
                return false;
            }
        }

        return true;
    };

    /**
     * @brief Add elementary term
     */
    void addElementaryTerm( const std::shared_ptr< ElementaryTerm< ValueType > > &elemTerm ) {
        this->checkMesh( elemTerm->getMesh() );
        _elemComp->addElementaryTerm( elemTerm->getName() );
        _elemTerm.push_back( elemTerm );
    };

    /**
     * @brief Add elementary term
     */
    void addElementaryTerm(
        const std::vector< std::shared_ptr< ElementaryTerm< ValueType > > > &elemTerm ) {
        for ( auto &elem : elemTerm ) {
            this->addElementaryTerm( elem );
        }
    };

    std::vector< std::shared_ptr< ElementaryTerm< ValueType > > > getElementaryTerms() const {
        return _elemTerm;
    }

    /**
     * @brief Get all FiniteElementDescriptors
     * @return vector of all FiniteElementDescriptors
     */
    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptors() const {
        std::vector< FiniteElementDescriptorPtr > FEDs;
        FEDs.reserve( _elemTerm.size() );

        for ( auto &elemTerm : _elemTerm ) {
            auto FED = elemTerm->getFiniteElementDescriptor();
            if ( FED ) {
                FEDs.push_back( FED );
            }
        }

        return FEDs;
    };

    bool hasElementaryTerms() const { return ( _elemTerm.size() != 0 ); };

    ASTERINTEGER getNumberOfElementaryTerms() const { return _elemTerm.size(); }

    /**
     * @brief TimesEqual overloading
     */
    ElementaryMatrix< ValueType, PhysicalQuantity > &operator*=( const ValueType &scal ) {

        for ( auto &et : _elemTerm ) {
            ( *et ) *= scal;
        }

        return *this;
    };
};

/** @typedef Elementary matrix for displacement-double */
template class ElementaryMatrix< ASTERDOUBLE, Displacement >;
using ElementaryMatrixDisplacementReal = ElementaryMatrix< ASTERDOUBLE, Displacement >;
using ElementaryMatrixDisplacementRealPtr = std::shared_ptr< ElementaryMatrixDisplacementReal >;

/** @typedef Elementary matrix for displacement-complex */
template class ElementaryMatrix< ASTERCOMPLEX, Displacement >;
using ElementaryMatrixDisplacementComplex = ElementaryMatrix< ASTERCOMPLEX, Displacement >;
using ElementaryMatrixDisplacementComplexPtr =
    std::shared_ptr< ElementaryMatrixDisplacementComplex >;

/** @typedef Elementary matrix for temperature-double */
template class ElementaryMatrix< ASTERDOUBLE, Temperature >;
using ElementaryMatrixTemperatureReal = ElementaryMatrix< ASTERDOUBLE, Temperature >;
using ElementaryMatrixTemperatureRealPtr = std::shared_ptr< ElementaryMatrixTemperatureReal >;

/** @typedef Elementary matrix for pressure-complex */
template class ElementaryMatrix< ASTERCOMPLEX, Pressure >;
using ElementaryMatrixPressureComplex = ElementaryMatrix< ASTERCOMPLEX, Pressure >;
using ElementaryMatrixPressureComplexPtr = std::shared_ptr< ElementaryMatrixPressureComplex >;

#endif /* ELEMENTARYMATRIX_H_ */
