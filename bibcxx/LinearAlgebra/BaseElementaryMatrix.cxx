/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2026 - EDF - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#include "LinearAlgebra/BaseElementaryMatrix.h"

BaseElementaryMatrix::BaseElementaryMatrix( const std::string name, const std::string type )
    : DataStructure( name, 19, type ),
      _isBuilt( true ),
      _model( nullptr ),
      _elemComp( std::make_shared< ElementaryCompute >( getName() ) ) {};

/** @brief Constructor with automatic name */
BaseElementaryMatrix::BaseElementaryMatrix( const std::string type )
    : BaseElementaryMatrix( ResultNaming::getNewResultName(), type ) {};

BaseElementaryMatrix::BaseElementaryMatrix( const ModelPtr model ) : BaseElementaryMatrix() {
    this->setModel( model );
};

/** @brief  Prepare compute */
void BaseElementaryMatrix::prepareCompute( const std::string option ) {
    if ( option != "WRAP_FORTRAN" ) {
        _elemComp->createDescriptor( _model, option );
    }
};

bool BaseElementaryMatrix::isSymmetric() const {
    const std::string typeco( "MATR_ELEM" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "F" );
    const std::string questi( "TYPE_MATRICE" );
    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );
    return strip( repk.toString() ) == "SYMETRI";
}

ASTERINTEGER BaseElementaryMatrix::numberOfSuperElement() const {
    const std::string typeco( "MATR_ELEM" );
    ASTERINTEGER repi = 0, ier = 0;
    JeveuxChar32 repk( " " );
    const std::string arret( "C" );
    const std::string questi( "NB_SS_ACTI" );

    CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );

    return repi;
};
