/**
 * @file PhysicalQuantityManager.cxx
 * @brief Implementation de PhysicalQuantityManager
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

#include "Modeling/PhysicalQuantityManager.h"

#include "aster_fort_calcul.h"
#include "aster_fort_utils.h"

NamesMapChar8 PhysicalQuantityManager::_nameOfPhysicalQuantity( "&CATA.GD.NOMGD" );
JeveuxCollectionChar8 PhysicalQuantityManager::_nameOfCmp( "&CATA.GD.NOMCMP" );

bool PhysicalQuantityManager::hasQuantityOfName( const std::string name ) {
    return _nameOfPhysicalQuantity->getIndexFromString( name ) != -1;
};

std::string PhysicalQuantityManager::getPhysicalQuantityName( const ASTERINTEGER quantityNumber ) {
    if ( quantityNumber <= 0 || quantityNumber > _nameOfPhysicalQuantity->size() )
        throw std::runtime_error( "Not a known physical quantity" );
    return strip( _nameOfPhysicalQuantity->getStringFromIndex( quantityNumber ) );
};

ASTERINTEGER PhysicalQuantityManager::getPhysicalQuantityNumber( const std::string name ) {
    return _nameOfPhysicalQuantity->getIndexFromString( name );
};

const VectorString PhysicalQuantityManager::getAllPhysicalQuantityNames() {
    VectorString result;
    auto sz = _nameOfPhysicalQuantity->size();
    result.reserve( sz );
    for ( int i = 0; i < sz; i++ )
        result.push_back( strip( _nameOfPhysicalQuantity->getStringFromIndex( i + 1 ) ) );
    return result;
};

ASTERINTEGER
PhysicalQuantityManager::getNumberOfEncodedInteger( const ASTERINTEGER quantityNumber ) {
    ASTERINTEGER toReturn = 0;
    toReturn = CALL_NBEC( &quantityNumber );
    return toReturn;
};

ASTERINTEGER PhysicalQuantityManager::getNumberOfComponents( const ASTERINTEGER quantityNumber ) {
    return _nameOfCmp->getObjectSize( quantityNumber );
};

const VectorString PhysicalQuantityManager::getComponentNames( const ASTERINTEGER quantityNumber ) {
    _nameOfCmp->build();
    auto vc = ( *_nameOfCmp )[quantityNumber]->toVector();
    VectorString result;
    result.reserve( vc.size() );
    for ( int i = 0; i < vc.size(); i++ )
        result.push_back( vc[i].toString() );
    return result;
};

ASTERINTEGER PhysicalQuantityManager::getPhysicalQuantityNumber( const std::string option,
                                                                 const std::string parameter ) {

    ASTERINTEGER gd = -1;
    std::string oprionStr = ljust( strip( option ), 16, ' ' );
    std::string parameterStr = ljust( strip( parameter ), 8, ' ' );

    CALL_GRANDEURID( oprionStr, parameterStr, &gd );

    return gd;
}

std::string PhysicalQuantityManager::getPhysicalQuantityName( const std::string option,
                                                              const std::string parameter ) {

    return getPhysicalQuantityName( getPhysicalQuantityNumber( option, parameter ) );
}
