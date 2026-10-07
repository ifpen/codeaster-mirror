/**
 * @file ListOfTables.cxx
 * @brief Implementation de ListOfTables
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

#include "DataFields/ListOfTables.h"

#include "aster_fort_jeveux.h"

#include "Utilities/Tools.h"

ListOfTables::ListOfTables( const std::string &name ) : _name( name ) {
    _name.resize( 19, ' ' );
    _dsId = JeveuxVectorChar16( _name + ".LTNT" );
    _dsName = JeveuxVectorChar24( _name + ".LTNS" );
}

bool ListOfTables::update_tables() {
    CALL_JEMARQ();

    if ( !_dsId.exists() ) {
        CALL_JEDEMA();
        return true;
    } else {
        _dsId->updateValuePointer();
    }

    if ( !_dsName.exists() ) {
        CALL_JEDEMA();
        return true;
    } else {
        _dsName->updateValuePointer();
    }

    const int size = _dsId->size();
    for ( int i = 0; i < size; i++ ) {
        const auto id = strip( ( *_dsId )[i].toString() );
        const auto name = strip( ( *_dsName )[i].toString() );
        if ( id == "" )
            continue;
        if ( _mapTables[id] == nullptr ) {
            _mapTables[id] = std::make_shared< Table >( name );
            _mapTables[id]->build();
        }
    }
    // todo: remove previously registered and not available anymore

    CALL_JEDEMA();
    return true;
}

TablePtr ListOfTables::getTable( const std::string id ) {
    this->update_tables();
    const auto id_ = strip( id );
    const auto curIter = _mapTables.find( id_ );
    if ( curIter == _mapTables.end() )
        return TablePtr( nullptr );
    return curIter->second;
}

void ListOfTables::setTable( const std::string id, const TablePtr table ) {
    this->update_tables();
    const auto id_ = strip( id );
    const auto curIter = _mapTables.find( id_ );
    if ( curIter != _mapTables.end() ) {
#ifdef ASTER_DEBUG_CXX
        std::cout << "Override existing table" << std::endl << std::flush;
#endif
        _mapTables[id_] = table;
        _dsId->updateValuePointer();
        const int size = _dsId->size();
        for ( int i = 0; i < size; i++ ) {
            if ( id_ == strip( ( *_dsId )[i].toString() ) )
                ( *_dsName )[i] = table->getName();
        }
    } else {
        const int i = _dsId->size();
        _dsId->resize( i + 1 );
        _dsName->resize( i + 1 );
        ( *_dsId )[i] = id_;
        ( *_dsName )[i] = table->getName();
    }
    _mapTables[id_] = table;
}
