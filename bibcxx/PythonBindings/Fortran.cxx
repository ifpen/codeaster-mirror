/**
 * @file Fortran.h
 * @brief Definition of interface functions between C++ and Fortran
 * @author Mathieu Courtois
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

#include "PythonBindings/Fortran.h"

#include "aster_fort_ds.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"
#include "aster_mpi.h"
#include "aster_pybind.h"
#include "aster_utils.h"
#include "shared_vars.h"

void jeveux_init( int fcomm ) {
    ASTERINTEGER dbg = 0;
    aster_mpi_init( (MPI_Fint)fcomm );
    CALL_IBMAIN();

    // now Jeveux is available
    register_sh_jeveux_status( 1 );
}

void jeveux_finalize( const ASTERINTEGER options ) {
    if ( get_sh_jeveux_status() != 1 ) {
        return;
    }
    ASTERINTEGER fopt = options;
    CALL_OP9999( &fopt );
    register_sh_jeveux_status( 0 );
}

void call_oper( py::object &syntax, int jxveri ) {
    ASTERINTEGER jxvrf = jxveri;

    // Add the new syntax object on the stack
    register_sh_etape( append_etape( syntax.ptr() ) );

    try {
        CALL_EXPASS( &jxvrf );

    } catch ( ... ) {
        // unstack the syntax object
        register_sh_etape( pop_etape() );
        throw;
    }
    // unstack the syntax object
    register_sh_etape( pop_etape() );
}

void call_oper_init() {
    // to execute before/after steps around operator
    ASTERINTEGER op = 8888;
    CALL_EXECOP( &op );
}

void call_cmd_ctxt_enter() { CALL_CMD_CTXT_ENTER(); }

void call_cmd_ctxt_exit() { CALL_CMD_CTXT_EXIT(); }

void call_affich( const std::string &code, const std::string &text ) { CALL_AFFICH( code, text ); }

void call_print( const std::string &text ) { call_affich( "MESSAGE", text ); }

void jeveux_delete( const std::string prefix ) {
    std::string type = " ";
    CALLO_DETRSD( type, prefix );
}

bool deleteTemporaryObjects( void ) {
    CALL_DELETE_TEMPORARY_OBJECTS();
    return true;
}

bool deleteCachedObjects( void ) {
    CALL_DELETE_CACHED_OBJECTS();
    return true;
}

std::string onFatalError( const std::string value ) {
    ASTERINTEGER lng = 16;
    char *tmp = MakeBlankFStr( lng );
    std::string blank( " " );

    if ( value == "ABORT" || value == "EXCEPTION" ||
         value == "EXCEPTION+VALID" | value == "INIT" ) {
        CALL_ONERRF( (char *)value.c_str(), tmp, &lng );
    } else {
        CALL_ONERRF( (char *)blank.c_str(), tmp, &lng );
    }
    std::string state( tmp, lng );
    FreeStr( tmp );
    return state;
}

void call_matfpe( const int value ) {
    ASTERINTEGER enable;

    if ( value == -1 || value == 1 ) {
        enable = (ASTERINTEGER)value;
        CALL_MATFPE( &enable );
    } else {
        AS_ABORT( "Invalid value: expecting -1 or 1" );
    }
}

int asmpi_get() {
    const std::string action( "GET" );
    MPI_Fint comm;
    CALLO_ASMPI_COMM( action, &comm );
    return comm;
}

void asmpi_set( const int comm ) {
    MPI_Fint fcomm;
    const std::string action( "SET" );
    fcomm = MPI_Fint( comm );
    CALLO_ASMPI_COMM( action, &fcomm );
}

void asmpi_free( const int comm ) {
    MPI_Fint fcomm;
    const std::string action( "FREE" );
    fcomm = MPI_Fint( comm );
    CALLO_ASMPI_COMM( action, &fcomm );
}

VectorInt asmpi_info( const int comm ) {
    MPI_Fint fcomm;
    MPI_Fint size, rank;
    fcomm = MPI_Fint( comm );
    CALL_ASMPI_INFO( &fcomm, &rank, &size );
    return VectorInt { rank, size };
}

int asmpi_split( const int parent, int color, std::string name ) {
    MPI_Fint fparent, newcomm;
    int key( 0 );
    fparent = MPI_Fint( parent );
    CALL_ASMPI_SPLIT_COMM( &fparent, &color, &key, name, &newcomm );
    return newcomm;
}

void set_option( const std::string &option, const ASTERINTEGER value ) {
    py::object method = py::cast< py::object >( GetJdcAttr( (char *)"set_option" ) );
    py::object none = method( option, value );
}
