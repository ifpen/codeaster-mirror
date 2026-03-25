#ifndef FORTRAN_H_
#define FORTRAN_H_

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

#include "aster_pybind.h"

void jeveux_init( int fcomm = 0 );

void jeveux_finalize( const ASTERINTEGER options = 0 );

void call_oper( py::object &syntax, int jxveri );

void call_oper_init();

void call_cmd_ctxt_enter();

void call_cmd_ctxt_exit();

void call_affich( const std::string &code, const std::string &text );

void call_print( const std::string &text );

void jeveux_delete( const std::string prefix );

bool deleteTemporaryObjects( void );

bool deleteCachedObjects( void );

std::string onFatalError( const std::string value = "" );

void call_matfpe( const int value );

int asmpi_get();

void asmpi_set( const int comm );

void asmpi_free( const int comm );

VectorInt asmpi_info( const int comm );

int asmpi_split( const int parent, int color, std::string name );

void set_option( const std::string &option, const ASTERINTEGER value );

// from aster_core_module
extern "C" {
void reset_tpmax();
extern PyObject *GetJdcAttr( _IN char * );
}

#endif
