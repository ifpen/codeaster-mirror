#ifndef ASTERCXX_H_
#define ASTERCXX_H_

/* ==================================================================== */
/* Copyright (C) 1991 - 2026  EDF www.code-aster.org */
/*                                                                      */
/* This file is part of Code_Aster.                                     */
/*                                                                      */
/* Code_Aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* Code_Aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* ==================================================================== */

#include "aster.h"

#include "asterc_config.h"

#include <algorithm>
#include <array>
#include <complex>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

using ASTERBOOL = bool;
using ASTERCOMPLEX = std::complex< ASTERDOUBLE >;

using VectorBool = std::vector< ASTERBOOL >;
using VectorInt = std::vector< ASTERINTEGER4 >;
using VectorLong = std::vector< ASTERINTEGER >;
using VectorReal = std::vector< ASTERDOUBLE >;
using VectorComplex = std::vector< ASTERCOMPLEX >;
using VectorString = std::vector< std::string >;

using VectorOfVectorsLong = std::vector< VectorLong >;
using VectorOfVectorsReal = std::vector< VectorReal >;

using SetInt = std::set< ASTERINTEGER4 >;
using SetLong = std::set< ASTERINTEGER >;
using SetString = std::set< std::string >;

using PairInt = std::pair< ASTERINTEGER4, ASTERINTEGER4 >;
using PairLong = std::pair< ASTERINTEGER, ASTERINTEGER >;
using PairReal = std::pair< ASTERDOUBLE, ASTERDOUBLE >;
using PairString = std::pair< std::string, std::string >;

using VectorPairInt = std::vector< PairInt >;
using VectorPairLong = std::vector< PairLong >;
using VectorPairReal = std::vector< PairReal >;
using VectorPairString = std::vector< PairString >;

using MapLong = std::map< ASTERINTEGER, ASTERINTEGER >;
using MapString = std::map< std::string, std::string >;
using MapVectorLong = std::map< ASTERINTEGER, VectorLong >;

#define AS_ABORT( message )                                                                        \
    DEBUG_LOC;                                                                                     \
    std::cout << message << std::endl;                                                             \
    INTERRUPT( 17 );

// Exceptions identifiers - keep consistency with asterf.h
#define ASTER_ERROR 1
#define ASTER_CONVERGENCE_ERROR 2
#define ASTER_INTEGRATION_ERROR 3
#define ASTER_SOLVER_ERROR 4
#define ASTER_CONTACT_ERROR 5
#define ASTER_TIMELIMIT_ERROR 6

#endif // ASTERCXX_H_
