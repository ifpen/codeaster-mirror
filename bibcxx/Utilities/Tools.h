#ifndef TOOLS_H_
#define TOOLS_H_

/**
 * @file Tools.h
 * @brief Fichier entete des outils
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

#include <algorithm>
#include <iomanip>
#include <sstream>

std::string strip( const std::string &str, const std::string &whitespace = " \t" );

std::string ljust( const std::string &str, const ASTERINTEGER &length, char fillchar = ' ' );

std::string toUpper( const std::string &in_str );

std::string toLower( const std::string &in_str );

std::string remove_brackets( const std::string &in_str );

VectorString split( const std::string &toSplit );

std::string strToupper( std::string s );

// Convert integer in string for name of object
std::string to_string( const int varInt, const int stringSize );

// wrapper arround dismoi;
std::tuple< bool, ASTERINTEGER, std::string > dismoi( const std::string &question,
                                                      const std::string &name,
                                                      const std::string &type, const bool stop );

/**
 * @brief irange Create a vector of integer from begin to end (included).
 *      for exemple {-1, 0, 1, 2, 3}
 */
VectorInt irange( const ASTERINTEGER4 begin, const ASTERINTEGER4 end );

VectorLong irange( const ASTERINTEGER begin, const ASTERINTEGER end );

/**
 * @brief vectorStringToFStrArray Fill an array of Fortran strings with a vector of strings.
 *      The caller must ensure that 'tabFStr' is big enough.
 */
void vectorStringToFStrArray( char *tabFStr, const int flen, const VectorString &vector );

/**
 * @brief vectorStringAsFStrArray Create an array of Fortran strings from a vector of strings.
 *      The output array must be freed by the caller.
 */
char *vectorStringAsFStrArray( const VectorString &vector, const int flen );

template < typename T >
std::vector< T > concatenate( const std::vector< std::vector< T > > &vec ) {

    ASTERINTEGER total_size = 0;
    for ( auto &v : vec ) {
        total_size += v.size();
    }

    std::vector< T > ret;
    ret.reserve( total_size );

    for ( auto &v : vec ) {
        for ( auto &elem : v ) {
            ret.push_back( elem );
        }
    }

    return ret;
}

// Set and sort a vector
template < typename T >
std::vector< T > unique( const std::vector< T > &vec ) {
    // make unique & sort
    std::set< T > s;
    std::copy( vec.begin(), vec.end(), std::inserter( s, s.end() ) );

    // recopy
    std::vector< T > r;
    r.resize( s.size() );
    std::copy( s.begin(), s.end(), r.begin() );

    return r;
}

// Get unique list of Aster Concept in a map (indexed by rank)
template < typename T, typename T2 >
std::vector< T > unique( const std::map< T2, T > &_map ) {
    std::map< std::string, T > unique_map;
    for ( auto it : _map ) {
        unique_map[it.second->getName()] = it.second;
    }

    std::vector< T > ret;
    ret.reserve( unique_map.size() );
    for ( auto it : unique_map ) {
        ret.push_back( it.second );
    }

    return ret;
}

/**
 * @brief Return all elements in common between the two vectors. Be carefull, input vectors are
 * modified inplace since a std::sort operation is performed. So create a copy before
 * if you don't want modifications of your inputs
 */
template < typename T >
std::vector< T > set_intersection( std::vector< T > &vec1, std::vector< T > &vec2 ) {
    std::vector< T > common;

    if ( vec1.empty() || vec2.empty() )
        return common;

    std::sort( vec1.begin(), vec1.end() );
    std::sort( vec2.begin(), vec2.end() );

    std::set_intersection( vec1.begin(), vec1.end(), vec2.begin(), vec2.end(),
                           std::back_inserter( common ) );

    return common;
}

template < typename T >
std::vector< T > set_union( std::vector< T > &vec1, std::vector< T > &vec2 ) {
    return unique( concatenate( std::vector< std::vector< T > >( { vec1, vec2 } ) ) );
}

/**
 * @brief Return all elements in vect1 that are not in vect2.
 */
template < typename T >
std::vector< T > set_difference( std::vector< T > &vec1, std::vector< T > &vec2 ) {
    std::vector< T > common;

    std::map< T, bool > remove;
    for ( auto &elem : vec2 ) {
        remove[elem] = true;
    }

    common.reserve( vec1.size() );
    for ( auto &elem : vec1 ) {
        if ( remove.count( elem ) == 0 ) {
            common.push_back( elem );
        }
    }

    return common;
}

template < typename T >
void print( const T &vec ) {
    std::cout << "Size: " << vec.size() << std::endl;
    for ( const auto &e : vec ) {
        std::cout << e << ", ";
    }
    std::cout << std::endl;
}

template < typename T >
std::set< T > toSet( const std::vector< T > &vec ) {
    std::set< T > set;
    std::copy( vec.begin(), vec.end(), std::inserter( set, set.end() ) );
    return set;
}

template < typename T >
std::vector< T > toVector( const std::set< T > &set ) {
    return std::vector< T >( set.begin(), set.end() );
}

#endif /* TOOLS_H_ */
