#pragma once

/**
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

/* FIXME aster.h should not include Python.h */
#include "astercxx.h"

#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;

/** @brief Factory for '__init__' constructor without 'DSTypePtr' */
template < typename DSType, typename... Args >
static std::shared_ptr< DSType > initFactoryPtr( Args... args ) {
    return std::make_shared< DSType >( args... );
};

/** @brief Defines pickling functions */
/* If the object is a DataStructure, it must inherit from DSWithCppPickling.
 */
template < typename DSType >
static auto define_pickling() {
    return py::pickle( []( const DSType &obj ) { return obj._getState(); },
                       []( const py::tuple &tup ) { return std::make_shared< DSType >( tup ); } );
};

/** @brief Explicitly marks that pickling is not supported by an object */
template < typename DSType >
static auto disable_pickling() {
    return []( const DSType &obj ) { return 1; };
};

namespace pybind11 {
namespace detail {

/** @brief Converter for JeveuxVector */
template < typename T >
struct type_caster< JeveuxVector< T > > {
  public:
    PYBIND11_TYPE_CASTER( JeveuxVector< T >, const_name( "JeveuxVector" ) );

    bool load( handle /* src */, bool ) { return false; }

    static handle cast( const JeveuxVector< T > &vect, return_value_policy /* policy */,
                        handle /* parent */ ) {
        py::list pylist;
        if ( vect.exists() ) {
            vect->updateValuePointer();
            auto size = vect->size();
            for ( int i = 0; i < size; ++i ) {
                pylist.append( ( *vect )[i] );
            }
        }
        return pylist.inc_ref();
    }
};

/** @brief Converter for JeveuxVector */
template < int lengthT >
struct type_caster< JeveuxVector< JeveuxString< lengthT > > > {
  public:
    PYBIND11_TYPE_CASTER( JeveuxVector< JeveuxString< lengthT > >, const_name( "JeveuxVector" ) );

    bool load( handle /* src */, bool ) { return false; }

    static handle cast( const JeveuxVector< JeveuxString< lengthT > > &vect,
                        return_value_policy /* policy */, handle /* parent */ ) {
        py::list pylist;
        if ( vect.exists() ) {
            vect->updateValuePointer();
            auto size = vect->size();
            for ( int i = 0; i < size; ++i ) {
                pylist.append( ( *vect )[i].rstrip() );
            }
        }
        return pylist.inc_ref();
    }
};

/** @brief Converter for JeveuxCollection */
template < typename T >
struct type_caster< JeveuxCollection< T > > {
  public:
    PYBIND11_TYPE_CASTER( JeveuxCollection< T >, const_name( "JeveuxCollection" ) );

    bool load( handle /* src */, bool ) { return false; }

    static handle cast( const JeveuxCollection< T > &coll, return_value_policy /* policy */,
                        handle /* parent */ ) {
        py::list pylist;
        if ( !coll->build() || coll->size() < 0 ) {
            return pylist.inc_ref();
        }
        for ( const auto &obj : coll ) {
            obj->updateValuePointer();
            py::list items;
            for ( const auto &val : obj ) {
                items.append( val );
            }
            items.inc_ref();
            pylist.append( items );
        }
        return pylist.inc_ref();
    }
};

/** @brief Converter for JeveuxContiguousCollection */
template < typename T >
struct type_caster< JeveuxCollection< T, ASTERINTEGER, Contiguous > > {
  public:
    typedef JeveuxCollection< T, ASTERINTEGER, Contiguous > ContiguousCollection;
    PYBIND11_TYPE_CASTER( ContiguousCollection, const_name( "JeveuxContiguousCollection" ) );

    bool load( handle /* src */, bool ) { return false; }

    static handle cast( const JeveuxCollection< T, ASTERINTEGER, Contiguous > &coll,
                        return_value_policy /* policy */, handle /* parent */ ) {
        py::list pylist;
        if ( !coll->build() || coll->size() < 0 ) {
            return pylist.inc_ref();
        }
        for ( const auto &obj : coll ) {
            obj->updateValuePointer();
            py::list items;
            for ( const auto &val : obj ) {
                items.append( val );
            }
            items.inc_ref();
            pylist.append( items );
        }
        return pylist.inc_ref();
    }
};

} // namespace detail
} // namespace pybind11

// ----------------------------------------------------------------------------------------
// This code comes from dolfinx https://github.com/FEniCS/dolfinx
// pybind11 casters for PETSc/petsc4py objects
#ifdef ASTER_HAVE_PETSC4PY
#include <petsc4py/petsc4py.h>
#include <petscdm.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscvec.h>

// Import petsc4py on demand
#define VERIFY_PETSC4PY( func )                                                                    \
    if ( !func ) {                                                                                 \
        if ( import_petsc4py() != 0 ) {                                                            \
            std::cout << "ERROR: could not import petsc4py!" << std::endl;                         \
            throw std::runtime_error( "Error when importing petsc4py" );                           \
        }                                                                                          \
    }

// Macro for casting between dolfin and petsc4py objects
#define PETSC_CASTER_MACRO( TYPE, NAME )                                                           \
    template <>                                                                                    \
    class type_caster< _p_##TYPE > {                                                               \
      public:                                                                                      \
        PYBIND11_TYPE_CASTER( TYPE, _( #NAME ) );                                                  \
        bool load( handle src, bool ) {                                                            \
            VERIFY_PETSC4PY( PyPetsc##TYPE##_Get );                                                \
            if ( PyObject_TypeCheck( src.ptr(), &PyPetsc##TYPE##_Type ) == 0 )                     \
                return false;                                                                      \
            value = PyPetsc##TYPE##_Get( src.ptr() );                                              \
            return true;                                                                           \
        }                                                                                          \
                                                                                                   \
        static handle cast( TYPE src, pybind11::return_value_policy policy, handle parent ) {      \
            VERIFY_PETSC4PY( PyPetsc##TYPE##_New );                                                \
            auto obj = PyPetsc##TYPE##_New( src );                                                 \
            if ( policy == py::return_value_policy::take_ownership )                               \
                PetscObjectDereference( (PetscObject)src );                                        \
            else if ( policy == py::return_value_policy::reference_internal )                      \
                keep_alive_impl( obj, parent );                                                    \
            return py::handle( obj );                                                              \
        }                                                                                          \
                                                                                                   \
        operator TYPE() { return value; }                                                          \
    }

namespace pybind11 {
namespace detail {
// PETSC_CASTER_MACRO(DM, dm);
PETSC_CASTER_MACRO( KSP, ksp );
PETSC_CASTER_MACRO( Mat, mat );
// PETSC_CASTER_MACRO(SNES, snes);
PETSC_CASTER_MACRO( Vec, vec );
} // namespace detail
} // namespace pybind11
#endif

#undef PETSC_CASTER_MACRO

// ----------------------------------------------------------------------------------------
