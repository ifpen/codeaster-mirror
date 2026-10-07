/**
 * @file ContactPairing.cxx
 * @brief Implementation de ParallelContactPairing
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

#include "Contact/ParallelContactPairing.h"

#include "aster_fort_ds.h"
#include "aster_fort_jeveux.h"
#include "aster_fort_mesh.h"
#include "aster_fort_utils.h"

#include "Meshes/ConnectionMesh.h"
#include "Messages/Messages.h"
#include "Modeling/Model.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

#ifdef ASTER_HAVE_MPI

ParallelContactPairing::ParallelContactPairing( const std::string name,
                                                const ParallelContactNewPtr cont )
    : ContactPairing( name, cont ), _pContNew( cont ) {};

void ParallelContactPairing::buildFiniteElementDescriptor() {
    ContactPairing::buildFiniteElementDescriptor();
    _SFEDesc = _fed;
    const auto cMesh = std::dynamic_pointer_cast< ConnectionMesh >( _pContNew->getMesh() );
    if ( cMesh ) {
        _PFEDesc = std::make_shared< ParallelContactFEDescriptor >( _fed, cMesh,
                                                                    _pContNew->getParallelModel() );
    } else {
        throw std::runtime_error( "Problem with ConnectionMesh while contact pairing" );
    }
}

#endif /* ASTER_HAVE_MPI */
