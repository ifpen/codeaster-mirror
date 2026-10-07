/**
 * @file ContactNew.cxx
 * @brief Implementation de ParallelContactNew
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
#include "Contact/ParallelContactNew.h"

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

using VectorLongIter = VectorLong::iterator;

ParallelContactNew::ParallelContactNew( const std::string name, const ModelPtr model,
                                        ParallelMeshPtr mesh, const std::string type )
    : ContactNew( name, nullptr, type ), _pMesh( mesh ), _pModel( model ) {
    if ( model->getMesh() != mesh ) {
        throw std::runtime_error( "Mesh and model inconsistent" );
    }
};

bool ParallelContactNew::build() {
    CALL_JEMARQ();

    VectorString mastersSlaves, empty, slaves, masters;
    for ( auto &zone_i : _zones ) {
        const auto &sl = zone_i->getSlaveGroupOfCells();
        const auto &ma = zone_i->getMasterGroupOfCells();
        mastersSlaves.push_back( sl );
        mastersSlaves.push_back( ma );
        slaves.push_back( sl );
        masters.push_back( ma );
    }

    auto connectionMesh = std::make_shared< ConnectionMesh >( _pMesh, empty, mastersSlaves );
    _model = std::make_shared< Model >( connectionMesh );
    _model->setFrom( _pModel );

    ContactNew::build();
    _SFEDesc = _FEDesc;
    _PFEDesc = std::make_shared< ParallelContactFEDescriptor >( _FEDesc, connectionMesh, _model,
                                                                _pModel, masters, slaves );

    CALL_JEDEMA();

    return true;
}

#endif /* ASTER_HAVE_MPI */
