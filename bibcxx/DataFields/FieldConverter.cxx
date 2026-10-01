/**
 *   Copyright (C) 1991 2026  EDF www.code-aster.org
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

#include "DataFields/FieldConverter.h"

#include "ParallelUtilities/AsterMPI.h"

FieldOnNodesReal toFieldOnNodes( const MeshCoordinatesField &field, const BaseMeshPtr mesh ) {
    FieldOnNodesReal chamno = FieldOnNodesReal();

    std::string type = "NOEU", base = "G";
    std::string prol = "NON", model = " ";

    CALLO_CHPCHD( field.getName(), type, mesh->getName(), prol, base, chamno.getName(), model );

    chamno.build( mesh );

    return chamno;
}

FieldOnNodesReal getRealPart( const FieldOnNodesComplex &field ) {

    auto newValues = JeveuxVectorReal( "&&TMP", field.size() );

    field.updateValuePointers();
    std::transform( field.getValues()->begin(), field.getValues()->end(), newValues.begin(),
                    []( ASTERCOMPLEX db ) { return db.real(); } );

    FieldOnNodesReal newField =
        FieldOnNodesReal( field.getEquationNumbering(), field.getReference(), newValues );

    return newField;
};

FieldOnNodesReal getImaginaryPart( const FieldOnNodesComplex &field ) {

    auto newValues = JeveuxVectorReal( "&&TMP", field.size() );

    field.updateValuePointers();
    std::transform( field.getValues()->begin(), field.getValues()->end(), newValues.begin(),
                    []( ASTERCOMPLEX db ) { return db.imag(); } );

    FieldOnNodesReal newField =
        FieldOnNodesReal( field.getEquationNumbering(), field.getReference(), newValues );

    return newField;
};

FieldOnNodesReal getRealPart( const FieldOnNodesReal &field ) { return field; };

FieldOnNodesReal getImaginaryPart( const FieldOnNodesReal &field ) {

    FieldOnNodesReal newField = FieldOnNodesReal( field );
    newField.setValues( 0. );

    return newField;
};

#ifdef ASTER_HAVE_MPI
FieldOnNodesRealPtr transferToConnectionMesh( const FieldOnNodesRealPtr toTransfer,
                                              const ConnectionMeshPtr cMesh ) {
    auto sFON = toSimpleFieldOnNodes( toTransfer );
    sFON->updateValuePointers();
    AS_ASSERT( cMesh );

    const auto &comp = sFON->getComponents();
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    VectorString uniqueCmp;
    std::map< std::string, int > namePosIn, namePosAll;
    VectorLong indirection( comp.size(), 0. );

    VectorString allComp;
    for ( int i = 0; i < comp.size(); ++i )
        namePosIn[comp[i]] = i;
    AsterMPI::all_gather( comp, allComp );
    SetString addedCmp;
    for ( int i = 0; i < allComp.size(); ++i ) {
        const auto &cmpName = allComp[i];
        if ( addedCmp.count( cmpName ) == 0 ) {
            uniqueCmp.push_back( cmpName );
            namePosAll[cmpName] = uniqueCmp.size() - 1;
        }
        addedCmp.insert( cmpName );
    }
    for ( const auto &cmpName : comp ) {
        const auto &pos1 = namePosIn[cmpName];
        const auto &pos2 = namePosAll[cmpName];
        indirection[pos1] = pos2;
    }

    SimpleFieldOnNodesRealPtr rSFON = std::make_shared< SimpleFieldOnNodesReal >(
        cMesh, sFON->getPhysicalQuantity(), uniqueCmp, true );
    rSFON->updateValuePointers();

    const auto &localNum = cMesh->getNodesLocalNumbering();
    const auto &owners = cMesh->getNodesOwner();
    const auto &nbNodes = localNum->size();
    const auto &nbCmpIn = sFON->getNumberOfComponents();
    for ( int iNode = 0; iNode < nbNodes; ++iNode ) {
        const auto &curOwner = ( *owners )[iNode];
        if ( curOwner == rank ) {
            const auto &curPos = ( *localNum )[iNode] - 1;
            for ( int j = 0; j < nbCmpIn; ++j ) {
                const auto &cmpOut = indirection[j];
                ( *rSFON )( iNode, cmpOut ) = ( *sFON )( curPos, j );
            }
        }
    }

    const auto &values = rSFON->getValues();
    JeveuxVectorReal res( "&&TMP", values->size() );
    AsterMPI::all_reduce( values, res, MPI_SUM );

    const auto &nbCmpOut = rSFON->getNumberOfComponents();
    for ( int iNode = 0; iNode < nbNodes; ++iNode ) {
        for ( int j = 0; j < nbCmpOut; ++j ) {
            ( *rSFON )( iNode, j ) = ( *res )[iNode * nbCmpOut + j];
        }
    }
    return toFieldOnNodes( rSFON );
};

FieldOnNodesRealPtr transferFromConnectionToParallelMesh( const FieldOnNodesRealPtr toTransfer,
                                                          const BaseMeshPtr mesh ) {
    const auto sFON = toSimpleFieldOnNodes( toTransfer );
    sFON->updateValuePointers();
    const auto inputMesh = toTransfer->getMesh();
    // Check if the input mesh is a ConnectionMesh
    const auto cMesh = std::dynamic_pointer_cast< const ConnectionMesh >( inputMesh );
    if ( !cMesh ) {
        throw std::runtime_error( "Only valid for ConnectionMesh" );
    }
    // Check the mesh is a ParallelMesh
    const auto pMesh = std::dynamic_pointer_cast< ParallelMesh >( mesh );
    if ( !pMesh ) {
        throw std::runtime_error( "Only valid for ParallelMesh" );
    }

    SimpleFieldOnNodesRealPtr rSFON = std::make_shared< SimpleFieldOnNodesReal >(
        pMesh, sFON->getPhysicalQuantity(), sFON->getComponents(), true );
    rSFON->updateValuePointers();

    const auto &cLNum = cMesh->getNodesLocalNumbering();
    const auto &nbNodes = cLNum->size();
    const auto &lAllocated = sFON->getLogicalValues();
    const auto &nbCmpIn = sFON->getNumberOfComponents();

    for ( int inode = 0; inode < nbNodes; ++inode ) {
        if ( ( *lAllocated )[inode] && ( *cLNum )[inode] != -1 ) {
            const auto &curPos = ( *cLNum )[inode] - 1;
            for ( int icmp = 0; icmp < nbCmpIn; ++icmp ) {
                ( *rSFON )( curPos, icmp ) = ( *sFON )[inode * nbCmpIn + icmp];
            }
        }
    }

    auto ret = toFieldOnNodes( rSFON );

    return ret;
};
#endif /* ASTER_HAVE_MPI */
