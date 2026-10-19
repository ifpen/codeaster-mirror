/**
 * @file MeshReader.cxx
 * @brief Implementation de MeshReader
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

#include "IOManager/MeshReader.h"

#include "Messages/Messages.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

#include <algorithm>

#ifdef ASTER_HAVE_MED

static const VectorInt asterTypeList = { 1,  2,  4,  6,  7,  9,  11, 12, 14, 16,
                                         18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };

static const std::map< int, med_int > asterMedMatching = {
    { 1, 1 },    { 2, 102 },  { 4, 103 },  { 6, 104 },  { 7, 203 },  { 9, 206 },  { 11, 207 },
    { 12, 204 }, { 14, 208 }, { 16, 209 }, { 18, 304 }, { 19, 310 }, { 20, 306 }, { 21, 315 },
    { 22, 318 }, { 23, 305 }, { 24, 313 }, { 25, 308 }, { 26, 320 }, { 27, 327 }
};

const std::set< med_int > medTypeToRenumber = { 304, 308, 305, 306, 310, 320, 313, 315, 318, 327 };

template < std::size_t N, const int indices[N] >
void applyPermutation( const med_int *in, med_int *out ) {
    for ( size_t i = 0; i < N; i++ ) {
        out[i] = in[indices[i]];
    }
};

std::vector< med_int > medToAsterRenumbering( const med_int &medType,
                                              const std::vector< med_int > &toRenumber,
                                              const int &nbElem ) {
    std::vector< med_int > out( toRenumber.size() );
    switch ( medType ) {
    case 304: {
        constexpr std::size_t N = 4;
        static constexpr int arr[N] { 0, 2, 1, 3 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 308: {
        constexpr std::size_t N = 8;
        static constexpr int arr[N] { 0, 3, 2, 1, 4, 7, 6, 5 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 305: {
        constexpr std::size_t N = 5;
        static constexpr int arr[N] { 0, 3, 2, 1, 4 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 306: {
        constexpr std::size_t N = 6;
        static constexpr int arr[N] { 0, 2, 1, 3, 5, 4 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 310: {
        constexpr std::size_t N = 10;
        static constexpr int arr[N] { 0, 2, 1, 3, 6, 5, 4, 7, 9, 8 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 320: {
        constexpr std::size_t N = 20;
        static constexpr int arr[N] { 0, 3, 2,  1,  4,  7,  6,  5,  11, 10,
                                      9, 8, 16, 19, 18, 17, 15, 14, 13, 12 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 313: {
        constexpr std::size_t N = 13;
        static constexpr int arr[N] { 0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 315: {
        constexpr std::size_t N = 15;
        static constexpr int arr[N] { 0, 2, 1, 3, 5, 4, 8, 7, 6, 12, 14, 13, 11, 10, 9 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 318: {
        constexpr std::size_t N = 18;
        static constexpr int arr[N] {
            0, 2, 1, 3, 5, 4, 8, 7, 6, 12, 14, 13, 11, 10, 9, 17, 16, 15
        };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 327: {
        constexpr std::size_t N = 27;
        static constexpr int arr[N] { 0,  3,  2,  1,  4,  7,  6,  5,  11, 10, 9,  8,  16, 19,
                                      18, 17, 15, 14, 13, 12, 20, 24, 23, 22, 21, 25, 26 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    }
    return out;
}

void MeshReader::readMeshFromMedFile( MeshPtr &toReturn, const std::filesystem::path &filename,
                                      const std::string &meshName, int verbosity ) {
    auto fr = MedFileReader();
    fr.open( filename, MedReadOnly );
    _readMesh( toReturn, fr, meshName, verbosity );
    toReturn->endDefinition();
}

#ifdef ASTER_HAVE_MPI
void MeshReader::readIncompleteMeshFromMedFile( IncompleteMeshPtr &toReturn,
                                                const std::filesystem::path &filename,
                                                const std::string &meshName, int verbosity ) {
    auto fr = MedFileReader();
    fr.openParallel( filename, MedReadOnly );
    _readMesh( toReturn, fr, meshName, verbosity );

    const auto curMesh = fr.getMesh( 0 );
    const auto &families = curMesh->getFamilies();
    for ( const auto &curFam : families ) {
        const auto &groups = curFam->getGroups();
        VectorString groupShort;
        for ( const auto &groupName : groups ) {
            if ( groupName.size() <= 24 ) {
                groupShort.push_back( groupName );
            }
        }
        toReturn->addFamily( curFam->getId(), groupShort );
    }
    const auto globNum = curMesh->getGlobalNodeNumberingAtSequence( -1, -1 );
    if ( globNum.size() != 0 ) {
        VectorLong arg( globNum.begin(), globNum.end() );
        toReturn->setLocalToGlobalNodeIds( arg );
    } else {
        VectorLong nodeGlobNum;
        const auto size = toReturn->getNumberOfNodes();
        const auto range = toReturn->getNodeRange();
        nodeGlobNum.reserve( size );
        for ( int i = 0; i < size; ++i ) {
            nodeGlobNum.push_back( i + range[0] );
        }
        toReturn->setLocalToGlobalNodeIds( nodeGlobNum );
    }
}

void MeshReader::readParallelMeshFromMedFile( ParallelMeshPtr &toReturn,
                                              const std::filesystem::path &filename,
                                              const std::string &meshName, int verbosity ) {
    auto fr = MedFileReader();
    fr.open( filename, MedReadOnly );
    _readMesh( toReturn, fr, meshName, verbosity );
    const auto rank = getMPIRank();

    // Read opposite domains, joint information and fill node owner vector
    const auto curMesh = fr.getMesh( 0 );
    const auto &joints = curMesh->getJoints();
    std::set< ASTERINTEGER > domainSet;
    std::map< std::string, VectorLong > allJointsMap;
    VectorLong nodeOwner( toReturn->getNumberOfNodes(), rank );
    for ( int i = 0; i < joints.size(); ++i ) {
        const auto &curJoint = joints[i];
        const auto &curName = curJoint->getName();
        if ( curJoint->getCorrespondenceNumber() != 1 || curJoint->getStepNumber() != 1 ) {
            throw std::runtime_error( "Unexpected joint in med file " + std::string( filename ) );
        }

        domainSet.insert( curJoint->getOppositeDomain() );

        const auto splitName = split( curName );
        AS_ASSERT( splitName.size() == 2 );
        const auto dIn = std::stoi( splitName[0] ), dOut = std::stoi( splitName[1] );
        AS_ASSERT( dIn == rank || dOut == rank );

        std::stringstream stream;
        stream << std::hex << curJoint->getOppositeDomain();
        std::string key = strToupper( ( dIn == rank ? "R" : "E" ) + std::string( stream.str() ) );
        const auto corresp = curJoint->getCorrespondence( 1, 1 );
        VectorLong corresp2;
        for ( const auto &tmp : corresp )
            corresp2.push_back( tmp );
        allJointsMap[key] = corresp2;

        if ( dIn == rank ) {
            for ( int j = 0; j < corresp.size() / 2; ++j ) {
                nodeOwner[corresp[2 * j] - 1] = -1;
            }
        }
    }
    auto domains = toVector( domainSet );
    std::sort( domains.begin(), domains.end() );

    // Read global node numbering
    const auto globNum = curMesh->getGlobalNodeNumberingAtSequence( -1, -1 );
    VectorLong globNum2;
    for ( const auto &num : globNum ) {
        globNum2.push_back( num );
    }

    VectorOfVectorsLong allJoints;
    for ( const auto &curDom : domains ) {
        std::stringstream stream;
        stream << std::hex << curDom;
        const std::string curDomStr( stream.str() );
        std::string eR( "ER" );
        for ( const auto &eOrR : eR ) {
            const std::string bis( 1, eOrR );
            const std::string key = strToupper( bis + curDomStr );
            allJoints.push_back( allJointsMap[key] );
        }
    }

    // Add parallel informations to mesh
    toReturn->create_joints( domains, globNum2, nodeOwner, {}, allJoints, 1 );
    toReturn->endDefinition();
}
#endif

void MeshReader::_readMesh( BaseMeshPtr toReturn, MedFileReader &fr, const std::string &meshName,
                            int verbosity ) {
#ifdef ASTER_HAVE_MPI
    const auto iM = std::dynamic_pointer_cast< IncompleteMesh >( toReturn );
    const bool incompleteMesh = ( iM ? true : false );
#else
    const bool incompleteMesh = false;
#endif
    auto coordsToFill = toReturn->getCoordinates();
    if ( coordsToFill->exists() ) {
        throw std::runtime_error( "not empty" );
    }
    auto coordValues = coordsToFill->getValues();

    // Read mesh from file
    auto curMeshId = -1;
    if ( meshName != "" ) {
        const auto meshNb = fr.getMeshNumber();
        for ( int meshId = 0; meshId < meshNb; ++meshId ) {
            const auto curMesh = fr.getMesh( meshId );
            if ( curMesh->getName() == meshName ) {
                curMeshId = meshId;
                break;
            }
        }
        if ( curMeshId == -1 ) {
            throw std::runtime_error( "Mesh " + meshName + " not found" );
        }
    } else {
        curMeshId = 0;
    }
    const auto curMesh = fr.getMesh( curMeshId );
    UtmessCore( "I", "MED_10", { curMesh->getName() } );
    const auto seq = curMesh->getSequence( 0 );
    const auto nodeNbAndStart = curMesh->getSplitNodeNumberAtSequence( seq[0], seq[1] );

    // Read node coordinates and copy in coordValues
    const auto nbNodes = nodeNbAndStart.first;
    const auto curCoords = curMesh->readCoordinates( seq[0], seq[1] );
    const auto dim = curMesh->getDimension();
    if ( dim == 3 ) {
        *( coordValues ) = curMesh->readCoordinates( seq[0], seq[1] );
    } else {
        const auto coordsToCopy = curMesh->readCoordinates( seq[0], seq[1] );
        coordValues->allocate( nbNodes * 3 );
        for ( int nodeId = 0; nodeId < nbNodes; ++nodeId ) {
            for ( int curDim = 0; curDim < dim; ++curDim ) {
                ( *coordValues )[nodeId * 3 + curDim] = coordsToCopy[nodeId * dim + curDim];
            }
        }
    }
    coordsToFill->buildDescriptor();

    // Get cell type and sort it according to aster sort
    auto cellTypes = curMesh->getGeometricTypesAtSequence( seq[0], seq[1] );
    std::set< med_int > typeSet( cellTypes.begin(), cellTypes.end() );
    std::vector< med_int > cellTypesSorted;
    VectorInt asterCellTypes;
    for ( const auto &asterType : asterTypeList ) {
        const auto &medType = asterMedMatching.at( asterType );
        if ( typeSet.count( medType ) != 0 ) {
            cellTypesSorted.push_back( medType );
            asterCellTypes.push_back( asterType );
        }
    }

    // Get cell informations by type
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    auto connectivity = toReturn->getConnectivity();
    auto cellType = toReturn->getCellTypeVector();
    int totalSize = 0, size = 0, cumCells = 0;
    VectorInt elemNbAndSizeVec;
    VectorOfVectorsLong cellRange;
    for ( const auto medType : cellTypesSorted ) {
        const auto totalCells =
            curMesh->getCellNumberForGeometricTypeAtSequence( seq[0], seq[1], medType );
        const auto cellNbAndStart =
            curMesh->getSplitCellNumberForGeometricTypeAtSequence( seq[0], seq[1], medType );
        if ( rank == nbProcs - 1 ) {
            cellRange.push_back( { cellNbAndStart.second - 1 + cumCells, totalCells + cumCells } );
        } else {
            cellRange.push_back( { cellNbAndStart.second - 1 + cumCells,
                                   cellNbAndStart.first * ( rank + 1 ) + cumCells } );
        }
        const auto cellNb = cellNbAndStart.first;
        const auto nbNodesForGeoT = curMesh->getNodeNumberForGeometricType( medType );
        totalSize += cellNb * nbNodesForGeoT;
        size += cellNb;
        cumCells += totalCells;
        elemNbAndSizeVec.push_back( cellNb );
        elemNbAndSizeVec.push_back( nbNodesForGeoT );
    }

#ifdef ASTER_HAVE_MPI
    if ( incompleteMesh ) {
        iM->setCellRange( cellRange );
    }
#endif
    // Get families in mesh
    const auto families = curMesh->getFamilies();
    med_int maxId = 0, minId = 0, nodeGrpCount = 0, cellGrpCount = 0;
    std::map< int, MedFamilyPtr > idToFamily;
    VectorString nodeGroupList, cellGroupList;
    std::map< std::string, int > nodeGroupNameToGroupId, cellGroupNameToGroupId;
    for ( const auto &fam : families ) {
        maxId = std::max( fam->getId(), maxId );
        minId = std::min( fam->getId(), minId );
        idToFamily[fam->getId()] = fam;
        const auto &curGroups = fam->getGroups();
        for ( const auto &grpName : curGroups ) {
            if ( fam->getId() > 0 ) {
                if ( nodeGroupNameToGroupId.count( grpName ) == 0 ) {
                    nodeGroupNameToGroupId[grpName] = nodeGrpCount;
                    nodeGroupList.push_back( grpName );
                    ++nodeGrpCount;
                }

            } else {
                if ( cellGroupNameToGroupId.count( grpName ) == 0 ) {
                    cellGroupNameToGroupId[grpName] = cellGrpCount;
                    cellGroupList.push_back( grpName );
                    ++cellGrpCount;
                }
            }
        }
    }
    const int familyOffset = -minId;
    maxId += familyOffset;
    VectorOfVectorsLong cellFamily( maxId + 1, VectorLong() );

    // Get node families
    auto curNFam = curMesh->getNodeFamilyAtSequence( seq[0], seq[1] );
    if ( incompleteMesh ) {
#ifdef ASTER_HAVE_MPI
        iM->setNodeRange(
            { nodeNbAndStart.second - 1, nodeNbAndStart.second - 1 + nodeNbAndStart.first } );
        VectorLong nodeFam( curNFam.begin(), curNFam.end() );
        iM->setNodeFamily( nodeFam );
#endif
    } else {
        // Group informations are only needed in non IncompleteMesh case
        auto index = 0;
        for ( const auto &cellFamId : curNFam ) {
            if ( cellFamId != 0 ) {
                cellFamily[cellFamId + familyOffset].push_back( index + 1 );
                if ( cellFamId < 0 ) {
                    throw std::runtime_error( "Error in MED families." );
                }
            }
            ++index;
        }
    }

    // Allocate connectivity
    connectivity->allocate( size, totalSize );
    VectorLong cellFam( size, 0 );
    auto cellFamStart = &cellFam[0];
    cellType->allocate( size );
    int count = 0, cumElem = 1, totalCount = 0, count2 = 0;
    ASTERINTEGER *connexPtr = nullptr;
    for ( const auto medType : cellTypesSorted ) {
        const auto &cellNb = elemNbAndSizeVec[count2 * 2];
        const auto &nbNodesForGeoT = elemNbAndSizeVec[count2 * 2 + 1];
        ++count2;
        if ( cellNb == 0 )
            continue;

        // Collection allocation
        for ( int cellId = 0; cellId < cellNb; ++cellId ) {
            if ( nbNodesForGeoT != 0 ) {
                connectivity->allocateObject( cumElem + cellId, nbNodesForGeoT );
            }
        }
        // Get contiguous collection start
        if ( count == 0 ) {
            connexPtr = connectivity->startPointer();
        }

        // Read connectivity from file
        auto curConn =
            curMesh->getConnectivityForGeometricTypeAtSequence( seq[0], seq[1], medType );
        // Apply renumbering
        if ( medTypeToRenumber.count( medType ) != 0 ) {
            curConn = medToAsterRenumbering( medType, curConn, cellNb );
        }
        // And copy to collection
        std::copy( curConn.begin(), curConn.end(), connexPtr + totalCount );
        curConn.clear();

        // Fill cell type
        auto cellTypePtr = &( ( *cellType )[cumElem - 1] );
        const auto &curAsterType = asterCellTypes[count2 - 1];
        std::fill( cellTypePtr, cellTypePtr + cellNb, curAsterType );

        // Get cell family
        auto curFam = curMesh->getCellFamilyForGeometricTypeAtSequence( seq[0], seq[1], medType );
        if ( incompleteMesh ) {
#ifdef ASTER_HAVE_MPI
            std::copy( curFam.begin(), curFam.end(), cellFamStart + cumElem - 1 );
#endif
        } else {
            auto index = 0;
            // Group informations are only needed in non IncompleteMesh case
            for ( const auto &cellFamId : curFam ) {
                if ( cellFamId != 0 ) {
                    if ( cellFamId > 0 ) {
                        throw std::runtime_error( "Error in MED families." );
                    }
                    cellFamily[cellFamId + familyOffset].push_back( index + cumElem );
                }
                ++index;
            }
        }

        cumElem += cellNb;
        totalCount += cellNb * nbNodesForGeoT;
        ++count;
    }

    if ( incompleteMesh ) {
#ifdef ASTER_HAVE_MPI
        iM->setCellFamily( cellFam );
#endif
    } else {
        int index = 0;
        VectorOfVectorsLong nodeIdGroupList( nodeGroupList.size(), VectorLong() );
        VectorOfVectorsLong cellIdGroupList( cellGroupList.size(), VectorLong() );
        // From families to groups
        for ( const auto &cellFamVector : cellFamily ) {
            const auto famId = index - familyOffset;
            if ( cellFamVector.size() != 0 ) {
                const auto &curFam = idToFamily.at( famId );
                const auto &curGroups = curFam->getGroups();
                for ( const auto &group : curGroups ) {
                    if ( famId > 0 ) {
                        const auto &grpId = nodeGroupNameToGroupId.at( group );
                        auto &curGrpToFill = nodeIdGroupList[grpId];
                        for ( const auto &id : cellFamVector ) {
                            curGrpToFill.push_back( id );
                        }
                    } else if ( famId < 0 ) {
                        const auto &grpId = cellGroupNameToGroupId.at( group );
                        auto &curGrpToFill = cellIdGroupList[grpId];
                        for ( const auto &id : cellFamVector ) {
                            curGrpToFill.push_back( id );
                        }
                    }
                }
            }
            ++index;
        }

        // Add non empty node groups
        VectorString nodeGroupList2;
        VectorOfVectorsLong nodeIdGroupList2;
        for ( int i = 0; i < nodeGroupList.size(); ++i ) {
            if ( nodeIdGroupList[i].size() != 0 ) {
                if ( nodeGroupList[i].size() <= 24 ) {
                    nodeGroupList2.push_back( nodeGroupList[i] );
                    std::sort( nodeIdGroupList[i].begin(), nodeIdGroupList[i].end() );
                    nodeIdGroupList2.push_back( nodeIdGroupList[i] );
                } else {
                    UtmessCore( "A", "MED_7", { nodeGroupList[i] } );
                }
            }
        }
        if ( nodeGroupList2.size() != 0 ) {
            toReturn->addGroupsOfNodes( nodeGroupList2, nodeIdGroupList2 );
        }

        // Add non empty cell groups
        VectorString cellGroupList2;
        VectorOfVectorsLong cellIdGroupList2;
        for ( int i = 0; i < cellGroupList.size(); ++i ) {
            if ( cellIdGroupList[i].size() != 0 ) {
                if ( cellGroupList[i].size() <= 24 ) {
                    cellGroupList2.push_back( cellGroupList[i] );
                    std::sort( cellIdGroupList[i].begin(), cellIdGroupList[i].end() );
                    cellIdGroupList2.push_back( cellIdGroupList[i] );
                } else {
                    UtmessCore( "A", "MED_7", { cellGroupList[i] } );
                }
            }
        }
        if ( cellGroupList2.size() != 0 ) {
            toReturn->addGroupsOfCells( cellGroupList2, cellIdGroupList2 );
        }
    }

    const auto dim2 = ( dim == 3 ? 3 : 2 );
    toReturn->buildInformations( dim2 );
    toReturn->buildNamesVectors();
    toReturn->show( verbosity );
}
#endif
