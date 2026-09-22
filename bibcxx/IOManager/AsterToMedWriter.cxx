/**
 * @file AsterToMedWriter.cxx
 * @brief Implementation de AsterToMedWriter
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

#include "IOManager/AsterToMedWriter.h"

#include "aster_fort_mpi.h"

#include "IOManager/MedUtilities.h"
#include "Meshes/ParallelMesh.h"
#include "Messages/Messages.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Utilities/Tools.h"

#include <algorithm>

#ifdef ASTER_HAVE_MED

static const std::map< int, int > asterToAsterAvailableTypeInMed = {
    { 1, 1 },
    { 2, 2 },
    { 4, 4 },
    { 6, 6 },
    { 7, 7 },
    { 9, 9 },
    { 11, 11 },
    { 12, 12 },
    { 14, 14 },
    { 16, 16 },
    { 18, 18 },
    { 19, 19 },
    { 20, 20 },
    { 21, 21 },
    { 22, 22 },
    { 23, 23 },
    { 24, 24 },
    { 25, 25 },
    { 26, 26 },
    { 27, 27 },
    // Here begins unavailable types in med format
    { 28, 19 },
    { 29, 22 },
    { 30, 24 },
    { 70, 25 },
    { 71, 20 },
    { 87, 7 },
    { 88, 12 },
    { 89, 18 },
};

template < typename U, typename T >
std::vector< U > vectorFilter( const VectorLong &indexes, T &coords, const int &nbCmp ) {
    std::vector< U > toReturn( indexes.size() * nbCmp, 0 );
    int count = 0;
    for ( const auto &index : indexes ) {
        for ( int i = 0; i < nbCmp; ++i ) {
            toReturn[count * nbCmp + i] = coords[index * nbCmp + i];
        }
        ++count;
    }
    return toReturn;
};

std::vector< med_int > asterToMedRenumbering( const med_int &medType,
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
                                      9, 8, 19, 18, 17, 16, 12, 15, 14, 13 };
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
        static constexpr int arr[N] { 0, 2, 1, 3, 5, 4, 8, 7, 6, 14, 13, 12, 9, 11, 10 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 318: {
        constexpr std::size_t N = 18;
        static constexpr int arr[N] {
            0, 2, 1, 3, 5, 4, 8, 7, 6, 14, 13, 12, 9, 11, 10, 17, 16, 15
        };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    case 327: {
        constexpr std::size_t N = 27;
        static constexpr int arr[N] { 0,  3,  2,  1,  4,  7,  6,  5,  11, 10, 9,  8,  19, 18,
                                      17, 16, 12, 15, 14, 13, 20, 24, 23, 22, 21, 25, 26 };
        for ( int i = 0; i < nbElem; ++i ) {
            const int pos = i * N;
            applyPermutation< N, arr >( &toRenumber[pos], &out[pos] );
        }
        break;
    }
    }
    return out;
};

bool AsterToMedWriter::printMesh( const Mesh &toPrint, const std::filesystem::path &filename,
                                  bool local, const std::string &meshName,
                                  const std::string &mode ) {
    if ( !local ) {
        throw std::runtime_error( "Parallel print not allowed with Mesh" );
    }
    VectorLong nodeList =
        irange( (ASTERINTEGER)0, (ASTERINTEGER)( toPrint.getNumberOfNodes() - 1 ) );
    VectorLong cellList =
        irange( (ASTERINTEGER)0, (ASTERINTEGER)( toPrint.getNumberOfCells() - 1 ) );
    return _printMeshFromList( toPrint, filename, nodeList, cellList, local, meshName, mode );
};
#ifdef ASTER_HAVE_MPI
bool AsterToMedWriter::printMesh( const ConnectionMesh &toPrint,
                                  const std::filesystem::path &filename, bool local,
                                  const std::string &meshName, const std::string &mode ) {
    VectorLong nodeList =
        irange( (ASTERINTEGER)0, (ASTERINTEGER)( toPrint.getNumberOfNodes() - 1 ) );
    VectorLong cellList =
        irange( (ASTERINTEGER)0, (ASTERINTEGER)( toPrint.getNumberOfCells() - 1 ) );
    return _printMeshFromList( toPrint, filename, nodeList, cellList, local, meshName, mode );
};

bool AsterToMedWriter::printMesh( const ParallelMesh &toPrint,
                                  const std::filesystem::path &filename, bool local,
                                  const std::string &meshName, const std::string &mode ) {
    VectorLong nodeList, cellList;
    if ( local ) {
        nodeList = irange( (ASTERINTEGER)0, (ASTERINTEGER)( toPrint.getNumberOfNodes() - 1 ) );
        cellList = irange( (ASTERINTEGER)0, (ASTERINTEGER)( toPrint.getNumberOfCells() - 1 ) );
    } else {
        // In parallel print case, only inner cells and nodes are print
        // by each procs
        nodeList = toPrint.getInnerNodes();
        cellList = toPrint.getInnerCells();
    }
    auto cret = _printMeshFromList( toPrint, filename, nodeList, cellList, local, meshName, mode );
    const auto &nodeGN = toPrint.getLocalToGlobalNodeIds();
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    std::vector< med_int > globNum( nodeGN->begin(), nodeGN->end() );

    auto fr = MedFileReader();

    if ( local ) {
        fr.open( filename, MedReadWrite );
        const auto name = ( meshName != "" ) ? meshName : toPrint.getName();
        auto mesh = fr.getMesh( name );

        // print global numbering
        mesh->printGlobalNodeNumberingAtSequence( MED_NO_DT, MED_NO_IT, globNum );

        const auto rankStr = std::to_string( rank );
        const auto oppDom = toPrint.getOppositeDomains();
        int count = 0;
        // print joints
        for ( const auto &domdis : *oppDom ) {
            const auto dStr = std::to_string( domdis );
            const auto jName = ( rank < domdis ) ? rankStr + " " + dStr : dStr + " " + rankStr;
            const auto jDesc = ( rank < domdis ) ? "Receive joint" : "Send joint";
            std::vector< med_int > medJoint;
            if ( rank < domdis ) {
                const auto v1 = toPrint.getReceiveJoint( count );
                medJoint.assign( v1.begin(), v1.end() );
            } else {
                const auto v1 = toPrint.getSendJoint( count );
                medJoint.assign( v1.begin(), v1.end() );
            }
            mesh->createJoint( jName, jDesc, domdis );

            mesh->printNodeJointAtSequence( jName, MED_NO_DT, MED_NO_IT, medJoint );
            const auto jName2 = ( rank < domdis ) ? dStr + " " + rankStr : rankStr + " " + dStr;
            const auto jDesc2 = ( rank < domdis ) ? "Send joint" : "Receive joint";
            if ( rank < domdis ) {
                const auto v1 = toPrint.getSendJoint( count );
                medJoint.assign( v1.begin(), v1.end() );
            } else {
                const auto v1 = toPrint.getReceiveJoint( count );
                medJoint.assign( v1.begin(), v1.end() );
            }
            mesh->createJoint( jName2, jDesc2, domdis );
            mesh->printNodeJointAtSequence( jName2, MED_NO_DT, MED_NO_IT, medJoint );
            ++count;
        }
    } else {
        fr.openParallel( filename, MedReadWrite );
        const auto name = ( meshName != "" ) ? meshName : toPrint.getName();
        auto mesh = fr.getMesh( name );
        auto innerGNum = vectorFilter< med_int >( nodeList, globNum, 1 );

        std::vector< med_int > allGlobIds;
        AsterMPI::all_gather( innerGNum, allGlobIds );
        mesh->printGlobalNodeNumberingAtSequence( MED_NO_DT, MED_NO_IT, allGlobIds );
    }
    return cret;
};
#endif
bool AsterToMedWriter::_printMeshFromList( const BaseMesh &toPrint,
                                           const std::filesystem::path &filename,
                                           const VectorLong &nodeList, const VectorLong &cellList,
                                           bool local, const std::string &meshName,
                                           const std::string &mode ) {
    // Open med file
    auto fr = MedFileReader();
    const MedFileAccessType accessType = ( mode == "w" ) ? MedCreate : MedReadWrite;
    if ( !local ) {
        fr.openParallel( filename, accessType );
    } else {
        fr.open( filename, accessType );
    }
    const auto name = ( meshName != "" ) ? meshName : toPrint.getName();
    if ( fr.getMesh( name ) != nullptr ) {
        return false;
    }

    // mesh creation
    auto mesh = fr.createMesh( name, 3, "Cree par code_aster" );

    VectorOfVectorsLong cellByType;
    const auto &cellTypeVec = toPrint.getCellTypeVector();
    cellTypeVec->updateValuePointer();

    std::array< med_int, 3 > cumNodeNb;
    std::map< med_geometry_type, std::array< med_int, 3 > > mapMedTypeCellsByProc;
    // sort cells by type and create med filter informations
    if ( !local ) {
#ifdef ASTER_HAVE_MPI
        _buildFilterInformations( nodeList, cellList, cellTypeVec, cumNodeNb, mapMedTypeCellsByProc,
                                  cellByType );
#endif
    } else {
        _sortCellsByType( cellList, cellTypeVec, cellByType );
    }
    VectorLong medGNum;
    // save med filter informations and build med global numbering
    if ( !local ) {
        _createMedGlobalNumbering( medGNum, nodeList, cumNodeNb[2], toPrint.getNumberOfNodes() );
        ASTERINTEGER size = medGNum.size();
        CALL_VECTOR_GHOSTS_COMM( toPrint.getName(), &size, &medGNum[0] );
    } else {
        medGNum = VectorLong( toPrint.getNumberOfNodes(), 0 );
        for ( int i = 0; i < toPrint.getNumberOfNodes(); ++i ) {
            medGNum[i] = i;
        }
    }

    // node coordinates
    const auto &meshCoords = toPrint.getCoordinates()->getValues();
    auto coords = vectorFilter< med_float >( nodeList, ( *meshCoords ), 3 );
    const auto nodeNb = nodeList.size();
    MedFilterPtr medCoordsFilter( nullptr );
    if ( !local ) {
        medCoordsFilter = MedFilterPtr( new MedFilter(
            fr.getFilePointer(), cumNodeNb[0], 1, 3, MED_ALL_CONSTITUENT, MED_FULL_INTERLACE,
            MED_COMPACT_STMODE, cumNodeNb[2] + 1, cumNodeNb[1], 1, cumNodeNb[1], 0 ) );
    }
    mesh->printCoordinatesAtSequence( MED_NO_IT, MED_NO_IT, nodeNb, coords, medCoordsFilter );
    medCoordsFilter = nullptr;

    // 0 family creation (mandatory)
    const std::string familleZeroName( "FAMILLE_ZERO" );
    mesh->addFamily( familleZeroName, 0, {} );

    // create families in med file
    std::vector< med_int > allNodeFamily;
    // node families
    _createGroups( toPrint, mesh, allNodeFamily, nodeList, Nodes, local );
    MedFilterPtr medNodeFamFilter( nullptr );
    if ( !local ) {
        medNodeFamFilter = MedFilterPtr( new MedFilter(
            fr.getFilePointer(), cumNodeNb[0], 1, 1, MED_ALL_CONSTITUENT, MED_FULL_INTERLACE,
            MED_COMPACT_STMODE, cumNodeNb[2] + 1, cumNodeNb[1], 1, cumNodeNb[1], 0 ) );
    }
    mesh->setNodeFamilyAtSequence( MED_NO_DT, MED_NO_IT, allNodeFamily, medNodeFamFilter );
    allNodeFamily = std::vector< med_int >();

    std::vector< med_int > cellListFamily, allCellFamily( toPrint.getNumberOfCells(), 0 );
    // cell families
    _createGroups( toPrint, mesh, cellListFamily, cellList, Cells, local );
    int count = 0;
    for ( const auto &id : cellList ) {
        allCellFamily[id] = cellListFamily[count];
        ++count;
    }

    // cell connectivity and families
    const auto connectivity = toPrint.getConnectivity();
    connectivity->updateValuePointer();
    int curType = 1;
    const auto endIter = asterMedMatching.end();
    for ( const auto &cellIds : cellByType ) {
        const auto &curIter = asterMedMatching.find( curType );
        if ( curIter == endIter ) {
            ++curType;
            continue;
        }
        const auto &geotype = asterMedMatching.at( curType );
        bool empty = false, localEmpty = false;
        if ( !local ) {
            const auto &partDesc = mapMedTypeCellsByProc.at( geotype );
            if ( partDesc[0] == 0 )
                empty = true;
            if ( cellIds.size() == 0 )
                localEmpty = true;
        } else if ( cellIds.size() == 0 ) {
            empty = true;
        }
        if ( empty ) {
            ++curType;
            continue;
        }
        std::vector< med_int > curConn, curFamilies;
        int nodeNbInCell = 0;
        const auto cellNb = cellIds.size();
        if ( !localEmpty ) {
            const auto &firstCell = ( *connectivity )[cellIds[0] + 1];
            nodeNbInCell = firstCell->size();
            curConn = std::vector< med_int >( nodeNbInCell * cellNb, 0 );
            curFamilies = std::vector< med_int >( cellNb, 0 );
            int curPos = 0, posInFamilies = 0;
            for ( const auto &j : cellIds ) {
                const auto &curCell = ( *connectivity )[j + 1];
                for ( int k = 0; k < nodeNbInCell; ++k ) {
                    curConn[curPos] = medGNum[( *curCell )[k] - 1] + 1;
                    ++curPos;
                }
                curFamilies[posInFamilies] = allCellFamily[j];
                ++posInFamilies;
            }
            if ( medTypeToRenumber.count( geotype ) != 0 ) {
                curConn = asterToMedRenumbering( geotype, curConn, cellNb );
            }
        }
        if ( !local ) {
#ifdef ASTER_HAVE_MPI
            nodeNbInCell = AsterMPI::max( nodeNbInCell );
#endif
        }
        MedFilterPtr medCellFilter( nullptr ), medCellFamFilter( nullptr );
        // create filters for parallel print
        if ( !local ) {
            const auto &partDesc = mapMedTypeCellsByProc.at( geotype );
            med_int nbblocs = 1;
            if ( partDesc[1] == 0 ) {
                nbblocs = 0;
            }
            medCellFilter = MedFilterPtr(
                new MedFilter( fr.getFilePointer(), partDesc[0], 1, nodeNbInCell,
                               MED_ALL_CONSTITUENT, MED_FULL_INTERLACE, MED_COMPACT_STMODE,
                               partDesc[2] + 1, partDesc[1], nbblocs, partDesc[1], 0 ) );
            medCellFamFilter = MedFilterPtr( new MedFilter(
                fr.getFilePointer(), partDesc[0], 1, 1, MED_ALL_CONSTITUENT, MED_FULL_INTERLACE,
                MED_COMPACT_STMODE, partDesc[2] + 1, partDesc[1], nbblocs, partDesc[1], 0 ) );
        }
        mesh->printConnectivityForGeometricTypeAtSequence( MED_NO_DT, MED_NO_IT, geotype, cellNb,
                                                           curConn, medCellFilter );
        mesh->setCellFamilyAtSequence( MED_NO_DT, MED_NO_IT, geotype, curFamilies,
                                       medCellFamFilter );
        ++curType;
    }
    return true;
};

#ifdef ASTER_HAVE_MPI
std::set< std::set< std::string > > gatherFamilies( std::set< std::set< std::string > > in ) {
    VectorString grpVec, allGrpVec;
    VectorLong grpNum, allGrpNum;
    for ( const auto &fam : in ) {
        grpNum.push_back( fam.size() );
        for ( const auto &grpName : fam ) {
            grpVec.push_back( grpName );
        }
    }
    AsterMPI::all_gather( grpVec, allGrpVec );
    AsterMPI::all_gather( grpNum, allGrpNum );
    std::set< std::set< std::string > > out;
    int position = 0;
    for ( const auto &num : allGrpNum ) {
        std::set< std::string > toAdd;
        for ( int i = 0; i < num; ++i ) {
            toAdd.insert( allGrpVec[position] );
            ++position;
        }
        out.insert( toAdd );
    }
    return out;
};
#endif

void AsterToMedWriter::_createGroups( const BaseMesh &toPrint, MedMeshPtr mesh,
                                      std::vector< med_int > &allEntityFamily,
                                      const VectorLong &indexes, entityType entType, bool local ) {
    // this function is for nodes and cells
    const auto &grpNames =
        ( entType == Nodes ) ? toPrint.getGroupsOfNodes() : toPrint.getGroupsOfCells();
    int nodeNb = ( entType == Nodes ) ? toPrint.getNumberOfNodes() : toPrint.getNumberOfCells();
    std::vector< std::set< std::string > > allNodeGrp( nodeNb, std::set< std::string >() );
    std::set< std::set< std::string > > families;
    std::map< std::set< std::string >, int > famToInt;

    // loop over groups to get all node group
    for ( const auto &grpName : grpNames ) {
        const auto nodeIdVector =
            ( entType == Nodes ) ? toPrint.getNodes( grpName, true ) : toPrint.getCells( grpName );
        for ( const auto &nodeId : nodeIdVector ) {
            allNodeGrp[nodeId].insert( grpName );
        }
    }

    // build family set
    for ( const auto &family : allNodeGrp ) {
        families.insert( family );
    }
    // gather all families over procs
#ifdef ASTER_HAVE_MPI
    auto allFamilies = local ? families : gatherFamilies( families );
#else
    auto allFamilies = families;
#endif
    // add families in med mesh
    // and build famToInt: from family (set of groups) to family id
    int familyCount = 1;
    for ( const auto &family : allFamilies ) {
        if ( family.empty() ) {
            continue;
        }
        if ( entType == Nodes ) {
            famToInt[family] = familyCount;
        } else {
            famToInt[family] = -familyCount;
        }
        VectorString grpInFamily( family.begin(), family.end() );
        if ( entType == Nodes ) {
            const std::string familyname( "FAM_" + std::to_string( familyCount ) );
            mesh->addFamily( familyname, familyCount, grpInFamily );
        } else {
            const std::string familyname( "FAM_" + std::to_string( -familyCount ) );
            mesh->addFamily( familyname, -familyCount, grpInFamily );
        }
        ++familyCount;
    }
    // in med, empty family id must be 0
    famToInt[{}] = 0;

    // build allEntityFamily
    allEntityFamily = std::vector< med_int >( indexes.size(), 0 );
    for ( int i = 0; i < indexes.size(); ++i ) {
        const auto &nodeId = indexes[i];
        const auto &grps = allNodeGrp[nodeId];
        const auto famId = famToInt[grps];
        allEntityFamily[i] = famId;
    }
};

void AsterToMedWriter::_sortCellsByType( const VectorLong &cellVector, JeveuxVectorLong cellTypeVec,
                                         VectorOfVectorsLong &cellIdByType ) const {
    const auto maxCellType = *std::max_element( asterTypeList.begin(), asterTypeList.end() );
    cellIdByType = VectorOfVectorsLong( maxCellType, VectorLong() );
    // count cell by type
    for ( const auto &cellId : cellVector ) {
        const auto &cellType = ( *cellTypeVec )[cellId];
        auto cellTypeMod = asterToAsterAvailableTypeInMed.at( cellType );
        cellIdByType[cellTypeMod - 1].push_back( cellId );
    }
};
#ifdef ASTER_HAVE_MPI
void AsterToMedWriter::_buildFilterInformations(
    const VectorLong &nodeVector, const VectorLong &cellVector, JeveuxVectorLong cellTypeVec,
    std::array< med_int, 3 > &cumNodeNb,
    std::map< med_geometry_type, std::array< med_int, 3 > > &mapMedTypeCellsByProc,
    VectorOfVectorsLong &cellIdByType ) {
    const auto nbInnerNode = nodeVector.size();
    const auto rank = getMPIRank();
    const auto nbProcs = getMPISize();
    VectorInt nbNodes = { (int)nbInnerNode }, allNodeNb;
    std::vector< med_int > cumNodeNb0( nbProcs + 1, 0 );
    // compute all node number for all procs

    AsterMPI::all_gather( nbNodes, allNodeNb );
    for ( int i = 1; i <= nbProcs; ++i ) {
        cumNodeNb0[i] = cumNodeNb0[i - 1] + allNodeNb[i - 1];
    }
    cumNodeNb = { cumNodeNb0[nbProcs], cumNodeNb0[rank + 1] - cumNodeNb0[rank], cumNodeNb0[rank] };

    const auto maxCellType = *std::max_element( asterTypeList.begin(), asterTypeList.end() );
    _sortCellsByType( cellVector, cellTypeVec, cellIdByType );
    VectorInt cellNbByType, allCellNbByType;
    typedef std::vector< std::vector< med_int > > VectorOfVectorsMedInt;
    typedef std::vector< med_int > VectorMedInt;
    VectorOfVectorsMedInt cumCellNbByType( maxCellType, VectorMedInt( nbProcs + 1, 0 ) );
    for ( const auto &vec : cellIdByType ) {
        cellNbByType.push_back( vec.size() );
    }

    // gather cell by type number for all procs
    AsterMPI::all_gather( cellNbByType, allCellNbByType );

    // compute all cell number by type for all procs
    for ( int i = 0; i < maxCellType; ++i ) {
        for ( int j = 1; j <= nbProcs; ++j ) {
            const auto &refAllCellNbByType = allCellNbByType[( j - 1 ) * maxCellType + i];
            const auto &lastcumCellNbByType = cumCellNbByType[i][j - 1];
            auto &curcumCellNbByType = cumCellNbByType[i][j];
            curcumCellNbByType = lastcumCellNbByType + refAllCellNbByType;
        }
    }
    // store it in a map
    const auto &end = asterMedMatching.end();
    for ( int i = 0; i < maxCellType; ++i ) {
        const auto &curIter = asterMedMatching.find( i + 1 );
        if ( curIter != end ) {
            const auto &toCopy = cumCellNbByType[i];
            const auto medType = asterMedMatching.at( i + 1 );
            mapMedTypeCellsByProc[medType] = { toCopy[nbProcs], toCopy[rank + 1] - toCopy[rank],
                                               toCopy[rank] };
        }
    }
};
#endif

void AsterToMedWriter::_createMedGlobalNumbering( VectorLong &globNum, const VectorLong &innerNodes,
                                                  int startNum, int localNodeNumber ) {
    globNum = VectorLong( localNodeNumber, 0 );
    auto curNum = startNum;
    for ( const auto &nodeId : innerNodes ) {
        globNum[nodeId] = curNum;
        ++curNum;
    }
    return;
};

bool AsterToMedWriter::printResult( const ResultPtr &resu, const std::filesystem::path &filename,
                                    bool local, const std::string &mode ) {
    // first, print mesh if needed
    const auto meshInResu = resu->getMesh();
    if ( meshInResu->isParallel() ) {
#ifdef ASTER_HAVE_MPI
        const auto pMesh0 = std::dynamic_pointer_cast< ParallelMesh >( meshInResu );
        printMesh( pMesh0, filename, local, mode );
#endif
    } else {
        const auto mesh0 = std::dynamic_pointer_cast< Mesh >( meshInResu );
        printMesh( mesh0, filename, local, mode );
    }

    // open med file
    auto fr = MedFileReader();

    VectorLong nodeList, cellList;
    VectorOfVectorsLong cellByType;
    if ( !local ) {
        fr.openParallel( filename, MedReadWrite );
    } else {
        fr.open( filename, MedReadWrite );
    }
    const auto &medMesh = fr.getMesh( meshInResu->getName() );
    std::array< med_int, 3 > cumNodeNb;
    std::map< med_geometry_type, std::array< med_int, 3 > > mapMedTypeCellsByProc;
    if ( !local ) {
        const auto &cellTypeVec = meshInResu->getCellTypeVector();
        cellTypeVec->updateValuePointer();
        nodeList = meshInResu->getInnerNodes();
        cellList = meshInResu->getInnerCells();
        // sort cells by type and create med filter informations
#ifdef ASTER_HAVE_MPI
        _buildFilterInformations( nodeList, cellList, cellTypeVec, cumNodeNb, mapMedTypeCellsByProc,
                                  cellByType );
#endif
    }

    const auto indexes = resu->getIndexes();
    const auto userName = resu->getUserName();

    MedFieldPtr medField( nullptr );
    const auto nodeFieldNameList = resu->getFieldsOnNodesRealNames();
    const auto nodeNb = meshInResu->getNumberOfNodes();
    for ( const auto &fieldName : nodeFieldNameList ) {
        bool first = true;
        const std::string fName =
            std::string( "________" ).replace( 0, userName.size(), userName ) + fieldName;
        for ( const auto &index : indexes ) {
            const auto &curField = resu->getFieldOnNodesReal( fieldName, index );
            const auto sFON = toSimpleFieldOnNodes( curField );
            const auto cmpsVector = sFON->getComponents();
            if ( first ) {
                medField = fr.createField( fName, medMesh, cmpsVector );
                first = false;
            }
            medField->addSequence( index, index, resu->getTime( index ), MED_NO_DT, MED_NO_IT );
            auto valuesDescPair = sFON->getValuesWithDescription( cmpsVector, ( VectorString ) {} );
            if ( !local ) {
                auto values = vectorFilter< ASTERDOUBLE >( nodeList, valuesDescPair.first,
                                                           cmpsVector.size() );
                auto valuesFilter = MedFilterPtr(
                    new MedFilter( fr.getFilePointer(), cumNodeNb[0], 1, cmpsVector.size(),
                                   MED_ALL_CONSTITUENT, MED_FULL_INTERLACE, MED_COMPACT_STMODE,
                                   cumNodeNb[2] + 1, cumNodeNb[1], 1, cumNodeNb[1], 0 ) );
                medField->printRealValuesAtSequenceOnNodes( index, index, nodeNb, values,
                                                            valuesFilter );
            } else {
                medField->printRealValuesAtSequenceOnNodes( index, index, nodeNb,
                                                            valuesDescPair.first );
            }
        }
    }
    return true;
};

bool AsterToMedWriter::printMesh( const MeshPtr &mesh, const std::filesystem::path &filename,
                                  bool local, const std::string &meshName,
                                  const std::string &mode ) {
    return printMesh( *mesh, filename, local, meshName );
};

#ifdef ASTER_HAVE_MPI
bool AsterToMedWriter::printMesh( const ParallelMeshPtr &mesh,
                                  const std::filesystem::path &filename, bool local,
                                  const std::string &meshName, const std::string &mode ) {
    return printMesh( *mesh, filename, local, meshName );
};

bool AsterToMedWriter::printMesh( const ConnectionMeshPtr &mesh,
                                  const std::filesystem::path &filename, bool local,
                                  const std::string &meshName, const std::string &mode ) {
    return printMesh( *mesh, filename, local, meshName );
};
#endif

#endif
