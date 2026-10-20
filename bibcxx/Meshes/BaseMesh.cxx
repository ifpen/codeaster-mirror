/**
 * @file BaseMesh.cxx
 * @brief Implementation de BaseMesh
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

#include "Meshes/BaseMesh.h"

#include "aster_fort_mesh.h"
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"

#include "DataFields/ConstantFieldOnCells.h"
#include "ParallelUtilities/AsterMPI.h"
#include "PythonBindings/LogicalUnitManager.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/CapyConvertibleValue.h"
#include "Utilities/Tools.h"

BaseMesh::BaseMesh( const std::string &name, const std::string &type )
    : DataStructure( name, 8, type ),
      ListOfTables( name ),
      _dimensionInformations( JeveuxVectorLong( getName() + ".DIME      " ) ),
      _coordinates( new MeshCoordinatesField( getName() + ".COORDO    " ) ),
      _nameOfGrpNodes( NamesMapChar24( getName() + ".PTRNOMNOE " ) ),
      _groupsOfNodes( JeveuxCollectionLongNamePtr( getName() + ".GROUPENO  ", _nameOfGrpNodes ) ),
      _connectivity( JeveuxContiguousCollectionLong( getName() + ".CONNEX    " ) ),
      _cellsType( JeveuxVectorLong( getName() + ".TYPMAIL   " ) ),
      _nameOfGrpCells( NamesMapChar24( getName() + ".PTRNOMMAI " ) ),
      _groupsOfCells( JeveuxCollectionLongNamePtr( getName() + ".GROUPEMA  ", _nameOfGrpCells ) ),
      _adapt( JeveuxVectorLong( getName() + ".ADAPTATION" ) ),
      _oriMeshName( JeveuxVectorChar8( getName() + ".MAOR" ) ),
      _resMeshCells( JeveuxVectorLong( getName() + ".CRMA" ) ),
      _resMeshNodes( JeveuxVectorLong( getName() + ".CRNO" ) ),
      _patch( JeveuxCollectionLong( getName() + ".PATCH" ) ),
      _nodePatchConnectivity( JeveuxVectorLong( getName() + ".CONOPA" ) ),
      _cellPatchConnectivity( JeveuxVectorLong( getName() + ".COMAPA" ) ),
      _namePatch( JeveuxVectorChar24( getName() + ".PTRNOMPAT" ) ),
      // use BaseMeshPtr(NULL) instead of this to avoid cross destruction
      _curvAbsc( new ConstantFieldOnCellsReal( getName().substr( 0, 8 ) + ".ABSC_CURV ",
                                               BaseMeshPtr( NULL ) ) ),
      _explorer( ConnectivityMeshExplorer( _connectivity, _cellsType ) ) {};

ASTERINTEGER BaseMesh::getNumberOfNodes() const {
    if ( isEmpty() )
        return 0;
    _dimensionInformations->updateValuePointer();
    return ( *_dimensionInformations )[0];
}

ASTERINTEGER BaseMesh::getNumberOfCells() const {
    if ( isEmpty() )
        return 0;
    _dimensionInformations->updateValuePointer();
    return ( *_dimensionInformations )[2];
}

ASTERINTEGER BaseMesh::getDimension() const {
    if ( isEmpty() )
        return 0;
    _dimensionInformations->updateValuePointer();
    const auto dimGeom = ( *_dimensionInformations )[5];

    if ( dimGeom == 3 ) {
        const std::string typeco( "MAILLAGE" );
        ASTERINTEGER repi = 0, ier = 0;
        JeveuxChar32 repk( " " );
        const std::string arret( "F" );
        const std::string questi( "DIM_GEOM" );

        CALLO_DISMOI( questi, getName(), typeco, &repi, repk, arret, &ier );

        return repi;
    }
    return dimGeom;
}

bool BaseMesh::readMeshFile( const std::filesystem::path &fileName, const std::string &format,
                             const int verbosity ) {
    FileType type = Ascii;
    if ( format == "MED" )
        type = Binary;
    LogicalUnitFile file1( fileName, type, Old );

    SyntaxMapContainer syntax;
    syntax.container["INFO"] = ASTERINTEGER( verbosity + 1 );

    if ( format == "GIBI" || format == "GMSH" ) {
        // Fichier temporaire
        LogicalUnitFile file2( "", Ascii, New );
        std::string preCmd = "PRE_" + format;
        ASTERINTEGER op2 = 47;
        if ( format == "GIBI" )
            op2 = 49;

        CommandSyntax *cmdSt2 = new CommandSyntax( preCmd );
        SyntaxMapContainer syntax2;
        syntax2.container["UNITE_" + format] = file1.getLogicalUnit();
        syntax2.container["UNITE_MAILLAGE"] = file2.getLogicalUnit();
        cmdSt2->define( syntax2 );

        CALL_EXECOP( &op2 );

        delete cmdSt2;
        syntax.container["FORMAT"] = "ASTER";
        syntax.container["UNITE"] = file2.getLogicalUnit();

        CommandSyntax cmdSt( "LIRE_MAILLAGE" );
        cmdSt.setResult( ResultNaming::getCurrentName(), "MAILLAGE" );

        cmdSt.define( syntax );

        ASTERINTEGER op = 1;
        CALL_EXECOP( &op );
    } else {
        syntax.container["FORMAT"] = format;
        syntax.container["UNITE"] = file1.getLogicalUnit();

        CommandSyntax cmdSt( "LIRE_MAILLAGE" );
        cmdSt.setResult( ResultNaming::getCurrentName(), "MAILLAGE" );

        cmdSt.define( syntax );

        ASTERINTEGER op = 1;
        CALL_EXECOP( &op );
    }

    return build();
}

const JeveuxCollectionLong BaseMesh::getMedConnectivity() const {
    JeveuxChar24 objv( " " );
    CALLO_GET_MED_CONNECTIVITY( getName(), objv );
    JeveuxCollectionLong result( objv.toString() );
    return result;
}

const JeveuxVectorLong BaseMesh::getMedCellsTypes() const {
    JeveuxChar24 objv( " " );
    CALLO_GET_MED_TYPES( getName(), objv );
    JeveuxVectorLong result( objv.toString() );
    return result;
}

bool BaseMesh::printMedFile( const std::filesystem::path &fileName, bool local ) const {
    const auto rank = getMPIRank();
    LogicalUnitFile a;
    ASTERINTEGER retour = -1;
    // In case that the print file (single and absolute path) is unique between processors,
    // it must only be created on proc 0.
    if ( isParallel() || ( !isParallel() && rank == 0 ) ) {
        a.openFile( fileName, Binary, New );
        retour = a.getLogicalUnit();
    }
    CommandSyntax cmdSt( "IMPR_RESU" );

    SyntaxMapContainer dict;

    if ( isParallel() || isConnection() )
        dict.container["PROC0"] = "NON";
    else
        dict.container["PROC0"] = "OUI";

    dict.container["FORMAT"] = "MED";
    dict.container["UNITE"] = retour;

    ListSyntaxMapContainer listeResu;
    SyntaxMapContainer dict2;
    dict2.container["MAILLAGE"] = getName();
    listeResu.push_back( dict2 );
    dict.container["RESU"] = listeResu;

    if ( !local && isParallel() )
        dict.container["FICHIER_UNIQUE"] = "OUI";

    cmdSt.define( dict );

    ASTERINTEGER op = 39;
    CALL_EXECOP( &op );

    return true;
};

const JeveuxCollectionLong BaseMesh::getInverseConnectivity() const {
    JeveuxChar24 objv( DataStructureNaming::getNewName( 24 ) );
    std::string base( "G" );

    ASTERINTEGER listCell;
    ASTERINTEGER nbCell = 0;

    CALLO_CNCINV( getName(), &listCell, &nbCell, base, objv );
    JeveuxCollectionLong result( objv.toString() );
    return result;
};

const std::vector< VectorLong > BaseMesh::getConnectivityZeroBased() const {
    std::vector< VectorLong > result;
    if ( !_connectivity->build() || _connectivity->size() < 0 ) {
        return result;
    }
    result.reserve( _connectivity->size() );
    for ( const auto &obj : _connectivity ) {
        obj->updateValuePointer();
        VectorLong items;
        items.reserve( obj->size() );
        for ( const auto &val : obj )
            items.push_back( val - 1 );
        result.push_back( items );
    }
    return result;
}

const std::vector< VectorLong > BaseMesh::getMedConnectivityZeroBased() const {
    std::vector< VectorLong > result;
    auto connectivity = getMedConnectivity();
    if ( !connectivity->build() || connectivity->size() < 0 ) {
        return result;
    }
    result.reserve( connectivity->size() );
    for ( const auto &obj : connectivity ) {
        obj->updateValuePointer();
        VectorLong items;
        items.reserve( obj->size() );
        for ( const auto &val : obj )
            items.push_back( val - 1 );
        result.push_back( items );
    }
    return result;
}

ASTERINTEGER BaseMesh::getCellType( const ASTERINTEGER &index ) const {
    if ( isEmpty() )
        return 0;
    if ( !_cellsType.exists() )
        return 0;
    _cellsType->updateValuePointer();
    return ( *_cellsType )[index];
};

JeveuxVectorLong BaseMesh::getCellsType() const {
    if ( _cellsType->exists() )
        _cellsType->updateValuePointer();
    return _cellsType;
};

ASTERINTEGER BaseMesh::getCellDime( const ASTERINTEGER &index ) const {
    ASTERINTEGER returnValue;

    auto cellType = getCellType( index );
    const std::string cata = "&CATA.TM.TMDIM";
    JeveuxChar32 objName, charName;
    CALLO_JEXNUM( objName, cata, &cellType );
    CALLO_JENONU( objName, &returnValue );

    return returnValue;
};

std::string BaseMesh::getCellTypeName( const ASTERINTEGER &index ) const {
    auto cellType = getCellType( index );
    const std::string cata = "&CATA.TM.NOMTM";
    JeveuxChar32 objName, charName;

    CALLO_JEXNUM( objName, cata, &cellType );
    CALLO_JENUNO( objName, charName );
    return strip( charName.toString() );
};

bool BaseMesh::hasCellsOfType( const std::string typma ) const {

    if ( isEmpty() )
        return false;
    if ( !_cellsType->exists() )
        return false;
    _cellsType->updateValuePointer();

    ASTERINTEGER typv;
    JeveuxChar32 objName( " " );
    std::string name = "&CATA.TM.NOMTM";
    CALLO_JEXNOM( objName, name, typma );
    CALLO_JENONU( objName, &typv );
    if ( typv <= 0 )
        return false;

    auto nbCells = _cellsType->size();
    for ( ASTERINTEGER i = 0; i < nbCells; i++ ) {
        if ( ( *_cellsType )[i] == typv )
            return true;
    }
    return false;
}

bool BaseMesh::isSkin( const std::string groupName ) const {

    if ( isParallel() ) {
        auto meshDime = getDimension();
        return true;
    } else {
        auto meshDime = getDimension();
        std::cout << "Mesh dime: " << meshDime << std::endl;
        if ( hasGroupOfCells( groupName ) ) {
            const VectorLong cells = getCells( groupName );
            const auto cellsType = getCellsType();
            for ( auto &cellId : cells ) {
                auto cellType = ( *cellsType )[cellId];
                auto cellDime = getCellDime( cellId );
                if ( cellDime != ( meshDime - 1 ) ) {
                    return false;
                }
            }

        } else {
            throw std::runtime_error( "The given group " + groupName + " doesn't exist in mesh" );
        }
    }
    return true;
}

bool BaseMesh::build() {
    _groupsOfNodes->build( true );
    _groupsOfCells->build( true );
    _patch->build();
    _connectivity->build();
    return update_tables();
}

void BaseMesh::initDefinition( const int &dim, const VectorReal &coord,
                               const VectorOfVectorsLong &connectivity, const VectorLong &types,
                               const int &nbGrpCells, const int &nbGrpNodes ) {
    int nbNodes = coord.size() / 3;
    int nbCells = types.size();

    AS_ASSERT( !_dimensionInformations.exists() );
    _dimensionInformations->allocate( 6 );
    ( *_dimensionInformations )[0] = nbNodes;
    ( *_dimensionInformations )[2] = nbCells;
    ( *_dimensionInformations )[5] = (ASTERINTEGER)dim;

    AS_ASSERT( !_cellsType.exists() );
    ( *_cellsType ) = types;

    AS_ASSERT( !_connectivity.exists() );
    _connectivity->allocate( connectivity );

    const JeveuxVectorReal values( "&&TMP", coord );
    _coordinates->assign( values );

    // created to the max capacity, groups may be added by several calls to addGroupsOfxxx
    if ( nbGrpCells > 0 ) {
        _groupsOfCells->allocate( nbGrpCells );
    }
    if ( nbGrpNodes > 0 ) {
        _groupsOfNodes->allocate( nbGrpNodes );
    }
}

bool BaseMesh::buildInformations( const int &dim ) {
    if ( _dimensionInformations->exists() )
        return false;
    if ( !_coordinates->exists() )
        throw std::runtime_error( "Coordinates vector must exist" );
    if ( !_cellsType->exists() )
        throw std::runtime_error( "Cells type vector must exist" );
    int nbNodes = _coordinates->size() / 3;
    int nbCells = _cellsType->size();

    _dimensionInformations->allocate( 6 );
    ( *_dimensionInformations )[0] = nbNodes;
    ( *_dimensionInformations )[2] = nbCells;
    ( *_dimensionInformations )[5] = (ASTERINTEGER)dim;
    return true;
}

bool BaseMesh::buildNamesVectors() { return true; }

void BaseMesh::show( const int verbosity ) const {
    ASTERINTEGER level( verbosity );
    CALLO_INFOMA( getName(), &level );
}

void BaseMesh::addGroupsOfNodes( const VectorString &names,
                                 const VectorOfVectorsLong &groupsOfNodes ) {
    int nbGroups = names.size();
    AS_ASSERT( nbGroups == groupsOfNodes.size() );
    if ( !_groupsOfNodes->exists() )
        _groupsOfNodes->allocate( nbGroups );
    AS_ASSERT( _groupsOfNodes->capacity() >= _groupsOfNodes->size() + nbGroups );

    for ( auto i = 0; i < nbGroups; ++i ) {
        _groupsOfNodes->push_back( names[i], groupsOfNodes[i] );
    }
}

void BaseMesh::addGroupsOfCells( const VectorString &names,
                                 const VectorOfVectorsLong &groupsOfCells ) {
    int nbGroups = names.size();
    AS_ASSERT( nbGroups == groupsOfCells.size() );
    if ( !_groupsOfCells->exists() )
        _groupsOfCells->allocate( nbGroups );
    AS_ASSERT( _groupsOfCells->capacity() >= _groupsOfCells->size() + nbGroups );

    for ( auto i = 0; i < nbGroups; ++i ) {
        _groupsOfCells->push_back( names[i], groupsOfCells[i] );
    }
}

void BaseMesh::endDefinition() {
    CALLO_CARGEO( getName() );
    AS_ASSERT( build() );
}
void BaseMesh::check( const ASTERDOUBLE tolerance ) {
    ASTERDOUBLE value = tolerance;
    CALLO_CHCKMA( getName(), &value );
}

void add_automatic_names( NamesMapChar8 &map, int size, std::string prefix ) {
    map->allocate( size );
    if ( size > 10000000 ) {
        for ( auto i = 0; i < size; ++i ) {
            std::ostringstream oss;
            oss << std::hex << i + 1;
            std::string name = prefix + toUpper( oss.str() );
            map->add( i + 1, name );
        }
    } else {
        for ( auto i = 0; i < size; ++i ) {
            map->add( i + 1, prefix + std::to_string( i + 1 ) );
        }
    }
}

const std::map< int, std::set< int > > &BaseMesh::buildReverseConnectivity() {
    if ( _bReverseConnex )
        return _reverseConnex;
    auto &meshConn = getConnectivity();
    meshConn->build();
    meshConn->updateValuePointer();
    const auto size = meshConn->size();
    for ( int elemId = 1; elemId <= size; ++elemId ) {
        auto &curCell = *( ( *meshConn )[elemId] );
        const auto nbNodes = curCell.size();
        for ( int j = 0; j < nbNodes; ++j ) {
            const auto &nodeId = curCell[j];
            _reverseConnex[nodeId - 1].insert( elemId - 1 );
        }
    }
    _bReverseConnex = true;
    return _reverseConnex;
};

void BaseMesh::deleteReverseConnectivity() {
    // free memory
    _reverseConnex = std::map< int, std::set< int > >();
    _bReverseConnex = false;
};

VectorLong BaseMesh::getRestrictedToOriginalNodesIds() const {

    if ( !_resMeshNodes->exists() ) {
        return VectorLong();
    }

    auto v_nodes = _resMeshNodes->toVector();
    for ( auto &val : v_nodes ) {
        // 0-based
        val -= 1;
    }

    return v_nodes;
}

VectorLong BaseMesh::getRestrictedToOriginalCellsIds() const {

    if ( !_resMeshCells->exists() ) {
        return VectorLong();
    }

    auto v_cells = _resMeshCells->toVector();
    for ( auto &val : v_cells ) {
        // 0-based
        val -= 1;
    }

    return v_cells;
}

MapLong BaseMesh::getOriginalToRestrictedNodesIds() const {

    MapLong nodes;

    if ( _resMeshNodes->exists() ) {
        _resMeshNodes->updateValuePointer();
        const int size = _resMeshNodes->size();
        for ( int i = 0; i < size; i++ ) {
            // 0-based
            nodes[( *_resMeshNodes )[i] - 1] = i;
        }
    }

    return nodes;
}

MapLong BaseMesh::getOriginalToRestrictedCellsIds() const {

    MapLong cells;

    if ( _resMeshCells->exists() ) {
        _resMeshCells->updateValuePointer();
        const int size = _resMeshCells->size();
        for ( int i = 0; i < size; i++ ) {
            // 0-based
            cells[( *_resMeshCells )[i] - 1] = i;
        }
    }

    return cells;
}

std::pair< ASTERDOUBLE, ASTERDOUBLE >
BaseMesh::getMinMaxEdgeSizes( const std::string cellGroupName ) {
    ASTERDOUBLE edgeMin, edgeMax;

    CALLO_MINMAXEDGESWRAP( this->getName(), cellGroupName, &edgeMin, &edgeMax );

    std::pair< ASTERDOUBLE, ASTERDOUBLE > returnValue = std::make_pair( edgeMin, edgeMax );

    return returnValue;
}
