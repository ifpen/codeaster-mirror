#pragma once

/**
 * @file BaseMesh.h
 * @brief Fichier entete de la classe BaseMesh
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

#include "DataFields/ListOfTables.h"
#include "DataFields/MeshCoordinatesField.h"
#include "DataStructures/DataStructure.h"
#include "MemoryManager/NamesMap.h"
#include "Meshes/MeshExplorer.h"
#include "Utilities/GenericEnum.h"

#include <filesystem>

/** @brief Forward declaration of ConstantFieldOnCells */
template < class ValueType >
class ConstantFieldOnCells;
using ConstantFieldOnCellsReal = ConstantFieldOnCells< ASTERDOUBLE >;
using ConstantFieldOnCellsRealPtr = std::shared_ptr< ConstantFieldOnCellsReal >;

/**
 * @class BaseMesh
 * @brief This object is the base class for all meshes variants
 */
class BaseMesh : public DataStructure, public ListOfTables {
  private:
    static unsigned long int _cell_idx;
    static unsigned long int _node_idx;

  public:
    using ConnectivityMeshExplorer =
        MeshExplorer< CellsIteratorFromConnectivity, const JeveuxContiguousCollectionLong &,
                      const JeveuxVectorLong & >;

  protected:
    using JeveuxCollectionLongNamePtr = JeveuxCollection< ASTERINTEGER, NamesMapChar24 >;
    /** @brief Objet Jeveux '.DIME' */
    JeveuxVectorLong _dimensionInformations;
    /** @brief Champ aux noeuds '.COORDO' */
    MeshCoordinatesFieldPtr _coordinates;
    /** @brief Pointeur de nom Jeveux '.PTRNOMNOE' */
    NamesMapChar24 _nameOfGrpNodes;
    /** @brief Collection Jeveux '.GROUPENO' */
    JeveuxCollectionLongNamePtr _groupsOfNodes;
    /** @brief Collection Jeveux '.CONNEX' */
    JeveuxContiguousCollectionLong _connectivity;
    /** @brief Objet Jeveux '.TYPMAIL' */
    JeveuxVectorLong _cellsType;
    /** @brief Pointeur de nom Jeveux '.PTRNOMMAI' */
    NamesMapChar24 _nameOfGrpCells;
    /** @brief Objet Jeveux '.GROUPEMA' */
    JeveuxCollectionLongNamePtr _groupsOfCells;
    /** @brief jeveux vector '.ADAPTATION' */
    JeveuxVectorLong _adapt;
    /** @brief jeveux vector '.MAOR' */
    JeveuxVectorChar8 _oriMeshName;
    /** @brief jeveux vector '.CRMA' */
    JeveuxVectorLong _resMeshCells;
    /** @brief jeveux vector '.CRNO' */
    JeveuxVectorLong _resMeshNodes;
    /** @brief Collection jeveux '.PATCH' */
    JeveuxCollectionLong _patch;
    /** @brief jeveux vector '.CONOPA' */
    JeveuxVectorLong _nodePatchConnectivity;
    /** @brief jeveux vector '.COMAPA' */
    JeveuxVectorLong _cellPatchConnectivity;
    /** @brief jeveux vector '.PTRNOMPAT' */
    JeveuxVectorChar24 _namePatch;
    /** @brief card '.ABSC_CURV' */
    ConstantFieldOnCellsRealPtr _curvAbsc;
    /** @brief Object to allow loop over connectivity */
    const ConnectivityMeshExplorer _explorer;
    /** @brief Reverse connectivity */
    std::map< int, std::set< int > > _reverseConnex;
    bool _bReverseConnex = false;

    /**
     * @brief Constructeur
     * @param name nom jeveux de l'objet
     * @param type jeveux de l'objet
     */
    BaseMesh( const std::string &name, const std::string &type );

    /**
     * @brief Read a Mesh file
     * @return return true if it succeeds, false otherwise
     */
    bool readMeshFile( const std::filesystem::path &fileName, const std::string &format,
                       const int verbosity );

  public:
    using BaseMeshPtr = std::shared_ptr< BaseMesh >;

    virtual void addFamily( int id, VectorString groups ) {};

    /**
     * @brief Get the connectivity
     */
    const ConnectivityMeshExplorer &getConnectivityExplorer() const {
        _cellsType->updateValuePointer();
        _connectivity->build();
        return _explorer;
    };

    /**
     * @brief Return the connectivity
     */
    const JeveuxContiguousCollectionLong getConnectivity() const { return _connectivity; }

    /**
     * @brief Return cell type vector
     */
    const JeveuxVectorLong getCellTypeVector() const { return _cellsType; }

    /**
     * @brief Return the connectivity, node zero based
     */
    const std::vector< VectorLong > getConnectivityZeroBased() const;

    virtual std::vector< VectorString > getCellFamilyGroups() const {
        return std::vector< VectorString >();
    };

    virtual std::vector< VectorString > getNodeFamilyGroups() const {
        return std::vector< VectorString >();
    };

    const JeveuxCollectionLong getInverseConnectivity() const;

    /**
     * @brief Return the connectivity with MED numberings
     */
    const JeveuxContiguousCollectionLong getMedConnectivity() const;

    /**
     * @brief Return the connectivity with MED numberings, node zero based
     */
    const std::vector< VectorLong > getMedConnectivityZeroBased() const;

    /**
     * @brief Return the MED type for each cell
     */
    const JeveuxVectorLong getMedCellsTypes() const;

    /**
     * @brief Recuperation des coordonnees du maillage
     * @return champ aux noeuds contenant les coordonnees des noeuds du maillage
     */
    MeshCoordinatesFieldPtr getCoordinates() const {
        _coordinates->updateValuePointers();
        return _coordinates;
    };

    ConstantFieldOnCellsRealPtr getCurvilinearAbscissa() const { return _curvAbsc; }

    /**
     * @brief Get all the names of group of cells
     * @return NamesMapChar24 _nameOfGrpCells
     */
    const NamesMapChar24 &getGroupsOfNodesMap() const { return _nameOfGrpCells; };

    /**
     * @brief Returns the number of nodes
     */
    ASTERINTEGER getNumberOfNodes() const;

    /**
     * @brief Returns the number of cells
     */
    ASTERINTEGER getNumberOfCells() const;

    ASTERINTEGER getCellType( const ASTERINTEGER &index ) const;

    virtual VectorLong getCellFamily() const { return VectorLong(); };

    virtual VectorLong getNodeFamily() const { return VectorLong(); };

    JeveuxVectorLong getCellsType() const;

    ASTERINTEGER getCellDime( const ASTERINTEGER &index ) const;

    std::string getCellTypeName( const ASTERINTEGER &index ) const;

    bool hasCellsOfType( const std::string ) const;

    bool isSkin( const std::string groupName ) const;

    /**
     * @brief Recuperation de la dimension du maillage
     */
    virtual ASTERINTEGER getDimension() const;

    /**
     * @brief Teste l'existence d'un groupe de mailles dans le maillage
     * @return true si le groupe existe
     */
    virtual bool hasGroupOfCells( const std::string &name, const bool local = false ) const {
        AS_ASSERT( false );
        return false;
    };

    /**
     * @brief Teste l'existence d'un groupe de noeuds dans le maillage
     * @return true si le groupe existe
     */
    virtual bool hasGroupOfNodes( const std::string &name, const bool local = false ) const {
        AS_ASSERT( false );
        return false;
    };

    /**
     * @brief Returns the names of the groups of cells
     * @return VectorString
     */
    virtual VectorString getGroupsOfCells( const bool local = false ) const {
        AS_ASSERT( false );
        return {};
    };

    /**
     * @brief Returns the names of the groups of nodes
     * @return VectorString
     */
    virtual VectorString getGroupsOfNodes( const bool local = false ) const {
        AS_ASSERT( false );
        return {};
    };

    /**
     * @brief Create a group of cells
     */
    virtual void setGroupOfCells( const std::string &name, const VectorLong &cell_ids ) {
        AS_ASSERT( false );
    };

    /**
     * @brief Create a group of nodes
     */
    virtual void setGroupOfNodes( const std::string &name, const VectorLong &node_ids,
                                  const bool localNumbering = false ) {
        AS_ASSERT( false );
    };
    /**
     * @brief Steal from input vector ids in local numbering of ghost nodes on the last layer
     */
    virtual void setLastGhostsLayer( const VectorInt &node_ids ) { AS_ASSERT( false ); };
    /**
     * @brief Returns the cells indexes of a group of cells
     * @return VectorLong
     */
    virtual VectorLong getCells( const std::string name ) const {
        AS_ASSERT( false );
        return {};
    }

    virtual VectorLong getCells( const VectorString &names = {} ) const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Return ids in local numbering of ghost nodes on the last layer
     */
    virtual JeveuxVectorShort getLastGhostsLayer() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Returns the nodes indexes of a group of nodes
     * @return VectorLong
     */

    virtual VectorLong getNodes( const VectorString &names, const bool localNumbering = false,
                                 const ASTERINTEGER same_rank = PythonBool::None ) const;

    virtual VectorLong getNodes( const std::string name = std::string(),
                                 const bool localNumbering = false,
                                 const ASTERINTEGER same_rank = PythonBool::None ) const;

    /**
     * @brief Returns the nodes indexes of a group of cells
     * @return VectorLong
     */
    virtual VectorLong getNodesFromCells( const std::string name, const bool localNumbering = false,
                                          const ASTERINTEGER same_rank = PythonBool::None ) const {
        AS_ASSERT( false );
        return {};
    }

    virtual VectorLong getNodesFromCells( const VectorString &names,
                                          const bool localNumbering = false,
                                          const ASTERINTEGER same_rank = PythonBool::None ) const {
        AS_ASSERT( false );
        return {};
    }

    virtual VectorLong getNodesFromCells( const VectorLong &cells,
                                          const bool localNumbering = false,
                                          const ASTERINTEGER same_rank = PythonBool::None ) const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Returns the nodes indexes of inner nodes
     * @return VectorLong
     */
    virtual VectorLong getInnerNodes() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Returns the nodes indexes of inner nodes
     * @return VectorLong
     */
    virtual VectorLong getInnerCells() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Returns the nodes indexes of outer nodes
     * @return VectorLong
     */
    virtual VectorLong getOuterNodes() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Get the JeveuxVector for outer subdomain nodes
     * @return VectorLong
     */
    virtual const JeveuxVectorLong getNodesOwner() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Get the JeveuxVector for outer subdomain cells
     * @return VectorLong
     */
    virtual const JeveuxVectorLong getCellsOwner() const {
        AS_ASSERT( false );
        return {};
    }

    VectorOfVectorsLong getNodesRanks() const {
        AS_ASSERT( false );
        return {};
    }

    VectorOfVectorsLong getCellsRanks() const {
        AS_ASSERT( false );
        return {};
    }

    /**
     * @brief Fonction permettant de savoir si un maillage est vide (non relu par exemple)
     * @return retourne true si le maillage est vide
     */
    bool isEmpty() const { return !_dimensionInformations.exists(); };

    /**
     * @brief Fonction permettant de savoir si un maillage est parallel
     * @return retourne true si le maillage est parallel
     */
    virtual bool isParallel() const { return false; };

    /**
     * @brief Fonction permettant de savoir si un maillage est complet
     * @return retourne true si le maillage est complet
     */
    virtual bool isIncomplete() const { return false; };

    /**
     * @brief Fonction permettant de savoir si un maillage est partiel
     * @return retourne true si le maillage est partiel
     */
    virtual bool isConnection() const { return false; };

    /**
     * @brief Tester le maillage a des cells quadratiques
     * @return true si quadratique
     */
    virtual bool isQuadratic() const {
        AS_ASSERT( false );
        return false;
    };

    /**
     * @brief Impression du maillage au format MED
     * @param fileName Nom du fichier MED à imprimer
     * @param local booleen pour demander une impression locale
     * @param version version de med de l'impression
     * @param mode mode d'ecriture 'a' ou 'w'
     * @return true
     */
    virtual bool printMedFile( const std::filesystem::path &fileName, bool local = true,
                               std::array< int, 3 > version = { 0, 0, 0 },
                               const std::string &mode = "a" ) const;

    /**
     * @brief Get the mapping between local and global numbering of nodes
     * @return JeveuxVector of the indirection
     */
    virtual const JeveuxVectorLong getLocalToGlobalNodeIds() const {
        AS_ASSERT( false );
        return {};
    }

    virtual ASTERINTEGER getGlobalToLocalNodeId( const ASTERINTEGER &nodeId ) {
        AS_ASSERT( false );
        return -1;
    }

    virtual std::shared_ptr< const MapLong > getGlobalToLocalNodeIds() const {
        AS_ASSERT( false );
        return {};
    };

    virtual const JeveuxVectorLong getLocalToGlobalCellIds() const {
        AS_ASSERT( false );
        return {};
    }

    VectorLong getRestrictedToOriginalNodesIds() const;
    VectorLong getRestrictedToOriginalCellsIds() const;

    MapLong getOriginalToRestrictedNodesIds() const;
    MapLong getOriginalToRestrictedCellsIds() const;

    /**
     * @brief Build the mesh
     * @return true if success
     */
    bool build();

    const std::map< int, std::set< int > > &buildReverseConnectivity();

    void deleteReverseConnectivity();

    /* Mesh builder functions */
    void initDefinition( const int &dim, const VectorReal &coord,
                         const VectorOfVectorsLong &connectivity, const VectorLong &types,
                         const int &nbGrpCells, const int &nbGrpNodes );

    bool buildInformations( const int &dim );

    bool buildNamesVectors();

    void addGroupsOfNodes( const VectorString &names, const VectorOfVectorsLong &groupsOfNodes );

    void addGroupsOfCells( const VectorString &names, const VectorOfVectorsLong &groupsOfCells );

    void endDefinition();

    void show( const int verbosity = 1 ) const;

    void check( const ASTERDOUBLE tolerance );

    std::pair< ASTERDOUBLE, ASTERDOUBLE > getMinMaxEdgeSizes( const std::string cellGroupName );
};

using BaseMeshPtr = std::shared_ptr< BaseMesh >;

/**
 * @brief Helper function to automatically name nodes and cells
 */
void add_automatic_names( NamesMapChar8 &map, int size, std::string prefix );
