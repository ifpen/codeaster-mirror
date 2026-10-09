#ifndef FINITEELEMENTDESCRIPTOR_H_
#define FINITEELEMENTDESCRIPTOR_H_

/**
 * @file FiniteElementDescriptor.h
 * @brief Fichier entete de la classe FiniteElementDescriptor
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

#include "DataStructures/DataStructure.h"
#include "MemoryManager/JeveuxCollection.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/MeshExplorer.h"
#include "Modeling/Partition.h"

/**
 * @class FiniteElementDescriptor
 * @brief Class which describes the finite elements
 * @author Nicolas Sellenet
 */
class FiniteElementDescriptor : public DataStructure {
  public:
    using ConnectivityVirtualCellsExplorer = MeshExplorer< CellsIteratorFromFiniteElementDescriptor,
                                                           const JeveuxContiguousCollectionLong & >;

  protected:
    /** @brief Vecteur Jeveux '.NBNO' */
    JeveuxVectorLong _numberOfDelayedNumberedConstraintNodes;
    /** @brief Vecteur Jeveux '.LGRF' */
    JeveuxVectorChar8 _parameters;
    /** @brief Vecteur Jeveux '.PRNM' */
    JeveuxVectorLong _dofDescriptor;
    /** @brief Collection '.LIEL' */
    JeveuxContiguousCollectionLong _listOfGroupsOfElements;
    /** @brief Vecteur Jeveux '.REPE' */
    JeveuxVectorLong _groupsOfCellsNumberByElement;
    /** @brief Collection '.NEMA' */
    JeveuxContiguousCollectionLong _virtualCellsDescriptor;
    /** @brief Vecteur Jeveux '.PRNS' */
    JeveuxVectorLong _dofOfDelayedNumberedConstraintNodes;
    /** @brief Vecteur Jeveux '.LGNS' */
    JeveuxVectorLong _virtualNodesNumbering;
    /** @brief Vecteur Jeveux '.SSSA' */
    JeveuxVectorLong _superElementsDescriptor;
    /** @brief Vecteur Jeveux '.NVGE' */
    JeveuxVectorChar16 _nameOfNeighborhoodStructure;
    /** @brief Vecteur Jeveux '.TYFE' */
    JeveuxVectorLong _typeOfCells;
    /** @brief Base mesh */
    BaseMeshPtr _mesh;
    /** @brief Object to loop over connectivity of delayed numbered cells */
    const ConnectivityVirtualCellsExplorer _explorer;
    /** @brief Object to loop over list of group of cells */
    const ConnectivityVirtualCellsExplorer _explorer2;
    PartitionPtr _partition;

    /**
     * @brief Constructeur
     */
    FiniteElementDescriptor( const std::string &name, const std::string &type,
                             const BaseMeshPtr mesh );

  public:
    /**
     * @typedef FiniteElementDescriptorPtr
     * @brief Pointeur intelligent vers un FiniteElementDescriptor
     */
    using FiniteElementDescriptorPtr = std::shared_ptr< FiniteElementDescriptor >;

    /**
     * @brief Constructeur
     */
    FiniteElementDescriptor( const std::string &name, const BaseMeshPtr mesh );

    FiniteElementDescriptor( const std::string &name, const FiniteElementDescriptor &fed );

    FiniteElementDescriptor( const BaseMeshPtr mesh );

    FiniteElementDescriptor( const FiniteElementDescriptorPtr FEDesc,
                             const VectorString &groupOfCells );

    /**
     * @brief Destructor
     */
    ~FiniteElementDescriptor() {};

    FiniteElementDescriptor &operator=( const FiniteElementDescriptor &toCopy );

    void addVirtualCells( const FiniteElementDescriptor &other );

    const ConnectivityVirtualCellsExplorer &getVirtualCellsExplorer() const;

    const JeveuxVectorLong &getVirtualNodesComponentDescriptor() const;

    const JeveuxVectorLong &getVirtualNodesNumbering() const;

    const ConnectivityVirtualCellsExplorer &getListOfGroupsOfElementsExplorer() const;

    const JeveuxContiguousCollectionLong &getListOfGroupsOfElements() const;

    const JeveuxContiguousCollectionLong &getVirtualCellsDescriptor() const;

    ASTERINTEGER getNumberOfVirtualNodes() const;

    ASTERINTEGER getNumberOfGroupsOfElements() const;

    ASTERINTEGER getNumberOfElementsOfGroupOfElements( const ASTERINTEGER grel_idx ) const;

    void setNumberOfVirtualNodes( const ASTERINTEGER nbNodes );

    JeveuxVectorLong getNumberOfVirtualNodesDescriptor();

    JeveuxVectorChar8 getParameters() const;

    const JeveuxVectorLong &getPhysicalNodesComponentDescriptor() const;

    const JeveuxVectorLong &getListOfGroupsOfElementsbyElement() const;

    const BaseMeshPtr getMesh() const;

    void setMesh( const BaseMeshPtr &currentMesh );

    int getPhysics( void ) const;

    bool exists() const;

    bool build();

    /**@brief patition method */
    const std::string getPartitionMethod() const;

    ASTERINTEGER getNumberOfCells() const;

    FiniteElementDescriptorPtr restrict( const VectorString &groupsOfCells ) const;

    FiniteElementDescriptorPtr restrict( const VectorLong &cells ) const;

    /** @brief Get index of elem type */
    ASTERINTEGER getElemTypeNume( const std::string &elemTypeName ) const;

    std::string getElemTypeName( const ASTERINTEGER &elemTypeNume ) const;

    JeveuxVectorLong getFiniteElementType() const;

    ASTERINTEGER getFiniteElementTypeIdOfGroupOfElements( const ASTERINTEGER grel_idx ) const;

    SetString getFiniteElementNames() const;

    /**
     * @brief Number of super-elements in model
     * @return Number of super elements in model
     */
    ASTERINTEGER numberOfSuperElement() const;

    /**@brief Has super elements in model ? */
    bool existsSuperElement() const;

    /** @brief Has finite element in model ? */
    bool existsFiniteElement() const;

    bool existsHHO() const;

#ifdef ASTER_HAVE_MPI
    /** @brief Transert .PRNM from other FiniteElementDescriptor.
     * this should be associated to a ConnectionMesh,
     * other should be associated to the parallelMesh of the ConnectionMesh */
    void transferDofDescriptorFrom( FiniteElementDescriptorPtr & );

    void setFrom( FiniteElementDescriptorPtr & );

    void transferListOfGroupOfCellFrom( FiniteElementDescriptorPtr & );

#endif /* ASTER_HAVE_MPI */
};

/**
 * @typedef FiniteElementDescriptor
 * @brief Pointeur intelligent vers un FiniteElementDescriptor
 */
using FiniteElementDescriptorPtr = std::shared_ptr< FiniteElementDescriptor >;

#endif /* FINITEELEMENTDESCRIPTOR_H_ */
