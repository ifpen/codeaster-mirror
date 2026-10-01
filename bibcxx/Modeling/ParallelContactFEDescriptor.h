
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#ifndef PARALLELCONTACTFEDESCRIPTOR_H_
#define PARALLELCONTACTFEDESCRIPTOR_H_

/**
 * @file ParallelContactFEDescriptor.h
 * @brief Fichier entete de la classe ParallelContactFEDescriptor
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

#include "Meshes/ConnectionMesh.h"
#include "Meshes/Joints.h"
#include "Modeling/Model.h"
#include "Modeling/ParallelFiniteElementDescriptor.h"

/**
 * @class ParallelContactFEDescriptor
 * @brief Classe definissant un ligrel parallèle
 * @author Nicolas Sellenet
 */
class ParallelContactFEDescriptor : public FiniteElementDescriptor {
  protected:
    /** @brief Matching numbering between keeped delayed elements and base elements */
    VectorLong _virtualCellToKeep;
    /** @brief All joints */
    JointsPtr _joints;
    /** @brief Delayed nodes owner */
    JeveuxVectorLong _owner;
    /** @brief Number of elements in which a given node is located */
    JeveuxVectorLong _multiplicity;
    /** @brief Number of non local elements in which a given node is located */
    JeveuxVectorLong _outerMultiplicity;
    /** @brief Global numbering for delayed nodes */
    JeveuxVectorLong _globalNumberingVirtualNodes;
    FiniteElementDescriptorPtr _FEDesc;
    /** @brief Element matching in element group list */
    VectorOfVectorsLong _cellMatching;
    /** @brief LIEL matching in element group list */
    VectorLong _lielMatching;
    /** @brief True if delayed node has component */
    JeveuxVectorLong _localVirtualIdToGlobalNodeId;
    JeveuxVectorLong _globalNodeIdToLocalVirtual;
    /** @brief connection mesh name */
    JeveuxVectorChar24 _cMeshName;

  public:
    /**
     * @brief Constructeur
     */
    ParallelContactFEDescriptor( const std::string &name, const FiniteElementDescriptorPtr &FEDesc,
                                 const ConnectionMeshPtr &mesh, const ModelPtr &connectionModel,
                                 const ModelPtr &model, const VectorString &masters,
                                 const VectorString &slaves );

    ParallelContactFEDescriptor( const std::string &name, const FiniteElementDescriptorPtr &FEDesc,
                                 const ConnectionMeshPtr &mesh, const ModelPtr &model );

    ParallelContactFEDescriptor( const FiniteElementDescriptorPtr &FEDesc,
                                 const ConnectionMeshPtr &mesh, const ModelPtr &connectionModel,
                                 const ModelPtr &model, const VectorString &masters,
                                 const VectorString &slaves );

    ParallelContactFEDescriptor( const FiniteElementDescriptorPtr &FEDesc,
                                 const ConnectionMeshPtr &mesh, const ModelPtr &model );

    /**
     * @brief Get vector of delayed elements keeped from the base FiniteElementDescriptor
     * @return reference on VectorLong
     */
    const VectorLong &getVirtualCellsToKeep() const { return _virtualCellToKeep; };

    /**
     * @brief Get vector of joints between subdomains
     * @return reference on VectorLong
     */
    const JeveuxVectorLong getJoints() const { return _joints->getOppositeDomains(); };
    std::string getJointObjectName() const { return _joints->getName(); };

    /**
     * @brief Get the mapping between local and global numbering of nodes
     * @return JeveuxVector of the indirection
     */
    const JeveuxVectorLong getLocalToGlobalMapping() const { return _globalNumberingVirtualNodes; }

    /**
     * @brief Get the mapping between local and global numbering of cells
     * @return JeveuxVector of the indirection
     */
    const VectorOfVectorsLong &getCellMatching() const { return _cellMatching; }

    const auto &getListOfGroupsOfElementsMatching() const { return _lielMatching; }

    /**
     * @typedef ParallelContactFEDescriptorPtr
     * @brief Pointeur intelligent vers un ParallelContactFEDescriptor
     */
    typedef std::shared_ptr< ParallelContactFEDescriptor > ParallelContactFEDescriptorPtr;

    FiniteElementDescriptorPtr getSupportFiniteElementDescriptor() const { return _FEDesc; }
};

/**
 * @typedef ParallelContactFEDescriptorPtr
 * @brief Pointeur intelligent vers un ParallelContactFEDescriptor
 */
typedef std::shared_ptr< ParallelContactFEDescriptor > ParallelContactFEDescriptorPtr;

#endif /* PARALLELCONTACTFEDESCRIPTOR_H_ */

#endif /* ASTER_HAVE_MPI */
