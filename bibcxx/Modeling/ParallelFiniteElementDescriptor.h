
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#ifndef PARALLELFINITEELEMENTDESCRIPTOR_H_
#define PARALLELFINITEELEMENTDESCRIPTOR_H_

/**
 * @file ParallelFiniteElementDescriptor.h
 * @brief Fichier entete de la classe ParallelFiniteElementDescriptor
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
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"

/**
 * @class ParallelFiniteElementDescriptor
 * @brief Classe definissant un ligrel parallèle
 * @author Nicolas Sellenet
 */
class ParallelFiniteElementDescriptor : public FiniteElementDescriptor {
  protected:
    /** @brief Matching numbering between keeped delayed elements and base elements */
    VectorLong _virtualCellToKeep;
    /** @brief Join to send */
    std::vector< JeveuxVectorLong > _joinToSend;
    /** @brief Join to receive */
    std::vector< JeveuxVectorLong > _joinToReceive;
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

  public:
    /**
     * @brief Constructeur
     */
    ParallelFiniteElementDescriptor( const std::string &name,
                                     const FiniteElementDescriptorPtr &FEDesc,
                                     const ConnectionMeshPtr &mesh, const ModelPtr &model );

    /**
     * @brief Constructeur
     */
    ParallelFiniteElementDescriptor( const std::string &name, const std::string &jName,
                                     const BaseMeshPtr &mesh );

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
     * @typedef ParallelFiniteElementDescriptorPtr
     * @brief Pointeur intelligent vers un ParallelFiniteElementDescriptor
     */
    typedef std::shared_ptr< ParallelFiniteElementDescriptor > ParallelFiniteElementDescriptorPtr;
};

/**
 * @typedef ParallelFiniteElementDescriptorPtr
 * @brief Pointeur intelligent vers un ParallelFiniteElementDescriptor
 */
typedef std::shared_ptr< ParallelFiniteElementDescriptor > ParallelFiniteElementDescriptorPtr;

#endif /* PARALLELFINITEELEMENTDESCRIPTOR_H_ */

#endif /* ASTER_HAVE_MPI */
