#ifndef BASEELEMENTARYMATRIX_H_
#define BASEELEMENTARYMATRIX_H_

/**
 * @file ElementaryMatrix.h
 * @brief Definition of elementary matrices
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

#include "DataFields/ElementaryTerm.h"
#include "DataStructures/DataStructure.h"
#include "Discretization/ElementaryCompute.h"
#include "Loads/MechanicalLoad.h"
#include "Loads/PhysicalQuantity.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Supervis/ResultNaming.h"

/**
 * @class BaseElementaryMatrix
 * @brief Generic class for sd_matr_elem
 */
class BaseElementaryMatrix : public DataStructure {
  protected:
    /** @brief Flag for empty datastructure (either empty or built)*/
    bool _isBuilt;

    /** @brief Model */
    ModelPtr _model;

    /** @brief Elementary compute */
    ElementaryComputePtr _elemComp;

  public:
    /** @brief Constructor with a name */
    BaseElementaryMatrix( const std::string name, const std::string type = "MATR_ELEM" );

    /** @brief Constructor with automatic name */
    BaseElementaryMatrix( const std::string type = "MATR_ELEM" );

    BaseElementaryMatrix( const ModelPtr model );

    /** @brief Get the model */
    ModelPtr getModel() const { return _model; };

    /** @brief Get option */
    std::string getOption() const { return _elemComp->getOption(); };

    /**
     * @brief Detect state of datastructure
     * @return true if datastructure is not empty
     */
    bool isBuilt() const { return _isBuilt; };

    /**
     * @brief Set state of datastructure
     * @param bBuilt flag for state of datastructure
     */
    void isBuilt( bool bBuilt ) { _isBuilt = bBuilt; };

    /**
     * @brief Set the model
     * @param currModel pointer to model
     */
    void setModel( const ModelPtr &currModel ) { _model = currModel; };

    /** @brief  Prepare compute */
    void prepareCompute( const std::string option );

    bool isSymmetric() const;

    ASTERINTEGER numberOfSuperElement() const;

    bool existsSuperElement() const { return this->numberOfSuperElement() > 0; }
};

/** @typedef BaseElementaryMatrixPtr */
using BaseElementaryMatrixPtr = std::shared_ptr< BaseElementaryMatrix >;

#endif /* BASEELEMENTARYMATRIX_H_ */
