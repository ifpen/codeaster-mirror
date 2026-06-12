#ifndef BASEASSEMBLYMATRIX_H_
#define BASEASSEMBLYMATRIX_H_

/**
 * @file AssemblyMatrix.h
 * @brief Fichier entete de la classe AssemblyMatrix
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

#include "aster_fort_calcul.h"
#include "aster_fort_ds.h"
#include "aster_fort_petsc.h"

#include "DataFields/FieldOnNodes.h"
#include "DataStructures/DataStructure.h"
#include "LinearAlgebra/ElementaryMatrix.h"
#include "Loads/ListOfLoads.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Studies/PhysicalProblem.h"
#include "Supervis/Exceptions.h"
#include "Utilities/Tools.h"

#ifdef ASTER_HAVE_PETSC4PY
#include <petscmat.h>
#endif

#include "Modeling/Model.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/ParallelDOFNumbering.h"

/**
 * @class BaseAssemblyMatrix
 * @brief Classe template definissant la base d'une sd_matr_asse.
 */

class BaseAssemblyMatrix : public DataStructure {
  protected:
    /** @brief Objet Jeveux '.REFA' */
    JeveuxVectorChar24 _description;
    /** @brief Objet '.CONL' */
    JeveuxVectorReal _scaleFactorLagrangian;
    /** @brief Objet Jeveux '.PERM' */
    JeveuxVectorLong _perm;

    /** @brief Objet Jeveux '.CCID' */
    JeveuxVectorLong _ccid;
    /** @brief Objet Jeveux '.CCLL' */
    JeveuxVectorLong _ccll;
    /** @brief Objet Jeveux '.CCII' */
    JeveuxVectorLong _ccii;

    /** @brief Objet nume_ddl */
    BaseDOFNumberingPtr _dofNum;
    /** @brief La matrice a-t-elle été construite ? */
    bool _isBuilt;
    /** @brief La matrice est elle vide ? */
    bool _isFactorized;
    /** @brief Solver name (MUMPS or PETSc) */
    std::string _solverName;

  public:
    /**
     * @typedef BaseAssemblyMatrixPtr
     * @brief Pointeur intelligent vers un BaseAssemblyMatrix
     */
    typedef std::shared_ptr< BaseAssemblyMatrix > BaseAssemblyMatrixPtr;

    /**
     * @brief Constructeur
     */
    BaseAssemblyMatrix() = delete;

    /**
     * @brief Constructeur
     */
    BaseAssemblyMatrix( const std::string &name, const std::string &type );

    /**
     * @brief Constructeur
     */
    BaseAssemblyMatrix( const std::string &name, const std::string &type,
                        const BaseAssemblyMatrix &toCopy );

    /**
     * @brief Constructeur
     */
    BaseAssemblyMatrix( const std::string &type )
        : BaseAssemblyMatrix( ResultNaming::getNewResultName(), type ) {};

    /**
     * @brief Constructeur
     */
    BaseAssemblyMatrix( const PhysicalProblemPtr phys_prob, const std::string &type );

    /**
     * @brief Constructeur
     */

    BaseAssemblyMatrix( BaseAssemblyMatrix &&other );

    /**
     * @brief Make the matrix symmetric
     */
    void symmetrize();

    /**
     * @brief Tell if the matrix is fully filled
     */
    bool isMPIFull();

    /**
     * @brief Tell if the matrix is symmetric
     */
    bool isSymmetric();

    /**
     * @brief Tell the calcul option
     */
    std::string getCalculOption();

    /**
     * @brief Get the size of the matrix
     * @param local [bool] local or global size
     */
    VectorLong size( const bool local = true ) const {
        if ( !_dofNum )
            raiseAsterError( "Sizes not available" );
        ASTERINTEGER shape = _dofNum->getNumberOfDOFs( local );
        return VectorLong( 2, shape );
    };

    /**
     * @brief Get the internal DOFNumbering
     * @return Internal DOFNumbering
     */
    BaseDOFNumberingPtr getDOFNumbering() const { return _dofNum; };

    /**
     * @brief Get mesh
     * @return Internal mesh
     */
    BaseMeshPtr getMesh() const {
        if ( _dofNum ) {
            return _dofNum->getMesh();
        }
        return nullptr;
    };

    /**
     * @brief Set new values
     */
    virtual void setValues( const VectorLong &idx, const VectorLong &jdx,
                            const VectorReal &values ) {
        // Template class raises error. It must be specialized in each instanciated class.
        throw std::runtime_error( "Not implemented" );
    };

    /**
     * @brief Scale the matrix using right and left multiplication
     * by diagonal matrices stored as vectors
     */
    virtual void scale( const VectorReal &lvect, const VectorReal &rvect ) {
        // Template class raises error. It must be specialized in each instanciated class.
        throw std::runtime_error( "Not implemented" );
    };

    /**
     * @brief Transpose
     */
    void transpose() { CALLO_MATR_ASSE_TRANSPOSE( getName() ); };

    /**
     * @brief Transpose and conjugate
     */
    void transposeConjugate() { CALLO_MATR_ASSE_TRANSPOSE_CONJUGATE( getName() ); };

    /**
     * @brief Print the matrix in code_aster or matlab format in given logical unit
     */
    void print( const std::string format = "ASTER", const ASTERINTEGER unit = 6 ) const {
        CALLO_MATR_ASSE_PRINT( getName(), &unit, format );
    };

#ifdef ASTER_HAVE_PETSC4PY
    /**
     * @brief Conversion to petsc4py
     * @return converted matrix
     */
    Mat toPetsc( const bool &local ) {
        Mat myMat;
        PetscErrorCode ierr;
        const ASTERINTEGER local2 = local ? 1 : 0;

        if ( !_isBuilt )
            throw std::runtime_error( "Assembly matrix is empty" );
        if ( getType() != "MATR_ASSE_DEPL_R" && getType() != "MATR_ASSE_TEMP_R" &&
             getType() != "MATR_ASSE_ELIM_R" )
            throw std::runtime_error( "Not yet implemented" );
        if ( not getMesh()->isParallel() && local )
            throw std::runtime_error( "local export is only usable on a ParallelMesh" );

        CALLO_MATASS2PETSC( getName(), &local2, &myMat, &ierr );

        return myMat;
    };
#endif

    /**
     * @brief Methode permettant de savoir si la matrice a été construite
     * @return true si construite
     */
    bool isBuilt() const { return _isBuilt; };

    /**
     * @brief Methode permettant de savoir si la matrice est factorisée
     * @return true si factorisée
     */
    bool isFactorized() const { return _isFactorized; };

    void setFactorized( const bool &facto ) { _isFactorized = facto; };

    /**
     * @brief Methode permettant de definir la numerotation
     * @param currentNum objet DOFNumbering
     */
    void setDOFNumbering( const BaseDOFNumberingPtr currentNum ) {
        _dofNum = currentNum;
        _isBuilt = true;
    };

    /** @brief update _dofNum using DOFNumbering name created in Fortran */
    void updateDOFNumbering();

    /**
     * @brief Function to set the solver name (MUMS or PETSc)
     * @param sName name of solver ("MUMPS" or "PETSC")
     * @todo delete this function and the attribute _solverName
     */
    void setSolverName( const std::string &sName ) { _solverName = sName; };

    /**
     * @brief Delete the factorized matrix used by MUMPS or PETSc if it exist
     * @param sName name of solver ("MUMPS" or "PETSC")
     * @todo delete this function and the attribute _solverName
     */
    bool deleteFactorizedMatrix( void ) {
        if ( _description.exists() && get_sh_jeveux_status() == 1 ) {
            CALLO_DELETE_MATRIX( getName(), _solverName );
        }

        _isFactorized = false;

        return true;
    };

    /**
     * @brief Return True if CCID object exists for DirichletElimination
     */
    bool hasDirichletEliminationDOFs() const { return _ccid.exists(); }

    /**
     * @brief Return CCID object if exists for DirichletElimination
     */
    VectorLong getDirichletBCDOFs( void ) const {
        if ( _ccid.exists() ) {
            _ccid->updateValuePointer();
            return _ccid->toVector();
        } else {
            ASTERINTEGER shape = _dofNum->getNumberOfDOFs( true );
            return VectorLong( shape + 1, 0 );
        }
    }

    /**
     * @brief Zero the Dirichlet BC DOFs in the given field
     */
    void zeroDirichletBCDOFs( FieldOnNodesReal &field ) const {
        if ( !_ccid.exists() )
            return; // No Dirichlet BCs to apply
        VectorLong dirBC = getDirichletBCDOFs();
        if ( dirBC.size() != field.size() + 1 )
            raiseAsterError( "Field has incompatible size" );
        field.updateValuePointers();
        for ( auto i = 0; i < field.size(); i++ ) {
            if ( dirBC[i] == 1 ) {
                field[i] = 0.;
            }
        }
    }

    /**
     * @brief Return the scaling factor of Lagrange multipliers
     */
    ASTERDOUBLE getLagrangeScaling() const;

    virtual BaseAssemblyMatrixPtr getEmptyMatrix( const std::string &name ) const {
        AS_ABORT( "Not allowed" );
    }
};

typedef std::shared_ptr< BaseAssemblyMatrix > BaseAssemblyMatrixPtr;

#endif /* BASEASSEMBLYMATRIX_H_ */
