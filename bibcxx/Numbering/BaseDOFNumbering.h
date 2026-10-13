/**
 * @file BaseDOFNumbering.h
 * @brief Fichier entete de la classe BaseDOFNumbering
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
#include "LinearAlgebra/ElementaryMatrix.h"
#include "LinearAlgebra/MatrixStorage.h"
#include "Loads/DirichletBC.h"
#include "Loads/ListOfLoads.h"
#include "Loads/MechanicalLoad.h"
#include "MemoryManager/JeveuxVector.h"
#include "Meshes/BaseMesh.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Numbering/EquationNumbering.h"

#pragma once

// Forward declaration
template < typename ValueType, PhysicalQuantityEnum PhysicalQuantity >
class ElementaryMatrix;

class DOFNumbering;
class ParallelDOFNumbering;

using ElementaryMatrixDisplacementRealPtr =
    std::shared_ptr< ElementaryMatrix< ASTERDOUBLE, Displacement > >;
using ElementaryMatrixDisplacementComplexPtr =
    std::shared_ptr< ElementaryMatrix< ASTERCOMPLEX, Displacement > >;
using ElementaryMatrixTemperatureRealPtr =
    std::shared_ptr< ElementaryMatrix< ASTERDOUBLE, Temperature > >;
using ElementaryMatrixPressureComplexPtr =
    std::shared_ptr< ElementaryMatrix< ASTERCOMPLEX, Pressure > >;

/**
 * @class BaseDOFNumbering
 * @brief Class definissant un nume_ddl
 * @author Nicolas Sellenet
 */
class BaseDOFNumbering : public DataStructure {
  public:
    typedef std::variant< ElementaryMatrixDisplacementRealPtr,
                          ElementaryMatrixDisplacementComplexPtr,
                          ElementaryMatrixTemperatureRealPtr, ElementaryMatrixPressureComplexPtr >
        MatrElem;

  private:
    class ElementaryMatrixGetModel {
      public:
        template < typename T >
        ModelPtr operator()( const T &operand ) const {
            return operand->getModel();
        };
    };

    class ElementaryMatrixGetName {
      public:
        template < typename T >
        std::string operator()( const T &operand ) const {
            return operand->getName();
        };
    };

    class ElementaryMatrixGetFEDescrp {
      public:
        template < typename T >
        std::vector< FiniteElementDescriptorPtr > operator()( const T &operand ) const {
            return operand->getFiniteElementDescriptors();
        };
    };

    class ElementaryMatrixGetMesh {
      public:
        template < typename T >
        BaseMeshPtr operator()( const T &operand ) const {
            return operand->getMesh();
        };
    };

  private:
    class MultFrontGarbage {
        /** @brief Objet Jeveux '.ADNT' */
        JeveuxVectorShort _adnt;
        /** @brief Objet Jeveux '.GLOB' */
        JeveuxVectorShort _glob;
        /** @brief Objet Jeveux '.LOCL' */
        JeveuxVectorShort _locl;
        /** @brief Objet Jeveux '.PNTI' */
        JeveuxVectorShort _pnti;
        /** @brief Objet Jeveux '.RENU' */
        JeveuxVectorChar8 _renu;
        /** @brief Objet Jeveux '.ADPI' */
        JeveuxVectorLong _adpi;
        /** @brief Objet Jeveux '.ADRE' */
        JeveuxVectorLong _adre;
        /** @brief Objet Jeveux '.ANCI' */
        JeveuxVectorLong _anci;
        /** @brief Objet Jeveux '.LFRN' */
        JeveuxVectorLong _debf;
        /** @brief Objet Jeveux '.DECA' */
        JeveuxVectorLong _deca;
        /** @brief Objet Jeveux '.DEFS' */
        JeveuxVectorLong _defs;
        /** @brief Objet Jeveux '.DESC' */
        JeveuxVectorLong _desc;
        /** @brief Objet Jeveux '.DIAG' */
        JeveuxVectorLong _diag;
        /** @brief Objet Jeveux '.FILS' */
        JeveuxVectorLong _fils;
        /** @brief Objet Jeveux '.FRER' */
        JeveuxVectorLong _frer;
        /** @brief Objet Jeveux '.LGBL' */
        JeveuxVectorLong _lgbl;
        /** @brief Objet Jeveux '.LGSN' */
        JeveuxVectorLong _lgsn;
        /** @brief Objet Jeveux '.NBAS' */
        JeveuxVectorLong _nbas;
        /** @brief Objet Jeveux '.NBLI' */
        JeveuxVectorLong _nbli;
        /** @brief Objet Jeveux '.NCBL' */
        JeveuxVectorLong _nbcl;
        /** @brief Objet Jeveux '.NOUV' */
        JeveuxVectorLong _nouv;
        /** @brief Objet Jeveux '.PARE' */
        JeveuxVectorLong _pare;
        /** @brief Objet Jeveux '.SEQU' */
        JeveuxVectorLong _sequ;
        /** @brief Objet Jeveux '.SUPN' */
        JeveuxVectorLong _supn;

        MultFrontGarbage( const std::string &DOFNumName )
            : _adnt( DOFNumName + ".ADNT" ),
              _glob( DOFNumName + ".GLOB" ),
              _locl( DOFNumName + ".LOCL" ),
              _pnti( DOFNumName + ".PNTI" ),
              _renu( DOFNumName + ".RENU" ),
              _adpi( DOFNumName + ".ADPI" ),
              _adre( DOFNumName + ".ADRE" ),
              _anci( DOFNumName + ".ANCI" ),
              _debf( DOFNumName + ".LFRN" ),
              _deca( DOFNumName + ".DECA" ),
              _defs( DOFNumName + ".DEFS" ),
              _desc( DOFNumName + ".DESC" ),
              _diag( DOFNumName + ".DIAG" ),
              _fils( DOFNumName + ".FILS" ),
              _frer( DOFNumName + ".FRER" ),
              _lgbl( DOFNumName + ".LGBL" ),
              _lgsn( DOFNumName + ".LGSN" ),
              _nbas( DOFNumName + ".NBAS" ),
              _nbli( DOFNumName + ".NBLI" ),
              _nbcl( DOFNumName + ".NCBL" ),
              _nouv( DOFNumName + ".NOUV" ),
              _pare( DOFNumName + ".PARE" ),
              _sequ( DOFNumName + ".SEQU" ),
              _supn( DOFNumName + ".SUPN" ) {};
        friend class BaseDOFNumbering;
    };
    typedef std::shared_ptr< MultFrontGarbage > MultFrontGarbagePtr;

  protected:
    class LocalEquationNumbering : public DataStructure {
      protected:
        /** @brief Objet Jeveux '.NEQU' */
        JeveuxVectorLong _numberOfEquations;
        /** @brief Objet Jeveux '.DELG' */
        JeveuxVectorLong _lagrangianInformations;
        /** @brief Objet Jeveux '.PRNO' */
        JeveuxCollectionLong _componentsOnNodes;
        /** @brief Objet Jeveux '.NUEQ' */
        JeveuxVectorLong _indexationVector;
        /** @brief Objet Jeveux '.NULG' */
        JeveuxVectorLong _localToGlobal;
        /** @brief Objet Jeveux '.NUGL' */
        JeveuxVectorLong _globalToLocal;
        /*protected* @brief Objet Jeveux '.PDDL' */
        JeveuxVectorLong _localToRank;

      public:
        LocalEquationNumbering( const std::string &baseName )
            : DataStructure( baseName + ".NUML", 19, "NUML_EQUA" ),
              _numberOfEquations( getName() + ".NEQU" ),
              _lagrangianInformations( getName() + ".DELG" ),
              _componentsOnNodes( getName() + ".PRNO" ),
              _indexationVector( getName() + ".NUEQ" ),
              _localToGlobal( getName() + ".NULG" ),
              _globalToLocal( getName() + ".NUGL" ),
              _localToRank( getName() + ".PDDL" ) {};

        /**
         * @brief Returns the vector of local to global numbering
         */
        const JeveuxVectorLong getLocalToGlobal() const { return _localToGlobal; }
        /**
         * @brief Returns the vector of global to local numbering
         */
        const JeveuxVectorLong getGlobalToLocal() const { return _globalToLocal; }
        /**
         * @brief Returns the vector of the rank owning the local dof number
         */
        const JeveuxVectorLong getLocalToRank() const { return _localToRank; }

        bool exists() const {
            return _numberOfEquations->exists() && _lagrangianInformations->exists() &&
                   _componentsOnNodes->exists();
        }

        friend class BaseDOFNumbering;
        friend class DOFNumbering;
        friend class ParallelDOFNumbering;
    };
    typedef std::shared_ptr< LocalEquationNumbering > LocalEquationNumberingPtr;

  protected:
    /** @brief Objet Jeveux '.NSLV' */
    JeveuxVectorChar24 _nameOfSolverDataStructure;
    /** @brief Objet Jeveux '.SMOS' */
    MorseStoragePtr _smos;
    /** @brief Objet Jeveux '.SLCS' */
    LigneDeCielPtr _slcs;
    /** @brief Objet Jeveux '.MLTF' */
    MultFrontGarbagePtr _mltf;

    /** @brief Objet '.NUML' */
    LocalEquationNumberingPtr _localNumbering;

    /** @brief Vectors of FiniteElementDescriptor */
    std::vector< FiniteElementDescriptorPtr > _FEDVector;
    std::set< std::string > _FEDNames;

    /**
     * @brief Add a FiniteElementDescriptor to elementary matrix
     * @param FiniteElementDescriptorPtr FiniteElementDescriptor
     */
    bool addFiniteElementDescriptor( const FiniteElementDescriptorPtr &curFED );

    bool addFiniteElementDescriptor( const std::vector< FiniteElementDescriptorPtr > &curFED );

  protected:
    /**
     * @brief Constructeur
     * @param name nom souhaité de la sd (utile pour le BaseDOFNumbering d'une sd_resu)
     */

    BaseDOFNumbering( const std::string name, const std::string &type );

  public:
    virtual ~BaseDOFNumbering() {};

    /**
     * @typedef BaseDOFNumberingPtr
     * @brief Pointeur intelligent vers un BaseDOFNumbering
     */
    typedef std::shared_ptr< BaseDOFNumbering > BaseDOFNumberingPtr;

    /**
     * @brief Returns the EquationNumberingPtr
     */
    virtual EquationNumberingPtr getEquationNumbering() const = 0;

    virtual std::string getPhysicalQuantity() const = 0;

    /**
     * @brief Build the Numbering of DOFs
     */
    virtual bool computeNumbering( const ModelPtr model, const ListOfLoadsPtr listOfLoads,
                                   bool verbose = true );

    virtual bool computeNumbering( const ModelPtr model, const ListOfLoadsPtr listOfLoads,
                                   const FiniteElementDescriptorPtr defiCont, bool verbose = true );

    /**
     * @brief Build the Numbering of DOFs
     */
    virtual bool computeNumbering( const std::vector< MatrElem > matrix, bool verbose = true );

    /**
     * @brief renumbering of DOFs
     */
    virtual bool computeRenumbering( const ModelPtr model, const ListOfLoadsPtr listOfLoads,
                                     const FiniteElementDescriptorPtr defiCont,
                                     const FiniteElementDescriptorPtr virtContElem,
                                     bool verbose = true );

    /**
     * @brief Build the Numbering of DOFs
     */
    virtual bool computeNumbering( const std::vector< FiniteElementDescriptorPtr > &Feds,
                                   const std::string &locatedComponent, bool verbose = true );

    /**
     * @brief Are Lagrange Multipliers used for BC or MPC
     */
    virtual bool useLagrangeDOF() const { AS_ABORT( "Not allowed" ); };

    /**
     * @brief Are Single Lagrange Multipliers used for BC or MPC
     */
    virtual bool useSingleLagrangeDOF() const { AS_ABORT( "Not allowed" ); };

    /**
     * @brief Get The Component Associated To A Given Row
     */
    virtual std::string getComponentFromDOF( const ASTERINTEGER dof,
                                             const bool local = false ) const {
        AS_ABORT( "Not allowed" );
    };

    /**
     * @brief Get The Components Associated To A Given Node
     */
    virtual VectorString getComponentFromNode( const ASTERINTEGER node,
                                               const bool local = false ) const {
        AS_ABORT( "Not allowed" );
    };

    /**
     * @brief Get The Node Id Associated To A Given Row
     */
    virtual ASTERINTEGER getNodeFromDOF( const ASTERINTEGER dof, const bool local = false ) const {
        AS_ABORT( "Not allowed" );
    };

    /**
     * @brief Return true if a physical dof is Associated To A Given Row
     */
    virtual bool isPhysicalDOF( const ASTERINTEGER dof, const bool local = false ) const {
        AS_ABORT( "Not allowed" );
    };

    /**
     * @brief get the Row index Associated To the Component of a Node
     */
    virtual ASTERINTEGER getDOFFromNodeAndComponent( const ASTERINTEGER &node,
                                                     const std::string &comp,
                                                     const bool local = false ) const {
        AS_ABORT( "Not allowed" );
    }

    /**
     * @brief Get The total number of Dofs
     */
    virtual ASTERINTEGER getNumberOfDOFs( const bool local = false ) const {
        AS_ABORT( "Not allowed" );
    };

    /**
     * @brief Get Rows Associated to all Physical Dof
     */
    virtual VectorLong getPhysicalDOFs( const bool local = false ) const {
        AS_ABORT( "Not allowed" );
    }

    /**
     * @brief Get Rows Associated to Lagrange Multipliers Dof
     */
    virtual VectorLong getLagrangeDOFs( const bool local = false ) const {
        AS_ABORT( "Not allowed" );
    }

    /**
     * @brief Get Rows Associated to first and second Lagrange Multipliers Dof
     */

    virtual std::map< ASTERINTEGER, VectorLong >
    getDictOfLagrangeDOFs( const bool local = false ) const {
        AS_ABORT( "Not allowed" );
    }

    /**
     * @brief Get Assigned Components
     */
    virtual VectorString getComponents() const { AS_ABORT( "Not allowed" ); }

    /**
     * @brief Get all FiniteElementDescriptors
     * @return vector of all FiniteElementDescriptors
     */
    std::vector< FiniteElementDescriptorPtr > getFiniteElementDescriptors() { return _FEDVector; };

    /**
     * @brief Assign FiniteElementDescriptors
     */
    void setFiniteElementDescriptors( std::vector< FiniteElementDescriptorPtr > &descr ) {
        _FEDVector = descr;
    };

    /**
     * @brief Get mesh
     * @return Internal mesh
     */
    virtual BaseMeshPtr getMesh() const = 0;

    virtual void setMesh( const BaseMeshPtr mesh ) const = 0;

    /**
     * @brief Methode permettant de savoir si la numerotation est vide
     * @return true si la numerotation est vide
     */

    virtual bool exists() const { AS_ABORT( "Not allowed" ); };

    /**
     * @brief Methode permettant de savoir si l'objet est parallel
     * @return false
     */
    virtual bool isParallel() { return false; };

    /**
     * @brief Get model
     */
    virtual ModelPtr getModel() const = 0;

    /**
     * @brief Set model
     */
    virtual void setModel( const ModelPtr & ) = 0;

    MorseStoragePtr getMorseStorage() const { return _smos; };
};

/**
 * @typedef BaseDOFNumberingPtr
 * @brief Enveloppe d'un pointeur intelligent vers un BaseDOFNumbering
 * @author Nicolas Sellenet
 */
typedef std::shared_ptr< BaseDOFNumbering > BaseDOFNumberingPtr;
