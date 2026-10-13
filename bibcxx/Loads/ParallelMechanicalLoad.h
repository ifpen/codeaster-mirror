
#include "astercxx.h"

#ifdef ASTER_HAVE_MPI

#ifndef PARALLELMECHANICALLOAD_H_
#define PARALLELMECHANICALLOAD_H_

/**
 * @file ParallelMechanicalLoad.h
 * @brief Fichier entete de la classe ParallelMechanicalLoad
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

#include "Coupling/CouplingPairingData.h"
#include "DataFields/ConstantFieldOnCells.h"
#include "DataStructures/DataStructure.h"
#include "Loads/MechanicalLoad.h"
#include "Meshes/ConnectionMesh.h"
#include "Meshes/MeshExplorer.h"
#include "Modeling/FiniteElementDescriptor.h"
#include "Modeling/Model.h"
#include "Modeling/ParallelContactFEDescriptor.h"
#include "Modeling/ParallelFiniteElementDescriptor.h"
#include "ParallelUtilities/AsterMPI.h"
#include "Supervis/ResultNaming.h"
#include "Utilities/SyntaxSaver.h"

#include <fstream>

/**
 * @class ParallelMechanicalLoad
 * @brief Classe definissant une charge dualisée parallèle
 */
template < typename ConstantFieldOnCellsType >
class ParallelMechanicalLoad : public DataStructure {
  public:
    using ConstantFieldOnCellsTypePtr = std::shared_ptr< ConstantFieldOnCellsType >;

    using ParallelMechanicalLoadPtr = std::shared_ptr< ParallelMechanicalLoad >;

  private:
    template < typename ConstantFieldOnCellsType2Ptr >
    void transferConstantFieldOnCells( const ConstantFieldOnCellsType2Ptr &fieldIn,
                                       ConstantFieldOnCellsType2Ptr &fieldOut ) {
        const auto &toKeep = _FEDesc->getVirtualCellsToKeep();

        std::string savedName( "" );
        fieldIn->build();
        fieldOut->allocate( fieldIn );
        const auto sizeFieldIn = ( *fieldIn ).size();
        const auto vect_resu = fieldIn->getValues();

        fieldIn->updateValuePointers();
        for ( int pos = 0; pos < sizeFieldIn; ++pos ) {
            const auto &zone = fieldIn->getZoneDescription( pos );
            const auto &curFEDesc = zone.getFiniteElementDescriptor();
            if ( curFEDesc->getName() != savedName && savedName != "" ) {
                std::string a( "Different FiniteElementDescriptor in one ConstantFieldOnCells is "
                               "not allowed" );
                throw std::runtime_error( a );
            }
            savedName = curFEDesc->getName();

            const auto &listCells = zone.getListOfCells();
            VectorLong toCopy;
            toCopy.reserve( listCells.size() );
            for ( const auto &num : listCells ) {
                if ( toKeep[-num - 1] != 1 )
                    toCopy.push_back( toKeep[-num - 1] );
            }

            if ( toCopy.size() != 0 ) {
                const auto newZone =
                    ConstantFieldOnZone( zone.getFiniteElementDescriptor(), toCopy );
                const auto resu = vect_resu[pos];
                fieldOut->setValueOnZone( newZone, resu );
            }
        }

        fieldOut->build();
    };

    void allocateType( const MechanicalLoadPtr< ConstantFieldOnCellsType > &load ) {
        auto typeLoadStd = load->getType();
        typeLoadStd->updateValuePointer();

        _type->allocate( 1 );
        ( *_type )[0] = ( *typeLoadStd )[0];
        _modelName->allocate( 1 );
        ( *_modelName )[0] = _model->getName();
    }

    void allocateFields( const MechanicalLoadPtr< ConstantFieldOnCellsType > &load ) {

        transferConstantFieldOnCells( load->getMechanicalLoadDescription()->getImposedField(),
                                      _cimpo );
        transferConstantFieldOnCells(
            load->getMechanicalLoadDescription()->getMultiplicativeField(), _cmult );
    };

    void buildParallelCouplingFED( const MechanicalLoadPtr< ConstantFieldOnCellsType > &load ) {
        AS_ASSERT( load->getPairingField() && load->getPairingField()->exists() );
        _pairingPts = std::make_shared< FieldOnCellsReal >( getName() + ".CHME.PAIRS",
                                                            *( load->getPairingField() ) );
        const auto hhoProj = load->getHHOProjector();
        if ( hhoProj && hhoProj->exists() ) {
            _hhoProj = std::make_shared< FieldOnCellsReal >( getName() + ".CHME.MPHHO", *hhoProj );
        }

        const auto cMesh =
            std::dynamic_pointer_cast< ConnectionMesh >( load->getModel()->getMesh() );

        if ( cMesh ) {
            // as ParallelContactPairing - same objects.
            _FEDesc = nullptr;
            _cplFED = std::make_shared< ParallelContactFEDescriptor >(
                getName() + ".CHME.LIGRE", _pairingPts->getDescription(), cMesh, _model );
        } else {
            throw std::runtime_error( "Problem with ConnectionMesh while coupling pairing" );
        }
    };

    void buildParallelFED( const MechanicalLoadPtr< ConstantFieldOnCellsType > &load ) {
        _FEDesc = std::make_shared< ParallelFiniteElementDescriptor >(
            getName() + ".CHME.LIGRE", load->getFiniteElementDescriptor(),
            load->getModel()->getConnectionMesh(), _model );
        _cimpo = std::make_shared< ConstantFieldOnCellsType >( getName() + ".CHME.CIMPO", _FEDesc );
        _cmult = std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CHME.CMULT", _FEDesc );
    }

  protected:
    /** @brief Modele */
    ModelPtr _model;
    /** @brief Vecteur Jeveux '.LIGRE' */
    ParallelFiniteElementDescriptorPtr _FEDesc;
    /** @brief Carte '.CIMPO' */
    ConstantFieldOnCellsTypePtr _cimpo;
    /** @brief Carte '.CMULT' */
    ConstantFieldOnCellsRealPtr _cmult;
    /** @brief Vecteur Jeveux '.TYPE' */
    JeveuxVectorChar8 _type;
    /** @brief Vecteur Jeveux '.MODEL.NOMO' */
    JeveuxVectorChar8 _modelName;
    /** @brief Syntax used to build object */
    SyntaxSaverPtr _syntax;
    /** @brief node and cell support groups of ParallelMechnicalLoad */
    VectorString _grpNo, _grpMa;
    /** @brief Field '.PAIRS' for LIAISON_MASSIF */
    FieldOnCellsRealPtr _pairingPts;
    /** @brief Field '.MPHHO' for LIAISON_MASSIF */
    FieldOnCellsRealPtr _hhoProj;
    /** @brief LIGREL for LIAISON_MASSIF - Numbering and Compute */
    ParallelContactFEDescriptorPtr _cplFED;

  public:
    /**
     * @brief Constructeur
     */
    ParallelMechanicalLoad( void ) = delete;

    /**
     *
     * @brief Constructeur
     */
    ParallelMechanicalLoad( const MechanicalLoadPtr< ConstantFieldOnCellsType > &load,
                            const ModelPtr &model )
        : ParallelMechanicalLoad( ResultNaming::getNewResultName(), load, model ) {};

    /** @brief Constructor */
    ParallelMechanicalLoad( const ParallelMechanicalLoadPtr load, const ModelPtr &model )
        : DataStructure( ResultNaming::getNewResultName(), 8, "CHAR_MECA" ),
          _type( getName() + ".TYPE" ),
          _modelName( getName() + ".CHME.MODEL.NOMO" ),
          _model( model ) {
        _syntax = load->_syntax;
        _grpNo = load->_grpNo;
        _grpMa = load->_grpMa;
        auto pMesh = std::static_pointer_cast< ParallelMesh >( model->getMesh() );
        auto connectionMesh = std::make_shared< ConnectionMesh >( pMesh, _grpNo, _grpMa );
        auto connectionModel = std::make_shared< Model >( connectionMesh );
        connectionModel->setFrom( model );
        auto partialMechanicalLoad =
            std::make_shared< MechanicalLoad< ConstantFieldOnCellsType > >( connectionModel );
        partialMechanicalLoad->buildFromSyntax( _syntax );

        buildParallelFED( partialMechanicalLoad );
        allocateType( partialMechanicalLoad );
        allocateFields( partialMechanicalLoad );
    };

    /**
     * @brief Constructeur
     */
    ParallelMechanicalLoad( const std::string &name,
                            const MechanicalLoadPtr< ConstantFieldOnCellsType > &load,
                            const ModelPtr &model )
        : DataStructure( name, 8, "CHAR_MECA" ),
          _type( getName() + ".TYPE" ),
          _modelName( getName() + ".CHME.MODEL.NOMO" ),
          _model( model ) {

        allocateType( load );

        const auto pairing = load->getPairingField();

        if ( pairing && pairing->exists() ) {
            // Can not mixed LIAISON_MASSIF with loads
            AS_ASSERT( !load->getMechanicalLoadDescription()->getImposedField()->exists() );
            AS_ASSERT( !load->getMechanicalLoadDescription()->getMultiplicativeField()->exists() );

            buildParallelCouplingFED( load );
        } else {
            buildParallelFED( load );
            allocateFields( load );
        }
    };

    /**
     * @brief Constructeur
     */
    ParallelMechanicalLoad( const std::string &name,
                            const ParallelFiniteElementDescriptorPtr &fEDesc,
                            const ModelPtr &model )
        : DataStructure( name, 8, "CHAR_MECA" ),
          _FEDesc( fEDesc ),
          _cimpo(
              std::make_shared< ConstantFieldOnCellsType >( getName() + ".CHME.CIMPO", _FEDesc ) ),
          _cmult(
              std::make_shared< ConstantFieldOnCellsReal >( getName() + ".CHME.CMULT", _FEDesc ) ),
          _type( getName() + ".TYPE" ),
          _modelName( getName() + ".CHME.MODEL.NOMO" ),
          _model( model ) {};

    /**
     * @brief Function membre debugPrint
     * @param logicalUnit Unite logique d'impression
     */
    void debugPrint( const int logicalUnit = 6 ) const {
        this->DataStructure::debugPrint( logicalUnit );
        const auto &explorer = _FEDesc->getVirtualCellsExplorer();
        const auto &mesh = _FEDesc->getMesh();
        auto &LToGmapMesh = mesh->getLocalToGlobalNodeIds();
        auto &LToGmapFE = _FEDesc->getLocalToGlobalMapping();
        auto e1 = LToGmapMesh.exists();
        auto e2 = LToGmapFE.exists();
        if ( !e1 || !e2 )
            return;
        LToGmapMesh->updateValuePointer();
        LToGmapFE->updateValuePointer();

        std::ofstream outFile;
        outFile.open( "fort." + std::to_string( logicalUnit ), std::ios::app );
        outFile << "\nEcriture de la connectivité des mailles fantômes en numérotation globale\n";
        for ( const auto meshElem : explorer ) {
            const auto &numElem = meshElem.getId() + 1;
            outFile << numElem << " : ";
            for ( auto numNode : meshElem ) {
                if ( numNode > 0 ) {
                    outFile << ( *LToGmapMesh )[numNode - 1] << " ";
                } else {
                    outFile << ( *LToGmapFE )[-numNode - 1] << " ";
                }
            }
            outFile << std::endl;
        }
        outFile.close();
    };

    /**
     * @brief Get the finite element descriptor
     */
    FiniteElementDescriptorPtr getFiniteElementDescriptor() const {
        if ( _cplFED ) {
            return _cplFED;
        }

        return _FEDesc;
    };

    /**
     * @brief Get the model
     */
    ModelPtr getModel() const {
        if ( ( !_model ) || _model->isEmpty() )
            throw std::runtime_error( "Model of current load is empty" );
        return _model;
    };

    ConstantFieldOnCellsRealPtr getMultiplicativeField() const { return _cmult; };

    ConstantFieldOnCellsTypePtr getImposedField() const { return _cimpo; };

    /**
     * @brief Function to set rebuild parameters (used for balancing)
     * @param syntax aster syntax
     * @param grpNo cell groups
     * @param syntax node groups
     */
    void setRebuildParameters( SyntaxSaverPtr syntax, const VectorString &grpNo,
                               const VectorString &grpMa ) {
        _syntax = syntax;
        _grpNo = grpNo;
        _grpMa = grpMa;
    };

    auto getPairingField() const { return _pairingPts; }
    auto getHHOProjector() const { return _hhoProj; }

    auto getParallelCouplingFEDescriptor() const { return _cplFED; }
};

/**
 * @typedef ParallelMechanicalLoadPtr
 * @brief Pointeur intelligent vers un ParallelMechanicalLoad
 */
using ParallelMechanicalLoadReal = ParallelMechanicalLoad< ConstantFieldOnCellsReal >;

using ParallelMechanicalLoadFunction = ParallelMechanicalLoad< ConstantFieldOnCellsChar24 >;

using ParallelMechanicalLoadRealPtr = std::shared_ptr< ParallelMechanicalLoadReal >;

using ParallelMechanicalLoadFunctionPtr = std::shared_ptr< ParallelMechanicalLoadFunction >;

/** @typedef std::list de ParallelMechanicalLoad */
using ListParaMecaLoadReal = std::list< ParallelMechanicalLoadRealPtr >;
/** @typedef Iterateur sur une std::list de ParallelMechanicalLoad */
using ListParaMecaLoadRealIter = ListParaMecaLoadReal::iterator;
/** @typedef Iterateur constant sur une std::list de ParallelMechanicalLoad */
using ListParaMecaLoadRealCIter = ListParaMecaLoadReal::const_iterator;

/** @typedef std::list de ParallelMechanicalLoad */
using ListParaMecaLoadFunction = std::list< ParallelMechanicalLoadFunctionPtr >;
/** @typedef Iterateur sur une std::list de ParallelMechanicalLoad */
using ListParaMecaLoadFunctionIter = ListParaMecaLoadFunction::iterator;
/** @typedef Iterateur constant sur une std::list de ParallelMechanicalLoad */
using ListParaMecaLoadFunctionCIter = ListParaMecaLoadFunction::const_iterator;

#endif /* PARALLELMECHANICALLOAD_H_ */

#endif /* ASTER_HAVE_MPI */
