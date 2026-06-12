/**
 * @file PhysicalProblem.cxx
 * @brief Implementation of class PhysicalProblem
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

#include "Studies/PhysicalProblem.h"

#include "aster_pybind.h"

#include "Messages/Messages.h"
#include "Numbering/DOFNumbering.h"
#include "Numbering/ParallelDOFNumbering.h"
#include "Supervis/CommandSyntax.h"
#include "Supervis/Exceptions.h"

PhysicalProblem::PhysicalProblem( const ModelPtr curModel, const MaterialFieldPtr curMat,
                                  const ElementaryCharacteristicsPtr cara )
    : _model( curModel ),
      _mesh( curModel->getMesh() ),
      _materialField( curMat ),
      _elemChara( cara ),
      _virtualSlavCell( nullptr ),
      _virtualCell( nullptr ),
      _listOfLoads( std::make_shared< ListOfLoads >( _model ) ),
      _dofNume( nullptr ),
      _behavProp( nullptr ),
      _externVarRefe( nullptr ) {

    // Add checks
    if ( _elemChara ) {
        if ( _elemChara->containsFieldOnCells() ) {
            if ( _model != _elemChara->getModel() ) {
                UTMESS( "A", "MECANONLINE5_74" );
            }
        } else {
            if ( _model->getMesh() != _elemChara->getModel()->getMesh() ) {
                UTMESS( "A", "MECANONLINE5_75" );
            }
        }
    }

    if ( _materialField ) {
        if ( _mesh != _materialField->getMesh() ) {
            const std::string msg = "Inconsistent meshes: " + _mesh->getName() + " vs " +
                                    _materialField->getMesh()->getName();
            AS_ABORT( msg );
        }
    }
};

PhysicalProblem::PhysicalProblem( const BaseDOFNumberingPtr dofNume )
    : _model( dofNume->getModel() ),
      _mesh( dofNume->getMesh() ),
      _dofNume( dofNume ),
      _listOfLoads( std::make_shared< ListOfLoads >( _model ) ) {};

py::tuple PhysicalProblem::_getState() const {
    return py::make_tuple( _model, _materialField, _elemChara, _dofNume, _behavProp, _externVarRefe,
                           _listOfLoads );
}

PhysicalProblem::PhysicalProblem( const py::tuple &tup )
    : PhysicalProblem( tup[0].cast< ModelPtr >(), tup[1].cast< MaterialFieldPtr >(),
                       tup[2].cast< ElementaryCharacteristicsPtr >() ) {
    if ( tup.size() != 7 ) {
        throw std::runtime_error( "Invalid state!" );
    }
    _dofNume = tup[3].cast< BaseDOFNumberingPtr >();
    _behavProp = tup[4].cast< BehaviourPropertyPtr >();
    _externVarRefe = tup[5].cast< FieldOnCellsRealPtr >();
    _listOfLoads = tup[6].cast< ListOfLoadsPtr >();
}

CodedMaterialPtr PhysicalProblem::getCodedMaterial() const {

    if ( _materialField ) {
        auto codedMater = std::make_shared< CodedMaterial >( _materialField, _model );
        codedMater->allocate( true );
        return codedMater;
    }

    return nullptr;
};

void PhysicalProblem::setListOfLoads( const ListOfLoadsPtr loads ) {
    if ( loads ) {
        _listOfLoads = loads;
    };
};

void PhysicalProblem::setDOFNumbering( const BaseDOFNumberingPtr dofNume ) {
    if ( dofNume ) {
        _dofNume = dofNume;
    }
    computeListOfLoads( "" );
};

void PhysicalProblem::computeBehaviourProperty( py::object &keywords,
                                                const std::string &initialState,
                                                const ASTERINTEGER verbosity ) {
    // Create object for behaviour
    _behavProp = std::make_shared< BehaviourProperty >( _model, _materialField );
    _behavProp->setInitialState( initialState == "OUI" );
    _behavProp->setVerbosity( verbosity > 1 );

    // Check input PyObject
    if ( !PyDict_Check( keywords.ptr() ) && !PyList_Check( keywords.ptr() ) &&
         !PyTuple_Check( keywords.ptr() ) )
        throw std::runtime_error( "Unexpected value for 'COMPORTEMENT'." );

    // Create syntax
    std::string cata;

    if ( _model->isMechanical() ) {
        cata = "code_aster.Cata.Commons.c_comportement.C_COMPORTEMENT_MNL";
    } else if ( _model->isThermal() ) {
        cata = "code_aster.Cata.Commons.c_comportement.C_COMPORTEMENT_TNL";

    } else {
        AS_ABORT( "Should not be here" );
    }

    CommandSyntax cmdSt( cata );
    py::dict kwfact( py::arg( "COMPORTEMENT" ) = keywords );
    cmdSt.define( kwfact );

    // Build objects
    AS_ASSERT( _behavProp->build() );
};

FieldOnCellsRealPtr PhysicalProblem::getExternalStateVariables( const ASTERDOUBLE &time ) const {

    // Create field
    auto FEDesc = getModel()->getFiniteElementDescriptor();
    auto field = std::make_shared< FieldOnCellsReal >( FEDesc );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( getModel()->getName(), 24 );
    std::string materialFieldName = ljust( getMaterialField()->getName(), 24 );
    auto currElemChara = getElementaryCharacteristics();
    std::string elemCharaName( " " );
    if ( currElemChara )
        elemCharaName = currElemChara->getName();
    elemCharaName.resize( 24, ' ' );
    std::string fieldName = ljust( field->getName(), 19 );

    // Output
    std::string out( ' ', 2 );
    std::string base( "G" );

    // Call Fortran WRAPPER
    CALLO_VRCINS_WRAP( modelName, materialFieldName, elemCharaName, &time, fieldName, out, base );

    return field;
}

void PhysicalProblem::computeReferenceExternalStateVariables() {

    // Create field
    auto FEDesc = getModel()->getFiniteElementDescriptor();
    _externVarRefe = std::make_shared< FieldOnCellsReal >( FEDesc );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( getModel()->getName(), 8 );
    std::string materialFieldName = ljust( getMaterialField()->getName(), 8 );
    auto currElemChara = getElementaryCharacteristics();
    std::string elemCharaName( 8, ' ' );
    if ( currElemChara )
        elemCharaName = std::string( currElemChara->getName(), 0, 8 );
    std::string fieldName = ljust( _externVarRefe->getName(), 19 );
    std::string base( "G" );

    // Call Fortran WRAPPER
    CALLO_VRCREF( modelName, materialFieldName, elemCharaName, fieldName, base );
}

bool PhysicalProblem::computeDOFNumbering() {
    // create dofNume
#ifdef ASTER_HAVE_MPI
    if ( getMesh()->isParallel() )
        _dofNume = std::make_shared< ParallelDOFNumbering >();
    else
#endif /* ASTER_HAVE_MPI */
        _dofNume = std::make_shared< DOFNumbering >();

    if ( _virtualSlavCell ) {
        return _dofNume->computeNumbering( getModel(), getListOfLoads(), getVirtualSlavCell() );
    } else {
        return _dofNume->computeNumbering( getModel(), getListOfLoads() );
    }
};

bool PhysicalProblem::computeListOfLoads( std::string command_name ) {
    return _listOfLoads->build( _model, command_name );
};

VectorLong PhysicalProblem::getDirichletBCDOFs( void ) const {
    JeveuxVectorLong ccid( "&&NUME_CCID" );
    std::string base( "V" );

    if ( !_dofNume )
        raiseAsterError( "DOFNumbering not available ; call computeDOFNumering first." );
    // Il faudrait eventuellement rajouter une liste de charge en plus donné par le user
    CALLO_NUMCIMA( getListOfLoads()->getName(), _dofNume->getName(), ccid->getName(), base );

    ccid->updateValuePointer();
    return ccid->toVector();
};

void PhysicalProblem::zeroDirichletBCDOFs( FieldOnNodesReal &field ) const {
    VectorLong dirBC = getDirichletBCDOFs();
    if ( dirBC.size() != field.size() )
        raiseAsterError( "Field has incompatible size" );
    field.updateValuePointers();
    for ( auto i = 0; i < field.size(); i++ ) {
        if ( dirBC[i] == 1 ) {
            field[i] = 0.;
        }
    }
};

bool PhysicalProblem::isMechanical( void ) const { return this->getModel()->isMechanical(); };

bool PhysicalProblem::isThermal( void ) const { return this->getModel()->isThermal(); };

bool PhysicalProblem::isAcoustic( void ) const { return this->getModel()->isAcoustic(); };
