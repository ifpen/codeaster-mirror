/**
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

#include "LinearAlgebra/BaseAssemblyMatrix.h"

#include "aster_fort_calcul.h"

BaseAssemblyMatrix::BaseAssemblyMatrix( const std::string &name, const std::string &type )
    : DataStructure( name, 19, type ),
      _description( JeveuxVectorChar24( getName() + ".REFA" ) ),
      _scaleFactorLagrangian( JeveuxVectorReal( getName() + ".CONL" ) ),
      _perm( JeveuxVectorLong( getName() + ".PERM" ) ),
      _ccid( JeveuxVectorLong( getName() + ".CCID" ) ),
      _ccll( JeveuxVectorLong( getName() + ".CCLL" ) ),
      _ccii( JeveuxVectorLong( getName() + ".CCII" ) ),
      _isBuilt( false ),
      _isFactorized( false ),
      _dofNum( nullptr ) {};

BaseAssemblyMatrix::BaseAssemblyMatrix( const PhysicalProblemPtr phys_prob,
                                        const std::string &type )
    : BaseAssemblyMatrix( type ) {
    _dofNum = phys_prob->getDOFNumbering();
};

BaseAssemblyMatrix::BaseAssemblyMatrix( const std::string &name, const std::string &type,
                                        const BaseAssemblyMatrix &toCopy )
    : BaseAssemblyMatrix( name, type ) {
    // Jeveux Pointer
    ( *_description ) = ( *toCopy._description );
    ( *_scaleFactorLagrangian ) = ( *toCopy._scaleFactorLagrangian );
    ( *_perm ) = ( *toCopy._perm );
    ( *_ccid ) = ( *toCopy._ccid );
    ( *_ccll ) = ( *toCopy._ccll );
    ( *_ccii ) = ( *toCopy._ccii );
    // Objects
    _dofNum = toCopy._dofNum;
    _isBuilt = toCopy._isBuilt;
    _isFactorized = toCopy._isFactorized;
}

BaseAssemblyMatrix::BaseAssemblyMatrix( BaseAssemblyMatrix &&other )
    : DataStructure( std::move( other ) ) {
    // Jeveux Pointer
    _description = other._description;
    _scaleFactorLagrangian = other._scaleFactorLagrangian;
    _perm = other._perm;
    _ccid = other._ccid;
    _ccll = other._ccll;
    _ccii = other._ccii;
    // Objects
    _dofNum = other._dofNum;
    _isBuilt = other._isBuilt;
    _isFactorized = other._isFactorized;
}

ASTERDOUBLE BaseAssemblyMatrix::getLagrangeScaling() const {
    // Special case of a matrix on a ParallelMesh
    if ( getMesh() && getMesh()->isParallel() ) {
        ASTERDOUBLE scaling0, scaling;
#ifdef ASTER_HAVE_MPI
        CALLO_CONLAG( getName(), &scaling0 );
        scaling = AsterMPI::max( scaling0 );
#endif
        return scaling;
    } else {
        // Other cases
        ASTERDOUBLE scaling = 1.;
        if ( _scaleFactorLagrangian.exists() ) {
            CALLO_CONLAG( getName(), &scaling );
        }
        return scaling;
    }
}

void BaseAssemblyMatrix::updateDOFNumbering() {
    if ( _description.exists() ) {
        _description->updateValuePointer();
        std::string dofName = ( *_description )[1];
        _dofNum = std::make_shared< DOFNumbering >( dofName );
        _isBuilt = true;
    }
};

void BaseAssemblyMatrix::symmetrize() { CALL_MATR_ASSE_SYME( getName() ); };

bool BaseAssemblyMatrix::isMPIFull() const {
    _description->updateValuePointer();
    return strip( ( *_description )[10].toString() ) == "MPI_COMPLET";
};

bool BaseAssemblyMatrix::isSymmetric() const {
    _description->updateValuePointer();
    return strip( ( *_description )[8].toString() ).substr( 0, 2 ) == "MS";
};

std::string BaseAssemblyMatrix::getCalculOption() const {
    _description->updateValuePointer();
    return strip( ( *_description )[3].toString() );
}
