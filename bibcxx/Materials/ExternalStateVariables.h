#ifndef EXTERNALSTATEVARIABLES_H_
#define EXTERNALSTATEVARIABLES_H_

/**
 * @file ExternalStateVariable.h
 * @brief Header of ExternalStateVariables
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

#include "DataFields/DataField.h"
#include "Functions/Formula.h"
#include "Functions/Function.h"
#include "Meshes/BaseMesh.h"
#include "Meshes/Mesh.h"
#include "Meshes/MeshEntities.h"
#include "Meshes/ParallelMesh.h"
#include "Meshes/Skeleton.h"

class TransientResult;
using TransientResultPtr = std::shared_ptr< TransientResult >;

/**
 * @class EvolutionParameter
 * @brief Class EvolutionParameter to be used when ExternalStateVariable is time dependant
 */
class EvolutionParameter {
  private:
    TransientResultPtr _transientResult;
    std::string _fieldName;
    std::string _leftExtension;
    std::string _rightExtension;
    FunctionPtr _timeFunction;
    FormulaPtr _timeFormula;

  public:
    EvolutionParameter( const TransientResultPtr &result, const std::string fieldName );

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    // EvolutionParameter( const py::tuple &tup );
    // py::tuple _getState() const;

    std::string getFieldName() const { return _fieldName; };

    std::string getLeftExtension() const { return _leftExtension; };

    std::string getRightExtension() const { return _rightExtension; };

    TransientResultPtr getTransientResult() const { return _transientResult; };

    FormulaPtr getTimeFormula() const { return _timeFormula; };

    FunctionPtr getTimeFunction() const { return _timeFunction; };

    void setLeftExtension( const std::string typeExtension );

    void setRightExtension( const std::string typeExtension );

    void setTimeFunction( const FormulaPtr &func ) {
        _timeFormula = func;
        _timeFunction = nullptr;
    };

    void setTimeFunction( const FunctionPtr &func ) {
        _timeFormula = nullptr;
        _timeFunction = func;
    };
};

using EvolutionParameterPtr = std::shared_ptr< EvolutionParameter >;

/** @brief Enum to identify external state variable */
enum class externVarEnumInt : int {
    Unknown = -1,
    Temperature,
    Geometry,
    Corrosion,
    IrreversibleStrain,
    ConcreteHydration,
    Irradiation,
    SteelPhases,
    ZircaloyPhases,
    Neutral1,
    Neutral2,
    Neutral3,
    ConcreteDrying,
    TotalFluidPressure,
    VolumetricStrain,
    SteelTemperPhases,
    NumberOfExternVarTypes
};

class ExternalVariableTraits {
  public:
    static std::string getExternVarTypeStr( const externVarEnumInt e ) {
        if ( e == externVarEnumInt::Unknown ) {
            return std::string( "XXXX" );
        } else if ( e == externVarEnumInt::Temperature ) {
            return std::string( "TEMP" );
        } else if ( e == externVarEnumInt::Geometry ) {
            return std::string( "GEOM" );
        } else if ( e == externVarEnumInt::Corrosion ) {
            return std::string( "CORR" );
        } else if ( e == externVarEnumInt::IrreversibleStrain ) {
            return std::string( "EPSA" );
        } else if ( e == externVarEnumInt::ConcreteHydration ) {
            return std::string( "HYDR" );
        } else if ( e == externVarEnumInt::Irradiation ) {
            return std::string( "IRRA" );
        } else if ( e == externVarEnumInt::SteelPhases ) {
            return std::string( "M_ACIER" );
        } else if ( e == externVarEnumInt::ZircaloyPhases ) {
            return std::string( "M_ZIRC" );
        } else if ( e == externVarEnumInt::Neutral1 ) {
            return std::string( "NEUT1" );
        } else if ( e == externVarEnumInt::Neutral2 ) {
            return std::string( "NEUT2" );
        } else if ( e == externVarEnumInt::Neutral3 ) {
            return std::string( "NEUT3" );
        } else if ( e == externVarEnumInt::ConcreteDrying ) {
            return std::string( "SECH" );
        } else if ( e == externVarEnumInt::TotalFluidPressure ) {
            return std::string( "PTOT" );
        } else if ( e == externVarEnumInt::VolumetricStrain ) {
            return std::string( "DIVU" );
        } else if ( e == externVarEnumInt::SteelTemperPhases ) {
            return std::string( "M_REVENU" );
        } else {
            AS_ABORT( "Unknown external state variables" );
        }
    }

    static bool externVarHasRefeValue( const externVarEnumInt e ) {
        if ( e == externVarEnumInt::Temperature || e == externVarEnumInt::ConcreteDrying ) {
            return true;
        } else {
            return false;
        }
    }

    static externVarEnumInt getExternVarTypeInt( const std::string e ) {
        if ( e == "TEMP" ) {
            return externVarEnumInt::Temperature;
        } else if ( e == "GEOM" ) {
            return externVarEnumInt::Geometry;
        } else if ( e == "CORR" ) {
            return externVarEnumInt::Corrosion;
        } else if ( e == "EPSA" ) {
            return externVarEnumInt::IrreversibleStrain;
        } else if ( e == "HYDR" ) {
            return externVarEnumInt::ConcreteHydration;
        } else if ( e == "IRRA" ) {
            return externVarEnumInt::Irradiation;
        } else if ( e == "M_ACIER" ) {
            return externVarEnumInt::SteelPhases;
        } else if ( e == "M_ZIRC" ) {
            return externVarEnumInt::ZircaloyPhases;
        } else if ( e == "NEUT1" ) {
            return externVarEnumInt::Neutral1;
        } else if ( e == "NEUT2" ) {
            return externVarEnumInt::Neutral2;
        } else if ( e == "NEUT3" ) {
            return externVarEnumInt::Neutral3;
        } else if ( e == "SECH" ) {
            return externVarEnumInt::ConcreteDrying;
        } else if ( e == "PTOT" ) {
            return externVarEnumInt::TotalFluidPressure;
        } else if ( e == "DIVU" ) {
            return externVarEnumInt::VolumetricStrain;
        } else if ( e == "M_REVENU" ) {
            return externVarEnumInt::SteelTemperPhases;
        } else {
            return externVarEnumInt::Unknown;
        }
    }
    static bool externVarHasStrain( const externVarEnumInt e ) {
        if ( e == externVarEnumInt::Temperature ) {
            return true;
        } else if ( e == externVarEnumInt::Geometry ) {
            return false;
        } else if ( e == externVarEnumInt::Corrosion ) {
            return false;
        } else if ( e == externVarEnumInt::IrreversibleStrain ) {
            return true;
        } else if ( e == externVarEnumInt::ConcreteHydration ) {
            return true;
        } else if ( e == externVarEnumInt::Irradiation ) {
            return false;
        } else if ( e == externVarEnumInt::SteelPhases ) {
            return true;
        } else if ( e == externVarEnumInt::ZircaloyPhases ) {
            return true;
        } else if ( e == externVarEnumInt::Neutral1 ) {
            return false;
        } else if ( e == externVarEnumInt::Neutral2 ) {
            return false;
        } else if ( e == externVarEnumInt::Neutral3 ) {
            return false;
        } else if ( e == externVarEnumInt::ConcreteDrying ) {
            return true;
        } else if ( e == externVarEnumInt::TotalFluidPressure ) {
            return true;
        } else if ( e == externVarEnumInt::VolumetricStrain ) {
            return false;
        } else if ( e == externVarEnumInt::SteelTemperPhases ) {
            return true;
        } else {
            return false;
        }
    }
    static std::string getExternVarOption( const externVarEnumInt e ) {
        if ( e == externVarEnumInt::Temperature ) {
            return std::string( "CHAR_MECA_TEMP_R" );
        } else if ( e == externVarEnumInt::IrreversibleStrain ) {
            return std::string( "CHAR_MECA_EPSA_R" );
        } else if ( e == externVarEnumInt::ConcreteHydration ) {
            return std::string( "CHAR_MECA_HYDR_R" );
        } else if ( e == externVarEnumInt::SteelPhases ) {
            return std::string( "CHAR_MECA_META_Z" );
        } else if ( e == externVarEnumInt::ZircaloyPhases ) {
            return std::string( "CHAR_MECA_META_Z" );
        } else if ( e == externVarEnumInt::ConcreteDrying ) {
            return std::string( "CHAR_MECA_SECH_R" );
        } else if ( e == externVarEnumInt::TotalFluidPressure ) {
            return std::string( "CHAR_MECA_PTOT_R" );
        } else if ( e == externVarEnumInt::SteelTemperPhases ) {
            return std::string( "CHAR_MECA_META_Z" );
        } else {
            AS_ABORT( "No strain option for this external state variable." );
        }
    }
};

/**
 * @class ExternalStateVariable
 * @brief External state variable
 */
class ExternalStateVariable {
  private:
    BaseMeshPtr _mesh;
    MeshEntityPtr _localization;
    ASTERDOUBLE _refValue;
    DataFieldPtr _field;
    EvolutionParameterPtr _evolParameter;
    externVarEnumInt _type;

  public:
    /** @brief Constructor on all mesh */
    ExternalStateVariable( const externVarEnumInt _currType, const BaseMeshPtr mesh )
        : _type( _currType ),
          _mesh( mesh ),
          _localization( new AllMeshEntities() ),
          _refValue( 0. ),
          _evolParameter( nullptr ) {};

    /** @brief Constructor on all mesh */
    ExternalStateVariable( const std::string _currType, const BaseMeshPtr mesh )
        : _type( ExternalVariableTraits::getExternVarTypeInt( _currType ) ),
          _mesh( mesh ),
          _localization( new AllMeshEntities() ),
          _refValue( 0. ),
          _evolParameter( nullptr ) {};

    /** @brief Constructor on group of cells */
    ExternalStateVariable( const externVarEnumInt _currType, const BaseMeshPtr mesh,
                           const std::string nameOfGroup )
        : _type( _currType ),
          _mesh( mesh ),
          _localization( new GroupOfCells( nameOfGroup ) ),
          _refValue( 0. ),
          _evolParameter( nullptr ) {
        if ( !_mesh->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );
    }

    /** @brief Constructor on group of cells */
    ExternalStateVariable( const std::string _currType, const BaseMeshPtr mesh,
                           const std::string nameOfGroup )
        : _type( ExternalVariableTraits::getExternVarTypeInt( _currType ) ),
          _mesh( mesh ),
          _localization( new GroupOfCells( nameOfGroup ) ),
          _refValue( 0. ),
          _evolParameter( nullptr ) {
        if ( !_mesh->hasGroupOfCells( nameOfGroup ) )
            throw std::runtime_error( nameOfGroup + " not in mesh" );
    }

    /** @brief Destructor */
    ~ExternalStateVariable() {};

    /** @brief restricted constructor (Set) and method (Get) to support pickling */
    ExternalStateVariable( const py::tuple &tup );
    py::tuple _getState() const;

    /** @brief Function to know if a reference value exists */
    bool isSetRefe() const { return ExternalVariableTraits::externVarHasRefeValue( _type ); };

    /** @brief Get evolution parameter */
    EvolutionParameterPtr getEvolutionParameter() const { return _evolParameter; };

    /** @brief Get transient result */
    TransientResultPtr getTransientResult() const;

    /** @brief Get the field of values */
    DataFieldPtr getField() const { return _field; };

    /** @brief Get the reference value */
    ASTERDOUBLE getReferenceValue() const {
        if ( !isSetRefe() )
            throw std::runtime_error( "Reference value not set" );
        return _refValue;
    };

    /** @brief Function to set the evolution of external state variable */
    void setEvolutionParameter( const EvolutionParameterPtr &evolParameter ) {
        _field = nullptr;
        _evolParameter = evolParameter;
    };

    /**  @brief Function to set the field of external state variable */
    void setField( const DataFieldPtr &field ) {
        _field = field;
        _evolParameter = nullptr;
    };

    /** @brief Function to set the reference value */
    void setReferenceValue( const ASTERDOUBLE &value );

    /** @brief Get type of external state variable */
    externVarEnumInt getType() const { return _type; }

    /** @brief Get mesh */
    BaseMeshPtr getMesh() { return _mesh; };

    /** @brief Get localization */
    MeshEntityPtr getLocalization() { return _localization; };
};

using ExternalStateVariablePtr = std::shared_ptr< ExternalStateVariable >;

#endif /* EXTERNALSTATEVARIABLES_H_ */
