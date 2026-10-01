/**
 * @file DiscreteComputation.cxx
 * @brief Implementation of class DiscreteComputation
 * @section LICENCE
 *   Copyright (C) 1991 2026  EDF www.code-aster.org
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
#include "aster_fort_superv.h"
#include "aster_fort_utils.h"

#include "DataFields/FieldOnCellsBuilder.h"
#include "Discretization/Calcul.h"
#include "Discretization/DiscreteComputation.h"
#include "LinearAlgebra/ElementaryMatrixConverter.h"
#include "Loads/DirichletBC.h"
#include "Loads/MechanicalLoad.h"
#include "Materials/MaterialField.h"
#include "MemoryManager/JeveuxVector.h"
#include "Modeling/Model.h"
#include "Modeling/ParallelContactFEDescriptor.h"
#include "Modeling/XfemModel.h"
#include "Utilities/Tools.h"

#include <limits>

FieldOnNodesRealPtr DiscreteComputation::getDifferentialDualDisplacement() const {
    // reimplement nmdidi.F90 & vecdid.F90
    auto chalph = std::make_shared< ConstantFieldOnCellsReal >( _phys_problem->getMesh() );
    const std::string physicalName( "NEUT_R" );
    chalph->allocate( physicalName );
    ConstantFieldOnZone a( _phys_problem->getMesh() );
    ConstantFieldValues< ASTERDOUBLE > b( { "X1" }, { 1.0 } );
    chalph->setValueOnZone( a, b );

    auto elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_MECA" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();

    FieldOnNodesRealPtr disp_didi = listOfLoads->getDifferentialDisplacement();

    auto calcul = std::make_unique< Calcul >( calcul_option );
    calcul->setOption( "MECA_BU_R" );

    auto impl = [&]( auto &load, const ASTERINTEGER &iload ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();
        auto mult_field = load->getMultiplicativeField();
        if ( mult_field && mult_field->exists() && load_FEDesc ) {
            calcul->clearInputs();
            calcul->clearOutputs();
            calcul->setFiniteElementDescriptor( load_FEDesc );
            calcul->addInputField( "PDDLMUR", mult_field );
            calcul->addInputField( "PDDLIMR", disp_didi );
            calcul->addInputField( "PALPHAR", chalph );

            calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ),
                                             iload );
            }
        }
    };

    auto types = listOfLoads->getListOfMechaTyp();
    for ( const auto &load : listOfLoads->getMechanicalLoadsReal() ) {
        if ( types[iload - 1] == "DIDI" )
            impl( load, iload );
        iload++;
    }

    auto nloads = types.size();
    types = listOfLoads->getListOfMechaFuncTyp();
    for ( const auto &load : listOfLoads->getMechanicalLoadsFunction() ) {
        if ( types[iload - 1 - nloads] == "DIDI" )
            impl( load, iload );
        iload++;
    }

    elemVect->build();

    if ( elemVect->hasElementaryTerm() ) {
        return elemVect->assemble( _phys_problem->getDOFNumbering() );
    } else {
        FieldOnNodesRealPtr vectAsse =
            std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
        vectAsse->setValues( 0.0 );
        vectAsse->build();
        return vectAsse;
    }
}

FieldOnNodesRealPtr DiscreteComputation::getDualDisplacement( FieldOnNodesRealPtr disp_curr,
                                                              ASTERDOUBLE scaling ) const {

    auto elemVect = std::make_shared< ElementaryVectorReal >( _phys_problem->getModel() );

    if ( !_phys_problem->isMechanical() ) {
        AS_ABORT( "Not implemented" );
    };

    // Prepare loads
    auto listOfLoads = _phys_problem->getListOfLoads();
    std::string listLoadsName = ljust( listOfLoads->getName(), 19 );

    // Get JEVEUX names of objects to call Fortran
    std::string modelName = ljust( _phys_problem->getModel()->getName(), 24 );
    std::string dispName = ljust( disp_curr->getName(), 24 );
    std::string vectElemName = ljust( elemVect->getName(), 24 );
    const std::string base( "G" );
    const ASTERDOUBLE const_scaling = scaling;

    // Wrapper FORTRAN
    CALLO_VEBUME( modelName, dispName, listLoadsName, vectElemName, &const_scaling, base );

    // Construct vect_elem object
    auto FEDs = _phys_problem->getListOfLoads()->getFiniteElementDescriptors();
    FEDs.push_back( _phys_problem->getModel()->getFiniteElementDescriptor() );
    elemVect->build( FEDs );

    // Assemble
    FieldOnNodesRealPtr bume = elemVect->assemble( _phys_problem->getDOFNumbering() );

    if ( _phys_problem->getMesh()->isParallel() )
        CALLO_AP_ASSEMBLY_VECTOR( bume->getName() );

    if ( listOfLoads->hasDifferentialDisplacement() )
        bume->operator-=( *getDifferentialDualDisplacement() );

    return bume;
};

std::variant< ElementaryVectorDisplacementRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getMechanicalNeumannForces( const ASTERDOUBLE time_curr,
                                                 const ASTERDOUBLE time_step,
                                                 const ASTERDOUBLE theta,
                                                 const ASTERINTEGER modeFourier,
                                                 const FieldOnCellsRealPtr varc_curr,
                                                 const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    auto elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_MECA" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );

    if ( currMater && currMater->hasExternalStateVariable() ) {
        if ( !varc_curr || !varc_curr->exists() ) {
            raiseAsterError( "External state variables are needed but not given" );
        }
    }

    auto calcul = std::make_unique< Calcul >( calcul_option );

    auto impl = [&]( auto load, const ASTERINTEGER &load_i, const std::string &option,
                     const std::string &name, const std::string &param,
                     const FiniteElementDescriptorPtr FED,
                     std::vector< std::pair< std::string, DataFieldPtr > > field_in =
                         std::vector< std::pair< std::string, DataFieldPtr > >() ) {
        if ( load->hasLoadField( name ) ) {
            calcul->setOption( option );
            calcul->setFiniteElementDescriptor( FED );

            calcul->clearInputs();
            calcul->addTimeField( "PINSTR", time_curr, time_step, theta );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            if ( currMater ) {
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
                calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
            }
            if ( varc_curr ) {
                calcul->addInputField( "PVARCPR", varc_curr );
            }
            if ( currElemChara ) {
                calcul->addElementaryCharacteristicsField( currElemChara );
                calcul->addOrientationField( currElemChara );
            }

            calcul->addXFEMField( currModel );
            calcul->addHHOField( currModel );

            if ( !param.empty() ) {
                calcul->addInputField( param, load->getConstantLoadField( name ) );
            }
            calcul->addFourierModeField( modeFourier );

            for ( auto &[param_in, field] : field_in ) {
                if ( field && field->exists() ) {
                    calcul->addInputField( param_in, field );
                }
            }

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ),
                                             load_i );
            }
        }
    };

    auto mecaLoadReal = listOfLoads->getMechanicalLoadsReal();
    for ( const auto &load : mecaLoadReal ) {

        impl( load, iload, "CHAR_MECA_FRCO2D", "FCO2D", "PFRCO2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FRCO3D", "FCO3D", "PFRCO3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FR2D3D", "F2D3D", "PFR2D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FR1D3D", "F1D3D", "PFR1D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FR1D2D", "F1D2D", "PFR1D2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_PRES_R", "PRESS", "PPRESSR", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FLUX_R", "FLUX", "PFLUXR", model_FEDesc );
        impl( load, iload, "CHAR_MECA_EFON_R", "EFOND", "PEFOND", model_FEDesc,
              { { "PPREFFR", load->getConstantLoadField( "PREFF" ) } } );
        iload++;
    }

    auto mecaLoadFunc = listOfLoads->getMechanicalLoadsFunction();
    for ( const auto &load : mecaLoadFunc ) {

        impl( load, iload, "CHAR_MECA_FFCO2D", "FCO2D", "PFFCO2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FFCO3D", "FCO3D", "PFFCO3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FF2D3D", "F2D3D", "PFF2D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FF1D3D", "F1D3D", "PFF1D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FF1D2D", "F1D2D", "PFF1D2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_PRES_F", "PRESS", "PPRESSF", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FLUX_F", "FLUX", "PFLUXF", model_FEDesc );
        impl( load, iload, "CHAR_MECA_EFON_F", "EFOND", "PEFOND", model_FEDesc,
              { { "PPREFFF", load->getConstantLoadField( "PREFF" ) } } );

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions(
                _phys_problem->getDOFNumbering(), _phys_problem->getListOfLoads(), time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

std::variant< ElementaryVectorDisplacementRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getMechanicalVolumetricForces( const ASTERDOUBLE time_curr,
                                                    const ASTERDOUBLE time_step,
                                                    const ASTERDOUBLE theta,
                                                    const ASTERINTEGER modeFourier,
                                                    const FieldOnCellsRealPtr varc_curr,
                                                    const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    auto elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_MECA" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto listOfLoads = _phys_problem->getListOfLoads();
    auto model_FEDesc = currModel->getFiniteElementDescriptor();
    AS_ASSERT( model_FEDesc );

    if ( currMater && currMater->hasExternalStateVariable() ) {
        if ( !varc_curr || !varc_curr->exists() ) {
            raiseAsterError( "External state variables are needed but not given" );
        }
    }

    auto calcul = std::make_unique< Calcul >( calcul_option );

    auto impl = [&]( auto load, const ASTERINTEGER &load_i, const std::string &option,
                     const std::string &name, const std::string &param,
                     const FiniteElementDescriptorPtr FED,
                     std::vector< std::pair< std::string, DataFieldPtr > > field_in =
                         std::vector< std::pair< std::string, DataFieldPtr > >() ) {
        if ( load->hasLoadField( name ) ) {
            calcul->setOption( option );
            calcul->setFiniteElementDescriptor( FED );

            calcul->clearInputs();
            calcul->addTimeField( "PINSTR", time_curr, time_step, theta );
            calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
            if ( currMater ) {
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
                calcul->addInputField( "PCOMPOR", currMater->getBehaviourField() );
            }
            if ( varc_curr ) {
                calcul->addInputField( "PVARCPR", varc_curr );
            }
            if ( currElemChara ) {
                calcul->addElementaryCharacteristicsField( currElemChara );
                calcul->addOrientationField( currElemChara );
            }

            calcul->addXFEMField( currModel );
            calcul->addHHOField( currModel );

            if ( !param.empty() ) {
                calcul->addInputField( param, load->getConstantLoadField( name ) );
            }
            calcul->addFourierModeField( modeFourier );

            for ( auto &[param_in, field] : field_in ) {
                if ( field && field->exists() ) {
                    calcul->addInputField( param_in, field );
                }
            }

            calcul->clearOutputs();
            calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );
            calcul->compute();
            if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ),
                                             load_i );
            }
        }
    };

    auto mecaLoadReal = listOfLoads->getMechanicalLoadsReal();
    for ( const auto &load : mecaLoadReal ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();

        impl( load, iload, "CHAR_MECA_FORC_R", "FORNO", "PFORNOR", load_FEDesc );
        impl( load, iload, "CHAR_MECA_FR3D3D", "F3D3D", "PFR3D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FR2D2D", "F2D2D", "PFR2D2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FR1D1D", "F1D1D", "PFR1D1D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_PESA_R", "PESAN", "PPESANR", model_FEDesc );
        impl( load, iload, "CHAR_MECA_ROTA_R", "ROTAT", "PROTATR", model_FEDesc );
        impl( load, iload, "CHAR_MECA_EPSI_R", "EPSIN", "PEPSINR", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FRELEC", "FELEC", "PFRELEC", model_FEDesc );
        impl( load, iload, "CHAR_MECA_ONDE", "ONDE", "PONDECR", model_FEDesc );

        // PRE_SIGM
        if ( load->hasLoadField( "SIINT" ) ) {
            // TODO: Do not create a copy - circular inclusion
            auto pre_sgm_name =
                ( load->getConstantLoadFieldChar8( "SIINT" )->getValues( 0 ).getValues()[0] )
                    .toString();

            const std::string typeco( "CHAMP" );
            auto [ok, repi, retour] = dismoi( "DOCU", pre_sgm_name, typeco, true );

            if ( retour == "CHML" ) {
                auto pre_sigm = std::make_shared< FieldOnCellsReal >();
                std::string base = "G";
                CALLO_COPISD( typeco, base, pre_sgm_name, pre_sigm->getName() );
                impl( load, iload, "FORC_NODA", "SIINT", "", model_FEDesc,
                      { { "PSIEFR", pre_sigm } } );
            } else if ( retour == "CART" ) {
                auto pre_sigm =
                    std::make_shared< ConstantFieldOnCellsReal >( currModel->getMesh() );
                std::string base = "G";
                CALLO_COPISD( typeco, base, pre_sgm_name, pre_sigm->getName() );
                impl( load, iload, "FORC_NODA", "SIINT", "", model_FEDesc,
                      { { "PSIEFR", pre_sigm } } );
            } else {
                AS_ABORT( "Error: " + retour );
            }
        }

        // CHAR_MECA_EVOL
        if ( load->hasLoadResult() ) {
            std::string cara_elem = " ", mater = " ", mateco = " ";
            if ( currElemChara ) {
                cara_elem = currElemChara->getName();
            }
            if ( currMater ) {
                mater = currMater->getName();
                mateco = currCodedMater->getCodedMaterialField()->getName();
            }
            ASTERDOUBLE inst_prev = time_curr - time_step, inst_curr = time_curr,
                        inst_theta = theta;
            auto resu_elem = std::make_shared< ElementaryTermReal >();
            ASTERINTEGER nharm = modeFourier;
            std::string base = "G";

            CALLO_ME2MME_EVOL( currModel->getName(), cara_elem, mater, mateco, &nharm, base, &iload,
                               load->getName(), model_FEDesc->getName(), &inst_prev, &inst_curr,
                               &inst_theta, resu_elem->getName(), elemVect->getName() );
        }

        // VECT_ASSE
        if ( load->hasLoadVectAsse() ) {
            // TODO: Do not create a copy - circular inclusion
            auto veass = std::make_shared< FieldOnNodesReal >();
            std::string type = "CHAMP_GD", base = "G";
            CALLO_COPISD( type, base, load->getLoadVectAsseName(), veass->getName() );
            veass->build( currModel->getMesh() );
            elemVect->setVeass( veass, iload );
        }
        iload++;
    }

    auto mecaLoadFunc = listOfLoads->getMechanicalLoadsFunction();
    for ( const auto &load : mecaLoadFunc ) {
        auto load_FEDesc = load->getFiniteElementDescriptor();

        impl( load, iload, "CHAR_MECA_FORC_F", "FORNO", "PFORNOF", load_FEDesc );
        impl( load, iload, "CHAR_MECA_FF3D3D", "F3D3D", "PFF3D3D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FF2D2D", "F2D2D", "PFF2D2D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_FF1D1D", "F1D1D", "PFF1D1D", model_FEDesc );
        impl( load, iload, "CHAR_MECA_EPSI_F", "EPSIN", "PEPSINF", model_FEDesc );
        impl( load, iload, "ONDE_PLAN", "ONDPL", "PONDPLA", model_FEDesc,
              { { "PONDPLR", load->getConstantLoadField( "ONDPR" ) } } );

        iload++;
    }

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions(
                _phys_problem->getDOFNumbering(), _phys_problem->getListOfLoads(), time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

/** @brief Compute internal forces, stress and internal state variables */
std::tuple< FieldOnCellsLongPtr, ASTERINTEGER, FieldOnCellsRealPtr, FieldOnCellsRealPtr,
            FieldOnNodesRealPtr >
DiscreteComputation::getInternalMechanicalForces(
    const FieldOnNodesRealPtr displ_prev, const FieldOnNodesRealPtr displ_step,
    const FieldOnCellsRealPtr stress, const FieldOnCellsRealPtr internVar,
    const FieldOnCellsRealPtr internVarIter, const ASTERDOUBLE &time_prev,
    const ASTERDOUBLE &time_step, const FieldOnCellsRealPtr &varc_prev,
    const FieldOnCellsRealPtr &varc_curr, const VectorString &groupOfCells ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehaviour = _phys_problem->getBehaviourProperty();

    // Select option to compute
    std::string option = "RAPH_MECA";

    // Prepare computing:
    CalculPtr calcul = createCalculForNonLinear( option, time_prev, time_prev + time_step,
                                                 varc_prev, varc_curr, groupOfCells );
    FiniteElementDescriptorPtr FEDesc = calcul->getFiniteElementDescriptor();

    // Set current physical state
    calcul->addInputField( "PDEPLMR", displ_prev );
    calcul->addInputField( "PDEPLPR", displ_step );
    calcul->addInputField( "PCONTMR", stress );
    calcul->addInputField( "PVARIMR", internVar );

    // Coded Material
    auto currCodedMater = _phys_problem->getCodedMaterial();
    calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

    // Provisoire: pour TANGENTE=VERIFICATION, nécessité de variables internes à chaque itération
    // Nécessaire également pour Deborst
    calcul->addInputField( "PVARIMP", internVarIter );

    calcul->addHHOField( currModel );

    // Create output vector
    auto elemVect = std::make_shared< ElementaryVectorReal >( _phys_problem->getModel() );

    // Create output fields
    auto stress_curr =
        FieldOnCellsPtrBuilder< ASTERDOUBLE >( FEDesc, "ELGA", "SIEF_R", currElemChara );
    FieldOnCellsLongPtr exitField = std::make_shared< FieldOnCellsLong >( FEDesc );
    auto vari_curr = FieldOnCellsPtrBuilder< ASTERDOUBLE >( FEDesc, "ELGA", "VARI_R", currBehaviour,
                                                            currElemChara );

    // Add output fields
    calcul->addOutputField( "PVARIPR", vari_curr );
    calcul->addOutputField( "PCONTPR", stress_curr );
    calcul->addOutputField( "PCODRET", exitField );

    // Add output elementary
    calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );

    // Compute and assemble vector
    FieldOnNodesRealPtr internalForces;
    if ( currModel->existsFiniteElement() ) {
        calcul->compute();
        if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) )
            elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ) );
        elemVect->build();
        internalForces = elemVect->assemble( _phys_problem->getDOFNumbering() );
    };

    std::string exitFieldName = ljust( exitField->getName(), 19 );
    ASTERINTEGER exitCode = 0;
    CALLO_GETERRORCODE( exitFieldName, &exitCode );
#ifdef ASTER_HAVE_MPI
    exitCode = AsterMPI::max( exitCode );
#endif

    return std::make_tuple( exitField, exitCode, vari_curr, stress_curr, internalForces );
}

/** @brief Compute AFFE_CHAR_MECA DDL_IMPO */
std::variant< ElementaryVectorDisplacementRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getMechanicalImposedDualBC( const ASTERDOUBLE time_curr,
                                                 const bool assembly ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    auto elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >( _phys_problem->getModel() );

    // Init
    ASTERINTEGER iload = 1;

    // Setup
    const std::string calcul_option( "CHAR_MECA" );

    // Main parameters
    auto currModel = _phys_problem->getModel();
    auto listOfLoads = _phys_problem->getListOfLoads();

    auto calcul = std::make_unique< Calcul >( calcul_option );

    auto impl_disp = [&]( auto loads, bool real ) {
        std::string name;
        if ( real ) {
            calcul->setOption( "MECA_DDLI_R" );
            name = "PDDLIMR";
        } else {
            calcul->setOption( "MECA_DDLI_F" );
            name = "PDDLIMF";
        }
        for ( const auto &load : loads ) {
            auto load_FEDesc = load->getFiniteElementDescriptor();
            auto impo_field = load->getImposedField();
            if ( impo_field && impo_field->exists() && load_FEDesc ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( load_FEDesc );
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                calcul->addInputField( name, impo_field );
                if ( !real ) {
                    calcul->addTimeField( "PINSTR", time_curr );
                }
                calcul->addOutputElementaryTerm( "PVECTUR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ),
                                                 iload );
                }
            }
            iload++;
        }
    };

    impl_disp( listOfLoads->getMechanicalLoadsReal(), true );
    impl_disp( listOfLoads->getMechanicalLoadsFunction(), false );

#ifdef ASTER_HAVE_MPI
    impl_disp( listOfLoads->getParallelMechanicalLoadsReal(), true );
    impl_disp( listOfLoads->getParallelMechanicalLoadsFunction(), false );
#endif

    auto impl_vite = [&]( auto loads, bool real ) {
        std::string name;
        if ( real ) {
            calcul->setOption( "CHAR_MECA_VFAC" );
            name = "PVITEFR";
        } else {
            calcul->setOption( "CHAR_MECA_VFAC_F" );
            name = "PVITEFF";
        }

        for ( const auto &load : loads ) {
            if ( load->hasLoadField( "VFACE" ) ) {
                calcul->clearInputs();
                calcul->clearOutputs();
                calcul->setFiniteElementDescriptor( currModel->getFiniteElementDescriptor() );
                calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
                auto currCodedMater = _phys_problem->getCodedMaterial();
                calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
                calcul->addInputField( name, load->getConstantLoadField( "VFACE" ) );
                if ( !real ) {
                    calcul->addTimeField( "PINSTR", time_curr );
                }
                calcul->addOutputElementaryTerm( "PVECTUR",
                                                 std::make_shared< ElementaryTermReal >() );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ),
                                                 iload );
                }
            }
            iload++;
        }
    };

    impl_vite( listOfLoads->getMechanicalLoadsReal(), true );
    impl_vite( listOfLoads->getMechanicalLoadsFunction(), false );

    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assembleWithLoadFunctions(
                _phys_problem->getDOFNumbering(), _phys_problem->getListOfLoads(), time_curr );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
}

FieldOnNodesRealPtr DiscreteComputation::getContactForces(
    const MeshCoordinatesFieldPtr geom, const FieldOnNodesRealPtr displ_prev,
    const FieldOnNodesRealPtr displ_step, const ASTERDOUBLE &time_prev,
    const ASTERDOUBLE &time_step, const FieldOnCellsRealPtr data,
    const FieldOnNodesRealPtr coef_cont, const FieldOnNodesRealPtr coef_frot ) const {

    // Select option for matrix
    std::string option = "CHAR_MECA_CONT";

    auto Fed_pair = _phys_problem->getContactFED();
#ifdef ASTER_HAVE_MPI
    const auto pCFED = std::dynamic_pointer_cast< ParallelContactFEDescriptor >( Fed_pair );
    if ( pCFED ) {
        Fed_pair = pCFED->getSupportFiniteElementDescriptor();
    }
#endif /* ASTER_HAVE_MPI */

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );
    calcul->setFiniteElementDescriptor( Fed_pair );

    // Set input field
#ifdef ASTER_HAVE_MPI
    if ( pCFED ) {
        calcul->addInputField( "PGEOMER", Fed_pair->getMesh()->getCoordinates() );
    } else {
#endif /* ASTER_HAVE_MPI */
        calcul->addInputField( "PGEOMER", _phys_problem->getMesh()->getCoordinates() );
#ifdef ASTER_HAVE_MPI
    }
#endif /* ASTER_HAVE_MPI */
    calcul->addInputField( "PGEOMCR", geom );
    calcul->addInputField( "PDEPL_M", displ_prev );
    calcul->addInputField( "PDEPL_P", displ_step );
    calcul->addInputField( "PCONFR", data );
    calcul->addInputField( "PCCONTR", coef_cont );
    calcul->addInputField( "PCFROTR", coef_frot );

    // Add time fields
    calcul->addTimeField( "PINSTMR", time_prev );
    calcul->addTimeField( "PINSTPR", time_prev + time_step );

    // Create output vector
    auto elemVect =
        std::make_shared< ElementaryVectorDisplacementReal >( _phys_problem->getModel() );

    // Add output elementary
    calcul->addOutputElementaryTerm( "PVECTCR", std::make_shared< ElementaryTermReal >() );
    calcul->addOutputElementaryTerm( "PVECTFR", std::make_shared< ElementaryTermReal >() );

    // Computation
    calcul->compute();
    if ( calcul->hasOutputElementaryTerm( "PVECTCR" ) ) {
        elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTCR" ) );
    }
    if ( calcul->hasOutputElementaryTerm( "PVECTFR" ) ) {
        elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTFR" ) );
    }
    elemVect->build();
#ifdef ASTER_HAVE_MPI
    if ( pCFED ) {
        auto resu = transfertToParallelFEDesc( elemVect, pCFED );
        resu->build();
        return resu->assemble( _phys_problem->getDOFNumbering() );
    } else {
#endif /* ASTER_HAVE_MPI */
        return elemVect->assemble( _phys_problem->getDOFNumbering() );
#ifdef ASTER_HAVE_MPI
    }
#endif /* ASTER_HAVE_MPI */
}

/**  @brief Compute nodal forces  */
std::variant< ElementaryVectorDisplacementRealPtr, FieldOnNodesRealPtr >
DiscreteComputation::getMechanicalNodalForces( const FieldOnCellsRealPtr stress,
                                               const FieldOnNodesRealPtr disp,
                                               const ASTERINTEGER modeFourier,
                                               const FieldOnCellsRealPtr varc_curr,
                                               const ConstantFieldOnCellsChar16Ptr behaviourMap,
                                               const VectorString &groupOfCells,
                                               const bool assembly ) const {

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currBehav = _phys_problem->getBehaviourProperty();

    // Some checks
    AS_ASSERT( currModel->isMechanical() );
    AS_ASSERT( currMater );
    if ( currMater && currMater->hasExternalStateVariable() ) {
        if ( !varc_curr || !varc_curr->exists() ) {
            raiseAsterError( "External state variables are needed but not given" );
        }
    }
    if ( currModel->existsSTRX() ) {
        throw std::runtime_error( "Beams not yet implemented" );
    }

    // Preparation of elementary vector
    const std::string vect_option( "CHAR_MECA" );
    auto elemVect = std::make_shared< ElementaryVectorDisplacementReal >( currModel );

    // Create object for calcul
    const std::string calcul_option( "FORC_NODA" );
    auto calcul = std::make_unique< Calcul >( calcul_option );

    // Prepare calcul
    if ( groupOfCells.empty() ) {
        calcul->setModel( currModel );
    } else {
        calcul->setGroupsOfCells( currModel, groupOfCells );
    }

    // Add input fields: geometry
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );

    // Add input fields: material parameters
    calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );

    // Add input fields: elementary characteristics
    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
        calcul->addOrientationField( currElemChara );
    }

    // Add input fields: for XFEM
    calcul->addXFEMField( currModel );

    calcul->addHHOField( currModel );

    // Add Fourier field
    calcul->addFourierModeField( modeFourier );

    // Add behaviour field (non-linear)
    if ( currBehav && !behaviourMap ) {
        calcul->addInputField( "PCOMPOR", currBehav->getBehaviourField() );
    };
    if ( behaviourMap ) {
        calcul->addInputField( "PCOMPOR", behaviourMap );
    };

    // Add external state variables
    if ( varc_curr ) {
        calcul->addInputField( "PVARCPR", varc_curr );
    }

    // Set current physical state
    if ( disp != nullptr ) {
        calcul->addInputField( "PDEPLAR", disp );
    }
    calcul->addInputField( "PSIEFR", stress );
    // calcul->addInputField( "PSTRXMR", strx );

    // Set output fields
    calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );

    // // Compute
    if ( calcul->hasComplexInputFields() ) {
        AS_ABORT( "Complex fields no supported" );
    } else {
        calcul->compute();
    }

    // Add RESU_ELEM to VECT_ELEM
    if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
        elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ) );
    }

    // Prepare output
    elemVect->build();

    if ( assembly ) {
        if ( elemVect->hasElementaryTerm() ) {
            return elemVect->assemble( _phys_problem->getDOFNumbering() );
        } else {
            FieldOnNodesRealPtr vectAsse =
                std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
            vectAsse->setValues( 0.0 );
            vectAsse->build();
            return vectAsse;
        }
    }

    return elemVect;
};

FieldOnNodesRealPtr
DiscreteComputation::getMechanicalForces( const ASTERDOUBLE time_curr, const ASTERDOUBLE time_step,
                                          const ASTERDOUBLE theta, const ASTERINTEGER modeFourier,
                                          const FieldOnCellsRealPtr varc_curr ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    const FieldOnNodesRealPtr neumannForces =
        std::get< FieldOnNodesRealPtr >( DiscreteComputation::getMechanicalNeumannForces(
            time_curr, time_step, theta, modeFourier, varc_curr ) );

    const FieldOnNodesRealPtr volumetricForces =
        std::get< FieldOnNodesRealPtr >( DiscreteComputation::getMechanicalVolumetricForces(
            time_curr, time_step, theta, modeFourier, varc_curr ) );

    FieldOnNodesRealPtr vectAsse =
        std::make_shared< FieldOnNodesReal >( *neumannForces + *volumetricForces );

    vectAsse->build();
    return vectAsse;
};

/** @brief Compute reaction forces */
FieldOnNodesRealPtr DiscreteComputation::getMechanicalReactionForces(
    const FieldOnCellsRealPtr stress, const FieldOnNodesRealPtr disp, const ASTERDOUBLE time_prev,
    const ASTERDOUBLE time_curr, const ASTERDOUBLE theta, const ASTERINTEGER modeFourier,
    const FieldOnCellsRealPtr varc_curr, const ConstantFieldOnCellsChar16Ptr behaviourMap ) const {

    const FieldOnNodesRealPtr nodalForces =
        std::get< FieldOnNodesRealPtr >( DiscreteComputation::getMechanicalNodalForces(
            stress, disp, modeFourier, varc_curr, behaviourMap ) );

    const ASTERDOUBLE time_step = time_curr - time_prev;

    const FieldOnNodesRealPtr mechanicalForces = DiscreteComputation::getMechanicalForces(
        time_curr, time_step, theta, modeFourier, varc_curr );

    FieldOnNodesRealPtr vectAsse =
        std::make_shared< FieldOnNodesReal >( *nodalForces - *mechanicalForces );

    vectAsse->build();
    return vectAsse;
};

FieldOnNodesRealPtr
DiscreteComputation::dualMechanicalVector( FieldOnNodesRealPtr lagr_curr ) const {

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    // Init
    ASTERINTEGER iload = 1;

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currMater = _phys_problem->getMaterialField();
    auto currListOfLoads = _phys_problem->getListOfLoads();

    // Select option to compute
    std::string option = "MECA_BTLA_R";

    // Create elementary vector
    auto elemVect = std::make_shared< ElementaryVectorReal >( currModel );

    // Setup
    const std::string calcul_option( "CHAR_MECA" );

    // Prepare computing
    CalculPtr calcul = std::make_unique< Calcul >( option );

    auto impl_load = [&]( auto loads ) {
        for ( const auto &load : loads ) {
            auto load_FEDesc = load->getFiniteElementDescriptor();
            auto mult_field = load->getMultiplicativeField();
            if ( mult_field && mult_field->exists() && load_FEDesc ) {

                // Add input fields
                calcul->clearInputs();
                calcul->addInputField( "PLAGRAR", lagr_curr );
                calcul->addInputField( "PDDLMUR", mult_field );

                // Add output fields
                calcul->clearOutputs();
                calcul->addOutputElementaryTerm( "PVECTUR",
                                                 std::make_shared< ElementaryTermReal >() );

                // Compute
                calcul->setFiniteElementDescriptor( load_FEDesc );
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ),
                                                 iload );
                }
            }
            iload++;
        }
    };

    impl_load( currListOfLoads->getMechanicalLoadsReal() );
    impl_load( currListOfLoads->getMechanicalLoadsFunction() );

#ifdef ASTER_HAVE_MPI
    impl_load( currListOfLoads->getParallelMechanicalLoadsReal() );
    impl_load( currListOfLoads->getParallelMechanicalLoadsFunction() );
#endif

    // Construct elementary vector
    auto FEDs = _phys_problem->getListOfLoads()->getFiniteElementDescriptors();
    FEDs.push_back( _phys_problem->getModel()->getFiniteElementDescriptor() );
    elemVect->build( FEDs );

    // Assemble
    return elemVect->assemble( _phys_problem->getDOFNumbering() );
};

FieldOnNodesRealPtr DiscreteComputation::getResidualReference(
    const std::map< std::string, ASTERDOUBLE > &vale_by_name ) const {
    // reimplement nmrefe.F90

    AS_ASSERT( _phys_problem->getModel()->isMechanical() );

    VectorString list_cmp = { "SIGM",   "EPSI", "FLUXTHER", "FLUXHYD1", "FLUXHYD2", "EFFORT",
                              "MOMENT", "VARI", "DEPL",     "LAGR",     "PI" };
    VectorString list_name = { "SIGM_REFE",      "EPSI_REFE",   "FLUX_THER_REFE", "FLUX_HYD1_REFE",
                               "FLUX_HYD2_REFE", "EFFORT_REFE", "MOMENT_REFE",    "VARI_REFE",
                               "DEPL_REFE",      "LAGR_REFE",   "PI_REFE" };
    VectorString list_neut = { "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11" };

    const ASTERDOUBLE r8vide = CALL_R8VIDE();
    VectorReal list_vale;
    list_vale.reserve( list_name.size() );
    for ( int i = 0; i < list_name.size(); ++i ) {
        auto it = vale_by_name.find( list_name[i] );
        if ( it != vale_by_name.end() )
            list_vale.push_back( it->second );
        else
            list_vale.push_back( r8vide );
    }

    ConstantFieldOnZone a( _phys_problem->getMesh() );

    auto chrefe_val = std::make_shared< ConstantFieldOnCellsReal >( _phys_problem->getMesh() );
    const std::string physicalName_val( "RESIREF" );
    chrefe_val->allocate( physicalName_val );
    AS_ASSERT( list_cmp.size() == list_vale.size() );
    ConstantFieldValues< ASTERDOUBLE > b_val( list_cmp, list_vale );
    chrefe_val->setValueOnZone( a, b_val );

    std::vector< JeveuxChar8 > list_cmp_k8;
    for ( const auto &s : list_cmp ) {
        list_cmp_k8.emplace_back( s );
    }
    auto chrefe_cmp = std::make_shared< ConstantFieldOnCellsChar8 >( _phys_problem->getMesh() );
    const std::string physicalName_cmp( "RESICMP" );
    chrefe_cmp->allocate( physicalName_cmp );
    ConstantFieldValues< JeveuxChar8 > b_cmp( list_neut, list_cmp_k8 );
    chrefe_cmp->setValueOnZone( a, b_cmp );

    // Get main parameters
    auto currModel = _phys_problem->getModel();
    auto currMater = _phys_problem->getMaterialField();
    auto currElemChara = _phys_problem->getElementaryCharacteristics();
    auto currListOfLoads = _phys_problem->getListOfLoads();
    auto currCodedMater = _phys_problem->getCodedMaterial();
    auto currBehaviour = _phys_problem->getBehaviourProperty();

    // Setup
    auto elemVect = std::make_shared< ElementaryVectorReal >( currModel );

    // Prepare computing
    auto calcul = std::make_unique< Calcul >( "REFE_FORC_NODA" );
    calcul->setModel( currModel );
    calcul->addBehaviourField( currBehaviour );
    if ( currElemChara ) {
        calcul->addElementaryCharacteristicsField( currElemChara );
        calcul->addOrientationField( currElemChara );
    }
    calcul->addInputField( "PGEOMER", currModel->getMesh()->getCoordinates() );
    calcul->addInputField( "PMATERC", currCodedMater->getCodedMaterialField() );
    calcul->addInputField( "PRESICMP", chrefe_cmp );
    calcul->addInputField( "PRESIREF", chrefe_val );
    calcul->addXFEMField( currModel );
    calcul->addHHOField( currModel );
    calcul->addOutputElementaryTerm( "PVECTUR", std::make_shared< ElementaryTermReal >() );
    calcul->compute();

    elemVect->addElementaryTerm( calcul->getOutputElementaryTermReal( "PVECTUR" ) );
    elemVect->build();
    FieldOnNodesRealPtr vectAsse = elemVect->assemble( _phys_problem->getDOFNumbering(), true );

    // for (ASTERINTEGER i=0; i<values->size(); i++)
    //     std::cout << ( *values )[i] << std::endl;

    return vectAsse;
};

FieldOnNodesRealPtr DiscreteComputation::getMechanicalCouplingForces(
    const FieldOnNodesRealPtr displ_prev, const FieldOnNodesRealPtr displ_step,
    const ASTERDOUBLE &time_prev, const ASTERDOUBLE &time_step ) const {

    // Select option for vector
    std::string option = "CHAR_MECA_CPL";

    const auto model = _phys_problem->getModel();

    AS_ASSERT( model->isMechanical() );

    auto elemVect = std::make_shared< ElementaryVectorDisplacementReal >( model );

    const auto isHPC = model->getMesh()->isParallel();

    auto getFED = []( const auto &load ) {
#ifdef ASTER_HAVE_MPI
        if constexpr ( std::is_same_v< std::decay_t< decltype( load ) >,
                                       ParallelMechanicalLoadRealPtr > ) {
            return load->getParallelCouplingFEDescriptor()->getSupportFiniteElementDescriptor();
        } else
#endif
        {
            return load->getPairingField()->getDescription();
        }
    };

    auto addET = []( const auto &load, const ElementaryTermRealPtr &ET ) {
#ifdef ASTER_HAVE_MPI
        if constexpr ( std::is_same_v< std::decay_t< decltype( load ) >,
                                       ParallelMechanicalLoadRealPtr > ) {
            return transfertToParallelFEDesc( ET, load->getParallelCouplingFEDescriptor() );
        }
#endif
        return ET;
    };

    auto impl = [&]( auto loads ) {
        for ( const auto &load : loads ) {

            FieldOnNodesRealPtr mate_noeu = nullptr;

            const auto pairing = load->getPairingField();
            if ( pairing && pairing->exists() ) {

                const auto fed = getFED( load );

                // Prepare computing
                CalculPtr calcul = std::make_unique< Calcul >( option );
                calcul->setFiniteElementDescriptor( fed );

                // Set input field
                calcul->addInputField( "PGEOMER", fed->getMesh()->getCoordinates() );
                calcul->addInputField( "PDEPLMR", displ_prev );
                calcul->addInputField( "PDEPLPR", displ_step );
                calcul->addInputField( "PPAIRR", pairing );
                std::cout << "PAIR: " << pairing->getName() << std::endl;
                calcul->addHHOField( model );

                // Add time fields
                calcul->addTimeField( "PINSTMR", time_prev );
                calcul->addTimeField( "PINSTPR", time_prev + time_step );

                // Add output elementary
                calcul->addOutputElementaryTerm( "PVECTUR",
                                                 std::make_shared< ElementaryTermReal >() );

                // Computation
                calcul->compute();
                if ( calcul->hasOutputElementaryTerm( "PVECTUR" ) ) {
                    elemVect->addElementaryTerm(
                        addET( load, calcul->getOutputElementaryTermReal( "PVECTUR" ) ) );
                }
            }
        }
    };

    const auto listOfLoads = _phys_problem->getListOfLoads();

    if ( isHPC ) {
#ifdef ASTER_HAVE_MPI
        impl( listOfLoads->getParallelMechanicalLoadsReal() );
#endif
    } else {
        impl( listOfLoads->getMechanicalLoadsReal() );
    }

    elemVect->build();

    if ( elemVect->hasElementaryTerm() ) {
        return elemVect->assemble( _phys_problem->getDOFNumbering() );
    }

    FieldOnNodesRealPtr vectAsse =
        std::make_shared< FieldOnNodesReal >( _phys_problem->getDOFNumbering() );
    vectAsse->setValues( 0.0 );
    vectAsse->build();
    return vectAsse;
}
