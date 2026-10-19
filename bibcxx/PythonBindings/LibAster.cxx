/**
 * @file LibAster.cxx
 * @brief Création de LibAster
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

#define CODEASTER_IMPORT_ARRAY 1

#include "astercxx.h"

#include "aster_init.h"
#include "aster_numpy.h"
#include "aster_pybind.h"

// Please keep '*Interface.h' files in alphabetical order to ease merging
#include "PythonBindings/AcousticLoadInterface.h"
#include "PythonBindings/AcousticModeResultInterface.h"
#include "PythonBindings/AssemblyMatrixInterface.h"
#include "PythonBindings/AsterToMedWriterInterface.h"
#include "PythonBindings/BaseAssemblyMatrixInterface.h"
#include "PythonBindings/BaseDOFNumberingInterface.h"
#include "PythonBindings/BaseMeshInterface.h"
#include "PythonBindings/BehaviourDefinitionInterface.h"
#include "PythonBindings/BehaviourPropertyInterface.h"
#include "PythonBindings/BucklingModeResultInterface.h"
#include "PythonBindings/CodedMaterialInterface.h"
#include "PythonBindings/CombinedFourierResultInterface.h"
#include "PythonBindings/CommGraphInterface.h"
#include "PythonBindings/ConnectionMeshInterface.h"
#include "PythonBindings/ConstantFieldOnCellsInterface.h"
#include "PythonBindings/ContactComputationInterface.h"
#include "PythonBindings/ContactEnumInterface.h"
#include "PythonBindings/ContactInterface.h"
#include "PythonBindings/ContactNewInterface.h"
#include "PythonBindings/ContactPairingInterface.h"
#include "PythonBindings/ContactParameterInterface.h"
#include "PythonBindings/ContactZoneInterface.h"
#include "PythonBindings/CouplingPairingInterface.h"
#include "PythonBindings/CouplingZonePairingInterface.h"
#include "PythonBindings/CppToFortranGlossaryInterface.h"
#include "PythonBindings/CrackInterface.h"
#include "PythonBindings/CrackShapeInterface.h"
#include "PythonBindings/CreateEnthalpyInterface.h"
#include "PythonBindings/CyclicSymmetryModeInterface.h"
#include "PythonBindings/DOFNumberingInterface.h"
#include "PythonBindings/DataFieldInterface.h"
#include "PythonBindings/DataStructureInterface.h"
#include "PythonBindings/DebugInterface.h"
#include "PythonBindings/DirichletBCInterface.h"
#include "PythonBindings/DiscreteComputationInterface.h"
#include "PythonBindings/DryingResultInterface.h"
#include "PythonBindings/DynamicMacroElementInterface.h"
#include "PythonBindings/ElasticFourierResultInterface.h"
#include "PythonBindings/ElasticResultInterface.h"
#include "PythonBindings/ElementaryCharacteristicsInterface.h"
#include "PythonBindings/ElementaryMatrixInterface.h"
#include "PythonBindings/ElementaryTermInterface.h"
#include "PythonBindings/ElementaryVectorInterface.h"
#include "PythonBindings/EmpiricalModeResultInterface.h"
#include "PythonBindings/EquationNumberingInterface.h"
#include "PythonBindings/ExternalStateVariablesInterface.h"
#include "PythonBindings/ExternalStateVariablesResultInterface.h"
#include "PythonBindings/FiberGeometryInterface.h"
#include "PythonBindings/FieldCharacteristicsInterface.h"
#include "PythonBindings/FieldOnCellsInterface.h"
#include "PythonBindings/FieldOnNodesInterface.h"
#include "PythonBindings/FiniteElementDescriptorInterface.h"
#include "PythonBindings/FluidStructureInteractionInterface.h"
#include "PythonBindings/FluidStructureModalBasisInterface.h"
#include "PythonBindings/FormulaInterface.h"
#include "PythonBindings/FortranInterface.h"
#include "PythonBindings/FullHarmonicAcousticResultInterface.h"
#include "PythonBindings/FullHarmonicResultInterface.h"
#include "PythonBindings/FullResultInterface.h"
#include "PythonBindings/FullTransientResultInterface.h"
#include "PythonBindings/Function2DInterface.h"
#include "PythonBindings/FunctionInterface.h"
#include "PythonBindings/GeneralizedAssemblyMatrixInterface.h"
#include "PythonBindings/GeneralizedAssemblyVectorInterface.h"
#include "PythonBindings/GeneralizedDOFNumberingInterface.h"
#include "PythonBindings/GeneralizedLoadInterface.h"
#include "PythonBindings/GeneralizedModeResultInterface.h"
#include "PythonBindings/GeneralizedModelInterface.h"
#include "PythonBindings/GeneralizedResultInterface.h"
#include "PythonBindings/GenericEnumInterface.h"
#include "PythonBindings/GenericFunctionInterface.h"
#include "PythonBindings/GridInterface.h"
#include "PythonBindings/HHOInterface.h"
#include "PythonBindings/IncompleteMeshInterface.h"
#include "PythonBindings/InterspectralMatrixInterface.h"
#include "PythonBindings/LinearSolverInterface.h"
#include "PythonBindings/ListOfFloatsInterface.h"
#include "PythonBindings/ListOfIntegersInterface.h"
#include "PythonBindings/ListOfLoadsInterface.h"
#include "PythonBindings/LoadResultInterface.h"
#include "PythonBindings/MGISBehaviourInterface.h"
#include "PythonBindings/MaterialFieldInterface.h"
#include "PythonBindings/MaterialInterface.h"
#include "PythonBindings/MatrixStorageInterface.h"
#include "PythonBindings/MatrixToPetscInterface.h"
#include "PythonBindings/MechanicalLoadDescriptionInterface.h"
#include "PythonBindings/MechanicalLoadInterface.h"
#include "PythonBindings/MedCouplingConversionInterface.h"
#include "PythonBindings/MedFamilyInterface.h"
#include "PythonBindings/MedFieldInterface.h"
#include "PythonBindings/MedFileReaderInterface.h"
#include "PythonBindings/MedMeshInterface.h"
#include "PythonBindings/MedToAsterReaderInterface.h"
#include "PythonBindings/MedVectorInterface.h"
#include "PythonBindings/MeshBalancerInterface.h"
#include "PythonBindings/MeshConnectionGraphInterface.h"
#include "PythonBindings/MeshCoordinatesFieldInterface.h"
#include "PythonBindings/MeshEntitiesInterface.h"
#include "PythonBindings/MeshInterface.h"
#include "PythonBindings/MeshPairingInterface.h"
#include "PythonBindings/MeshUtilsInterface.h"
#include "PythonBindings/MeshesMappingInterface.h"
#include "PythonBindings/ModalBasisInterface.h"
#include "PythonBindings/ModeResultInterface.h"
#include "PythonBindings/ModelInterface.h"
#include "PythonBindings/ModelingUtilitiesInterface.h"
#include "PythonBindings/MultipleElasticResultInterface.h"
#include "PythonBindings/NodeInterface.h"
#include "PythonBindings/NonLinearResultInterface.h"
#include "PythonBindings/NumpyInterface.h"
#include "PythonBindings/ObjectBalancerInterface.h"
#include "PythonBindings/ParMetisPartitionerInterface.h"
#include "PythonBindings/ParallelContactFEDescriptorInterface.h"
#include "PythonBindings/ParallelContactNewInterface.h"
#include "PythonBindings/ParallelContactPairingInterface.h"
#include "PythonBindings/ParallelDOFNumberingInterface.h"
#include "PythonBindings/ParallelEquationNumberingInterface.h"
#include "PythonBindings/ParallelFiniteElementDescriptorInterface.h"
#include "PythonBindings/ParallelMechanicalLoadInterface.h"
#include "PythonBindings/ParallelMeshInterface.h"
#include "PythonBindings/ParallelThermalLoadInterface.h"
#include "PythonBindings/PetscApplyFactorOnSubBlocksInterface.h"
#include "PythonBindings/PetscRedistributeInterface.h"
#include "PythonBindings/PhysicalProblemInterface.h"
#include "PythonBindings/PhysicalQuantityInterface.h"
#include "PythonBindings/PhysicalQuantityManagerInterface.h"
#include "PythonBindings/PhysicalSolutionRestitutorInterface.h"
#include "PythonBindings/PhysicsAndModelingsInterface.h"
#include "PythonBindings/PostProcessingInterface.h"
#include "PythonBindings/PrestressingCableInterface.h"
#include "PythonBindings/PtScotchPartitionerInterface.h"
#include "PythonBindings/ResultBalancerInterface.h"
#include "PythonBindings/ResultInterface.h"
#include "PythonBindings/ResultManagerInterface.h"
#include "PythonBindings/ResultNamingInterface.h"
#include "PythonBindings/SetLoggingLevelInterface.h"
#include "PythonBindings/SimpleFieldOnCellsInterface.h"
#include "PythonBindings/SimpleFieldOnNodesInterface.h"
#include "PythonBindings/SkeletonInterface.h"
#include "PythonBindings/StaticMacroElementInterface.h"
#include "PythonBindings/StructureInterfaceInterface.h"
#include "PythonBindings/SuperMeshInterface.h"
#include "PythonBindings/SyntaxSaverInterface.h"
#include "PythonBindings/TableContainerInterface.h"
#include "PythonBindings/TableInterface.h"
#include "PythonBindings/ThermalFourierResultInterface.h"
#include "PythonBindings/ThermalLoadDescriptionInterface.h"
#include "PythonBindings/ThermalLoadInterface.h"
#include "PythonBindings/ThermalResultInterface.h"
#include "PythonBindings/TimesListInterface.h"
#include "PythonBindings/TransientResultInterface.h"
#include "PythonBindings/TurbulentSpectrumInterface.h"
#include "PythonBindings/UnitaryMechanicalLoadInterface.h"
#include "PythonBindings/XfemCrackInterface.h"
#include "Supervis/Exceptions.h"
// Please keep '*Interface.h' files in alphabetical order to ease merging

void *numpyInitialize() {
    import_array();
    return NULL;
}

PYBIND11_MODULE( libaster, mod ) {
    numpyInitialize();
    initAsterModules();

    // hide c++ signatures
    py::options options;
    // options.disable_function_signatures();
    options.disable_enum_members_docstring();

    auto cleanup_callback = []() { jeveux_finalize(); };
    mod.add_object( "_cleanup", py::capsule( cleanup_callback ) );

    // Definition of exceptions, thrown from 'Exceptions.cxx'/uexcep
    createExceptions( mod );

    mod.def( "raiseAsterError", &raiseAsterError, py::arg( "idmess" ) = "VIDE_1",
             py::arg( "valk" ) = VectorString(), py::arg( "vali" ) = VectorLong(),
             py::arg( "valr" ) = VectorReal() );

    // do not sort (compilation error)
    exportGenericEnumToPython( mod );
    exportDataStructureToPython( mod );
    exportDebugToPython( mod );
    exportNumpyToPython( mod );
    exportPhysicalQuantityManagerToPython( mod );
    exportNodeToPython( mod );
    exportMeshEntitiesToPython( mod );
    exportBaseMeshToPython( mod );
    exportMeshToPython( mod );
    exportMedCouplingConversionToPython( mod );
    exportDiscreteComputationToPython( mod );
    exportEquationNumberingToPython( mod );
    exportMatrixStorageToPython( mod );
    exportBaseDOFNumberingToPython( mod );
    exportDOFNumberingToPython( mod );
    exportElementaryCharacteristicsToPython( mod );
    exportFiniteElementDescriptorToPython( mod );
    exportFiberGeometryToPython( mod );
    exportDataFieldToPython( mod );
    exportFieldOnCellsToPython( mod );
    exportFieldOnNodesToPython( mod );
    exportConstantFieldOnCellsToPython( mod );
    exportSimpleFieldOnCellsToPython( mod );
    exportSimpleFieldOnNodesToPython( mod );
    exportTableToPython( mod );
    exportTableContainerToPython( mod );
    exportTimesListToPython( mod );
    exportGeneralizedDOFNumberingToPython( mod );
    exportFluidStructureInteractionToPython( mod );
    exportTurbulentSpectrumToPython( mod );
    exportGenericFunctionToPython( mod );
    exportListOfLoadsToPython( mod );
    exportFunctionToPython( mod );
    exportFormulaToPython( mod );
    exportFortranToPython( mod );
    exportFunction2DToPython( mod );
    exportContactToPython( mod );
    exportContactEnumToPython( mod );
    exportContactParameterToPython( mod );
    exportContactNewToPython( mod );
    exportContactZoneToPython( mod );
    exportMeshPairingToPython( mod );
    exportContactPairingToPython( mod );
    exportContactComputationToPython( mod );
    exportCouplingZonePairingToPython( mod );
    exportCouplingPairingToPython( mod );
    exportBaseAssemblyMatrixToPython( mod );
    exportAssemblyMatrixToPython( mod );
    exportElementaryTermToPython( mod );
    exportElementaryMatrixToPython( mod );
    exportElementaryVectorToPython( mod );
    exportGeneralizedAssemblyMatrixToPython( mod );
    exportGeneralizedAssemblyVectorToPython( mod );
    exportGeneralizedLoadToPython( mod );
    exportInterspectralMatrixToPython( mod );
    exportLinearSolverToPython( mod );
    exportModalBasisToPython( mod );
    exportStructureInterfaceToPython( mod );
    exportAcousticLoadToPython( mod );
    exportDirichletBCToPython( mod );
    exportMechanicalLoadToPython( mod );
    exportMechanicalLoadDescriptionToPython( mod );
    exportUnitaryMechanicalLoadToPython( mod );
    exportPhysicalQuantityToPython( mod );
    exportThermalLoadToPython( mod );
    exportThermalLoadDescriptionToPython( mod );
    exportBehaviourDefinitionToPython( mod );
    exportMaterialToPython( mod );
    exportMaterialFieldToPython( mod );
    exportGridToPython( mod );
    exportMeshesMappingToPython( mod );
    exportSkeletonToPython( mod );
    exportDynamicMacroElementToPython( mod );
    exportStaticMacroElementToPython( mod );
    exportSuperMeshToPython( mod );
    exportCrackShapeToPython( mod );
    exportCrackToPython( mod );
    exportGeneralizedModelToPython( mod );
    exportPhysicsAndModelingsToPython( mod );
    exportModelToPython( mod );
    exportPrestressingCableToPython( mod );
    exportXfemCrackToPython( mod );
    exportResultToPython( mod );
    exportTransientResultToPython( mod );
    exportLoadResultToPython( mod );
    exportThermalResultToPython( mod );
    exportDryingResultToPython( mod );
    exportCombinedFourierResultToPython( mod );
    exportElasticFourierResultToPython( mod );
    exportThermalFourierResultToPython( mod );
    exportMultipleElasticResultToPython( mod );
    exportNonLinearResultToPython( mod );
    exportPhysicalProblemToPython( mod );
    exportCppToFortranGlossaryToPython( mod );
    exportCyclicSymmetryModeToPython( mod );
    exportFullResultToPython( mod );
    exportModeResultToPython( mod );
    exportModeResultComplexToPython( mod );
    exportAcousticModeResultToPython( mod );
    exportBucklingModeResultToPython( mod );
    exportGeneralizedResultToPython( mod );
    exportElasticResultToPython( mod );
    exportMeshCoordinatesFieldToPython( mod );
    exportFullTransientResultToPython( mod );
    exportFullHarmonicResultToPython( mod );
    exportFullHarmonicAcousticResultToPython( mod );
    exportFluidStructureModalBasisToPython( mod );
    exportGeneralizedModeResultToPython( mod );
    exportMGISBehaviourToPython( mod );

#ifdef ASTER_HAVE_MPI
    /* These objects must be declared in ObjectsExt
       as an 'UnavailableObject' for sequential version. */
    exportParallelMeshToPython( mod );
    exportParallelEquationNumberingToPython( mod );
    exportParallelDOFNumberingToPython( mod );
    exportParallelMechanicalLoadToPython( mod );
    exportParallelThermalLoadToPython( mod );
    exportParallelFiniteElementDescriptorToPython( mod );
    exportParallelContactFEDescriptorToPython( mod );
    exportParallelContactNewToPython( mod );
    exportParallelContactPairingToPython( mod );
#endif /* ASTER_HAVE_MPI */
    exportApplyFactorOnSubBlocksToPython( mod );
    exportConnectionMeshToPython( mod );
    exportResultNamingToPython( mod );
    exportListOfFloatsToPython( mod );
    exportListOfIntegersToPython( mod );
    exportEmpiricalModeResultToPython( mod );
    exportExternalStateVariablesToPython( mod );
    exportExternalStateVariablesResultToPython( mod );
    exportCreateEnthalpyToPython( mod );
    exportMatrixToPetscToPython( mod );
    exportBehaviourPropertyToPython( mod );
    exportCodedMaterialToPython( mod );
    exportSetLoggingLevelToPython( mod );
    exportPhysicalSolutionRestitutorToPython( mod );
    exportPostProcessingToPython( mod );
    exportHHOToPython( mod );
#ifdef ASTER_HAVE_MPI
    exportCommGraphToPython( mod );
    exportObjectBalancerToPython( mod );
    exportMeshBalancerToPython( mod );
    exportIncompleteMeshToPython( mod );
#ifdef ASTER_HAVE_SCOTCH
    exportPtScotchPartitionerToPython( mod );
#endif /* ASTER_HAVE_SCOTCH */
#ifdef ASTER_HAVE_PARMETIS
    exportParMetisPartitionerToPython( mod );
#endif /* ASTER_HAVE_PARMETIS */
    exportMeshConnectionGraphToPython( mod );
    exportResultBalancerToPython( mod );
#endif /* ASTER_HAVE_MPI */
    exportPetscRedistributeToPython( mod );
#ifdef ASTER_HAVE_MED
    exportMedFileReaderToPython( mod );
    exportMedFieldToPython( mod );
    exportMedMeshToPython( mod );
    exportMedFamilyToPython( mod );
    exportMedVectorToPython( mod );
#endif /* ASTER_HAVE_MED */
    exportMedToAsterReaderToPython( mod );
    exportAsterToMedWriterToPython( mod );
    exportFieldCharacteristicsToPython( mod );
    exportModelingUtilitiesToPython( mod );
    exportSyntaxSaverToPython( mod );
    exportMeshUtilsToPython( mod );
    exportResultManagerToPython( mod );
};
