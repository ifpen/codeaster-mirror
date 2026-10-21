# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------


from code_aster.Commands import *
from code_aster import CA
import libaster
import tempfile
import numpy as np
from code_aster.Utilities import PETSc

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

refinement = 1
monMaillage = CA.Mesh.buildCube(refine=refinement)

coords = monMaillage.getCoordinates()
npcoords = coords.toNumpy()
c0val = coords[0][0]
npcoords[0][0] = 1.2345
test.assertAlmostEqual(coords[0][0], 1.2345, msg="coords.toNumpy")
npcoords[0][0] = c0val

monModel = CA.Model(monMaillage)
monModel.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)
monModel.build()
test.assertEqual(monModel.getType(), "MODELE_SDASTER")

YOUNG = 200000.0
POISSON = 0.3

Kinv = 3.2841e-4
Kv = 1.0 / Kinv
SY = 437.0
Rinf = 758.0
Qzer = 758.0 - 437.0
Qinf = Qzer + 100.0
b = 2.3
C1inf = 63767.0 / 2.0
C2inf = 63767.0 / 2.0
Gam1 = 341.0
Gam2 = 341.0
C_Pa = 1.0e6

acier = DEFI_MATERIAU(
    ELAS=_F(E=YOUNG, NU=POISSON),
    VISCOCHAB=_F(
        K=SY * C_Pa,
        B=b,
        MU=10,
        Q_M=Qinf * C_Pa,
        Q_0=Qzer * C_Pa,
        C1=C1inf * C_Pa,
        C2=C2inf * C_Pa,
        G1_0=Gam1,
        G2_0=Gam2,
        K_0=Kv * C_Pa,
        N=11,
        A_K=1.0,
    ),
)
# acier.debugPrint(6)
test.assertEqual(acier.getType(), "MATER_SDASTER")

affectMat = CA.MaterialField(monMaillage)
affectMat.addMaterialOnMesh(acier)
affectMat.addMaterialOnGroupOfCells(acier, ["TOP", "BOTTOM"])
affectMat.build()
test.assertEqual(affectMat.getType(), "CHAM_MATER")

imposedDof1 = CA.DisplacementReal()
imposedDof1.setValue(CA.PhysicalQuantityComponent.Dx, 0.0)
imposedDof1.setValue(CA.PhysicalQuantityComponent.Dy, 0.0)
imposedDof1.setValue(CA.PhysicalQuantityComponent.Dz, 0.0)
CharMeca1 = CA.ImposedDisplacementReal(monModel)
CharMeca1.setValue(imposedDof1, "BOTTOM")
CharMeca1.build()
test.assertEqual(CharMeca1.getType(), "CHAR_MECA")

imposedPres1 = CA.PressureReal()
imposedPres1.setValue(CA.PhysicalQuantityComponent.Pres, 1000.0)
CharMeca2 = CA.DistributedPressureReal(monModel)
CharMeca2.setValue(imposedPres1, "TOP")
CharMeca2.build()
test.assertEqual(CharMeca2.getType(), "CHAR_MECA")


study = CA.PhysicalProblem(monModel, affectMat)
study.addLoad(CharMeca1)
study.addLoad(CharMeca2)
listLoads = study.getListOfLoads()
study.computeDOFNumbering()
dComputation = CA.DiscreteComputation(study)
# compute Neumann
retour = dComputation.getNeumannForces(1)
matr_elem = dComputation.getLinearStiffnessMatrix(with_dual=False)
matr_elem_dual = dComputation.getDualStiffnessMatrix()

test.assertEqual(matr_elem.getType(), "MATR_ELEM_DEPL_R")

monSolver = CA.MumpsSolver()

numeDDL = CA.DOFNumbering()
numeDDL.computeNumbering([matr_elem, matr_elem_dual])
test.assertEqual(numeDDL.getType(), "NUME_DDL_SDASTER")
# test.assertFalse(numeDDL.hasDirichletBC())

matrAsse = CA.AssemblyMatrixDisplacementReal()
matrAsse.setDOFNumbering(numeDDL)
test.assertFalse(listLoads.hasDirichletBC())
matrAsse.assemble([matr_elem, matr_elem_dual])

ccid = matrAsse.getDirichletBCDOFs()
test.assertEqual(sum(ccid), 0)
test.assertEqual(len(ccid), numeDDL.getNumberOfDOFs() + 1)

x = matrAsse.getValuesWithDescription()
test.assertTrue("numpy" in str(type(x[0])))

# test setValues
# -----------------
values, idx, jdx, neq = matrAsse.getValuesWithDescription()
K1 = matrAsse.toNumpy()

neq = K1.shape[0]
matrAsse.setValues([0, 1], [0, 1], [1.0, 1.0])
K2 = matrAsse.toNumpy()
test.assertAlmostEqual(np.linalg.norm(K2), np.sqrt(2))

matrAsse.setValues(idx.tolist(), jdx.tolist(), values.tolist())
K3 = matrAsse.toNumpy()
test.assertEqual(np.linalg.norm(K1 - K3), 0)

A = matrAsse.toPetsc()

v = PETSc.Viewer()
A.view(v)
v = PETSc.Viewer().createASCII("test.txt")
v.pushFormat(PETSc.Viewer.Format.ASCII_MATLAB)
A.view(v)

test.assertEqual(matrAsse.getType(), "MATR_ASSE_DEPL_R")
test.assertTrue(isinstance(matrAsse, CA.AssemblyMatrixDisplacementReal))
monSolver.factorize(matrAsse)
matrfact = monSolver.getMatrix()
test.assertTrue(matrAsse == matrfact)
test.assertEqual(matrfact.getType(), "MATR_ASSE_DEPL_R")
test.assertTrue(isinstance(matrfact, CA.AssemblyMatrixDisplacementReal))
precond = monSolver.getPrecondMatrix()
test.assertEqual(precond.getType(), "MATR_ASSE_DEPL_R")
test.assertTrue(isinstance(precond, CA.AssemblyMatrixDisplacementReal))

vcine = dComputation.getDirichletBC(0.0)
resu = monSolver.solve(retour, vcine)

y, _ = resu.getValuesWithDescription()
test.assertEqual(len(y), len(resu.getValues()))

resu2 = resu.toSimpleFieldOnNodes()
resu2.updateValuePointers()
# point at (1., 1., 1.)
corner = monMaillage.getNumberOfNodes() - 1
test.assertAlmostEqual(resu2[corner, "DX"], 0.000757555469653289)

resu.printMedFile("fort.med")

# test setValues + solve
# ----------------------
matrAsse.setValues(idx.tolist(), jdx.tolist(), [10 * v for v in values])
monSolver.factorize(matrAsse)
resu = monSolver.solve(retour)
resu2 = resu.toSimpleFieldOnNodes()
resu2.updateValuePointers()
test.assertAlmostEqual(resu2[corner, 0], 0.000757555469653289 / 10.0)

# Test conversions
fnodes0 = CA.SimpleFieldOnNodesReal(monMaillage, "DEPL_R", ["DX", "DY", "DZ"], True)
nval, nma = fnodes0.toNumpy()
nma[:, :] = True
nval[:, 1] = 3.6
full_nodes0 = fnodes0.toFieldOnNodes()
test.assertAlmostEqual(
    full_nodes0.norm("NORM_2"), 18.706148721743872, places=6, msg="Conversion SimpleFieldOnNodes"
)

fcells0 = CA.SimpleFieldOnCellsReal(monMaillage, "ELEM", "DEPL_R", ["DX", "DY", "DZ"], 1, 1, True)
cval, cma, _, _ = fcells0.toNumpy()
cma[:, :] = True
cval[:, 1] = 2.5
full_nodes1 = fcells0.toFieldOnNodes()
test.assertAlmostEqual(
    full_nodes1.norm("NORM_2"), 12.99038105676658, places=6, msg="Conversion SimpleFieldOnCells"
)

### medcoupling conversion

medmesh = resu.getMesh().createMedCouplingMesh()
medfield = resu2.toMedFileField1TS(medmesh, "DEPL")
medcfield = resu2.toMedCouplingField(medmesh.getMeshAtLevel(0), "DEPL")

resu3 = CA.SimpleFieldOnNodesReal.fromMedCouplingField(medcfield, monMaillage)

with tempfile.NamedTemporaryFile(prefix="test_", suffix=".rmed", mode="w", delete=True) as f:
    medmesh.write(f.name, 2)
    medfield.write(f.name, 0)

test.assertEqual(medfield.getName(), "DEPL")
test.assertEqual(medfield.getNumberOfComponents(), 3)

# To be sure that vcine is Permanent #30689
libaster.deleteTemporaryObjects()
vcine.updateValuePointers()

# check other solvers attributes
test.assertEqual(monSolver.getSolverName(), "MUMPS")
test.assertTrue(monSolver.supportParallelMesh(), "support of ParallelMesh")

monSolver = CA.LdltSolver()
test.assertEqual(monSolver.getSolverName(), "LDLT")
test.assertFalse(monSolver.supportParallelMesh(), "support of ParallelMesh")

monSolver = CA.MultFrontSolver()
test.assertEqual(monSolver.getSolverName(), "MULT_FRONT")
test.assertFalse(monSolver.supportParallelMesh(), "support of ParallelMesh")

monSolver = CA.PetscSolver()
test.assertEqual(monSolver.getSolverName(), "PETSC")
test.assertTrue(monSolver.supportParallelMesh(), "support of ParallelMesh")

monSolver = CA.GcpcSolver()
test.assertEqual(monSolver.getSolverName(), "GCPC")
test.assertFalse(monSolver.supportParallelMesh(), "support of ParallelMesh")

# Test Component
elas = CA.ElasticResult()
elas.allocate(1)
elas.setField(resu, "DEPL", para="NUME_ORDRE", value=1)
elas.setTime(1.0, 1)
elas.setModel(monModel)
elas.setMaterialField(affectMat)

elas = CALC_CHAMP(reuse=elas, RESULTAT=elas, CONTRAINTE="SIEF_ELGA")
chcoor = CALC_CHAM_ELEM(MODELE=monModel, OPTION="COOR_ELGA")
coor = chcoor.toSimpleFieldOnCells()

sief = elas.getField("SIEF_ELGA", 1).toSimpleFieldOnCells()
sixx = sief.SIXX

# where sixx is defined (=3D)
cells = sixx.getCells()
cells = cells[(sixx._idx >= 0).sum(axis=1) != 0]
test.assertTrue(np.all(cells == np.array(monMaillage.getCells("VOLUME"))))

weight = coor.W.onSupportOf(sixx)
integr = (sixx * weight).getValuesByCells().sum(axis=1)
print(integr[cells])

test.printSummary()

FIN()
