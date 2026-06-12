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
from code_aster.Utilities import petscInitialize
from code_aster.CA import MPI

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

# force PETSc to start before solves for testing purpose only - no need in regular study
petscInitialize("-ksp_view -log_view -ksp_monitor")

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
print("Nb procs", MPI.ASTER_COMM_WORLD.Get_size())
print("Rank", MPI.ASTER_COMM_WORLD.Get_rank())

pMesh2 = CA.ParallelMesh()
pMesh2.readMedFile("mesh004a/%d.med" % rank, partitioned=True)
pMesh2 = DEFI_GROUP(reuse=pMesh2, MAILLAGE=pMesh2, CREA_GROUP_NO=_F(TOUT_GROUP_MA="OUI"))
del pMesh2

pMesh = CA.ParallelMesh()
pMesh.readMedFile("mesh004a/%d.med" % rank, partitioned=True)
pMesh.debugPrint(rank + 30)

model = CA.Model(pMesh)
test.assertEqual(model.getType(), "MODELE_SDASTER")
model.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)
model.build()

testMesh = model.getMesh()
test.assertEqual(testMesh.getType(), "MAILLAGE_P")

model.debugPrint(rank + 30)

acier = DEFI_MATERIAU(ELAS=_F(E=2.0e11, NU=0.3))
acier.debugPrint(8)

affectMat = CA.MaterialField(pMesh)

testMesh2 = affectMat.getMesh()
test.assertEqual(testMesh2.getType(), "MAILLAGE_P")

affectMat.addMaterialOnMesh(acier)
affectMat.build()

charCine = CA.MechanicalDirichletBC(model)
charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Dx, 0.0, "COTE_B")
charCine.build()

study = CA.PhysicalProblem(model, affectMat)
study.addDirichletBC(charCine)
nume = NUME_DDL(MODELE=model)
study.setDOFNumbering(nume)
# Without calling computeListOfLoads, the Dirichlet BC DOFs should still be available
if rank == 0 or rank == 1:
    test.assertTrue(1 in study.getDirichletBCDOFs(), "No Dirichlet BC DOFs found")
# Calling computeListOfLoads should not affect the availability of Dirichlet BC DOFs
if rank == 0 or rank == 1:
    study.computeListOfLoads()
    test.assertTrue(1 in study.getDirichletBCDOFs(), "No Dirichlet BC DOFs found")
#
dComputation = CA.DiscreteComputation(study)
matr_elem = dComputation.getElasticStiffnessMatrix()

listLoads = study.getListOfLoads()

numeDDL = CA.ParallelDOFNumbering()
numeDDL.computeNumbering([matr_elem])
numeDDL.debugPrint(rank + 30)

matrAsse = CA.AssemblyMatrixDisplacementReal()
matrAsse.setDOFNumbering(numeDDL)
matrAsse.assemble(matr_elem, charCine)
matrAsse.debugPrint(rank + 30)

ccid = matrAsse.getDirichletBCDOFs()
bcNb = {0: 142, 1: 148, 2: 0, 3: 0}
test.assertEqual(sum(ccid), bcNb[rank])
test.assertEqual(len(ccid), numeDDL.getNumberOfDOFs(local=True) + 1)

vec = CA.FieldOnNodesReal(model)
vec.setValues(1.0)
study.zeroDirichletBCDOFs(vec)
test.assertEqual(sum(vec.getValues()), numeDDL.getNumberOfDOFs(local=True) - sum(ccid[:-1]))
vec.setValues(1.0)
matrAsse.zeroDirichletBCDOFs(vec)
test.assertEqual(sum(vec.getValues()), numeDDL.getNumberOfDOFs(local=True) - sum(ccid[:-1]))

retour = matrAsse.getDOFNumbering()
test.assertEqual(retour.isParallel(), True)

# tests on numbering of DOFs
physicalRows = numeDDL.getPhysicalDOFs(local=True)
test.assertListEqual(physicalRows, list(range(3 * len(pMesh.getNodes(localNumbering=True)))))
multipliersRows = numeDDL.getLagrangeDOFs(local=True)
test.assertListEqual(multipliersRows, [])
test.assertFalse(numeDDL.useLagrangeDOF())
test.assertFalse(numeDDL.useSingleLagrangeDOF())
test.assertEqual(numeDDL.getComponents(), ["DX", "DY", "DZ"])
test.assertEqual(numeDDL.getComponentFromNode(0, local=True), ["DX", "DY", "DZ"])
test.assertEqual(numeDDL.getNodeFromDOF(0, local=True), 0)
test.assertTrue(numeDDL.isPhysicalDOF(0, local=True))
test.assertEqual(numeDDL.getNumberOfDOFs(local=True), 3 * len(pMesh.getNodes(localNumbering=True)))
test.assertEqual(numeDDL.getNumberOfDOFs(local=False), 3993)
test.assertEqual(numeDDL.getPhysicalQuantity(), "DEPL_R")
nnodes = pMesh.getNumberOfNodes()
test.assertEqual(
    numeDDL.getNodeAndComponentFromDOF(local=True)[::3], [(i, "DX") for i in range(nnodes)]
)
l2G = pMesh.getLocalToGlobalNodeIds()
test.assertEqual(
    numeDDL.getNodeAndComponentFromDOF(local=False)[::3], [(l2G[i], "DX") for i in range(nnodes)]
)

dof_glo = [3003, 1002, 1014, 15]
node_glo = [12, 7, 47, 11]
dof_loc = [18, 3, 6, 9]
node_loc = [6, 1, 2, 3]
test.assertEqual(
    numeDDL.getNodeAndComponentFromDOF(dof_glo[rank], local=False), (node_glo[rank], "DX")
)
test.assertEqual(
    numeDDL.getNodeAndComponentFromDOF(dof_loc[rank], local=True), (node_loc[rank], "DX")
)

test.printSummary()

CA.close()
