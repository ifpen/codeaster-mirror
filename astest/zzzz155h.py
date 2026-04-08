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

import numpy as np

from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()


def buildCompleteFieldOnCells(field):
    """
    Build complete (over processes) field on cells
    Arguments:
        field (SimpleFieldOnCells): field to complete

    Returns:
        list: list of lists containing all field values over processes on global numbering
    """
    field.updateValuePointers()
    mesh = field.getMesh()
    lTGC = mesh.getLocalToGlobalCellIds()
    maxCells = max(lTGC)
    maxCells = MPI.ASTER_COMM_WORLD.allreduce(maxCells, MPI.MAX)

    innerCellsSet = set(mesh.getInnerCells())
    values = field.toNumpy()[0]
    nbCmp = field.getNumberOfComponents()
    nbPt = field.getMaxNumberOfPoints()
    nbSPt = field.getMaxNumberOfSubPoints()
    cWComps = field.getCellsWithValues()

    completeSief = np.zeros((maxCells + 1, nbCmp * nbPt * nbSPt))
    cmpt = 0
    for idCell in cWComps:
        if idCell in innerCellsSet:
            globCellId = lTGC[idCell]
            toAdd = []
            for iCmp in range(nbCmp):
                toAdd.extend(list(values[cmpt * nbPt * nbSPt : (cmpt + 1) * nbPt * nbSPt, iCmp]))
            completeSief[globCellId] = np.array(toAdd)
        cmpt += 1

    return MPI.ASTER_COMM_WORLD.allreduce(completeSief, MPI.SUM)


def buildCompleteFieldOnNodes(field):
    """
    Build complete (over processes) field on nodes
    Arguments:
        field (SimpleFieldOnNodes): field to complete

    Returns:
        list: list of lists containing all field values over processes on global numbering
    """
    field.updateValuePointers()
    mesh = field.getMesh()
    lTGN = mesh.getLocalToGlobalNodeIds()
    maxNodes = max(lTGN)
    maxNodes = MPI.ASTER_COMM_WORLD.allreduce(maxNodes, MPI.MAX)

    innerNodesSet = set(mesh.getInnerNodes())
    values, mask = field.getValues()
    nbNode = field.getNumberOfNodes()
    nbCmp = field.getNumberOfComponents()

    completeField = np.zeros((maxNodes + 1, nbCmp))
    cmpt = 0
    for idNode in range(nbNode):
        if idNode in innerNodesSet:
            globNodeId = lTGN[idNode]
            toAdd = []
            for iCmp in range(nbCmp):
                toAdd.append(values[idNode, iCmp])
            completeField[globNodeId] = np.array(toAdd)
        cmpt += 1

    return MPI.ASTER_COMM_WORLD.allreduce(completeField, MPI.SUM)


rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()

filename = "fort.20"

mesh = CA.ParallelMesh()
mesh.readMedFile(filename)

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))

cl = AFFE_CHAR_CINE(
    MODELE=model,
    SYNTAXE="OUI",
    MECA_IMPO=(_F(GROUP_MA="BAS", DX=0.0, DY=0.0, DZ=0.0), _F(GROUP_MA="HAUT", DZ=1.0)),
)


###########################################
# First MECA_STATIQUE with AFFE_CHAR_CINE #
###########################################
resu = MECA_STATIQUE(
    MODELE=model, CHAM_MATER=mater, EXCIT=_F(CHARGE=cl), INST=0.0, SOLVEUR=_F(METHODE="MUMPS")
)

# Balancing ElasticResult
bMesh = CA.MeshBalancer()
bMesh.buildFromBaseMesh(mesh)
outResu = None
if rank == 0:
    outResu = CA.applyBalancingStrategy(
        resu, [1, 2, 3, 4, 9, 10, 11, 12, 13, 14, 15, 16, 20, 26, 33, 34, 35, 36, 41, 43]
    )
elif rank == 1:
    outResu = CA.applyBalancingStrategy(
        resu, [17, 18, 29, 30, 37, 38, 39, 45, 46, 47, 49, 50, 51, 52, 57, 59, 61, 62, 63, 64]
    )
elif rank == 2:
    outResu = CA.applyBalancingStrategy(
        resu, [6, 7, 19, 21, 22, 23, 24, 25, 27, 28, 40, 42, 44, 48, 53, 54, 55, 56, 58, 60]
    )
elif rank == 3:
    outResu = CA.applyBalancingStrategy(
        resu, [5, 8, 31, 32, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80]
    )

del bMesh

oldSiefElga = resu.getField("SIEF_ELGA", 1).toSimpleFieldOnCells()
initSief = buildCompleteFieldOnCells(oldSiefElga)

oldDepl = resu.getField("DEPL", 1).toSimpleFieldOnNodes()
initDepl = buildCompleteFieldOnNodes(oldDepl)

newSiefElga = outResu.getField("SIEF_ELGA", 1).toSimpleFieldOnCells()
newSief = buildCompleteFieldOnCells(newSiefElga)

depl = outResu.getField("DEPL", 1).toSimpleFieldOnNodes()
newDepl = buildCompleteFieldOnNodes(depl)

# Compare SIEF_ELGA (reference and balanced)
for initArray, newArray in zip(initSief, newSief):
    for initVal, newVal in zip(initArray, newArray):
        test.assertEqual(initVal, newVal)

# Compare DEPL (reference and balanced)
for initArray, newArray in zip(initDepl, newDepl):
    for initVal, newVal in zip(initArray, newArray):
        test.assertEqual(initVal, newVal)

test.assertEqual(outResu.getType(), "EVOL_ELAS")
test.assertTrue(isinstance(outResu, CA.ElasticResult))


#######################################
# MECA_STATIQUE with balanced objects #
#######################################
bModel = outResu.getModel()
bMater = outResu.getMaterialField()
bLoads = outResu.getListOfLoads(1)
bBCLoad = bLoads.getDirichletBCs()
assert len(bBCLoad) == 1

resuBis = MECA_STATIQUE(
    MODELE=bModel,
    CHAM_MATER=bMater,
    EXCIT=_F(CHARGE=bBCLoad[0]),
    INST=0.0,
    SOLVEUR=_F(METHODE="MUMPS"),
)

newSiefElga = resuBis.getField("SIEF_ELGA", 1).toSimpleFieldOnCells()
newSief = buildCompleteFieldOnCells(newSiefElga)

depl = resuBis.getField("DEPL", 1).toSimpleFieldOnNodes()
newDepl = buildCompleteFieldOnNodes(depl)

# Compare SIEF_ELGA (reference and with balanced objects)
for initArray, newArray in zip(initSief, newSief):
    for initVal, newVal in zip(initArray, newArray):
        test.assertAlmostEqual(initVal, newVal, delta=5.0e-10)

# Compare DEPL (reference and balanced)
for initArray, newArray in zip(initDepl, newDepl):
    for initVal, newVal in zip(initArray, newArray):
        test.assertAlmostEqual(initVal, newVal, delta=1.0e-15)


########################################################
# MECA_STATIQUE with AFFE_CHAR_CINE and AFFE_CHAR_MECA #
########################################################
clCine = AFFE_CHAR_CINE(
    MODELE=model, MECA_IMPO=(_F(GROUP_MA="BAS", DX=0.0, DY=0.0, DZ=0.0),), SYNTAXE="OUI"
)

clMeca = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="HAUT", DZ=1.0),), SYNTAXE="OUI")

resu2 = MECA_STATIQUE(
    MODELE=model,
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=clCine), _F(CHARGE=clMeca)),
    INST=0.0,
    SOLVEUR=_F(METHODE="MUMPS"),
)

# Balancing ElasticResult
bMesh = CA.MeshBalancer()
bMesh.buildFromBaseMesh(mesh)
outResu = None
if rank == 0:
    outResu = CA.applyBalancingStrategy(
        resu2, [1, 2, 3, 4, 9, 10, 11, 12, 13, 14, 15, 16, 20, 26, 33, 34, 35, 36, 41, 43]
    )
elif rank == 1:
    outResu = CA.applyBalancingStrategy(
        resu2, [17, 18, 29, 30, 37, 38, 39, 45, 46, 47, 49, 50, 51, 52, 57, 59, 61, 62, 63, 64]
    )
elif rank == 2:
    outResu = CA.applyBalancingStrategy(
        resu2, [6, 7, 19, 21, 22, 23, 24, 25, 27, 28, 40, 42, 44, 48, 53, 54, 55, 56, 58, 60]
    )
elif rank == 3:
    outResu = CA.applyBalancingStrategy(
        resu2, [5, 8, 31, 32, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80]
    )

del bMesh

bModel = outResu.getModel()
bMater = outResu.getMaterialField()
bLoads = outResu.getListOfLoads(1)
bBCLoad = bLoads.getDirichletBCs()
assert len(bBCLoad) == 1
bMLLoad = bLoads.getParallelMechanicalLoadsReal()
assert len(bMLLoad) == 1

resuBis = MECA_STATIQUE(
    MODELE=bModel,
    CHAM_MATER=bMater,
    EXCIT=(_F(CHARGE=bBCLoad[0]), _F(CHARGE=bMLLoad[0])),
    INST=0.0,
    SOLVEUR=_F(METHODE="MUMPS"),
)

oldSiefElga = resu2.getField("SIEF_ELGA", 1).toSimpleFieldOnCells()
initSief = buildCompleteFieldOnCells(oldSiefElga)

oldDepl = resu2.getField("DEPL", 1).toSimpleFieldOnNodes()
initDepl = buildCompleteFieldOnNodes(oldDepl)

newSiefElga = resuBis.getField("SIEF_ELGA", 1).toSimpleFieldOnCells()
newSief = buildCompleteFieldOnCells(newSiefElga)

depl = resuBis.getField("DEPL", 1).toSimpleFieldOnNodes()
newDepl = buildCompleteFieldOnNodes(depl)

# Compare SIEF_ELGA (reference and with balanced objects)
for initArray, newArray in zip(initSief, newSief):
    for initVal, newVal in zip(initArray, newArray):
        test.assertAlmostEqual(initVal, newVal, delta=5.0e-10)

# Compare DEPL (reference and balanced)
for initArray, newArray in zip(initDepl, newDepl):
    for initVal, newVal in zip(initArray, newArray):
        test.assertAlmostEqual(initVal, newVal, delta=1.0e-15)

mesh2 = CREA_MAILLAGE(MAILLAGE=mesh, CREA_POI1=_F(NOM_GROUP_MA="POI1", GROUP_NO="Noeuds_poi1"))
globalCellIds = mesh2.getLocalToGlobalCellIds()
if rank == 0:
    test.assertEqual(
        globalCellIds,
        [
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
            31,
            32,
            33,
            34,
            35,
            36,
            37,
            38,
            39,
            40,
            41,
            42,
            43,
            44,
            45,
            46,
            47,
            48,
            49,
            50,
            51,
            52,
            53,
            54,
            -9223372036854775807,
            -9223372036854775807,
            55,
            -9223372036854775807,
            -9223372036854775807,
            56,
            -9223372036854775807,
            -9223372036854775807,
        ],
    )
elif rank == 1:
    test.assertEqual(
        globalCellIds,
        [
            57,
            -9223372036854775807,
            58,
            59,
            60,
            61,
            -9223372036854775807,
            62,
            63,
            64,
            65,
            66,
            67,
            68,
            -9223372036854775807,
            -9223372036854775807,
            69,
            70,
            71,
            72,
            73,
            74,
            75,
            76,
            -9223372036854775807,
            77,
            -9223372036854775807,
            78,
            -9223372036854775807,
            79,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            80,
            81,
            82,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            83,
            84,
            85,
            -9223372036854775807,
            -9223372036854775807,
            86,
            87,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            88,
            89,
            90,
            -9223372036854775807,
            91,
            92,
            93,
            94,
            -9223372036854775807,
            95,
            96,
            97,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
        ],
    )
elif rank == 2:
    test.assertEqual(
        globalCellIds,
        [
            98,
            99,
            100,
            101,
            102,
            103,
            104,
            105,
            106,
            107,
            108,
            -9223372036854775807,
            109,
            -9223372036854775807,
            110,
            111,
            112,
            113,
            114,
            -9223372036854775807,
            -9223372036854775807,
            115,
            116,
            117,
            118,
            119,
            120,
            121,
            122,
            123,
            124,
            125,
            126,
            127,
            128,
            129,
            130,
            131,
            132,
            133,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            134,
            -9223372036854775807,
            135,
            -9223372036854775807,
            136,
            -9223372036854775807,
            -9223372036854775807,
            137,
            138,
            139,
            140,
            141,
            142,
            -9223372036854775807,
            143,
            -9223372036854775807,
            144,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            145,
            146,
            147,
            148,
            -9223372036854775807,
            149,
            150,
            -9223372036854775807,
            151,
            152,
            -9223372036854775807,
            153,
            -9223372036854775807,
        ],
    )
elif rank == 3:
    test.assertEqual(
        globalCellIds,
        [
            -9223372036854775807,
            154,
            155,
            156,
            157,
            -9223372036854775807,
            158,
            159,
            160,
            161,
            162,
            163,
            -9223372036854775807,
            -9223372036854775807,
            164,
            -9223372036854775807,
            165,
            166,
            -9223372036854775807,
            -9223372036854775807,
            167,
            -9223372036854775807,
            168,
            -9223372036854775807,
            169,
            170,
            -9223372036854775807,
            171,
            172,
            173,
            -9223372036854775807,
            174,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            175,
            -9223372036854775807,
            176,
            -9223372036854775807,
            177,
            -9223372036854775807,
            178,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            179,
            -9223372036854775807,
            180,
            -9223372036854775807,
            181,
            -9223372036854775807,
            182,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            183,
            184,
            185,
            -9223372036854775807,
            186,
            -9223372036854775807,
            187,
            -9223372036854775807,
            188,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            -9223372036854775807,
            189,
            190,
            -9223372036854775807,
            191,
        ],
    )

FIN()
