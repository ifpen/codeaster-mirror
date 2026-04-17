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


import os.path as osp

from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI
from code_aster.Utilities import SharedTmpdir


CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

# petscInitialize("-start_in_debugger ")

# check ParallelMesh object API
test = CA.TestCase()

# MPI test
rank = MPI.ASTER_COMM_WORLD.Get_rank()
nbproc = MPI.ASTER_COMM_WORLD.Get_size()

test.assertEqual(nbproc, 3)

# from MED format (only this one a ParallelMesh)
mesh = CA.ParallelMesh()
mesh.readMedFile("mesh001b/%d.med" % rank, partitioned=True)


mesh = DEFI_GROUP(
    reuse=mesh,
    MAILLAGE=mesh,
    CREA_GROUP_NO=(_F(NOM="GN" + str(rank), GROUP_NO="EXT_0"),),
    CREA_GROUP_MA=(_F(NOM="GC" + str(rank), GROUP_MA="Cable0"),),
)

test.assertTrue(mesh.isParallel())
test.assertEqual(mesh.getDimension(), 3)

nbNodes = [89, 90, 109]
nbCells = [59, 54, 75]

test.assertEqual(mesh.getNumberOfNodes(), nbNodes[rank])
test.assertEqual(mesh.getNumberOfCells(), nbCells[rank])

# test groups
globalGroupsOfCells = [
    "Cable0",
    "Cable1",
    "Cable2",
    "Cable3",
    "Cable4",
    "Cable5",
    "Cable6",
    "Cable7",
    "Cable8",
    "Press",
    "Encast",
    "Beton",
    "GC0",
    "GC1",
    "GC2",
]
groupsOfCells = [
    [
        "Cable0",
        "Cable1",
        "Cable2",
        "Cable3",
        "Cable4",
        "Cable5",
        "Cable6",
        "Cable7",
        "Cable8",
        "Press",
        "Encast",
        "Beton",
        "GC0",
    ],
    [
        "Cable0",
        "Cable1",
        "Cable2",
        "Cable3",
        "Cable4",
        "Cable5",
        "Cable6",
        "Cable7",
        "Cable8",
        "Press",
        "Beton",
        "GC1",
    ],
    [
        "Cable0",
        "Cable1",
        "Cable2",
        "Cable3",
        "Cable4",
        "Cable5",
        "Cable6",
        "Cable7",
        "Cable8",
        "Press",
        "Encast",
        "Beton",
        "GC2",
    ],
]

test.assertSequenceEqual(sorted(mesh.getGroupsOfCells()), sorted(globalGroupsOfCells))
test.assertSequenceEqual(sorted(mesh.getGroupsOfCells(False)), sorted(globalGroupsOfCells))
test.assertSequenceEqual(sorted(mesh.getGroupsOfCells(True)), sorted(groupsOfCells[rank]))

test.assertTrue(mesh.hasGroupOfCells("Beton"))
test.assertTrue(mesh.hasGroupOfCells("Beton", False))
test.assertTrue(mesh.hasGroupOfCells("Beton", True))
test.assertTrue(mesh.hasGroupOfCells("GC1", False))
test.assertTrue(mesh.hasGroupOfCells("GC" + str(rank), True))
test.assertFalse(mesh.hasGroupOfCells("GC4", True))
test.assertFalse(mesh.hasGroupOfCells("GC4", False))

globalGroupsOfNodes = ["EXT_0", "EXT_1", "EXT_2", "GN0", "GN1", "GN2"]
groupsOfNodes = [
    ["EXT_0", "EXT_1", "EXT_2", "GN0"],
    ["EXT_0", "EXT_1", "EXT_2", "GN1"],
    ["EXT_0", "EXT_1", "EXT_2", "GN2"],
]

test.assertSequenceEqual(sorted(mesh.getGroupsOfNodes()), sorted(globalGroupsOfNodes))
test.assertSequenceEqual(sorted(mesh.getGroupsOfNodes(False)), sorted(globalGroupsOfNodes))
test.assertSequenceEqual(sorted(mesh.getGroupsOfNodes(True)), sorted(groupsOfNodes[rank]))

test.assertTrue(mesh.hasGroupOfNodes("EXT_0"))
test.assertTrue(mesh.hasGroupOfNodes("EXT_0", False))
test.assertTrue(mesh.hasGroupOfNodes("EXT_0", True))
test.assertTrue(mesh.hasGroupOfNodes("GN1", False))
test.assertTrue(mesh.hasGroupOfNodes("GN" + str(rank), True))
test.assertFalse(mesh.hasGroupOfNodes("GN4", True))
test.assertFalse(mesh.hasGroupOfNodes("GN4", False))

# Link between local and global numbering
globalNodesNum = mesh.getNodes(localNumbering=False)
nodesGlobFirst = [4, 44, 0]
test.assertEqual(globalNodesNum[0], nodesGlobFirst[rank])
nodesGlobLast = [97, 97, 95]
test.assertEqual(globalNodesNum[-1], nodesGlobLast[rank])

nbrgField = mesh.getNumberingField()
test.assertEqual(nbrgField.getValues(), globalNodesNum)
with SharedTmpdir("foo") as tmpdir:
    nbrgField.printMedFile(osp.join(tmpdir.path, f"nbrg_{rank}.resu.med"))

# Owner of Nodes
nodesOwner = mesh.getNodesOwner()
test.assertEqual(nodesOwner[1], rank)

ownerField = mesh.getOwnerField()
test.assertEqual(ownerField.getValues(), nodesOwner)
with SharedTmpdir("foo") as tmpdir:
    ownerField.printMedFile(osp.join(tmpdir.path, f"owner_{rank}.resu.med"))

# Node 92 (index is 91) is shared by all meshes (owner is 1)
node92 = [85, 84, 106]
test.assertEqual(globalNodesNum[node92[rank] - 1], 92 - 1)
test.assertEqual(nodesOwner[node92[rank] - 1], 1)

nodesRanks = mesh.getNodesRanks()
print(nodesRanks, flush=True)
test.assertTrue(len(nodesRanks[node92[rank] - 1]) > 1)

# Cells rank
cellsRankRef = [
    [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ],
    [
        0,
        1,
        1,
        0,
        1,
        1,
        0,
        0,
        1,
        1,
        0,
        1,
        1,
        0,
        0,
        1,
        0,
        0,
        1,
        1,
        0,
        1,
        1,
        0,
        1,
        1,
        0,
        1,
        1,
        0,
        0,
        0,
        1,
        1,
        0,
        1,
        1,
        0,
        0,
        0,
        1,
        1,
        0,
        1,
        1,
        0,
        0,
        0,
        1,
        0,
        0,
        1,
        1,
        0,
    ],
    [
        2,
        2,
        0,
        2,
        1,
        1,
        2,
        1,
        2,
        2,
        0,
        2,
        1,
        1,
        2,
        1,
        2,
        2,
        0,
        2,
        1,
        0,
        2,
        1,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        1,
        2,
        2,
        1,
        2,
        2,
        1,
        2,
        2,
        2,
        0,
        2,
        2,
        1,
        1,
        2,
        2,
        1,
        2,
        2,
        2,
        0,
        2,
        2,
        1,
        1,
        2,
        2,
        1,
        2,
        2,
        2,
        0,
        2,
        2,
        1,
        0,
        2,
        2,
        1,
    ],
]
test.assertEqual(len(mesh.getCellsOwner()), mesh.getNumberOfCells())
test.assertSequenceEqual(mesh.getCellsOwner(), cellsRankRef[rank])

cellsRanks = mesh.getCellsRanks()
test.assertSequenceEqual([x[0] for x in cellsRanks], cellsRankRef[rank])

local_map = mesh.getLocalToGlobalCellIds()
test.assertEqual(len(local_map), [59, 54, 75][rank])


def inter(list1, list2):
    return list(set(list1).intersection(list2))


test.assertEqual(mesh.getNumberOfNodes(), len(mesh.getNodes(localNumbering=True)))
test.assertEqual(mesh.getNumberOfNodes(), len(mesh.getNodes(localNumbering=True)))
test.assertEqual(mesh.getNumberOfNodes(), len(mesh.getNodes(localNumbering=False)))

test.assertEqual(mesh.getNumberOfCells(), len(mesh.getCells()))

test.assertSequenceEqual(mesh.getNodes(localNumbering=True), range(mesh.getNumberOfNodes()))
test.assertSequenceEqual(mesh.getNodes(localNumbering=True), range(mesh.getNumberOfNodes()))
test.assertSequenceEqual(mesh.getCells(), range(mesh.getNumberOfCells()))

allNodes = []
innerNodes = []
outerNodes = []

for i in range(mesh.getNumberOfNodes()):
    allNodes.append(i)
    if nodesOwner[i] == rank:
        innerNodes.append(i)
    else:
        outerNodes.append(i)

test.assertTrue(len(inter(innerNodes, outerNodes)) == 0)


test.assertSequenceEqual(sorted(mesh.getNodes(localNumbering=True)), sorted(allNodes))
test.assertSequenceEqual(sorted(mesh.getNodes(localNumbering=True)), sorted(allNodes))
test.assertSequenceEqual(
    sorted(mesh.getNodes(localNumbering=False)), sorted([globalNodesNum[i] for i in allNodes])
)
test.assertSequenceEqual(
    sorted(mesh.getNodes(localNumbering=True, same_rank=False)), sorted(outerNodes)
)
test.assertSequenceEqual(
    sorted(mesh.getNodes(localNumbering=False, same_rank=False)),
    sorted([globalNodesNum[i] for i in outerNodes]),
)
test.assertSequenceEqual(
    sorted(mesh.getNodes(localNumbering=True, same_rank=True)), sorted(innerNodes)
)
test.assertSequenceEqual(
    sorted(mesh.getNodes(localNumbering=False, same_rank=True)),
    sorted([globalNodesNum[i] for i in innerNodes]),
)

test.assertSequenceEqual(
    sorted(mesh.getNodes("Beton", localNumbering=True)), sorted(mesh.getNodes("Beton", True, None))
)
test.assertSequenceEqual(
    sorted(mesh.getNodes("Beton", localNumbering=True)),
    sorted(inter(mesh.getNodes("Beton", localNumbering=True), mesh.getNodes(localNumbering=True))),
)
test.assertSequenceEqual(
    sorted(mesh.getNodes("Beton", localNumbering=True)),
    sorted(inter(mesh.getNodes("Beton", True), mesh.getNodes(localNumbering=True))),
)
test.assertSequenceEqual(
    sorted(mesh.getNodes("Beton", localNumbering=True)),
    sorted(inter(mesh.getNodes("Beton", localNumbering=True), mesh.getNodes(localNumbering=True))),
)
test.assertSequenceEqual(
    sorted(mesh.getNodes("Beton", localNumbering=True)),
    sorted(inter(mesh.getNodes("Beton", True), mesh.getNodes(localNumbering=True))),
)
test.assertSequenceEqual(
    sorted(mesh.getNodes("Beton", False)),
    sorted(inter(mesh.getNodes("Beton", False), mesh.getNodes(localNumbering=False))),
)
test.assertSequenceEqual(
    sorted(mesh.getNodes("Beton", False)),
    sorted(inter(mesh.getNodes("Beton", False), mesh.getNodes(localNumbering=False))),
)
test.assertSequenceEqual(
    sorted(mesh.getNodes("Beton", False)), sorted(inter(mesh.getNodes("Beton", False), allNodes))
)

# rank-owned nodes in local numbering
nodeLoc = [[58, 68, 69], [], [76, 85, 86]]
test.assertSequenceEqual(mesh.getNodesFromCells(["Cable0"], True, True), nodeLoc[rank])
test.assertSequenceEqual(
    mesh.getNodesFromCells("Cable0", False, True),
    sorted([globalNodesNum[i] for i in nodeLoc[rank]]),
)

# nodes in local numbering
nodeLocWithGhosts = [[58, 67, 68, 69], [60, 61], [76, 85, 86, 87]]
test.assertSequenceEqual(mesh.getNodesFromCells("Cable0", True, None), nodeLocWithGhosts[rank])
test.assertSequenceEqual(
    mesh.getNodesFromCells("Cable0", False, None),
    sorted([globalNodesNum[i] for i in nodeLocWithGhosts[rank]]),
)

# test add group
mesh.setGroupOfNodes("TEST_GNO", mesh.getNodesFromCells("Cable0", True, None), localNumbering=True)
test.assertSequenceEqual(
    mesh.getNodesFromCells("Cable0", True, None), mesh.getNodes("TEST_GNO", True)
)
mesh.setGroupOfNodes("TEST_GN1", [0, 12, 48, 18, 87], localNumbering=False)
nodeglob = [[87], [], [0, 12, 18, 48, 87]]
test.assertSequenceEqual(nodeglob[rank], mesh.getNodes("TEST_GN1", False))

# test global to local mapping
local_map = mesh.getLocalToGlobalNodeIds()
global_map = mesh.getGlobalToLocalNodeIds()

for i in range(len(local_map)):
    test.assertEqual(i, global_map[local_map[i]])

new_mesh = mesh.refine(2)

# test mesh builder
square = CA.ParallelMesh.buildSquare(refine=4, deterministic=True)
square2 = square.restrict(["TOP", "LEFT"])
test.assertEqual(square.getNumberOfNodes(), [114, 115, 120][rank])
test.assertEqual(square2.getNumberOfNodes(), [22, 13, 0][rank])
test.assertEqual(mesh.getLocalToGlobalNodeIds()[0], [4, 44, 0][rank])
test.assertEqual(
    CA.ParallelMesh.buildCube(refine=4, deterministic=True).getNumberOfNodes(),
    [1955, 2151, 2205][rank],
)
test.assertEqual(
    CA.ParallelMesh.buildDisk(refine=5, deterministic=True).getNumberOfNodes(),
    [5546, 5426, 6027][rank],
)
test.assertEqual(
    CA.ParallelMesh.buildCylinder(refine=3, deterministic=True).getNumberOfNodes(),
    [3321, 3321, 3825][rank],
)

field_r = CREA_CHAMP(
    TYPE_CHAM="NOEU_DEPL_R",
    OPERATION="AFFE",
    MAILLAGE=square2,
    AFFE=_F(TOUT="OUI", NOM_CMP=("DX", "DY"), VALE=(1.0, 2.0)),
)


# test ghosts
square = CA.ParallelMesh.buildSquare(refine=2, deterministic=True, ghost=2)

test.assertEqual(
    square.getNodes(localNumbering=True, same_rank=False),
    [
        [2, 3, 6, 7, 10, 11, 13, 14, 15, 17, 18, 19],
        [2, 3, 4, 7, 8, 9, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23],
        [4, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24],
    ][rank],
)

test.assertEqual(
    square.getCellsOwner(),
    [
        [0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2, 2],
        [1, 1, 1, 0, 0, 1, 1, 2, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 2, 2, 2],
        [0, 0, 2, 2, 2, 2, 0, 1, 0, 1, 1, 1, 0, 2, 2, 2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
    ][rank],
)

test.assertEqual(
    square.getLocalToGlobalCellIds(),
    [
        [0, 1, 2, 17, 3, 4, 5, 6, 25, 7, 8, 21, 9, 10, 19, 11, 12, 23, 13, 29, 30],
        [14, 15, 16, 0, 2, 17, 18, 28, 10, 19, 20, 11, 12, 7, 8, 21, 22, 9, 23, 24, 29, 30, 31],
        [
            5,
            6,
            25,
            26,
            27,
            28,
            3,
            15,
            4,
            16,
            23,
            24,
            13,
            29,
            30,
            31,
            7,
            8,
            21,
            22,
            9,
            10,
            19,
            20,
            11,
            12,
        ],
    ][rank],
)

test.printSummary()

CA.close()
