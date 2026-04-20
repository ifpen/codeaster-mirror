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
import re
from glob import glob

import medcoupling as mc

from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI

CA.init("--test")

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
nbproc = MPI.ASTER_COMM_WORLD.Get_size()

if nbproc > 1:
    is_parallel = True
else:
    is_parallel = False

# list of meshes with error
failure = [
    # for mesh002s
    "sdll06a.mmed",
    # for mesh002t
    "ssls134a.med",
    "ssna107c.mmed",
    "ssnl121b.mmed",
    "ssnp05a.mmed",
    "ssnp107a.mmed",
    # for mesh002u
    "ssnv214a.med",
]

# get list of meshes from .export
export = glob("*.export")
assert len(export) == 1, "expecting exactly one export file!"

re_mesh = re.compile("^F +nom +(.*med) +D +0", re.M)
with open(export[0], "r") as fexp:
    content = fexp.read()

meshes = re_mesh.findall(content)

nb_mesh = 0
nb_mesh_converted = 0
conversion_error = []
# loop on mesh
for mesh_file in meshes:
    MPI.ASTER_COMM_WORLD.Barrier()
    mesh_name = osp.basename(mesh_file)
    nb_mesh += 1
    if mesh_name in failure:
        print("SKIP:", mesh_name)
        continue

    mesh_vers = mc.MEDFileVersionOfFileStr(mesh_name)
    vers_num = int(mesh_vers.replace(".", ""))

    print("MESHNAME %d: %s (version: %s)" % (nb_mesh, mesh_name, mesh_vers), flush=True)

    # convert old med mesh < 3.0.0
    if vers_num < 300:
        mfd = mc.MEDFileData(mesh_name)
        mesh_name = mesh_name.split(".")[0] + "_tmp.med"
        mfd.write(mesh_name, 2)
        print(
            "Mesh converted: %s (version: %s)" % (mesh_name, mc.MEDFileVersionOfFileStr(mesh_name)),
            flush=True,
        )

    # read parallel mesh and partitioning
    pmesh = CA.ParallelMesh()
    try:
        pmesh.readMedFile(mesh_name)
    except Exception as exc:
        test.assertIsNone(exc, "partitioning failed (%s)" % mesh_name)
        print("ERROR:", str(exc))
        conversion_error.append(mesh_name)
        continue

    # read std mesh
    mesh = CA.Mesh()
    mesh.readMedFile(mesh_name)

    if not pmesh.checkConsistency(mesh_name):
        print("ERROR in consistency check:", mesh_name)
        conversion_error.append(mesh_name)

    test.assertTrue(pmesh.getNumberOfNodes() > 0)
    test.assertTrue(pmesh.getNumberOfCells() > 0)

    nb_mesh_converted += 1

    MPI.ASTER_COMM_WORLD.Barrier()


list_nb_mesh = {2: 516, 3: 486, 4: 484}
list_nb_mesh_conv = {2: 515, 3: 481, 4: 483}
list_nb_conv_error = {2: 0, 3: 0, 4: 0}

print("Number of mesh: %s" % (nb_mesh), flush=True)
print("Number of mesh converted: %s" % (nb_mesh_converted), flush=True)
print("Number of conversion error: %s " % (len(conversion_error)), flush=True)
print("List files: ", conversion_error)

# all the mesh are partitioned
test.assertTrue(is_parallel)
test.assertEqual(nb_mesh, list_nb_mesh[nbproc])
test.assertEqual(nb_mesh_converted, list_nb_mesh_conv[nbproc])
test.assertEqual(len(conversion_error), list_nb_conv_error[nbproc])

test.printSummary()

CA.close()
