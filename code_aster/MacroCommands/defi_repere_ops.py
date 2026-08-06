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
from ..Cata.Syntax import _F
from ..CodeCommands import CREA_CHAMP

from math import cos, sin, atan2, pi
from typing import Sequence
import numpy as np


def can_perform_atan2(a: float, b: float, tol: float):
    arr = np.array([a, b])
    return np.any(np.abs(arr) > tol)


def matrot_eul(angles: Sequence[float]):
    "Rotation matrix from EULER angles"
    cos0 = cos(angles[0])
    cos1 = cos(angles[1])
    cos2 = cos(angles[2])
    sin0 = sin(angles[0])
    sin1 = sin(angles[1])
    sin2 = sin(angles[2])

    # Rotation matrix
    mat = np.zeros((3, 3))
    mat[0, 0] = cos2 * cos0 - sin2 * cos1 * sin0
    mat[1, 0] = cos2 * sin0 + sin2 * cos1 * cos0
    mat[2, 0] = sin2 * sin1
    mat[0, 1] = -sin2 * cos0 - cos2 * cos1 * sin0
    mat[1, 1] = -sin2 * sin0 + cos2 * cos1 * cos0
    mat[2, 1] = cos2 * sin1
    mat[0, 2] = sin1 * sin0
    mat[1, 2] = -sin1 * cos0
    mat[2, 2] = cos1
    return mat


def matrot_nau(angles: Sequence[float]):
    "Rotation matrix from NAUTIQUE angles"
    cos0 = cos(angles[0])
    cos1 = cos(angles[1])
    cos2 = cos(angles[2])
    sin0 = sin(angles[0])
    sin1 = sin(angles[1])
    sin2 = sin(angles[2])

    mat = np.zeros((3, 3))
    mat[0, 0] = cos1 * cos0
    mat[1, 0] = sin2 * sin1 * cos0 - cos2 * sin0
    mat[2, 0] = sin2 * sin0 + cos2 * sin1 * cos0
    mat[0, 1] = cos1 * sin0
    mat[1, 1] = cos2 * cos0 + sin2 * sin1 * sin0
    mat[2, 1] = cos2 * sin1 * sin0 - cos0 * sin2
    mat[0, 2] = -sin1
    mat[1, 2] = sin2 * cos1
    mat[2, 2] = cos2 * cos1
    return mat


def eulnau(angeul: Sequence[float]):
    "Transform EULER angles to NAUTIQUE angles"
    # Compute rotation matrix
    mat = matrot_eul(np.deg2rad(angeul))

    # Local vectors written in the local system
    xg = mat[:, 0]
    yg = mat[:, 1]
    zg = mat[:, 2]

    # Compute new angles
    angnau = angvxy(xg, yg)
    angnau = np.where(angnau < 0.0, angnau + 2 * pi, angnau)
    return np.rad2deg(angnau)


def angvxy(gx: Sequence[float], gn: Sequence[float], tol: float = 1e-14):
    "Compute 3 NAUTIQUE angles from vectors GX and GN"
    alpha, beta = angvx(gx, tol=tol)
    angle = np.zeros(3)
    angle[0] = alpha
    angle[1] = beta

    mat = matrot_nau(angle)
    gy = mat @ gn
    if can_perform_atan2(gy[2], gy[1], tol):
        angle[2] = atan2(gy[2], gy[1])
    return angle


def angvx(gx: Sequence[float], tol: float = 1e-14):
    "Compute 2 NAUTIQUE angles from vector GX"
    alpha = 0.0
    if can_perform_atan2(gx[1], gx[0], tol):
        alpha = atan2(gx[1], gx[0])

    p = np.sqrt(gx[0] ** 2 + gx[1] ** 2)
    beta = 0.0
    if can_perform_atan2(gx[2], p, tol):
        beta = -atan2(gx[2], p)
    return alpha, beta


def defi_repere_ops(self, **kwargs):

    # Recover model
    model = kwargs["MODELE"]

    # TODO: add verification with geometrical and topological dimension
    ndim = model.getGeometricDimension()

    # Define arguments for keyword AFFE in CREA_CHAMP
    AFFE_ARGS = {}

    # 1. Set TOUT or GROUP_MA
    DEFI_ARGS = kwargs["DEFI"]
    if "TOUT" in DEFI_ARGS:
        AFFE_ARGS["TOUT"] = DEFI_ARGS["TOUT"]
    elif "GROUP_MA" in DEFI_ARGS:
        AFFE_ARGS["GROUP_MA"] = DEFI_ARGS["GROUP_MA"]
    else:
        raise NotImplementedError

    # 2. Set NOM_CMP
    # TODO: replace X1 ... X7 for ("C", "ALPHA", "BETA", "KAPPA", "X", "Y", "Z")
    AFFE_ARGS["NOM_CMP"] = ("X1", "X2", "X3", "X4", "X5", "X6", "X7")

    # 3. Set VALE for different cases
    SYST_COOR = DEFI_ARGS["SYST_COOR"]
    if SYST_COOR == "CYLINDRIQUE":
        dtype = -1
        orig = DEFI_ARGS["ORIG_AXE"]
        angl = DEFI_ARGS["ANGL_AXE"]
    elif SYST_COOR == "CARTESIEN":
        dtype = 1
        orig = (0.0, 0.0, 0.0)
        if "ANGL_NAUT" in DEFI_ARGS:
            angl = DEFI_ARGS["ANGL_NAUT"]
        elif "ANGL_EULER" in DEFI_ARGS:
            eule = DEFI_ARGS["ANGL_EULER"]
            angl = eulnau([eule[0], eule[1], eule[2]])
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError

    VALE = (dtype, angl[0], angl[1], angl[2], orig[0], orig[1], orig[2])
    AFFE_ARGS["VALE"] = VALE

    # Create CREA_CHAMP datastructure
    champ = CREA_CHAMP(
        OPERATION="AFFE",
        TYPE_CHAM="ELEM_CAMA_R",  # or ELEM_NEUT_R
        MODELE=model,
        PROL_ZERO="OUI",
        AFFE=_F(**AFFE_ARGS),
    )

    return champ
