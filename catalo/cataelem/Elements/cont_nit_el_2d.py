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


from cataelem.Tools.base_objects import LocatedComponents, ArrayOfComponents, SetOfNodes, ElrefeLoc
from cataelem.Tools.base_objects import Calcul, Element
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.located_components as LC
import cataelem.Commons.parameters as SP
import cataelem.Commons.mesh_types as MT
from cataelem.Options.options import OP


# ELEMENTARY TREATMENT OF 2D FRICTIONLESS ELEMENT WITH DEFI_CONTACT OPERATOR
# NITSCHE FORMULATION

# Convention:
# EN1: Slave nodes
# EN2: Master nodes

# ----------------
# Modes locaux :
# ----------------


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        # Slave nodes
        ("EN1", ("DX", "DY")),
        # Master nodes
        ("EN2", ("DX", "DY")),
    ),
)

ECCONT = LocatedComponents(
    phys=PHY.CONT_R, type="ELNO", diff=True, components=(("EN1", ("COEF_C",)), ("EN2", ()))
)

NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))

CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class CNT3S2(Element):
    """
    THE CNT3S2 CLASS ELEMENT : TRIA3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2-3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.TR3SE2
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5)))
    calculs = (
        OP.CHAR_MECA_CONT(
            te=356,
            para_in=(
                (SP.PCONFR, LC.CCONFR),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PGEOMER, NGEOMER),
                (SP.PGEOMCR, NGEOMER),
                (SP.PCCONTR, ECCONT),
            ),
            para_out=((SP.PVECTCR, MVECTUR),),
        ),
        OP.RIGI_CONT(
            te=356,
            para_in=(
                (SP.PCONFR, LC.CCONFR),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PGEOMER, NGEOMER),
                (SP.PGEOMCR, NGEOMER),
                (SP.PCCONTR, ECCONT),
            ),
            para_out=((SP.PMATUUR, MMATUUR), (SP.PMATUNS, MMATUNS)),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D), (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER)),
        ),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
    )


# ------------------------------------------------------------
class CNT3S3(CNT3S2):
    """
    CNT3S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT  : SEG3/SEG3
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 4-5-6 (DX,DY,LAGS_C)
        SEG3 MASTER ELEMENT : 1-2-3 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.TR3SE3
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6)))


# ------------------------------------------------------------
class CNT6S2(CNT3S2):
    """
    CNT6S2 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG2/SEG3
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C)
        SEG3 MASTER ELEMENT : 3-4-5 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.TR6SE2
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)), SetOfNodes("EN2", (7, 8)))


# ------------------------------------------------------------
class CNT6S3(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.TR6SE3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)), SetOfNodes("EN2", (7, 8, 9)))


# ------------------------------------------------------------


class CNQ4S2(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU4SE2
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6)))


# ------------------------------------------------------------
class CNQ4S3(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU4SE3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7)))


# ------------------------------------------------------------
class CNQ8S2(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU8SE2
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)), SetOfNodes("EN2", (9, 10)))


# ------------------------------------------------------------


class CNQ8S3(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU8SE3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)), SetOfNodes("EN2", (9, 10, 11)))


# ------------------------------------------------------------


class CNQ9S2(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU9SE2
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)), SetOfNodes("EN2", (10, 11)))


# ------------------------------------------------------------


class CNQ9S3(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU9SE3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)), SetOfNodes("EN2", (10, 11, 12)))


# ------------------------------------------------------------
class CNT3S2A(CNT3S2):
    """
    THE CNT3S2 CLASS ELEMENT : TRIA3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2-3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.TR3SE2
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5)))


# ------------------------------------------------------------
class CNT3S3A(CNT3S2):
    """
    CNT3S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT  : SEG3/SEG3
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 4-5-6 (DX,DY,LAGS_C)
        SEG3 MASTER ELEMENT : 1-2-3 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.TR3SE3
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6)))


# ------------------------------------------------------------
class CNT6S2A(CNT3S2):
    """
    CNT6S2 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG2/SEG3
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C)
        SEG3 MASTER ELEMENT : 3-4-5 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.TR6SE2
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)), SetOfNodes("EN2", (7, 8)))


# ------------------------------------------------------------
class CNT6S3A(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.TR6SE3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)), SetOfNodes("EN2", (7, 8, 9)))


# ------------------------------------------------------------


class CNQ4S2A(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU4SE2
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6)))


# ------------------------------------------------------------
class CNQ4S3A(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU4SE3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7)))


# ------------------------------------------------------------
class CNQ8S2A(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU8SE2
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)), SetOfNodes("EN2", (9, 10)))


# ------------------------------------------------------------


class CNQ8S3A(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU8SE3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)), SetOfNodes("EN2", (9, 10, 11)))


# ------------------------------------------------------------


class CNQ9S2A(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU9SE2
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)), SetOfNodes("EN2", (10, 11)))


# ------------------------------------------------------------


class CNQ9S3A(CNT3S2):
    """
    CNT6S3 DERIVED FROM THE CNT3S2 CLASS ELEMENT : SEG3/SEG2
    DEFI_CONTACT / NITSCHE / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY,LAGS_C) / 3 (DX,DY)
        SEG2 MASTER ELEMENT : 4-5     (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=356)
        PMMATUR : SYMMETRIC MATRIX (te=356)
        PMMATUR : VECTOR OF CONTACT LOAD (te=356)
    """

    meshType = MT.QU9SE3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)), SetOfNodes("EN2", (10, 11, 12)))
