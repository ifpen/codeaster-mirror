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
    phys=PHY.CONT_R, type="ELNO", diff=True, components=(("EN1", ("COEF_C")), ("EN2", ()))
)

ECFROT = LocatedComponents(
    phys=PHY.CONT_R, type="ELNO", diff=True, components=(("EN1", ("COEF_F",)), ("EN2", ()))
)

NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))

CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class CPS2S2(Element):
    """
    THE CPS2S2 CLASS ELEMENT : SEG2/SEG2
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numbering :
        SEG2 SLAVE  ELEMENT : 1-2 (DX,DY)
        SEG2 MASTER ELEMENT : 3-4 (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.SEG22
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4)))
    calculs = (
        OP.CHAR_MECA_CONT(
            te=352,
            para_in=(
                (SP.PCONFR, LC.CCONFR),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PGEOMER, NGEOMER),
                (SP.PGEOMCR, NGEOMER),
                (SP.PCCONTR, ECCONT),
                (SP.PCFROTR, ECFROT),
            ),
            para_out=((SP.PVECTCR, MVECTUR), (SP.PVECTFR, MVECTUR)),
        ),
        OP.RIGI_CONT(
            te=352,
            para_in=(
                (SP.PCONFR, LC.CCONFR),
                (SP.PDEPL_M, DDL_MECA),
                (SP.PDEPL_P, DDL_MECA),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PGEOMER, NGEOMER),
                (SP.PGEOMCR, NGEOMER),
                (SP.PCCONTR, ECCONT),
                (SP.PCFROTR, ECFROT),
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
class CPS3S3(CPS2S2):
    """
    CPS2S2A DERIVED FROM THE CPS2S2 CLASS ELEMENT : SEG2/SEG2 (AXIS)
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2  (DX,DY)
        SEG2 MASTER ELEMENT : 3-4  (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.SEG33
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6)))


# ------------------------------------------------------------
class CPS2S3(CPS2S2):
    """
    CPS2S2A DERIVED FROM THE CPS2S2 CLASS ELEMENT : SEG2/SEG2 (AXIS)
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2  (DX,DY)
        SEG2 MASTER ELEMENT : 3-4  (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.SEG23
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4, 5)))


# ------------------------------------------------------------
class CPS3S2(CPS2S2):
    """
    CPS2S2A DERIVED FROM THE CPS2S2 CLASS ELEMENT : SEG2/SEG2 (AXIS)
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2  (DX,DY)
        SEG2 MASTER ELEMENT : 3-4  (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.SEG32
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5)))


# ------------------------------------------------------------
class CPS2S2A(CPS2S2):
    """
    CPS2S2A DERIVED FROM THE CPS2S2 CLASS ELEMENT : SEG2/SEG2 (AXIS)
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2  (DX,DY)
        SEG2 MASTER ELEMENT : 3-4  (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.SEG22
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4)))


# ------------------------------------------------------------
class CPS3S3A(CPS2S2):
    """
    CPS3S3 DERIVED FROM THE CPS2S2 CLASS ELEMENT  : SEG3/SEG3 (AXIS)
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 4-5-6 (DX,DY)
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
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.SEG33
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6)))


# ------------------------------------------------------------
class CPS2S3A(CPS2S2):
    """
    CPS2S3 DERIVED FROM THE CPS2S2 CLASS ELEMENT : SEG2/SEG3 (AXIS)
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG2 SLAVE  ELEMENT : 1-2   (DX,DY)
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
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.SEG23
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4, 5)))


# ------------------------------------------------------------
class CPS3S2A(CPS2S2):
    """
    CPS3S2 DERIVED FROM THE CPS2S2 CLASS ELEMENT : SEG3/SEG2 (AXIS)
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        SEG3 SLAVE  ELEMENT : 1-2   (DX,DY) / 3 - (DX,DY)
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
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.SEG32
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5)))


# ------------------------------------------------------------


class CPP1L2(CPS2S2):
    """
    CPS3S2 DERIVED FROM THE CPS2S2 CLASS ELEMENT : POINT1.
    This element is for nodes that have Lagrange and are not pairing
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        POI1 SLAVE  ELEMENT : 1   (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.POI1
    nodes = (SetOfNodes("EN1", (1,)), SetOfNodes("EN2", ()))


# ------------------------------------------------------------


class CPP1N2(CPS2S2):
    """
    CPS3S2 DERIVED FROM THE CPS2S2 CLASS ELEMENT : POINT1.
    This element is for nodes that have no Lagrange and are not pairing
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        POI1 SLAVE  ELEMENT : 1   (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.POI1
    nodes = (SetOfNodes("EN1", ()), SetOfNodes("EN2", (1,)))


# ------------------------------------------------------------
class CPP1L2A(CPS2S2):
    """
    CPS3S2 DERIVED FROM THE CPS2S2 CLASS ELEMENT : POINT1.
    This element is for nodes that have Lagrange and are not pairing
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        POI1 SLAVE  ELEMENT : 1   (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.POI1
    nodes = (SetOfNodes("EN1", (1,)), SetOfNodes("EN2", ()))


# ------------------------------------------------------------


class CPP1N2A(CPS2S2):
    """
    CPS3S2 DERIVED FROM THE CPS2S2 CLASS ELEMENT : POINT1.
    This element is for nodes that have no Lagrange and are not pairing
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave Penalized Frictional Contact Element in 2D : elementary treatments
    Local Numerotation :
        POI1 SLAVE  ELEMENT : 1   (DX,DY)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.POI1
    nodes = (SetOfNodes("EN1", ()), SetOfNodes("EN2", (1,)))
