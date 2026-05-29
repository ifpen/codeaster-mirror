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

# ELEMENTARY TREATMENT OF 3D FRICTIONLESS ELEMENT WITH DEFI_CONTACT OPERATOR

# ----------------
# Modes locaux :
# ----------------


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        # Slave nodes with LAG_C
        ("EN1", ("DX", "DY", "DZ", "LAGS_C", "LAGS_F1", "LAGS_F2")),
        # Slave nodes without LAG_C
        ("EN2", ("DX", "DY", "DZ")),
        # Master nodes
        ("EN3", ("DX", "DY", "DZ")),
    ),
)


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))

CGAPR = LocatedComponents(
    phys=PHY.NEUT_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("X1",)), ("EN2", ("X1",)), ("EN3", ())),
)

CSTATR = LocatedComponents(
    phys=PHY.NEUT_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("X1",)), ("EN2", ("X1",)), ("EN3", ())),
)

ECCONT = LocatedComponents(
    phys=PHY.CONT_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("COEF_C",)), ("EN2", ()), ("EN3", ())),
)

ECFROT = LocatedComponents(
    phys=PHY.CONT_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("COEF_F",)), ("EN2", ()), ("EN3", ())),
)

CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))

MVECGAP = ArrayOfComponents(phys=PHY.VNEU_R, locatedComponents=CGAPR)

MVEIGAP = ArrayOfComponents(phys=PHY.VNEU_R, locatedComponents=CSTATR)

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class FMQ4Q4(Element):
    """
    THE FMQ4Q4 CLASS ELEMENT : SEG2/SEG2 (3D Edge / 3D edge )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD4 SLAVE  ELEMENT : 1-2-3-4 (DX,DY,DZ,LAGS_C)
        QUAD4 MASTER ELEMENT : 5-6-7-8 (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QUAD44
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", ()),
        SetOfNodes("EN3", (5, 6, 7, 8)),
    )
    calculs = (
        OP.CHAR_MECA_CONT(
            te=355,
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
            te=355,
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
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D), (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER)),
        ),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
    )


# ------------------------------------------------------------
class FMT3T3(FMQ4Q4):
    """
    THE COT3T3 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : TRIA3/TRIA3 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        TRIA3 MASTER  ELEMENT : 4-5-6 (DX,DY,DZ)
        TRIA3 SLAVE ELEMENT : 1-2-3 (DX,DY,DZ,LAGS_C)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.TRIA33
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", ()), SetOfNodes("EN3", (4, 5, 6)))


# ------------------------------------------------------------
class FMQ4T3(FMQ4Q4):
    """
    THE COQ4T3 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD4/TRIA3 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD4 SLAVE  ELEMENT : 1-2-3-4 (DX,DY,DZ,LAGS_C)
        TRIA3 MASTER ELEMENT : 5-6-7   (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU4TR3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", ()), SetOfNodes("EN3", (5, 6, 7)))


# ------------------------------------------------------------
class FMT3Q4(FMQ4Q4):
    """
    THE COT3Q4 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : TRIA3/QUAD4 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        TRIA3 SLAVE  ELEMENT : 1-2-3     (DX,DY,DZ,LAGS_C)
        QUAD4 MASTER ELEMENT : 4-5-6-7   (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.TR3QU4
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", ()), SetOfNodes("EN3", (4, 5, 6, 7)))


# ------------------------------------------------------------
class FMT6T3(FMQ4Q4):
    """
    THE COT6T3 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : TRIA6/TRIA3 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        TRIA6 SLAVE  ELEMENT : 1-2-3-4-5-6     (DX,DY,DZ,LAGS_C)
        TRIA3 MASTER ELEMENT : 7-8-9           (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.TR6TR3
    nodes = (
        SetOfNodes("EN1", (1, 2, 3)),
        SetOfNodes("EN2", (4, 5, 6)),
        SetOfNodes("EN3", (7, 8, 9)),
    )


# ------------------------------------------------------------
class FMT3T6(FMQ4Q4):
    """
    THE COT3T6 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : TRIA3/TRIA6 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        TRIA3 SLAVE  ELEMENT : 1-2-3          (DX,DY,DZ,LAGS_C)
        TRIA6 MASTER ELEMENT : 4-5-6-7-8-9    (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.TR3TR6
    nodes = (
        SetOfNodes("EN1", (1, 2, 3)),
        SetOfNodes("EN2", ()),
        SetOfNodes("EN3", (4, 5, 6, 7, 8, 9)),
    )


# ------------------------------------------------------------
class FMT6Q4(FMQ4Q4):
    """
    THE COT6Q4 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : TRIA6/QUAD4 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        TRIA6 SLAVE  ELEMENT : 1-2-3-4-5-6       (DX,DY,DZ,LAGS_C)
        QUAD4 MASTER ELEMENT : 7-8-9-10          (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.TR6QU4
    nodes = (
        SetOfNodes("EN1", (1, 2, 3)),
        SetOfNodes("EN2", (4, 5, 6)),
        SetOfNodes("EN3", (7, 8, 9, 10)),
    )


# ------------------------------------------------------------
class FMQ4T6(FMQ4Q4):
    """
    THE COQ4T6 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD4/TRIA6 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD4 SLAVE  ELEMENT : 1-2-3-4       (DX,DY,DZ,LAGS_C)
        TRIA6 MASTER ELEMENT : 5-6-7-8-9-10  (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU4TR6
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", ()),
        SetOfNodes("EN3", (5, 6, 7, 8, 9, 10)),
    )


# ------------------------------------------------------------
class FMT6Q8(FMQ4Q4):
    """
    THE COT6Q8 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : TRIA6/QUAD8 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        TRIA6 SLAVE  ELEMENT : 1-2-3-4-5-6           (DX,DY,DZ,LAGS_C)
        QUAD8 MASTER ELEMENT : 7-8-9-10-11-12-13-14  (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.TR6QU8
    nodes = (
        SetOfNodes("EN1", (1, 2, 3)),
        SetOfNodes("EN2", (4, 5, 6)),
        SetOfNodes("EN3", (7, 8, 9, 10, 11, 12, 13, 14)),
    )


# ------------------------------------------------------------
class FMQ8T6(FMQ4Q4):
    """
    THE COQ8T6 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD8/TRIA6 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD8 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8   (DX,DY,DZ,LAGS_C)
        TRIA6 MASTER ELEMENT : 9-10-11-12-13-14  (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU8TR6
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8)),
        SetOfNodes("EN3", (9, 10, 11, 12, 13, 14)),
    )


# ------------------------------------------------------------
class FMT6Q9(FMQ4Q4):
    """
    THE COT6Q9 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : TRIA6/QUAD9 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        TRIA6 SLAVE  ELEMENT : 1-2-3-4-5-6              (DX,DY,DZ,LAGS_C)
        QUAD9 MASTER ELEMENT : 7-8-9-10-11-12-13-14-15  (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.TR6QU9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3)),
        SetOfNodes("EN2", (4, 5, 6)),
        SetOfNodes("EN3", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
    )


# ------------------------------------------------------------
class FMQ9T6(FMQ4Q4):
    """
    THE COQ9T6 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD9/TRIA6 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9  (DX,DY,DZ,LAGS_C)
        TRIA6 MASTER ELEMENT : 10-11-12-13-14-15  (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU9TR6
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8, 9)),
        SetOfNodes("EN3", (10, 11, 12, 13, 14, 15)),
    )


# ------------------------------------------------------------
class FMQ8T3(FMQ4Q4):
    """
    THE COQ8T3 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD8/TRIA3 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD8 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8  (DX,DY,DZ,LAGS_C)
        TRIA3 MASTER ELEMENT : 9-10-11          (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU8TR3
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8)),
        SetOfNodes("EN3", (9, 10, 11)),
    )


# ------------------------------------------------------------
class FMT3Q8(FMQ4Q4):
    """
    THE COT3Q8 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : TRIA3/QUAD8 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        TRIA3 SLAVE  ELEMENT : 1-2-3                (DX,DY,DZ,LAGS_C)
        QUAD8 MASTER ELEMENT : 4-5-6-7-8-9-10-11    (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.TR3QU8
    nodes = (
        SetOfNodes("EN1", (1, 2, 3)),
        SetOfNodes("EN2", ()),
        SetOfNodes("EN3", (4, 5, 6, 7, 8, 9, 10, 11)),
    )


# ------------------------------------------------------------
class FMQ8Q4(FMQ4Q4):
    """
    THE COQ8Q4 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD8/QUAD4 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD8 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8   (DX,DY,DZ,LAGS_C)
        QUAD4 MASTER ELEMENT : 9-10-11-12        (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU8QU4
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8)),
        SetOfNodes("EN3", (9, 10, 11, 12)),
    )


# ------------------------------------------------------------
class FMQ4Q8(FMQ4Q4):
    """
    THE COQ4Q8 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD4/QUAD8 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD4 SLAVE  ELEMENT : 1-2-3-4              (DX,DY,DZ,LAGS_C)
        QUAD8 MASTER ELEMENT : 5-6-7-8-9-10-11-12   (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU4QU8
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", ()),
        SetOfNodes("EN3", (5, 6, 7, 8, 9, 10, 11, 12)),
    )


# ------------------------------------------------------------
class FMQ8Q9(FMQ4Q4):
    """
    THE COQ8Q9 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD8/QUAD9 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD8 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8              (DX,DY,DZ,LAGS_C)
        QUAD9 MASTER ELEMENT : 9-10-11-12-13-14-15-16-17    (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU8QU9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8)),
        SetOfNodes("EN3", (9, 10, 11, 12, 13, 14, 15, 16, 17)),
    )


# ------------------------------------------------------------
class FMQ9Q8(FMQ4Q4):
    """
    THE COQ9Q8 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD9/QUAD8 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9              (DX,DY,DZ,LAGS_C)
        QUAD8 MASTER ELEMENT : 10-11-12-13-14-15-16-17        (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU9QU8
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8, 9)),
        SetOfNodes("EN3", (10, 11, 12, 13, 14, 15, 16, 17)),
    )


# ------------------------------------------------------------
class FMQ9Q4(FMQ4Q4):
    """
    THE COQ9Q4 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD9/QUAD4 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9      (DX,DY,DZ,LAGS_C)
        QUAD4 MASTER ELEMENT : 10-11-12-13            (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU9QU4
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8, 9)),
        SetOfNodes("EN3", (10, 11, 12, 13)),
    )


# ------------------------------------------------------------
class FMQ4Q9(FMQ4Q4):
    """
    THE COQ4Q9 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD4/QUAD9 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD4 SLAVE  ELEMENT : 1-2-3-4                    (DX,DY,DZ,LAGS_C)
        QUAD9 MASTER ELEMENT : 5-6-7-8-9-10-11-12-13      (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU4QU9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", ()),
        SetOfNodes("EN3", (5, 6, 7, 8, 9, 10, 11, 12, 13)),
    )


# ------------------------------------------------------------
class FMQ9T3(FMQ4Q4):
    """
    THE COQ9T3 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD9/TRIA3 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9   (DX,DY,DZ,LAGS_C)
        TRIA3 MASTER ELEMENT : 10-11-12         (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QU9TR3
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8, 9)),
        SetOfNodes("EN3", (10, 11, 12)),
    )


# ------------------------------------------------------------
class FMT3Q9(FMQ4Q4):
    """
    THE COT3Q9 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : TRIA3/QUAD9 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        TRIA3 SLAVE  ELEMENT : 1-2-3                     (DX,DY,DZ,LAGS_C)
        QUAD9 MASTER ELEMENT : 4-5-6-7-8-9-10-11-12      (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.TR3QU9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3)),
        SetOfNodes("EN2", ()),
        SetOfNodes("EN3", (4, 5, 6, 7, 8, 9, 10, 11, 12)),
    )


# ------------------------------------------------------------
class FMQ8Q8(FMQ4Q4):
    """
    THE COQ8Q8 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : TRIA3/QUAD9 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD8 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8                (DX,DY,DZ,LAGS_C)
        QUAD8 MASTER ELEMENT : 9-10-11-12-13-14-14-15-16      (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QUAD88
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8)),
        SetOfNodes("EN3", (9, 10, 11, 12, 13, 14, 15, 16)),
    )


# ------------------------------------------------------------
class FMQ9Q9(FMQ4Q4):
    """
    THE COQ9Q9 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD9/QUAD9 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9                  (DX,DY,DZ,LAGS_C)
        QUAD9 MASTER ELEMENT : 10-11-12-13-14-14-15-16-17-18      (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.QUAD99
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8, 9)),
        SetOfNodes("EN3", (10, 11, 12, 13, 14, 15, 16, 17, 18)),
    )


# ------------------------------------------------------------
class FMT6T6(FMQ4Q4):
    """
    THE COQ9Q9 DERIVED FROM  FMQ4Q4 CLASS ELEMENT  : QUAD9/QUAD9 (3D Face / 3D Face )
    DEFI_CONTACT / LAGRANGIAN / SURFACE-TO-SURFACE
        Slave frictionless Contact Element in 3D  : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9                  (DX,DY,DZ,LAGS_C)
        QUAD9 MASTER ELEMENT : 10-11-12-13-14-14-15-16-17-18      (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=364)
        PMMATUR : SYMMETRIC MATRIX (te=364)
        PMMATUR : VECTOR OF CONTACT LOAD (te=365)
    """

    meshType = MT.TRIA66
    nodes = (
        SetOfNodes("EN1", (1, 2, 3)),
        SetOfNodes("EN2", (4, 5, 6)),
        SetOfNodes("EN3", (7, 8, 9, 10, 11, 12)),
    )


# ------------------------------------------------------------
class FMP1L3(FMQ4Q4):
    """
    FMS3S2 DERIVED FROM THE FMQ4Q4 CLASS ELEMENT : POINT1.
    This element is for nodes that have Lagrange and are not pairing
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        POI1 SLAVE  ELEMENT : 1   (DX,DY,DZ,LAGS_C)
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
        PMMATUR : VECTOR OF CONTACT LOAD (te=355)
    """

    meshType = MT.POI1
    nodes = (SetOfNodes("EN1", (1,)), SetOfNodes("EN2", ()), SetOfNodes("EN3", ()))


# ------------------------------------------------------------


class FMP1N3(FMQ4Q4):
    """
    FMS3S2 DERIVED FROM THE FMQ4Q4 CLASS ELEMENT : POINT1.
    This element is for nodes that have no Lagrange and are not pairing
    DEFI_CONTACT / MORTAR / SEGMENT-TO-SEGMENT
        Slave frictionless Contact Element in 2D : elementary treatments
    Local Numerotation :
        POI1 SLAVE  ELEMENT : 1   (DX,DY,DZ)
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
        PMMATUR : VECTOR OF CONTACT LOAD (te=355)
    """

    meshType = MT.POI1
    nodes = (SetOfNodes("EN1", ()), SetOfNodes("EN2", (1,)), SetOfNodes("EN3", ()))
