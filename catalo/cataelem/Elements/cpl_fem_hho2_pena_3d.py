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

# ELEMENTARY TREATMENT OF 3D PENALIZATION FRICTIONLESS ELEMENT WITH DEFI_CONTACT OPERATOR

# ----------------
# Modes locaux :
# ----------------


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        # Slave nodes - FEM
        ("EN1", ("DX", "DY", "DZ")),
        # Master nodes - HHO
        ("EN2", ()),
        ("EN3", ("HHO_FX[6]", "HHO_FY[6]", "HHO_FZ[6]")),
    ),
)

CHHOBS = LocatedComponents(
    phys=PHY.N3600R,
    type="ELNO",
    diff=True,
    components=(("EN1", ()), ("EN2", ()), ("EN3", ("X[21]",))),
)


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class CP_Q4Q9_HHO2(Element):
    """
    THE CMQ4Q4 CLASS ELEMENT : SEG2/SEG2 (3D Edge / 3D edge )
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.QU4QU9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8, 9, 10, 11, 12)),
        SetOfNodes("EN3", (13,)),
    )
    calculs = (
        OP.CHAR_MECA_CPL(
            te=520,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPAIRR, LC.CPAIRR),
                (OP.CHAR_MECA_CPL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.HHO_CALC_PROJ(
            te=527,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPAIRR, LC.CPAIRR),
                (OP.HHO_CALC_PROJ.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.RIGI_CPL(
            te=520,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPAIRR, LC.CPAIRR),
                (OP.RIGI_CPL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_ELAS_CPL(
            te=520,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPAIRR, LC.CPAIRR),
                (OP.RIGI_ELAS_CPL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
                (SP.PNEUT_F, LC.CNTINIF),
                (SP.PNEUT_R, LC.CNTINIR),
                (SP.PNEUT_I, LC.CNTINII),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM3D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENTINIF),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENTINIR),
                (OP.TOU_INI_ELNO.PNEUT_I, LC.ENTINII),
            ),
        ),
    )


# ------------------------------------------------------------
class CP_Q8Q9_HHO2(CP_Q4Q9_HHO2):
    """
    THE COT3T3 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA3/TRIA3 (3D Face / 3D Face )
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.QU8QU9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16)),
        SetOfNodes("EN3", (17,)),
    )


# ------------------------------------------------------------
class CP_Q9Q9_HHO2(CP_Q4Q9_HHO2):
    """
    THE COQ4T3 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD4/TRIA3 (3D Face / 3D Face )
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.QUAD99
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)),
        SetOfNodes("EN2", (10, 11, 12, 13, 14, 15, 16, 17)),
        SetOfNodes("EN3", (18,)),
    )


# ------------------------------------------------------------
class CP_T3Q9_HHO2(CP_Q4Q9_HHO2):
    """
    THE COT3Q4 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA3/QUAD4 (3D Face / 3D Face )
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.TR3QU9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3)),
        SetOfNodes("EN2", (4, 5, 6, 7, 8, 9, 10, 11)),
        SetOfNodes("EN3", (12,)),
    )


# ------------------------------------------------------------
class CP_T6Q9_HHO2(CP_Q4Q9_HHO2):
    """
    THE COT6T3 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA6/TRIA3 (3D Face / 3D Face )
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.TR6QU9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12, 13, 14)),
        SetOfNodes("EN3", (15,)),
    )


# ------------------------------------------------------------
class CP_Q4T7_HHO2(CP_Q4Q9_HHO2):
    """
    THE COT3T6 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA3/TRIA6 (3D Face / 3D Face )
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.QU4TR7
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8, 9, 10)),
        SetOfNodes("EN3", (11,)),
    )


# ------------------------------------------------------------
class CP_Q8T7_HHO2(CP_Q4Q9_HHO2):
    """
    THE COT6Q4 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA6/QUAD4 (3D Face / 3D Face )
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.QU8TR7
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14)),
        SetOfNodes("EN3", (15,)),
    )


# ------------------------------------------------------------
class CP_Q9T7_HHO2(CP_Q4Q9_HHO2):
    """
    THE COQ4T6 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD4/TRIA6 (3D Face / 3D Face )
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.QU9TR7
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)),
        SetOfNodes("EN2", (10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN3", (16,)),
    )


# ------------------------------------------------------------
class CP_T3T7_HHO2(CP_Q4Q9_HHO2):
    """
    THE COT6Q8 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA6/QUAD8 (3D Face / 3D Face )
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.TR3TR7
    nodes = (
        SetOfNodes("EN1", (1, 2, 3)),
        SetOfNodes("EN2", (4, 5, 6, 7, 8, 9)),
        SetOfNodes("EN3", (10,)),
    )


# ------------------------------------------------------------
class CP_T6T9_HHO2(CP_Q4Q9_HHO2):
    """
    THE COQ8T6 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD8/TRIA6 (3D Face / 3D Face )
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.TR6TR7
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12)),
        SetOfNodes("EN3", (13,)),
    )
