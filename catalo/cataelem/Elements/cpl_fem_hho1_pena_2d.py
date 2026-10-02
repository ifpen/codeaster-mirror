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


# ELEMENTARY TREATMENT OF 2D COUPLING FEM-HHO
# PENALISATION FORMULATION

# Convention:
# EN1: Slave nodes - FEM
# EN2: Master nodes - HHO

# ----------------
# Modes locaux :
# ----------------


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        # Slave nodes - FEM
        ("EN1", ("DX", "DY")),
        # Master nodes - HHO
        ("EN2", ()),
        ("EN3", ("HHO_FX[2]", "HHO_FY[2]")),
    ),
)


CHHOBS = LocatedComponents(
    phys=PHY.N3600R,
    type="ELNO",
    diff=True,
    components=(("EN1", ()), ("EN2", ()), ("EN3", ("X[3]",))),
)


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class CP_S2S3_HHO1(Element):
    """
    THE CP_S2S3_HHO1 CLASS ELEMENT : TRIA3/SEG2
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 2D : elementary treatments
    """

    meshType = MT.SEG23
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4)), SetOfNodes("EN3", (5,)))
    calculs = (
        OP.CHAR_MECA_CPL(
            te=520,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PPAIRR, LC.CPAIRR),
                (OP.CHAR_MECA_CPL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.HHO_CALC_PROJ(
            te=527,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (OP.HHO_CALC_PROJ.PCHHOBS, CHHOBS)),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.RIGI_CPL(
            te=520,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PPAIRR, LC.CPAIRR),
                (OP.RIGI_CPL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_ELAS_CPL(
            te=520,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PPAIRR, LC.CPAIRR),
                (OP.RIGI_ELAS_CPL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
                (SP.PNEUT_F, LC.CNTINIF),
                (SP.PNEUT_R, LC.CNTINIR),
                (SP.PNEUT_I, LC.CNTINII),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM2D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENTINIF),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENTINIR),
                (OP.TOU_INI_ELNO.PNEUT_I, LC.ENTINII),
            ),
        ),
    )


# ------------------------------------------------------------
class CP_S3S3_HHO1(CP_S2S3_HHO1):
    """
    CP_S3S3_HHO1 DERIVED FROM THE CP_S2S3_HHO1 CLASS ELEMENT : SEG3/SEG2
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 2D : elementary treatments
    """

    meshType = MT.SEG33
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5)), SetOfNodes("EN3", (6,)))


# ------------------------------------------------------------
class CP_S4S3_HHO1(CP_S2S3_HHO1):
    """
    CP_S3S3_HHO2 DERIVED FROM THE CP_S2S3_HHO2 CLASS ELEMENT : SEG3/SEG2
    LIAISON_ELEM / PENALISATION / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 2D : elementary treatments
    """

    meshType = MT.SEG43
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6)), SetOfNodes("EN3", (7,)))
