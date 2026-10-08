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


# ELEMENTARY TREATMENT OF 3D COUPLING FEM-HHO
# NITSCHE FORMULATION

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
        ("EN1", ("DX", "DY", "DZ")),
        # Master nodes - HHO
        ("EN2", ()),
        ("EN3", ("HHO_FX[3]", "HHO_FY[3]", "HHO_FZ[3]")),
    ),
)


CHHOBS = LocatedComponents(
    phys=PHY.N3600R,
    type="ELNO",
    diff=True,
    components=(("EN1", ()), ("EN2", ()), ("EN3", ("X[6]",))),
)


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class CN_TE4TR7_HHO1(Element):
    """
    THE CN_TE4TR7_HHO1 CLASS ELEMENT : TETRA4/TRIA7
    LIAISON_MASSIF / NITSCHE / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.TE4TR7
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8, 9, 10)),
        SetOfNodes("EN3", (11,)),
    )
    calculs = (
        OP.HHO_CALC_PROJ(
            te=527,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPAIRR, LC.CPAIRR),
                (OP.HHO_CALC_PROJ.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.HHO_EXTR_PROJ(
            te=547,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PMAELNS1, MMATUNS), (SP.PPAIRR, LC.CPAIRR)),
            para_out=((OP.HHO_EXTR_PROJ.PMPHHOR, LC.CMPHHOR),),
        ),
        OP.RIGI_ELAS_CPL(
            te=519,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPAIRR, LC.CPAIRR),
                (OP.RIGI_ELAS_CPL.PCHHOBS, CHHOBS),
                (OP.RIGI_ELAS_CPL.PMATERR, LC.ENMATE_R),
                (OP.RIGI_ELAS_CPL.PMPHHOR, LC.CMPHHOR),
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
class CN_T10TR7_HHO1(CN_TE4TR7_HHO1):
    """
    CN_T10TR7_HHO1 DERIVED FROM THE CN_TE4TR7_HHO1 CLASS ELEMENT : SEG3/SEG2
    LIAISON_MASSIF / NITSCHE / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.T10TR7
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
        SetOfNodes("EN2", (11, 12, 13, 14, 15, 16)),
        SetOfNodes("EN3", (17,)),
    )


# ------------------------------------------------------------
class CN_T20TR7_HHO1(CN_TE4TR7_HHO1):
    """
    CN_T10TR7_HHO1 DERIVED FROM THE CN_TE4TR7_HHO1 CLASS ELEMENT : SEG3/SEG2
    LIAISON_MASSIF / NITSCHE / SEGMENT-TO-SEGMENT
        Coupling FEM/HHO Element in 3D : elementary treatments
    """

    meshType = MT.T20TR7
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN2", (21, 22, 23, 24, 25, 26)),
        SetOfNodes("EN3", (27,)),
    )
