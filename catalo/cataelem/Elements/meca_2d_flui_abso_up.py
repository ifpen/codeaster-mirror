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

# ----------------
# Modes locaux :
# ----------------


DDL_MECA = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("PRES"))

NDEPLAC = LocatedComponents(phys=PHY.DEPL_C, type="ELNO", components=("PRES",))


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y"))


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "W")
)

EAMORFL = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X[2]"))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)


# ------------------------------------------------------------
class MEFASE2UP(Element):
    """Please document this element"""

    meshType = MT.SEG2
    elrefe = (ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2", "FPG1=FPG1"), mater=("RIGI", "FPG1")),)
    calculs = (
        OP.COOR_ELGA(
            te=478, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.FORC_NODA(
            te=167,
            para_in=((SP.PDEPLAR, DDL_MECA), (SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=167,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PCODRET, LC.ECODRET), (SP.PMATUUR, MMATUUR), (SP.PVECTUR, MVECTUR)),
        ),
        OP.AMOR_MECA(
            te=258,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.AMOR_MECA.PAMORFL, EAMORFL),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.CHAR_MECA_VFAC(
            te=255,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC), (SP.PVITEFR, LC.EVITEFR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_VFAC_F(
            te=255,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVITEFF, LC.EVITEFF),
                (SP.PINSTR, LC.MTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.MASS_MECA(
            te=184,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RAPH_MECA(
            te=167,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PCODRET, LC.ECODRET), (SP.PVECTUR, MVECTUR)),
        ),
        OP.RIGI_MECA(
            te=167,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_HYST(
            te=167,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PMATUUC, MMATUUC),),
        ),
        OP.RIGI_MECA_TANG(
            te=167,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=(
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.TOU_INI_ELGA(te=99, para_out=((OP.TOU_INI_ELGA.PGEOM_R, EGGEOP_R),)),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D), (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER)),
        ),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
    )


# ------------------------------------------------------------
class MEFASE3UP(MEFASE2UP):
    """Please document this element"""

    meshType = MT.SEG3
    elrefe = (ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("RIGI", "FPG1")),)
