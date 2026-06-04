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
import cataelem.Commons.attributes as AT

# --------------------------------------------------------------------------------------------------
# Located components
# --------------------------------------------------------------------------------------------------
DDL_THER = LocatedComponents(
    phys=PHY.TEMP_R, type="ELNO", diff=True, components=(("EN1", ("HHO_FT[6]",)), ("EN2", ()))
)

CHHOBS = LocatedComponents(
    phys=PHY.N3600R, type="ELNO", diff=True, components=(("EN1", ("X[21]",)), ("EN2", ()))
)

MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)


# --------------------------------------------------------------------------------------------------
class THER3DQU9_HHO2_F(Element):
    """Thermics - Skin element 3D_HHO_QUAD - QUAD"""

    meshType = MT.QUAD9
    nodes = (SetOfNodes("EN1", (9,)), SetOfNodes("EN2", (1, 2, 3, 4, 5, 6, 7, 8)))
    attrs = ((AT.BORD_ISO, "OUI"),)
    elrefe = (ElrefeLoc(MT.QU9, gauss=("RIGI=FPG9",), mater=("RIGI",)),)
    calculs = (
        OP.CHAR_THER_FLUN_F(
            te=461,
            para_in=(
                (SP.PFLUXNF, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_FLUN_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_R(
            te=461,
            para_in=(
                (SP.PFLUXNR, LC.CFLUXNR),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_FLUN_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUNL(
            te=461,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_FLUNL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ECHA_F(
            te=461,
            para_in=(
                (SP.PCOEFHF, LC.CHECHPF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (SP.PT_EXTF, LC.CTEMPEF),
                (OP.CHAR_THER_ECHA_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ECHA_R(
            te=461,
            para_in=(
                (SP.PCOEFHR, LC.CHECHPR),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (SP.PT_EXTR, LC.ET_EXTR),
                (OP.CHAR_THER_ECHA_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_RAYO_F(
            te=461,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PRAYONF, LC.CRAYONF),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_RAYO_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_RAYO_R(
            te=461,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PRAYONR, LC.CRAYONR),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_RAYO_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.RIGI_THER_ECHA_F(
            te=457,
            para_in=(
                (SP.PCOEFHF, LC.CHECHPF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTR, LC.CTIMETR),
                (OP.RIGI_THER_ECHA_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_ECHA_R(
            te=457,
            para_in=(
                (SP.PCOEFHR, LC.CHECHPR),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTR, LC.CTIMETR),
                (OP.RIGI_THER_ECHA_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_FLUXNL(
            te=457,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MTAN_THER_FLUXNL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_RAYO_F(
            te=457,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PRAYONF, LC.CRAYONF),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MTAN_THER_RAYO_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_RAYO_R(
            te=457,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PRAYONR, LC.CRAYONR),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MTAN_THER_RAYO_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PGEOM_R, LC.EGGEO3D),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
            ),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),
                (OP.TOU_INI_ELEM.PNEU1_R, LC.CNEUTR1),
                (OP.TOU_INI_ELEM.PCOEH_R, LC.CHECHPR),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM3D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
            ),
        ),
    )


# --------------------------------------------------------------------------------------------------
class THER3DTR7_HHO2_F(THER3DQU9_HHO2_F):
    """Thermics - Skin element 3D_HHO_QUAD - TRIA"""

    meshType = MT.TRIA7
    nodes = (SetOfNodes("EN1", (7,)), SetOfNodes("EN2", (1, 2, 3, 4, 5, 6)))
    attrs = ((AT.BORD_ISO, "OUI"),)
    elrefe = (ElrefeLoc(MT.TR7, gauss=("RIGI=FPG6",), mater=("RIGI",)),)
