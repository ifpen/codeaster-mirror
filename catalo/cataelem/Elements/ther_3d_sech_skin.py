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

# --------------------------------------------------------------------------------------------------
# Located components
# --------------------------------------------------------------------------------------------------
DDL_THER = LocatedComponents(phys=PHY.TEMP_R, type="ELNO", components=("SECH",))

NACCELR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))

MVECTAR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=NACCELR)

MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)


# --------------------------------------------------------------------------------------------------
class SECH_FACE3(Element):
    """Drying - Skin element 3D - TRIA3"""

    meshType = MT.TRIA3
    elrefe = (ElrefeLoc(MT.TR3, gauss=("RIGI=COT3", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),)
    calculs = (
        OP.ACCEPTANCE(
            te=329,
            para_in=((SP.PACCELR, NACCELR), (SP.PGEOMER, LC.EGEOM3D), (SP.PNUMMOD, LC.CNUMMOD)),
            para_out=((SP.PVECTUR, MVECTAR),),
        ),
        OP.AMOR_AJOU(
            te=327,
            para_in=((SP.PACCELR, NACCELR), (SP.PGEOMER, LC.EGEOM3D)),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.CHAR_THER_ACCE_R(
            te=325,
            para_in=((SP.PACCELR, NACCELR), (SP.PGEOMER, LC.EGEOM3D), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ACCE_X(
            te=325,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PMATERC, LC.CMATERC), (SP.PTEMPER, DDL_THER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ACCE_Y(
            te=325,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PMATERC, LC.CMATERC), (SP.PTEMPER, DDL_THER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ACCE_Z(
            te=325,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PMATERC, LC.CMATERC), (SP.PTEMPER, DDL_THER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUNL(
            te=274,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_F(
            te=75,
            para_in=((SP.PFLUXNF, LC.CFLUXNF), (SP.PGEOMER, LC.EGEOM3D), (SP.PINSTR, LC.CTIMETR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_R(
            te=75,
            para_in=((SP.PFLUXNR, LC.CFLUXNR), (SP.PGEOMER, LC.EGEOM3D)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUTNL(
            te=506,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.CHAR_THER_FLUX_F(
            te=75,
            para_in=((SP.PFLUXVF, LC.CFLUX3F), (SP.PGEOMER, LC.EGEOM3D), (SP.PINSTR, LC.CTIMETR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_PHID_R(
            te=326,
            para_in=(
                (SP.PACCELR, NACCELR),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ECHA_F(
            te=75,
            para_in=(
                (SP.PCOEFHF, LC.CHECHPF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (SP.PT_EXTF, LC.CTEMPEF),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ECHA_R(
            te=75,
            para_in=(
                (SP.PCOEFHR, LC.CHECHPR),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (SP.PT_EXTR, LC.ET_EXTR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.COOR_ELGA(
            te=488,
            para_in=((SP.PGEOMER, LC.EGEOM3D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU3D),),
        ),
        OP.MTAN_THER_FLUXNL(
            te=251,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.NORME_L2(
            te=563,
            para_in=(
                (SP.PCALCI, LC.EMNEUT_I),
                (SP.PCHAMPG, LC.EGTINIR),
                (SP.PCOEFR, LC.CNORMCF),
                (OP.NORME_L2.PCOORPG, LC.EGGEO3D),
            ),
            para_out=((SP.PNORME, LC.ENORME),),
        ),
        OP.RESI_THER_COEF_F(
            te=128,
            para_in=(
                (SP.PCOEFHF, LC.CHECHPF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_COEF_R(
            te=127,
            para_in=(
                (SP.PCOEFHR, LC.CHECHPR),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_FLUXNL(
            te=129,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RIGI_THER_ECHA_F(
            te=251,
            para_in=((SP.PCOEFHF, LC.CHECHPF), (SP.PGEOMER, LC.EGEOM3D), (SP.PINSTR, LC.CTIMETR)),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_ECHA_R(
            te=251,
            para_in=((SP.PCOEFHR, LC.CHECHPR), (SP.PGEOMER, LC.EGEOM3D), (SP.PINSTR, LC.CTIMETR)),
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
class SECH_FACE4(SECH_FACE3):
    """Drying - Skin element 3D - QUAD4"""

    meshType = MT.QUAD4
    elrefe = (ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),)


# --------------------------------------------------------------------------------------------------
class SECH_FACE6(SECH_FACE3):
    """Drying - Skin element 3D - TRIA6"""

    meshType = MT.TRIA6
    elrefe = (ElrefeLoc(MT.TR6, gauss=("RIGI=FPG6", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),)


# --------------------------------------------------------------------------------------------------
class SECH_FACE8(SECH_FACE3):
    """Drying - Skin element 3D - QUAD8"""

    meshType = MT.QUAD8
    elrefe = (ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),)


# --------------------------------------------------------------------------------------------------
class SECH_FACE9(SECH_FACE3):
    """Drying - Skin element 3D - QUAD9"""

    meshType = MT.QUAD9
    elrefe = (ElrefeLoc(MT.QU9, gauss=("RIGI=FPG9", "NOEU=NOEU", "FPG1=FPG1"), mater=("FPG1",)),)
