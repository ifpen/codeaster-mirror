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

NACCELR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY"))

MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)


# --------------------------------------------------------------------------------------------------
class SEAXSE2(Element):
    """Drying - Skin element AXIS - SEG2"""

    meshType = MT.SEG2
    elrefe = (ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("FPG1",)),)
    calculs = (
        OP.CHAR_THER_ACCE_R(
            te=315,
            para_in=((SP.PACCELR, NACCELR), (SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ACCE_X(
            te=315,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC), (SP.PTEMPER, DDL_THER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ACCE_Y(
            te=315,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC), (SP.PTEMPER, DDL_THER)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUNL(
            te=274,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_F(
            te=75,
            para_in=((SP.PFLUXNF, LC.CFLUXNF), (SP.PGEOMER, LC.EGEOM2D), (SP.PINSTR, LC.CTIMETR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUN_R(
            te=75,
            para_in=((SP.PFLUXNR, LC.CFLUXNR), (SP.PGEOMER, LC.EGEOM2D)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_FLUTNL(
            te=506,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.CHAR_THER_FLUX_F(
            te=75,
            para_in=((SP.PFLUXVF, LC.CFLUX2F), (SP.PGEOMER, LC.EGEOM2D), (SP.PINSTR, LC.CTIMETR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_ECHA_F(
            te=75,
            para_in=(
                (SP.PCOEFHF, LC.CHECHPF),
                (SP.PGEOMER, LC.EGEOM2D),
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
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (SP.PT_EXTR, LC.ET_EXTR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.COOR_ELGA(
            te=478,
            para_in=((SP.PGEOMER, LC.EGEOM2D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU2D),),
        ),
        OP.MTAN_THER_FLUXNL(
            te=251,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((OP.MTAN_THER_FLUXNL.PMATTTR, MMATTTR),),
        ),
        OP.RESI_THER_COEF_F(
            te=137,
            para_in=(
                (SP.PCOEFHF, LC.CHECHPF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_COEF_R(
            te=136,
            para_in=(
                (SP.PCOEFHR, LC.CHECHPR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_FLUXNL(
            te=138,
            para_in=(
                (SP.PFLUXNL, LC.CFLUXNF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RIGI_THER_ECHA_F(
            te=251,
            para_in=((SP.PCOEFHF, LC.CHECHPF), (SP.PGEOMER, LC.EGEOM2D), (SP.PINSTR, LC.CTIMETR)),
            para_out=((OP.RIGI_THER_ECHA_F.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_ECHA_R(
            te=251,
            para_in=((SP.PCOEFHR, LC.CHECHPR), (SP.PGEOMER, LC.EGEOM2D), (SP.PINSTR, LC.CTIMETR)),
            para_out=((OP.RIGI_THER_ECHA_R.PMATTTR, MMATTTR),),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PFLUX_R, LC.EFLUX2R),
                (OP.TOU_INI_ELGA.PGEOM_R, LC.EGGEO2D),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (SP.PTEMP_R, LC.ESECHPG),
            ),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D), (OP.TOU_INI_ELEM.PFLUN_R, LC.CFLUXNR)),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PFLUX_R, LC.NFLUX2R),
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM2D),
                (OP.TOU_INI_ELNO.PHYDR_R, LC.EHYDRNO),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PVARI_R, LC.EPHASES),
            ),
        ),
    )


# --------------------------------------------------------------------------------------------------
class SEAXSE3(SEAXSE2):
    """Drying - Skin element AXIS - SEG3"""

    meshType = MT.SEG3
    elrefe = (ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("FPG1",)),)


# --------------------------------------------------------------------------------------------------
class SEAXSL2(SEAXSE2):
    """Drying - Skin (diag) element AXIS - SEG2"""

    meshType = MT.SEG2
    elrefe = (ElrefeLoc(MT.SE2, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("FPG1",)),)
