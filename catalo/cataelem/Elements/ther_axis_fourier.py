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
DDL_THER = LocatedComponents(phys=PHY.TEMP_R, type="ELNO", components=("TEMP",))

MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)


# --------------------------------------------------------------------------------------------------
class THFOQU4(Element):
    """Thermics - AXIS_FOURIER - QUAD4"""

    meshType = MT.QUAD4
    elrefe = (
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "FPG1=FPG1", "MASS=FPG4"), mater=("FPG1",)),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )
    calculs = (
        OP.CARA_CISA(te=-1),
        OP.CARA_GAUCHI(te=-1),
        OP.CARA_TORSION(te=-1),
        OP.CHAR_THER_SOUR_F(
            te=264,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PSOURCF, LC.CSOURCF), (SP.PINSTR, LC.CTIMETR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOUR_R(
            te=263,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PSOURCR, LC.ESOURCR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.COOR_ELGA(
            te=479,
            para_in=((SP.PGEOMER, LC.EGEOM2D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU2D),),
        ),
        OP.DURT_ELNO(
            te=551,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.DURT_ELNO.PPHASIN, LC.EPHASES)),
            para_out=((SP.PDURT_R, LC.EDURTNO),),
        ),
        OP.FLUX_ELGA(
            te=266,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PHARMON, LC.CHARMON),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((OP.FLUX_ELGA.PFLUXPG, LC.EFLUX3R),),
        ),
        OP.FLUX_ELNO(
            te=4,
            para_in=((OP.FLUX_ELNO.PFLUXPG, LC.EFLUX3R),),
            para_out=((SP.PFLUXNO, LC.NFLUX3R),),
        ),
        OP.META_ELNO(
            te=64,
            para_in=(
                (SP.PCOMPME, LC.CCOMPOT),
                (SP.PCOMPMT, LC.CCOMPOT),
                (SP.PFTRC, LC.CFTRC),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPHASIN, LC.EPHASES),
                (SP.PTEMPAR, DDL_THER),
                (SP.PTEMPER, DDL_THER),
                (SP.PTEMPIR, DDL_THER),
                (SP.PTIMMTR, LC.CTIMMTR),
                (SP.PPHASEP, LC.EPHASES),
            ),
            para_out=((OP.META_ELNO.PPHASOUT, LC.EPHASES),),
        ),
        OP.META_INIT_ELNO(
            te=320,
            para_in=(
                (SP.PCOMPME, LC.CCOMPOT),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPHASII, LC.CPHASES),
                (SP.PTEMPER, DDL_THER),
            ),
            para_out=((OP.META_INIT_ELNO.PPHASOUT, LC.EPHASES),),
        ),
        OP.META_TRAN_ELNO(
            te=319,
            para_in=((OP.META_TRAN_ELNO.PPHASIN, LC.EPHASES),),
            para_out=((OP.META_TRAN_ELNO.PPHASOUT, LC.EPHASES),),
        ),
        OP.RIGI_THER(
            te=260,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PHARMON, LC.CHARMON),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.TOU_INI_ELEM(te=99, para_out=((OP.TOU_INI_ELEM.PSOUR_R, LC.CSOURCR),)),
        OP.TOU_INI_ELGA(te=99, para_out=((OP.TOU_INI_ELGA.PGEOM_R, LC.EGGEO2D),)),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM2D),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PVARI_R, LC.EPHASES),
            ),
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, LC.EGEOM2D),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# --------------------------------------------------------------------------------------------------
class THFOQU8(THFOQU4):
    """Thermics - AXIS_FOURIER - QUAD8"""

    meshType = MT.QUAD8
    elrefe = (
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG9", "FPG1=FPG1", "MASS=FPG9"), mater=("FPG1",)),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )


# --------------------------------------------------------------------------------------------------
class THFOQU9(THFOQU4):
    """Thermics - AXIS_FOURIER - QUAD9"""

    meshType = MT.QUAD9
    elrefe = (
        ElrefeLoc(MT.QU9, gauss=("RIGI=FPG9", "FPG1=FPG1", "MASS=FPG9"), mater=("FPG1",)),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )


# --------------------------------------------------------------------------------------------------
class THFOTR3(THFOQU4):
    """Thermics - AXIS_FOURIER - TRIA3"""

    meshType = MT.TRIA3
    elrefe = (
        ElrefeLoc(MT.TR3, gauss=("RIGI=FPG1", "FPG1=FPG1", "MASS=FPG3"), mater=("FPG1",)),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )


# --------------------------------------------------------------------------------------------------
class THFOTR6(THFOQU4):
    """Thermics - AXIS_FOURIER - TRIA6"""

    meshType = MT.TRIA6
    elrefe = (
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG3", "FPG1=FPG1", "MASS=FPG6"), mater=("FPG1",)),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )
