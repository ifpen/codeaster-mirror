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
DDL_THER = LocatedComponents(
    phys=PHY.TEMP_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("HHO_FT[10]",)), ("EN2", ()), ("EN3", ("HHO_CT[20]"))),
)

CPPOST = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X[1]"))
TEMPHHO = LocatedComponents(phys=PHY.TEMP_R, type="ELNO", components=("TEMP",))

# too large to be saved so used size 1
CHHOGT = LocatedComponents(phys=PHY.N1920R, type="ELEM", components=("X[1]",))

CHHOST = LocatedComponents(phys=PHY.N1360R, type="ELEM", components=("X[1]",))

CHHOBS = LocatedComponents(
    phys=PHY.N3600R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("X[55]",)), ("EN2", ()), ("EN3", ("X[630]"))),
)

PFONC = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z[7]",))

PFONCR = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z[2]",))

MVECTTR = ArrayOfComponents(phys=PHY.VTEM_R, locatedComponents=DDL_THER)

MMATTTR = ArrayOfComponents(phys=PHY.MTEM_R, locatedComponents=DDL_THER)

MMATTSR = ArrayOfComponents(phys=PHY.MTNS_R, locatedComponents=DDL_THER)


# --------------------------------------------------------------------------------------------------
class THER3DH27_HHO333(Element):
    """Thermics - HHO_QUAD - 3D - HEXA"""

    meshType = MT.HEXA27
    nodes = (
        SetOfNodes("EN1", (21, 22, 23, 24, 25, 26)),
        SetOfNodes("EN2", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN3", (27,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.H27,
            gauss=("RIGI=FPG64", "FPG1=FPG1", "MASS=FPG64", "MTGA=FPG64"),
            mater=("RIGI", "FPG1", "MTGA", "MASS"),
        ),
    )
    calculs = (
        OP.CHAR_THER_EVOL(
            te=445,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_EVOL.PVARCPR, LC.ZVARCPG),
                (OP.CHAR_THER_EVOL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOUR_F(
            te=465,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PSOURCF, LC.CSOURCF),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_SOUR_F.PVARCPR, LC.ZVARCPG),
                (OP.CHAR_THER_SOUR_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOURNL(
            te=465,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PSOURNL, LC.CSOURCF),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_SOURNL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOUR_R(
            te=465,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PSOURCR, LC.ESOURCR),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_SOUR_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.COOR_ELGA(
            te=488,
            para_in=((SP.PGEOMER, LC.EGEOM3D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU3D),),
        ),
        OP.FLUX_ELGA(
            te=487,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.FLUX_ELGA.PVARCPR, LC.ZVARCPG),
                (OP.FLUX_ELGA.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.FLUX_ELGA.PFLUXPG, LC.EFLUX3R),),
        ),
        OP.FLUX_ELNO(
            te=4,
            para_in=((OP.FLUX_ELNO.PFLUXPG, LC.EFLUX3R),),
            para_out=((SP.PFLUXNO, LC.NFLUX3R),),
        ),
        OP.HHO_PRECALC_BS(
            te=494,
            para_in=((SP.PGEOMER, LC.EGEOM3D),),
            para_out=((OP.HHO_PRECALC_BS.PCHHOBO, CHHOBS),),
        ),
        OP.HHO_PRECALC_OP(
            te=460,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (OP.HHO_PRECALC_OP.PCHHOBS, CHHOBS)),
            para_out=((OP.HHO_PRECALC_OP.PCHHOGT, CHHOGT), (OP.HHO_PRECALC_OP.PCHHOST, CHHOST)),
        ),
        OP.HHO_PROJ_THER(
            te=473,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.HHO_PROJ_THER.PFUNC_R, PFONCR),
                (SP.PINSTPR, LC.MTEMPSR),
                (OP.HHO_PROJ_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_PROJ_THER.PTEMP_R, DDL_THER),),
        ),
        OP.HHO_PROJ2_THER(
            te=484,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.HHO_PROJ2_THER.PH1TP_R, TEMPHHO),
                (OP.HHO_PROJ2_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_PROJ_THER.PTEMP_R, DDL_THER),),
        ),
        OP.HHO_PROJ3_THER(
            te=484,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.HHO_PROJ3_THER.PQPTP_R, LC.ETEMPPG),
                (OP.HHO_PROJ3_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_PROJ3_THER.PTEMP_R, DDL_THER),),
        ),
        OP.HHO_TEMP_THER(
            te=456,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTMPCHF, DDL_THER),
                (OP.HHO_TEMP_THER.PCHHOBS, CHHOBS),
                (SP.POPPOST, CPPOST),
            ),
            para_out=((OP.HHO_TEMP_THER.PTEMP_R, TEMPHHO),),
        ),
        OP.HHO_CINE_F_THER(
            te=492,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTPR, LC.MTEMPSR),
                (OP.HHO_CINE_F_THER.PFONC, PFONC),
                (OP.HHO_CINE_F_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_CINE_F_THER.PCINE, DDL_THER),),
        ),
        OP.HHO_CINE_R_THER(
            te=492,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.HHO_CINE_R_THER.PCMPVALE, TEMPHHO),
                (OP.HHO_CINE_R_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_CINE_R_THER.PCINE, DDL_THER),),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_THER(
            te=449,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MASS_THER.PVARCPR, LC.ZVARCPG),
                (OP.MASS_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MASS_THER_TANG(
            te=429,
            para_in=(
                (OP.MASS_THER_TANG.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.MASS_THER_TANG.PVARCPR, LC.ZVARCPG),
                (OP.MASS_THER_TANG.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MASS_THER_RESI(
            te=429,
            para_in=(
                (OP.MASS_THER_RESI.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.MASS_THER_RESI.PHYDRPR, LC.EHYDRR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.MASS_THER_RESI.PVARCPR, LC.ZVARCPG),
                (OP.MASS_THER_RESI.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.MTAN_THER_SOURNL(
            te=437,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PSOURNL, LC.CSOURCF),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MTAN_THER_SOURNL.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.RAPH_THER(
            te=477,
            para_in=(
                (OP.RAPH_THER.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.RAPH_THER.PVARCPR, LC.ZVARCPG),
                (OP.RAPH_THER.PCHHOGT, CHHOGT),
                (OP.RAPH_THER.PCHHOST, CHHOST),
                (OP.RAPH_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PRESIDU, MVECTTR), (OP.RAPH_THER.PFLUXPR, LC.EFLUX3R)),
        ),
        OP.RIGI_THER(
            te=454,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.CTIMETR),
                (OP.RIGI_THER.PVARCPR, LC.ZVARCPG),
                (OP.RIGI_THER.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_TANG(
            te=454,
            para_in=(
                (OP.RIGI_THER_TANG.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.RIGI_THER_TANG.PVARCPR, LC.ZVARCPG),
                (OP.RIGI_THER_TANG.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.RIGI_THER_TANG.PMATTSR, MMATTSR),),
        ),
        OP.TEMP_ELGA(
            te=456,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPER, DDL_THER),
                (OP.TEMP_ELGA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PTEMP_R, LC.ETEMPPG),),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D), (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER)),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PCOEH_R, LC.EHECHPR),
                (OP.TOU_INI_ELEM.PSOUR_R, LC.CSOURCR),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
            ),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PFLUX_R, LC.EFLUX3R),
                (OP.TOU_INI_ELGA.PGEOM_R, LC.EGGEO3D),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGNEUT1R),
                (OP.TOU_INI_ELGA.PSOUR_R, LC.ESOURCR),
                (OP.TOU_INI_ELGA.PVARI_R, LC.ZVARIPG),
                (SP.PTEMP_R, LC.ETEMPPG),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PFLUX_R, LC.NFLUX3R),
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM3D),
                (OP.TOU_INI_ELNO.PHYDR_R, LC.EHYDRNO),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PVARI_R, LC.EPHASES),
                (OP.TOU_INI_ELNO.PSOUR_R, LC.NSOURCR),
            ),
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, LC.EGEOM3D),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# --------------------------------------------------------------------------------------------------
class THER3DT15_HHO333(THER3DH27_HHO333):
    """Thermics - HHO_QUAD - 3D - TETRA"""

    meshType = MT.TETRA15
    nodes = (
        SetOfNodes("EN1", (11, 12, 13, 14)),
        SetOfNodes("EN2", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
        SetOfNodes("EN3", (15,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.T15,
            gauss=("RIGI=FPG24", "FPG1=FPG1", "MASS=FPG24", "MTGA=FPG24"),
            mater=("RIGI", "FPG1", "MTGA", "MASS"),
        ),
    )


# --------------------------------------------------------------------------------------------------
class THER3DP21_HHO333(THER3DH27_HHO333):
    """Thermics - HHO_QUAD - 3D - PENTA"""

    meshType = MT.PENTA21
    nodes = (
        SetOfNodes("EN1", (16, 17, 18, 19, 20)),
        SetOfNodes("EN2", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN3", (21,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.P21,
            gauss=("RIGI=FPG29", "FPG1=FPG1", "MTGA=FPG29", "MASS=FPG29"),
            mater=("RIGI", "FPG1", "MTGA", "MASS"),
        ),
    )


# --------------------------------------------------------------------------------------------------
class THER3DP19_HHO333(THER3DH27_HHO333):
    """Thermics - HHO_QUAD - 3D - PYRAM"""

    meshType = MT.PYRAM19
    nodes = (
        SetOfNodes("EN1", (14, 15, 16, 17, 18)),
        SetOfNodes("EN2", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)),
        SetOfNodes("EN3", (19,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.P19,
            gauss=("RIGI=FPG24", "FPG1=FPG1", "MTGA=FPG24", "MASS=FPG24"),
            mater=("RIGI", "FPG1", "MTGA", "MASS"),
        ),
    )
