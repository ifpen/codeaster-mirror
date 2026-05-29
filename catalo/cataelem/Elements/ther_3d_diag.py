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

MMATTSR = ArrayOfComponents(phys=PHY.MTNS_R, locatedComponents=DDL_THER)


# --------------------------------------------------------------------------------------------------
class THER_HEXA8_D(Element):
    """Thermics - 3D_DIAG - HEXA8"""

    meshType = MT.HEXA8
    elrefe = (
        ElrefeLoc(
            MT.HE8,
            gauss=("RIGI=FPG8", "FPG1=FPG1", "MASS=NOEU_S", "NOEU=NOEU_S"),
            mater=("FPG1", "RIGI", "MASS"),
        ),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG4")),
    )
    calculs = (
        OP.CHAR_THER_EVOL(
            te=78,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_EVOL.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_EVOLNI(
            te=244,
            para_in=(
                (OP.CHAR_THER_EVOLNI.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PCAMASS, LC.CCAMA3D),
                (OP.CHAR_THER_EVOLNI.PHYDRPM, LC.EHYDRR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_EVOLNI.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTTI, MVECTTR), (SP.PVECTTR, MVECTTR)),
        ),
        OP.CHAR_THER_GRAI_F(
            te=217,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PGRAINF, LC.CFLUX3F),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_GRAI_F.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_GRAI_R(
            te=217,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PGRAINR, LC.CFLUX3R),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_THER_GRAI_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOURNL(
            te=354,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PSOURNL, LC.CSOURCF),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOUR_F(
            te=80,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PSOURCF, LC.CSOURCF),
                (SP.PINSTR, LC.CTIMETR),
                (OP.CHAR_THER_SOUR_F.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_SOUR_R(
            te=80,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PSOURCR, LC.ESOURCR)),
            para_out=((SP.PVECTTR, MVECTTR),),
        ),
        OP.CHAR_THER_TNL(
            te=525,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PLAGRM, LC.EGNEUT1R),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (SP.PVITESR, LC.NVITE3R),
            ),
            para_out=((SP.PLAGRP, LC.EGNEUT1R), (SP.PRESIDU, MVECTTR), (SP.PVECTTR, MVECTTR)),
        ),
        OP.COOR_ELGA(
            te=488,
            para_in=((SP.PGEOMER, LC.EGEOM3D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU3D),),
        ),
        OP.DURT_ELNO(
            te=551,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.DURT_ELNO.PPHASIN, LC.EPHASES)),
            para_out=((SP.PDURT_R, LC.EDURTNO),),
        ),
        OP.ERTH_ELEM(
            te=3,
            para_in=(
                (SP.PCHARG, LC.CREFERK),
                (SP.PFLUX_M, LC.NFLUX3R),
                (SP.PFLUX_P, LC.NFLUX3R),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PSOURCF, LC.CSOURCF),
                (SP.PSOURCR, LC.ESOURCR),
                (SP.PTEMP_M, DDL_THER),
                (SP.PTEMP_P, DDL_THER),
                (OP.ERTH_ELEM.PVOISIN, LC.EVOISIN),
            ),
            para_out=((OP.ERTH_ELEM.PERREUR, LC.EERREURT),),
        ),
        OP.ERTH_ELNO(
            te=379,
            para_in=((OP.ERTH_ELNO.PERREUR, LC.EERREURT),),
            para_out=((SP.PERRENO, LC.EERRENOT),),
        ),
        OP.ETHE_ELEM(
            te=66,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PINSTR, LC.CTIMETR),
                (OP.ETHE_ELEM.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.ETHE_ELEM.PENERDR, LC.CENEISO),),
        ),
        OP.FLUX_ELGA(
            te=69,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.FLUX_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.FLUX_ELGA.PFLUXPG, LC.EFLUX3R),),
        ),
        OP.FLUX_ELNO(
            te=4,
            para_in=((OP.FLUX_ELNO.PFLUXPG, LC.EFLUX3R),),
            para_out=((SP.PFLUXNO, LC.NFLUX3R),),
        ),
        OP.GRAT_ELGA(
            te=52,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PTEMPER, DDL_THER),
                (OP.GRAT_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.GRAT_ELGA.PGRATPG, LC.EGRAT3R),),
        ),
        OP.GRAT_ELNO(
            te=4,
            para_in=((OP.GRAT_ELNO.PGRATPG, LC.EGRAT3R),),
            para_out=((SP.PGRATNO, LC.NGRAT3R),),
        ),
        OP.HYDR_ELGA(
            te=385,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.HYDR_ELGA.PHYDRMR, LC.EHYDRR),
                (OP.HYDR_ELGA.PCOMPOR, LC.CCOMPOT),
                (SP.PMATERC, LC.CMATERC),
                (OP.HYDR_ELGA.PTEMPMR, DDL_THER),
                (OP.HYDR_ELGA.PTEMPPR, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((OP.HYDR_ELGA.PHYDRPR, LC.EHYDRR),),
        ),
        OP.HYDR_ELNO(
            te=4, para_in=((OP.HYDR_ELNO.PHYDRPG, LC.EHYDRR),), para_out=((SP.PHYDRNO, LC.EHYDRNO),)
        ),
        OP.INIT_MAIL_VOIS(te=99, para_out=((OP.INIT_MAIL_VOIS.PVOISIN, LC.EVOISIN),)),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_THER(
            te=77,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.CTIMETR),
                (OP.MASS_THER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
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
            te=321,
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
        OP.RIGI_THER_TANG(
            te=243,
            para_in=(
                (OP.RIGI_THER_TANG.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.RIGI_THER_TANG.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATTTR, MMATTTR), (OP.RIGI_THER_TANG.PMATTSR, MMATTSR)),
        ),
        OP.MASS_THER_TANG(
            te=246,
            para_in=(
                (OP.MASS_THER_TANG.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.MASS_THER_TANG.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.MTAN_THER_SOURNL(
            te=354,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PSOURNL, LC.CSOURCF),
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
                (SP.PCOEFR, LC.CNTINIR),
                (OP.NORME_L2.PCOORPG, LC.EGGAU3D),
            ),
            para_out=((SP.PNORME, LC.ENORME),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.REPERE_LOCAL(
            te=133,
            para_in=((SP.PCAMASS, LC.CCAMA3D), (SP.PGEOMER, LC.EGEOM3D)),
            para_out=((SP.PREPLO1, LC.CGEOM3D), (SP.PREPLO2, LC.CGEOM3D), (SP.PREPLO3, LC.CGEOM3D)),
        ),
        OP.RAPH_THER(
            te=243,
            para_in=(
                (OP.RAPH_THER.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.RAPH_THER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PRESIDU, MVECTTR), (OP.RAPH_THER.PFLUXPR, LC.EFLUX3R)),
        ),
        OP.MASS_THER_RESI(
            te=252,
            para_in=(
                (OP.MASS_THER_RESI.PCOMPOR, LC.CCOMPOT),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.MASS_THER_RESI.PHYDRPR, LC.EHYDRR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (OP.MASS_THER_RESI.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RESI_THER_SOURNL(
            te=354,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PSOURNL, LC.CSOURCF),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
            ),
            para_out=((SP.PRESIDU, MVECTTR),),
        ),
        OP.RIGI_THER(
            te=76,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.CTIMETR),
                (OP.RIGI_THER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.RIGI_THER_CONV(
            te=522,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (SP.PVITESR, LC.NVITE3R),
            ),
            para_out=((OP.RIGI_THER_CONV.PMATTTR, MMATTSR),),
        ),
        OP.RIGI_THER_TRANS(
            te=501,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPEI, DDL_THER),
                (SP.PTEMPER, DDL_THER),
            ),
            para_out=((SP.PMATTTR, MMATTTR),),
        ),
        OP.SOUR_ELGA(
            te=318,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTEMPER, DDL_THER),
                (SP.PINSTR, LC.CTIMETR),
                (OP.SOUR_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.SOUR_ELGA.PSOUR_R, LC.ESOURCR),),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),
                (OP.TOU_INI_ELEM.PCOEH_R, LC.CHECHPR),
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
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (OP.TOU_INI_ELGA.PSOUR_R, LC.ESOURCR),
                (OP.TOU_INI_ELGA.PVARI_R, LC.ZVARIPG),
                (OP.TOU_INI_ELGA.PHYDR_R, LC.EHYDRR),
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
class THER_PENTA6_D(THER_HEXA8_D):
    """Thermics - 3D_DIAG - PENTA6"""

    meshType = MT.PENTA6
    elrefe = (
        ElrefeLoc(
            MT.PE6,
            gauss=("RIGI=FPG6", "FPG1=FPG1", "MASS=NOEU_S", "NOEU=NOEU_S"),
            mater=("FPG1", "RIGI", "MASS"),
        ),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG4", "NOEU=NOEU")),
        ElrefeLoc(MT.TR3, gauss=("RIGI=COT3", "MASS=COT3", "NOEU=NOEU")),
    )


# --------------------------------------------------------------------------------------------------
class THER_TETRA4_D(THER_HEXA8_D):
    """Thermics - 3D_DIAG - TETRA4"""

    meshType = MT.TETRA4
    elrefe = (
        ElrefeLoc(
            MT.TE4,
            gauss=("RIGI=FPG4", "FPG1=FPG1", "MASS=NOEU_S", "NOEU=NOEU_S"),
            mater=("FPG1", "RIGI", "MASS"),
        ),
        ElrefeLoc(MT.TR3, gauss=("RIGI=COT3", "MASS=COT3")),
    )
