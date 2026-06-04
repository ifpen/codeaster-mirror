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


# Elementary characteristics for pipes: radius and thickness
CCAGEPO = LocatedComponents(phys=PHY.CAGEPO_R, type="ELEM", components=("R1", "EP1"))

# Sub-points for pipes: layers and sectors
ENBSP_I = LocatedComponents(phys=PHY.NBSP_I, type="ELEM", components=("TUY_NCOU", "TUY_NSEC"))

# Strains for pipes
EDEFGNO = LocatedComponents(
    phys=PHY.EPSI_R, type="ELNO", components=("EPXX", "GAXY", "GAXZ", "GAT", "KY", "KZ")
)

EDEFGPG = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "GAXY", "GAXZ", "GAT", "KY", "KZ"),
)

# Curvilinear coordinates (for SEG3)
CABSCUR = LocatedComponents(phys=PHY.ABSC_R, type="ELEM", components=("ABSC[3]",))


# Local orientation for pipes
CCAORIE = LocatedComponents(
    phys=PHY.CAORIE_R,
    type="ELEM",
    components=(
        "ALPHA",
        "BETA",
        "GAMMA",
        "ALPHA2",
        "BETA2",
        "GAMMA2",
        "ALPHA3",
        "BETA3",
        "GAMMA3",
        "ICOUDE",
        "DN1N2",
        "RCOURB",
        "ANGCOU",
        "ANGZZK",
    ),
)

NDEPLAC = LocatedComponents(
    phys=PHY.DEPL_C,
    type="ELNO",
    components=(
        "DX",
        "DY",
        "DZ",
        "DRX",
        "DRY",
        "DRZ",
        "UI2",
        "VI2",
        "WI2",
        "UI3",
        "VI3",
        "WI3",
        "UO2",
        "VO2",
        "WO2",
        "UO3",
        "VO3",
        "WO3",
        "WO",
        "WI1",
        "WO1",
    ),
)

DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    components=(
        "DX",
        "DY",
        "DZ",
        "DRX",
        "DRY",
        "DRZ",
        "UI2",
        "VI2",
        "WI2",
        "UI3",
        "VI3",
        "WI3",
        "UO2",
        "VO2",
        "WO2",
        "UO3",
        "VO3",
        "WO3",
        "WO",
        "WI1",
        "WO1",
    ),
)

MVECTUC = ArrayOfComponents(phys=PHY.VDEP_C, locatedComponents=NDEPLAC)
MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)
MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)
MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

CABSCUR.setName("CABSCUR")
CCAORIE.setName("CCAORIE")
NDEPLAC.setName("NDEPLAC")
DDL_MECA.setName("DDL_MECA")
MVECTUC.setName("MVECTUC")
MVECTUR.setName("MVECTUR")
MMATUUR.setName("MMATUUR")
MMATUUC.setName("MMATUUC")


class MET3SEG3(Element):
    """Mechanics - TUYAU_3M - SEG3"""

    meshType = MT.SEG3
    elrefe = (
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG2", "MASS=FPG3", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
    )

    calculs = (
        OP.AMOR_MECA(
            te=121,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMASSEL, MMATUUR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PRIGIEL, MMATUUR),
                (OP.AMOR_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.CHAR_MECA_FC1D1D(
            te=583,
            para_in=((SP.PCAORIE, CCAORIE), (SP.PFC1D1D, LC.CFORCEC), (SP.PGEOMER, LC.EGEOM3D)),
            para_out=((SP.PVECTUC, MVECTUC),),
        ),
        OP.CHAR_MECA_FF1D1D(
            te=583,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PFF1D1D, LC.CFORPIF),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.CHAR_MECA_FF1D1D.PNBSP_I, ENBSP_I),
                (SP.PINSTR, LC.MTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FR1D1D(
            te=583,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PFR1D1D, LC.CFORPIR),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.CHAR_MECA_FR1D1D.PNBSP_I, ENBSP_I),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=583,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_PESA_R.PNBSP_I, ENBSP_I),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_F(
            te=583,
            para_in=(
                (SP.PABSCUR, CABSCUR),
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.CHAR_MECA_PRES_F.PNBSP_I, ENBSP_I),
                (SP.PPRESSF, LC.CPRE3DF),
                (SP.PINSTR, LC.MTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=583,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PGEOMER, LC.CGEOM3D),
                (OP.CHAR_MECA_PRES_R.PNBSP_I, ENBSP_I),
                (SP.PPRESSR, LC.EPRE3DR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=589,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_TEMP_R.PNBSP_I, ENBSP_I),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.COOR_ELGA(
            te=478,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.COOR_ELGA.PNBSP_I, ENBSP_I),
            ),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU3D), (OP.COOR_ELGA.PCOORSU, LC.EGGAU3D)),
        ),
        OP.DEGE_ELGA(
            te=584,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.DEGE_ELGA.PNBSP_I, ENBSP_I),
                (OP.DEGE_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.DEGE_ELGA.PDEFOPG, EDEFGPG),),
        ),
        OP.DEGE_ELNO(
            te=584,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.DEGE_ELNO.PNBSP_I, ENBSP_I),
                (OP.DEGE_ELNO.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PDEFOGR, EDEFGNO),),
        ),
        OP.EFGE_ELGA(
            te=587,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (OP.EFGE_ELGA.PNBSP_I, ENBSP_I),
                (SP.PSIEFR, LC.EGIG3DR),
            ),
            para_out=((SP.PEFGEC, LC.EEFGEBC), (SP.PEFGER, LC.EEFGEBR)),
        ),
        OP.EFGE_ELNO(
            te=185,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (OP.EFGE_ELNO.PCOMPOR, LC.CCOMPOR),
                (OP.EFGE_ELNO.PCONTRR, LC.EGIG3DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.EFGE_ELNO.PNBSP_I, ENBSP_I),
                (SP.PNONLIN, LC.ENONLIN),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EFGE_ELNO.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PEFFORC, LC.NEFGEBC), (OP.EFGE_ELNO.PEFFORR, LC.NEFGEBR)),
        ),
        OP.EFEQ_ELNO(
            te=83, para_in=((SP.PEFFONR, LC.NEFGEBR),), para_out=((SP.PEFFOENR, LC.EEFGENOQ),)
        ),
        OP.EPEQ_ELGA(
            te=335,
            para_in=((OP.EPEQ_ELGA.PDEFORR, LC.EGPS3DR),),
            para_out=((OP.EPEQ_ELGA.PDEFOEQ, LC.EDFEQPG),),
        ),
        OP.EPEQ_ELNO(
            te=335,
            para_in=((OP.EPEQ_ELNO.PDEFORR, LC.EEPS3DR),),
            para_out=((OP.EPEQ_ELNO.PDEFOEQ, LC.EDFEQNO),),
        ),
        OP.EPME_ELGA(
            te=531,
            para_in=(
                (OP.EPME_ELGA.PDEFORR, LC.EGPS3DR),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPME_ELGA.PNBSP_I, ENBSP_I),
                (OP.EPME_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPME_ELGA.PDEFOPG, LC.EGPS3DR),),
        ),
        OP.EPME_ELNO(
            te=4,
            para_in=((OP.EPME_ELNO.PDEFOPG, LC.EGPS3DR),),
            para_out=((SP.PDEFONO, LC.EEPS3DR),),
        ),
        OP.EPSI_ELGA(
            te=584,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PDEPLAR, DDL_MECA),
                (OP.EPSI_ELGA.PNBSP_I, ENBSP_I),
            ),
            para_out=((SP.PDEFOPC, LC.EGPS3DC), (OP.EPSI_ELGA.PDEFOPG, LC.EGPS3DR)),
        ),
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, LC.EGPS3DR),),
            para_out=((SP.PDEFONC, LC.EEPS3DC), (SP.PDEFONO, LC.EEPS3DR)),
        ),
        OP.EPSP_ELGA(
            te=531,
            para_in=(
                (OP.EPSP_ELGA.PCONTRR, LC.EGIG3DR),
                (OP.EPSP_ELGA.PDEFORR, LC.EGPS3DR),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPSP_ELGA.PNBSP_I, ENBSP_I),
                (OP.EPSP_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPSP_ELGA.PDEFOPG, LC.EGPS3DR),),
        ),
        OP.EPSP_ELNO(
            te=4,
            para_in=((OP.EPSP_ELNO.PDEFOPG, LC.EGPS3DR),),
            para_out=((SP.PDEFONO, LC.EEPS3DR),),
        ),
        OP.EPVC_ELGA(
            te=531,
            para_in=(
                (OP.EPVC_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPVC_ELGA.PNBSP_I, ENBSP_I),
                (OP.EPVC_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPVC_ELGA.PDEFOPG, LC.EGVARC1D),),
        ),
        OP.EPVC_ELNO(
            te=4,
            para_in=((OP.EPVC_ELNO.PDEFOPG, LC.EGVARC1D),),
            para_out=((SP.PDEFONO, LC.NVARC1D),),
        ),
        OP.FORC_NODA(
            te=585,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, LC.EGIG3DR),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.FORC_NODA.PNBSP_I, ENBSP_I),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=586,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, LC.EGIG3DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (OP.FULL_MECA.PNBSP_I, ENBSP_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, LC.ZVARIPG),
                (OP.FULL_MECA.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, LC.EGIG3DR),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_INER(
            te=38,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_INER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMASSINE, LC.EMASSINE),),
        ),
        OP.MASS_MECA(
            te=582,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PNBSP_I, ENBSP_I),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MINMAX_SP(te=99, para_out=((SP.PGAMIMA, LC.EGMINMAX), (SP.PNOMIMA, LC.NMINMAX))),
        OP.M_GAMMA(
            te=582,
            para_in=(
                (SP.PACCELR, DDL_MECA),
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.M_GAMMA.PNBSP_I, ENBSP_I),
                (OP.M_GAMMA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2), (OP.NSPG_NBVA.PNBSP_I, ENBSP_I)),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.PAS_COURANT(
            te=404,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.PAS_COURANT.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PCOURAN, LC.ECOURAN),),
        ),
        OP.RAPH_MECA(
            te=586,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, LC.EGIG3DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (OP.RAPH_MECA.PNBSP_I, ENBSP_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, LC.ZVARIPG),
                (OP.RAPH_MECA.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, LC.EGIG3DR),
                (OP.RAPH_MECA.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.REFE_FORC_NODA(
            te=585,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (OP.REFE_FORC_NODA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.REFE_FORC_NODA.PNBSP_I, ENBSP_I),
                (SP.PREFCO, LC.CRESSIG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.REPERE_LOCAL(
            te=135,
            para_in=((SP.PCAORIE, CCAORIE),),
            para_out=((SP.PREPLO1, LC.CGEOM3D), (SP.PREPLO2, LC.CGEOM3D), (SP.PREPLO3, LC.CGEOM3D)),
        ),
        OP.RIGI_MECA(
            te=582,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PNBSP_I, ENBSP_I),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_HYST(
            te=50,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PRIGIEL, MMATUUR),
                (OP.RIGI_MECA_HYST.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUC, MMATUUC),),
        ),
        OP.RIGI_MECA_TANG(
            te=586,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PCAORIE, CCAORIE),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, LC.EGIG3DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA_TANG.PNBSP_I, ENBSP_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, LC.ZVARIPG),
                (OP.RIGI_MECA_TANG.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, LC.EGIG3DR),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.SIEF_ELGA(
            te=584,
            para_in=(
                (SP.PCAGEPO, CCAGEPO),
                (SP.PCAORIE, CCAORIE),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.SIEF_ELGA.PNBSP_I, ENBSP_I),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.SIEF_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PCONTRC, LC.EGIG3DC), (OP.SIEF_ELGA.PCONTRR, LC.EGIG3DR)),
        ),
        OP.SIEF_ELNO(
            te=4,
            para_in=((OP.SIEF_ELNO.PCONTRR, LC.EGIG3DR), (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PSIEFNOC, LC.ESIG3DC), (OP.SIEF_ELNO.PSIEFNOR, LC.ESIG3DR)),
        ),
        OP.SIEQ_ELGA(
            te=335,
            para_in=((OP.SIEQ_ELGA.PCONTRR, LC.EGIG3DR),),
            para_out=((OP.SIEQ_ELGA.PCONTEQ, LC.ECOEQPG),),
        ),
        OP.SIEQ_ELNO(
            te=335,
            para_in=((OP.SIEQ_ELNO.PCONTRR, LC.ESIG3DR),),
            para_out=((OP.SIEQ_ELNO.PCONTEQ, LC.ECOEQNO),),
        ),
        OP.SIGM_ELGA(
            te=546,
            para_in=((SP.PSIEFR, LC.EGIG3DR),),
            para_out=((SP.PSIGMC, LC.EGIG3DC), (SP.PSIGMR, LC.EGIG3DR)),
        ),
        OP.SIGM_ELNO(
            te=4,
            para_in=((OP.SIGM_ELNO.PCONTRR, LC.EGIG3DR),),
            para_out=((SP.PSIEFNOC, LC.ESIG3DC), (OP.SIGM_ELNO.PSIEFNOR, LC.ESIG3DR)),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),
                (OP.TOU_INI_ELEM.PNBSP_I, ENBSP_I),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
            ),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (OP.TOU_INI_ELGA.PSIEF_R, LC.EGIG3DR),
                (OP.TOU_INI_ELGA.PVARI_R, LC.ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM3D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PPRES_R, LC.EPRE3DR),
                (OP.TOU_INI_ELNO.PSIEF_R, LC.ESIG3DR),
            ),
        ),
        OP.VARI_ELNO(
            te=4,
            para_in=((SP.PVARIGR, LC.ZVARIPG),),
            para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),),
        ),
        OP.VERI_CARA_ELEM(
            te=119,
            para_in=((SP.PCAGEPO, CCAGEPO), (SP.PCAORIE, CCAORIE), (SP.PCHCKPR, LC.CCHCKPR)),
            para_out=((SP.PCODRET, LC.ECODRET), (SP.PINDICR, LC.CINDICR)),
        ),
    )


class MET3SEG4(MET3SEG3):
    """Mechanics - TUYAU_3M - SEG4"""

    meshType = MT.SEG4
    elrefe = (
        ElrefeLoc(MT.SE4, gauss=("RIGI=FPG3", "MASS=FPG3", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
    )

    def postInit(self):
        """Redefine components that differ from the parent element"""
        self.changeComponents("CABSCUR", ("ABSC[4]",))
        self.changeComponents(
            "CCAORIE",
            (
                "ALPHA",
                "BETA",
                "GAMMA",
                "ALPHA2",
                "BETA2",
                "GAMMA2",
                "ALPHA3",
                "BETA3",
                "GAMMA3",
                "ALPHA4",
                "BETA4",
                "GAMMA4",
                "ICOUDE",
                "DN1N2",
                "RCOURB",
                "ANGCOU",
                "ANGZZK",
            ),
        )
        self.changeComponents(
            "NDEPLAC",
            (
                "DX",
                "DY",
                "DZ",
                "DRX",
                "DRY",
                "DRZ",
                "UI2",
                "VI2",
                "WI2",
                "UI3",
                "VI3",
                "WI3",
                "UO2",
                "VO2",
                "WO2",
                "UO3",
                "VO3",
                "WO3",
                "WO",
                "WI1",
                "WO1",
            ),
        )
        self.changeComponents(
            "DDL_MECA",
            (
                "DX",
                "DY",
                "DZ",
                "DRX",
                "DRY",
                "DRZ",
                "UI2",
                "VI2",
                "WI2",
                "UI3",
                "VI3",
                "WI3",
                "UO2",
                "VO2",
                "WO2",
                "UO3",
                "VO3",
                "WO3",
                "WO",
                "WI1",
                "WO1",
            ),
        )


class MET6SEG3(MET3SEG3):
    """Mechanics - TUYAU_6M - SEG3"""

    meshType = MT.SEG3
    elrefe = (
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG2", "MASS=FPG3", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
    )

    def postInit(self):
        """Redefine components that differ from the parent element"""
        self.changeComponents("CABSCUR", ("ABSC[3]",))
        self.changeComponents(
            "CCAORIE",
            (
                "ALPHA",
                "BETA",
                "GAMMA",
                "ALPHA2",
                "BETA2",
                "GAMMA2",
                "ALPHA3",
                "BETA3",
                "GAMMA3",
                "ICOUDE",
                "DN1N2",
                "RCOURB",
                "ANGCOU",
                "ANGZZK",
            ),
        )
        self.changeComponents(
            "NDEPLAC",
            (
                "DX",
                "DY",
                "DZ",
                "DRX",
                "DRY",
                "DRZ",
                "UI2",
                "VI2",
                "WI2",
                "UI3",
                "VI3",
                "WI3",
                "UI4",
                "VI4",
                "WI4",
                "UI5",
                "VI5",
                "WI5",
                "UI6",
                "VI6",
                "WI6",
                "UO2",
                "VO2",
                "WO2",
                "UO3",
                "VO3",
                "WO3",
                "UO4",
                "VO4",
                "WO4",
                "UO5",
                "VO5",
                "WO5",
                "UO6",
                "VO6",
                "WO6",
                "WO",
                "WI1",
                "WO1",
            ),
        )
        self.changeComponents(
            "DDL_MECA",
            (
                "DX",
                "DY",
                "DZ",
                "DRX",
                "DRY",
                "DRZ",
                "UI2",
                "VI2",
                "WI2",
                "UI3",
                "VI3",
                "WI3",
                "UI4",
                "VI4",
                "WI4",
                "UI5",
                "VI5",
                "WI5",
                "UI6",
                "VI6",
                "WI6",
                "UO2",
                "VO2",
                "WO2",
                "UO3",
                "VO3",
                "WO3",
                "UO4",
                "VO4",
                "WO4",
                "UO5",
                "VO5",
                "WO5",
                "UO6",
                "VO6",
                "WO6",
                "WO",
                "WI1",
                "WO1",
            ),
        )
