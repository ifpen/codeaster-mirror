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

# ----------------------------------------------------------------------------------------------
# Located components
# ----------------------------------------------------------------------------------------------
NDEPLAC = LocatedComponents(
    phys=PHY.DEPL_C, type="ELNO", components=("DX", "DY", "DZ", "DRX", "DRY", "DRZ")
)

DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        ("EN1", ("DZ", "DRX", "DRY")),
        ("EN2", ("DRX", "DRY", "PHI", "LAG_GV")),
        ("EN3", ("DRX", "DRY")),
    ),
)

CCACOQU = LocatedComponents(
    phys=PHY.CACOQU_R, type="ELEM", components=("EP", "ALPHA", "BETA", "CTOR", "EXCENT", "INERTIE")
)

CCAORIE = LocatedComponents(
    phys=PHY.CAORIE_R,
    type="ELEM",
    components=("ALPHA", "BETA", "REP", "AXE_X", "AXE_Y", "AXE_Z", "O_X", "O_Y", "O_Z"),
)

ECHGREP = LocatedComponents(phys=PHY.CHGREPER, type="ELEM", components=("NATCHG", "CMAT[9]"))

NACCELR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))

EENERR = LocatedComponents(
    phys=PHY.ENER_R,
    type="ELEM",
    components=("TOTALE", "MEMBRANE", "FLEXION", "CISAILLE", "COUPL_MF"),
)

EENERPG = LocatedComponents(
    phys=PHY.ENER_R,
    type="ELGA",
    location="RIGI",
    components=("TOTALE", "MEMBRANE", "FLEXION", "CISAILLE", "COUPL_MF"),
)

EENERNO = LocatedComponents(
    phys=PHY.ENER_R,
    type="ELNO",
    components=("TOTALE", "MEMBRANE", "FLEXION", "CISAILLE", "COUPL_MF"),
)

EDEFOPC = LocatedComponents(
    phys=PHY.EPSI_C,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"),
)

EDEFGPC = LocatedComponents(
    phys=PHY.EPSI_C,
    type="ELGA",
    location="RIGI",
    components=("EXX", "EYY", "EXY", "KXX", "KYY", "KXY", "GAX", "GAY"),
)

CEPSINR = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EXX", "EYY", "EXY", "KXX", "KYY", "KXY"),
)

CEPSINF = LocatedComponents(
    phys=PHY.EPSI_F, type="ELEM", components=("EXX", "EYY", "EXY", "KXX", "KYY", "KXY")
)

EDEFOPG = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"),
)

EDEFONO = LocatedComponents(
    phys=PHY.EPSI_R, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)

EDFVCPG = LocatedComponents(phys=PHY.EPSI_R, type="ELGA", location="RIGI", components=("EPTHER_L",))

EDFVCNO = LocatedComponents(phys=PHY.EPSI_R, type="ELNO", components=("EPTHER_L",))

EDEFGNO = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELNO",
    components=("EXX", "EYY", "EXY", "KXX", "KYY", "KXY", "GAX", "GAY"),
)

EDEFGPG = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EXX", "EYY", "EXY", "KXX", "KYY", "KXY", "GAX", "GAY"),
)

CFORCEF = LocatedComponents(
    phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ", "MX", "MY", "MZ", "REP", "PLAN")
)

EFORCNO = LocatedComponents(
    phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ", "MX", "MY", "MZ", "REP", "PLAN")
)

ENBSP_I = LocatedComponents(phys=PHY.NBSP_I, type="ELEM", components=("COQ_NCOU",))

ENEU1_R = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X1",))

ELNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELEM", components=("X[30]",))

EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))

EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))

EMNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="MATER", components=("X1",))

ECASECT = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[16]",))

CPRESSF = LocatedComponents(phys=PHY.PRES_F, type="ELEM", components=("PRES",))

CPRES_R = LocatedComponents(phys=PHY.PRES_R, type="ELEM", components=("PRES",))

EPRESNO = LocatedComponents(phys=PHY.PRES_R, type="ELNO", components=("PRES",))

EEFGENOC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("NXX", "NYY", "NXY", "MXX", "MYY", "MXY", "QX", "QY")
)

EEFGEPGR = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=("NXX", "NYY", "NXY", "MXX", "MYY", "MXY", "QX", "QY"),
)

EEFGENOR = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("NXX", "NYY", "NXY", "MXX", "MYY", "MXY", "QX", "QY")
)

EGAMIMA = LocatedComponents(
    phys=PHY.SPMX_R,
    type="ELGA",
    location="RIGI",
    components=("VAL", "NUCOU", "NUSECT", "NUFIBR", "POSIC", "POSIS"),
)

ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))

MVECTAR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=NACCELR)

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


class PLAQ_MITC(Element):
    """Mechanics - Plate (Reissner -Mindlin Mixed Interpolation of Tensorial Components)"""

    meshType = MT.QUAD9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4)),
        SetOfNodes("EN2", (5, 6, 7, 8)),
        SetOfNodes("EN3", (9,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.QU9,
            gauss=(
                "RIGI=FPG9",
                "MASS=FPG9",
                "FPG1=FPG1",
                "NOEU_S=NOEU_S",
                "NOEU=NOEU",
                "MTGA=FPG9",
            ),
            mater=("RIGI", "NOEU", "FPG1", "MTGA"),
        ),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )
    calculs = (
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PERREUR, LC.CERROR),
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),
                (OP.TOU_INI_ELEM.PNEUT_F, LC.CNTINIF),
                (SP.PNEU1_R, LC.CNTINIR),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
            ),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PDEPL_R, LC.EGDEP2D),
                (OP.TOU_INI_ELGA.PDOMMAG, LC.EDOMGGA),
                (OP.TOU_INI_ELGA.PEPSI_R, LC.EGPS2DR),
                (OP.TOU_INI_ELGA.PGEOM_R, LC.EGGEO2D),
                (OP.TOU_INI_ELGA.PINST_R, LC.EGINST_R),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (OP.TOU_INI_ELGA.PSIEF_R, LC.EGIG2DR),
                (OP.TOU_INI_ELGA.PSOUR_R, LC.ESOURCR),
                (OP.TOU_INI_ELGA.PVARI_R, LC.ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PDOMMAG, LC.EDOMGNO),
                (OP.TOU_INI_ELNO.PEPSI_R, LC.EEPS2DR),
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM2D),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PSIEF_R, LC.ESIG2DR),
                (OP.TOU_INI_ELNO.PVARI_R, LC.ZVARINO),
            ),
        ),
        OP.SIEF_ELGA(
            te=33,
            para_in=(
                (SP.PCACOQU, CCACOQU),
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
        OP.VERI_PLAN(
            te=51,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PCHCKPR, LC.CCHCKPR)),
            para_out=((OP.VERI_PLAN.PCODRET, LC.ECODRET), (OP.VERI_PLAN.PINDICR, LC.CINDICR)),
        ),
        OP.RIGI_MECA(
            te=28,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PNBSP_I, ENBSP_I),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=29,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PMATERC, LC.CMATERC),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPRESSR, EPRESNO),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_F(
            te=29,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PMATERC, LC.CMATERC),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPRESSF, CPRESSF),
                (SP.PINSTR, LC.MTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.VERI_CARA_ELEM(
            te=119,
            para_in=((SP.PCACOQU, CCACOQU),),
            para_out=((SP.PCODRET, LC.ECODRET), (SP.PINDICR, LC.CINDICR)),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2), (OP.NSPG_NBVA.PNBSP_I, ENBSP_I)),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
    )
