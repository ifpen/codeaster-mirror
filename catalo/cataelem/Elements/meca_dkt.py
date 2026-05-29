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
    phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ", "DRX", "DRY", "DRZ")
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

ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))

MVECTAR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=NACCELR)

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MEDKQU4(Element):
    """Mechanics - Plate (DKT) - QUAD4"""

    meshType = MT.QUAD4
    elrefe = (
        ElrefeLoc(
            MT.QU4,
            gauss=("RIGI=FPG4", "MASS=FPG9", "NOEU=NOEU", "FPG1=FPG1"),
            mater=("RIGI", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )
    calculs = (
        OP.ACCEPTANCE(
            te=329,
            para_in=((SP.PACCELR, NACCELR), (SP.PGEOMER, LC.EGEOM3D), (SP.PNUMMOD, LC.CNUMMOD)),
            para_out=((SP.PVECTUR, MVECTAR),),
        ),
        OP.ADD_SIGM(
            te=581,
            para_in=((SP.PEPCON1, LC.EGIG3DR), (SP.PEPCON2, LC.EGIG3DR)),
            para_out=((SP.PEPCON3, LC.EGIG3DR),),
        ),
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
        OP.CARA_SECT_POUT3(
            te=513, para_in=((SP.PGEOMER, LC.EGEOM3D),), para_out=((SP.PCASECT, ECASECT),)
        ),
        OP.CARA_SECT_POUT4(
            te=513,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PORIGIN, LC.CGEOM3D)),
            para_out=((SP.PVECTU1, MVECTUR), (SP.PVECTU2, MVECTUR)),
        ),
        OP.CHAR_MECA_EPSA_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_EPSA_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_EPSI_R(
            te=35,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PEPSINR, CEPSINR),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_EPSI_R.PNBSP_I, ENBSP_I),
                (OP.CHAR_MECA_EPSI_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_EPSI_F(
            te=35,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PEPSINF, CEPSINF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_EPSI_F.PNBSP_I, ENBSP_I),
                (OP.CHAR_MECA_EPSI_F.PVARCPR, LC.ZVARCPG),
                (SP.PINSTR, LC.MTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FFCO3D(
            te=32,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PFFCO3D, CFORCEF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTR, LC.MTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FRCO3D(
            te=32,
            para_in=((SP.PCACOQU, CCACOQU), (SP.PFRCO3D, EFORCNO), (SP.PGEOMER, LC.EGEOM3D)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_HYDR_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_HYDR_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=32,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_F(
            te=32,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPRESSF, CPRESSF),
                (SP.PINSTR, LC.MTEMPSR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=32,
            para_in=((SP.PCACOQU, CCACOQU), (SP.PGEOMER, LC.EGEOM3D), (SP.PPRESSR, EPRESNO)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRSU_F(
            te=580, para_in=((SP.PPRESSF, CPRESSF),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.CHAR_MECA_PRSU_R(
            te=580, para_in=((SP.PPRESSR, EPRESNO),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.CHAR_MECA_SECH_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_SECH_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SFCO3D(
            te=580, para_in=((SP.PFFCO3D, CFORCEF),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.CHAR_MECA_SRCO3D(
            te=580, para_in=((SP.PFRCO3D, EFORCNO),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.CHAR_MECA_TEMP_R(
            te=35,
            para_in=(
                (SP.PCACOQU, CCACOQU),
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
            te=488,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.COOR_ELGA.PNBSP_I, ENBSP_I),
            ),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU3D), (OP.COOR_ELGA.PCOORSU, LC.EGGAU3D)),
        ),
        OP.COOR_ELGA_MATER(
            te=462,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.COOR_ELGA_MATER.PNBSP_I, ENBSP_I),
            ),
            para_out=((SP.PCOOPGM, LC.EGGEMA_R),),
        ),
        OP.DEGE_ELGA(
            te=33,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.DEGE_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.DEGE_ELGA.PDEFOPG, EDEFGPG),),
        ),
        OP.DEGE_ELNO(
            te=33,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.DEGE_ELNO.PNBSP_I, ENBSP_I),
                (OP.DEGE_ELNO.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PDEFOGR, EDEFGNO),),
        ),
        OP.DEPL_ELGA(
            te=552,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.DEPL_ELGA.PNBSP_I, ENBSP_I),
                (SP.PDEPLAR, DDL_MECA),
            ),
            para_out=((SP.PDEPLGA, LC.EGDEP3D),),
        ),
        OP.ECIN_ELEM(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.ECIN_ELEM.PNBSP_I, ENBSP_I),
                (SP.POMEGA2, LC.COMEG2R),
                (OP.ECIN_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVITESR, DDL_MECA),
            ),
            para_out=((SP.PENERCR, EENERR),),
        ),
        OP.EFGE_ELGA(
            te=451,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PMATERC, LC.CMATERC),
                (OP.EFGE_ELGA.PNBSP_I, ENBSP_I),
                (SP.PSIEFR, LC.EGIG3DR),
            ),
            para_out=((SP.PEFGEC, LC.EGFG3DC), (SP.PEFGER, EEFGEPGR)),
        ),
        OP.EFGE_ELNO(
            te=185,
            para_in=(
                (SP.PCACOQU, CCACOQU),
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
            para_out=((SP.PEFFORC, EEFGENOC), (OP.EFGE_ELNO.PEFFORR, EEFGENOR)),
        ),
        OP.EFGE_EXCENT(
            te=452,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PEFFOGC, LC.EGFG3DC),
                (SP.PEFFOGR, EEFGEPGR),
                (SP.PEFFONC, EEFGENOC),
                (SP.PEFFONR, EEFGENOR),
            ),
            para_out=(
                (SP.PEFFOEGC, LC.EGFG3DC),
                (SP.PEFFOEGR, EEFGEPGR),
                (SP.PEFFOENC, EEFGENOC),
                (SP.PEFFOENR, EEFGENOR),
            ),
        ),
        OP.ENEL_ELEM(
            te=412,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.ENEL_ELEM.PCOMPOR, LC.CCOMPOR),
                (OP.ENEL_ELEM.PCONTPR, LC.EGIG3DR),
                (SP.PDEPLR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.ENEL_ELEM.PNBSP_I, ENBSP_I),
                (OP.ENEL_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.ENEL_ELEM.PVARIPR, ZVARIPG),
            ),
            para_out=((SP.PENERD1, EENERR),),
        ),
        OP.ENEL_ELGA(
            te=412,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.ENEL_ELGA.PCOMPOR, LC.CCOMPOR),
                (OP.ENEL_ELGA.PCONTRR, LC.EGIG3DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.ENEL_ELGA.PNBSP_I, ENBSP_I),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.ENEL_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARIGR, ZVARIPG),
            ),
            para_out=((OP.ENEL_ELGA.PENERDR, EENERPG),),
        ),
        OP.ENEL_ELNO(
            te=4, para_in=((OP.ENEL_ELNO.PENERPG, EENERPG),), para_out=((SP.PENERNO, EENERNO),)
        ),
        OP.EPEQ_ELGA(
            te=335,
            para_in=((OP.EPEQ_ELGA.PDEFORR, EDEFOPG),),
            para_out=((OP.EPEQ_ELGA.PDEFOEQ, LC.EDFEQPG),),
        ),
        OP.EPEQ_ELNO(
            te=335,
            para_in=((OP.EPEQ_ELNO.PDEFORR, EDEFONO),),
            para_out=((OP.EPEQ_ELNO.PDEFOEQ, LC.EDFEQNO),),
        ),
        OP.EPME_ELGA(
            te=531,
            para_in=(
                (OP.EPME_ELGA.PDEFORR, EDEFOPG),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPME_ELGA.PNBSP_I, ENBSP_I),
                (OP.EPME_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPME_ELGA.PDEFOPG, EDEFOPG),),
        ),
        OP.EPME_ELNO(
            te=4, para_in=((OP.EPME_ELNO.PDEFOPG, EDEFOPG),), para_out=((SP.PDEFONO, EDEFONO),)
        ),
        OP.EPOT_ELEM(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPOT_ELEM.PNBSP_I, ENBSP_I),
                (OP.EPOT_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPOT_ELEM.PENERDR, EENERR),),
        ),
        OP.CHAR_MECA_EPSA_R(
            te=312,
            para_in=((SP.PMATERC, LC.CMATERC), (OP.CHAR_MECA_EPSA_R.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.EPSI_ELGA(
            te=33,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPSI_ELGA.PNBSP_I, ENBSP_I),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EPSI_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PDEFOPC, EDEFOPC), (OP.EPSI_ELGA.PDEFOPG, EDEFOPG)),
        ),
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, EDEFOPG),),
            para_out=((SP.PDEFONC, LC.EEPS3DC), (SP.PDEFONO, EDEFONO)),
        ),
        OP.EPSP_ELGA(
            te=531,
            para_in=(
                (OP.EPSP_ELGA.PCONTRR, LC.EGIG3DR),
                (OP.EPSP_ELGA.PDEFORR, EDEFOPG),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPSP_ELGA.PNBSP_I, ENBSP_I),
                (OP.EPSP_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPSP_ELGA.PDEFOPG, EDEFOPG),),
        ),
        OP.EPSP_ELNO(
            te=4, para_in=((OP.EPSP_ELNO.PDEFOPG, EDEFOPG),), para_out=((SP.PDEFONO, EDEFONO),)
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
            para_out=((OP.EPVC_ELGA.PDEFOPG, EDFVCPG),),
        ),
        OP.EPVC_ELNO(
            te=4, para_in=((OP.EPVC_ELNO.PDEFOPG, EDFVCPG),), para_out=((SP.PDEFONO, EDFVCNO),)
        ),
        OP.FERR_ELEM(
            te=146,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.FERR_ELEM.PEFFORR, EEFGENOR),
                (SP.PFERRA1, LC.CFER1_R),
            ),
            para_out=((SP.PFERRA2, LC.CFER2_R),),
        ),
        OP.MARG_ELEM(
            te=53,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.MARG_ELEM.PEFFORR, EEFGENOR),
                (OP.MARG_ELEM.PEFFOR0, EEFGENOR),
                (SP.PVFER1, LC.CVFER1_R),
                (SP.PVFER0, LC.CVFER0_R),
            ),
            para_out=((SP.PVFER2, LC.CVFER2_R),),
        ),
        OP.FORC_NODA(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, LC.EGIG3DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.FORC_NODA.PNBSP_I, ENBSP_I),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, LC.EGIG3DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PMATERC, LC.CMATERC),
                (OP.FULL_MECA.PNBSP_I, ENBSP_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.FULL_MECA.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, LC.EGIG3DR),
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.FULL_MECA_ELAS(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.FULL_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA_ELAS.PCONTMR, LC.EGIG3DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (OP.FULL_MECA_ELAS.PNBSP_I, ENBSP_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.FULL_MECA_ELAS.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA_ELAS.PCONTPR, LC.EGIG3DR),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA_ELAS.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.INI_SP_MATER(
            te=99,
            para_out=((SP.PHYDMAT, LC.EHYDRMA), (SP.PNEUMAT, EMNEUT_R), (SP.PTEMMAT, LC.ETEMPMA)),
        ),
        OP.INI_SP_RIGI(te=99, para_out=((OP.INI_SP_RIGI.PCONTPR, LC.EGIG3DR),)),
        OP.MASS_INER(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_INER.PNBSP_I, ENBSP_I),
                (OP.MASS_INER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMASSINE, LC.EMASSINE),),
        ),
        OP.MASS_MECA(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PNBSP_I, ENBSP_I),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_MECA_EXPLI(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA_EXPLI.PNBSP_I, ENBSP_I),
                (OP.MASS_MECA_EXPLI.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MINMAX_SP(te=99, para_out=((SP.PGAMIMA, LC.EGMINMAX), (SP.PNOMIMA, LC.NMINMAX))),
        OP.MODI_REPERE(
            te=63,
            para_in=((OP.MODI_REPERE.PMATPASS, ECHGREP), (SP.PSIEFR, LC.EGIG3DR)),
            para_out=((OP.MODI_REPERE.PCONTPR, LC.EGIG3DR),),
        ),
        OP.M_GAMMA(
            te=31,
            para_in=(
                (SP.PACCELR, DDL_MECA),
                (SP.PCACOQU, CCACOQU),
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
                (SP.PCACOQU, CCACOQU),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (OP.PAS_COURANT.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PCOURAN, LC.ECOURAN),),
        ),
        OP.PREP_VRC(
            te=408,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.PREP_VRC.PINST_R, LC.MTEMPSR),
                (OP.PREP_VRC.PNBSP_I, ENBSP_I),
                (SP.PTEMPEF, LC.CTEMPEF),
                (SP.PTEMPER, LC.NTEMPER),
                (SP.PGEOMER, LC.EGEOM3D),
            ),
            para_out=((SP.PTEMPCR, LC.CTEREFE),),
        ),
        OP.RAPH_MECA(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, LC.EGIG3DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PMATERC, LC.CMATERC),
                (OP.RAPH_MECA.PNBSP_I, ENBSP_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.RAPH_MECA.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, LC.EGIG3DR),
                (OP.RAPH_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.REFE_FORC_NODA(
            te=31,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PREFCO, LC.CRESEFM)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.REPERE_LOCAL(
            te=134,
            para_in=((SP.PCACOQU, CCACOQU), (SP.PGEOMER, LC.EGEOM3D)),
            para_out=(
                (OP.REPERE_LOCAL.PMATPASS, ECHGREP),
                (SP.PREPLO1, LC.CGEOM3D),
                (SP.PREPLO2, LC.CGEOM3D),
                (SP.PREPLO3, LC.CGEOM3D),
            ),
        ),
        OP.REPE_GENE(
            te=442,
            para_in=(
                (SP.PANGREP, CCAORIE),
                (SP.PCACOQU, CCACOQU),
                (SP.PDGGAIN, EDEFGPG),
                (SP.PDGGAINC, EDEFGPC),
                (SP.PDGNOIN, EDEFGNO),
                (SP.PDGNOINC, LC.EDEFGNC),
                (SP.PEFGAIN, EEFGEPGR),
                (SP.PEFGAINC, LC.EGFG3DC),
                (SP.PEFNOIN, EEFGENOR),
                (SP.PEFNOINC, EEFGENOC),
                (SP.PGEOMER, LC.EGEOM3D),
            ),
            para_out=(
                (SP.PDGGAOUC, EDEFGPC),
                (SP.PDGGAOUT, EDEFGPG),
                (SP.PDGNOOUC, LC.EDEFGNC),
                (SP.PDGNOOUT, EDEFGNO),
                (SP.PEFGAOUC, LC.EGFG3DC),
                (SP.PEFGAOUT, EEFGEPGR),
                (SP.PEFNOOUC, EEFGENOC),
                (SP.PEFNOOUT, EEFGENOR),
            ),
        ),
        OP.REPE_TENS(
            te=442,
            para_in=(
                (SP.PANGREP, CCAORIE),
                (SP.PCACOQU, CCACOQU),
                (SP.PCOGAIN, LC.EGIG3DR),
                (SP.PCONOIN, LC.ESIG3DR),
                (SP.PDEGAIN, EDEFOPG),
                (SP.PDENOIN, EDEFONO),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.REPE_TENS.PNBSP_I, ENBSP_I),
            ),
            para_out=(
                (SP.PCOGAOUT, LC.EGIG3DR),
                (SP.PCONOOUT, LC.ESIG3DR),
                (SP.PDEGAOUT, EDEFOPG),
                (SP.PDENOOUT, EDEFONO),
            ),
        ),
        OP.RIGI_MECA(
            te=31,
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
        OP.RIGI_MECA_ELAS(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_ELAS.PCONTMR, LC.EGIG3DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA_ELAS.PNBSP_I, ENBSP_I),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.RIGI_MECA_ELAS.PVARIMR, ZVARIPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_GEOM(
            te=428,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.RIGI_GEOM.PCONTRR, LC.EGIG3DR),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.RIGI_GEOM.PNBSP_I, ENBSP_I),
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
        #       -- te0580 : ne resout que les cas triviaux : 0.
        OP.RIGI_MECA_PRSU_F(
            te=580, para_in=((SP.PPRESSF, CPRESSF),), para_out=((SP.PMATUNS, MMATUNS),)
        ),
        OP.RIGI_MECA_PRSU_R(
            te=580, para_in=((SP.PPRESSR, EPRESNO),), para_out=((SP.PMATUNS, MMATUNS),)
        ),
        OP.RIGI_MECA_SFCO3D(
            te=580, para_in=((SP.PFFCO3D, CFORCEF),), para_out=((SP.PMATUNS, MMATUNS),)
        ),
        OP.RIGI_MECA_SRCO3D(
            te=580, para_in=((SP.PFRCO3D, EFORCNO),), para_out=((SP.PMATUNS, MMATUNS),)
        ),
        OP.RIGI_MECA_TANG(
            te=31,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (SP.PCARCRI, LC.CCARCRI),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, LC.EGIG3DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PNBSP_I, ENBSP_I),
                (OP.RIGI_MECA_TANG.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, LC.EGIG3DR),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
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
        OP.TEMP_ELGA(
            te=126,
            para_in=(
                (SP.PCACOQU, CCACOQU),
                (OP.TEMP_ELGA.PNBSP_I, ENBSP_I),
                (OP.TEMP_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PTEMP_R, LC.ETEMPPG),),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PNEUT_F, ELNEUT_F),
                (SP.PNEU1_R, ENEU1_R),
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),
                (OP.TOU_INI_ELEM.PNBSP_I, ENBSP_I),
                (OP.TOU_INI_ELEM.PPRES_R, CPRES_R),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
            ),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
                (OP.TOU_INI_ELGA.PSIEF_R, LC.EGIG3DR),
                (OP.TOU_INI_ELGA.PVARI_R, ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM3D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (SP.PTEMPN_R, LC.ETEMPNO),
                (OP.TOU_INI_ELNO.PPRES_R, EPRESNO),
                (OP.TOU_INI_ELNO.PSIEF_R, EEFGENOR),
            ),
        ),
        OP.VARC_ELGA(
            te=530,
            para_in=((OP.VARC_ELGA.PVARCPR, LC.ZVARCPG),),
            para_out=((SP.PVARC_R, LC.EVARC_R),),
        ),
        OP.VARC_ELNO(
            te=4, para_in=((SP.PVARCGR, LC.EVARC_R),), para_out=((SP.PVARCNR, LC.EVARCNR),)
        ),
        OP.VARI_ELNO(
            te=4, para_in=((SP.PVARIGR, ZVARIPG),), para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),)
        ),
        OP.VERI_CARA_ELEM(
            te=119,
            para_in=((SP.PCACOQU, CCACOQU),),
            para_out=((SP.PCODRET, LC.ECODRET), (SP.PINDICR, LC.CINDICR)),
        ),
        OP.VERI_PLAN(
            te=51,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PCHCKPR, LC.CCHCKPR)),
            para_out=((SP.PCODRET, LC.ECODRET), (SP.PINDICR, LC.CINDICR)),
        ),
    )


# ------------------------------------------------------------
class MEDKTR3(MEDKQU4):
    """Mechanics - Plate (DKT) - TRIA3"""

    meshType = MT.TRIA3
    elrefe = (
        ElrefeLoc(
            MT.TR3,
            gauss=("RIGI=FPG3", "MASS=FPG4", "NOEU=NOEU", "FPG1=FPG1"),
            mater=("RIGI", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )

    calculs = (OP.ACCEPTANCE(te=-1),)
