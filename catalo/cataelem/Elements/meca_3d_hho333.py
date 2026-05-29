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


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        ("EN1", ("HHO_FX[10]", "HHO_FY[10]", "HHO_FZ[10]")),
        ("EN2", ()),
        ("EN3", ("HHO_CX[20]", "HHO_CY[20]", "HHO_CZ[20]")),
    ),
)

EDEPLPG = LocatedComponents(
    phys=PHY.DEPL_R, type="ELGA", location="RIGI", components=("DX", "DY", "DZ")
)


EENERR = LocatedComponents(phys=PHY.ENER_R, type="ELEM", components=("TOTALE",))


EENERPG = LocatedComponents(phys=PHY.ENER_R, type="ELGA", location="RIGI", components=("TOTALE",))


EENERNO = LocatedComponents(phys=PHY.ENER_R, type="ELNO", components=("TOTALE",))


EDEFOPC = LocatedComponents(
    phys=PHY.EPSI_C,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"),
)


EDEFONC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)


CEPSINF = LocatedComponents(
    phys=PHY.EPSI_F, type="ELEM", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
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


CEPSINR = LocatedComponents(
    phys=PHY.EPSI_R, type="ELEM", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
)


EDFVCPG = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EPTHER_L", "EPTHER_T", "EPTHER_N", "EPSECH", "EPHYDR", "EPPTOT"),
)


EDFVCNO = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELNO",
    components=("EPTHER_L", "EPTHER_T", "EPTHER_N", "EPSECH", "EPHYDR", "EPPTOT"),
)


EERREUR = LocatedComponents(
    phys=PHY.ERRE_R,
    type="ELEM",
    components=(
        "ERREST",
        "NUEST",
        "SIGCAL",
        "TERMRE",
        "TERMR2",
        "TERMNO",
        "TERMN2",
        "TERMSA",
        "TERMS2",
        "TAILLE",
    ),
)


EERRENO = LocatedComponents(
    phys=PHY.ERRE_R,
    type="ELNO",
    components=(
        "ERREST",
        "NUEST",
        "SIGCAL",
        "TERMRE",
        "TERMR2",
        "TERMNO",
        "TERMN2",
        "TERMSA",
        "TERMS2",
        "TAILLE",
    ),
)


EFACY_R = LocatedComponents(
    phys=PHY.FACY_R,
    type="ELGA",
    location="RIGI",
    components=(
        "DTAUM1",
        "VNM1X",
        "VNM1Y",
        "VNM1Z",
        "SINMAX1",
        "SINMOY1",
        "EPNMAX1",
        "EPNMOY1",
        "SIGEQ1",
        "NBRUP1",
        "ENDO1",
        "DTAUM2",
        "VNM2X",
        "VNM2Y",
        "VNM2Z",
        "SINMAX2",
        "SINMOY2",
        "EPNMAX2",
        "EPNMOY2",
        "SIGEQ2",
        "NBRUP2",
        "ENDO2",
        "VMIS",
        "TRESCA",
    ),
)


CFORCEF = LocatedComponents(phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ"))


NFORCER = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ"))


EFORCER = LocatedComponents(
    phys=PHY.FORC_R, type="ELGA", location="RIGI", components=("FX", "FY", "FZ")
)

NREACR = LocatedComponents(phys=PHY.REAC_R, type="ELNO", components=("DX", "DY", "DZ"))

NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGGEOM_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z")
)


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


ENGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


EMNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELEM", components=("X[30]",))


ECOPILO = LocatedComponents(
    phys=PHY.PILO_R, type="ELGA", location="RIGI", components=("A0", "A[3]", "ETA")
)


ECONTNC = LocatedComponents(
    phys=PHY.SIEF_C, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
)


ECONTPC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
)


ECONTNO = LocatedComponents(
    phys=PHY.SIEF_R, type="ELNO", components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
)


ECONTPG = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"),
)


ECOEQPG = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=(
        "VMIS",
        "TRESCA",
        "PRIN_[3]",
        "VMIS_SG",
        "VECT_1_X",
        "VECT_1_Y",
        "VECT_1_Z",
        "VECT_2_X",
        "VECT_2_Y",
        "VECT_2_Z",
        "VECT_3_X",
        "VECT_3_Y",
        "VECT_3_Z",
        "TRSIG",
        "TRIAX",
    ),
)


ESOURCR = LocatedComponents(phys=PHY.SOUR_R, type="ELGA", location="RIGI", components=("SOUR",))


ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))

# too large to be saved so used size 1
CHHOGT = LocatedComponents(phys=PHY.N1920R, type="ELEM", components=("X[1]",))

CHHOST = LocatedComponents(phys=PHY.N1360R, type="ELEM", components=("X[1]",))

CHHOBS = LocatedComponents(
    phys=PHY.N3600R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("X[55]",)), ("EN2", ()), ("EN3", ("X[630]"))),
)

CPPOST = LocatedComponents(phys=PHY.NEUT_I, type="ELEM", components=("X[1]"))
DEPLHHO = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))

PFONC = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z[21]",))

PFONCR = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z[4]",))

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)

# ------------------------------------------------------------


class MECA3DH27_HHO333(Element):
    """Please document this element"""

    meshType = MT.HEXA27
    nodes = (
        SetOfNodes("EN1", (21, 22, 23, 24, 25, 26)),
        SetOfNodes("EN2", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN3", (27,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.H27,
            gauss=("RIGI=FPG64", "FPG1=FPG1", "MTGA=FPG64", "MASS=FPG64"),
            mater=("RIGI", "FPG1", "MTGA", "MASS"),
        ),
    )
    calculs = (
        OP.AMOR_MECA(
            te=121,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMASSEL, MMATUUR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PRIGIEL, MMATUUR),
                (OP.AMOR_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.CHAR_MECA_FF3D3D(
            te=476,
            para_in=(
                (SP.PFF3D3D, CFORCEF),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTR, CTEMPSR),
                (OP.CHAR_MECA_FF3D3D.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FR3D3D(
            te=476,
            para_in=(
                (SP.PFR3D3D, NFORCER),
                (SP.PGEOMER, NGEOMER),
                (OP.CHAR_MECA_FR3D3D.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
        ),
        OP.DEPL_ELGA(
            te=456,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PDEPLAR, DDL_MECA), (OP.DEPL_ELGA.PCHHOBS, CHHOBS)),
            para_out=((SP.PDEPLGA, LC.EGDEP3D),),
        ),
        OP.EPSI_ELGA(
            te=448,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (OP.EPSI_ELGA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PDEFOPC, EDEFOPC), (OP.EPSI_ELGA.PDEFOPG, EDEFOPG)),
        ),
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, EDEFOPG),),
            para_out=((SP.PDEFONC, EDEFONC), (SP.PDEFONO, EDEFONO)),
        ),
        OP.EPVC_ELGA(
            te=529,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPVC_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPVC_ELGA.PDEFOPG, EDFVCPG),),
        ),
        OP.EPVC_ELNO(
            te=4, para_in=((OP.EPVC_ELNO.PDEFOPG, EDFVCPG),), para_out=((SP.PDEFONO, EDFVCNO),)
        ),
        OP.FORC_NODA(
            te=450,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, ECONTPG),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCPR, LC.ZVARCPG),
                (OP.FORC_NODA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=455,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, ZVARIPG),
                (OP.FULL_MECA.PVARIMR, ZVARIPG),
                (OP.FULL_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, ECONTPG),
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.HHO_PRECALC_BS(
            te=494,
            para_in=((SP.PGEOMER, NGEOMER),),
            para_out=((OP.HHO_PRECALC_BS.PCHHOBO, CHHOBS),),
        ),
        OP.HHO_PRECALC_OP(
            te=460,
            para_in=((SP.PGEOMER, NGEOMER), (OP.HHO_PRECALC_OP.PCHHOBS, CHHOBS)),
            para_out=((OP.HHO_PRECALC_OP.PCHHOGT, CHHOGT), (OP.HHO_PRECALC_OP.PCHHOST, CHHOST)),
        ),
        OP.HHO_PROJ_MECA(
            te=473,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (OP.HHO_PROJ_MECA.PFUNC_R, PFONCR),
                (SP.PINSTPR, CTEMPSR),
                (OP.HHO_PROJ_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_PROJ_MECA.PDEPL_R, DDL_MECA),),
        ),
        OP.HHO_PROJ3_MECA(
            te=484,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (OP.HHO_PROJ3_MECA.PQPTP_R, LC.EGDEP3D),
                (OP.HHO_PROJ3_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_PROJ3_MECA.PDEPL_R, DDL_MECA),),
        ),
        OP.HHO_DEPL_MECA(
            te=456,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PDEPLPR, DDL_MECA),
                (OP.HHO_DEPL_MECA.PCHHOBS, CHHOBS),
                (SP.POPPOST, CPPOST),
            ),
            para_out=((OP.HHO_DEPL_MECA.PDEPL_R, DEPLHHO),),
        ),
        OP.HHO_CINE_F_MECA(
            te=458,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTPR, CTEMPSR),
                (OP.HHO_CINE_F_MECA.PFONC, PFONC),
                (OP.HHO_CINE_F_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_CINE_F_MECA.PCINE, DDL_MECA),),
        ),
        OP.HHO_CINE_R_MECA(
            te=458,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (OP.HHO_CINE_R_MECA.PCMPVALE, DEPLHHO),
                (OP.HHO_CINE_R_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_CINE_R_MECA.PCINE, DDL_MECA),),
        ),
        OP.INIT_VARC(
            te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG), (OP.INIT_VARC.PVARCNO, LC.ZVARCNO))
        ),
        OP.MASS_MECA(
            te=474,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
                (OP.MASS_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MATE_ELGA(
            te=142,
            para_in=(
                (SP.PMATERC, LC.CMATERC),
                (SP.PGEOMER, NGEOMER),
                (OP.MATE_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.MATE_ELGA.PMATERR, LC.EGMATE_R),),
        ),
        OP.MATE_ELEM(
            te=142,
            para_in=(
                (SP.PMATERC, LC.CMATERC),
                (SP.PGEOMER, NGEOMER),
                (OP.MATE_ELEM.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.MATE_ELEM.PMATERR, LC.EEMATE_R),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.RAPH_MECA(
            te=455,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARIMR, ZVARIPG),
                (OP.RAPH_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, ECONTPG),
                (OP.RAPH_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.REFE_FORC_NODA(
            te=450,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PREFCO, LC.CRESSIG),
                (OP.REFE_FORC_NODA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_MECA(
            te=455,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PINSTR, CTEMPSR),
                (OP.RIGI_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_ELAS(
            te=455,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (OP.RIGI_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_ELAS.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RIGI_MECA_ELAS.PVARIMR, ZVARIPG),
                (OP.RIGI_MECA_ELAS.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.RIGI_MECA_TANG(
            te=455,
            para_in=(
                (SP.PCAMASS, LC.CCAMA3D),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, ECONTPG),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, NGEOMER),
                (SP.PINSTMR, CTEMPSR),
                (SP.PINSTPR, CTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARIMR, ZVARIPG),
                (OP.RIGI_MECA_TANG.PCHHOBS, CHHOBS),
            ),
            para_out=(
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, ECONTPG),
                (SP.PCODRET, LC.ECODRET),
                (SP.PCOPRED, LC.ECODRET),
            ),
        ),
        OP.SIEF_ELNO(
            te=4,
            para_in=((OP.SIEF_ELNO.PCONTRR, ECONTPG), (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PSIEFNOC, ECONTNC), (OP.SIEF_ELNO.PSIEFNOR, ECONTNO)),
        ),
        OP.SIEQ_ELGA(
            te=335,
            para_in=((OP.SIEQ_ELGA.PCONTRR, ECONTPG),),
            para_out=((OP.SIEQ_ELGA.PCONTEQ, ECOEQPG),),
        ),
        OP.SIEQ_ELNO(
            te=335,
            para_in=((OP.SIEQ_ELNO.PCONTRR, ECONTNO),),
            para_out=((OP.SIEQ_ELNO.PCONTEQ, LC.ECOEQNO),),
        ),
        OP.SIGM_ELGA(
            te=546,
            para_in=((SP.PSIEFR, ECONTPG),),
            para_out=((SP.PSIGMC, ECONTPC), (SP.PSIGMR, ECONTPG)),
        ),
        OP.SIGM_ELNO(
            te=4,
            para_in=((OP.SIGM_ELNO.PCONTRR, ECONTPG),),
            para_out=((SP.PSIEFNOC, ECONTNC), (OP.SIGM_ELNO.PSIEFNOR, ECONTNO)),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D), (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER)),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PDEPL_R, EDEPLPG),
                (OP.TOU_INI_ELGA.PDOMMAG, LC.EDOMGGA),
                (OP.TOU_INI_ELGA.PEPSI_R, EDEFOPG),
                (SP.PFACY_R, EFACY_R),
                (OP.TOU_INI_ELGA.PGEOM_R, EGGEOM_R),
                (OP.TOU_INI_ELGA.PINST_R, LC.EGINST_R),
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
                (OP.TOU_INI_ELGA.PSIEF_R, ECONTPG),
                (OP.TOU_INI_ELGA.PSOUR_R, ESOURCR),
                (OP.TOU_INI_ELGA.PVARI_R, ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PDOMMAG, LC.EDOMGNO),
                (OP.TOU_INI_ELNO.PEPSI_R, EDEFONO),
                (OP.TOU_INI_ELNO.PGEOM_R, ENGEOM_R),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PSIEF_R, ECONTNO),
                (OP.TOU_INI_ELNO.PVARI_R, LC.ZVARINO),
            ),
        ),
        OP.VARC_ELGA(
            te=530,
            para_in=((OP.VARC_ELGA.PVARCPR, LC.ZVARCPG),),
            para_out=((SP.PVARC_R, LC.EVARC_R),),
        ),
        OP.VARI_ELNO(
            te=4, para_in=((SP.PVARIGR, ZVARIPG),), para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),)
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# ------------------------------------------------------------


class MECA3DT15_HHO333(MECA3DH27_HHO333):
    """Please document this element"""

    meshType = MT.TETRA15
    nodes = (
        SetOfNodes("EN1", (11, 12, 13, 14)),
        SetOfNodes("EN2", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
        SetOfNodes("EN3", (15,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.T15,
            gauss=("RIGI=FPG24", "FPG1=FPG1", "MTGA=FPG24", "MASS=FPG24"),
            mater=("RIGI", "FPG1", "MTGA", "MASS"),
        ),
    )


class MECA3DP21_HHO333(MECA3DH27_HHO333):
    """Maille pentaèdre / prisme à 21 noeuds pour HHO en quadratique"""

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


class MECA3DP19_HHO333(MECA3DH27_HHO333):
    """Maille pyramide à 19 noeuds pour HHO en quadratique"""

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
