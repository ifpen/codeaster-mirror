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
CCAMASS = LocatedComponents(
    phys=PHY.CAMA_R, type="ELEM", components=("C", "ALPHA", "BETA", "KAPPA", "X", "Y", "Z")
)


CTEMPSR = LocatedComponents(phys=PHY.INST_R, type="ELEM", components=("INST",))

NDEPLAR = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))

DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("DX", "DY", "DZ", "PRES", "GONF")), ("EN2", ("DX", "DY", "DZ"))),
)


EDEFOPC = LocatedComponents(
    phys=PHY.EPSI_C,
    type="ELGA",
    location="RIGI",
    components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ"),
)


EDEFONC = LocatedComponents(
    phys=PHY.EPSI_C, type="ELNO", components=("EPXX", "EPYY", "EPZZ", "EPXY", "EPXZ", "EPYZ")
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


NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))

NFORCER = LocatedComponents(phys=PHY.FORC_R, type="ELNO", components=("FX", "FY", "FZ"))

EGGEOM_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z")
)


EGGEOP_R = LocatedComponents(
    phys=PHY.GEOM_R, type="ELGA", location="RIGI", components=("X", "Y", "Z", "W")
)


ENGEOM_R = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))


EGNEUT_F = LocatedComponents(phys=PHY.NEUT_F, type="ELGA", location="RIGI", components=("X[30]",))


EGNEUT_R = LocatedComponents(phys=PHY.NEUT_R, type="ELGA", location="RIGI", components=("X[30]",))


EREFCO = LocatedComponents(phys=PHY.PREC_R, type="ELEM", components=("SIGM", "EPSI", "LAG_GV"))


ECONTNC = LocatedComponents(
    phys=PHY.SIEF_C,
    type="ELNO",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ", "PRES_M", "PRES", "SIG[3]", "DEPV"),
)


ECONTPG = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELGA",
    location="RIGI",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ", "PRES_M", "PRES", "SIG[3]", "DEPV"),
)


ECONTNO = LocatedComponents(
    phys=PHY.SIEF_R,
    type="ELNO",
    components=("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ", "PRES_M", "PRES", "SIG[3]", "DEPV"),
)

CFORCEF = LocatedComponents(phys=PHY.FORC_F, type="ELEM", components=("FX", "FY", "FZ"))

ZVARIPG = LocatedComponents(phys=PHY.VARI_R, type="ELGA", location="RIGI", components=("VARI",))


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MVECTDR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=NDEPLAR)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class T10_3D(Element):
    """Please document this element"""

    meshType = MT.TETRA10
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8, 9, 10)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(MT.T10, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.TE4, gauss=("RIGI=FPG4",)),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG3",)),
    )
    calculs = (
        OP.CHAR_MECA_FF3D3D(
            te=17,
            para_in=((SP.PFF3D3D, CFORCEF), (SP.PGEOMER, NGEOMER), (SP.PINSTR, CTEMPSR)),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_FR3D3D(
            te=16,
            para_in=((SP.PFR3D3D, NFORCER), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=15,
            para_in=(
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTDR),),
        ),
        OP.COOR_ELGA(
            te=488, para_in=((SP.PGEOMER, NGEOMER),), para_out=((OP.COOR_ELGA.PCOORPG, EGGEOP_R),)
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
        OP.EPSI_ELGA(
            te=25,
            para_in=(
                (SP.PCAMASS, CCAMASS),
                (SP.PDEPLAR, NDEPLAR),
                (SP.PGEOMER, NGEOMER),
                (SP.PMATERC, LC.CMATERC),
            ),
            para_out=((SP.PDEFOPC, EDEFOPC), (OP.EPSI_ELGA.PDEFOPG, EDEFOPG)),
        ),
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, EDEFOPG),),
            para_out=((SP.PDEFONC, EDEFONC), (SP.PDEFONO, EDEFONO)),
        ),
        OP.FORC_NODA(
            te=5,
            para_in=((SP.PCOMPOR, LC.CCOMPOR), (SP.PSIEFR, ECONTPG), (SP.PGEOMER, NGEOMER)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=5,
            para_in=(
                (SP.PCAMASS, CCAMASS),
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
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.RAPH_MECA(
            te=5,
            para_in=(
                (SP.PCAMASS, CCAMASS),
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
                (SP.PVARIMP, ZVARIPG),
                (OP.RAPH_MECA.PVARIMR, ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, ECONTPG),
                (OP.RAPH_MECA.PVARIPR, ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.REFE_FORC_NODA(
            te=5,
            para_in=(
                (OP.REFE_FORC_NODA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, NGEOMER),
                (SP.PREFCO, EREFCO),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_MECA_ELAS(
            te=5,
            para_in=(
                (SP.PCAMASS, CCAMASS),
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
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.RIGI_MECA_TANG(
            te=5,
            para_in=(
                (SP.PCAMASS, CCAMASS),
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
            ),
            para_out=(
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, ECONTPG),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.SIEF_ELNO(
            te=4,
            para_in=((OP.SIEF_ELNO.PCONTRR, ECONTPG), (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PSIEFNOC, ECONTNC), (OP.SIEF_ELNO.PSIEFNOR, ECONTNO)),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PGEOM_R, EGGEOM_R),
                (OP.TOU_INI_ELGA.PINST_R, LC.EGINST_R),
                (OP.TOU_INI_ELGA.PNEUT_F, EGNEUT_F),
                (OP.TOU_INI_ELGA.PNEUT_R, EGNEUT_R),
                (OP.TOU_INI_ELGA.PSIEF_R, ECONTPG),
                (OP.TOU_INI_ELGA.PVARI_R, ZVARIPG),
            ),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D), (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER)),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, ENGEOM_R),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
            ),
        ),
        OP.VARI_ELNO(
            te=4, para_in=((SP.PVARIGR, ZVARIPG),), para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),)
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, NGEOMER),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# ------------------------------------------------------------
class P15_3D(T10_3D):
    """Please document this element"""

    meshType = MT.PENTA15
    nodes = (
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
    )
    elrefe = (
        ElrefeLoc(MT.P15, gauss=("RIGI=FPG6", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.PE6, gauss=("RIGI=FPG6",)),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG4",)),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG3",)),
    )


# ------------------------------------------------------------
class H20_3D(T10_3D):
    """Please document this element"""

    meshType = MT.HEXA20
    nodes = (
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
    )
    elrefe = (
        ElrefeLoc(MT.H20, gauss=("RIGI=FPG8", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.HE8, gauss=("RIGI=FPG8",)),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG4",)),
    )


class T10_3DI(T10_3D):
    """Please document this element"""

    meshType = MT.TETRA10
    nodes = (SetOfNodes("EN2", (5, 6, 7, 8, 9, 10)), SetOfNodes("EN1", (1, 2, 3, 4)))
    elrefe = (
        ElrefeLoc(MT.T10, gauss=("RIGI=FPG4", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.TE4, gauss=("RIGI=FPG4",)),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG3",)),
    )


class P15_3DI(T10_3D):
    """Please document this element"""

    meshType = MT.PENTA15
    nodes = (
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
    )
    elrefe = (
        ElrefeLoc(MT.P15, gauss=("RIGI=FPG6", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.PE6, gauss=("RIGI=FPG6",)),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG4",)),
        ElrefeLoc(MT.TR6, gauss=("RIGI=FPG3",)),
    )


# ------------------------------------------------------------
class H20_3DI(T10_3D):
    """Please document this element"""

    meshType = MT.HEXA20
    nodes = (
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)),
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
    )
    elrefe = (
        ElrefeLoc(MT.H20, gauss=("RIGI=FPG8", "FPG1=FPG1"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.HE8, gauss=("RIGI=FPG8",)),
        ElrefeLoc(MT.QU8, gauss=("RIGI=FPG4",)),
    )
