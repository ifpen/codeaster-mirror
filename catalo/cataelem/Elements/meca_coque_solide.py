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

DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("DX", "DY", "DZ")), ("EN2", ("PINCH",))),
)

# For EPVC_ELGA
EDFVCPG = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELGA",
    location="RIGI",
    components=("EPTHER_L", "EPTHER_T", "EPTHER_N", "EPSECH", "EPHYDR", "EPPTOT"),
)

# For EPVC_ELNO
EDFVCNO = LocatedComponents(
    phys=PHY.EPSI_R,
    type="ELNO",
    components=("EPTHER_L", "EPTHER_T", "EPTHER_N", "EPSECH", "EPHYDR", "EPPTOT"),
)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)


# ----------------------------------------------------------------------------------------------
class MESSHELL_SB9(Element):
    """Solid-shell element on HEXA9 geometric support"""

    meshType = MT.HEXA9
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)), SetOfNodes("EN2", (9,)))
    elrefe = (
        ElrefeLoc(MT.HE9, gauss=("RIGI=LOB5", "FPG1=FPG1", "MASS=FPG8"), mater=("RIGI", "FPG1")),
        ElrefeLoc(MT.QU4, gauss=("RIGI=FPG4", "MASS=FPG4")),
    )
    calculs = (
        OP.CHAR_MECA_EFSU_R(
            te=-1,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PPRESSR, LC.CPRESBR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_EPSA_R(
            te=125,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FF3D3D(
            te=125,
            para_in=((SP.PFF3D3D, LC.CFOR3DF), (SP.PGEOMER, LC.EGEOM3D), (SP.PINSTR, LC.MTEMPSR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FR3D3D(
            te=125,
            para_in=((SP.PFR3D3D, LC.NFOR3DR), (SP.PGEOMER, LC.EGEOM3D)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_HYDR_R(
            te=125,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=125,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=125,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PPRESSR, LC.CPRESBR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=125,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PINSTR, LC.MTEMPSR),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SECH_R(
            te=125,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.COOR_ELGA(
            te=488,
            para_in=((SP.PGEOMER, LC.EGEOM3D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU3D),),
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
        OP.EPSI_ELGA(
            te=125,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PMATERC, LC.CMATERC), (SP.PDEPLAR, DDL_MECA)),
            para_out=((OP.EPSI_ELGA.PDEFOPG, LC.EGPS3DR),),
        ),
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, LC.EGPS3DR),),
            para_out=((SP.PDEFONO, LC.EEPS3DR),),
        ),
        OP.EPSL_ELGA(
            te=125,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPSL_ELGA.PDEFOPG, LC.EGPS3DR),),
        ),
        OP.EPSL_ELNO(
            te=4,
            para_in=((OP.EPSL_ELNO.PDEFOPG, LC.EGPS3DR),),
            para_out=((SP.PDEFONO, LC.EEPS3DR),),
        ),
        OP.EPVC_ELGA(
            te=125,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPVC_ELGA.PDEFOPG, EDFVCPG),),
        ),
        OP.EPVC_ELNO(
            te=4, para_in=((OP.EPVC_ELNO.PDEFOPG, EDFVCPG),), para_out=((SP.PDEFONO, EDFVCNO),)
        ),
        OP.FORC_NODA(
            te=125,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, LC.EGIG3DR),
                (SP.PDEPLAR, DDL_MECA),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=124,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PVARCMR, LC.ZVARCPG),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PCONTMR, LC.EGIG3DR),
                (SP.PVARIMP, LC.ZVARIPG),
                (SP.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (SP.PCONTPR, LC.EGIG3DR),
                (SP.PVARIPR, LC.ZVARIPG),
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INIT_VARC(te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG),)),
        OP.MASS_INER(
            te=65,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PMATERC, LC.CMATERC), (SP.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PMASSINE, LC.EMASSINE),),
        ),
        OP.MASS_MECA(
            te=125,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PMATERC, LC.CMATERC), (SP.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.NSPG_NBVA(
            te=496, para_in=((SP.PCOMPOR, LC.CCOMPO2),), para_out=((SP.PDCEL_I, LC.EDCEL_I),)
        ),
        OP.REFE_FORC_NODA(
            te=125,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PREFCO, LC.CRESSIG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_GEOM(
            te=125,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (OP.RIGI_GEOM.PCONTRR, LC.EGIG3DR)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA(
            te=125,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PMATERC, LC.CMATERC), (SP.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_TANG(
            te=124,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PVARCMR, LC.ZVARCPG),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PCONTMR, LC.EGIG3DR),
                (SP.PVARIMP, LC.ZVARIPG),
                (SP.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCONTPR, LC.EGIG3DR),
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.RAPH_MECA(
            te=124,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PVARCMR, LC.ZVARCPG),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PCONTMR, LC.EGIG3DR),
                (SP.PVARIMP, LC.ZVARIPG),
                (SP.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (SP.PCONTPR, LC.EGIG3DR),
                (SP.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.SIEF_ELGA(
            te=125,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PDEPLAR, DDL_MECA),
            ),
            para_out=((OP.SIEF_ELGA.PCONTRR, LC.EGIG3DR),),
        ),
        OP.SIEF_ELNO(
            te=4,
            para_in=((SP.PVARCPR, LC.ZVARCPG), (OP.SIEF_ELNO.PCONTRR, LC.EGIG3DR)),
            para_out=((OP.SIEF_ELNO.PSIEFNOR, LC.ESIG3DR),),
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
            te=546, para_in=((SP.PSIEFR, LC.EGIG3DR),), para_out=((SP.PSIGMR, LC.EGIG3DR),)
        ),
        OP.SIGM_ELNO(
            te=4,
            para_in=((OP.SIGM_ELNO.PCONTRR, LC.EGIG3DR),),
            para_out=((OP.SIGM_ELNO.PSIEFNOR, LC.ESIG3DR),),
        ),
        OP.TOU_INI_ELEM(
            te=99, para_out=((SP.PGEOM_R, LC.CGEOM3D), (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER))
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (SP.PGEOM_R, LC.EGGEO3D),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (OP.TOU_INI_ELGA.PSIEF_R, LC.EGIG3DR),
                (OP.TOU_INI_ELGA.PVARI_R, LC.ZVARIPG),
                (OP.TOU_INI_ELGA.PEPSI_R, LC.EGPS3DR),
                (OP.TOU_INI_ELGA.PTEMP_R, LC.ETEMPPG),
                (OP.TOU_INI_ELGA.PDOMMAG, LC.EDOMGGA),
                (SP.PFACY_R, LC.EGFC3DR),
                (OP.TOU_INI_ELGA.PINST_R, LC.EGINST_R),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (SP.PGEOM_R, LC.EGEOM3D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENNEUT_F),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENNEUT_R),
                (OP.TOU_INI_ELNO.PSIEF_R, LC.ESIG3DR),
                (OP.TOU_INI_ELNO.PVARI_R, LC.ZVARINO),
                (OP.TOU_INI_ELNO.PEPSI_R, LC.EEPS3DR),
                (OP.TOU_INI_ELNO.PDOMMAG, LC.EDOMGNO),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
            ),
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, LC.EGEOM3D),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )
