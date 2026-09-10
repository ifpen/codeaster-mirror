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
        ("EN1", ("HHO_FX[4]", "HHO_FY[4]")),
        ("EN2", ()),
        ("EN3", ("HHO_CX[10]", "HHO_CY[10]")),
    ),
)

CHHOGT = LocatedComponents(phys=PHY.N1920R, type="ELEM", components=("X[520]",))

CHHOST = LocatedComponents(phys=PHY.N1360R, type="ELEM", components=("X[351]",))

CHHOBS = LocatedComponents(
    phys=PHY.N3600R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("X[10]",)), ("EN2", ()), ("EN3", ("X[120]"))),
)


DEPLHHO = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY"))

PFONC = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z[10]",))

PFONCR = LocatedComponents(phys=PHY.NEUT_K8, type="ELEM", components=("Z[3]",))

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)

# ------------------------------------------------------------


class MECA_DPQ9_HHO333(Element):
    """Please document this element"""

    meshType = MT.QUAD9
    nodes = (
        SetOfNodes("EN1", (5, 6, 7, 8)),
        SetOfNodes("EN2", (1, 2, 3, 4)),
        SetOfNodes("EN3", (9,)),
    )
    elrefe = (
        ElrefeLoc(
            MT.QU9,
            gauss=("RIGI=FPG16", "FPG1=FPG1", "MTGA=FPG16", "MASS=FPG16"),
            mater=("RIGI", "FPG1", "MTGA", "MASS"),
        ),
    )
    calculs = (
        OP.AMOR_MECA(
            te=121,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMASSEL, MMATUUR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PRIGIEL, MMATUUR),
                (OP.AMOR_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.CHAR_MECA_FF2D2D(
            te=476,
            para_in=(
                (SP.PFF2D2D, LC.CFOR2DF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_FF2D2D.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FR2D2D(
            te=476,
            para_in=(
                (SP.PFR2D2D, LC.NFOR2DR),
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.CHAR_MECA_FR2D2D.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_HYDR_R(
            te=503,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_HYDR_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.CHAR_MECA_HYDR_R.PCHHOGT, CHHOGT),
                (OP.CHAR_MECA_HYDR_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PTOT_R(
            te=503,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_PTOT_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.CHAR_MECA_PTOT_R.PCHHOGT, CHHOGT),
                (OP.CHAR_MECA_PTOT_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SECH_R(
            te=503,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_SECH_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.CHAR_MECA_SECH_R.PCHHOGT, CHHOGT),
                (OP.CHAR_MECA_SECH_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=503,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.CHAR_MECA_TEMP_R.PCHHOGT, CHHOGT),
                (OP.CHAR_MECA_TEMP_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.COOR_ELGA(
            te=479,
            para_in=((SP.PGEOMER, LC.EGEOM2D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU2D),),
        ),
        OP.EPSI_ELGA(
            te=448,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.EPSI_ELGA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PDEFOPC, LC.EGPS2DC), (OP.EPSI_ELGA.PDEFOPG, LC.EGPS2DR)),
        ),
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONC, LC.EEPS2DC), (SP.PDEFONO, LC.EEPS2DR)),
        ),
        OP.FORC_NODA(
            te=450,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, LC.EGIG2DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.FORC_NODA.PCHHOGT, CHHOGT),
                (OP.FORC_NODA.PCHHOST, CHHOST),
                (OP.FORC_NODA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=455,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (OP.FULL_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA.PCONTMR, LC.EGIG2DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, LC.ZVARIPG),
                (OP.FULL_MECA.PVARIMR, LC.ZVARIPG),
                (OP.FULL_MECA.PCHHOGT, CHHOGT),
                (OP.FULL_MECA.PCHHOST, CHHOST),
                (OP.FULL_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA.PCONTPR, LC.EGIG2DR),
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.HHO_COND_MECA(
            te=507,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PMAELS1, MMATUUR), (SP.PVEELE1, MVECTUR)),
            para_out=(
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (SP.PMATUND, MMATUNS),
                (SP.PVECTUD, MVECTUR),
            ),
        ),
        OP.HHO_PRECALC_BS(
            te=494,
            para_in=((SP.PGEOMER, LC.EGEOM2D),),
            para_out=((OP.HHO_PRECALC_BS.PCHHOBO, CHHOBS),),
        ),
        OP.HHO_PRECALC_OP(
            te=460,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (OP.HHO_PRECALC_OP.PCHHOBS, CHHOBS)),
            para_out=((OP.HHO_PRECALC_OP.PCHHOGT, CHHOGT), (OP.HHO_PRECALC_OP.PCHHOST, CHHOST)),
        ),
        OP.HHO_DEPL_MECA(
            te=456,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PDEPLPR, DDL_MECA),
                (OP.HHO_DEPL_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_DEPL_MECA.PDEPL_R, DEPLHHO),),
        ),
        OP.HHO_CINE_F_MECA(
            te=458,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTPR, LC.MTEMPSR),
                (OP.HHO_CINE_F_MECA.PFONC, PFONC),
                (OP.HHO_CINE_F_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_CINE_F_MECA.PCINE, DDL_MECA),),
        ),
        OP.HHO_CINE_R_MECA(
            te=458,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.HHO_CINE_R_MECA.PCMPVALE, DEPLHHO),
                (OP.HHO_CINE_R_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_CINE_R_MECA.PCINE, DDL_MECA),),
        ),
        OP.HHO_PROJ_MECA(
            te=473,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.HHO_PROJ_MECA.PFUNC_R, PFONCR),
                (SP.PINSTPR, LC.MTEMPSR),
                (OP.HHO_PROJ_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_PROJ_MECA.PDEPL_R, DDL_MECA),),
        ),
        OP.INIT_VARC(
            te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG), (OP.INIT_VARC.PVARCNO, LC.ZVARCNO))
        ),
        OP.MASS_MECA(
            te=474,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
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
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.MATE_ELGA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.MATE_ELGA.PMATERR, LC.EGMATE_R),),
        ),
        OP.MATE_ELEM(
            te=142,
            para_in=(
                (SP.PMATERC, LC.CMATERC),
                (SP.PGEOMER, LC.EGEOM2D),
                (OP.MATE_ELEM.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.MATE_ELEM.PMATERR, LC.EEMATE_R),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.PILO_PRED_DEFO(
            te=495,
            para_in=(
                (OP.PILO_PRED_DEFO.PCOMPOR, LC.CCOMPOR),
                (OP.PILO_PRED_DEFO.PCONTMR, LC.EGIG2DR),
                (SP.PDDEPLR, DDL_MECA),
                (SP.PDEPL0R, DDL_MECA),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PDEPL1R, DDL_MECA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTYPEPI, LC.CTYPEPI),
                (OP.PILO_PRED_DEFO.PVARIMR, LC.ZVARIPG),
                (OP.PILO_PRED_DEFO.PCHHOGT, CHHOGT),
                (OP.PILO_PRED_DEFO.PCHHOBS, CHHOBS),
                (SP.PCDTAU, LC.CCDTAU),
            ),
            para_out=((OP.PILO_PRED_DEFO.PCOPILO, LC.ECOPILO),),
        ),
        OP.PILO_PRED_ELAS(
            te=495,
            para_in=(
                (SP.PBORNPI, LC.CBORNPI),
                (SP.PCDTAU, LC.CCDTAU),
                (OP.PILO_PRED_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.PILO_PRED_ELAS.PCONTMR, LC.EGIG2DR),
                (SP.PDDEPLR, DDL_MECA),
                (SP.PDEPL0R, DDL_MECA),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PDEPL1R, DDL_MECA),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTYPEPI, LC.CTYPEPI),
                (OP.PILO_PRED_ELAS.PVARIMR, LC.ZVARIPG),
                (OP.PILO_PRED_ELAS.PCHHOGT, CHHOGT),
                (OP.PILO_PRED_ELAS.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.PILO_PRED_ELAS.PCOPILO, LC.ECOPILO),),
        ),
        OP.RAPH_MECA(
            te=455,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (OP.RAPH_MECA.PCOMPOR, LC.CCOMPOR),
                (OP.RAPH_MECA.PCONTMR, LC.EGIG2DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RAPH_MECA.PVARIMR, LC.ZVARIPG),
                (OP.RAPH_MECA.PCHHOGT, CHHOGT),
                (OP.RAPH_MECA.PCHHOST, CHHOST),
                (OP.RAPH_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, LC.EGIG2DR),
                (OP.RAPH_MECA.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.REFE_FORC_NODA(
            te=450,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PRESIREF, LC.CRESIREF),
                (SP.PRESICMP, LC.CRESICMP),
                (OP.REFE_FORC_NODA.PCHHOGT, CHHOGT),
                (OP.REFE_FORC_NODA.PCHHOST, CHHOST),
                (OP.REFE_FORC_NODA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_MECA(
            te=455,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.RIGI_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_ELAS(
            te=455,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (OP.RIGI_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_ELAS.PCONTMR, LC.EGIG2DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RIGI_MECA_ELAS.PVARIMR, LC.ZVARIPG),
                (OP.RIGI_MECA_ELAS.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.RIGI_MECA_TANG(
            te=455,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (OP.RIGI_MECA_TANG.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_TANG.PCONTMR, LC.EGIG2DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PITERAT, LC.CITERAT),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RIGI_MECA_TANG.PVARIMR, LC.ZVARIPG),
                (OP.RIGI_MECA_TANG.PCHHOGT, CHHOGT),
                (OP.RIGI_MECA_TANG.PCHHOST, CHHOST),
                (OP.RIGI_MECA_TANG.PCHHOBS, CHHOBS),
            ),
            para_out=(
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, LC.EGIG2DR),
                (SP.PCODRET, LC.ECODRET),
                (SP.PCOPRED, LC.ECODRET),
            ),
        ),
        OP.SIEF_ELGA(
            te=504,
            para_in=(
                (SP.PCAMASS, LC.CCAMASS),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.SIEF_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.SIEF_ELGA.PCHHOGT, CHHOGT),
                (OP.SIEF_ELGA.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PCONTRC, LC.EGIG2DC), (OP.SIEF_ELGA.PCONTRR, LC.EGIG2DR)),
        ),
        OP.SIEF_ELNO(
            te=4,
            para_in=((OP.SIEF_ELNO.PCONTRR, LC.EGIG2DR), (OP.SIEF_ELNO.PVARCPR, LC.ZVARCPG)),
            para_out=((SP.PSIEFNOC, LC.ESIG2DC), (OP.SIEF_ELNO.PSIEFNOR, LC.ESIG2DR)),
        ),
        OP.SIEQ_ELGA(
            te=335,
            para_in=((OP.SIEQ_ELGA.PCONTRR, LC.EGIG2DR),),
            para_out=((OP.SIEQ_ELGA.PCONTEQ, LC.ECOEQPG),),
        ),
        OP.SIEQ_ELNO(
            te=335,
            para_in=((OP.SIEQ_ELNO.PCONTRR, LC.ESIG2DR),),
            para_out=((OP.SIEQ_ELNO.PCONTEQ, LC.ECOEQNO),),
        ),
        OP.SIGM_ELGA(
            te=546,
            para_in=((SP.PSIEFR, LC.EGIG2DR),),
            para_out=((SP.PSIGMC, LC.EGIG2DC), (SP.PSIGMR, LC.EGIG2DR)),
        ),
        OP.SIGM_ELNO(
            te=4,
            para_in=((OP.SIGM_ELNO.PCONTRR, LC.EGIG2DR),),
            para_out=((SP.PSIEFNOC, LC.ESIG2DC), (OP.SIGM_ELNO.PSIEFNOR, LC.ESIG2DR)),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PERREUR, LC.CERROR),
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
                (SP.PNEUT_F, LC.CNTINIF),
                (SP.PNEUT_R, LC.CNTINIR),
                (SP.PNEUT_I, LC.CNTINII),
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
                (OP.TOU_INI_ELGA.PSIEF_R, LC.EGIG2DR),
                (OP.TOU_INI_ELGA.PSOUR_R, LC.ESOURCR),
                (OP.TOU_INI_ELGA.PVARI_R, LC.ZVARIPG),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (OP.TOU_INI_ELGA.PNEUT_I, LC.EGTINII),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PDOMMAG, LC.EDOMGNO),
                (OP.TOU_INI_ELNO.PEPSI_R, LC.EEPS2DR),
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM2D),
                (OP.TOU_INI_ELNO.PINST_R, LC.ENINST_R),
                (OP.TOU_INI_ELNO.PSIEF_R, LC.ESIG2DR),
                (OP.TOU_INI_ELNO.PVARI_R, LC.ZVARINO),
                (OP.TOU_INI_ELNO.PMATE_R, LC.ENMATE_R),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENTINIF),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENTINIR),
                (OP.TOU_INI_ELNO.PNEUT_I, LC.ENTINII),
            ),
        ),
        OP.VARC_ELGA(
            te=530,
            para_in=((OP.VARC_ELGA.PVARCPR, LC.ZVARCPG),),
            para_out=((SP.PVARC_R, LC.EVARC_R),),
        ),
        OP.VARI_ELNO(
            te=4,
            para_in=((SP.PVARIGR, LC.ZVARIPG),),
            para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),),
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, LC.EGEOM2D),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
    )


# ------------------------------------------------------------


class MECA_DPT7_HHO333(MECA_DPQ9_HHO333):
    """Please document this element"""

    meshType = MT.TRIA7
    nodes = (SetOfNodes("EN1", (4, 5, 6)), SetOfNodes("EN2", (1, 2, 3)), SetOfNodes("EN3", (7,)))
    elrefe = (
        ElrefeLoc(
            MT.TR7,
            gauss=("RIGI=FPG12", "FPG1=FPG1", "MTGA=FPG12", "MASS=FPG12"),
            mater=("RIGI", "FPG1", "MTGA", "MASS"),
        ),
    )
