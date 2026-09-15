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
import cataelem.Commons.attributes as AT

# ----------------
# Modes locaux :
# ----------------


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(("EN1", ("HHO_FX[3]", "HHO_FY[3]", "HHO_FZ[3]")), ("EN2", ())),
)

CHHOBS = LocatedComponents(
    phys=PHY.N3600R, type="ELNO", diff=True, components=(("EN1", ("X[6]",)), ("EN2", ()))
)


DEPLHHO = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY", "DZ"))

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class MECA3DQU9_HHO1_F(Element):
    """Please document this element"""

    meshType = MT.QUAD9
    nodes = (SetOfNodes("EN1", (9,)), SetOfNodes("EN2", (1, 2, 3, 4, 5, 6, 7, 8)))
    attrs = ((AT.BORD_ISO, "OUI"),)
    elrefe = (ElrefeLoc(MT.QU9, gauss=("RIGI=FPG4",), mater=("RIGI",)),)
    calculs = (
        OP.CHAR_MECA_PRES_F(
            te=459,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPRESSF, LC.CPRE3DF),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_PRES_F.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRES_R(
            te=459,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPRESSR, LC.EPRE3DR),
                (OP.CHAR_MECA_PRES_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PRSU_R(
            te=459,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPRESSR, LC.EPRE3DR),
                (OP.CHAR_MECA_PRES_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FF2D3D(
            te=459,
            para_in=(
                (SP.PFF2D3D, LC.CFOR3DF),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_FF2D3D.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FR2D3D(
            te=459,
            para_in=(
                (SP.PFR2D3D, LC.NFOR3DR),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.CHAR_MECA_FR2D3D.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FRSU23(
            te=459,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PFR2D3D, LC.NFOR3DR),
                (SP.PGEOMER, LC.EGEOM3D),
                (OP.CHAR_MECA_FRSU23.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.COOR_ELGA(
            te=488,
            para_in=((SP.PGEOMER, LC.EGEOM3D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU3D),),
        ),
        OP.HHO_COND_MECA(
            te=574,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PMAELS1, MMATUUR), (SP.PVEELE1, MVECTUR)),
            para_out=(
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (SP.PMATUND, MMATUNS),
                (SP.PVECTUD, MVECTUR),
            ),
        ),
        OP.HHO_DEPL_MECA(
            te=427,
            para_in=(
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PDEPLPR, DDL_MECA),
                (OP.HHO_DEPL_MECA.PCHHOBS, CHHOBS),
            ),
            para_out=((OP.HHO_DEPL_MECA.PDEPL_R, DEPLHHO),),
        ),
        OP.INIT_VARC(
            te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG), (OP.INIT_VARC.PVARCNO, LC.ZVARCNO))
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.RIGI_MECA_PRSU_R(
            te=561,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPRESSR, LC.EPRE3DR),
                (OP.RIGI_MECA_PRSU_R.PCHHOBS, CHHOBS),
            ),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (SP.PFORC_R, LC.CFOR3DR),
                (OP.TOU_INI_ELEM.PPRES_R, LC.CPRE3DR),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
                (SP.PNEUT_F, LC.CNTINIF),
                (SP.PNEUT_R, LC.CNTINIR),
                (SP.PNEUT_I, LC.CNTINII),
            ),
        ),
        OP.TOU_INI_ELGA(
            te=99,
            para_out=(
                (OP.TOU_INI_ELGA.PDEPL_R, LC.EGDEP3D),
                (OP.TOU_INI_ELGA.PGEOM_R, LC.EGGEO3D),
                (OP.TOU_INI_ELGA.PPRES_R, LC.EPRESGA),
                (OP.TOU_INI_ELGA.PNEUT_F, LC.EGTINIF),
                (OP.TOU_INI_ELGA.PNEUT_R, LC.EGTINIR),
                (OP.TOU_INI_ELGA.PNEUT_I, LC.EGTINII),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM3D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENTINIF),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENTINIR),
                (OP.TOU_INI_ELNO.PNEUT_I, LC.ENTINII),
                (OP.TOU_INI_ELNO.PPRES_R, LC.EPRE3DR),
                (OP.TOU_INI_ELNO.PSIEF_R, LC.ESIG3DR),
                (OP.TOU_INI_ELNO.PFORC_R, LC.ENFOR3R),
            ),
        ),
    )


# ------------------------------------------------------------
class MECA3DTR7_HHO1_F(MECA3DQU9_HHO1_F):
    """Please document this element"""

    meshType = MT.TRIA7
    nodes = (SetOfNodes("EN1", (7,)), SetOfNodes("EN2", (1, 2, 3, 4, 5, 6)))
    attrs = ((AT.BORD_ISO, "OUI"),)
    elrefe = (ElrefeLoc(MT.TR7, gauss=("RIGI=FPG3",), mater=("RIGI",)),)


# ------------------------------------------------------------
class MECA3DGVT_HHO1_F(MECA3DTR7_HHO1_F):
    """Please document this element"""

    calculs = (OP.HHO_COND_MECA(te=-1),)


# ------------------------------------------------------------
class MECA3DGVQ_HHO1_F(MECA3DQU9_HHO1_F):
    """Please document this element"""

    calculs = (OP.HHO_COND_MECA(te=-1),)
