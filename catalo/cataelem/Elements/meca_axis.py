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
NDEPLAC = LocatedComponents(phys=PHY.DEPL_C, type="ELNO", components=("DX", "DY"))

DDL_MECA = LocatedComponents(phys=PHY.DEPL_R, type="ELNO", components=("DX", "DY"))

MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MVECZZR = ArrayOfComponents(phys=PHY.VSIZ_R, locatedComponents=LC.DDL_NOZ1)

MMATUUC = ArrayOfComponents(phys=PHY.MDEP_C, locatedComponents=NDEPLAC)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)

MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)

MMATZZR = ArrayOfComponents(phys=PHY.MSIZ_R, locatedComponents=LC.DDL_NOZ1)


# ------------------------------------------------------------
class MEAXQU4(Element):
    """Mechanics - Axisymmetric - QUAD4"""

    meshType = MT.QUAD4
    elrefe = (
        ElrefeLoc(
            MT.QU4,
            gauss=("RIGI=FPG4", "FPG1=FPG1", "MASS=FPG4", "NOEU_S=NOEU_S", "NOEU=NOEU"),
            mater=("RIGI", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )
    calculs = (
        OP.ADD_SIGM(
            te=581,
            para_in=((SP.PEPCON1, LC.EGIG2DR), (SP.PEPCON2, LC.EGIG2DR)),
            para_out=((SP.PEPCON3, LC.EGIG2DR),),
        ),
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
        OP.CALC_ESTI_ERRE(
            te=291,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.CALC_ESTI_ERRE.PSIEF_R, LC.EGIG2DR),
                (SP.PSIGMA, LC.ESIG2DR),
                (OP.CALC_ESTI_ERRE.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.CALC_ESTI_ERRE.PERREUR, LC.CERROR),),
        ),
        OP.CALC_G(
            te=222,
            para_in=(
                (SP.PACCELE, DDL_MECA),
                (OP.CALC_G.PCOMPOR, LC.CCOMPOR),
                (SP.PCONTGR, LC.EGIG2DR),
                (OP.CALC_G.PCONTRR, LC.EGIG2DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PEPSINR, LC.EEPS2DR),
                (SP.PFRVOLU, LC.NFOR2DR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, LC.ESIG2DR),
                (OP.CALC_G.PTHETAR, LC.ETHETA),
                (OP.CALC_G.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVITESS, DDL_MECA),
            ),
            para_out=((SP.PGTHETA, LC.CKGTHET),),
        ),
        OP.CALC_G_F(
            te=222,
            para_in=(
                (SP.PACCELE, DDL_MECA),
                (OP.CALC_G_F.PCOMPOR, LC.CCOMPOR),
                (SP.PCONTGR, LC.EGIG2DR),
                (OP.CALC_G_F.PCONTRR, LC.EGIG2DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PEPSINF, LC.CEPS2DF),
                (SP.PFFVOLU, LC.CFOR2DF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, LC.ESIG2DR),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CALC_G_F.PTHETAR, LC.ETHETA),
                (OP.CALC_G_F.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVITESS, DDL_MECA),
            ),
            para_out=((SP.PGTHETA, LC.CKGTHET),),
        ),
        OP.CALC_K_G(
            te=222,
            para_in=(
                (OP.CALC_K_G.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                #                     (SP.PFISSR, LC.CFISSR),
                (SP.PEPSINR, LC.EEPS2DR),
                (SP.PFRVOLU, LC.NFOR2DR),
                (OP.CALC_K_G.PBASLOR, LC.N6NEUT_R),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PPULPRO, LC.CFREQR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, LC.ESIG2DR),
                (OP.CALC_K_G.PTHETAR, LC.ETHETA),
                (OP.CALC_K_G.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PGTHETA, LC.CKGTHET),),
        ),
        OP.CALC_K_G_F(
            te=222,
            para_in=(
                (OP.CALC_K_G_F.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PEPSINF, LC.CEPS2DF),
                (SP.PFFVOLU, LC.CFOR2DF),
                (OP.CALC_K_G_F.PBASLOR, LC.N6NEUT_R),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PPULPRO, LC.CFREQR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, LC.ESIG2DR),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CALC_K_G_F.PTHETAR, LC.ETHETA),
                (OP.CALC_K_G_F.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PGTHETA, LC.CKGTHET),),
        ),
        OP.CALC_G_XFEM(
            te=96,
            para_in=(
                (SP.PACCELE, DDL_MECA),
                (OP.CALC_G_XFEM.PCOMPOR, LC.CCOMPOR),
                (SP.PCONTGR, LC.EGIG2DR),
                (OP.CALC_G_XFEM.PCONTRR, LC.EGIG2DR),
                (SP.PDEFOPL, LC.EEPS2DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PEPSINR, LC.EEPS2DR),
                (SP.PFRVOLU, LC.NFOR2DR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, LC.ESIG2DR),
                (SP.PTHETAR, DDL_MECA),
                (OP.CALC_G_XFEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.CALC_G_XFEM.PVARIPR, LC.ZVARINO),
                (SP.PVITESS, DDL_MECA),
            ),
            para_out=((SP.PGTHETA, LC.CGTHETA),),
        ),
        OP.CALC_G_XFEM_F(
            te=96,
            para_in=(
                (SP.PACCELE, DDL_MECA),
                (OP.CALC_G_XFEM_F.PCOMPOR, LC.CCOMPOR),
                (SP.PCONTGR, LC.EGIG2DR),
                (OP.CALC_G_XFEM_F.PCONTRR, LC.EGIG2DR),
                (SP.PDEFOPL, LC.EEPS2DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PEPSINF, LC.CEPS2DF),
                (SP.PFFVOLU, LC.CFOR2DF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, LC.ESIG2DR),
                (SP.PINSTR, LC.MTEMPSR),
                (SP.PTHETAR, DDL_MECA),
                (OP.CALC_G_XFEM_F.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.CALC_G_XFEM_F.PVARIPR, LC.ZVARINO),
                (SP.PVITESS, DDL_MECA),
            ),
            para_out=((SP.PGTHETA, LC.CGTHETA),),
        ),
        OP.CALC_K_G_XFEM(
            te=299,
            para_in=(
                (OP.CALC_K_G_XFEM.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PFISSR, LC.CFISSR),
                (SP.PFRVOLU, LC.NFOR2DR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PPULPRO, LC.CFREQR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, LC.ESIG2DR),
                (SP.PTHETAR, DDL_MECA),
                (OP.CALC_K_G_XFEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PGTHETA, LC.CKGTX2D),),
        ),
        OP.CALC_K_G_XFEM_F(
            te=299,
            para_in=(
                (OP.CALC_K_G_XFEM_F.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PFFVOLU, LC.CFOR2DF),
                (SP.PFISSR, LC.CFISSR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PPULPRO, LC.CFREQR),
                (SP.PROTATR, LC.CROTATR),
                (SP.PSIGINR, LC.ESIG2DR),
                (SP.PINSTR, LC.MTEMPSR),
                (SP.PTHETAR, DDL_MECA),
                (OP.CALC_K_G_XFEM_F.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PGTHETA, LC.CKGTX2D),),
        ),
        OP.CALC_NOEU_BORD(
            te=290, para_in=((SP.PGEOMER, LC.EGEOM2D),), para_out=((SP.PVECTUR, MVECTUR),)
        ),
        OP.CHAR_MECA_EPSA_R(
            te=421,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_EPSA_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_EPSI_F(
            te=284,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PEPSINF, LC.CEPS2DF),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_EPSI_F.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_EPSI_R(
            te=284,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PEPSINR, LC.EGPS2DR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.CHAR_MECA_EPSI_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FF2D2D(
            te=94,
            para_in=((SP.PFF2D2D, LC.CFOR2DF), (SP.PGEOMER, LC.EGEOM2D), (SP.PINSTR, LC.MTEMPSR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_FR2D2D(
            te=93,
            para_in=((SP.PFR2D2D, LC.NFOR2DR), (SP.PGEOMER, LC.EGEOM2D)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_HYDR_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_HYDR_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_META_Z(
            te=353,
            para_in=(
                (OP.CHAR_MECA_META_Z.PCOMPOR, LC.CCOMPOR),
                (OP.CHAR_MECA_META_Z.PCONTMR, LC.EGIG2DR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.CHAR_MECA_META_Z.PVARCPR, LC.ZVARCPG),
                (OP.CHAR_MECA_META_Z.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_PESA_R(
            te=85,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (OP.CHAR_MECA_PESA_R.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_ROTA_R(
            te=84,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PMATERC, LC.CMATERC), (SP.PROTATR, LC.CROTATR)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_SECH_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_SECH_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.CHAR_MECA_TEMP_R(
            te=13,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.CHAR_MECA_TEMP_R.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.COOR_ELGA(
            te=479,
            para_in=((SP.PGEOMER, LC.EGEOM2D),),
            para_out=((OP.COOR_ELGA.PCOORPG, LC.EGGAU2D),),
        ),
        OP.DERA_ELGA(
            te=489,
            para_in=(
                (OP.DERA_ELGA.PCOMPOR, LC.CCOMPOR),
                (OP.DERA_ELGA.PCONTMR, LC.EGIG2DR),
                (OP.DERA_ELGA.PCONTPR, LC.EGIG2DR),
                (SP.PDERAMG, LC.EDERAPG),
                (SP.PMATERC, LC.CMATERC),
                (OP.DERA_ELGA.PVARCPR, LC.ZVARCPG),
                (OP.DERA_ELGA.PVARIMR, LC.ZVARIPG),
                (OP.DERA_ELGA.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((OP.DERA_ELGA.PDERAPG, LC.EDERAPG),),
        ),
        OP.DERA_ELNO(
            te=4,
            para_in=((OP.DERA_ELNO.PDERAPG, LC.EDERAPG),),
            para_out=((SP.PDERANO, LC.EDERANO),),
        ),
        OP.ECIN_ELEM(
            te=82,
            para_in=(
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.POMEGA2, LC.COMEG2R),
                (OP.ECIN_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVITESR, DDL_MECA),
            ),
            para_out=((SP.PENERCR, LC.CENEISO),),
        ),
        OP.ENDO_ELGA(
            te=512,
            para_in=(
                (OP.ENDO_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PCONTGP, LC.EGIG2DR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PTRIAGM, LC.ETRIAPG),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.ENDO_ELGA.PVARCPR, LC.ZVARCPG),
                (OP.ENDO_ELGA.PVARIMR, LC.ZVARIPG),
                (OP.ENDO_ELGA.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((OP.ENDO_ELGA.PTRIAPG, LC.ETRIAPG),),
        ),
        OP.ENDO_ELNO(
            te=4,
            para_in=((OP.ENDO_ELNO.PTRIAPG, LC.ETRIAPG),),
            para_out=((SP.PTRIANO, LC.ETRIANO),),
        ),
        OP.ENEL_ELEM(
            te=490,
            para_in=(
                (OP.ENEL_ELEM.PCOMPOR, LC.CCOMPOR),
                (OP.ENEL_ELEM.PCONTPR, LC.EGIG2DR),
                (SP.PDEPLR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.ENEL_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.ENEL_ELEM.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((SP.PENERD1, LC.CENEISO),),
        ),
        OP.ENTR_ELEM(
            te=490,
            para_in=(
                (OP.ENTR_ELEM.PCOMPOR, LC.CCOMPOR),
                (OP.ENTR_ELEM.PCONTPR, LC.EGIG2DR),
                (SP.PDEPLR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.ENTR_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.ENTR_ELEM.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((SP.PENTRD1, LC.CENEISO),),
        ),
        OP.ENEL_ELGA(
            te=575,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (OP.ENEL_ELGA.PCOMPOR, LC.CCOMPOR),
                (OP.ENEL_ELGA.PCONTRR, LC.EGIG2DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.ENEL_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIGR, LC.ZVARIPG),
            ),
            para_out=((OP.ENEL_ELGA.PENERDR, LC.EENEISO),),
        ),
        OP.ENEL_ELNO(
            te=4,
            para_in=((OP.ENEL_ELNO.PENERPG, LC.EENEISO),),
            para_out=((SP.PENERNO, LC.NENEISO),),
        ),
        OP.ENER_TOTALE(
            te=490,
            para_in=(
                (OP.ENER_TOTALE.PCOMPOR, LC.CCOMPOR),
                (OP.ENER_TOTALE.PCONTMR, LC.EGIG2DR),
                (OP.ENER_TOTALE.PCONTPR, LC.EGIG2DR),
                (SP.PDEPLM, DDL_MECA),
                (SP.PDEPLR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.ENER_TOTALE.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.ENER_TOTALE.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((SP.PENERD1, LC.CENEISO),),
        ),
        OP.EPEQ_ELGA(
            te=335,
            para_in=((OP.EPEQ_ELGA.PDEFORR, LC.EGPS2DR),),
            para_out=((OP.EPEQ_ELGA.PDEFOEQ, LC.EDFEQPG),),
        ),
        OP.EPEQ_ELNO(
            te=335,
            para_in=((OP.EPEQ_ELNO.PDEFORR, LC.EEPS2DR),),
            para_out=((OP.EPEQ_ELNO.PDEFOEQ, LC.EDFEQNO),),
        ),
        OP.EPGQ_ELGA(
            te=335,
            para_in=((OP.EPGQ_ELGA.PDEFORR, LC.EGPS2DR),),
            para_out=((OP.EPGQ_ELGA.PDEFOEQ, LC.EDFEQPG),),
        ),
        OP.EPGQ_ELNO(
            te=335,
            para_in=((OP.EPGQ_ELNO.PDEFORR, LC.EEPS2DR),),
            para_out=((OP.EPGQ_ELNO.PDEFOEQ, LC.EDFEQNO),),
        ),
        OP.EPFD_ELGA(
            te=528,
            para_in=(
                (OP.EPFD_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PVARIGR, LC.ZVARIPG),
            ),
            para_out=((OP.EPFD_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPFD_ELNO(
            te=4,
            para_in=((OP.EPFD_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPFP_ELGA(
            te=528,
            para_in=(
                (OP.EPFP_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EPFP_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIGR, LC.ZVARIPG),
            ),
            para_out=((OP.EPFP_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPFP_ELNO(
            te=4,
            para_in=((OP.EPFP_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPME_ELGA(
            te=87,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (OP.EPME_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EPME_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPME_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPME_ELNO(
            te=4,
            para_in=((OP.EPME_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPMG_ELGA(
            te=87,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (OP.EPMG_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EPMG_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPMG_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPMG_ELNO(
            te=4,
            para_in=((OP.EPMG_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPMQ_ELGA(
            te=335,
            para_in=((OP.EPMQ_ELGA.PDEFORR, LC.EGPS2DR),),
            para_out=((OP.EPMQ_ELGA.PDEFOEQ, LC.EDFEQPG),),
        ),
        OP.EPMQ_ELNO(
            te=335,
            para_in=((OP.EPMQ_ELNO.PDEFORR, LC.EEPS2DR),),
            para_out=((OP.EPMQ_ELNO.PDEFOEQ, LC.EDFEQNO),),
        ),
        OP.EPOT_ELEM(
            te=286,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPOT_ELEM.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPOT_ELEM.PENERDR, LC.CENEISO),),
        ),
        OP.EPSG_ELGA(
            te=87,
            para_in=(
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EPSG_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPSG_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPSG_ELNO(
            te=4,
            para_in=((OP.EPSG_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPSL_ELGA(
            te=87,
            para_in=(
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EPSL_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPSL_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPSL_ELNO(
            te=4,
            para_in=((OP.EPSL_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPSI_ELGA(
            te=87,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (OP.EPSI_ELGA.PCOMPOR, LC.CCOMPOR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
            ),
            para_out=((SP.PDEFOPC, LC.EGPS2DC), (OP.EPSI_ELGA.PDEFOPG, LC.EGPS2DR)),
        ),
        OP.EPSI_ELNO(
            te=4,
            para_in=((OP.EPSI_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONC, LC.EEPS2DC), (SP.PDEFONO, LC.EEPS2DR)),
        ),
        OP.EPSP_ELGA(
            te=334,
            para_in=(
                (OP.EPSP_ELGA.PCOMPOR, LC.CCOMPOR),
                (OP.EPSP_ELGA.PCONTRR, LC.EGIG2DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.EPSP_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIGR, LC.ZVARIPG),
            ),
            para_out=((OP.EPSP_ELGA.PDEFOPG, LC.EGPS2DR),),
        ),
        OP.EPSP_ELNO(
            te=4,
            para_in=((OP.EPSP_ELNO.PDEFOPG, LC.EGPS2DR),),
            para_out=((SP.PDEFONO, LC.EEPS2DR),),
        ),
        OP.EPVC_ELGA(
            te=529,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.EPVC_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
            ),
            para_out=((OP.EPVC_ELGA.PDEFOPG, LC.EGVARC3D),),
        ),
        OP.EPVC_ELNO(
            te=4,
            para_in=((OP.EPVC_ELNO.PDEFOPG, LC.EGVARC3D),),
            para_out=((SP.PDEFONO, LC.NVARC3D),),
        ),
        OP.ERME_ELEM(
            te=377,
            para_in=(
                (SP.PCONTNO, LC.ESIG2DR),
                (SP.PFFVOLU, LC.CFOR2DF),
                (SP.PFORCE, LC.CREFERI),
                (SP.PFRVOLU, LC.EFOR2DR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PPESANR, LC.CPESANR),
                (SP.PPRESS, LC.CREFERI),
                (SP.PROTATR, LC.CROTATR),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.ERME_ELEM.PVOISIN, LC.EVOISIN),
            ),
            para_out=((OP.ERME_ELEM.PERREUR, LC.CERROR),),
        ),
        OP.ERME_ELNO(
            te=379,
            para_in=((OP.ERME_ELNO.PERREUR, LC.CERROR),),
            para_out=((SP.PERRENO, LC.NERROR),),
        ),
        OP.ERRE_QIZZ(
            te=292,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PSIEFD_R, LC.EGIG2DR),
                (SP.PSIEFP_R, LC.EGIG2DR),
                (SP.PSIGMAD, LC.ESIG2DR),
                (SP.PSIGMAP, LC.ESIG2DR),
                (OP.ERRE_QIZZ.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((OP.ERRE_QIZZ.PERREUR, LC.CERROR),),
        ),
        OP.ETOT_ELEM(
            te=575,
            para_in=(
                (OP.ETOT_ELEM.PCONTMR, LC.EGIG2DR),
                (OP.ETOT_ELEM.PCONTPR, LC.EGIG2DR),
                (SP.PDEPLM, DDL_MECA),
                (SP.PDEPLR, DDL_MECA),
                (OP.ETOT_ELEM.PENERDM, LC.CENEISO),
                (SP.PGEOMER, LC.EGEOM2D),
            ),
            para_out=((OP.ETOT_ELEM.PENERDR, LC.CENEISO),),
        ),
        OP.ETOT_ELGA(
            te=575,
            para_in=(
                (OP.ETOT_ELGA.PCONTMR, LC.EGIG2DR),
                (OP.ETOT_ELGA.PCONTPR, LC.EGIG2DR),
                (SP.PDEPLM, DDL_MECA),
                (SP.PDEPLR, DDL_MECA),
                (OP.ETOT_ELGA.PENERDM, LC.EENEISO),
                (SP.PGEOMER, LC.EGEOM2D),
            ),
            para_out=((OP.ETOT_ELGA.PENERDR, LC.EENEISO),),
        ),
        OP.ETOT_ELNO(
            te=4,
            para_in=((OP.ETOT_ELNO.PENERPG, LC.EENEISO),),
            para_out=((SP.PENERNO, LC.NENEISO),),
        ),
        OP.FORC_NODA(
            te=8,
            para_in=(
                (SP.PCOMPOR, LC.CCOMPOR),
                (SP.PSIEFR, LC.EGIG2DR),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.FULL_MECA(
            te=139,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
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
        OP.FULL_MECA_ELAS(
            te=139,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (OP.FULL_MECA_ELAS.PCOMPOR, LC.CCOMPOR),
                (OP.FULL_MECA_ELAS.PCONTMR, LC.EGIG2DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.FULL_MECA_ELAS.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (SP.PVARIMP, LC.ZVARIPG),
                (OP.FULL_MECA_ELAS.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.FULL_MECA_ELAS.PCONTPR, LC.EGIG2DR),
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (OP.FULL_MECA_ELAS.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.INDIC_ENER(
            te=490,
            para_in=(
                (OP.INDIC_ENER.PCOMPOR, LC.CCOMPOR),
                (OP.INDIC_ENER.PCONTPR, LC.EGIG2DR),
                (SP.PDEPLR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.INDIC_ENER.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.INDIC_ENER.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((SP.PENERD1, LC.CENEISO), (SP.PENERD2, LC.CENEISO)),
        ),
        OP.INDIC_SEUIL(
            te=490,
            para_in=(
                (OP.INDIC_SEUIL.PCOMPOR, LC.CCOMPOR),
                (OP.INDIC_SEUIL.PCONTPR, LC.EGIG2DR),
                (SP.PDEPLR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.INDIC_SEUIL.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.INDIC_SEUIL.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((SP.PENERD1, LC.CENEISO), (SP.PENERD2, LC.CENEISO)),
        ),
        OP.INDL_ELGA(
            te=30,
            para_in=(
                (OP.INDL_ELGA.PCOMPOR, LC.CCOMPOR),
                (OP.INDL_ELGA.PCONTPR, LC.EGIG2DR),
                (SP.PMATERC, LC.CMATERC),
                (OP.INDL_ELGA.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((SP.PINDLOC, LC.EGINDLO),),
        ),
        OP.INIT_MAIL_VOIS(te=99, para_out=((OP.INIT_MAIL_VOIS.PVOISIN, LC.EVOISIN),)),
        OP.INIT_VARC(
            te=99, para_out=((OP.INIT_VARC.PVARCPR, LC.ZVARCPG), (OP.INIT_VARC.PVARCNO, LC.ZVARCNO))
        ),
        OP.MASS_INER(
            te=285,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_INER.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMASSINE, LC.EMASSINE),),
        ),
        OP.MASS_MECA(
            te=82,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.MASS_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.MASS_ZZ1(te=293, para_in=((SP.PGEOMER, LC.EGEOM2D),), para_out=((SP.PMATZZR, MMATZZR),)),
        OP.M_GAMMA(
            te=82,
            para_in=(
                (SP.PACCELR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.M_GAMMA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.NORME_FROB(
            te=563,
            para_in=(
                (SP.PCALCI, LC.EMNEUT_I),
                (SP.PCHAMPG, LC.EGTINIR),
                (SP.PCOEFR, LC.CNORMCF),
                (OP.NORME_FROB.PCOORPG, LC.EGGAU2D),
            ),
            para_out=((SP.PNORME, LC.ENORME),),
        ),
        OP.NORME_L2(
            te=563,
            para_in=(
                (SP.PCALCI, LC.EMNEUT_I),
                (SP.PCHAMPG, LC.EGTINIR),
                (SP.PCOEFR, LC.CNORMCF),
                (OP.NORME_L2.PCOORPG, LC.EGGAU2D),
            ),
            para_out=((SP.PNORME, LC.ENORME),),
        ),
        OP.NSPG_NBVA(
            te=496,
            para_in=((OP.NSPG_NBVA.PCOMPOR, LC.CCOMPO2),),
            para_out=((SP.PDCEL_I, LC.EDCEL_I),),
        ),
        OP.PAS_COURANT(
            te=404,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.PAS_COURANT.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PCOURAN, LC.ECOURAN),),
        ),
        OP.PILO_PRED_DEFO(
            te=543,
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
            ),
            para_out=((OP.PILO_PRED_DEFO.PCOPILO, LC.ECOPILO),),
        ),
        OP.PILO_PRED_ELAS(
            te=543,
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
            ),
            para_out=((OP.PILO_PRED_ELAS.PCOPILO, LC.ECOPILO),),
        ),
        OP.QIRE_ELEM(
            te=378,
            para_in=(
                (SP.PCONSTR, LC.CCONSTR),
                (SP.PCONTNOD, LC.ESIG2DR),
                (SP.PCONTNOP, LC.ESIG2DR),
                (SP.PFFVOLUD, LC.CFOR2DF),
                (SP.PFFVOLUP, LC.CFOR2DF),
                (SP.PFORCED, LC.CREFERI),
                (SP.PFORCEP, LC.CREFERI),
                (SP.PFRVOLUD, LC.EFOR2DR),
                (SP.PFRVOLUP, LC.EFOR2DR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PPESANRD, LC.CPESANR),
                (SP.PPESANRP, LC.CPESANR),
                (SP.PPRESSD, LC.CREFERI),
                (SP.PPRESSP, LC.CREFERI),
                (SP.PROTATRD, LC.CROTATR),
                (SP.PROTATRP, LC.CROTATR),
                (SP.PINSTR, LC.MTEMPSR),
                (OP.QIRE_ELEM.PVOISIN, LC.EVOISIN),
            ),
            para_out=((OP.QIRE_ELEM.PERREUR, LC.CERROR),),
        ),
        OP.QIRE_ELNO(
            te=379,
            para_in=((OP.QIRE_ELNO.PERREUR, LC.CERROR),),
            para_out=((SP.PERRENO, LC.NERROR),),
        ),
        OP.RAPH_MECA(
            te=139,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
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
                (SP.PVARIMP, LC.ZVARIPG),
                (OP.RAPH_MECA.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, LC.EGIG2DR),
                (OP.RAPH_MECA.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.RAPH_MECA_IMPLEX(
            te=139,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
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
                (SP.PVARIMP, LC.ZVARIPG),
                (OP.RAPH_MECA.PVARIMR, LC.ZVARIPG),
            ),
            para_out=(
                (SP.PCODRET, LC.ECODRET),
                (OP.RAPH_MECA.PCONTPR, LC.EGIG2DR),
                (OP.RAPH_MECA.PVARIPR, LC.ZVARIPG),
                (SP.PVECTUR, MVECTUR),
            ),
        ),
        OP.REFE_FORC_NODA(
            te=8,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (SP.PREFCO, LC.CRESSIG)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.REPERE_LOCAL(
            te=133,
            para_in=((SP.PCAMASS, LC.CCAMA2D), (SP.PGEOMER, LC.EGEOM2D)),
            para_out=((SP.PREPLO1, LC.CGEOM2D), (SP.PREPLO2, LC.CGEOM2D)),
        ),
        OP.REST_ECRO(
            te=116,
            para_in=(
                (OP.REST_ECRO.PCOMPOR, LC.CCOMPOR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.REST_ECRO.PVARCPR, LC.ZVARCPG),
                (OP.REST_ECRO.PVARIMR, LC.ZVARIPG),
            ),
            para_out=((OP.REST_ECRO.PVARIPR, LC.ZVARIPG),),
        ),
        OP.RICE_TRACEY(
            te=332,
            para_in=(
                (OP.RICE_TRACEY.PCOMPOR, LC.CCOMPOR),
                (OP.RICE_TRACEY.PCONTPR, LC.EGIG2DR),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PSDRMR, LC.EGNEUT1R),
                (SP.PSOUSOP, LC.CSOUSOP),
                (OP.RICE_TRACEY.PVARIMR, LC.ZVARIPG),
                (OP.RICE_TRACEY.PVARIPR, LC.ZVARIPG),
            ),
            para_out=((SP.PRICTRA, LC.ERICTRA), (SP.PSDRPR, LC.EGNEUT1R)),
        ),
        OP.RIGI_MECA(
            te=81,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.RIGI_MECA.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_ELAS(
            te=139,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
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
            ),
            para_out=((SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.RIGI_GEOM(
            te=92,
            para_in=((OP.RIGI_GEOM.PCONTRR, LC.EGIG2DR), (SP.PGEOMER, LC.EGEOM2D)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_MECA_HYST(
            te=230,
            para_in=(
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PRIGIEL, MMATUUR),
                (OP.RIGI_MECA_HYST.PVARCPR, LC.ZVARCPG),
            ),
            para_out=((SP.PMATUUC, MMATUUC),),
        ),
        OP.RIGI_MECA_IMPLEX(
            te=139,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PCARCRI, LC.CCARCRI),
                (SP.PMULCOM, LC.CMLCOMP),
                (OP.RIGI_MECA_IMPLEX.PCOMPOR, LC.CCOMPOR),
                (OP.RIGI_MECA_IMPLEX.PCONTMR, LC.EGIG2DR),
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PINSTMR, LC.MTEMPSR),
                (SP.PINSTPR, LC.MTEMPSR),
                (SP.PMATERC, LC.CMATERC),
                (SP.PVARCMR, LC.ZVARCPG),
                (OP.RIGI_MECA_IMPLEX.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
                (OP.RIGI_MECA_IMPLEX.PVARIMR, LC.ZVARIPG),
            ),
            para_out=((SP.PCONTXR, LC.EGIG2DR), (SP.PMATUNS, MMATUNS), (SP.PMATUUR, MMATUUR)),
        ),
        OP.RIGI_MECA_TANG(
            te=139,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
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
            ),
            para_out=(
                (SP.PMATUNS, MMATUNS),
                (SP.PMATUUR, MMATUUR),
                (SP.PVECTUR, MVECTUR),
                (OP.RIGI_MECA_TANG.PCONTPR, LC.EGIG2DR),
                (SP.PCOPRED, LC.ECODRET),
                (SP.PCODRET, LC.ECODRET),
            ),
        ),
        OP.SECM_ZZ1(
            te=294,
            para_in=((SP.PGEOMER, LC.EGEOM2D), (OP.SECM_ZZ1.PSIEF_R, LC.EGIG2DR)),
            para_out=(
                (SP.PVECTR1, MVECZZR),
                (SP.PVECTR2, MVECZZR),
                (SP.PVECTR3, MVECZZR),
                (SP.PVECTR4, MVECZZR),
            ),
        ),
        OP.SIEF_ELGA(
            te=22,
            para_in=(
                (SP.PCAMASS, LC.CCAMA2D),
                (SP.PDEPLAR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (OP.SIEF_ELGA.PVARCPR, LC.ZVARCPG),
                (SP.PVARCRR, LC.ZVARCPG),
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
        OP.SIMY_ELGA(
            te=6,
            para_in=((OP.SIMY_ELGA.PCONTRR, LC.EGIG2DR), (SP.PGEOMER, LC.EGEOM2D)),
            para_out=((OP.SIMY_ELGA.PSIEFNOR, LC.EGIG2DR),),
        ),
        OP.SING_ELEM(te=99, para_out=((SP.PSING_R, LC.ESINGUL),)),
        OP.SING_ELNO(te=99, para_out=((SP.PSINGNO, LC.ESINGNO),)),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM2D), (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER)),
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
                (SP.PTEMPN_R, LC.ETEMPNO),
                (OP.TOU_INI_ELNO.PSIEF_R, LC.ESIG2DR),
                (OP.TOU_INI_ELNO.PVARI_R, LC.ZVARINO),
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
            te=4,
            para_in=((SP.PVARIGR, LC.ZVARIPG),),
            para_out=((OP.VARI_ELNO.PVARINR, LC.ZVARINO),),
        ),
        OP.VERI_JACOBIEN(
            te=328, para_in=((SP.PGEOMER, LC.EGEOM2D),), para_out=((SP.PCODRET, LC.ECODRET),)
        ),
        OP.WEIBULL(
            te=331,
            para_in=(
                (OP.WEIBULL.PCOMPOR, LC.CCOMPOR),
                (SP.PCONTRG, LC.EGIG2DR),
                (OP.WEIBULL.PDEFORR, LC.EGPS2DR),
                (OP.WEIBULL.PDOMMAG, LC.EDOMGGA),
                (SP.PGEOMER, LC.EGEOM2D),
                (SP.PMATERC, LC.CMATERC),
                (SP.PSOUSOP, LC.CSOUSOP),
                (OP.WEIBULL.PVARCPR, LC.ZVARCPG),
                (SP.PVARIPG, LC.ZVARIPG),
            ),
            para_out=((SP.PSIGISG, LC.EDOMGGA), (SP.PWEIBUL, LC.EWEIBUL)),
        ),
    )


# ------------------------------------------------------------
class MEAXQU8(MEAXQU4):
    """Mechanics - Axisymmetric - QUAD4"""

    meshType = MT.QUAD8
    elrefe = (
        ElrefeLoc(
            MT.QU8,
            gauss=("RIGI=FPG9", "FPG1=FPG1", "MASS=FPG9", "NOEU_S=NOEU_S", "NOEU=NOEU"),
            mater=("RIGI", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )


# ------------------------------------------------------------
class MEAXQU9(MEAXQU4):
    """Mechanics - Axisymmetric - QUAD9"""

    meshType = MT.QUAD9
    elrefe = (
        ElrefeLoc(
            MT.QU9,
            gauss=("RIGI=FPG9", "FPG1=FPG1", "MASS=FPG9", "NOEU_S=NOEU_S", "NOEU=NOEU"),
            mater=("RIGI", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )


# ------------------------------------------------------------
class MEAXTR3(MEAXQU4):
    """Mechanics - Axisymmetric - TRIA3"""

    meshType = MT.TRIA3
    elrefe = (
        ElrefeLoc(
            MT.TR3,
            gauss=("RIGI=FPG1", "FPG1=FPG1", "MASS=FPG3", "NOEU_S=NOEU_S", "NOEU=NOEU"),
            mater=("RIGI", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.SE2, gauss=("RIGI=FPG2",)),
    )


# ------------------------------------------------------------
class MEAXTR6(MEAXQU4):
    """Mechanics - Axisymmetric - TRIA6"""

    meshType = MT.TRIA6
    elrefe = (
        ElrefeLoc(
            MT.TR6,
            gauss=("RIGI=FPG6", "FPG1=FPG1", "MASS=FPG6", "NOEU_S=NOEU_S", "NOEU=NOEU"),
            mater=("RIGI", "NOEU", "FPG1"),
        ),
        ElrefeLoc(MT.SE3, gauss=("RIGI=FPG4",)),
    )
