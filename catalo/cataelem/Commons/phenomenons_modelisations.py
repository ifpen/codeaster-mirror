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


from cataelem.Tools.base_objects import Phenomenon, Modelisation, objects_from_context
import cataelem.Commons.mesh_types as MT
from cataelem.Elements.elements import EL
import cataelem.Commons.attributes as AT

PhenMod = {}


############################################################
# Pour le phenomene : MECANIQUE :
############################################################
MECANIQUE = Phenomenon(code="ME")
phen = MECANIQUE

phen.add(
    "PLAQ_MITC",
    Modelisation(
        dim=(2, 3),
        code="MIT",
        attrs=(
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.COQUE, "OUI"),
            (AT.PLAQUE, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.SOUS_POINT, "OUI"),
        ),
        elements=((MT.QUAD9, EL.PLAQ_MITC_QUAD9), (MT.TRIA6, EL.PLAQ_MITC_TRIA6)),
    ),
)

phen.add(
    "2D_BARRE",
    Modelisation(
        dim=(1, 2),
        code="2DB",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.TYPMOD, "1D"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
            (AT.STRX, "NON"),
        ),
        elements=((MT.SEG2, EL.MECA_2D_BARRE),),
    ),
)

phen.add(
    "2D_DIS_T",
    Modelisation(
        dim=(-1, 2),
        code="2DT",
        attrs=((AT.TYPMOD, "0D"), (AT.EFGE, "OUI"), (AT.SIGM, "NON")),
        elements=((MT.SEG2, EL.MECA_2D_DIS_T_L), (MT.POI1, EL.MECA_2D_DIS_T_N)),
    ),
)

phen.add(
    "2D_DIS_TR",
    Modelisation(
        dim=(-1, 2),
        code="2TR",
        attrs=((AT.TYPMOD, "0D"), (AT.EFGE, "OUI"), (AT.SIGM, "NON")),
        elements=((MT.SEG2, EL.MECA_2D_DIS_TR_L), (MT.POI1, EL.MECA_2D_DIS_TR_N)),
    ),
)

phen.add(
    "2D_FLUIDE#1",
    Modelisation(
        dim=(2, 2),
        code="2FL",
        attrs=((AT.TYPMOD, "PLAN"), (AT.FLUIDE, "OUI"), (AT.FORMULATION, "U_P_PHI")),
        elements=(
            (MT.TRIA3, EL.MEFLTR3),
            (MT.QUAD4, EL.MEFLQU4),
            (MT.TRIA6, EL.MEFLTR6),
            (MT.QUAD8, EL.MEFLQU8),
            (MT.QUAD9, EL.MEFLQU9),
            (MT.SEG2, EL.MEFLSE2),
            (MT.SEG3, EL.MEFLSE3),
        ),
    ),
)

phen.add(
    "2D_FLUIDE#2",
    Modelisation(
        dim=(2, 2),
        code="2FP",
        attrs=((AT.TYPMOD, "PLAN"), (AT.FLUIDE, "OUI"), (AT.FORMULATION, "U_P")),
        elements=(
            (MT.TRIA3, EL.MEFLTR3P),
            (MT.QUAD4, EL.MEFLQU4P),
            (MT.TRIA6, EL.MEFLTR6P),
            (MT.QUAD8, EL.MEFLQU8P),
            (MT.QUAD9, EL.MEFLQU9P),
            (MT.SEG2, EL.MEFLSE2P),
            (MT.SEG3, EL.MEFLSE3P),
        ),
    ),
)

phen.add(
    "2D_FLUIDE#3",
    Modelisation(
        dim=(2, 2),
        code="2FI",
        attrs=((AT.TYPMOD, "PLAN"), (AT.FLUIDE, "OUI"), (AT.FORMULATION, "U_PSI")),
        elements=(
            (MT.TRIA3, EL.MEFLTR3PSI),
            (MT.QUAD4, EL.MEFLQU4PSI),
            (MT.TRIA6, EL.MEFLTR6PSI),
            (MT.QUAD8, EL.MEFLQU8PSI),
            (MT.QUAD9, EL.MEFLQU9PSI),
            (MT.SEG2, EL.MEFLSE2PSI),
            (MT.SEG3, EL.MEFLSE3PSI),
        ),
    ),
)

phen.add(
    "2D_FLUI_ABSO#1",
    Modelisation(
        dim=(1, 2),
        code="2FA",
        attrs=(
            (AT.TYPMOD, "PLAN"),
            (AT.FLUIDE, "OUI"),
            (AT.ABSO, "OUI"),
            (AT.FORMULATION, "U_P_PHI"),
        ),
        elements=((MT.SEG2, EL.MEFASE2), (MT.SEG3, EL.MEFASE3)),
    ),
)

phen.add(
    "AXIS_FLUI_ABSO#1",
    Modelisation(
        dim=(1, 2),
        code="AFA",
        attrs=(
            (AT.TYPMOD, "PLAN"),
            (AT.FLUIDE, "OUI"),
            (AT.ABSO, "OUI"),
            (AT.FORMULATION, "U_P_PHI"),
            (AT.AXIS, "OUI"),
        ),
        elements=((MT.SEG2, EL.MEFAXSE2), (MT.SEG3, EL.MEFAXSE3)),
    ),
)

phen.add(
    "2D_FLUI_ABSO#2",
    Modelisation(
        dim=(1, 2),
        code="2FX",
        attrs=((AT.TYPMOD, "PLAN"), (AT.FLUIDE, "OUI"), (AT.ABSO, "OUI"), (AT.FORMULATION, "U_P")),
        elements=((MT.SEG2, EL.MEFASE2UP), (MT.SEG3, EL.MEFASE3UP)),
    ),
)

phen.add(
    "AXIS_FLUI_ABSO#2",
    Modelisation(
        dim=(1, 2),
        code="AFX",
        attrs=(
            (AT.TYPMOD, "PLAN"),
            (AT.FLUIDE, "OUI"),
            (AT.ABSO, "OUI"),
            (AT.FORMULATION, "U_P"),
            (AT.AXIS, "OUI"),
        ),
        elements=((MT.SEG2, EL.MEFAXSE2UP), (MT.SEG3, EL.MEFAXSE3UP)),
    ),
)

phen.add(
    "2D_FLUI_ABSO#3",
    Modelisation(
        dim=(1, 2),
        code="2FY",
        attrs=(
            (AT.TYPMOD, "PLAN"),
            (AT.FLUIDE, "OUI"),
            (AT.ABSO, "OUI"),
            (AT.FORMULATION, "U_PSI"),
        ),
        elements=((MT.SEG2, EL.MEFASE2UPSI), (MT.SEG3, EL.MEFASE3UPSI)),
    ),
)

phen.add(
    "AXIS_FLUI_ABSO#3",
    Modelisation(
        dim=(1, 2),
        code="AFY",
        attrs=(
            (AT.TYPMOD, "PLAN"),
            (AT.FLUIDE, "OUI"),
            (AT.ABSO, "OUI"),
            (AT.FORMULATION, "U_PSI"),
            (AT.AXIS, "OUI"),
        ),
        elements=((MT.SEG2, EL.MEFAXSE2UPSI), (MT.SEG3, EL.MEFAXSE3UPSI)),
    ),
)

phen.add(
    "2D_FLUI_PESA#1",
    Modelisation(
        dim=(2, 2),
        code="2FP",
        attrs=(
            (AT.TYPMOD, "PLAN"),
            (AT.FLUIDE, "OUI"),
            (AT.PESA, "OUI"),
            (AT.FORMULATION, "U_P_PHI"),
        ),
        elements=(
            (MT.TRIA3, EL.MEFP_FACE3),
            (MT.QUAD4, EL.MEFP_FACE4),
            (MT.TRIA6, EL.MEFP_FACE6),
            (MT.QUAD8, EL.MEFP_FACE8),
            (MT.QUAD9, EL.MEFP_FACE9),
        ),
    ),
)

phen.add(
    "2D_FLUI_STRU#1",
    Modelisation(
        dim=(1, 2),
        code="FS2",
        attrs=(
            (AT.TYPMOD, "PLAN"),
            (AT.FLUIDE, "OUI"),
            (AT.FSI, "OUI"),
            (AT.FORMULATION, "U_P_PHI"),
        ),
        elements=((MT.SEG2, EL.MEFSSE2), (MT.SEG3, EL.MEFSSE3)),
    ),
)

phen.add(
    "2D_FLUI_STRU#2",
    Modelisation(
        dim=(1, 2),
        code="FP2",
        attrs=((AT.TYPMOD, "PLAN"), (AT.FLUIDE, "OUI"), (AT.FSI, "OUI"), (AT.FORMULATION, "U_P")),
        elements=((MT.SEG2, EL.MEFSSE2P), (MT.SEG3, EL.MEFSSE3P)),
    ),
)

phen.add(
    "2D_FLUI_STRU#3",
    Modelisation(
        dim=(1, 2),
        code="FI2",
        attrs=((AT.TYPMOD, "PLAN"), (AT.FLUIDE, "OUI"), (AT.FSI, "OUI"), (AT.FORMULATION, "U_PSI")),
        elements=((MT.SEG2, EL.MEFSSE2PSI), (MT.SEG3, EL.MEFSSE3PSI)),
    ),
)

phen.add(
    "3D",
    Modelisation(
        dim=(3, 3),
        code="3D_",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D")),
        elements=(
            (MT.HEXA8, EL.MECA_HEXA8),
            (MT.PENTA6, EL.MECA_PENTA6),
            (MT.TETRA4, EL.MECA_TETRA4),
            (MT.QUAD4, EL.MECA_FACE4),
            (MT.TRIA3, EL.MECA_FACE3),
            (MT.SEG2, EL.MECA_ARETE2),
            (MT.HEXA27, EL.MECA_HEXA27),
            (MT.HEXA20, EL.MECA_HEXA20),
            (MT.PENTA15, EL.MECA_PENTA15),
            (MT.PENTA18, EL.MECA_PENTA18),
            (MT.TETRA10, EL.MECA_TETRA10),
            (MT.TETRA20, EL.MECA_TETRA20),
            (MT.QUAD9, EL.MECA_FACE9),
            (MT.QUAD8, EL.MECA_FACE8),
            (MT.TRIA6, EL.MECA_FACE6),
            (MT.TRIA10, EL.MECA_FACE10),
            (MT.SEG3, EL.MECA_ARETE3),
            (MT.SEG4, EL.MECA_ARETE4),
            (MT.PYRAM5, EL.MECA_PYRAM5),
            (MT.PYRAM13, EL.MECA_PYRAM13),
        ),
    ),
)

phen.add(
    "3D1XH",
    Modelisation(
        dim=(3, 3),
        code="3X1",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.LXFEM, "OUI"), (AT.XFEM, "XH")),
        elements=(
            (MT.HEXA20, EL.MECA_XH_HEXA20),
            (MT.HEXA8, EL.MECA_XH_HEXA8),
            (MT.PENTA15, EL.MECA_XH_PENTA15),
            (MT.PENTA6, EL.MECA_XH_PENTA6),
            (MT.PYRAM13, EL.MECA_XH_PYRAM13),
            (MT.PYRAM5, EL.MECA_XH_PYRAM5),
            (MT.TETRA10, EL.MECA_XH_TETRA10),
            (MT.TETRA4, EL.MECA_XH_TETRA4),
            (MT.TRIA3, EL.MECA_XH_FACE3),
            (MT.TRIA6, EL.MECA_XH_FACE6),
            (MT.QUAD4, EL.MECA_XH_FACE4),
            (MT.QUAD8, EL.MECA_XH_FACE8),
        ),
    ),
)

phen.add(
    "3D1XHT",
    Modelisation(
        dim=(3, 3),
        code="3X3",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.LXFEM, "OUI"), (AT.XFEM, "XHT")),
        elements=(
            (MT.HEXA8, EL.MECA_XHT_HEXA8),
            (MT.PENTA6, EL.MECA_XHT_PENTA6),
            (MT.PYRAM5, EL.MECA_XHT_PYRAM5),
            (MT.TETRA4, EL.MECA_XHT_TETRA4),
            (MT.QUAD4, EL.MECA_XHT_FACE4),
            (MT.TRIA3, EL.MECA_XHT_FACE3),
        ),
    ),
)

phen.add(
    "3D1XT",
    Modelisation(
        dim=(3, 3),
        code="3X2",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.LXFEM, "OUI"), (AT.XFEM, "XT")),
        elements=(
            (MT.HEXA8, EL.MECA_XT_HEXA8),
            (MT.PENTA6, EL.MECA_XT_PENTA6),
            (MT.PYRAM5, EL.MECA_XT_PYRAM5),
            (MT.TETRA4, EL.MECA_XT_TETRA4),
            (MT.QUAD4, EL.MECA_XT_FACE4),
            (MT.TRIA3, EL.MECA_XT_FACE3),
        ),
    ),
)

phen.add(
    "3D2XHT",
    Modelisation(
        dim=(3, 3),
        code="3X6",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.LXFEM, "OUI"), (AT.XFEM, "XHT")),
        elements=(
            (MT.HEXA20, EL.MECA_XHT_HEXA20),
            (MT.PENTA15, EL.MECA_XHT_PENTA15),
            (MT.PYRAM13, EL.MECA_XHT_PYRAM13),
            (MT.TETRA10, EL.MECA_XHT_TETRA10),
            (MT.QUAD8, EL.MECA_XHT_FACE8),
            (MT.TRIA6, EL.MECA_XHT_FACE6),
        ),
    ),
)

phen.add(
    "3D2XT",
    Modelisation(
        dim=(3, 3),
        code="3X5",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.LXFEM, "OUI"), (AT.XFEM, "XT")),
        elements=(
            (MT.HEXA20, EL.MECA_XT_HEXA20),
            (MT.PENTA15, EL.MECA_XT_PENTA15),
            (MT.PYRAM13, EL.MECA_XT_PYRAM13),
            (MT.TETRA10, EL.MECA_XT_TETRA10),
            (MT.QUAD8, EL.MECA_XT_FACE8),
            (MT.TRIA6, EL.MECA_XT_FACE6),
        ),
    ),
)

phen.add(
    "3DXH1",
    Modelisation(
        dim=(3, 3),
        code="3XA",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.LXFEM, "OUI"), (AT.XFEM, "XH1")),
        elements=(
            (MT.HEXA8, EL.MECA_XH1_HEXA8),
            (MT.PENTA6, EL.MECA_XH1_PENTA6),
            (MT.PYRAM5, EL.MECA_XH1_PYRAM5),
            (MT.TETRA4, EL.MECA_XH1_TETRA4),
            (MT.QUAD4, EL.MECA_XH1_FACE4),
            (MT.TRIA3, EL.MECA_XH1_FACE3),
        ),
    ),
)

phen.add(
    "3DXH2",
    Modelisation(
        dim=(3, 3),
        code="3XB",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.LXFEM, "OUI"), (AT.XFEM, "XH2")),
        elements=(
            (MT.HEXA8, EL.MECA_XH2_HEXA8),
            (MT.PENTA6, EL.MECA_XH2_PENTA6),
            (MT.PYRAM5, EL.MECA_XH2_PYRAM5),
            (MT.TETRA4, EL.MECA_XH2_TETRA4),
            (MT.QUAD4, EL.MECA_XH2_FACE4),
            (MT.TRIA3, EL.MECA_XH2_FACE3),
        ),
    ),
)

phen.add(
    "3DXH3",
    Modelisation(
        dim=(3, 3),
        code="3XC",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.LXFEM, "OUI"), (AT.XFEM, "XH3")),
        elements=(
            (MT.HEXA8, EL.MECA_XH3_HEXA8),
            (MT.PENTA6, EL.MECA_XH3_PENTA6),
            (MT.PYRAM5, EL.MECA_XH3_PYRAM5),
            (MT.TETRA4, EL.MECA_XH3_TETRA4),
            (MT.QUAD4, EL.MECA_XH3_FACE4),
            (MT.TRIA3, EL.MECA_XH3_FACE3),
        ),
    ),
)

phen.add(
    "3DXH4",
    Modelisation(
        dim=(3, 3),
        code="3XD",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.LXFEM, "OUI"), (AT.XFEM, "XH4")),
        elements=(
            (MT.HEXA8, EL.MECA_XH4_HEXA8),
            (MT.PENTA6, EL.MECA_XH4_PENTA6),
            (MT.PYRAM5, EL.MECA_XH4_PYRAM5),
            (MT.TETRA4, EL.MECA_XH4_TETRA4),
            (MT.QUAD4, EL.MECA_XH4_FACE4),
            (MT.TRIA3, EL.MECA_XH4_FACE3),
        ),
    ),
)

phen.add(
    "3D_ABSO",
    Modelisation(
        dim=(2, 3),
        code="3DA",
        attrs=((AT.TYPMOD, "3D"), (AT.FLUIDE, "NON"), (AT.ABSO, "OUI")),
        elements=(
            (MT.TRIA3, EL.MEAB_FACE3),
            (MT.QUAD4, EL.MEAB_FACE4),
            (MT.TRIA6, EL.MEAB_FACE6),
            (MT.QUAD8, EL.MEAB_FACE8),
            (MT.QUAD9, EL.MEAB_FACE9),
        ),
    ),
)

phen.add(
    "3D_DIL#1",
    Modelisation(
        dim=(3, 3),
        code="D3D",
        attrs=((AT.TYPMOD, "3D"), (AT.NBSIGM, "6"), (AT.FORMULATION, "DIL")),
        elements=(
            (MT.TETRA10, EL.T10_3D),
            (MT.PENTA15, EL.P15_3D),
            (MT.HEXA20, EL.H20_3D),
            (MT.QUAD8, EL.MECA_FACE8),
            (MT.TRIA6, EL.MECA_FACE6),
            (MT.SEG3, EL.MECA_ARETE3),
        ),
    ),
)

phen.add(
    "3D_DIL#2",
    Modelisation(
        dim=(3, 3),
        code="D3I",
        attrs=((AT.TYPMOD, "3D"), (AT.NBSIGM, "6"), (AT.FORMULATION, "DIL_INCO")),
        elements=(
            (MT.TETRA10, EL.T10_3DI),
            (MT.PENTA15, EL.P15_3DI),
            (MT.HEXA20, EL.H20_3DI),
            (MT.QUAD8, EL.MECA_FACE8),
            (MT.TRIA6, EL.MECA_FACE6),
            (MT.SEG3, EL.MECA_ARETE3),
        ),
    ),
)

phen.add(
    "3D_FAISCEAU",
    Modelisation(
        dim=(3, 3),
        code="3DF",
        elements=((MT.HEXA8, EL.MECA_POHO_HEXA8), (MT.HEXA20, EL.MECA_POHO_HEXA20)),
    ),
)

phen.add(
    "3D_FLUIDE#1",
    Modelisation(
        dim=(3, 3),
        code="3FL",
        attrs=((AT.TYPMOD, "3D"), (AT.FLUIDE, "OUI"), (AT.FORMULATION, "U_P_PHI")),
        elements=(
            (MT.HEXA8, EL.MEFL_HEXA8),
            (MT.HEXA20, EL.MEFL_HEXA20),
            (MT.HEXA27, EL.MEFL_HEXA27),
            (MT.PENTA6, EL.MEFL_PENTA6),
            (MT.PENTA15, EL.MEFL_PENTA15),
            (MT.PYRAM5, EL.MEFL_PYRAM5),
            (MT.PYRAM13, EL.MEFL_PYRAM13),
            (MT.TETRA4, EL.MEFL_TETRA4),
            (MT.TETRA10, EL.MEFL_TETRA10),
            (MT.TRIA3, EL.MEFL_FACE3),
            (MT.QUAD4, EL.MEFL_FACE4),
            (MT.TRIA6, EL.MEFL_FACE6),
            (MT.QUAD8, EL.MEFL_FACE8),
            (MT.QUAD9, EL.MEFL_FACE9),
        ),
    ),
)

phen.add(
    "3D_FLUIDE#2",
    Modelisation(
        dim=(3, 3),
        code="3FP",
        attrs=((AT.TYPMOD, "3D"), (AT.FLUIDE, "OUI"), (AT.FORMULATION, "U_P")),
        elements=(
            (MT.HEXA8, EL.MEFL_HEXA8P),
            (MT.HEXA20, EL.MEFL_HEXA20P),
            (MT.HEXA27, EL.MEFL_HEXA27P),
            (MT.PENTA6, EL.MEFL_PENTA6P),
            (MT.PENTA15, EL.MEFL_PENTA15P),
            (MT.PYRAM5, EL.MEFL_PYRAM5P),
            (MT.PYRAM13, EL.MEFL_PYRAM13P),
            (MT.TETRA4, EL.MEFL_TETRA4P),
            (MT.TETRA10, EL.MEFL_TETRA10P),
            (MT.TRIA3, EL.MEFL_FACE3P),
            (MT.QUAD4, EL.MEFL_FACE4P),
            (MT.TRIA6, EL.MEFL_FACE6P),
            (MT.QUAD8, EL.MEFL_FACE8P),
            (MT.QUAD9, EL.MEFL_FACE9P),
        ),
    ),
)

phen.add(
    "3D_FLUIDE#3",
    Modelisation(
        dim=(3, 3),
        code="3FI",
        attrs=((AT.TYPMOD, "3D"), (AT.FLUIDE, "OUI"), (AT.FORMULATION, "U_PSI")),
        elements=(
            (MT.HEXA8, EL.MEFL_HEXA8PSI),
            (MT.HEXA20, EL.MEFL_HEXA20PSI),
            (MT.HEXA27, EL.MEFL_HEXA27PSI),
            (MT.PENTA6, EL.MEFL_PENTA6PSI),
            (MT.PENTA15, EL.MEFL_PENTA15PSI),
            (MT.PYRAM5, EL.MEFL_PYRAM5PSI),
            (MT.PYRAM13, EL.MEFL_PYRAM13PSI),
            (MT.TETRA4, EL.MEFL_TETRA4PSI),
            (MT.TETRA10, EL.MEFL_TETRA10PSI),
            (MT.TRIA3, EL.MEFL_FACE3PSI),
            (MT.QUAD4, EL.MEFL_FACE4PSI),
            (MT.TRIA6, EL.MEFL_FACE6PSI),
            (MT.QUAD8, EL.MEFL_FACE8PSI),
            (MT.QUAD9, EL.MEFL_FACE9PSI),
        ),
    ),
)

phen.add(
    "3D_FLUI_ABSO#1",
    Modelisation(
        dim=(2, 3),
        code="3FA",
        attrs=((AT.FLUIDE, "OUI"), (AT.ABSO, "OUI"), (AT.FORMULATION, "U_P_PHI")),
        elements=(
            (MT.TRIA3, EL.MEFA_FACE3),
            (MT.QUAD4, EL.MEFA_FACE4),
            (MT.TRIA6, EL.MEFA_FACE6),
            (MT.QUAD8, EL.MEFA_FACE8),
            (MT.QUAD9, EL.MEFA_FACE9),
        ),
    ),
)

phen.add(
    "3D_FLUI_ABSO#2",
    Modelisation(
        dim=(2, 3),
        code="3FX",
        attrs=((AT.FLUIDE, "OUI"), (AT.ABSO, "OUI"), (AT.FORMULATION, "U_P")),
        elements=(
            (MT.TRIA3, EL.MEFA_FACE3UP),
            (MT.QUAD4, EL.MEFA_FACE4UP),
            (MT.TRIA6, EL.MEFA_FACE6UP),
            (MT.QUAD8, EL.MEFA_FACE8UP),
            (MT.QUAD9, EL.MEFA_FACE9UP),
        ),
    ),
)

phen.add(
    "3D_FLUI_ABSO#3",
    Modelisation(
        dim=(2, 3),
        code="3FY",
        attrs=((AT.FLUIDE, "OUI"), (AT.ABSO, "OUI"), (AT.FORMULATION, "U_PSI")),
        elements=(
            (MT.TRIA3, EL.MEFA_FACE3UPSI),
            (MT.QUAD4, EL.MEFA_FACE4UPSI),
            (MT.TRIA6, EL.MEFA_FACE6UPSI),
            (MT.QUAD8, EL.MEFA_FACE8UPSI),
            (MT.QUAD9, EL.MEFA_FACE9UPSI),
        ),
    ),
)

phen.add(
    "3D_GRAD_VARI",
    Modelisation(
        dim=(3, 3),
        code="3DV",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.TYPMOD2, "GRADVARI")),
        elements=(
            (MT.HEXA20, EL.MVCA_HEXA20),
            (MT.PENTA15, EL.MVCA_PENTA15),
            (MT.PYRAM13, EL.MVCA_PYRAM13),
            (MT.TETRA10, EL.MVCA_TETRA10),
            (MT.QUAD8, EL.MECA_FACE8),
            (MT.TRIA6, EL.MECA_FACE6),
            (MT.SEG3, EL.MECA_ARETE3),
        ),
    ),
)

phen.add(
    "3D_GRAD_INCO",
    Modelisation(
        dim=(3, 3),
        code="3GI",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.INCO, "C5GV"), (AT.TYPMOD2, "GRADVARI")),
        elements=(
            (MT.HEXA20, EL.GVI_3D_HE20),
            (MT.PENTA15, EL.GVI_3D_PE15),
            (MT.PYRAM13, EL.GVI_3D_PY13),
            (MT.TETRA10, EL.GVI_3D_TE10),
            (MT.QUAD8, EL.MECA_FACE8),
            (MT.TRIA6, EL.MECA_FACE6),
            (MT.SEG3, EL.MECA_ARETE3),
        ),
    ),
)

phen.add(
    "3D_HH2D",
    Modelisation(
        dim=(3, 3),
        code="3Z4",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.TETRA10, EL.HH2_TETRA10D),
            (MT.PENTA15, EL.HH2_PENTA15D),
            (MT.HEXA20, EL.HH2_HEXA20D),
            (MT.QUAD8, EL.HH2_FACE8),
            (MT.TRIA6, EL.HH2_FACE6),
        ),
    ),
)

phen.add(
    "3D_HH2MD",
    Modelisation(
        dim=(3, 3),
        code="3A3",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.TETRA10, EL.HH2M_TETRA10D),
            (MT.PENTA15, EL.HH2M_PENTA15D),
            (MT.HEXA20, EL.HH2M_HEXA20D),
            (MT.QUAD8, EL.HH2M_FACE8),
            (MT.TRIA6, EL.HH2M_FACE6),
        ),
    ),
)

phen.add(
    "3D_HH2MS",
    Modelisation(
        dim=(3, 3),
        code="3R9",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.HH2M_TETRA10S),
            (MT.PENTA15, EL.HH2M_PENTA15S),
            (MT.HEXA20, EL.HH2M_HEXA20S),
            (MT.QUAD8, EL.HH2M_FACE8),
            (MT.TRIA6, EL.HH2M_FACE6),
        ),
    ),
)

phen.add(
    "3D_HH2MS_DIL",
    Modelisation(
        dim=(3, 3),
        code="3D9",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.INTTHM, "RED"),
            (AT.DIL, "OUI"),
            (AT.FORMULATION, "DIL"),
        ),
        elements=(
            (MT.TETRA10, EL.HH2M_TETRA10S_DI),
            (MT.PENTA15, EL.HH2M_PENTA15S_DI),
            (MT.HEXA20, EL.HH2M_HEXA20S_DIL),
            (MT.QUAD8, EL.HH2M_FACE8),
            (MT.TRIA6, EL.HH2M_FACE6),
        ),
    ),
)

phen.add(
    "3D_HH2M_SI",
    Modelisation(
        dim=(3, 3),
        code="3M1",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
        ),
        elements=(
            (MT.TETRA10, EL.HH2M_TETRA10M),
            (MT.PENTA15, EL.HH2M_PENTA15M),
            (MT.HEXA20, EL.HH2M_HEXA20M),
            (MT.QUAD8, EL.HH2M_FACE8),
            (MT.TRIA6, EL.HH2M_FACE6),
        ),
    ),
)

phen.add(
    "3D_HH2S",
    Modelisation(
        dim=(3, 3),
        code="3Z3",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.HH2_TETRA10S),
            (MT.PENTA15, EL.HH2_PENTA15S),
            (MT.HEXA20, EL.HH2_HEXA20S),
            (MT.QUAD8, EL.HH2_FACE8),
            (MT.TRIA6, EL.HH2_FACE6),
        ),
    ),
)

phen.add(
    "3D_HH2SUDA",
    Modelisation(
        dim=(3, 3),
        code="3AD",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.TYPMOD3, "SUSHI"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
        ),
        elements=((MT.HEXA27, EL.ZHH2_HEXA27_SUDA), (MT.QUAD9, EL.ZHH2_FACE9_SUDA)),
    ),
)

phen.add(
    "3D_HHD",
    Modelisation(
        dim=(3, 3),
        code="3Z2",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.TETRA10, EL.HH_TETRA10D),
            (MT.PENTA15, EL.HH_PENTA15D),
            (MT.HEXA20, EL.HH_HEXA20D),
            (MT.QUAD8, EL.HH_FACE8),
            (MT.TRIA6, EL.HH_FACE6),
        ),
    ),
)

phen.add(
    "3D_HHM",
    Modelisation(
        dim=(3, 3),
        code="3H1",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
        ),
        elements=(
            (MT.TETRA10, EL.HHM_TETRA10),
            (MT.PENTA15, EL.HHM_PENTA15),
            (MT.HEXA20, EL.HHM_HEXA20),
            (MT.QUAD8, EL.HHM_FACE8),
            (MT.TRIA6, EL.HHM_FACE6),
        ),
    ),
)

phen.add(
    "3D_HHMD",
    Modelisation(
        dim=(3, 3),
        code="3H6",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.TETRA10, EL.HHM_TETRA10D),
            (MT.PENTA15, EL.HHM_PENTA15D),
            (MT.HEXA20, EL.HHM_HEXA20D),
            (MT.QUAD8, EL.HHM_FACE8),
            (MT.TRIA6, EL.HHM_FACE6),
        ),
    ),
)

phen.add(
    "3D_HHMS",
    Modelisation(
        dim=(3, 3),
        code="3R1",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.HHM_TETRA10S),
            (MT.PENTA15, EL.HHM_PENTA15S),
            (MT.HEXA20, EL.HHM_HEXA20S),
            (MT.QUAD8, EL.HHM_FACE8),
            (MT.TRIA6, EL.HHM_FACE6),
        ),
    ),
)

phen.add(
    "3D_HHS",
    Modelisation(
        dim=(3, 3),
        code="3Z1",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.HH_TETRA10S),
            (MT.PENTA15, EL.HH_PENTA15S),
            (MT.HEXA20, EL.HH_HEXA20S),
            (MT.QUAD8, EL.HH_FACE8),
            (MT.TRIA6, EL.HH_FACE6),
        ),
    ),
)

phen.add(
    "3D_HM",
    Modelisation(
        dim=(3, 3),
        code="3H2",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10),
            (MT.PYRAM13, EL.HM_PYRAM13),
            (MT.PENTA15, EL.HM_PENTA15),
            (MT.HEXA20, EL.HM_HEXA20),
            (MT.QUAD8, EL.HM_FACE8),
            (MT.TRIA6, EL.HM_FACE6),
        ),
    ),
)

phen.add(
    "3D_HMD",
    Modelisation(
        dim=(3, 3),
        code="3H7",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10D),
            (MT.PYRAM13, EL.HM_PYRAM13D),
            (MT.PENTA15, EL.HM_PENTA15D),
            (MT.HEXA20, EL.HM_HEXA20D),
            (MT.QUAD8, EL.HM_FACE8),
            (MT.TRIA6, EL.HM_FACE6),
        ),
    ),
)

phen.add(
    "3D_HMS",
    Modelisation(
        dim=(3, 3),
        code="3R2",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10S),
            (MT.PYRAM13, EL.HM_PYRAM13S),
            (MT.PENTA15, EL.HM_PENTA15S),
            (MT.HEXA20, EL.HM_HEXA20S),
            (MT.QUAD8, EL.HM_FACE8),
            (MT.TRIA6, EL.HM_FACE6),
        ),
    ),
)

phen.add(
    "3D_HMS_DIL",
    Modelisation(
        dim=(3, 3),
        code="3D8",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.INTTHM, "RED"),
            (AT.DIL, "OUI"),
            (AT.FORMULATION, "DIL"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10S_DIL),
            (MT.PYRAM13, EL.HM_PYRAM13S_DIL),
            (MT.PENTA15, EL.HM_PENTA15S_DIL),
            (MT.HEXA20, EL.HM_HEXA20S_DIL),
            (MT.QUAD8, EL.HM_FACE8),
            (MT.TRIA6, EL.HM_FACE6),
        ),
    ),
)

phen.add(
    "3D_HM_SI",
    Modelisation(
        dim=(3, 3),
        code="3M2",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10M),
            (MT.PYRAM13, EL.HM_PYRAM13M),
            (MT.PENTA15, EL.HM_PENTA15M),
            (MT.HEXA20, EL.HM_HEXA20M),
            (MT.QUAD8, EL.HM_FACE8),
            (MT.TRIA6, EL.HM_FACE6),
        ),
    ),
)

phen.add(
    "3D_HM_SI_DIL",
    Modelisation(
        dim=(3, 3),
        code="3D7",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.DIL, "OUI"),
            (AT.FORMULATION, "DIL"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10M_DIL),
            (MT.PYRAM13, EL.HM_PYRAM13M_DIL),
            (MT.PENTA15, EL.HM_PENTA15M_DIL),
            (MT.HEXA20, EL.HM_HEXA20M_DIL),
            (MT.QUAD8, EL.HM_FACE8),
            (MT.TRIA6, EL.HM_FACE6),
        ),
    ),
)

phen.add(
    "3D_HM_XH",
    Modelisation(
        dim=(3, 3),
        code="3XL",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10_XH),
            (MT.PYRAM13, EL.HM_PYRAM13_XH),
            (MT.PENTA15, EL.HM_PENTA15_XH),
            (MT.HEXA20, EL.HM_HEXA20_XH),
            (MT.QUAD8, EL.HM_FACE8_XH),
            (MT.TRIA6, EL.HM_FACE6_XH),
        ),
    ),
)

phen.add(
    "3D_HM_XH1",
    Modelisation(
        dim=(3, 3),
        code="3XM",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH1"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10_XH1),
            (MT.PYRAM13, EL.HM_PYRAM13_XH1),
            (MT.PENTA15, EL.HM_PENTA15_XH1),
            (MT.HEXA20, EL.HM_HEXA20_XH1),
            (MT.QUAD8, EL.HM_FACE8_XH1),
            (MT.TRIA6, EL.HM_FACE6_XH1),
        ),
    ),
)

phen.add(
    "3D_HM_XH2",
    Modelisation(
        dim=(3, 3),
        code="3XN",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH2"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10_XH2),
            (MT.PYRAM13, EL.HM_PYRAM13_XH2),
            (MT.PENTA15, EL.HM_PENTA15_XH2),
            (MT.HEXA20, EL.HM_HEXA20_XH2),
            (MT.QUAD8, EL.HM_FACE8_XH2),
            (MT.TRIA6, EL.HM_FACE6_XH2),
        ),
    ),
)

phen.add(
    "3D_HM_XH3",
    Modelisation(
        dim=(3, 3),
        code="3XO",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH3"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10_XH3),
            (MT.PYRAM13, EL.HM_PYRAM13_XH3),
            (MT.PENTA15, EL.HM_PENTA15_XH3),
            (MT.HEXA20, EL.HM_HEXA20_XH3),
            (MT.QUAD8, EL.HM_FACE8_XH3),
            (MT.TRIA6, EL.HM_FACE6_XH3),
        ),
    ),
)

phen.add(
    "3D_HM_XH_D",
    Modelisation(
        dim=(3, 3),
        code="3XJ",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10D_XH),
            (MT.PENTA15, EL.HM_PENTA15D_XH),
            (MT.HEXA20, EL.HM_HEXA20D_XH),
            (MT.QUAD8, EL.HM_FACE8_XH),
            (MT.TRIA6, EL.HM_FACE6_XH),
        ),
    ),
)

phen.add(
    "3D_HM_XH_S",
    Modelisation(
        dim=(3, 3),
        code="3XK",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10S_XH),
            (MT.PENTA15, EL.HM_PENTA15S_XH),
            (MT.HEXA20, EL.HM_HEXA20S_XH),
            (MT.QUAD8, EL.HM_FACE8_XH),
            (MT.TRIA6, EL.HM_FACE6_XH),
        ),
    ),
)

phen.add(
    "3D_HM_XH_SI",
    Modelisation(
        dim=(3, 3),
        code="3XI",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=(
            (MT.TETRA10, EL.HM_TETRA10M_XH),
            (MT.PENTA15, EL.HM_PENTA15M_XH),
            (MT.HEXA20, EL.HM_HEXA20M_XH),
            (MT.QUAD8, EL.HM_FACE8_XH),
            (MT.TRIA6, EL.HM_FACE6_XH),
        ),
    ),
)

phen.add(
    "3D_HS",
    Modelisation(
        dim=(3, 3),
        code="EH3",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.H_TETRA10S),
            (MT.PENTA15, EL.H_PENTA15S),
            (MT.HEXA20, EL.H_HEXA20S),
            (MT.QUAD8, EL.H_FACE8),
            (MT.TRIA6, EL.H_FACE6),
        ),
    ),
)

phen.add(
    "3D_INCO_UP",
    Modelisation(
        dim=(3, 3),
        code="3UP",
        attrs=((AT.NBSIGM, "6"), (AT.INCO, "C2"), (AT.TYPMOD, "3D")),
        elements=(
            (MT.HEXA20, EL.MIUP_HEXA20),
            (MT.PENTA15, EL.MIUP_PENTA15),
            (MT.TETRA10, EL.MIUP_TETRA10),
            (MT.TETRA4, EL.MIUP_TETRA4),
            (MT.QUAD8, EL.MECA_FACE8),
            (MT.TRIA6, EL.MECA_FACE6),
            (MT.TRIA3, EL.MECA_FACE3),
            (MT.SEG3, EL.MECA_ARETE3),
            (MT.SEG2, EL.MECA_ARETE2),
        ),
    ),
)

phen.add(
    "3D_INCO_UPG",
    Modelisation(
        dim=(3, 3),
        code="3DI",
        attrs=((AT.NBSIGM, "6"), (AT.INCO, "C3"), (AT.TYPMOD, "3D")),
        elements=(
            (MT.HEXA20, EL.MINC_HEXA20),
            (MT.PENTA15, EL.MINC_PENTA15),
            (MT.TETRA10, EL.MINC_TETRA10),
            (MT.QUAD8, EL.MECA_FACE8),
            (MT.TRIA6, EL.MECA_FACE6),
            (MT.SEG3, EL.MECA_ARETE3),
        ),
    ),
)

phen.add(
    "3D_INCO_UPO",
    Modelisation(
        dim=(3, 3),
        code="3OS",
        attrs=((AT.NBSIGM, "6"), (AT.INCO, "C2O"), (AT.TYPMOD, "3D")),
        elements=(
            (MT.HEXA8, EL.MINCOS_HEXA8),
            (MT.PENTA6, EL.MINCOS_PENTA6),
            (MT.PYRAM5, EL.MINCOS_PYRAM5),
            (MT.TETRA4, EL.MINCOS_TETRA4),
            (MT.QUAD4, EL.MECA_FACE4),
            (MT.TRIA3, EL.MECA_FACE3),
            (MT.SEG3, EL.MECA_ARETE3),
            (MT.SEG2, EL.MECA_ARETE2),
        ),
    ),
)

phen.add(
    "3D_INTERFACE",
    Modelisation(
        dim=(3, 3),
        code="3EI",
        attrs=((AT.TYPMOD, "3D"), (AT.TYPMOD2, "INTERFAC"), (AT.INTERFACE, "OUI")),
        elements=((MT.HEXA20, EL.MEEI_HEXA20), (MT.PENTA15, EL.MEEI_PENTA15)),
    ),
)

phen.add(
    "3D_INTERFACE_S",
    Modelisation(
        dim=(3, 3),
        code="3IS",
        attrs=((AT.TYPMOD, "3D"), (AT.TYPMOD2, "INTERFAC"), (AT.INTERFACE, "OUI")),
        elements=((MT.HEXA20, EL.MEEI_HEXS20), (MT.PENTA15, EL.MEEI_PENTS15)),
    ),
)

phen.add(
    "3D_INTERF_POU",
    Modelisation(
        dim=(3, 3),
        code="3SP",
        attrs=((AT.TYPMOD, "3D"), (AT.TYPMOD2, "INTSOLPI"), (AT.INTERFACE, "OUI")),
        elements=((MT.HEXA27, EL.MECA_S8P2),),
    ),
)

phen.add(
    "3D_JOINT",
    Modelisation(
        dim=(3, 3),
        code="3FI",
        attrs=((AT.TYPMOD, "3D"), (AT.TYPMOD2, "ELEMJOIN"), (AT.INTERFACE, "OUI")),
        elements=(
            (MT.HEXA8, EL.MEFI_HEXA8),
            (MT.PENTA6, EL.MEFI_PENTA6),
            (MT.HEXA20, EL.MEFI_HEXA20),
            (MT.PENTA15, EL.MEFI_PENTA15),
        ),
    ),
)

phen.add(
    "3D_JOINT_HYME",
    Modelisation(
        dim=(3, 3),
        code="3FH",
        attrs=((AT.TYPMOD, "3D"), (AT.TYPMOD2, "EJ_HYME"), (AT.INTERFACE, "OUI")),
        elements=((MT.HEXA20, EL.EJHYME_HEXA20), (MT.PENTA15, EL.EJHYME_PENTA15)),
    ),
)

phen.add(
    "3D_MIX_STA#1",
    Modelisation(
        dim=(3, 3),
        code="3MV",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.FORMULATION, "STA")),
        elements=(
            (MT.HEXA8, EL.MEMS_HEXA8),
            (MT.PENTA6, EL.MEMS_PENTA6),
            (MT.TETRA4, EL.MEMS_TETRA4),
            (MT.PYRAM5, EL.MEMS_PYRAM5),
            (MT.QUAD4, EL.MECA_FACE4),
            (MT.TRIA3, EL.MECA_FACE3),
            (MT.SEG2, EL.MECA_ARETE2),
        ),
    ),
)

phen.add(
    "3D_MIX_STA#2",
    Modelisation(
        dim=(3, 3),
        code="3SI",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D"), (AT.FORMULATION, "STA_INCO")),
        elements=((MT.HEXA8, EL.MEMI_HEXA8), (MT.QUAD4, EL.MECA_FACE4), (MT.SEG2, EL.MECA_ARETE2)),
    ),
)

phen.add(
    "3D_SI",
    Modelisation(
        dim=(3, 3),
        code="3DS",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D")),
        elements=(
            (MT.HEXA20, EL.MECA_HEXS20),
            (MT.HEXA8, EL.MECA_HEXS8),
            (MT.QUAD8, EL.MECA_FACE8),
            (MT.QUAD4, EL.MECA_FACE4),
            (MT.TETRA10, EL.MECA_TETRS10),
            (MT.TRIA6, EL.MECA_FACE6),
        ),
    ),
)

phen.add(
    "3D_THH2D",
    Modelisation(
        dim=(3, 3),
        code="3A1",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.TETRA10, EL.THH2_TETRA10D),
            (MT.PENTA15, EL.THH2_PENTA15D),
            (MT.HEXA20, EL.THH2_HEXA20D),
            (MT.QUAD8, EL.THH2_FACE8),
            (MT.TRIA6, EL.THH2_FACE6),
        ),
    ),
)

phen.add(
    "3D_THH2MD",
    Modelisation(
        dim=(3, 3),
        code="3A2",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.TETRA10, EL.THH2M_TETRA10D),
            (MT.PENTA15, EL.THH2M_PENTA15D),
            (MT.HEXA20, EL.THH2M_HEXA20D),
            (MT.QUAD8, EL.THH2M_FACE8),
            (MT.TRIA6, EL.THH2M_FACE6),
        ),
    ),
)

phen.add(
    "3D_THH2MS",
    Modelisation(
        dim=(3, 3),
        code="3R8",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.THH2M_TETRA10S),
            (MT.PENTA15, EL.THH2M_PENTA15S),
            (MT.HEXA20, EL.THH2M_HEXA20S),
            (MT.QUAD8, EL.THH2M_FACE8),
            (MT.TRIA6, EL.THH2M_FACE6),
        ),
    ),
)

phen.add(
    "3D_THH2S",
    Modelisation(
        dim=(3, 3),
        code="3R7",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.THH2_TETRA10S),
            (MT.PENTA15, EL.THH2_PENTA15S),
            (MT.HEXA20, EL.THH2_HEXA20S),
            (MT.QUAD8, EL.THH2_FACE8),
            (MT.TRIA6, EL.THH2_FACE6),
        ),
    ),
)

phen.add(
    "3D_THHD",
    Modelisation(
        dim=(3, 3),
        code="3H8",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.TETRA10, EL.THH_TETRA10D),
            (MT.PENTA15, EL.THH_PENTA15D),
            (MT.HEXA20, EL.THH_HEXA20D),
            (MT.QUAD8, EL.THH_FACE8),
            (MT.TRIA6, EL.THH_FACE6),
        ),
    ),
)

phen.add(
    "3D_THHM",
    Modelisation(
        dim=(3, 3),
        code="3H4",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
        ),
        elements=(
            (MT.HEXA20, EL.THHM_HEXA20),
            (MT.TETRA10, EL.THHM_TETRA10),
            (MT.PENTA15, EL.THHM_PENTA15),
            (MT.TRIA6, EL.THHM_FACE6),
            (MT.QUAD8, EL.THHM_FACE8),
        ),
    ),
)

phen.add(
    "3D_THHMD",
    Modelisation(
        dim=(3, 3),
        code="3H9",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.HEXA20, EL.THHM_HEXA20D),
            (MT.TETRA10, EL.THHM_TETRA10D),
            (MT.PENTA15, EL.THHM_PENTA15D),
            (MT.TRIA6, EL.THHM_FACE6),
            (MT.QUAD8, EL.THHM_FACE8),
        ),
    ),
)

phen.add(
    "3D_THHMS",
    Modelisation(
        dim=(3, 3),
        code="3R5",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.HEXA20, EL.THHM_HEXA20S),
            (MT.TETRA10, EL.THHM_TETRA10S),
            (MT.PENTA15, EL.THHM_PENTA15S),
            (MT.TRIA6, EL.THHM_FACE6),
            (MT.QUAD8, EL.THHM_FACE8),
        ),
    ),
)

phen.add(
    "3D_THHS",
    Modelisation(
        dim=(3, 3),
        code="3R4",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.THH_TETRA10S),
            (MT.PENTA15, EL.THH_PENTA15S),
            (MT.HEXA20, EL.THH_HEXA20S),
            (MT.QUAD8, EL.THH_FACE8),
            (MT.TRIA6, EL.THH_FACE6),
        ),
    ),
)

phen.add(
    "3D_THM",
    Modelisation(
        dim=(3, 3),
        code="3H5",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
        ),
        elements=(
            (MT.TETRA10, EL.THM_TETRA10),
            (MT.PENTA15, EL.THM_PENTA15),
            (MT.HEXA20, EL.THM_HEXA20),
            (MT.QUAD8, EL.THM_FACE8),
            (MT.TRIA6, EL.THM_FACE6),
        ),
    ),
)

phen.add(
    "3D_THMD",
    Modelisation(
        dim=(3, 3),
        code="3H0",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.TETRA10, EL.THM_TETRA10D),
            (MT.PENTA15, EL.THM_PENTA15D),
            (MT.HEXA20, EL.THM_HEXA20D),
            (MT.QUAD8, EL.THM_FACE8),
            (MT.TRIA6, EL.THM_FACE6),
        ),
    ),
)

phen.add(
    "3D_THMS",
    Modelisation(
        dim=(3, 3),
        code="3R6",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.THM_TETRA10S),
            (MT.PENTA15, EL.THM_PENTA15S),
            (MT.HEXA20, EL.THM_HEXA20S),
            (MT.QUAD8, EL.THM_FACE8),
            (MT.TRIA6, EL.THM_FACE6),
        ),
    ),
)

phen.add(
    "3D_THMS_DIL",
    Modelisation(
        dim=(3, 3),
        code="3D6",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.INTTHM, "RED"),
            (AT.DIL, "OUI"),
            (AT.FORMULATION, "DIL"),
        ),
        elements=(
            (MT.TETRA10, EL.THM_TETRA10S_DIL),
            (MT.PENTA15, EL.THM_PENTA15S_DIL),
            (MT.HEXA20, EL.THM_HEXA20S_DIL),
            (MT.QUAD8, EL.THM_FACE8),
            (MT.TRIA6, EL.THM_FACE6),
        ),
    ),
)

phen.add(
    "3D_THVD",
    Modelisation(
        dim=(3, 3),
        code="3I3",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "0"),
        ),
        elements=(
            (MT.TETRA10, EL.THV_TETRA10D),
            (MT.PENTA15, EL.THV_PENTA15D),
            (MT.HEXA20, EL.THV_HEXA20D),
            (MT.QUAD8, EL.THV_FACE8),
            (MT.TRIA6, EL.THV_FACE6),
        ),
    ),
)

phen.add(
    "3D_THVS",
    Modelisation(
        dim=(3, 3),
        code="3R3",
        attrs=(
            (AT.TYPMOD, "3D"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "0"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.TETRA10, EL.THV_TETRA10S),
            (MT.PENTA15, EL.THV_PENTA15S),
            (MT.HEXA20, EL.THV_HEXA20S),
            (MT.QUAD8, EL.THV_FACE8),
            (MT.TRIA6, EL.THV_FACE6),
        ),
    ),
)

phen.add(
    "AXIS",
    Modelisation(
        dim=(2, 2),
        code="AX_",
        attrs=((AT.AXIS, "OUI"), (AT.NBSIGM, "4"), (AT.TYPMOD, "AXIS")),
        elements=(
            (MT.TRIA3, EL.MEAXTR3),
            (MT.QUAD4, EL.MEAXQU4),
            (MT.SEG2, EL.MEAXSE2),
            (MT.TRIA6, EL.MEAXTR6),
            (MT.QUAD8, EL.MEAXQU8),
            (MT.QUAD9, EL.MEAXQU9),
            (MT.SEG3, EL.MEAXSE3),
        ),
    ),
)

phen.add(
    "AXIS_FLUIDE#1",
    Modelisation(
        dim=(2, 2),
        code="AXF",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.FLUIDE, "OUI"),
            (AT.AXIS, "OUI"),
            (AT.FORMULATION, "U_P_PHI"),
        ),
        elements=(
            (MT.TRIA3, EL.MEAXFLT3),
            (MT.QUAD4, EL.MEAXFLQ4),
            (MT.TRIA6, EL.MEAXFLT6),
            (MT.QUAD8, EL.MEAXFLQ8),
            (MT.QUAD9, EL.MEAXFLQ9),
            (MT.SEG2, EL.MEAXFLS2),
            (MT.SEG3, EL.MEAXFLS3),
        ),
    ),
)

phen.add(
    "AXIS_FLUIDE#2",
    Modelisation(
        dim=(2, 2),
        code="AX2",
        attrs=((AT.TYPMOD, "AXIS"), (AT.FLUIDE, "OUI"), (AT.AXIS, "OUI"), (AT.FORMULATION, "U_P")),
        elements=(
            (MT.TRIA3, EL.MEAXFLT3P),
            (MT.QUAD4, EL.MEAXFLQ4P),
            (MT.TRIA6, EL.MEAXFLT6P),
            (MT.QUAD8, EL.MEAXFLQ8P),
            (MT.QUAD9, EL.MEAXFLQ9P),
            (MT.SEG2, EL.MEAXFLS2P),
            (MT.SEG3, EL.MEAXFLS3P),
        ),
    ),
)

phen.add(
    "AXIS_FLUIDE#3",
    Modelisation(
        dim=(2, 2),
        code="AX3",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.FLUIDE, "OUI"),
            (AT.AXIS, "OUI"),
            (AT.FORMULATION, "U_PSI"),
        ),
        elements=(
            (MT.TRIA3, EL.MEAXFLT3PSI),
            (MT.QUAD4, EL.MEAXFLQ4PSI),
            (MT.TRIA6, EL.MEAXFLT6PSI),
            (MT.QUAD8, EL.MEAXFLQ8PSI),
            (MT.QUAD9, EL.MEAXFLQ9PSI),
            (MT.SEG2, EL.MEAXFLS2PSI),
            (MT.SEG3, EL.MEAXFLS3PSI),
        ),
    ),
)

phen.add(
    "AXIS_FLUI_STRU#1",
    Modelisation(
        dim=(1, 2),
        code="FSA",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.FLUIDE, "OUI"),
            (AT.AXIS, "OUI"),
            (AT.FSI, "OUI"),
            (AT.FORMULATION, "U_P_PHI"),
        ),
        elements=((MT.SEG2, EL.MEAXFSS2), (MT.SEG3, EL.MEAXFSS3)),
    ),
)

phen.add(
    "AXIS_FLUI_STRU#2",
    Modelisation(
        dim=(1, 2),
        code="FS2",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.FLUIDE, "OUI"),
            (AT.AXIS, "OUI"),
            (AT.FSI, "OUI"),
            (AT.FORMULATION, "U_P"),
        ),
        elements=((MT.SEG2, EL.MEAXFSS2P), (MT.SEG3, EL.MEAXFSS3P)),
    ),
)

phen.add(
    "AXIS_FLUI_STRU#3",
    Modelisation(
        dim=(1, 2),
        code="FS3",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.FLUIDE, "OUI"),
            (AT.AXIS, "OUI"),
            (AT.FSI, "OUI"),
            (AT.FORMULATION, "U_PSI"),
        ),
        elements=((MT.SEG2, EL.MEAXFSS2PSI), (MT.SEG3, EL.MEAXFSS3PSI)),
    ),
)


phen.add(
    "AXIS_FOURIER",
    Modelisation(
        dim=(2, 2),
        code="AFO",
        attrs=((AT.AXIS, "OUI"), (AT.NBSIGM, "6"), (AT.FOURIER, "OUI"), (AT.TYPMOD, "AXIS")),
        elements=(
            (MT.TRIA3, EL.MEFOTR3),
            (MT.QUAD4, EL.MEFOQU4),
            (MT.SEG2, EL.MEFOSE2),
            (MT.TRIA6, EL.MEFOTR6),
            (MT.QUAD8, EL.MEFOQU8),
            (MT.QUAD9, EL.MEFOQU9),
            (MT.SEG3, EL.MEFOSE3),
        ),
    ),
)

phen.add(
    "AXIS_GRAD_VARI",
    Modelisation(
        dim=(2, 2),
        code="AXV",
        attrs=((AT.AXIS, "OUI"), (AT.NBSIGM, "4"), (AT.TYPMOD, "AXIS"), (AT.TYPMOD2, "GRADVARI")),
        elements=((MT.SEG3, EL.MEAXSE3), (MT.TRIA6, EL.MVAXTR6), (MT.QUAD8, EL.MVAXQS8)),
    ),
)

phen.add(
    "AXIS_GRAD_INCO",
    Modelisation(
        dim=(2, 2),
        code="AGI",
        attrs=(
            (AT.AXIS, "OUI"),
            (AT.NBSIGM, "4"),
            (AT.TYPMOD, "AXIS"),
            (AT.INCO, "C5GV"),
            (AT.TYPMOD2, "GRADVARI"),
        ),
        elements=((MT.SEG3, EL.MEAXSE3), (MT.TRIA6, EL.GVI_AX_TR6), (MT.QUAD8, EL.GVI_AX_QU8)),
    ),
)

phen.add(
    "AXIS_HH2D",
    Modelisation(
        dim=(2, 2),
        code="AZ4",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.HH2_AXIS_QU8D),
            (MT.TRIA6, EL.HH2_AXIS_TR6D),
            (MT.SEG3, EL.HH2_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HH2MD",
    Modelisation(
        dim=(2, 2),
        code="AA3",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.HH2M_AXIS_QU8D),
            (MT.TRIA6, EL.HH2M_AXIS_TR6D),
            (MT.SEG3, EL.HH2M_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HH2MS",
    Modelisation(
        dim=(2, 2),
        code="AR2",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.HH2M_AXIS_QU8S),
            (MT.TRIA6, EL.HH2M_AXIS_TR6S),
            (MT.SEG3, EL.HH2M_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HH2S",
    Modelisation(
        dim=(2, 2),
        code="AZ3",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.HH2_AXIS_QU8S),
            (MT.TRIA6, EL.HH2_AXIS_TR6S),
            (MT.SEG3, EL.HH2_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HHD",
    Modelisation(
        dim=(2, 2),
        code="AZ2",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.HH_AXIS_QU8D),
            (MT.TRIA6, EL.HH_AXIS_TR6D),
            (MT.SEG3, EL.HH_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HHM",
    Modelisation(
        dim=(2, 2),
        code="AH1",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.AXIS, "OUI"),
        ),
        elements=(
            (MT.QUAD8, EL.HHM_AXIS_QU8),
            (MT.TRIA6, EL.HHM_AXIS_TR6),
            (MT.SEG3, EL.HHM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HHMD",
    Modelisation(
        dim=(2, 2),
        code="AH6",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.HHM_AXIS_QU8D),
            (MT.TRIA6, EL.HHM_AXIS_TR6D),
            (MT.SEG3, EL.HHM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HHMS",
    Modelisation(
        dim=(2, 2),
        code="AR1",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.HHM_AXIS_QU8S),
            (MT.TRIA6, EL.HHM_AXIS_TR6S),
            (MT.SEG3, EL.HHM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HHS",
    Modelisation(
        dim=(2, 2),
        code="AZ1",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.HH_AXIS_QU8S),
            (MT.TRIA6, EL.HH_AXIS_TR6S),
            (MT.SEG3, EL.HH_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HM",
    Modelisation(
        dim=(2, 2),
        code="AH2",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.AXIS, "OUI"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_AXIS_QU8),
            (MT.TRIA6, EL.HM_AXIS_TR6),
            (MT.SEG3, EL.HM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HMD",
    Modelisation(
        dim=(2, 2),
        code="AH7",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_AXIS_QU8D),
            (MT.TRIA6, EL.HM_AXIS_TR6D),
            (MT.SEG3, EL.HM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_HMS",
    Modelisation(
        dim=(2, 2),
        code="AR3",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_AXIS_QU8S),
            (MT.TRIA6, EL.HM_AXIS_TR6S),
            (MT.SEG3, EL.HM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_INCO_UP",
    Modelisation(
        dim=(2, 2),
        code="AUP",
        attrs=((AT.AXIS, "OUI"), (AT.NBSIGM, "4"), (AT.TYPMOD, "AXIS"), (AT.INCO, "C2")),
        elements=(
            (MT.TRIA6, EL.MUAXTR6),
            (MT.TRIA3, EL.MUAXTR3),
            (MT.QUAD8, EL.MUAXQU8),
            (MT.SEG3, EL.MEAXSE3),
            (MT.SEG2, EL.MEAXSE2),
        ),
    ),
)

phen.add(
    "AXIS_INCO_UPG",
    Modelisation(
        dim=(2, 2),
        code="AXC",
        attrs=((AT.AXIS, "OUI"), (AT.NBSIGM, "4"), (AT.TYPMOD, "AXIS"), (AT.INCO, "C3")),
        elements=((MT.TRIA6, EL.MIAXTR6), (MT.QUAD8, EL.MIAXQU8), (MT.SEG3, EL.MEAXSE3)),
    ),
)

phen.add(
    "AXIS_INCO_UPO",
    Modelisation(
        dim=(2, 2),
        code="AOS",
        attrs=((AT.AXIS, "OUI"), (AT.NBSIGM, "4"), (AT.TYPMOD, "AXIS"), (AT.INCO, "C2O")),
        elements=(
            (MT.TRIA3, EL.MIAXOSTR3),
            (MT.QUAD4, EL.MIAXOSQU4),
            (MT.SEG3, EL.MEAXSE3),
            (MT.SEG2, EL.MEAXSE2),
        ),
    ),
)

phen.add(
    "AXIS_INTERFACE",
    Modelisation(
        dim=(2, 2),
        code="AEI",
        attrs=(
            (AT.AXIS, "OUI"),
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "INTERFAC"),
            (AT.INTERFACE, "OUI"),
        ),
        elements=((MT.QUAD8, EL.EIAXQU8),),
    ),
)

phen.add(
    "AXIS_INTERFACE_S",
    Modelisation(
        dim=(2, 2),
        code="AIS",
        attrs=(
            (AT.AXIS, "OUI"),
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "INTERFAC"),
            (AT.INTERFACE, "OUI"),
        ),
        elements=((MT.QUAD8, EL.EIAXQS8),),
    ),
)

phen.add(
    "AXIS_JHMS",
    Modelisation(
        dim=(2, 2),
        code="JH2",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "JHMS"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.HM_J_AXQ8S), (MT.SEG3, EL.HM_J_AXSE3)),
    ),
)

phen.add(
    "AXIS_JOINT",
    Modelisation(
        dim=(2, 2),
        code="AFI",
        attrs=(
            (AT.AXIS, "OUI"),
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "ELEMJOIN"),
            (AT.INTERFACE, "OUI"),
        ),
        elements=((MT.QUAD4, EL.MFAXQU4),),
    ),
)

phen.add(
    "AXIS_SI",
    Modelisation(
        dim=(2, 2),
        code="AXS",
        attrs=((AT.AXIS, "OUI"), (AT.NBSIGM, "4"), (AT.TYPMOD, "AXIS")),
        elements=((MT.QUAD8, EL.MEAXQS8), (MT.SEG3, EL.MEAXSE3)),
    ),
)

phen.add(
    "AXIS_THH2D",
    Modelisation(
        dim=(2, 2),
        code="AA1",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.THH2_AXIS_QU8D),
            (MT.TRIA6, EL.THH2_AXIS_TR6D),
            (MT.SEG3, EL.THH2_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THH2MD",
    Modelisation(
        dim=(2, 2),
        code="AA2",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.THH2M_AXIS_QU8D),
            (MT.TRIA6, EL.THH2M_AXIS_TR6D),
            (MT.SEG3, EL.THH2M_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THH2MS",
    Modelisation(
        dim=(2, 2),
        code="AR7",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.THH2M_AXIS_QU8S),
            (MT.TRIA6, EL.THH2M_AXIS_TR6S),
            (MT.SEG3, EL.THH2M_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THH2S",
    Modelisation(
        dim=(2, 2),
        code="AR5",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.THH2_AXIS_QU8S),
            (MT.TRIA6, EL.THH2_AXIS_TR6S),
            (MT.SEG3, EL.THH2_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THHD",
    Modelisation(
        dim=(2, 2),
        code="AH8",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.THH_AXIS_QU8D),
            (MT.TRIA6, EL.THH_AXIS_TR6D),
            (MT.SEG3, EL.THH_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THHMD",
    Modelisation(
        dim=(2, 2),
        code="AH9",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.THHM_AXIS_QU8D),
            (MT.TRIA6, EL.THHM_AXIS_TR6D),
            (MT.SEG3, EL.THHM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THHMS",
    Modelisation(
        dim=(2, 2),
        code="AR6",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.THHM_AXIS_QU8S),
            (MT.TRIA6, EL.THHM_AXIS_TR6S),
            (MT.SEG3, EL.THHM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THHS",
    Modelisation(
        dim=(2, 2),
        code="AR4",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.THH_AXIS_QU8S),
            (MT.TRIA6, EL.THH_AXIS_TR6S),
            (MT.SEG3, EL.THH_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THM",
    Modelisation(
        dim=(2, 2),
        code="AH5",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.AXIS, "OUI"),
        ),
        elements=(
            (MT.QUAD8, EL.THM_AXIS_QU8),
            (MT.TRIA6, EL.THM_AXIS_TR6),
            (MT.SEG3, EL.THM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THMD",
    Modelisation(
        dim=(2, 2),
        code="AH0",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.THM_AXIS_QU8D),
            (MT.TRIA6, EL.THM_AXIS_TR6D),
            (MT.SEG3, EL.THM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THMS",
    Modelisation(
        dim=(2, 2),
        code="AR8",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.THM_AXIS_QU8S),
            (MT.TRIA6, EL.THM_AXIS_TR6S),
            (MT.SEG3, EL.THM_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THVD",
    Modelisation(
        dim=(2, 2),
        code="AG3",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "0"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.THV_AXIS_QU8D),
            (MT.TRIA6, EL.THV_AXIS_TR6D),
            (MT.SEG3, EL.THV_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_THVS",
    Modelisation(
        dim=(2, 2),
        code="AR9",
        attrs=(
            (AT.TYPMOD, "AXIS"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "0"),
            (AT.AXIS, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.THV_AXIS_QU8S),
            (MT.TRIA6, EL.THV_AXIS_TR6S),
            (MT.SEG3, EL.THV_AXIS_SE3),
        ),
    ),
)

phen.add(
    "AXIS_XH",
    Modelisation(
        dim=(2, 2),
        code="AX1",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.AXIS, "OUI"),
            (AT.TYPMOD, "AXIS"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=(
            (MT.TRIA6, EL.MEAXTR6_XH),
            (MT.TRIA3, EL.MEAXTR3_XH),
            (MT.QUAD8, EL.MEAXQU8_XH),
            (MT.QUAD4, EL.MEAXQU4_XH),
            (MT.SEG3, EL.MEAXSE3_XH),
            (MT.SEG2, EL.MEAXSE2_XH),
        ),
    ),
)

phen.add(
    "AXIS_XHT",
    Modelisation(
        dim=(2, 2),
        code="AX3",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.AXIS, "OUI"),
            (AT.TYPMOD, "AXIS"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XHT"),
        ),
        elements=(
            (MT.TRIA6, EL.MEAXTR6_XHT),
            (MT.TRIA3, EL.MEAXTR3_XHT),
            (MT.QUAD8, EL.MEAXQU8_XHT),
            (MT.QUAD4, EL.MEAXQU4_XHT),
            (MT.SEG3, EL.MEAXSE3_XHT),
            (MT.SEG2, EL.MEAXSE2_XHT),
        ),
    ),
)

phen.add(
    "AXIS_XT",
    Modelisation(
        dim=(2, 2),
        code="AX2",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.AXIS, "OUI"),
            (AT.TYPMOD, "AXIS"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XT"),
        ),
        elements=(
            (MT.TRIA6, EL.MEAXTR6_XT),
            (MT.TRIA3, EL.MEAXTR3_XT),
            (MT.QUAD8, EL.MEAXQU8_XT),
            (MT.QUAD4, EL.MEAXQU4_XT),
            (MT.SEG3, EL.MEAXSE3_XT),
            (MT.SEG2, EL.MEAXSE2_XT),
        ),
    ),
)

phen.add(
    "BARRE",
    Modelisation(
        dim=(1, 3),
        code="BAR",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.TYPMOD, "1D"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
            (AT.STRX, "NON"),
        ),
        elements=((MT.SEG2, EL.MECA_BARRE),),
    ),
)

phen.add(
    "BARRE_2D_NL",
    Modelisation(
        dim=(1, 2),
        code="B2D",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.TYPMOD, "1D"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
            (AT.STRX, "NON"),
        ),
        elements=((MT.SEG2, EL.MEBA2DSE2), (MT.SEG3, EL.MEBA2DSE3)),
    ),
)

phen.add(
    "BARRE_3D_NL",
    Modelisation(
        dim=(1, 3),
        code="B3D",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.TYPMOD, "1D"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
            (AT.STRX, "NON"),
        ),
        elements=((MT.SEG2, EL.MEBA3DSE2), (MT.SEG3, EL.MEBA3DSE3)),
    ),
)

phen.add(
    "CABLE",
    Modelisation(
        dim=(1, 3),
        code="CAB",
        attrs=(
            (AT.TYPMOD, "1D"),
            (AT.POUTRE, "OUI"),
            (AT.CABLE, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
            (AT.STRX, "NON"),
        ),
        elements=((MT.SEG2, EL.MECABL2),),
    ),
)

phen.add(
    "CABLE_GAINE",
    Modelisation(
        dim=(1, 3),
        code="1GC",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.TYPMOD, "1D"),
            (AT.TYPMOD2, "CABLE_GA"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
            (AT.STRX, "NON"),
        ),
        elements=((MT.SEG3, EL.MECGSEG3),),
    ),
)

phen.add(
    "CABLE_POULIE",
    Modelisation(
        dim=(1, 3),
        code="CAP",
        attrs=(
            (AT.TYPMOD, "1D"),
            (AT.POUTRE, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
            (AT.STRX, "NON"),
        ),
        elements=((MT.SEG3, EL.MEPOULI),),
    ),
)

# ------------------------------------------------------------------------------------
# Modelisations sous-terraines pour les raccords
# ------------------------------------------------------------------------------------


phen.add(
    "RACC_COQ_3D",
    Modelisation(
        dim=(2, 3),
        code="RC3",
        attrs=((AT.CONTACT, "OUI"),),
        elements=(
            (MT.SE2TR3, EL.RACS2T3),
            (MT.SE2TR6, EL.RACS2T6),
            (MT.SE2QU4, EL.RACS2Q4),
            (MT.SE2QU8, EL.RACS2Q8),
            (MT.SE3TR3, EL.RACS3T3),
            (MT.SE3TR6, EL.RACS3T6),
            (MT.SE3QU4, EL.RACS3Q4),
            (MT.SE3QU8, EL.RACS3Q8),
        ),
    ),
)


# ------------------------------------------------------------------------------------
# Modelisations sous-terraines pour le contact/frottement en methode "continue" :
# ------------------------------------------------------------------------------------

# -- Define SLAVE elements for CONTINUE method (in DEFI_CONTACT) - Contact

phen.add(
    "CONT_SL_2D",
    Modelisation(dim=(1, 2), code="CT2", elements=((MT.SEG2, EL.COS2E2D), (MT.SEG3, EL.COS3E2D))),
)

phen.add(
    "CONT_SL_3D",
    Modelisation(
        dim=(2, 3),
        code="CT3",
        elements=(
            (MT.TRIA3, EL.COT3E3D),
            (MT.QUAD4, EL.COQ4E3D),
            (MT.QUAD8, EL.COQ8E3D),
            (MT.TRIA6, EL.COT6E3D),
            (MT.QUAD9, EL.COQ9E3D),
            (MT.SEG2, EL.COP2E3D),
            (MT.SEG3, EL.COP3E3D),
        ),
    ),
)

# -- Define SLAVE elements for LAC method (in DEFI_CONTACT) - Contact

phen.add(
    "CONT_LAC_SL_2D",
    Modelisation(
        dim=(1, 2),
        code="LC2",
        attrs=((AT.CONTACT, "OUI"),),
        elements=((MT.SEG2, EL.LACS22D), (MT.SEG3, EL.LACS32D)),
    ),
)

phen.add(
    "CONT_LAC_SL_2DB",
    Modelisation(
        dim=(1, 2), code="L2B", attrs=((AT.CONTACT, "OUI"),), elements=((MT.SEG2, EL.LACS22DB),)
    ),
)

phen.add(
    "CONT_LAC_SL_2DT",
    Modelisation(
        dim=(1, 2), code="L2T", attrs=((AT.CONTACT, "OUI"),), elements=((MT.SEG2, EL.LACS22DT),)
    ),
)

phen.add(
    "CONT_LAC_SL_3D",
    Modelisation(
        dim=(2, 3),
        code="LC3",
        attrs=((AT.CONTACT, "OUI"),),
        elements=(
            (MT.TRIA3, EL.LACT33D),
            (MT.TRIA6, EL.LACT63D),
            (MT.QUAD9, EL.LACQ93D),
            (MT.QUAD8, EL.LACQ83D),
            (MT.QUAD4, EL.LACQ43D),
        ),
    ),
)

phen.add(
    "CONT_LAC_SL_3DB",
    Modelisation(
        dim=(2, 3),
        code="L3B",
        attrs=((AT.CONTACT, "OUI"),),
        elements=((MT.QUAD8, EL.LACQ83DB), (MT.QUAD4, EL.LACQ43DB)),
    ),
)

# -- Define SLAVE elements for CONTINUE method (in DEFI_CONTACT) - Friction

phen.add(
    "FRIC_SL_2D",
    Modelisation(dim=(1, 2), code="CF2", elements=((MT.SEG2, EL.CFS2E2D), (MT.SEG3, EL.CFS3E2D))),
)

phen.add(
    "FRIC_SL_3D",
    Modelisation(
        dim=(2, 3),
        code="CF3",
        elements=(
            (MT.TRIA3, EL.CFT3E3D),
            (MT.QUAD4, EL.CFQ4E3D),
            (MT.QUAD8, EL.CFQ8E3D),
            (MT.TRIA6, EL.CFT6E3D),
            (MT.QUAD9, EL.CFQ9E3D),
            (MT.SEG2, EL.CFP2E3D),
        ),
    ),
)

# -- Define SLAVE elements for LAGRANGIAN method (in DEFI_CONTACT) - Contact

phen.add(
    "CONT_LAG_SL_2D",
    Modelisation(dim=(1, 2), code="CM2", elements=((MT.SEG2, EL.CMS22D), (MT.SEG3, EL.CMS32D))),
)

phen.add(
    "CONT_LAG_SL_3D",
    Modelisation(
        dim=(2, 3),
        code="CM3",
        elements=(
            (MT.TRIA3, EL.CMT33D),
            (MT.TRIA6, EL.CMT63D),
            (MT.QUAD4, EL.CMQ43D),
            (MT.QUAD8, EL.CMQ83D),
            (MT.QUAD9, EL.CMQ93D),
        ),
    ),
)


phen.add(
    "FRIC_LAG_SL_2D",
    Modelisation(dim=(1, 2), code="FM2", elements=((MT.SEG2, EL.FMS22D), (MT.SEG3, EL.FMS32D))),
)

phen.add(
    "FRIC_LAG_SL_3D",
    Modelisation(
        dim=(2, 3),
        code="FM3",
        elements=(
            (MT.TRIA3, EL.FMT33D),
            (MT.TRIA6, EL.FMT63D),
            (MT.QUAD4, EL.FMQ43D),
            (MT.QUAD8, EL.FMQ83D),
            (MT.QUAD9, EL.FMQ93D),
        ),
    ),
)


# -- Define SLAVE elements for NITSCHE method (in DEFI_CONTACT) - Contact

phen.add(
    "CONT_NIT_SL_2D",
    Modelisation(dim=(1, 2), code="CN2", elements=((MT.SEG2, EL.CNS22D), (MT.SEG3, EL.CNS32D))),
)

phen.add(
    "CONT_NIT_SL_3D",
    Modelisation(
        dim=(2, 3),
        code="CN3",
        elements=(
            (MT.TRIA3, EL.CNT33D),
            (MT.TRIA6, EL.CNT63D),
            (MT.QUAD4, EL.CNQ43D),
            (MT.QUAD8, EL.CNQ83D),
            (MT.QUAD9, EL.CNQ93D),
        ),
    ),
)


# -- Define SLAVE elements for PENALISATION method (in DEFI_CONTACT) - Contact

phen.add(
    "CONT_PENA_SL_2D",
    Modelisation(dim=(1, 2), code="CP2", elements=((MT.SEG2, EL.CPS22D), (MT.SEG3, EL.CPS32D))),
)

phen.add(
    "CONT_PENA_SL_3D",
    Modelisation(
        dim=(2, 3),
        code="CP3",
        elements=(
            (MT.TRIA3, EL.CPT33D),
            (MT.TRIA6, EL.CPT63D),
            (MT.QUAD4, EL.CPQ43D),
            (MT.QUAD8, EL.CPQ83D),
            (MT.QUAD9, EL.CPQ93D),
        ),
    ),
)

# -- Define CONTACT elements for CONTINUE method (in STAT_NON_LINE) - Contact

phen.add(
    "CONT_EL_3D1",
    Modelisation(
        dim=(2, 3),
        code="CC1",
        attrs=((AT.CONTACT, "OUI"),),
        elements=(
            (MT.QU4QU8, EL.COQ4Q8),
            (MT.QU4QU9, EL.COQ4Q9),
            (MT.QU4TR3, EL.COQ4T3),
            (MT.QU4TR6, EL.COQ4T6),
            (MT.QU8QU4, EL.COQ8Q4),
            (MT.QU8QU9, EL.COQ8Q9),
            (MT.QU8TR3, EL.COQ8T3),
            (MT.QU8TR6, EL.COQ8T6),
            (MT.QU9QU4, EL.COQ9Q4),
            (MT.QU9QU8, EL.COQ9Q8),
            (MT.QU9TR3, EL.COQ9T3),
            (MT.QU9TR6, EL.COQ9T6),
            (MT.QUAD44, EL.COQ4Q4),
            (MT.QUAD88, EL.COQ8Q8),
            (MT.QUAD99, EL.COQ9Q9),
            (MT.SE2QU4, EL.COS2Q4),
            (MT.SE2QU8, EL.COS2Q8),
            (MT.SE2QU9, EL.COS2Q9),
            (MT.SE2TR3, EL.COS2T3),
            (MT.SE2TR6, EL.COS2T6),
            (MT.SE3QU4, EL.COS3Q4),
            (MT.SE3QU8, EL.COS3Q8),
            (MT.SE3QU9, EL.COS3Q9),
            (MT.SE3TR3, EL.COS3T3),
            (MT.SE3TR6, EL.COS3T6),
            (MT.TR3QU4, EL.COT3Q4),
            (MT.TR3QU8, EL.COT3Q8),
            (MT.TR3QU9, EL.COT3Q9),
            (MT.TR3TR6, EL.COT3T6),
            (MT.TR6QU4, EL.COT6Q4),
            (MT.TR6QU8, EL.COT6Q8),
            (MT.TR6QU9, EL.COT6Q9),
            (MT.TR6TR3, EL.COT6T3),
            (MT.TRIA33, EL.COT3T3),
            (MT.TRIA66, EL.COT6T6),
        ),
    ),
)

phen.add(
    "CONT_EL_3D2",
    Modelisation(
        dim=(1, 2), code="CC4", attrs=((AT.CONTACT, "OUI"),), elements=((MT.SEG22, EL.COP2P2),)
    ),
)

phen.add(
    "CONT_EL_2D1",
    Modelisation(
        dim=(1, 2),
        code="CC2",
        attrs=((AT.CONTACT, "OUI"),),
        elements=(
            (MT.SEG22, EL.COS2S2),
            (MT.SEG23, EL.COS2S3),
            (MT.SEG32, EL.COS3S2),
            (MT.SEG33, EL.COS3S3),
        ),
    ),
)

phen.add(
    "CONT_EL_2D2",
    Modelisation(
        dim=(1, 2),
        code="CC3",
        attrs=((AT.CONTACT, "OUI"), (AT.AXIS, "OUI")),
        elements=(
            (MT.SEG22, EL.COS2S2A),
            (MT.SEG23, EL.COS2S3A),
            (MT.SEG32, EL.COS3S2A),
            (MT.SEG33, EL.COS3S3A),
        ),
    ),
)

# -- Define CONTACT elements for CONTINUE method (in STAT_NON_LINE) - Friction

phen.add(
    "FRIC_EL_3D1",
    Modelisation(
        dim=(2, 3),
        code="CC5",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI")),
        elements=(
            (MT.QU4QU8, EL.CFQ4Q8),
            (MT.QU4QU9, EL.CFQ4Q9),
            (MT.QU4TR3, EL.CFQ4T3),
            (MT.QU4TR6, EL.CFQ4T6),
            (MT.QU8QU4, EL.CFQ8Q4),
            (MT.QU8QU9, EL.CFQ8Q9),
            (MT.QU8TR3, EL.CFQ8T3),
            (MT.QU8TR6, EL.CFQ8T6),
            (MT.QU9QU4, EL.CFQ9Q4),
            (MT.QU9QU8, EL.CFQ9Q8),
            (MT.QU9TR3, EL.CFQ9T3),
            (MT.QU9TR6, EL.CFQ9T6),
            (MT.QUAD44, EL.CFQ4Q4),
            (MT.QUAD88, EL.CFQ8Q8),
            (MT.QUAD99, EL.CFQ9Q9),
            (MT.SE2QU4, EL.CFS2Q4),
            (MT.SE2QU8, EL.CFS2Q8),
            (MT.SE2QU9, EL.CFS2Q9),
            (MT.SE2TR3, EL.CFS2T3),
            (MT.SE2TR6, EL.CFS2T6),
            (MT.SE3QU4, EL.CFS3Q4),
            (MT.SE3QU8, EL.CFS3Q8),
            (MT.SE3QU9, EL.CFS3Q9),
            (MT.SE3TR3, EL.CFS3T3),
            (MT.SE3TR6, EL.CFS3T6),
            (MT.TR3QU4, EL.CFT3Q4),
            (MT.TR3QU8, EL.CFT3Q8),
            (MT.TR3QU9, EL.CFT3Q9),
            (MT.TR3TR6, EL.CFT3T6),
            (MT.TR6QU4, EL.CFT6Q4),
            (MT.TR6QU8, EL.CFT6Q8),
            (MT.TR6QU9, EL.CFT6Q9),
            (MT.TR6TR3, EL.CFT6T3),
            (MT.TRIA33, EL.CFT3T3),
            (MT.TRIA66, EL.CFT6T6),
        ),
    ),
)

phen.add(
    "FRIC_EL_3D2",
    Modelisation(
        dim=(1, 2),
        code="CC8",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI")),
        elements=((MT.SEG22, EL.CFP2P2),),
    ),
)

phen.add(
    "FRIC_EL_2D1",
    Modelisation(
        dim=(1, 2),
        code="CC6",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI")),
        elements=(
            (MT.SEG22, EL.CFS2S2),
            (MT.SEG23, EL.CFS2S3),
            (MT.SEG32, EL.CFS3S2),
            (MT.SEG33, EL.CFS3S3),
        ),
    ),
)

phen.add(
    "FRIC_EL_2D2",
    Modelisation(
        dim=(1, 2),
        code="CC7",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI"), (AT.AXIS, "OUI")),
        elements=(
            (MT.SEG22, EL.CFS2S2A),
            (MT.SEG23, EL.CFS2S3A),
            (MT.SEG32, EL.CFS3S2A),
            (MT.SEG33, EL.CFS3S3A),
        ),
    ),
)

# -- Define CONTACT elements for LAC method (in STAT_NON_LINE) - Contact

phen.add(
    "CONT_LAC_EL_2DC",
    Modelisation(
        dim=(1, 2),
        code="L2C",
        attrs=((AT.CONTACT, "OUI"),),
        elements=(
            (MT.SEG22, EL.LCS2S2C),
            (MT.SEG33, EL.LCS3S3C),
            (MT.SEG23, EL.LCS2S3C),
            (MT.SEG32, EL.LCS3S2C),
        ),
    ),
)


phen.add(
    "CONT_LAC_EL_2DD",
    Modelisation(
        dim=(1, 2),
        code="L2D",
        attrs=((AT.CONTACT, "OUI"),),
        elements=((MT.SEG22, EL.LCS2S2D), (MT.SEG23, EL.LCS2S3D)),
    ),
)


phen.add(
    "CONT_LAC_EL_2DE",
    Modelisation(
        dim=(1, 2),
        code="L2E",
        attrs=((AT.CONTACT, "OUI"),),
        elements=((MT.SEG22, EL.LCS2S2E), (MT.SEG23, EL.LCS2S3E)),
    ),
)

phen.add(
    "CONT_LAC_EL_2DF",
    Modelisation(
        dim=(1, 2),
        code="L2F",
        attrs=((AT.CONTACT, "OUI"), (AT.AXIS, "OUI")),
        elements=(
            (MT.SEG22, EL.LCS2S2CA),
            (MT.SEG33, EL.LCS3S3CA),
            (MT.SEG23, EL.LCS2S3CA),
            (MT.SEG32, EL.LCS3S2CA),
        ),
    ),
)


phen.add(
    "CONT_LAC_EL_2DG",
    Modelisation(
        dim=(1, 2),
        code="L2G",
        attrs=((AT.CONTACT, "OUI"), (AT.AXIS, "OUI")),
        elements=((MT.SEG22, EL.LCS2S2DA), (MT.SEG23, EL.LCS2S3DA)),
    ),
)


phen.add(
    "CONT_LAC_EL_2DH",
    Modelisation(
        dim=(1, 2),
        code="L2H",
        attrs=((AT.CONTACT, "OUI"), (AT.AXIS, "OUI")),
        elements=((MT.SEG22, EL.LCS2S2EA), (MT.SEG23, EL.LCS2S3EA)),
    ),
)

phen.add(
    "CONT_LAC_EL_3DD",
    Modelisation(
        dim=(2, 3),
        code="L3D",
        attrs=((AT.CONTACT, "OUI"),),
        elements=(
            (MT.QUAD44, EL.LACQ4Q4D),
            (MT.QU4TR3, EL.LACQ4T3D),
            (MT.QU4TR6, EL.LACQ4T6D),
            (MT.QU4QU8, EL.LACQ4Q8D),
            (MT.QU4QU9, EL.LACQ4Q9D),
            (MT.QUAD88, EL.LACQ8Q8D),
            (MT.QU8TR6, EL.LACQ8T6D),
            (MT.QU8TR3, EL.LACQ8T3D),
            (MT.QU8QU4, EL.LACQ8Q4D),
            (MT.QU8QU9, EL.LACQ8Q9D),
            (MT.TRIA33, EL.LACT3T3D),
            (MT.TR3TR6, EL.LACT3T6D),
            (MT.TR3QU4, EL.LACT3Q4D),
            (MT.TR3QU8, EL.LACT3Q8D),
            (MT.TR3QU9, EL.LACT3Q9D),
            (MT.TRIA66, EL.LACT6T6D),
            (MT.TR6TR3, EL.LACT6T3D),
            (MT.TR6QU4, EL.LACT6Q4D),
            (MT.TR6QU8, EL.LACT6Q8D),
            (MT.TR6QU9, EL.LACT6Q9D),
            (MT.QUAD99, EL.LACQ9Q9D),
            (MT.QU9TR3, EL.LACQ9T3D),
            (MT.QU9TR6, EL.LACQ9T6D),
            (MT.QU9QU4, EL.LACQ9Q4D),
            (MT.QU9QU8, EL.LACQ9Q8D),
        ),
    ),
)


phen.add(
    "CONT_LAC_EL_3DE",
    Modelisation(
        dim=(2, 3),
        code="L3E",
        attrs=((AT.CONTACT, "OUI"),),
        elements=(
            (MT.QUAD44, EL.LACQ4Q4E),
            (MT.QU4TR3, EL.LACQ4T3E),
            (MT.QU4TR6, EL.LACQ4T6E),
            (MT.QU4QU8, EL.LACQ4Q8E),
            (MT.QU4QU9, EL.LACQ4Q9E),
            (MT.QUAD88, EL.LACQ8Q8E),
            (MT.QU8TR6, EL.LACQ8T6E),
            (MT.QU8TR3, EL.LACQ8T3E),
            (MT.QU8QU4, EL.LACQ8Q4E),
            (MT.QU8QU9, EL.LACQ8Q9E),
        ),
    ),
)


# -- Define CONTACT elements for LAGRANGIAN AUGMENTED method (in MECA_NON_LINE) - Contact

phen.add(
    "CONT_LAG_EL_2D",
    Modelisation(
        dim=(1, 2),
        code="M2C",
        attrs=((AT.CONTACT, "OUI"),),
        elements=(
            (MT.SEG22, EL.CMS2S2),
            (MT.SEG33, EL.CMS3S3),
            (MT.SEG23, EL.CMS2S3),
            (MT.SEG32, EL.CMS3S2),
            (MT.POI1, EL.CMP1L2),
            (MT.POI1, EL.CMP1N2),
        ),
    ),
)

phen.add(
    "CONT_LAG_EL_2DA",
    Modelisation(
        dim=(1, 2),
        code="M2A",
        attrs=((AT.CONTACT, "OUI"), (AT.AXIS, "OUI")),
        elements=(
            (MT.SEG22, EL.CMS2S2A),
            (MT.SEG33, EL.CMS3S3A),
            (MT.SEG23, EL.CMS2S3A),
            (MT.SEG32, EL.CMS3S2A),
            (MT.POI1, EL.CMP1L2A),
            (MT.POI1, EL.CMP1N2A),
        ),
    ),
)

phen.add(
    "CONT_LAG_EL_3D",
    Modelisation(
        dim=(2, 3),
        code="M3C",
        attrs=((AT.CONTACT, "OUI"),),
        elements=(
            (MT.QU4QU8, EL.CMQ4Q8),
            (MT.QU4QU9, EL.CMQ4Q9),
            (MT.QU4TR3, EL.CMQ4T3),
            (MT.QU4TR6, EL.CMQ4T6),
            (MT.QU8QU4, EL.CMQ8Q4),
            (MT.QU8QU9, EL.CMQ8Q9),
            (MT.QU8TR3, EL.CMQ8T3),
            (MT.QU8TR6, EL.CMQ8T6),
            (MT.QU9QU4, EL.CMQ9Q4),
            (MT.QU9QU8, EL.CMQ9Q8),
            (MT.QU9TR3, EL.CMQ9T3),
            (MT.QU9TR6, EL.CMQ9T6),
            (MT.QUAD44, EL.CMQ4Q4),
            (MT.QUAD88, EL.CMQ8Q8),
            (MT.QUAD99, EL.CMQ9Q9),
            (MT.TR3QU4, EL.CMT3Q4),
            (MT.TR3QU8, EL.CMT3Q8),
            (MT.TR3QU9, EL.CMT3Q9),
            (MT.TR3TR6, EL.CMT3T6),
            (MT.TR6QU4, EL.CMT6Q4),
            (MT.TR6QU8, EL.CMT6Q8),
            (MT.TR6QU9, EL.CMT6Q9),
            (MT.TR6TR3, EL.CMT6T3),
            (MT.TRIA33, EL.CMT3T3),
            (MT.TRIA66, EL.CMT6T6),
            (MT.POI1, EL.CMP1L3),
            (MT.POI1, EL.CMP1N3),
        ),
    ),
)

phen.add(
    "FRIC_LAG_EL_2D",
    Modelisation(
        dim=(1, 2),
        code="M2F",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI")),
        elements=(
            (MT.SEG22, EL.FMS2S2),
            (MT.SEG33, EL.FMS3S3),
            (MT.SEG23, EL.FMS2S3),
            (MT.SEG32, EL.FMS3S2),
            (MT.POI1, EL.FMP1L2),
            (MT.POI1, EL.FMP1N2),
        ),
    ),
)

phen.add(
    "FRIC_LAG_EL_2DA",
    Modelisation(
        dim=(1, 2),
        code="MAF",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI"), (AT.AXIS, "OUI")),
        elements=(
            (MT.SEG22, EL.FMS2S2A),
            (MT.SEG33, EL.FMS3S3A),
            (MT.SEG23, EL.FMS2S3A),
            (MT.SEG32, EL.FMS3S2A),
            (MT.POI1, EL.FMP1L2A),
            (MT.POI1, EL.FMP1N2A),
        ),
    ),
)

phen.add(
    "FRIC_LAG_EL_3D",
    Modelisation(
        dim=(2, 3),
        code="M3F",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI")),
        elements=(
            (MT.QU4QU8, EL.FMQ4Q8),
            (MT.QU4QU9, EL.FMQ4Q9),
            (MT.QU4TR3, EL.FMQ4T3),
            (MT.QU4TR6, EL.FMQ4T6),
            (MT.QU8QU4, EL.FMQ8Q4),
            (MT.QU8QU9, EL.FMQ8Q9),
            (MT.QU8TR3, EL.FMQ8T3),
            (MT.QU8TR6, EL.FMQ8T6),
            (MT.QU9QU4, EL.FMQ9Q4),
            (MT.QU9QU8, EL.FMQ9Q8),
            (MT.QU9TR3, EL.FMQ9T3),
            (MT.QU9TR6, EL.FMQ9T6),
            (MT.QUAD44, EL.FMQ4Q4),
            (MT.QUAD88, EL.FMQ8Q8),
            (MT.QUAD99, EL.FMQ9Q9),
            (MT.TR3QU4, EL.FMT3Q4),
            (MT.TR3QU8, EL.FMT3Q8),
            (MT.TR3QU9, EL.FMT3Q9),
            (MT.TR3TR6, EL.FMT3T6),
            (MT.TR6QU4, EL.FMT6Q4),
            (MT.TR6QU8, EL.FMT6Q8),
            (MT.TR6QU9, EL.FMT6Q9),
            (MT.TR6TR3, EL.FMT6T3),
            (MT.TRIA33, EL.FMT3T3),
            (MT.TRIA66, EL.FMT6T6),
            (MT.POI1, EL.FMP1L3),
            (MT.POI1, EL.FMP1N3),
        ),
    ),
)

# -- Define CONTACT elements for PENALTY method (in MECA_NON_LINE) - Contact

phen.add(
    "CONT_PENA_EL_2D",
    Modelisation(
        dim=(1, 2),
        code="P2C",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI")),
        elements=(
            (MT.SEG22, EL.CPS2S2),
            (MT.SEG33, EL.CPS3S3),
            (MT.SEG23, EL.CPS2S3),
            (MT.SEG32, EL.CPS3S2),
            (MT.POI1, EL.CPP1L2),
            (MT.POI1, EL.CPP1N2),
        ),
    ),
)

phen.add(
    "CONT_PENA_EL_2DA",
    Modelisation(
        dim=(1, 2),
        code="P2A",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI"), (AT.AXIS, "OUI")),
        elements=(
            (MT.SEG22, EL.CPS2S2A),
            (MT.SEG33, EL.CPS3S3A),
            (MT.SEG23, EL.CPS2S3A),
            (MT.SEG32, EL.CPS3S2A),
            (MT.POI1, EL.CPP1L2A),
            (MT.POI1, EL.CPP1N2A),
        ),
    ),
)

phen.add(
    "CONT_PENA_EL_3D",
    Modelisation(
        dim=(2, 3),
        code="P3C",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI")),
        elements=(
            (MT.QU4QU8, EL.CPQ4Q8),
            (MT.QU4QU9, EL.CPQ4Q9),
            (MT.QU4TR3, EL.CPQ4T3),
            (MT.QU4TR6, EL.CPQ4T6),
            (MT.QU8QU4, EL.CPQ8Q4),
            (MT.QU8QU9, EL.CPQ8Q9),
            (MT.QU8TR3, EL.CPQ8T3),
            (MT.QU8TR6, EL.CPQ8T6),
            (MT.QU9QU4, EL.CPQ9Q4),
            (MT.QU9QU8, EL.CPQ9Q8),
            (MT.QU9TR3, EL.CPQ9T3),
            (MT.QU9TR6, EL.CPQ9T6),
            (MT.QUAD44, EL.CPQ4Q4),
            (MT.QUAD88, EL.CPQ8Q8),
            (MT.QUAD99, EL.CPQ9Q9),
            (MT.TR3QU4, EL.CPT3Q4),
            (MT.TR3QU8, EL.CPT3Q8),
            (MT.TR3QU9, EL.CPT3Q9),
            (MT.TR3TR6, EL.CPT3T6),
            (MT.TR6QU4, EL.CPT6Q4),
            (MT.TR6QU8, EL.CPT6Q8),
            (MT.TR6QU9, EL.CPT6Q9),
            (MT.TR6TR3, EL.CPT6T3),
            (MT.TRIA33, EL.CPT3T3),
            (MT.TRIA66, EL.CPT6T6),
            (MT.POI1, EL.CPP1L3),
            (MT.POI1, EL.CPP1N3),
        ),
    ),
)


phen.add(
    "FRIC_PENA_EL_2DA",
    Modelisation(
        dim=(1, 2),
        code="PAF",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI"), (AT.AXIS, "OUI")),
        elements=(
            (MT.SEG22, EL.CPS2S2A),
            (MT.SEG33, EL.CPS3S3A),
            (MT.SEG23, EL.CPS2S3A),
            (MT.SEG32, EL.CPS3S2A),
            (MT.POI1, EL.CPP1L2A),
            (MT.POI1, EL.CPP1N2A),
        ),
    ),
)

phen.add(
    "FRIC_PENA_EL_3D",
    Modelisation(
        dim=(2, 3),
        code="P3F",
        attrs=((AT.CONTACT, "OUI"), (AT.FROTTEMENT, "OUI")),
        elements=(
            (MT.QU4QU8, EL.CPQ4Q8),
            (MT.QU4QU9, EL.CPQ4Q9),
            (MT.QU4TR3, EL.CPQ4T3),
            (MT.QU4TR6, EL.CPQ4T6),
            (MT.QU8QU4, EL.CPQ8Q4),
            (MT.QU8QU9, EL.CPQ8Q9),
            (MT.QU8TR3, EL.CPQ8T3),
            (MT.QU8TR6, EL.CPQ8T6),
            (MT.QU9QU4, EL.CPQ9Q4),
            (MT.QU9QU8, EL.CPQ9Q8),
            (MT.QU9TR3, EL.CPQ9T3),
            (MT.QU9TR6, EL.CPQ9T6),
            (MT.QUAD44, EL.CPQ4Q4),
            (MT.QUAD88, EL.CPQ8Q8),
            (MT.QUAD99, EL.CPQ9Q9),
            (MT.TR3QU4, EL.CPT3Q4),
            (MT.TR3QU8, EL.CPT3Q8),
            (MT.TR3QU9, EL.CPT3Q9),
            (MT.TR3TR6, EL.CPT3T6),
            (MT.TR6QU4, EL.CPT6Q4),
            (MT.TR6QU8, EL.CPT6Q8),
            (MT.TR6QU9, EL.CPT6Q9),
            (MT.TR6TR3, EL.CPT6T3),
            (MT.TRIA33, EL.CPT3T3),
            (MT.TRIA66, EL.CPT6T6),
            (MT.POI1, EL.CPP1L3),
            (MT.POI1, EL.CPP1N3),
        ),
    ),
)

# -- Define CONTACT elements for NITSCHE method (in MECA_NON_LINE) - Contact

phen.add(
    "CONT_NIT_EL_2D",
    Modelisation(
        dim=(1, 2),
        code="N2C",
        attrs=((AT.CONTACT, "OUI"),),
        elements=(
            (MT.TR3SE2, EL.CNT3S2),
            (MT.TR3SE3, EL.CNT3S3),
            (MT.TR6SE2, EL.CNT6S2),
            (MT.TR6SE3, EL.CNT6S3),
            (MT.QU4SE2, EL.CNQ4S2),
            (MT.QU4SE3, EL.CNQ4S3),
            (MT.QU8SE2, EL.CNQ8S2),
            (MT.QU8SE3, EL.CNQ8S3),
            (MT.QU9SE2, EL.CNQ9S2),
            (MT.QU9SE3, EL.CNQ9S3),
        ),
    ),
)

phen.add(
    "CONT_NIT_EL_2DA",
    Modelisation(
        dim=(1, 2),
        code="N2A",
        attrs=((AT.CONTACT, "OUI"), (AT.AXIS, "OUI")),
        elements=(
            (MT.TR3SE2, EL.CNT3S2A),
            (MT.TR3SE3, EL.CNT3S3A),
            (MT.TR6SE2, EL.CNT6S2A),
            (MT.TR6SE3, EL.CNT6S3A),
            (MT.QU4SE2, EL.CNQ4S2A),
            (MT.QU4SE3, EL.CNQ4S3A),
            (MT.QU8SE2, EL.CNQ8S2A),
            (MT.QU8SE3, EL.CNQ8S3A),
            (MT.QU9SE2, EL.CNQ9S2A),
            (MT.QU9SE3, EL.CNQ9S3A),
        ),
    ),
)

# -- Define COUPLING FEM/HHO elements for PENALISATION method

phen.add(
    "CPL_PEN_H1_EL_2D",
    Modelisation(
        dim=(1, 2),
        code="HP2",
        attrs=(
            (AT.HHO, "OUI"),
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO"),
            (AT.RACCORD, "OUI"),
        ),
        elements=(
            (MT.SEG23, EL.CP_S2S3_HHO1),
            (MT.SEG33, EL.CP_S3S3_HHO1),
            (MT.SEG43, EL.CP_S4S3_HHO1),
        ),
    ),
)

phen.add(
    "CPL_PEN_H2_EL_2D",
    Modelisation(
        dim=(1, 2),
        code="HP2",
        attrs=(
            (AT.HHO, "OUI"),
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO"),
            (AT.RACCORD, "OUI"),
        ),
        elements=(
            (MT.SEG23, EL.CP_S2S3_HHO2),
            (MT.SEG33, EL.CP_S3S3_HHO2),
            (MT.SEG43, EL.CP_S4S3_HHO2),
        ),
    ),
)

phen.add(
    "CPL_PEN_H1_EL_3D",
    Modelisation(
        dim=(2, 3),
        code="H31",
        attrs=(
            (AT.HHO, "OUI"),
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO"),
            (AT.RACCORD, "OUI"),
        ),
        elements=(
            (MT.QU4QU9, EL.CP_Q4Q9_HHO1),
            (MT.QU8QU9, EL.CP_Q8Q9_HHO1),
            (MT.QUAD99, EL.CP_Q9Q9_HHO1),
            (MT.TR3QU9, EL.CP_T3Q9_HHO1),
            (MT.TR6QU9, EL.CP_T6Q9_HHO1),
            (MT.QU4TR7, EL.CP_Q4T7_HHO1),
            (MT.QU8TR7, EL.CP_Q8T7_HHO1),
            (MT.QU9TR7, EL.CP_Q9T7_HHO1),
            (MT.TR3TR7, EL.CP_T3T7_HHO1),
            (MT.TR6TR7, EL.CP_T6T9_HHO1),
        ),
    ),
)

phen.add(
    "CPL_PEN_H2_EL_3D",
    Modelisation(
        dim=(2, 3),
        code="H32",
        attrs=(
            (AT.HHO, "OUI"),
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO"),
            (AT.RACCORD, "OUI"),
        ),
        elements=(
            (MT.QU4QU9, EL.CP_Q4Q9_HHO2),
            (MT.QU8QU9, EL.CP_Q8Q9_HHO2),
            (MT.QUAD99, EL.CP_Q9Q9_HHO2),
            (MT.TR3QU9, EL.CP_T3Q9_HHO2),
            (MT.TR6QU9, EL.CP_T6Q9_HHO2),
            (MT.QU4TR7, EL.CP_Q4T7_HHO2),
            (MT.QU8TR7, EL.CP_Q8T7_HHO2),
            (MT.QU9TR7, EL.CP_Q9T7_HHO2),
            (MT.TR3TR7, EL.CP_T3T7_HHO2),
            (MT.TR6TR7, EL.CP_T6T9_HHO2),
        ),
    ),
)

# -- Define COUPLING FEM/FEM elements for PENALISATION method

phen.add(
    "CPL_PEN_EL_2D",
    Modelisation(
        dim=(1, 2),
        code="CP2",
        attrs=((AT.RACCORD, "OUI"),),
        elements=(
            (MT.SEG22, EL.CP_S2S2),
            (MT.SEG23, EL.CP_S2S3),
            (MT.SEG32, EL.CP_S3S2),
            (MT.SEG33, EL.CP_S3S3),
            (MT.SEG44, EL.CP_S4S4),
            (MT.SEG43, EL.CP_S4S3),
        ),
    ),
)

phen.add(
    "CPL_PEN_EL_3D",
    Modelisation(
        dim=(2, 3),
        code="CP3",
        attrs=((AT.RACCORD, "OUI"),),
        elements=(
            (MT.QU4QU8, EL.CP_Q4Q8),
            (MT.QU4QU9, EL.CP_Q4Q9),
            (MT.QU4TR3, EL.CP_Q4T3),
            (MT.QU4TR6, EL.CP_Q4T6),
            (MT.QU8QU4, EL.CP_Q8Q4),
            (MT.QU8QU9, EL.CP_Q8Q9),
            (MT.QU8TR3, EL.CP_Q8T3),
            (MT.QU8TR6, EL.CP_Q8T6),
            (MT.QU9QU4, EL.CP_Q9Q4),
            (MT.QU9QU8, EL.CP_Q9Q8),
            (MT.QU9TR3, EL.CP_Q9T3),
            (MT.QU9TR6, EL.CP_Q9T6),
            (MT.QUAD44, EL.CP_Q4Q4),
            (MT.QUAD88, EL.CP_Q8Q8),
            (MT.QUAD99, EL.CP_Q9Q9),
            (MT.TR3QU4, EL.CP_T3Q4),
            (MT.TR3QU8, EL.CP_T3Q8),
            (MT.TR3QU9, EL.CP_T3Q9),
            (MT.TR3TR6, EL.CP_T3T6),
            (MT.TR6QU4, EL.CP_T6Q4),
            (MT.TR6QU8, EL.CP_T6Q8),
            (MT.TR6QU9, EL.CP_T6Q9),
            (MT.TR6TR3, EL.CP_T6T3),
            (MT.TRIA33, EL.CP_T3T3),
            (MT.TRIA66, EL.CP_T6T6),
        ),
    ),
)

# -- Define COUPLING FEM/FEM elements for LAGRANGIAN method

phen.add(
    "CPL_LAG_EL_2D",
    Modelisation(
        dim=(1, 2),
        code="CL2",
        attrs=((AT.RACCORD, "OUI"),),
        elements=(
            (MT.SEG22, EL.CL_S2S2),
            (MT.SEG23, EL.CL_S2S3),
            (MT.SEG32, EL.CL_S3S2),
            (MT.SEG33, EL.CL_S3S3),
            (MT.SEG44, EL.CL_S4S4),
            (MT.SEG43, EL.CL_S4S3),
            (MT.POI1, EL.CL_POI2D),
        ),
    ),
)

phen.add(
    "CPL_LAG_EL_3D",
    Modelisation(
        dim=(2, 3),
        code="CL3",
        attrs=((AT.RACCORD, "OUI"),),
        elements=(
            (MT.QU4QU8, EL.CL_Q4Q8),
            (MT.QU4QU9, EL.CL_Q4Q9),
            (MT.QU4TR3, EL.CL_Q4T3),
            (MT.QU4TR6, EL.CL_Q4T6),
            (MT.QU8QU4, EL.CL_Q8Q4),
            (MT.QU8QU9, EL.CL_Q8Q9),
            (MT.QU8TR3, EL.CL_Q8T3),
            (MT.QU8TR6, EL.CL_Q8T6),
            (MT.QU9QU4, EL.CL_Q9Q4),
            (MT.QU9QU8, EL.CL_Q9Q8),
            (MT.QU9TR3, EL.CL_Q9T3),
            (MT.QU9TR6, EL.CL_Q9T6),
            (MT.QUAD44, EL.CL_Q4Q4),
            (MT.QUAD88, EL.CL_Q8Q8),
            (MT.QUAD99, EL.CL_Q9Q9),
            (MT.TR3QU4, EL.CL_T3Q4),
            (MT.TR3QU8, EL.CL_T3Q8),
            (MT.TR3QU9, EL.CL_T3Q9),
            (MT.TR3TR6, EL.CL_T3T6),
            (MT.TR6QU4, EL.CL_T6Q4),
            (MT.TR6QU8, EL.CL_T6Q8),
            (MT.TR6QU9, EL.CL_T6Q9),
            (MT.TR6TR3, EL.CL_T6T3),
            (MT.TRIA33, EL.CL_T3T3),
            (MT.TRIA66, EL.CL_T6T6),
            (MT.POI1, EL.CL_POI3D),
        ),
    ),
)

# -- Define COUPLING FEM/HHO elements for NITSCHE method

phen.add(
    "CPL_NIT_H1_EL_2D",
    Modelisation(
        dim=(1, 2),
        code="HN1",
        attrs=(
            (AT.HHO, "OUI"),
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO"),
            (AT.RACCORD, "OUI"),
        ),
        elements=(
            (MT.TR3SE3, EL.CN_T3S3_HHO1),
            (MT.TR6SE3, EL.CN_T6S3_HHO1),
            (MT.T10SE3, EL.CN_T10S3_HHO1),
            (MT.QU4SE3, EL.CN_Q4S3_HHO1),
            (MT.QU8SE3, EL.CN_Q8S3_HHO1),
            (MT.QU9SE3, EL.CN_Q9S3_HHO1),
            (MT.Q12SE3, EL.CN_Q12S3_HHO1),
        ),
    ),
)

phen.add(
    "CPL_NIT_H2_EL_2D",
    Modelisation(
        dim=(1, 2),
        code="HN2",
        attrs=(
            (AT.HHO, "OUI"),
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO"),
            (AT.RACCORD, "OUI"),
        ),
        elements=(
            (MT.TR3SE3, EL.CN_T3S3_HHO2),
            (MT.TR6SE3, EL.CN_T6S3_HHO2),
            (MT.T10SE3, EL.CN_T10S3_HHO2),
            (MT.QU4SE3, EL.CN_Q4S3_HHO2),
            (MT.QU8SE3, EL.CN_Q8S3_HHO2),
            (MT.QU9SE3, EL.CN_Q9S3_HHO2),
            (MT.Q12SE3, EL.CN_Q12S3_HHO2),
        ),
    ),
)

phen.add(
    "CPL_NIT_H1_EL_3D",
    Modelisation(
        dim=(2, 3),
        code="HN6",
        attrs=(
            (AT.HHO, "OUI"),
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO"),
            (AT.RACCORD, "OUI"),
        ),
        elements=(
            (MT.TE4TR7, EL.CN_TE4TR7_HHO1),
            (MT.T10TR7, EL.CN_T10TR7_HHO1),
            (MT.T20TR7, EL.CN_T20TR7_HHO1),
        ),
    ),
)

phen.add(
    "CPL_NIT_H2_EL_3D",
    Modelisation(
        dim=(2, 3),
        code="HN7",
        attrs=(
            (AT.HHO, "OUI"),
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO"),
            (AT.RACCORD, "OUI"),
        ),
        elements=(
            (MT.TE4TR7, EL.CN_TE4TR7_HHO2),
            (MT.T10TR7, EL.CN_T10TR7_HHO2),
            (MT.T20TR7, EL.CN_T20TR7_HHO2),
        ),
    ),
)

# ------------------------------------------------------------------------------------
# Modelisations sous-terraines pour :
#  * Forces nodales
# ------------------------------------------------------------------------------------

phen.add("CL_FNOD2", Modelisation(dim=(0, 2), code="CL3", elements=((MT.POI1, EL.FORCE_NOD_2DDL),)))

phen.add("CL_FNOD3", Modelisation(dim=(0, 3), code="CL4", elements=((MT.POI1, EL.FORCE_NOD_3DDL),)))

phen.add("CL_FNOD6", Modelisation(dim=(0, 3), code="CL5", elements=((MT.POI1, EL.FORCE_NOD_6DDL),)))

phen.add(
    "CL_FNODCQ2", Modelisation(dim=(0, 2), code="CL6", elements=((MT.POI1, EL.FORCE_NOD_COQ2D),))
)

# ------------------------------------------------------------------------------------

phen.add(
    "COQUE_3D",
    Modelisation(
        dim=(2, 3),
        code="CQ3",
        attrs=(
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.COQUE, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.SOUS_POINT, "OUI"),
        ),
        elements=((MT.QUAD9, EL.MEC3QU9H), (MT.TRIA7, EL.MEC3TR7H), (MT.SEG3, EL.MEBOCQ3)),
    ),
)

phen.add(
    "COQUE_AXIS",
    Modelisation(
        dim=(1, 2),
        code="CQA",
        attrs=((AT.AXIS, "OUI"), (AT.TYPMOD, "C_PLAN"), (AT.COQUE, "OUI"), (AT.EFGE, "OUI")),
        elements=((MT.SEG3, EL.MECXSE3),),
    ),
)

phen.add(
    "C_PLAN",
    Modelisation(
        dim=(2, 2),
        code="CPL",
        attrs=((AT.NBSIGM, "4"), (AT.C_PLAN, "OUI"), (AT.TYPMOD, "C_PLAN")),
        elements=(
            (MT.TRIA3, EL.MECPTR3),
            (MT.QUAD4, EL.MECPQU4),
            (MT.TRIA6, EL.MECPTR6),
            (MT.QUAD8, EL.MECPQU8),
            (MT.QUAD9, EL.MECPQU9),
            (MT.SEG2, EL.MEPLSE2),
            (MT.SEG3, EL.MEPLSE3),
        ),
    ),
)

phen.add(
    "C_PLAN2XH",
    Modelisation(
        dim=(2, 2),
        code="CX4",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.C_PLAN, "OUI"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=((MT.TRIA6, EL.MECPTR6_XH), (MT.QUAD8, EL.MECPQU8_XH), (MT.SEG3, EL.MEPLSE3_XH)),
    ),
)

phen.add(
    "C_PLAN2XHT",
    Modelisation(
        dim=(2, 2),
        code="CX6",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.C_PLAN, "OUI"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XHT"),
        ),
        elements=(
            (MT.TRIA6, EL.MECPTR6_XHT),
            (MT.QUAD8, EL.MECPQU8_XHT),
            (MT.SEG3, EL.MEPLSE3_XHT),
        ),
    ),
)

phen.add(
    "C_PLAN2XT",
    Modelisation(
        dim=(2, 2),
        code="CX5",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.C_PLAN, "OUI"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XT"),
        ),
        elements=((MT.TRIA6, EL.MECPTR6_XT), (MT.QUAD8, EL.MECPQU8_XT), (MT.SEG3, EL.MEPLSE3_XT)),
    ),
)

phen.add(
    "C_PLAN_SI",
    Modelisation(
        dim=(2, 2),
        code="CPS",
        attrs=((AT.NBSIGM, "4"), (AT.C_PLAN, "OUI"), (AT.TYPMOD, "C_PLAN")),
        elements=(
            (MT.QUAD8, EL.MECPQS8),
            (MT.QUAD4, EL.MECPQS4),
            (MT.SEG3, EL.MEPLSE3),
            (MT.SEG2, EL.MEPLSE2),
        ),
    ),
)

phen.add(
    "C_PLAN_XH",
    Modelisation(
        dim=(2, 2),
        code="CX1",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.C_PLAN, "OUI"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=(
            (MT.TRIA6, EL.MECPTR6_XH),
            (MT.TRIA3, EL.MECPTR3_XH),
            (MT.QUAD8, EL.MECPQU8_XH),
            (MT.QUAD4, EL.MECPQU4_XH),
            (MT.SEG3, EL.MEPLSE3_XH),
            (MT.SEG2, EL.MEPLSE2_XH),
        ),
    ),
)

phen.add(
    "C_PLAN_XH1",
    Modelisation(
        dim=(2, 2),
        code="CXA",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.C_PLAN, "OUI"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH1"),
        ),
        elements=((MT.TRIA3, EL.MECPTR3_XH1), (MT.QUAD4, EL.MECPQU4_XH1)),
    ),
)

phen.add(
    "C_PLAN_XH2",
    Modelisation(
        dim=(2, 2),
        code="CXB",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.C_PLAN, "OUI"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH2"),
        ),
        elements=((MT.TRIA3, EL.MECPTR3_XH2), (MT.QUAD4, EL.MECPQU4_XH2)),
    ),
)

phen.add(
    "C_PLAN_XH3",
    Modelisation(
        dim=(2, 2),
        code="CXC",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.C_PLAN, "OUI"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH3"),
        ),
        elements=((MT.TRIA3, EL.MECPTR3_XH3), (MT.QUAD4, EL.MECPQU4_XH3)),
    ),
)

phen.add(
    "C_PLAN_XH4",
    Modelisation(
        dim=(2, 2),
        code="CXD",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.C_PLAN, "OUI"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH4"),
        ),
        elements=((MT.TRIA3, EL.MECPTR3_XH4), (MT.QUAD4, EL.MECPQU4_XH4)),
    ),
)

phen.add(
    "C_PLAN_XHT",
    Modelisation(
        dim=(2, 2),
        code="CX3",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.C_PLAN, "OUI"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XHT"),
        ),
        elements=(
            (MT.TRIA3, EL.MECPTR3_XHT),
            (MT.QUAD4, EL.MECPQU4_XHT),
            (MT.SEG3, EL.MEPLSE2_XHT),
        ),
    ),
)

phen.add(
    "C_PLAN_XT",
    Modelisation(
        dim=(2, 2),
        code="CX2",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.C_PLAN, "OUI"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XT"),
        ),
        elements=((MT.TRIA3, EL.MECPTR3_XT), (MT.QUAD4, EL.MECPQU4_XT), (MT.SEG3, EL.MEPLSE2_XT)),
    ),
)

phen.add(
    "DIS_T",
    Modelisation(
        dim=(-1, 3),
        code="DIT",
        attrs=((AT.TYPMOD, "0D"), (AT.EFGE, "OUI"), (AT.SIGM, "NON")),
        elements=((MT.SEG2, EL.MECA_DIS_T_L), (MT.POI1, EL.MECA_DIS_T_N)),
    ),
)

phen.add(
    "DIS_TR",
    Modelisation(
        dim=(-1, 3),
        code="DTR",
        attrs=((AT.TYPMOD, "0D"), (AT.EFGE, "OUI"), (AT.SIGM, "NON")),
        elements=((MT.SEG2, EL.MECA_DIS_TR_L), (MT.POI1, EL.MECA_DIS_TR_N)),
    ),
)

phen.add(
    "DKT",
    Modelisation(
        dim=(2, 3),
        code="DKT",
        attrs=(
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.COQUE, "OUI"),
            (AT.PLAQUE, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.SOUS_POINT, "OUI"),
        ),
        elements=((MT.TRIA3, EL.MEDKTR3), (MT.QUAD4, EL.MEDKQU4), (MT.SEG2, EL.MEBODKT)),
    ),
)

phen.add(
    "DKTG",
    Modelisation(
        dim=(2, 3),
        code="DTG",
        attrs=(
            (AT.NBSIGM, "6"),
            (AT.COQUE, "OUI"),
            (AT.PLAQUE, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
        ),
        elements=((MT.TRIA3, EL.MEDKTG3), (MT.QUAD4, EL.MEDKQG4), (MT.SEG2, EL.MEBODKT)),
    ),
)

phen.add(
    "DST",
    Modelisation(
        dim=(2, 3),
        code="DST",
        attrs=(
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.COQUE, "OUI"),
            (AT.PLAQUE, "OUI"),
            (AT.EFGE, "OUI"),
        ),
        elements=((MT.TRIA3, EL.MEDSTR3), (MT.QUAD4, EL.MEDSQU4), (MT.SEG2, EL.MEBODST)),
    ),
)

phen.add(
    "D_PLAN",
    Modelisation(
        dim=(2, 2),
        code="DPL",
        attrs=((AT.NBSIGM, "4"), (AT.D_PLAN, "OUI"), (AT.TYPMOD, "D_PLAN")),
        elements=(
            (MT.TRIA3, EL.MEDPTR3),
            (MT.TRIA6, EL.MEDPTR6),
            (MT.TRIA10, EL.MEDPTR10),
            (MT.QUAD4, EL.MEDPQU4),
            (MT.QUAD8, EL.MEDPQU8),
            (MT.QUAD9, EL.MEDPQU9),
            (MT.QUAD12, EL.MEDPQU12),
            (MT.SEG2, EL.MEPLSE2),
            (MT.SEG3, EL.MEPLSE3),
            (MT.SEG4, EL.MEPLSE4),
        ),
    ),
)

phen.add(
    "D_PLAN2XH",
    Modelisation(
        dim=(2, 2),
        code="DX4",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=((MT.TRIA6, EL.MEDPTR6_XH), (MT.QUAD8, EL.MEDPQU8_XH), (MT.SEG3, EL.MEPLSE3_XH)),
    ),
)

phen.add(
    "D_PLAN2XHT",
    Modelisation(
        dim=(2, 2),
        code="DX6",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XHT"),
        ),
        elements=(
            (MT.TRIA6, EL.MEDPTR6_XHT),
            (MT.QUAD8, EL.MEDPQU8_XHT),
            (MT.SEG3, EL.MEPLSE3_XHT),
        ),
    ),
)

phen.add(
    "D_PLAN2XT",
    Modelisation(
        dim=(2, 2),
        code="DX5",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XT"),
        ),
        elements=((MT.TRIA6, EL.MEDPTR6_XT), (MT.QUAD8, EL.MEDPQU8_XT), (MT.SEG3, EL.MEPLSE3_XT)),
    ),
)

phen.add(
    "D_PLAN_ABSO",
    Modelisation(
        dim=(1, 2),
        code="DPA",
        attrs=((AT.TYPMOD, "D_PLAN"), (AT.FLUIDE, "NON"), (AT.ABSO, "OUI")),
        elements=((MT.SEG2, EL.MEPASE2), (MT.SEG3, EL.MEPASE3)),
    ),
)

phen.add(
    "D_PLAN_DIL#1",
    Modelisation(
        dim=(2, 2),
        code="D2D",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.FORMULATION, "DIL"),
        ),
        elements=((MT.TRIA6, EL.TR6_DP_2D), (MT.QUAD8, EL.QU8_DP_2D), (MT.SEG3, EL.MEPLSE3)),
    ),
)

phen.add(
    "D_PLAN_DIL#2",
    Modelisation(
        dim=(2, 2),
        code="D2I",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.FORMULATION, "DIL_INCO"),
        ),
        elements=((MT.TRIA6, EL.TR6_DP_2DI), (MT.QUAD8, EL.QU8_DP_2DI), (MT.SEG3, EL.MEPLSE3)),
    ),
)

phen.add(
    "D_PLAN_GRAD_SIGM",
    Modelisation(
        dim=(2, 2),
        code="DSG",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "GRADSIGM"),
        ),
        elements=(
            (MT.TRIA6, EL.MGSPTR6),
            (MT.QUAD8, EL.MGSPQU8),
            (MT.SEG2, EL.MEPLSE2),
            (MT.SEG3, EL.MEPLSE3),
        ),
    ),
)

phen.add(
    "D_PLAN_GRAD_VARI",
    Modelisation(
        dim=(2, 2),
        code="DPV",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "GRADVARI"),
        ),
        elements=((MT.SEG3, EL.MEPLSE3), (MT.TRIA6, EL.MVDPTR6), (MT.QUAD8, EL.MVDPQS8)),
    ),
)

phen.add(
    "D_PLAN_GRAD_INCO",
    Modelisation(
        dim=(2, 2),
        code="DGI",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.INCO, "C5GV"),
            (AT.TYPMOD2, "GRADVARI"),
        ),
        elements=((MT.SEG3, EL.MEPLSE3), (MT.TRIA6, EL.GVI_DP_TR6), (MT.QUAD8, EL.GVI_DP_QU8)),
    ),
)

phen.add(
    "D_PLAN_HH2D",
    Modelisation(
        dim=(2, 2),
        code="DZ4",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=((MT.QUAD8, EL.HH2_DPQ8D), (MT.TRIA6, EL.HH2_DPTR6D), (MT.SEG3, EL.HH2_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HH2MD",
    Modelisation(
        dim=(2, 2),
        code="DA3",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=((MT.QUAD8, EL.HH2M_DPQ8D), (MT.TRIA6, EL.HH2M_DPTR6D), (MT.SEG3, EL.HH2M_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HH2MS",
    Modelisation(
        dim=(2, 2),
        code="DR3",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.HH2M_DPQ8S), (MT.TRIA6, EL.HH2M_DPTR6S), (MT.SEG3, EL.HH2M_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HH2MS_DIL",
    Modelisation(
        dim=(2, 2),
        code="DD1",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.DIL, "OUI"),
            (AT.FORMULATION, "DIL"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.HH2M_DPQ8S_DIL),
            (MT.TRIA6, EL.HH2M_DPTR6S_DIL),
            (MT.SEG3, EL.HH2M_DPSE3),
        ),
    ),
)

phen.add(
    "D_PLAN_HH2M_SI",
    Modelisation(
        dim=(2, 2),
        code="DM2",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.D_PLAN, "OUI"),
        ),
        elements=((MT.QUAD8, EL.HH2M_DPQ8M), (MT.TRIA6, EL.HH2M_DPTR6M), (MT.SEG3, EL.HH2M_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HH2S",
    Modelisation(
        dim=(2, 2),
        code="DZ3",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.HH2_DPQ8S), (MT.TRIA6, EL.HH2_DPTR6S), (MT.SEG3, EL.HH2_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HH2SUDA",
    Modelisation(
        dim=(2, 2),
        code="2DA",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.TYPMOD3, "SUSHI"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.D_PLAN, "OUI"),
        ),
        elements=(
            (MT.TRIA7, EL.DHH2T7_SUDA),
            (MT.QUAD9, EL.DHH2Q9_SUDA),
            (MT.SEG3, EL.DHH2S3_SUDA),
        ),
    ),
)

phen.add(
    "D_PLAN_HHD",
    Modelisation(
        dim=(2, 2),
        code="DZ2",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=((MT.QUAD8, EL.HH_DPQ8D), (MT.TRIA6, EL.HH_DPTR6D), (MT.SEG3, EL.HH_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HHM",
    Modelisation(
        dim=(2, 2),
        code="DH1",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.D_PLAN, "OUI"),
        ),
        elements=((MT.QUAD8, EL.HHM_DPQ8), (MT.TRIA6, EL.HHM_DPTR6), (MT.SEG3, EL.HHM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HHMD",
    Modelisation(
        dim=(2, 2),
        code="DH6",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=((MT.QUAD8, EL.HHM_DPQ8D), (MT.TRIA6, EL.HHM_DPTR6D), (MT.SEG3, EL.HHM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HHMS",
    Modelisation(
        dim=(2, 2),
        code="DR2",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.HHM_DPQ8S), (MT.TRIA6, EL.HHM_DPTR6S), (MT.SEG3, EL.HHM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HHS",
    Modelisation(
        dim=(2, 2),
        code="DZ1",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.HH_DPQ8S), (MT.TRIA6, EL.HH_DPTR6S), (MT.SEG3, EL.HH_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HM",
    Modelisation(
        dim=(2, 2),
        code="DH2",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
        ),
        elements=((MT.QUAD8, EL.HM_DPQ8), (MT.TRIA6, EL.HM_DPTR6), (MT.SEG3, EL.HM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HMD",
    Modelisation(
        dim=(2, 2),
        code="DH7",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=((MT.QUAD8, EL.HM_DPQ8D), (MT.TRIA6, EL.HM_DPTR6D), (MT.SEG3, EL.HM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HMS",
    Modelisation(
        dim=(2, 2),
        code="DR1",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.HM_DPQ8S), (MT.TRIA6, EL.HM_DPTR6S), (MT.SEG3, EL.HM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HMS_DIL",
    Modelisation(
        dim=(2, 2),
        code="DD2",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.DIL, "OUI"),
            (AT.FORMULATION, "DIL"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_DPQ8S_DIL),
            (MT.TRIA6, EL.HM_DPTR6S_DIL),
            (MT.SEG3, EL.HM_DPSE3),
        ),
    ),
)


phen.add(
    "D_PLAN_HM_SI",
    Modelisation(
        dim=(2, 2),
        code="DM1",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
        ),
        elements=((MT.QUAD8, EL.HM_DPQ8M), (MT.TRIA6, EL.HM_DPTR6M), (MT.SEG3, EL.HM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_HM_SI_DIL",
    Modelisation(
        dim=(2, 2),
        code="DD3",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.DIL, "OUI"),
            (AT.FORMULATION, "DIL"),
            (AT.D_PLAN, "OUI"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_DPQ8M_DIL),
            (MT.TRIA6, EL.HM_DPTR6M_DIL),
            (MT.SEG3, EL.HM_DPSE3),
        ),
    ),
)

phen.add(
    "D_PLAN_HM_XH",
    Modelisation(
        dim=(2, 2),
        code="DXI",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=((MT.QUAD8, EL.HM_DPQ8_XH), (MT.TRIA6, EL.HM_DPTR6_XH), (MT.SEG3, EL.HM_DPSE3_XH)),
    ),
)

phen.add(
    "D_PLAN_HM_XH1",
    Modelisation(
        dim=(2, 2),
        code="DXM",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH1"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_DPQ8_XH1),
            (MT.TRIA6, EL.HM_DPTR6_XH1),
            (MT.SEG3, EL.HM_DPSE3_XH1),
        ),
    ),
)

phen.add(
    "D_PLAN_HM_XH2",
    Modelisation(
        dim=(2, 2),
        code="DXN",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH2"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_DPQ8_XH2),
            (MT.TRIA6, EL.HM_DPTR6_XH2),
            (MT.SEG3, EL.HM_DPSE3_XH2),
        ),
    ),
)

phen.add(
    "D_PLAN_HM_XH3",
    Modelisation(
        dim=(2, 2),
        code="DXO",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH3"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_DPQ8_XH3),
            (MT.TRIA6, EL.HM_DPTR6_XH3),
            (MT.SEG3, EL.HM_DPSE3_XH3),
        ),
    ),
)

phen.add(
    "D_PLAN_HM_XH_D",
    Modelisation(
        dim=(2, 2),
        code="DXK",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_DPQ8D_XH),
            (MT.TRIA6, EL.HM_DPTR6D_XH),
            (MT.SEG3, EL.HM_DPSE3_XH),
        ),
    ),
)

phen.add(
    "D_PLAN_HM_XH_S",
    Modelisation(
        dim=(2, 2),
        code="DXL",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_DPQ8S_XH),
            (MT.TRIA6, EL.HM_DPTR6S_XH),
            (MT.SEG3, EL.HM_DPSE3_XH),
        ),
    ),
)

phen.add(
    "D_PLAN_HM_XH_SI",
    Modelisation(
        dim=(2, 2),
        code="DXJ",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "XFEM_HM"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=(
            (MT.QUAD8, EL.HM_DPQ8M_XH),
            (MT.TRIA6, EL.HM_DPTR6M_XH),
            (MT.SEG3, EL.HM_DPSE3_XH),
        ),
    ),
)

phen.add(
    "D_PLAN_HS",
    Modelisation(
        dim=(2, 2),
        code="DHA",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.H_DPQ8S), (MT.TRIA6, EL.H_DPTR6S), (MT.SEG3, EL.H_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_INCO_UP",
    Modelisation(
        dim=(2, 2),
        code="PUP",
        attrs=((AT.NBSIGM, "4"), (AT.D_PLAN, "OUI"), (AT.TYPMOD, "D_PLAN"), (AT.INCO, "C2")),
        elements=(
            (MT.TRIA6, EL.MUPLTR6),
            (MT.TRIA3, EL.MUPLTR3),
            (MT.QUAD8, EL.MUPLQU8),
            (MT.SEG3, EL.MEPLSE3),
            (MT.SEG2, EL.MEPLSE2),
        ),
    ),
)

phen.add(
    "D_PLAN_INCO_UPG",
    Modelisation(
        dim=(2, 2),
        code="PLI",
        attrs=((AT.NBSIGM, "4"), (AT.D_PLAN, "OUI"), (AT.TYPMOD, "D_PLAN"), (AT.INCO, "C3")),
        elements=((MT.TRIA6, EL.MIPLTR6), (MT.QUAD8, EL.MIPLQU8), (MT.SEG3, EL.MEPLSE3)),
    ),
)

phen.add(
    "D_PLAN_INCO_UPO",
    Modelisation(
        dim=(2, 2),
        code="POS",
        attrs=((AT.NBSIGM, "4"), (AT.D_PLAN, "OUI"), (AT.TYPMOD, "D_PLAN"), (AT.INCO, "C2O")),
        elements=(
            (MT.TRIA3, EL.MIPLOSTR3),
            (MT.QUAD4, EL.MIPLOSQU4),
            (MT.SEG3, EL.MEPLSE3),
            (MT.SEG2, EL.MEPLSE2),
        ),
    ),
)

phen.add(
    "D_PLAN_SI",
    Modelisation(
        dim=(2, 2),
        code="DPS",
        attrs=((AT.NBSIGM, "4"), (AT.D_PLAN, "OUI"), (AT.TYPMOD, "D_PLAN")),
        elements=(
            (MT.QUAD8, EL.MEDPQS8),
            (MT.QUAD4, EL.MEDPQS4),
            (MT.SEG3, EL.MEPLSE3),
            (MT.SEG2, EL.MEPLSE2),
        ),
    ),
)

phen.add(
    "D_PLAN_THH2D",
    Modelisation(
        dim=(2, 2),
        code="DA1",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=((MT.QUAD8, EL.THH2_DPQ8D), (MT.TRIA6, EL.THH2_DPTR6D), (MT.SEG3, EL.THH2_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_THH2MD",
    Modelisation(
        dim=(2, 2),
        code="DA2",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=(
            (MT.QUAD8, EL.THH2M_DPQ8D),
            (MT.TRIA6, EL.THH2M_DPTR6D),
            (MT.SEG3, EL.THH2M_DPSE3),
        ),
    ),
)

phen.add(
    "D_PLAN_THH2MS",
    Modelisation(
        dim=(2, 2),
        code="DR7",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.THH2M_DPQ8S),
            (MT.TRIA6, EL.THH2M_DPTR6S),
            (MT.SEG3, EL.THH2M_DPSE3),
        ),
    ),
)

phen.add(
    "D_PLAN_THH2S",
    Modelisation(
        dim=(2, 2),
        code="DR5",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "2"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.THH2_DPQ8S), (MT.TRIA6, EL.THH2_DPTR6S), (MT.SEG3, EL.THH2_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_THHD",
    Modelisation(
        dim=(2, 2),
        code="DH8",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=((MT.QUAD8, EL.THH_DPQ8D), (MT.TRIA6, EL.THH_DPTR6D), (MT.SEG3, EL.THH_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_THHMD",
    Modelisation(
        dim=(2, 2),
        code="DH9",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=((MT.QUAD8, EL.THHM_DPQ8D), (MT.TRIA6, EL.THHM_DPTR6D), (MT.SEG3, EL.THHM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_THHMS",
    Modelisation(
        dim=(2, 2),
        code="DR6",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.THHM_DPQ8S), (MT.TRIA6, EL.THHM_DPTR6S), (MT.SEG3, EL.THHM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_THHS",
    Modelisation(
        dim=(2, 2),
        code="DR4",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "1"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.THH_DPQ8S), (MT.TRIA6, EL.THH_DPTR6S), (MT.SEG3, EL.THH_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_THM",
    Modelisation(
        dim=(2, 2),
        code="DH5",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
        ),
        elements=((MT.QUAD8, EL.THM_DPQ8), (MT.TRIA6, EL.THM_DPTR6), (MT.SEG3, EL.THM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_THMD",
    Modelisation(
        dim=(2, 2),
        code="DH0",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=((MT.QUAD8, EL.THM_DPQ8D), (MT.TRIA6, EL.THM_DPTR6D), (MT.SEG3, EL.THM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_THMS",
    Modelisation(
        dim=(2, 2),
        code="DR8",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.THM_DPQ8S), (MT.TRIA6, EL.THM_DPTR6S), (MT.SEG3, EL.THM_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_THVD",
    Modelisation(
        dim=(2, 2),
        code="DG3",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "LUM"),
        ),
        elements=((MT.QUAD8, EL.THV_DPQ8D), (MT.TRIA6, EL.THV_DPTR6D), (MT.SEG3, EL.THV_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_THMS_DIL",
    Modelisation(
        dim=(2, 2),
        code="DD4",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "OUI"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.DIL, "OUI"),
            (AT.FORMULATION, "DIL"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=(
            (MT.QUAD8, EL.THM_DPQ8S_DIL),
            (MT.TRIA6, EL.THM_DPTR6S_DIL),
            (MT.SEG3, EL.THM_DPSE3),
        ),
    ),
)

phen.add(
    "D_PLAN_THVS",
    Modelisation(
        dim=(2, 2),
        code="DR9",
        attrs=(
            (AT.TYPMOD, "D_PLAN"),
            (AT.TYPMOD2, "THM"),
            (AT.MECA, "NON"),
            (AT.THER, "OUI"),
            (AT.HYDR1, "2"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.THV_DPQ8S), (MT.TRIA6, EL.THV_DPTR6S), (MT.SEG3, EL.THV_DPSE3)),
    ),
)

phen.add(
    "D_PLAN_XH",
    Modelisation(
        dim=(2, 2),
        code="DX1",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH"),
        ),
        elements=(
            (MT.TRIA6, EL.MEDPTR6_XH),
            (MT.TRIA3, EL.MEDPTR3_XH),
            (MT.QUAD8, EL.MEDPQU8_XH),
            (MT.QUAD4, EL.MEDPQU4_XH),
            (MT.SEG3, EL.MEPLSE3_XH),
            (MT.SEG2, EL.MEPLSE2_XH),
        ),
    ),
)

phen.add(
    "D_PLAN_XH1",
    Modelisation(
        dim=(2, 2),
        code="DXA",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH1"),
        ),
        elements=(
            (MT.TRIA3, EL.MEDPTR3_XH1),
            (MT.QUAD4, EL.MEDPQU4_XH1),
            (MT.SEG2, EL.MEPLSE2_XH1),
        ),
    ),
)

phen.add(
    "D_PLAN_XH2",
    Modelisation(
        dim=(2, 2),
        code="DXB",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH2"),
        ),
        elements=(
            (MT.TRIA3, EL.MEDPTR3_XH2),
            (MT.QUAD4, EL.MEDPQU4_XH2),
            (MT.SEG2, EL.MEPLSE2_XH2),
        ),
    ),
)

phen.add(
    "D_PLAN_XH3",
    Modelisation(
        dim=(2, 2),
        code="DXC",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH3"),
        ),
        elements=(
            (MT.TRIA3, EL.MEDPTR3_XH3),
            (MT.QUAD4, EL.MEDPQU4_XH3),
            (MT.SEG2, EL.MEPLSE2_XH3),
        ),
    ),
)

phen.add(
    "D_PLAN_XH4",
    Modelisation(
        dim=(2, 2),
        code="DXD",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XH4"),
        ),
        elements=(
            (MT.TRIA3, EL.MEDPTR3_XH4),
            (MT.QUAD4, EL.MEDPQU4_XH4),
            (MT.SEG2, EL.MEPLSE2_XH4),
        ),
    ),
)

phen.add(
    "D_PLAN_XHT",
    Modelisation(
        dim=(2, 2),
        code="DX3",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XHT"),
        ),
        elements=(
            (MT.TRIA3, EL.MEDPTR3_XHT),
            (MT.QUAD4, EL.MEDPQU4_XHT),
            (MT.SEG3, EL.MEPLSE2_XHT),
        ),
    ),
)

phen.add(
    "D_PLAN_XT",
    Modelisation(
        dim=(2, 2),
        code="DX2",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.LXFEM, "OUI"),
            (AT.XFEM, "XT"),
        ),
        elements=((MT.TRIA3, EL.MEDPTR3_XT), (MT.QUAD4, EL.MEDPQU4_XT), (MT.SEG3, EL.MEPLSE2_XT)),
    ),
)

phen.add(
    "FLUI_STRU#1",
    Modelisation(
        dim=(3, 3),
        code="FLS",
        attrs=((AT.TYPMOD, "3D"), (AT.FLUIDE, "OUI"), (AT.FSI, "OUI"), (AT.FORMULATION, "U_P_PHI")),
        elements=(
            (MT.TRIA3, EL.MEFS_FACE3),
            (MT.QUAD4, EL.MEFS_FACE4),
            (MT.TRIA6, EL.MEFS_FACE6),
            (MT.QUAD8, EL.MEFS_FACE8),
            (MT.QUAD9, EL.MEFS_FACE9),
        ),
    ),
)

phen.add(
    "FLUI_STRU#2",
    Modelisation(
        dim=(3, 3),
        code="FLP",
        attrs=((AT.TYPMOD, "3D"), (AT.FLUIDE, "OUI"), (AT.FSI, "OUI"), (AT.FORMULATION, "U_P")),
        elements=(
            (MT.TRIA3, EL.MEFS_FACE3P),
            (MT.QUAD4, EL.MEFS_FACE4P),
            (MT.TRIA6, EL.MEFS_FACE6P),
            (MT.QUAD8, EL.MEFS_FACE8P),
            (MT.QUAD9, EL.MEFS_FACE9P),
        ),
    ),
)

phen.add(
    "FLUI_STRU#3",
    Modelisation(
        dim=(3, 3),
        code="FLI",
        attrs=((AT.TYPMOD, "3D"), (AT.FLUIDE, "OUI"), (AT.FSI, "OUI"), (AT.FORMULATION, "U_PSI")),
        elements=(
            (MT.TRIA3, EL.MEFS_FACE3PSI),
            (MT.QUAD4, EL.MEFS_FACE4PSI),
            (MT.TRIA6, EL.MEFS_FACE6PSI),
            (MT.QUAD8, EL.MEFS_FACE8PSI),
            (MT.QUAD9, EL.MEFS_FACE9PSI),
        ),
    ),
)

phen.add(
    "GRILLE_EXCENTRE",
    Modelisation(
        dim=(2, 3),
        code="GRC",
        attrs=((AT.GRILLE, "OUI"), (AT.TYPMOD, "1D"), (AT.COQUE, "OUI"), (AT.SOUS_POINT, "OUI")),
        elements=((MT.TRIA3, EL.MEGCTR3), (MT.QUAD4, EL.MEGCQU4)),
    ),
)

phen.add(
    "GRILLE_MEMBRANE",
    Modelisation(
        dim=(2, 3),
        code="GRM",
        attrs=((AT.GRILLE, "OUI"), (AT.TYPMOD, "1D"), (AT.COQUE, "OUI")),
        elements=(
            (MT.TRIA3, EL.MEGMTR3),
            (MT.TRIA6, EL.MEGMTR6),
            (MT.QUAD4, EL.MEGMQU4),
            (MT.QUAD8, EL.MEGMQU8),
        ),
    ),
)

phen.add(
    "MEMBRANE",
    Modelisation(
        dim=(2, 3),
        code="MMB",
        attrs=((AT.TYPMOD, "MEMBRANE"), (AT.COQUE, "OUI"), (AT.EFGE, "OUI"), (AT.SIGM, "NON")),
        elements=(
            (MT.TRIA3, EL.MEMBTR3),
            (MT.TRIA6, EL.MEMBTR6),
            (MT.TRIA7, EL.MEMBTR7),
            (MT.QUAD4, EL.MEMBQU4),
            (MT.QUAD8, EL.MEMBQU8),
            (MT.QUAD9, EL.MEMBQU9),
        ),
    ),
)

phen.add(
    "PLAN_INTERFACE",
    Modelisation(
        dim=(2, 2),
        code="PEI",
        attrs=((AT.TYPMOD, "PLAN"), (AT.TYPMOD2, "INTERFAC"), (AT.INTERFACE, "OUI")),
        elements=((MT.QUAD8, EL.EIPLQU8),),
    ),
)

phen.add(
    "PLAN_INTERFACE_S",
    Modelisation(
        dim=(2, 2),
        code="PIS",
        attrs=((AT.TYPMOD, "PLAN"), (AT.TYPMOD2, "INTERFAC"), (AT.INTERFACE, "OUI")),
        elements=((MT.QUAD8, EL.EIPLQS8),),
    ),
)

phen.add(
    "PLAN_JHMS",
    Modelisation(
        dim=(2, 2),
        code="JH1",
        attrs=(
            (AT.TYPMOD, "PLAN"),
            (AT.TYPMOD2, "JHMS"),
            (AT.MECA, "OUI"),
            (AT.THER, "NON"),
            (AT.HYDR1, "1"),
            (AT.HYDR2, "0"),
            (AT.D_PLAN, "OUI"),
            (AT.INTTHM, "RED"),
        ),
        elements=((MT.QUAD8, EL.HM_J_DPQ8S), (MT.SEG3, EL.HM_J_DPSE3)),
    ),
)

phen.add(
    "PLAN_JOINT",
    Modelisation(
        dim=(2, 2),
        code="PFI",
        attrs=((AT.TYPMOD, "PLAN"), (AT.TYPMOD2, "ELEMJOIN"), (AT.INTERFACE, "OUI")),
        elements=((MT.QUAD4, EL.MFPLQU4), (MT.QUAD8, EL.MFPLQU8)),
    ),
)

phen.add(
    "PLAN_JOINT_HYME",
    Modelisation(
        dim=(2, 2),
        code="PFH",
        attrs=((AT.TYPMOD, "PLAN"), (AT.TYPMOD2, "EJ_HYME"), (AT.INTERFACE, "OUI")),
        elements=((MT.QUAD8, EL.EJHYME_PLQU8),),
    ),
)

phen.add(
    "POU_D_E",
    Modelisation(
        dim=(1, 3),
        code="PDE",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.POUX, "OUI"),
            (AT.EULER, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
            (AT.STRX, "OUI"),
        ),
        elements=((MT.SEG2, EL.MECA_POU_D_E),),
    ),
)

phen.add(
    "POU_D_EM",
    Modelisation(
        dim=(1, 3),
        code="PFM",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.POUX, "OUI"),
            (AT.TYPMOD, "1D"),
            (AT.TYPMOD2, "PMF"),
            (AT.EULER, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.STRX, "OUI"),
            (AT.SOUS_POINT, "OUI"),
        ),
        elements=((MT.SEG2, EL.MECA_POU_D_EM),),
    ),
)

phen.add(
    "POU_D_T",
    Modelisation(
        dim=(1, 3),
        code="PDT",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.POUX, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
            (AT.STRX, "OUI"),
        ),
        elements=((MT.SEG2, EL.MECA_POU_D_T),),
    ),
)

phen.add(
    "POU_D_TG",
    Modelisation(
        dim=(1, 3),
        code="PDG",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.POUX, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
            (AT.STRX, "OUI"),
        ),
        elements=((MT.SEG2, EL.MECA_POU_D_TG),),
    ),
)

phen.add(
    "POU_D_TGM",
    Modelisation(
        dim=(1, 3),
        code="PGM",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.POUX, "OUI"),
            (AT.TYPMOD, "1D"),
            (AT.TYPMOD2, "PMF"),
            (AT.EFGE, "OUI"),
            (AT.SOUS_POINT, "OUI"),
            (AT.STRX, "OUI"),
        ),
        elements=((MT.SEG2, EL.MECA_POU_D_TGM),),
    ),
)

phen.add(
    "POU_D_SQUE",
    Modelisation(
        dim=(1, 3),
        code="PGS",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.TYPMOD, "1D"),
            (AT.TYPMOD2, "PMF"),
            (AT.EFGE, "OUI"),
            (AT.SOUS_POINT, "OUI"),
            (AT.STRX, "OUI"),
        ),
        elements=((MT.SEG2, EL.MECA_POU_D_SQUE),),
    ),
)

phen.add(
    "POU_D_T_GD",
    Modelisation(
        dim=(1, 3),
        code="PGD",
        attrs=((AT.POUTRE, "OUI"), (AT.EFGE, "OUI"), (AT.SIGM, "NON"), (AT.STRX, "NON")),
        elements=((MT.SEG2, EL.MECA_POU_D_T_GD),),
    ),
)

phen.add(
    "POU_FLUI_STRU",
    Modelisation(
        dim=(1, 3),
        code="FS1",
        attrs=((AT.POUTRE, "OUI"), (AT.FLUIDE, "OUI"), (AT.FSI, "OUI"), (AT.STRX, "NON")),
        elements=((MT.SEG2, EL.MEFS_POU_D_T),),
    ),
)

phen.add(
    "Q4G",
    Modelisation(
        dim=(2, 3),
        code="Q4G",
        attrs=(
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.COQUE, "OUI"),
            (AT.PLAQUE, "OUI"),
            (AT.EFGE, "OUI"),
        ),
        elements=((MT.QUAD4, EL.MEQ4QU4), (MT.TRIA3, EL.MET3TR3), (MT.SEG2, EL.MEBOQ4G)),
    ),
)

phen.add(
    "Q4GG",
    Modelisation(
        dim=(2, 3),
        code="Q4S",
        attrs=(
            (AT.NBSIGM, "6"),
            (AT.COQUE, "OUI"),
            (AT.PLAQUE, "OUI"),
            (AT.EFGE, "OUI"),
            (AT.SIGM, "NON"),
        ),
        elements=((MT.QUAD4, EL.MEQ4GG4), (MT.TRIA3, EL.MET3GG3), (MT.SEG2, EL.MEBOQ4G)),
    ),
)

phen.add(
    "COQUE_SOLIDE",
    Modelisation(
        dim=(3, 3),
        code="SSH",
        attrs=((AT.NBSIGM, "6"), (AT.TYPMOD, "3D")),
        elements=(
            (MT.HEXA9, EL.MESSHELL_SB9),
            (MT.QUAD4, EL.MECA_FACE4),
            (MT.SEG2, EL.MECA_ARETE2),
        ),
    ),
)

phen.add(
    "TUYAU_3M",
    Modelisation(
        dim=(1, 3),
        code="TU3",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.TUYAU, "OUI"),
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.EFGE, "OUI"),
            (AT.SOUS_POINT, "OUI"),
            (AT.STRX, "NON"),
        ),
        elements=((MT.SEG3, EL.MET3SEG3), (MT.SEG4, EL.MET3SEG4)),
    ),
)

phen.add(
    "TUYAU_6M",
    Modelisation(
        dim=(1, 3),
        code="TU6",
        attrs=(
            (AT.POUTRE, "OUI"),
            (AT.TUYAU, "OUI"),
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "C_PLAN"),
            (AT.EFGE, "OUI"),
            (AT.SOUS_POINT, "OUI"),
            (AT.STRX, "NON"),
        ),
        elements=((MT.SEG3, EL.MET6SEG3),),
    ),
)

phen.add(
    "3D_HHO#4",
    Modelisation(
        dim=(3, 3),
        code="3D_",
        attrs=(
            (AT.FORMULATION, "HHO_QUAR"),
            (AT.TYPMOD2, "HHO"),
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "3D"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.MECA3DH27_HHO444),
            (MT.TETRA15, EL.MECA3DT15_HHO444),
            (MT.PENTA21, EL.MECA3DP21_HHO444),
            (MT.PYRAM19, EL.MECA3DP19_HHO444),
            (MT.QUAD9, EL.MECA3DQU9_HHO4_F),
            (MT.TRIA7, EL.MECA3DTR7_HHO4_F),
        ),
    ),
)

phen.add(
    "3D_HHO#3",
    Modelisation(
        dim=(3, 3),
        code="3D_",
        attrs=(
            (AT.FORMULATION, "HHO_CUBI"),
            (AT.TYPMOD2, "HHO"),
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "3D"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.MECA3DH27_HHO333),
            (MT.TETRA15, EL.MECA3DT15_HHO333),
            (MT.PENTA21, EL.MECA3DP21_HHO333),
            (MT.PYRAM19, EL.MECA3DP19_HHO333),
            (MT.QUAD9, EL.MECA3DQU9_HHO3_F),
            (MT.TRIA7, EL.MECA3DTR7_HHO3_F),
        ),
    ),
)

phen.add(
    "3D_HHO#2",
    Modelisation(
        dim=(3, 3),
        code="3D_",
        attrs=(
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO"),
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "3D"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.MECA3DH27_HHO222),
            (MT.TETRA15, EL.MECA3DT15_HHO222),
            (MT.PENTA21, EL.MECA3DP21_HHO222),
            (MT.PYRAM19, EL.MECA3DP19_HHO222),
            (MT.QUAD9, EL.MECA3DQU9_HHO2_F),
            (MT.TRIA7, EL.MECA3DTR7_HHO2_F),
        ),
    ),
)


phen.add(
    "3D_HHO#1",
    Modelisation(
        dim=(3, 3),
        code="3D_",
        attrs=(
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO"),
            (AT.NBSIGM, "6"),
            (AT.TYPMOD, "3D"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.MECA3DH27_HHO111),
            (MT.TETRA15, EL.MECA3DT15_HHO111),
            (MT.PENTA21, EL.MECA3DP21_HHO111),
            (MT.PYRAM19, EL.MECA3DP19_HHO111),
            (MT.QUAD9, EL.MECA3DQU9_HHO1_F),
            (MT.TRIA7, EL.MECA3DTR7_HHO1_F),
        ),
    ),
)

phen.add(
    "D_PLAN_HHO#4",
    Modelisation(
        dim=(2, 2),
        code="DPL",
        attrs=(
            (AT.FORMULATION, "HHO_QUAR"),
            (AT.TYPMOD2, "HHO"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.NBSIGM, "4"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_DPQ9_HHO444),
            (MT.TRIA7, EL.MECA_DPT7_HHO444),
            (MT.SEG3, EL.MECA_2D_HHO4_F),
        ),
    ),
)

phen.add(
    "D_PLAN_HHO#3",
    Modelisation(
        dim=(2, 2),
        code="DPL",
        attrs=(
            (AT.FORMULATION, "HHO_CUBI"),
            (AT.TYPMOD2, "HHO"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.NBSIGM, "4"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_DPQ9_HHO333),
            (MT.TRIA7, EL.MECA_DPT7_HHO333),
            (MT.SEG3, EL.MECA_2D_HHO3_F),
        ),
    ),
)

phen.add(
    "D_PLAN_HHO#2",
    Modelisation(
        dim=(2, 2),
        code="DPL",
        attrs=(
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.NBSIGM, "4"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_DPQ9_HHO222),
            (MT.TRIA7, EL.MECA_DPT7_HHO222),
            (MT.SEG3, EL.MECA_2D_HHO2_F),
        ),
    ),
)

phen.add(
    "D_PLAN_HHO#1",
    Modelisation(
        dim=(2, 2),
        code="DPL",
        attrs=(
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.NBSIGM, "4"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_DPQ9_HHO111),
            (MT.TRIA7, EL.MECA_DPT7_HHO111),
            (MT.SEG3, EL.MECA_2D_HHO1_F),
        ),
    ),
)

phen.add(
    "D_PLAN_MIX_STA#1",
    Modelisation(
        dim=(2, 2),
        code="2MV",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.FORMULATION, "STA"),
        ),
        elements=((MT.QUAD4, EL.MEMS_QU4), (MT.TRIA3, EL.MEMS_TR3), (MT.SEG2, EL.MEPLSE2)),
    ),
)

phen.add(
    "D_PLAN_MIX_STA#2",
    Modelisation(
        dim=(2, 2),
        code="2SI",
        attrs=(
            (AT.NBSIGM, "4"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.FORMULATION, "STA_INCO"),
        ),
        elements=((MT.QUAD4, EL.MEMI_QU4), (MT.SEG2, EL.MEPLSE2)),
    ),
)

phen.add(
    "D_PLAN_GRAD_HH#1",
    Modelisation(
        dim=(2, 2),
        code="DGL",
        attrs=(
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO_GRAD"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.NBSIGM, "4"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_DGVQ_HHO111),
            (MT.TRIA7, EL.MECA_DGVT_HHO111),
            (MT.SEG3, EL.MECA_2DGV_HHO1_F),
        ),
    ),
)

phen.add(
    "D_PLAN_GRAD_HH#2",
    Modelisation(
        dim=(2, 2),
        code="DGQ",
        attrs=(
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO_GRAD"),
            (AT.D_PLAN, "OUI"),
            (AT.TYPMOD, "D_PLAN"),
            (AT.NBSIGM, "4"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_DGVQ_HHO222),
            (MT.TRIA7, EL.MECA_DGVT_HHO222),
            (MT.SEG3, EL.MECA_2DGV_HHO2_F),
        ),
    ),
)

phen.add(
    "AXIS_HHO#1",
    Modelisation(
        dim=(2, 2),
        code="HA1",
        attrs=(
            (AT.FORMULATION, "HHO_MLINE"),
            (AT.AXIS, "OUI"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
            (AT.NBSIGM, "4"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_AXQ9_HHO121),
            (MT.TRIA7, EL.MECA_AXT7_HHO121),
            (MT.SEG3, EL.MECA_2DAX_HHO1_F),
        ),
    ),
)


phen.add(
    "AXIS_HHO#2",
    Modelisation(
        dim=(2, 2),
        code="HA2",
        attrs=(
            (AT.FORMULATION, "HHO_MQUAD"),
            (AT.AXIS, "OUI"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
            (AT.NBSIGM, "4"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_AXQ9_HHO232),
            (MT.TRIA7, EL.MECA_AXT7_HHO232),
            (MT.SEG3, EL.MECA_2DAX_HHO2_F),
        ),
    ),
)

phen.add(
    "AXIS_HHO#3",
    Modelisation(
        dim=(2, 2),
        code="HA3",
        attrs=(
            (AT.FORMULATION, "HHO_MCUBI"),
            (AT.AXIS, "OUI"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
            (AT.NBSIGM, "4"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_AXQ9_HHO343),
            (MT.TRIA7, EL.MECA_AXT7_HHO343),
            (MT.SEG3, EL.MECA_2DAX_HHO3_F),
        ),
    ),
)

phen.add(
    "AXIS_HHO#4",
    Modelisation(
        dim=(2, 2),
        code="HA3",
        attrs=(
            (AT.FORMULATION, "HHO_MQUAR"),
            (AT.AXIS, "OUI"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
            (AT.NBSIGM, "4"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_AXQ9_HHO454),
            (MT.TRIA7, EL.MECA_AXT7_HHO454),
            (MT.SEG3, EL.MECA_2DAX_HHO4_F),
        ),
    ),
)

phen.add(
    "3D_GRAD_HHO#1",
    Modelisation(
        dim=(3, 3),
        code="3GL",
        attrs=(
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO_GRAD"),
            (AT.TYPMOD, "3D"),
            (AT.NBSIGM, "6"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.MECA3DGVH_HHO111),
            (MT.TETRA15, EL.MECA3DGVT_HHO111),
            (MT.PENTA21, EL.MECA3GVPE_HHO111),
            (MT.PYRAM19, EL.MECA3GVPY_HHO111),
            (MT.QUAD9, EL.MECA3DGVQ_HHO1_F),
            (MT.TRIA7, EL.MECA3DGVT_HHO1_F),
        ),
    ),
)

phen.add(
    "3D_GRAD_HHO#2",
    Modelisation(
        dim=(3, 3),
        code="3GQ",
        attrs=(
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO_GRAD"),
            (AT.TYPMOD, "3D"),
            (AT.NBSIGM, "6"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.MECA3DGVH_HHO222),
            (MT.TETRA15, EL.MECA3DGVT_HHO222),
            (MT.PENTA21, EL.MECA3GVPE_HHO222),
            (MT.PYRAM19, EL.MECA3GVPY_HHO222),
            (MT.QUAD9, EL.MECA3DGVQ_HHO2_F),
            (MT.TRIA7, EL.MECA3DGVT_HHO2_F),
        ),
    ),
)

phen.add(
    "AXIS_GRAD_HHO#1",
    Modelisation(
        dim=(2, 2),
        code="HA1",
        attrs=(
            (AT.FORMULATION, "HHO_LINE"),
            (AT.AXIS, "OUI"),
            (AT.TYPMOD2, "HHO_GRAD"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
            (AT.NBSIGM, "4"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_AGVQ_HHO111),
            (MT.TRIA7, EL.MECA_AGVT_HHO111),
            (MT.SEG3, EL.MECA_AXGV_HHO1_F),
        ),
    ),
)


phen.add(
    "AXIS_GRAD_HHO#2",
    Modelisation(
        dim=(2, 2),
        code="HA2",
        attrs=(
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.AXIS, "OUI"),
            (AT.TYPMOD2, "HHO_GRAD"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
            (AT.NBSIGM, "4"),
        ),
        elements=(
            (MT.QUAD9, EL.MECA_AGVQ_HHO222),
            (MT.TRIA7, EL.MECA_AGVT_HHO222),
            (MT.SEG3, EL.MECA_AXGV_HHO2_F),
        ),
    ),
)

############################################################
# Pour le phenomene : THERMIQUE :
############################################################
THERMIQUE = Phenomenon(code="TH")
phen = THERMIQUE


phen.add(
    "3D",
    Modelisation(
        dim=(3, 3),
        code="3D_",
        elements=(
            (MT.HEXA8, EL.THER_HEXA8),
            (MT.PENTA6, EL.THER_PENTA6),
            (MT.TETRA4, EL.THER_TETRA4),
            (MT.PYRAM5, EL.THER_PYRAM5),
            (MT.QUAD4, EL.THER_FACE4),
            (MT.TRIA3, EL.THER_FACE3),
            (MT.HEXA27, EL.THER_HEXA27),
            (MT.HEXA20, EL.THER_HEXA20),
            (MT.PENTA15, EL.THER_PENTA15),
            (MT.TETRA10, EL.THER_TETRA10),
            (MT.TETRA20, EL.THER_TETRA20),
            (MT.PYRAM13, EL.THER_PYRAM13),
            (MT.QUAD8, EL.THER_FACE8),
            (MT.QUAD9, EL.THER_FACE9),
            (MT.TRIA6, EL.THER_FACE6),
            (MT.TRIA10, EL.THER_FACE10),
        ),
    ),
)

phen.add(
    "3D_SECH",
    Modelisation(
        dim=(3, 3),
        code="3DH",
        attrs=((AT.SECH, "OUI"),),
        elements=(
            (MT.HEXA8, EL.SECH_HEXA8),
            (MT.PENTA6, EL.SECH_PENTA6),
            (MT.TETRA4, EL.SECH_TETRA4),
            (MT.PYRAM5, EL.SECH_PYRAM5),
            (MT.QUAD4, EL.SECH_FACE4),
            (MT.TRIA3, EL.SECH_FACE3),
            (MT.HEXA27, EL.SECH_HEXA27),
            (MT.HEXA20, EL.SECH_HEXA20),
            (MT.PENTA15, EL.SECH_PENTA15),
            (MT.TETRA10, EL.SECH_TETRA10),
            (MT.PYRAM13, EL.SECH_PYRAM13),
            (MT.QUAD9, EL.SECH_FACE9),
            (MT.QUAD8, EL.SECH_FACE8),
            (MT.TRIA6, EL.SECH_FACE6),
        ),
    ),
)

phen.add(
    "3D1XH",
    Modelisation(
        dim=(3, 3),
        code="3X1",
        attrs=((AT.LXFEM, "OUI"), (AT.XFEM, "XH")),
        elements=(
            (MT.HEXA8, EL.THER_XH_HEXA8),
            (MT.PENTA6, EL.THER_XH_PENTA6),
            (MT.PYRAM5, EL.THER_XH_PYRAM5),
            (MT.TETRA4, EL.THER_XH_TETRA4),
            (MT.QUAD4, EL.THER_XH_FACE4),
            (MT.TRIA3, EL.THER_XH_FACE3),
        ),
    ),
)

phen.add(
    "3D1XHT",
    Modelisation(
        dim=(3, 3),
        code="3X3",
        attrs=((AT.LXFEM, "OUI"), (AT.XFEM, "XHT")),
        elements=(
            (MT.HEXA8, EL.THER_XHT_HEXA8),
            (MT.PENTA6, EL.THER_XHT_PENTA6),
            (MT.PYRAM5, EL.THER_XHT_PYRAM5),
            (MT.TETRA4, EL.THER_XHT_TETRA4),
            (MT.QUAD4, EL.THER_XHT_FACE4),
            (MT.TRIA3, EL.THER_XHT_FACE3),
        ),
    ),
)

phen.add(
    "3D1XT",
    Modelisation(
        dim=(3, 3),
        code="3X2",
        attrs=((AT.LXFEM, "OUI"), (AT.XFEM, "XT")),
        elements=(
            (MT.HEXA8, EL.THER_XT_HEXA8),
            (MT.PENTA6, EL.THER_XT_PENTA6),
            (MT.PYRAM5, EL.THER_XT_PYRAM5),
            (MT.TETRA4, EL.THER_XT_TETRA4),
            (MT.QUAD4, EL.THER_XT_FACE4),
            (MT.TRIA3, EL.THER_XT_FACE3),
        ),
    ),
)

phen.add(
    "3D_DIAG",
    Modelisation(
        dim=(3, 3),
        code="3DD",
        attrs=((AT.LUMPE, "OUI"),),
        elements=(
            (MT.HEXA8, EL.THER_HEXA8_D),
            (MT.PENTA6, EL.THER_PENTA6_D),
            (MT.TETRA4, EL.THER_TETRA4_D),
            (MT.QUAD4, EL.THER_FACE4_D),
            (MT.TRIA3, EL.THER_FACE3_D),
        ),
    ),
)

phen.add(
    "3D_SECH_DIAG",
    Modelisation(
        dim=(3, 3),
        code="3DJ",
        attrs=((AT.LUMPE, "OUI"), (AT.SECH, "OUI")),
        elements=(
            (MT.HEXA8, EL.SECH_HEXA8_D),
            (MT.PENTA6, EL.SECH_PENTA6_D),
            (MT.TETRA4, EL.SECH_TETRA4_D),
            (MT.QUAD4, EL.SECH_FACE4_D),
            (MT.TRIA3, EL.SECH_FACE3_D),
        ),
    ),
)

phen.add(
    "AXIS",
    Modelisation(
        dim=(2, 2),
        code="AX_",
        attrs=((AT.AXIS, "OUI"),),
        elements=(
            (MT.TRIA3, EL.THAXTR3),
            (MT.QUAD4, EL.THAXQU4),
            (MT.TRIA6, EL.THAXTR6),
            (MT.QUAD8, EL.THAXQU8),
            (MT.QUAD9, EL.THAXQU9),
            (MT.SEG2, EL.THAXSE2),
            (MT.SEG3, EL.THAXSE3),
        ),
    ),
)

phen.add(
    "AXIS_DIAG",
    Modelisation(
        dim=(2, 2),
        code="AXD",
        attrs=((AT.AXIS, "OUI"), (AT.LUMPE, "OUI")),
        elements=((MT.TRIA3, EL.THAXTL3), (MT.QUAD4, EL.THAXQL4), (MT.SEG2, EL.THAXSL2)),
    ),
)

phen.add(
    "AXIS_SECH",
    Modelisation(
        dim=(2, 2),
        code="AXH",
        attrs=((AT.AXIS, "OUI"), (AT.SECH, "OUI")),
        elements=(
            (MT.TRIA3, EL.SEAXTR3),
            (MT.QUAD4, EL.SEAXQU4),
            (MT.TRIA6, EL.SEAXTR6),
            (MT.QUAD8, EL.SEAXQU8),
            (MT.QUAD9, EL.SEAXQU9),
            (MT.SEG2, EL.SEAXSE2),
            (MT.SEG3, EL.SEAXSE3),
        ),
    ),
)

phen.add(
    "AXIS_SECH_DIAG",
    Modelisation(
        dim=(2, 2),
        code="AXJ",
        attrs=((AT.AXIS, "OUI"), (AT.LUMPE, "OUI"), (AT.SECH, "OUI")),
        elements=((MT.TRIA3, EL.SEAXTL3), (MT.QUAD4, EL.SEAXQL4), (MT.SEG2, EL.SEAXSL2)),
    ),
)

phen.add(
    "AXIS_FOURIER",
    Modelisation(
        dim=(2, 2),
        code="AXF",
        attrs=((AT.AXIS, "OUI"),),
        elements=(
            (MT.TRIA3, EL.THFOTR3),
            (MT.QUAD4, EL.THFOQU4),
            (MT.SEG2, EL.THFOSE2),
            (MT.TRIA6, EL.THFOTR6),
            (MT.QUAD8, EL.THFOQU8),
            (MT.QUAD9, EL.THFOQU9),
            (MT.SEG3, EL.THFOSE3),
        ),
    ),
)

phen.add(
    "AXIS_XH",
    Modelisation(
        dim=(2, 2),
        code="AX1",
        attrs=((AT.AXIS, "OUI"), (AT.LXFEM, "OUI"), (AT.XFEM, "XH")),
        elements=((MT.TRIA3, EL.THAXTR3_XH), (MT.QUAD4, EL.THAXQU4_XH), (MT.SEG2, EL.THAXSE2_XH)),
    ),
)

phen.add(
    "AXIS_XHT",
    Modelisation(
        dim=(2, 2),
        code="AX3",
        attrs=((AT.AXIS, "OUI"), (AT.LXFEM, "OUI"), (AT.XFEM, "XHT")),
        elements=(
            (MT.TRIA3, EL.THAXTR3_XHT),
            (MT.QUAD4, EL.THAXQU4_XHT),
            (MT.SEG2, EL.THAXSE2_XHT),
        ),
    ),
)

phen.add(
    "AXIS_XT",
    Modelisation(
        dim=(2, 2),
        code="AX2",
        attrs=((AT.AXIS, "OUI"), (AT.LXFEM, "OUI"), (AT.XFEM, "XT")),
        elements=((MT.TRIA3, EL.THAXTR3_XT), (MT.QUAD4, EL.THAXQU4_XT), (MT.SEG2, EL.THAXSE2_XT)),
    ),
)

# -------------------------------------------------------------------
# Modelisations sous-terraines pour :
#  * conditions d'echange thermique entre deux parois

phen.add(
    "CL_ECHANGE1",
    Modelisation(
        dim=(2, 3),
        code="CL1",
        elements=(
            (MT.TRIA33, EL.THER_FACE33),
            (MT.QUAD44, EL.THER_FACE44),
            (MT.TRIA66, EL.THER_FACE66),
            (MT.QUAD88, EL.THER_FACE88),
            (MT.QUAD99, EL.THER_FACE99),
        ),
    ),
)

phen.add(
    "CL_ECHANGE2",
    Modelisation(
        dim=(1, 2),
        code="CL2",
        elements=(
            (MT.SEG22, EL.THAXSE22),
            (MT.SEG33, EL.THAXSE33),
            (MT.SEG22, EL.THPLSE22),
            (MT.SEG33, EL.THPLSE33),
        ),
    ),
)
# -------------------------------------------------------------------

phen.add(
    "COQUE",
    Modelisation(
        dim=(2, 3),
        code="CQ_",
        attrs=((AT.COQUE, "OUI"),),
        elements=(
            (MT.TRIA3, EL.THCOTR3),
            (MT.TRIA6, EL.THCOTR6),
            (MT.TRIA7, EL.THCOTR7),
            (MT.QUAD4, EL.THCOQU4),
            (MT.QUAD8, EL.THCOQU8),
            (MT.QUAD9, EL.THCOQU9),
            (MT.SEG3, EL.THCOSE3),
            (MT.SEG2, EL.THCOSE2),
        ),
    ),
)

phen.add(
    "COQUE_AXIS",
    Modelisation(
        dim=(1, 2), code="CTA", attrs=((AT.COQUE, "OUI"),), elements=((MT.SEG3, EL.THCASE3),)
    ),
)

phen.add(
    "COQUE_PLAN",
    Modelisation(
        dim=(1, 2), code="CQP", attrs=((AT.COQUE, "OUI"),), elements=((MT.SEG3, EL.THCPSE3),)
    ),
)

phen.add(
    "PLAN",
    Modelisation(
        dim=(2, 2),
        code="PL_",
        elements=(
            (MT.TRIA3, EL.THPLTR3),
            (MT.QUAD4, EL.THPLQU4),
            (MT.TRIA6, EL.THPLTR6),
            (MT.TRIA10, EL.THPLTR10),
            (MT.QUAD8, EL.THPLQU8),
            (MT.QUAD9, EL.THPLQU9),
            (MT.QUAD12, EL.THPLQU12),
            (MT.SEG2, EL.THPLSE2),
            (MT.SEG3, EL.THPLSE3),
            (MT.SEG4, EL.THPLSE4),
        ),
    ),
)

phen.add(
    "PLAN_DIAG",
    Modelisation(
        dim=(2, 2),
        code="PLD",
        attrs=((AT.LUMPE, "OUI"),),
        elements=((MT.TRIA3, EL.THPLTL3), (MT.QUAD4, EL.THPLQL4), (MT.SEG2, EL.THPLSL2)),
    ),
)

phen.add(
    "PLAN_XH",
    Modelisation(
        dim=(2, 2),
        code="PX1",
        attrs=((AT.LXFEM, "OUI"), (AT.XFEM, "XH")),
        elements=((MT.TRIA3, EL.THPLTR3_XH), (MT.QUAD4, EL.THPLQU4_XH), (MT.SEG2, EL.THPLSE2_XH)),
    ),
)

phen.add(
    "PLAN_XHT",
    Modelisation(
        dim=(2, 2),
        code="PX3",
        attrs=((AT.LXFEM, "OUI"), (AT.XFEM, "XHT")),
        elements=(
            (MT.TRIA3, EL.THPLTR3_XHT),
            (MT.QUAD4, EL.THPLQU4_XHT),
            (MT.SEG2, EL.THPLSE2_XHT),
        ),
    ),
)

phen.add(
    "PLAN_XT",
    Modelisation(
        dim=(2, 2),
        code="PX2",
        attrs=((AT.LXFEM, "OUI"), (AT.XFEM, "XT")),
        elements=((MT.TRIA3, EL.THPLTR3_XT), (MT.QUAD4, EL.THPLQU4_XT), (MT.SEG2, EL.THPLSE2_XT)),
    ),
)

phen.add(
    "3D_HHO#4",
    Modelisation(
        dim=(3, 3),
        code="HT1",
        attrs=(
            (AT.FORMULATION, "HHO_QUAR"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "3D"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.THER3DH27_HHO444),
            (MT.TETRA15, EL.THER3DT15_HHO444),
            (MT.PYRAM19, EL.THER3DP19_HHO444),
            (MT.PENTA21, EL.THER3DP21_HHO444),
            (MT.QUAD9, EL.THER3DQU9_HHO4_F),
            (MT.TRIA7, EL.THER3DTR7_HHO4_F),
        ),
    ),
)

phen.add(
    "3D_HHO#3",
    Modelisation(
        dim=(3, 3),
        code="HT1",
        attrs=(
            (AT.FORMULATION, "HHO_CUBI"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "3D"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.THER3DH27_HHO333),
            (MT.TETRA15, EL.THER3DT15_HHO333),
            (MT.PYRAM19, EL.THER3DP19_HHO333),
            (MT.PENTA21, EL.THER3DP21_HHO333),
            (MT.QUAD9, EL.THER3DQU9_HHO3_F),
            (MT.TRIA7, EL.THER3DTR7_HHO3_F),
        ),
    ),
)

phen.add(
    "3D_HHO#2",
    Modelisation(
        dim=(3, 3),
        code="HT1",
        attrs=(
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "3D"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.THER3DH27_HHO222),
            (MT.TETRA15, EL.THER3DT15_HHO222),
            (MT.PYRAM19, EL.THER3DP19_HHO222),
            (MT.PENTA21, EL.THER3DP21_HHO222),
            (MT.QUAD9, EL.THER3DQU9_HHO2_F),
            (MT.TRIA7, EL.THER3DTR7_HHO2_F),
        ),
    ),
)


phen.add(
    "3D_HHO#1",
    Modelisation(
        dim=(3, 3),
        code="HT2",
        attrs=(
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "3D"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.THER3DH27_HHO111),
            (MT.TETRA15, EL.THER3DT15_HHO111),
            (MT.PYRAM19, EL.THER3DP19_HHO111),
            (MT.PENTA21, EL.THER3DP21_HHO111),
            (MT.QUAD9, EL.THER3DQU9_HHO1_F),
            (MT.TRIA7, EL.THER3DTR7_HHO1_F),
        ),
    ),
)

phen.add(
    "3D_HHO#0",
    Modelisation(
        dim=(3, 3),
        code="HT0",
        attrs=(
            (AT.FORMULATION, "HHO_CSTE"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "3D"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.HEXA27, EL.THER3DH27_HHO000),
            (MT.TETRA15, EL.THER3DT15_HHO000),
            (MT.PYRAM19, EL.THER3DP19_HHO000),
            (MT.PENTA21, EL.THER3DP21_HHO000),
            (MT.QUAD9, EL.THER3DQU9_HHO0_F),
            (MT.TRIA7, EL.THER3DTR7_HHO0_F),
        ),
    ),
)

phen.add(
    "PLAN_HHO#4",
    Modelisation(
        dim=(2, 2),
        code="PT3",
        attrs=(
            (AT.FORMULATION, "HHO_QUAR"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "PLAN"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.THER2DQ9_HHO444),
            (MT.TRIA7, EL.THER2DT7_HHO444),
            (MT.SEG3, EL.THER_2D_HHO4_F),
        ),
    ),
)

phen.add(
    "PLAN_HHO#3",
    Modelisation(
        dim=(2, 2),
        code="PT3",
        attrs=(
            (AT.FORMULATION, "HHO_CUBI"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "PLAN"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.THER2DQ9_HHO333),
            (MT.TRIA7, EL.THER2DT7_HHO333),
            (MT.SEG3, EL.THER_2D_HHO3_F),
        ),
    ),
)

phen.add(
    "PLAN_HHO#2",
    Modelisation(
        dim=(2, 2),
        code="PT3",
        attrs=(
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "PLAN"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.THER2DQ9_HHO222),
            (MT.TRIA7, EL.THER2DT7_HHO222),
            (MT.SEG3, EL.THER_2D_HHO2_F),
        ),
    ),
)

phen.add(
    "PLAN_HHO#1",
    Modelisation(
        dim=(2, 2),
        code="PT4",
        attrs=(
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "PLAN"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.THER2DQ9_HHO111),
            (MT.TRIA7, EL.THER2DT7_HHO111),
            (MT.SEG3, EL.THER_2D_HHO1_F),
        ),
    ),
)

phen.add(
    "PLAN_HHO#0",
    Modelisation(
        dim=(2, 2),
        code="PT0",
        attrs=(
            (AT.FORMULATION, "HHO_CSTE"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "PLAN"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.THER2DQ9_HHO000),
            (MT.TRIA7, EL.THER2DT7_HHO000),
            (MT.SEG3, EL.THER_2D_HHO0_F),
        ),
    ),
)

phen.add(
    "AXIS_HHO#0",
    Modelisation(
        dim=(2, 2),
        code="PA0",
        attrs=(
            (AT.FORMULATION, "HHO_CSTE"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.THERAXQ9_HHO000),
            (MT.TRIA7, EL.THERAXT7_HHO000),
            (MT.SEG3, EL.THER_AX_HHO0_F),
        ),
    ),
)

phen.add(
    "AXIS_HHO#1",
    Modelisation(
        dim=(2, 2),
        code="PA4",
        attrs=(
            (AT.FORMULATION, "HHO_LINE"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.THERAXQ9_HHO111),
            (MT.TRIA7, EL.THERAXT7_HHO111),
            (MT.SEG3, EL.THER_AX_HHO1_F),
        ),
    ),
)


phen.add(
    "AXIS_HHO#2",
    Modelisation(
        dim=(2, 2),
        code="PA3",
        attrs=(
            (AT.FORMULATION, "HHO_QUAD"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.THERAXQ9_HHO222),
            (MT.TRIA7, EL.THERAXT7_HHO222),
            (MT.SEG3, EL.THER_AX_HHO2_F),
        ),
    ),
)

phen.add(
    "AXIS_HHO#3",
    Modelisation(
        dim=(2, 2),
        code="PA3",
        attrs=(
            (AT.FORMULATION, "HHO_CUBI"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.THERAXQ9_HHO333),
            (MT.TRIA7, EL.THERAXT7_HHO333),
            (MT.SEG3, EL.THER_AX_HHO3_F),
        ),
    ),
)

phen.add(
    "AXIS_HHO#4",
    Modelisation(
        dim=(2, 2),
        code="PA3",
        attrs=(
            (AT.FORMULATION, "HHO_QUAR"),
            (AT.TYPMOD2, "HHO"),
            (AT.TYPMOD, "AXIS"),
            (AT.HHO, "OUI"),
        ),
        elements=(
            (MT.QUAD9, EL.THERAXQ9_HHO444),
            (MT.TRIA7, EL.THERAXT7_HHO444),
            (MT.SEG3, EL.THER_AX_HHO4_F),
        ),
    ),
)

############################################################
# Pour le phenomene : ACOUSTIQUE :
############################################################
ACOUSTIQUE = Phenomenon(code="AC")
phen = ACOUSTIQUE

phen.add(
    "3D",
    Modelisation(
        dim=(3, 3),
        code="3D_",
        elements=(
            (MT.HEXA8, EL.ACOU_HEXA8),
            (MT.PENTA6, EL.ACOU_PENTA6),
            (MT.TETRA4, EL.ACOU_TETRA4),
            (MT.QUAD4, EL.ACOUNA_FACE4),
            (MT.TRIA3, EL.ACOUNA_FACE3),
            (MT.HEXA27, EL.ACOU_HEXA27),
            (MT.HEXA20, EL.ACOU_HEXA20),
            (MT.PENTA15, EL.ACOU_PENTA15),
            (MT.TETRA10, EL.ACOU_TETRA10),
            (MT.QUAD9, EL.ACOUNA_FACE9),
            (MT.QUAD8, EL.ACOUNA_FACE8),
            (MT.TRIA6, EL.ACOUNA_FACE6),
        ),
    ),
)


phen.add(
    "3D_ABSO",
    Modelisation(
        dim=(3, 3),
        code="3AA",
        attrs=((AT.TYPMOD, "3D"), (AT.ABSO, "OUI")),
        elements=(
            (MT.QUAD4, EL.ACOU_FACE4),
            (MT.TRIA3, EL.ACOU_FACE3),
            (MT.QUAD9, EL.ACOU_FACE9),
            (MT.QUAD8, EL.ACOU_FACE8),
            (MT.TRIA6, EL.ACOU_FACE6),
        ),
    ),
)

phen.add(
    "PLAN",
    Modelisation(
        dim=(2, 2),
        code="PLA",
        elements=(
            (MT.TRIA3, EL.ACPLTR3),
            (MT.TRIA6, EL.ACPLTR6),
            (MT.QUAD4, EL.ACPLQU4),
            (MT.QUAD8, EL.ACPLQU8),
            (MT.QUAD9, EL.ACPLQU9),
            (MT.SEG2, EL.ACPLNASE2),
            (MT.SEG3, EL.ACPLNASE3),
        ),
    ),
)

phen.add(
    "PLAN_ABSO",
    Modelisation(
        dim=(2, 2),
        code="PAA",
        attrs=((AT.TYPMOD, "PLAN"), (AT.ABSO, "OUI")),
        elements=((MT.SEG2, EL.ACPLSE2), (MT.SEG3, EL.ACPLSE3)),
    ),
)


############################################################
# Pour le phenomene : PRESENTATION :
############################################################
PRESENTATION = Phenomenon(code="PR")
phen = PRESENTATION

# Les deux modelisations suivantes servent pour des calculs purement
# geometriques lies a XFEM

phen.add(
    "2D_GEOM",
    Modelisation(
        dim=(2, 2),
        code="G2D",
        elements=(
            (MT.QUAD8, EL.PR_G_QUAD8),
            (MT.QUAD4, EL.PR_G_QUAD4),
            (MT.TRIA6, EL.PR_G_TRIA6),
            (MT.TRIA3, EL.PR_G_TRIA3),
        ),
    ),
)

phen.add(
    "3D_GEOM",
    Modelisation(
        dim=(3, 3),
        code="G3D",
        elements=(
            (MT.HEXA20, EL.PR_G_HEXA20),
            (MT.HEXA8, EL.PR_G_HEXA8),
            (MT.PENTA15, EL.PR_G_PENTA15),
            (MT.PENTA6, EL.PR_G_PENTA6),
            (MT.TETRA10, EL.PR_G_TETRA10),
            (MT.TETRA4, EL.PR_G_TETRA4),
            (MT.PYRAM13, EL.PR_G_PYRAM13),
            (MT.PYRAM5, EL.PR_G_PYRAM5),
        ),
    ),
)

# La modelisation 'TOUT' sert dans IMPR_RESU
phen.add(
    "TOUT",
    Modelisation(
        dim=(3, 3),
        code="TOU",
        elements=(
            (MT.HEXA27, EL.PR_HEXA27),
            (MT.HEXA20, EL.PR_HEXA20),
            (MT.HEXA8, EL.PR_HEXA8),
            (MT.PENTA18, EL.PR_PENTA18),
            (MT.PENTA15, EL.PR_PENTA15),
            (MT.PENTA6, EL.PR_PENTA6),
            (MT.TETRA10, EL.PR_TETRA10),
            (MT.TETRA4, EL.PR_TETRA4),
            (MT.PYRAM13, EL.PR_PYRAM13),
            (MT.PYRAM5, EL.PR_PYRAM5),
            (MT.QUAD9, EL.PR_QUAD9),
            (MT.QUAD8, EL.PR_QUAD8),
            (MT.QUAD4, EL.PR_QUAD4),
            (MT.TRIA7, EL.PR_TRIA7),
            (MT.TRIA6, EL.PR_TRIA6),
            (MT.TRIA3, EL.PR_TRIA3),
            (MT.SEG4, EL.PR_SEG4),
            (MT.SEG3, EL.PR_SEG3),
            (MT.SEG2, EL.PR_SEG2),
            (MT.POI1, EL.PR_POI1),
        ),
    ),
)


# store Phenomenon objects
PHEN = objects_from_context(globals(), Phenomenon, ignore_names=["phen"])
