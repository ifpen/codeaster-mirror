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

from code_aster.CA import MPI
from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities import SharedTmpdir

import os.path as osp

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

mesh = CA.ParallelMesh()
mesh.readMedFile("fort.20")
coords = mesh.getCoordinates().getValues()
corresp = {}
for i in range(int(len(coords) / 3)):
    if corresp.get(coords[i * 3]) is None:
        corresp[coords[i * 3]] = {}
    corresp[coords[i * 3]][coords[i * 3 + 1]] = int(i)

model = AFFE_MODELE(
    MAILLAGE=mesh,
    AFFE=(
        _F(GROUP_MA=("acier"), PHENOMENE="MECANIQUE", MODELISATION="D_PLAN"),
        _F(GROUP_MA=("acier"), PHENOMENE="MECANIQUE", MODELISATION="D_PLAN_SI"),
        _F(GROUP_MA=("beton"), PHENOMENE="MECANIQUE", MODELISATION="D_PLAN_GRAD_VARI"),
    ),
)

corniere = DEFI_MATERIAU(
    ELAS=_F(E=211000.0, NU=0.3),
    ECRO_NL=_F(R0=383, R1=564, GAMMA_1=4.41),
    NON_LOCAL=_F(C_GRAD_VARI=0, PENA_LAGR=7660),
    NORTON=_F(K=50, N=5),
    VISC_ELAS=_F(TAU=1, K=41.5),
)

liner = DEFI_MATERIAU(
    ELAS=_F(E=215000, NU=0.3),
    ECRO_NL=_F(R0=298, R1=561, GAMMA_1=2.61),
    NON_LOCAL=_F(C_GRAD_VARI=0, PENA_LAGR=5960),
    NORTON=_F(K=50, N=5),
    VISC_ELAS=_F(TAU=1, K=41.5),
)

fieldmat = AFFE_MATERIAU(
    MAILLAGE=mesh,
    AFFE=(
        _F(GROUP_MA="beton", MATER=corniere),
        _F(GROUP_MA="corniere", MATER=corniere),
        _F(GROUP_MA="liner", MATER=liner),
    ),
)

liaison = AFFE_CHAR_MECA(
    MODELE=model,
    DDL_IMPO=(
        _F(GROUP_MA="edge_poutre", DY=0),
        _F(GROUP_MA="beton_o", DX=0),
        _F(GROUP_MA="beton_s", DY=0),
    ),
)

liaison2 = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="beton_e", DX=0.01),))

evol_endo = STAT_NON_LINE(
    CHAM_MATER=fieldmat,
    MODELE=model,
    EXCIT=(_F(CHARGE=liaison), _F(CHARGE=liaison2)),
    COMPORTEMENT=(
        _F(GROUP_MA="acier", RELATION="VISC_ISOT_NL", DEFORMATION="PETIT", REGU_VISC="OUI"),
        _F(GROUP_MA="beton", RELATION="VISC_ISOT_NL", DEFORMATION="PETIT", REGU_VISC="OUI"),
    ),
    INCREMENT=_F(LIST_INST=DEFI_LIST_REEL(VALE=[0, 1])),
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=1),
    CONVERGENCE=_F(
        RESI_REFE_RELA=0.01,
        SIGM_REFE=4.6,
        VARI_REFE=1,
        LAGR_REFE=0.013083333333333334,
        EFFORT_REFE=1788,
        ITER_GLOB_MAXI=10,
    ),
    ARCHIVAGE=_F(PAS_ARCH=200),
)

with SharedTmpdir("zzzz_") as tmpdir:
    medfile = osp.join(tmpdir.path, "resu.zzzz.med")
    DEFI_FICHIER(FICHIER=medfile, TYPE="BINARY", UNITE=99)
    IMPR_RESU(
        FICHIER_UNIQUE="OUI",
        FORMAT="MED",
        UNITE=99,
        RESU=(
            _F(
                RESULTAT=evol_endo,
                NOM_CHAM=("VARI_ELGA", "DEPL", "SIEF_ELGA"),
                TOUT_ORDRE="OUI",
                NOM_RESU_MED="RESU_SEQ",
            ),
        ),
    )
    DEFI_FICHIER(ACTION="LIBERER", UNITE=99)
    DEFI_FICHIER(UNITE=99, ACCES="OLD", FICHIER=medfile)

    mesh_seq = LIRE_MAILLAGE(UNITE=99)
    coords = mesh_seq.getCoordinates().getValues()
    corresp2 = {}
    for i in range(int(len(coords) / 3)):
        if corresp2.get(coords[i * 3]) is None:
            corresp2[coords[i * 3]] = {}
        corresp2[coords[i * 3]][coords[i * 3 + 1]] = int(i)

    parSeqCorresp = {}
    for x, val in corresp.items():
        val2 = corresp2[x]
        for y, id1 in val.items():
            id2 = val2[y]
            parSeqCorresp[id1] = id2

    model_seq = AFFE_MODELE(
        MAILLAGE=mesh_seq,
        AFFE=(
            _F(GROUP_MA=("acier"), PHENOMENE="MECANIQUE", MODELISATION="D_PLAN"),
            _F(GROUP_MA=("acier"), PHENOMENE="MECANIQUE", MODELISATION="D_PLAN_SI"),
            _F(GROUP_MA=("beton"), PHENOMENE="MECANIQUE", MODELISATION="D_PLAN_GRAD_VARI"),
        ),
    )
    fieldmat_seq = AFFE_MATERIAU(
        MAILLAGE=mesh_seq,
        AFFE=(
            _F(GROUP_MA="beton", MATER=corniere),
            _F(GROUP_MA="corniere", MATER=corniere),
            _F(GROUP_MA="liner", MATER=liner),
        ),
    )

    seq_result = LIRE_RESU(
        MODELE=model_seq,
        FORMAT="MED",
        UNITE=99,
        TYPE_RESU="EVOL_NOLI",
        COMPORTEMENT=(
            _F(GROUP_MA="acier", RELATION="VISC_ISOT_NL", DEFORMATION="PETIT", REGU_VISC="OUI"),
            _F(GROUP_MA="beton", RELATION="VISC_ISOT_NL", DEFORMATION="PETIT", REGU_VISC="OUI"),
        ),
        CHAM_MATER=fieldmat_seq,
        FORMAT_MED=(
            _F(NOM_RESU="RESU_SEQ", NOM_CHAM="DEPL"),
            _F(NOM_RESU="RESU_SEQ", NOM_CHAM="SIEF_ELGA"),
            _F(NOM_RESU="RESU_SEQ", NOM_CHAM="VARI_ELGA"),
        ),
        TOUT_ORDRE="OUI",
    )
    DEFI_FICHIER(ACTION="LIBERER", UNITE=99)

    evol_endo.build()
    deplPar = evol_endo.getField("DEPL", 1).toSimpleFieldOnNodes()
    seq_result.build()
    deplSeq = seq_result.getField("DEPL", 1).toSimpleFieldOnNodes()
    nbCmp = deplPar.getNumberOfComponents()
    valPar = deplPar.getValues()[0]
    valSeq = deplSeq.getValues()[0]
    for idPar, idSeq in parSeqCorresp.items():
        test.assertAlmostEqual(valPar[idPar][0], valSeq[idSeq][0])

    DEFI_FICHIER(FICHIER=medfile, TYPE="BINARY", UNITE=99)
    IMPR_RESU(
        FICHIER_UNIQUE="OUI",
        FORMAT="MED",
        UNITE=99,
        RESU=(
            _F(
                RESULTAT=evol_endo,
                NOM_CHAM=("VARI_ELGA", "DEPL", "SIEF_ELGA"),
                TOUT_ORDRE="OUI",
                NOM_RESU_MED="beton",
                GROUP_MA="beton",
            ),
            _F(
                RESULTAT=evol_endo,
                NOM_CHAM=("VARI_ELGA", "DEPL", "SIEF_ELGA"),
                TOUT_ORDRE="OUI",
                NOM_RESU_MED="acier",
                GROUP_MA="acier",
            ),
        ),
    )
    DEFI_FICHIER(ACTION="LIBERER", UNITE=99)

    DEFI_FICHIER(UNITE=99, ACCES="OLD", FICHIER=medfile)
    mesh_seq = LIRE_MAILLAGE(UNITE=99)

    seq_result_beton = LIRE_RESU(
        MODELE=model_seq,
        FORMAT="MED",
        UNITE=99,
        TYPE_RESU="EVOL_NOLI",
        COMPORTEMENT=(
            _F(GROUP_MA="acier", RELATION="VISC_ISOT_NL", DEFORMATION="PETIT", REGU_VISC="OUI"),
            _F(GROUP_MA="beton", RELATION="VISC_ISOT_NL", DEFORMATION="PETIT", REGU_VISC="OUI"),
        ),
        CHAM_MATER=fieldmat_seq,
        FORMAT_MED=(
            _F(NOM_RESU="beton", NOM_CHAM="DEPL"),
            _F(NOM_RESU="beton", NOM_CHAM="SIEF_ELGA"),
            _F(NOM_RESU="beton", NOM_CHAM="VARI_ELGA"),
        ),
        TOUT_ORDRE="OUI",
    )

    seq_result_acier = LIRE_RESU(
        MODELE=model_seq,
        FORMAT="MED",
        UNITE=99,
        TYPE_RESU="EVOL_NOLI",
        COMPORTEMENT=(
            _F(GROUP_MA="acier", RELATION="VISC_ISOT_NL", DEFORMATION="PETIT", REGU_VISC="OUI"),
            _F(GROUP_MA="beton", RELATION="VISC_ISOT_NL", DEFORMATION="PETIT", REGU_VISC="OUI"),
        ),
        CHAM_MATER=fieldmat_seq,
        FORMAT_MED=(
            _F(NOM_RESU="acier", NOM_CHAM="DEPL"),
            _F(NOM_RESU="acier", NOM_CHAM="SIEF_ELGA"),
            _F(NOM_RESU="acier", NOM_CHAM="VARI_ELGA"),
        ),
        TOUT_ORDRE="OUI",
    )

FIN()
