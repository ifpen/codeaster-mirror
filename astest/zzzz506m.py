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

from code_aster.Utilities import petscFinalize

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

Mail = LIRE_MAILLAGE(FORMAT="MED", PARTITIONNEUR="PTSCOTCH", UNITE=20)

Mail = MODI_MAILLAGE(
    reuse=Mail, MAILLAGE=Mail, ORIE_PEAU=_F(GROUP_MA_PEAU=("Press", "Sym_x", "Sym_y", "Sym_z"))
)

MODI = AFFE_MODELE(AFFE=_F(MODELISATION="3D", PHENOMENE="MECANIQUE", TOUT="OUI"), MAILLAGE=Mail)

mat1 = DEFI_MATERIAU(
    ECRO_LINE=_F(D_SIGM_EPSI=0.0, SY=150.0), ELAS=_F(COEF_AMOR=1.0, E=200000.0, NU=0.3)
)

AFFE = AFFE_MATERIAU(AFFE=_F(MATER=mat1, TOUT="OUI"), MAILLAGE=Mail, MODELE=MODI)

lisi = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1, NOMBRE=2))

LINST = DEFI_LIST_INST(
    DEFI_LIST=_F(LIST_INST=lisi),
    ECHEC=_F(
        ACTION="DECOUPE",
        EVENEMENT="ERREUR",
        SUBD_METHODE="MANUEL",
        SUBD_NIVEAU=3,
        SUBD_PAS=4,
        SUBD_PAS_MINI=0.0,
    ),
)

RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0, 0, 1, 1))

CHAR2 = AFFE_CHAR_MECA(MODELE=MODI, PRES_REP=_F(GROUP_MA=("Press",), PRES=400.0))

CHAR1 = AFFE_CHAR_CINE(
    MECA_IMPO=(
        _F(DX=0.0, GROUP_MA=("Sym_x",)),
        _F(DY=0.0, GROUP_MA=("Sym_y",)),
        _F(DZ=0.0, GROUP_MA=("Sym_z",)),
    ),
    MODELE=MODI,
)

RES = STAT_NON_LINE(
    CHAM_MATER=AFFE,
    COMPORTEMENT=_F(DEFORMATION="PETIT", RELATION="VMIS_ISOT_LINE", TOUT="OUI"),
    CONVERGENCE=_F(RESI_GLOB_MAXI=1e-6, ITER_GLOB_ELAS=25, ITER_GLOB_MAXI=10),
    EXCIT=(
        _F(CHARGE=CHAR1, FONC_MULT=RAMPE, TYPE_CHARGE="FIXE_CSTE"),
        _F(CHARGE=CHAR2, FONC_MULT=RAMPE, TYPE_CHARGE="FIXE_CSTE"),
    ),
    INCREMENT=_F(LIST_INST=lisi),
    INFO=1,
    MODELE=MODI,
    NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_INCR=1, REAC_ITER=1),
    SOLVEUR=_F(
        METHODE="PETSC",
        PRE_COND="HPDDM",
        RESI_RELA=1e-10,
        TYPE_RESOL="HARMO",
        SEUIL=1e-1,
        OPTION_PETSC="-ksp_monitor_true_residual -log_view",
    ),
)

# force PETSc to dump the log
petscFinalize()

RES_NEW = MECA_NON_LINE(
    CHAM_MATER=AFFE,
    COMPORTEMENT=_F(DEFORMATION="PETIT", RELATION="VMIS_ISOT_LINE", TOUT="OUI"),
    CONVERGENCE=_F(RESI_GLOB_MAXI=1e-6, ITER_GLOB_ELAS=25, ITER_GLOB_MAXI=10),
    EXCIT=(
        _F(CHARGE=CHAR1, FONC_MULT=RAMPE, TYPE_CHARGE="FIXE_CSTE"),
        _F(CHARGE=CHAR2, FONC_MULT=RAMPE, TYPE_CHARGE="FIXE_CSTE"),
    ),
    INCREMENT=_F(LIST_INST=lisi),
    INFO=1,
    MODELE=MODI,
    NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_INCR=1, REAC_ITER=1),
    SOLVEUR=_F(
        METHODE="PETSC",
        PRE_COND="GAMG",
        RESI_RELA=1e-10,
        OPTION_PETSC="-ksp_monitor_true_residual -options_view ",
    ),
)

# =========================================================
#          DETERMINATION DE LA REFERENCE
# =========================================================

nbIndexes = RES.getNumberOfIndexes()
test.assertEqual(RES.getNumberOfIndexes(), RES_NEW.getNumberOfIndexes())


# =========================================================
#            REALISATION DES TESTS
# =========================================================

for rank in range(nbIndexes):
    # ON EXTRAIT LES CHAMPS A TESTER au dernier instant
    DEPL_REF = RES.getField("DEPL", rank)
    SIGMA_REF = RES.getField("SIEF_ELGA", rank)
    VARI_REF = RES.getField("VARI_ELGA", rank)

    DEPL = RES_NEW.getField("DEPL", rank)
    SIGMA = RES_NEW.getField("SIEF_ELGA", rank)
    VARI = RES_NEW.getField("VARI_ELGA", rank)

    DIF_DEPL = DEPL_REF - DEPL

    DIF_SIG = SIGMA_REF - SIGMA

    DIF_VAR = VARI_REF - VARI

    TEST_RESU(
        CHAM_ELEM=(
            _F(
                CRITERE="ABSOLU",
                REFERENCE="ANALYTIQUE",
                ORDRE_GRANDEUR=10e-3,
                TYPE_TEST="MIN",
                CHAM_GD=DIF_SIG,
                VALE_CALC=0.0,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="ANALYTIQUE",
                ORDRE_GRANDEUR=10e-3,
                TYPE_TEST="MAX",
                CHAM_GD=DIF_SIG,
                VALE_CALC=0.0,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="ANALYTIQUE",
                ORDRE_GRANDEUR=10e-3,
                TYPE_TEST="MIN",
                CHAM_GD=DIF_VAR,
                VALE_CALC=0.0,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="ANALYTIQUE",
                ORDRE_GRANDEUR=10e-3,
                TYPE_TEST="MAX",
                CHAM_GD=DIF_VAR,
                VALE_CALC=0.0,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
        )
    )

    TEST_RESU(
        CHAM_NO=(
            _F(
                CRITERE="ABSOLU",
                REFERENCE="ANALYTIQUE",
                ORDRE_GRANDEUR=10e-3,
                TYPE_TEST="MIN",
                CHAM_GD=DIF_DEPL,
                VALE_CALC=0.0,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
            _F(
                CRITERE="ABSOLU",
                REFERENCE="ANALYTIQUE",
                ORDRE_GRANDEUR=10e-3,
                TYPE_TEST="MAX",
                CHAM_GD=DIF_DEPL,
                VALE_CALC=0.0,
                VALE_REFE=0.0,
                VALE_ABS="OUI",
            ),
        )
    )

FIN()
