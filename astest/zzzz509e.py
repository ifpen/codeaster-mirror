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

# This testcase is a copy of ssnp140g, but using MECA_NON_LINE.

from code_aster.Commands import *

command = MECA_NON_LINE

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"))

# LECTURE DU MAILLAGE
MAIL = LIRE_MAILLAGE(FORMAT="MED")

MAIL = DEFI_GROUP(
    reuse=MAIL,
    MAILLAGE=MAIL,
    CREA_GROUP_MA=_F(NOM="TOUT_2D", TOUT="OUI", TYPE_MAILLE="2D"),
    CREA_GROUP_NO=(_F(GROUP_MA="34"), _F(GROUP_MA="TOUT_2D")),
)

# AFFECTATION DU MODELE SUR LE MAILLAGE
MODE = AFFE_MODELE(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="D_PLAN"))

# DEFINITION DES PARAMETRES DE LA LOI DE COMPORTEMENT
MA = DEFI_MATERIAU(ELAS=_F(E=70, NU=0.2), ECRO_LINE=_F(D_SIGM_EPSI=2.170542, SY=0.24))

# AFFECTATION DU MATERIAU SUR L ENSEMBLE DU MAILLAGE
MATE = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MA))

# CONDITIONS AUX LIMITES ET CHARGEMENTS
CHAR = AFFE_CHAR_MECA(MODELE=MODE, DDL_IMPO=(_F(GROUP_MA="12", DY=0.0), _F(GROUP_MA="45", DX=0.0)))

CHAR2 = AFFE_CHAR_MECA(MODELE=MODE, DDL_IMPO=_F(GROUP_MA="34", DY=0.3))

F_DEPL = DEFI_FONCTION(NOM_PARA="INST", VALE=(0, 0, 1, 1), PROL_DROITE="LINEAIRE")

# DISCRETISATION EN TEMPS
TFIN = 1.1

# Dans ssnp140g, il y a 2 calculs avec 2 list_inst différents. On prend le second ici.
# L_INST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=50))
# DEFLIST = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=L_INST))

L_INST = DEFI_LIST_REEL(
    DEBUT=0.0,
    INTERVALLE=(
        _F(JUSQU_A=0.005, NOMBRE=1),
        _F(JUSQU_A=0.12, NOMBRE=1),
        _F(JUSQU_A=0.32, NOMBRE=1),
        _F(JUSQU_A=0.6, NOMBRE=1),
        _F(JUSQU_A=1, NOMBRE=1),
        _F(JUSQU_A=1.1, NOMBRE=1),
    ),
)

DEFLIST = DEFI_LIST_INST(
    METHODE="AUTO",
    MODELE=MODE,
    DEFI_LIST=_F(LIST_INST=L_INST),
    ADAPTATION=(
        _F(
            EVENEMENT="TOUT_INST",
            MODE_CALCUL_TPLUS="DELTA_GRANDEUR",
            NOM_CHAM="DEPL",
            NOM_CMP="DY",
            GROUP_NO="TOUT_2D",
            VALE_REF=0.005e0,
        ),
        _F(EVENEMENT="TOUT_INST", MODE_CALCUL_TPLUS="FIXE", PCENT_AUGM=300),
    ),
)

L_ARCH = DEFI_LIST_REEL(VALE=(0.00, 0.32, 0.6, 1.0, 1.1))

# RESOLUTION AVEC LA METHODE ITERATIVE DE NEWTON-RAPHSON
RESUNL = command(
    MODELE=MODE,
    CHAM_MATER=MATE,
    EXCIT=(_F(CHARGE=CHAR), _F(CHARGE=CHAR2, FONC_MULT=F_DEPL)),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
    INCREMENT=_F(LIST_INST=DEFLIST, INST_FIN=1.0),
    NEWTON=_F(REAC_ITER=1),
    # ARCHIVAGE=_F(LIST_INST=L_ARCH),
)

# EXTRACTION DE LA COURBE FORCE APPLIQUEE EN FONCTION DU TEMPS
RESUNL = CALC_CHAMP(reuse=RESUNL, RESULTAT=RESUNL, FORCE="FORC_NODA")

FORC = POST_RELEVE_T(
    ACTION=_F(
        OPERATION="EXTRACTION",
        INTITULE="FORCE",
        RESULTAT=RESUNL,
        NOM_CHAM="FORC_NODA",
        GROUP_NO="34",
        RESULTANTE="DY",
    )
)

FX = RECU_FONCTION(
    TABLE=FORC,
    PARA_X="INST",
    PARA_Y="DY",
    INTERPOL="LIN",
    PROL_DROITE="CONSTANT",
    PROL_GAUCHE="CONSTANT",
)

TEST_FONCTION(
    VALEUR=(
        _F(
            VALE_PARA=0.12,
            VALE_REFE=0.99480191824930,
            VALE_CALC=0.99480191824930,
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            FONCTION=FX,
        ),
        _F(
            VALE_PARA=0.32,
            VALE_CALC=1.565351533864,  # SNL NOOK because the timestep is incorrect at t=0.12
            VALE_REFE=1.565351533864,
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            FONCTION=FX,
        ),
        _F(
            VALE_PARA=0.6,
            VALE_CALC=1.770157957901,
            VALE_REFE=1.770157957901,
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            FONCTION=FX,
        ),
        _F(
            VALE_PARA=1.0,
            VALE_CALC=1.999693075339,
            VALE_REFE=1.999693075339,
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            FONCTION=FX,
        ),
    )
)

FIN()
