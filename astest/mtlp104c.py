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

import math

####################################################################################################
#
#### Paramètres
#
####################################################################################################

# parametres du modele (cf note EC2M RAP-EDF-CQN02549CE2_05/12/09)
bainCoefB = 0.064380728515158
bainCoefN = 0.216028726635953
#
b_martensite = 0.22414091
n_martensite = 0.14349531
#
R = 2.0
delta_H = 100000.0
C0 = (math.log(10) * R) / delta_H

ac3 = 802.60
tempTemperingBegin = 700.0
tempTemperingEnd = 550.0
tempHold = 610.0

# instants cycle thermique
t0 = 0.0
t1 = 10.0
t2 = 100.0
t3 = 110.0
t4 = 200.0
t5 = t1 + 200.0
t6 = t2 + 200.0
t7 = t3 + 200.0
t8 = t4 + 200.0

####################################################################################################
#
#### Calcul
#
####################################################################################################

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"))

# Get mesh
PRE_GIBI(UNITE_GIBI=19, UNITE_MAILLAGE=21)
meshMeca = LIRE_MAILLAGE(UNITE=21, FORMAT="ASTER")

# Material
INCLUDE(UNITE=12)

materialField = AFFE_MATERIAU(MAILLAGE=meshMeca, AFFE=(_F(TOUT="OUI", MATER=steel),))

# FE model
thermalModel = AFFE_MODELE(
    MAILLAGE=meshMeca, AFFE=_F(TOUT="OUI", PHENOMENE="THERMIQUE", MODELISATION=("AXIS",))
)

# Thermal boundary conditions
therFunction = DEFI_FONCTION(
    NOM_PARA="INST",
    VALE=(
        t0,
        20.0,
        t1,
        2000.0,
        t2,
        20.0,
        t3,
        700.0,
        t4,
        20.0,
        t5,
        2000.0,
        t6,
        20.0,
        t7,
        700.0,
        t8,
        20.0,
    ),
    PROL_GAUCHE="CONSTANT",
    PROL_DROITE="CONSTANT",
)


thermalBCs = AFFE_CHAR_THER_F(MODELE=thermalModel, TEMP_IMPO=(_F(TOUT="OUI", TEMP=therFunction),))

# Time steps
timesteps = DEFI_LIST_REEL(
    DEBUT=t0,
    INTERVALLE=(
        # Cycle 1 ---------------------
        _F(JUSQU_A=t1, NOMBRE=50),
        _F(JUSQU_A=t2, NOMBRE=50),
        # Cycle 2 ---------------------
        _F(JUSQU_A=t3, NOMBRE=50),
        _F(JUSQU_A=t4, NOMBRE=50),
        # Cycle 3 ---------------------
        _F(JUSQU_A=t5, NOMBRE=50),
        _F(JUSQU_A=t6, NOMBRE=50),
        # Cycle 4 ---------------------
        _F(JUSQU_A=t7, NOMBRE=50),
        _F(JUSQU_A=t8, NOMBRE=50),
    ),
)

last_inst = len(timesteps.getValues()) - 1

# resolution thermique
# -----------------------------------------------------
result = THER_NON_LINE(
    MODELE=thermalModel,
    CHAM_MATER=materialField,
    EXCIT=(_F(CHARGE=thermalBCs),),
    INCREMENT=_F(LIST_INST=timesteps, NUME_INST_FIN=last_inst),
    ETAT_INIT=_F(VALE=20.0),
    NEWTON=_F(ITER_LINE_MAXI=3),
    CONVERGENCE=_F(RESI_GLOB_RELA=5.0e-5, ITER_GLOB_MAXI=50),
)

####################################################################################################
#
#### Calcul métallurgique avec revenu
#
####################################################################################################
PHASINIT = CREA_CHAMP(
    TYPE_CHAM="CART_VAR2_R",
    OPERATION="AFFE",
    MAILLAGE=meshMeca,
    AFFE=_F(
        TOUT="OUI",
        NOM_CMP=("V1", "V2", "V3", "V4", "V5", "V6", "V7"),
        VALE=(0.7, 0.0, 0.3, 0.0, 0.0, 1.0, 10.0),
    ),
)


# V1 : FERRITE
# V2 : PERLITE
# V3 : BAINITE
# V4 : BAIN_REVENU
# V5 : MARTENSITE
# V6 : MARTENSITE_REVEN
# V7 : AUSTENITE
# V8 : SUMCOLD
# V9 : TAILLE_GRAIN
# V10 : TEMP
# V11 : TEMP_MARTENSITE
# V12 : CYCL_THER

result = CALC_META(
    reuse=result,
    MODELE=thermalModel,
    CHAM_MATER=materialField,
    RESULTAT=result,
    ETAT_INIT=_F(META_INIT_ELNO=PHASINIT),
    COMPORTEMENT=_F(RELATION="ACIER", TOUT="OUI", LOI_META="WAECKEL"),
    REVENU=_F(RELATION="ACIER_REVENU", TOUT="OUI", LOI_META="JMA"),
    OPTION=("META_ELNO", "DURT_ELNO"),
)


####################################################################################################
#
#### Tests avec revenu
#
####################################################################################################
##########################################################
#
#### Premier cycle
#
##########################################################
# Test à t = t1
TEST_RESU(
    RESU=(
        # BAINITE
        _F(
            RESULTAT=result,
            INST=t1,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V3",
            VALE_CALC=0.0,
        ),
        # BAIN_REVENU
        _F(
            RESULTAT=result,
            INST=t1,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V4",
            VALE_CALC=0.0,
        ),
        # MARTENSITE
        _F(
            RESULTAT=result,
            INST=t1,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V5",
            VALE_CALC=0.0,
        ),
        # MARTENS_REVENU
        _F(
            RESULTAT=result,
            INST=t1,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V6",
            VALE_CALC=0.0,
        ),
    )
)

# Test à t = t2
TEST_RESU(
    RESU=(
        # BAINITE
        _F(
            RESULTAT=result,
            INST=t2,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V3",
            VALE_CALC=0.16373991362659995,
        ),
        # BAIN_REVENU
        _F(
            RESULTAT=result,
            INST=t2,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V4",
            VALE_CALC=0.0,
        ),
        # MARTENSITE
        _F(
            RESULTAT=result,
            INST=t2,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V5",
            VALE_CALC=0.8306074069170569,
        ),
        # MARTENS_REVENU
        _F(
            RESULTAT=result,
            INST=t2,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V6",
            VALE_CALC=0.0,
        ),
    )
)

# Test à t = t3
TEST_RESU(
    RESU=(
        # BAINITE
        _F(
            RESULTAT=result,
            INST=t3,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V3",
            VALE_CALC=0.1485852924694056,
        ),
        # BAIN_REVENU
        _F(
            RESULTAT=result,
            INST=t3,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V4",
            VALE_CALC=0.015154621157194334,
        ),
        # MARTENSITE
        _F(
            RESULTAT=result,
            INST=t3,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V5",
            VALE_CALC=0.6187066145204585,
        ),
        # MARTENS_REVENU
        _F(
            RESULTAT=result,
            INST=t3,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V6",
            VALE_CALC=0.21190079239659837,
        ),
    )
)

# Test à t = t4
TEST_RESU(
    RESU=(
        # BAINITE
        _F(
            RESULTAT=result,
            INST=t4,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V3",
            VALE_CALC=0.14098716459862456,
        ),
        # BAIN_REVENU
        _F(
            RESULTAT=result,
            INST=t4,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V4",
            VALE_CALC=0.022752749027975385,
        ),
        # MARTENSITE
        _F(
            RESULTAT=result,
            INST=t4,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V5",
            VALE_CALC=0.5609993835986875,
        ),
        # MARTENS_REVENU
        _F(
            RESULTAT=result,
            INST=t4,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V6",
            VALE_CALC=0.2696080233183694,
        ),
    )
)

##########################################################
#
#### Second cycle -> Résultats identiques au premier cycle
#
##########################################################
# Test à t = t5
TEST_RESU(
    RESU=(
        # BAINITE
        _F(
            RESULTAT=result,
            INST=t5,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V3",
            VALE_CALC=0.0,
        ),
        # BAIN_REVENU
        _F(
            RESULTAT=result,
            INST=t5,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V4",
            VALE_CALC=0.0,
        ),
        # MARTENSITE
        _F(
            RESULTAT=result,
            INST=t5,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V5",
            VALE_CALC=0.0,
        ),
        # MARTENS_REVENU
        _F(
            RESULTAT=result,
            INST=t5,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V6",
            VALE_CALC=0.0,
        ),
    )
)

# Test à t = t6
TEST_RESU(
    RESU=(
        # BAINITE
        _F(
            RESULTAT=result,
            INST=t6,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V3",
            VALE_CALC=0.16373991362659995,
        ),
        # BAIN_REVENU
        _F(
            RESULTAT=result,
            INST=t6,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V4",
            VALE_CALC=0.0,
        ),
        # MARTENSITE
        _F(
            RESULTAT=result,
            INST=t6,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V5",
            VALE_CALC=0.8306074069170569,
        ),
        # MARTENS_REVENU
        _F(
            RESULTAT=result,
            INST=t6,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V6",
            VALE_CALC=0.0,
        ),
    )
)

# Test à t = t7
TEST_RESU(
    RESU=(
        # BAINITE
        _F(
            RESULTAT=result,
            INST=t7,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V3",
            VALE_CALC=0.1485852924694056,
        ),
        # BAIN_REVENU
        _F(
            RESULTAT=result,
            INST=t7,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V4",
            VALE_CALC=0.015154621157194334,
        ),
        # MARTENSITE
        _F(
            RESULTAT=result,
            INST=t7,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V5",
            VALE_CALC=0.6187066145204585,
        ),
        # MARTENS_REVENU
        _F(
            RESULTAT=result,
            INST=t7,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V6",
            VALE_CALC=0.21190079239659837,
        ),
    )
)

# Test à t = t8
TEST_RESU(
    RESU=(
        # BAINITE
        _F(
            RESULTAT=result,
            INST=t8,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V3",
            VALE_CALC=0.14098716459862456,
        ),
        # BAIN_REVENU
        _F(
            RESULTAT=result,
            INST=t8,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V4",
            VALE_CALC=0.022752749027975385,
        ),
        # MARTENSITE
        _F(
            RESULTAT=result,
            INST=t8,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V5",
            VALE_CALC=0.5609993835986875,
        ),
        # MARTENS_REVENU
        _F(
            RESULTAT=result,
            INST=t8,
            NOM_CHAM="META_ELNO",
            GROUP_MA="MO",
            GROUP_NO="NO1",
            NOM_CMP="V6",
            VALE_CALC=0.2696080233183694,
        ),
    )
)

FIN()
