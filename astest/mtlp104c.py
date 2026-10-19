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

import numpy
import math
from code_aster.Commands import *
from code_aster import CA

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
#### Fonction Python principale pour le revenu
#
####################################################################################################


def calc_revenu_last_inst1(RESUin__, thermalModel):
    # Extraction de la liste d instant
    L_inst0 = RESUin__.getAccessParameters()["INST"]

    Temp__0 = RESUin__.getField("META_ELNO", para="INST", value=L_inst0[-1])

    # ==========================
    # creation d une table type
    # -------------------------------
    # Rappel : 9 variables internes
    #           V1 : FERRITE
    #           V2 : PERLITE
    #           V3 : BAINITE
    #           V4 : MARTENSITE
    #           V5 : AUSTENITE
    #           V6 : SUMCOLD
    #           V7 : TAILLE_GRAIN
    #           V8 : TEMP
    #           V9 : TEMP_MARTENSITE
    # --------------------------------
    Tab_temp = CREA_TABLE(
        RESU=_F(CHAM_GD=Temp__0, NOM_CMP=("V1", "V2", "V3", "V4", "V5"), TOUT="OUI")
    )
    dico_Z = Tab_temp.EXTR_TABLE().values()
    l_maille0 = dico_Z["MAILLE"]
    l_spoint0 = dico_Z["SOUS_POINT"]
    l_noeud0 = dico_Z["NOEUD"]
    # on peut recuperer les coordonnees des noeuds
    l_coorx = dico_Z["COOR_X"]
    l_coory = dico_Z["COOR_Y"]
    #
    dim = len(l_noeud0)

    # Initialisation des variables du modele de revenu
    l_V1__pred = numpy.zeros(dim)  #    V1  : FERRITE
    l_V2__pred = numpy.zeros(dim)  #    V2  : PERLITE
    l_V3__pred = numpy.zeros(dim)  #    V3  : BAINITE_BRUTE
    l_V4__pred = numpy.zeros(dim)  #    V4  : MARTENSITE_BRUTE
    l_V5__pred = numpy.zeros(dim)  #    V5  : AUSTENITE
    l_V8__pred = numpy.zeros(dim)  #    V8  : TEMP
    l_V10__pred = numpy.zeros(dim)  #    V10 : BAINITE_REVENUE
    l_V11__pred = numpy.zeros(dim)  #    V11 : MARTENSITE_REVENUE

    #
    l_pct_Zr_mart_pred = numpy.zeros(
        dim
    )  #    l_pct_Zr_mart_pred : % MARTENSITE REVENUE par rapport a la MARTENSITE TOTALE, a l inst precedent
    l_pct_Zr_bain_pred = numpy.zeros(
        dim
    )  #    l_pct_Zr_bain_pred : % BAINITE REVENUE par rapport a la BAINITE TOTALE, a l inst precedent
    l_V3__new = numpy.zeros(dim)  #    V3__new            : BAINITE_BRUTE
    l_V4__new = numpy.zeros(dim)  #    V4__new            : MARTENSITE_BRUTE
    #
    l_V10__new = numpy.zeros(dim)  #    V10__new           : BAINITE_REVENUE
    l_V11__new = numpy.zeros(dim)  #    V11__new           : MARTENSITE_REVENUE
    #
    v_refroid_700 = numpy.zeros(dim)  # vitesse de refroidissement a 700°
    #
    inst_pred = L_inst0[0]
    #
    #
    indic = []
    for i in range(0, dim):
        indic.append(0)
    #
    #

    for ii, inst in enumerate(L_inst0):
        #
        delta_inst = inst - inst_pred

        # creation d un champ type de temperature
        Temp__i = RESUin__.getField("META_ELNO", para="INST", value=inst)

        # =================================================================
        # Extraction des phases mettalurgiques a post-traiter
        # =================================================================
        #
        l_V1__i = Temp__i.getValuesWithDescription("V1")[0]
        # V1 : FERRITE
        l_V2__i = Temp__i.getValuesWithDescription("V2")[0]
        # V2 : PERLITE
        l_V3__i = Temp__i.getValuesWithDescription("V3")[0]
        # V3 : BAINITE
        l_V4__i = Temp__i.getValuesWithDescription("V4")[0]
        # V4 : MARTENSITE
        l_V5__i = Temp__i.getValuesWithDescription("V5")[0]
        # V5 : AUSTENITE
        l_V8__i = Temp__i.getValuesWithDescription("V8")[0]
        # V8 : TEMP

        #
        # ==================================================================
        # A programmer : modification des phases et ajout des phase revenue
        # ==================================================================
        #
        if delta_inst == 0.0:
            for i in range(0, dim):
                l_V3__new[i] = l_V2__pred[i]
                l_V4__new[i] = l_V4__pred[i]
                l_V10__new[i] = l_V10__pred[i]
                l_V11__new[i] = l_V11__pred[i]
        else:
            for m, maille in enumerate(l_maille0):
                #
                delta_temp = l_V8__i[m] - l_V8__pred[m]
                pct_Zr_bain = 0.0
                pct_Zr_mart = 0.0
                # - calcul vitesse refroidissement a 700°
                if (l_V8__i[m] < tempTemperingBegin) and (l_V8__pred[m] > tempTemperingBegin):
                    v_refroid_700[m] = delta_temp / delta_inst

                if l_V8__i[m] > ac3:  # AC3 = 846.0
                    indic[m] = 0
                    pct_Zr_mart = 0.0
                    pct_Zr_bain = 0.0

                elif l_V8__i[m] > tempTemperingEnd:

                    if indic[m] == 2:
                        # calcul phase bainite revenu
                        #  - calcul de tau_0_bain
                        tau_0_bain = (-1.0 * math.log(1 - l_pct_Zr_bain_pred[m]) / bainCoefB) ** (
                            1.0 / bainCoefN
                        )
                        ##
                        Pa = 1.0 / (1.0 / (l_V8__i[m] + 273.0) - C0 * math.log(delta_inst))
                        #  - calcul delta_inst_eq_610
                        delta_inst_eq_610 = math.exp(
                            (1.0 / C0) * ((1.0 / (tempHold + 273.0)) - (1.0 / Pa))
                        )
                        #  - calcul de la nouvelle valeur de la phase bainite revenu
                        pct_Zr_bain = 1.0 - math.exp(
                            -1.0 * bainCoefB * (tau_0_bain + delta_inst_eq_610) ** bainCoefN
                        )
                        #
                        # calcul phase martensite revenu
                        #  - calcul de tau_0_mart
                        tau_0_mart = (
                            -1.0 * math.log(1 - l_pct_Zr_mart_pred[m]) / b_martensite
                        ) ** (1.0 / n_martensite)
                        ##
                        Pa = 1.0 / (1.0 / (l_V8__i[m] + 273.0) - C0 * math.log(delta_inst))
                        #  - calcul delta_inst_eq_610
                        delta_inst_eq_610 = math.exp(
                            (1.0 / C0) * ((1.0 / (tempHold + 273.0)) - (1.0 / Pa))
                        )
                        #  - calcul de la nouvelle valeur de la phase bainite revenu
                        pct_Zr_mart = 1.0 - math.exp(
                            -1.0 * b_martensite * (tau_0_mart + delta_inst_eq_610) ** n_martensite
                        )

                else:

                    if indic[m] == 2:
                        pct_Zr_mart = l_pct_Zr_mart_pred[m]
                        pct_Zr_bain = l_pct_Zr_bain_pred[m]
                    elif indic[m] == 1:
                        pct_Zr_mart = 0.0
                        pct_Zr_bain = 0.0
                        if delta_temp > 0.0:
                            indic[m] = 2
                    elif indic[m] == 0:
                        pct_Zr_mart = 0.0
                        pct_Zr_bain = 0.0
                        if l_V4__i[m] > 0.0 or l_V3__i[m] > 0.0:  # A VERIFIER
                            indic[m] = 1

                l_V10__new[m] = l_V3__i[m] * pct_Zr_bain
                l_V3__new[m] = l_V3__i[m] - l_V10__new[m]
                l_V11__new[m] = l_V4__i[m] * pct_Zr_mart
                l_V4__new[m] = l_V4__i[m] - l_V11__new[m]
                # mise a jour des variables a l instant precedent
                l_V8__pred[m] = l_V8__i[m]
                l_V4__pred[m] = l_V4__new[m]
                l_V10__pred[m] = l_V10__new[m]
                l_V11__pred[m] = l_V11__new[m]
                l_pct_Zr_mart_pred[m] = pct_Zr_mart
                l_pct_Zr_bain_pred[m] = pct_Zr_bain

        inst_pred = inst
    #
    # =========================================================================
    # Sauvegarde dernier instant- etape 2 : construction d'une champ a partir de la table TAB
    # =========================================================================

    TAB = CREA_TABLE(
        LISTE=(
            _F(LISTE_K=l_maille0, PARA="MAILLE"),
            _F(LISTE_K=l_noeud0, PARA="NOEUD"),
            _F(LISTE_I=l_spoint0, PARA="SOUS_POINT"),
            _F(LISTE_R=l_V1__i, PARA="V1"),
            _F(LISTE_R=l_V2__i, PARA="V2"),
            _F(LISTE_R=l_V3__new, PARA="V3"),  # BAINITE (brute)
            _F(LISTE_R=l_V4__new, PARA="V4"),  # MARTENSITE (brute)
            _F(LISTE_R=l_V5__i, PARA="V5"),
            _F(LISTE_R=l_V10__new, PARA="V10"),  # BAIN_REVENU
            _F(LISTE_R=l_V11__new, PARA="V11"),  # MARTENS_REVENU
        )
    )

    CHC = CREA_CHAMP(
        TYPE_CHAM="ELNO_VARI_R",
        OPERATION="EXTR",
        TABLE=TAB,
        MODELE=thermalModel,
        OPTION="META_ELNO",
    )

    return CHC


####################################################################################################
#
#### Calcul
#
####################################################################################################

CA.init("--test")

test = CA.TestCase()

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
#### Calcul métallurgique sans revenu
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
# V4 : MARTENSITE
# V5 : AUSTENITE
# V6 : SUMCOLD
# V7 : TAILLE_GRAIN
# V8 : TEMP
# V9 : TEMP_MARTENSITE

result_ = CALC_META(
    reuse=result,
    MODELE=thermalModel,
    CHAM_MATER=materialField,
    RESULTAT=result,
    ETAT_INIT=_F(META_INIT_ELNO=PHASINIT),
    COMPORTEMENT=_F(RELATION="ACIER", TOUT="OUI", LOI_META="WAECKEL"),
    OPTION=("META_ELNO", "DURT_ELNO"),
)

metaWithTempering = calc_revenu_last_inst1(result_, thermalModel)

####################################################################################################
#
#### Calcul métallurgique avec revenu
#
####################################################################################################

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

# Test à t = t8 pour le script python
# en utilisant les mêmes valeurs de non-regression que le résultat issu de CALC_META
temp_t8 = result.getField("META_ELNO", para="INST", value=t8)
#
l_V3 = temp_t8.getValuesWithDescription("V3")[0]
value_V3 = l_V3[1]
#
l_V4 = temp_t8.getValuesWithDescription("V4")[0]
value_V4 = l_V4[1]
#
l_V5 = temp_t8.getValuesWithDescription("V5")[0]
value_V5 = l_V5[1]
#
l_V6 = temp_t8.getValuesWithDescription("V6")[0]
value_V6 = l_V6[1]
#
test.assertAlmostEqual(abs(metaWithTempering.getValues()[2] - value_V3), 0.0)
test.assertAlmostEqual(abs(metaWithTempering.getValues()[3] - value_V5), 0.0)
test.assertAlmostEqual(abs(metaWithTempering.getValues()[9] - value_V4), 0.0)
test.assertAlmostEqual(abs(metaWithTempering.getValues()[10] - value_V6), 0.0)
#
test.printSummary()

FIN()
