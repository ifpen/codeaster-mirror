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

from code_aster.Commands import *
from code_aster import CA

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

from zzzz154a_cmd import MACRO_TEST

test = CA.TestCase()

mesh = CA.Mesh.buildSquare()

ther1 = CREA_CHAMP(
    MAILLAGE=mesh,
    TYPE_CHAM="NOEU_TEMP_R",
    OPERATION="AFFE",
    AFFE=_F(NOM_CMP="TEMP", VALE=1.0, TOUT="OUI"),
)
ther2 = CREA_CHAMP(
    MAILLAGE=mesh,
    TYPE_CHAM="NOEU_TEMP_R",
    OPERATION="AFFE",
    AFFE=_F(NOM_CMP="TEMP", VALE=2.0, TOUT="OUI"),
)
ther3 = CREA_CHAMP(
    MAILLAGE=mesh,
    TYPE_CHAM="NOEU_TEMP_R",
    OPERATION="AFFE",
    AFFE=_F(NOM_CMP="TEMP", VALE=3.0, TOUT="OUI"),
)

ther_dict = CA.ThermalResultDict("ther_dict")
test.assertEqual(len(ther_dict), 0, msg="check len()")
test.assertEqual(ther_dict.getType(), "EVOL_THER_DICT", msg="check type")

# mix direct (global)...
r12 = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_THER",
    AFFE=(
        _F(NOM_CHAM="TEMP", CHAM_GD=ther1, INST=1.0),
        _F(NOM_CHAM="TEMP", CHAM_GD=ther2, INST=2.0),
    ),
)
ther_dict["12"] = r12

# ...and indirect references to DataStructures, see #32978
ther_dict["13"] = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_THER",
    AFFE=(
        _F(NOM_CHAM="TEMP", CHAM_GD=ther1, INST=1.0),
        _F(NOM_CHAM="TEMP", CHAM_GD=ther3, INST=3.0),
    ),
)

test.assertEqual(len(ther_dict), 2, msg="check len()")
test.assertIsNone(ther_dict.get("unknown"), msg="check invalid access")
with test.assertRaises(KeyError):
    _ = ther_dict["unknown"]
with test.assertRaises(TypeError):
    ther_dict["invalid"] = ther1

ch1 = ther_dict["12"].getField("TEMP", 1)
ch3 = ther_dict["13"].getField("TEMP", 2)
test.assertAlmostEqual(max(ch1.getValues()), 1.0, msg="check ther1")
test.assertAlmostEqual(max(ch3.getValues()), 3.0, msg="check ther3")

# test for usage of dict-like objects in commands
dict_test = MACRO_TEST(
    AFFE=(_F(NOM_CAS="ab", RESULTAT=ther_dict["12"]), _F(NOM_CAS="ac", RESULTAT=ther_dict["13"]))
)

# test printing of dict-like objects

IMPR_RESU(RESU=_F(RESULTAT=ther_dict), UNITE=80)

test.printSummary()

CA.close()
