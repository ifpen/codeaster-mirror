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
from code_aster.CA import MPI

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()

mesh = CA.ParallelMesh()
mesh.readMedFile("zzzz155n.med")

MAT = DEFI_MATERIAU(ELAS=_F(E=9.84e4, NU=0.3))

MODE = AFFE_MODELE(
    MAILLAGE=mesh, AFFE=(_F(GROUP_MA=("vol1", "vol2"), PHENOMENE="MECANIQUE", MODELISATION="3D"),)
)

MATE = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=(_F(GROUP_MA=("vol1", "vol2"), MATER=MAT),))

BLOC = AFFE_CHAR_CINE(MODELE=MODE, MECA_IMPO=(_F(GROUP_MA=("gauche"), DX=0, DY=0, DZ=0.0),))

CHAR = AFFE_CHAR_CINE(MODELE=MODE, MECA_IMPO=(_F(GROUP_MA=("droite"), DX=0.1, DY=0.2, DZ=0),))

LIAISON = AFFE_CHAR_MECA(
    MODELE=MODE,
    INFO=2,
    LIAISON_GROUP=(
        _F(
            GROUP_NO_1="collage",
            GROUP_NO_2="collage2",
            DDL_1="DX",
            DDL_2="DX",
            COEF_MULT_1=1,
            COEF_MULT_2=-1,
            COEF_IMPO=0.3,
        ),
    ),
    DOUBLE_LAGRANGE="NON",
)

LINST = DEFI_LIST_REEL(DEBUT=0, INTERVALLE=(_F(JUSQU_A=1, NOMBRE=1), _F(JUSQU_A=2, NOMBRE=1)))

LINST2 = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=LINST), ECHEC=_F(SUBD_NIVEAU=5, SUBD_PAS=10))

LINEDEPL = DEFI_FONCTION(NOM_PARA="INST", ABSCISSE=(0, 1, 2), ORDONNEE=(0, 8, 10))

RESU = STAT_NON_LINE(
    MODELE=MODE,
    CHAM_MATER=MATE,
    EXCIT=(_F(CHARGE=BLOC), _F(CHARGE=LIAISON), _F(CHARGE=CHAR)),
    INCREMENT=_F(LIST_INST=LINST2, INST_FIN=1),
    COMPORTEMENT=(_F(DEFORMATION="PETIT", RELATION="ELAS", GROUP_MA=("vol1", "vol2")),),
    NEWTON=_F(),
    CONVERGENCE=_F(ITER_GLOB_MAXI=100, RESI_GLOB_MAXI=1),
    SOLVEUR=_F(ELIM_LAGR="NON"),
)

TEST_RESU(
    RESU=(
        _F(NUME_ORDRE=1, PARA="INST", RESULTAT=RESU, VALE_CALC=1.0),
        _F(NUME_ORDRE=1, PARA="ITER_GLOB", VALE_CALC_I=1, RESULTAT=RESU, CRITERE="ABSOLU"),
    )
)


mesh = CA.Mesh()
mesh.readMedFile("zzzz155n.med")

MAT = DEFI_MATERIAU(ELAS=_F(E=9.84e4, NU=0.3))

MODE = AFFE_MODELE(
    MAILLAGE=mesh, AFFE=(_F(GROUP_MA=("vol1", "vol2"), PHENOMENE="MECANIQUE", MODELISATION="3D"),)
)

MATE = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=(_F(GROUP_MA=("vol1", "vol2"), MATER=MAT),))

BLOC = AFFE_CHAR_CINE(MODELE=MODE, MECA_IMPO=(_F(GROUP_MA=("gauche"), DX=0, DY=0, DZ=0.0),))

CHAR = AFFE_CHAR_CINE(MODELE=MODE, MECA_IMPO=(_F(GROUP_MA=("droite"), DX=0.1, DY=0.2, DZ=0),))

LIAISON = AFFE_CHAR_MECA(
    MODELE=MODE,
    INFO=2,
    LIAISON_GROUP=(
        _F(
            GROUP_NO_1="collage",
            GROUP_NO_2="collage2",
            DDL_1="DX",
            DDL_2="DX",
            COEF_MULT_1=1,
            COEF_MULT_2=-1,
            COEF_IMPO=0.3,
        ),
    ),
    DOUBLE_LAGRANGE="NON",
)

LINST = DEFI_LIST_REEL(DEBUT=0, INTERVALLE=(_F(JUSQU_A=1, NOMBRE=1), _F(JUSQU_A=2, NOMBRE=1)))

LINST2 = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=LINST), ECHEC=_F(SUBD_NIVEAU=5, SUBD_PAS=10))

LINEDEPL = DEFI_FONCTION(NOM_PARA="INST", ABSCISSE=(0, 1, 2), ORDONNEE=(0, 8, 10))

RESU = STAT_NON_LINE(
    MODELE=MODE,
    CHAM_MATER=MATE,
    EXCIT=(_F(CHARGE=BLOC), _F(CHARGE=LIAISON), _F(CHARGE=CHAR)),
    INCREMENT=_F(LIST_INST=LINST2, INST_FIN=1),
    COMPORTEMENT=(_F(DEFORMATION="PETIT", RELATION="ELAS", GROUP_MA=("vol1", "vol2")),),
    NEWTON=_F(),
    CONVERGENCE=_F(ITER_GLOB_MAXI=100, RESI_GLOB_MAXI=1),
    SOLVEUR=_F(ELIM_LAGR="NON"),
)

TEST_RESU(
    RESU=(
        _F(NUME_ORDRE=1, PARA="INST", RESULTAT=RESU, VALE_CALC=1.0),
        _F(NUME_ORDRE=1, PARA="ITER_GLOB", VALE_CALC_I=1, RESULTAT=RESU, CRITERE="ABSOLU"),
    )
)


FIN()
