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

from cataelem.Tools.base_objects import MeshType, Elrefe, objects_from_context

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#      ATTENTION   ATTENTION   ATTENTION   ATTENTION   ATTENTION       %
#                                                                      %
#   IL NE FAUT PAS MODIFIER LE NOMBRE OU L'ORDRE DES TYPE_MAILLE       %
#   SANS METTRE A JOUR LES ROUTINES LRMTYP, CMLQLQ, IRMMMA, IRMMFA,    %
#   OP0150, LRMMMA, LRMMFA, LRMMDI, LRMHDF, LRCAME, IRMHDF, IRCMPR,    %
#   IRCMPE, IRCAME                                                     %
#                                                                      %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# ------------------------------------------------------------
POI1 = MeshType(nbno=1, dim=0, code="POI")

PO1 = Elrefe()
PO1.addLocation("NOEU", 1)
PO1.addLocation("NOEU_S", 1)
PO1.addLocation("FPG1", 1)
POI1.addElrefe(PO1)


# ------------------------------------------------------------
SEG2 = MeshType(nbno=2, dim=1, code="SE2")

SE2 = Elrefe()
SE2.addLocation("NOEU", 2)
SE2.addLocation("NOEU_S", 2)
SE2.addLocation("FPG1", 1)
SE2.addLocation("FPG2", 2)
SE2.addLocation("FPG2NOS", 4)
SE2.addLocation("FPG3", 3)
SE2.addLocation("FPG3NOS", 5)
SE2.addLocation("FPG4", 4)
SE2.addLocation("SIMP", 3)
SE2.addLocation("SIMP1", 5)
SE2.addLocation("COTES", 4)
SE2.addLocation("COTES1", 5)
SE2.addLocation("COTES2", 10)
SEG2.addElrefe(SE2)

CABPOU = Elrefe()
CABPOU.addLocation("FPG1", 1)
SEG2.addElrefe(CABPOU)

THCOSE2 = Elrefe()
THCOSE2.addLocation("FPG1", 3)
SEG2.addElrefe(THCOSE2)


# ------------------------------------------------------------
SEG22 = MeshType(nbno=4, dim=1, code="S22")


# ------------------------------------------------------------
SEG3 = MeshType(nbno=3, dim=1, code="SE3")

SE3 = Elrefe()
SE3.addLocation("NOEU", 3)
SE3.addLocation("NOEU_S", 2)
SE3.addLocation("FPG1", 1)
SE3.addLocation("FPG2", 2)
SE3.addLocation("FPG2NOS", 4)
SE3.addLocation("FPG3", 3)
SE3.addLocation("FPG3NOS", 5)
SE3.addLocation("FPG4", 4)
SE3.addLocation("SIMP", 3)
SE3.addLocation("SIMP1", 5)
SE3.addLocation("COTES", 4)
SE3.addLocation("COTES1", 5)
SE3.addLocation("COTES2", 10)
SE3.addLocation("FPG5", 5)
SE3.addLocation("FPG6", 6)
SE3.addLocation("FPG7", 7)
SE3.addLocation("FPG8", 8)
SEG3.addElrefe(SE3)

THCOSE3 = Elrefe()
THCOSE3.addLocation("FPG1", 3)
SEG3.addElrefe(THCOSE3)


# ------------------------------------------------------------
SEG33 = MeshType(nbno=6, dim=1, code="S33")


# ------------------------------------------------------------
SEG4 = MeshType(nbno=4, dim=1, code="SE4")

SE4 = Elrefe()
SE4.addLocation("NOEU", 4)
SE4.addLocation("NOEU_S", 2)
SE4.addLocation("FPG1", 1)
SE4.addLocation("FPG2", 2)
SE4.addLocation("FPG3", 3)
SE4.addLocation("FPG4", 4)
SE4.addLocation("FPG5", 5)
SE4.addLocation("FPG6", 6)
SEG4.addElrefe(SE4)


# ------------------------------------------------------------
TRIA3 = MeshType(nbno=3, dim=2, code="TR3")

TR3 = Elrefe()
TR3.addLocation("NOEU", 3)
TR3.addLocation("NOEU_S", 3)
TR3.addLocation("FPG1", 1)
TR3.addLocation("FPG3", 3)
TR3.addLocation("FPG3NOS", 6)
TR3.addLocation("COT3", 3)
TR3.addLocation("FPG4", 4)
TR3.addLocation("FPG6", 6)
TR3.addLocation("FPG7", 7)
TR3.addLocation("FPG12", 12)
TR3.addLocation("FPG13", 13)
TR3.addLocation("FPG16", 16)
TR3.addLocation("XFEM36", 36)
TR3.addLocation("XFEM72", 72)
TR3.addLocation("SIMP", 6)
TRIA3.addElrefe(TR3)


# ------------------------------------------------------------
TRIA33 = MeshType(nbno=6, dim=2, code="T33")


# ------------------------------------------------------------
TRIA6 = MeshType(nbno=6, dim=2, code="TR6")

TR6 = Elrefe()
TR6.addLocation("NOEU", 6)
TR6.addLocation("NOEU_S", 3)
TR6.addLocation("FPG1", 1)
TR6.addLocation("FPG3", 3)
TR6.addLocation("FPG3NOS", 6)
TR6.addLocation("FPG4", 4)
TR6.addLocation("FPG6", 6)
TR6.addLocation("FPG7", 7)
TR6.addLocation("FPG12", 12)
TR6.addLocation("FPG13", 13)
TR6.addLocation("FPG16", 16)
TR6.addLocation("XFEM36", 36)
TR6.addLocation("XFEM72", 72)
TR6.addLocation("SIMP", 6)
TRIA6.addElrefe(TR6)


# ------------------------------------------------------------
TRIA66 = MeshType(nbno=12, dim=2, code="T66")


# ------------------------------------------------------------
TRIA7 = MeshType(nbno=7, dim=2, code="TR7")

MEC3TR7H = Elrefe()
MEC3TR7H.addLocation("FPG1", 7)
TRIA7.addElrefe(MEC3TR7H)

TR7 = Elrefe()
TR7.addLocation("NOEU", 7)
TR7.addLocation("NOEU_S", 3)
TR7.addLocation("FPG1", 1)
TR7.addLocation("FPG3", 3)
TR7.addLocation("FPG4", 4)
TR7.addLocation("FPG6", 6)
TR7.addLocation("FPG7", 7)
TR7.addLocation("FPG12", 12)
TR7.addLocation("FPG13", 13)
TR7.addLocation("FPG16", 16)
TR7.addLocation("FPG19", 19)
TR7.addLocation("FPG25", 25)
TR7.addLocation("FPG28", 28)
TR7.addLocation("FPG33", 33)
TR7.addLocation("FPG37", 37)
TR7.addLocation("FPG42", 42)
TRIA7.addElrefe(TR7)


# ------------------------------------------------------------
QUAD4 = MeshType(nbno=4, dim=2, code="QU4")

QU4 = Elrefe()
QU4.addLocation("NOEU", 4)
QU4.addLocation("NOEU_S", 4)
QU4.addLocation("FPG1", 1)
QU4.addLocation("FIS2", 2)
QU4.addLocation("FPG4", 4)
QU4.addLocation("FPG4NOS", 8)
QU4.addLocation("FPG9", 9)
QU4.addLocation("FPG9COQ", 9)
QU4.addLocation("XFEM72", 72)
QU4.addLocation("XFEM144", 144)
QU4.addLocation("SIMP", 9)
QUAD4.addElrefe(QU4)


# ------------------------------------------------------------
QUAD44 = MeshType(nbno=8, dim=2, code="Q44")


# ------------------------------------------------------------
QUAD8 = MeshType(nbno=8, dim=2, code="QU8")

QU8 = Elrefe()
QU8.addLocation("NOEU", 8)
QU8.addLocation("NOEU_S", 4)
QU8.addLocation("FPG1", 1)
QU8.addLocation("FPG4", 4)
QU8.addLocation("FPG4NOS", 8)
QU8.addLocation("FPG9", 9)
QU8.addLocation("FPG16", 16)
QU8.addLocation("FPG9COQ", 9)
QU8.addLocation("XFEM72", 72)
QU8.addLocation("XFEM144", 144)
QUAD8.addElrefe(QU8)


# ------------------------------------------------------------
QUAD88 = MeshType(nbno=16, dim=2, code="Q88")


# ------------------------------------------------------------
QUAD9 = MeshType(nbno=9, dim=2, code="QU9")

QU9 = Elrefe()
QU9.addLocation("NOEU", 9)
QU9.addLocation("NOEU_S", 4)
QU9.addLocation("FPG1", 1)
QU9.addLocation("FPG4", 4)
QU9.addLocation("FPG5", 5)
QU9.addLocation("FPG9", 9)
QU9.addLocation("FPG9COQ", 9)
QU9.addLocation("FPG16", 16)
QU9.addLocation("FPG25", 25)
QU9.addLocation("FPG36", 36)
QU9.addLocation("FPG49", 49)
QU9.addLocation("FPG64", 64)
QUAD9.addElrefe(QU9)

MEC3QU9H = Elrefe()
MEC3QU9H.addLocation("FPG1", 9)
QUAD9.addElrefe(MEC3QU9H)


# ------------------------------------------------------------
QUAD99 = MeshType(nbno=18, dim=2, code="Q99")


# ------------------------------------------------------------
TETRA4 = MeshType(nbno=4, dim=3, code="TE4")

TE4 = Elrefe()
TE4.addLocation("NOEU", 4)
TE4.addLocation("NOEU_S", 4)
TE4.addLocation("FPG1", 1)
TE4.addLocation("FPG4", 4)
TE4.addLocation("FPG4NOS", 8)
TE4.addLocation("FPG5", 5)
TE4.addLocation("FPG11", 11)
TE4.addLocation("FPG15", 15)
TE4.addLocation("FPG24", 24)
TE4.addLocation("XFEM90", 90)
TE4.addLocation("XFEM180", 180)
TE4.addLocation("XFEM270", 270)
TE4.addLocation("XFEM360", 360)
TETRA4.addElrefe(TE4)


# ------------------------------------------------------------
TETRA10 = MeshType(nbno=10, dim=3, code="T10")

T10 = Elrefe()
T10.addLocation("NOEU", 10)
T10.addLocation("NOEU_S", 4)
T10.addLocation("FPG1", 1)
T10.addLocation("FPG4", 4)
T10.addLocation("FPG4NOS", 8)
T10.addLocation("FPG5", 5)
T10.addLocation("FPG15", 15)
T10.addLocation("XFEM90", 90)
T10.addLocation("XFEM270", 270)
T10.addLocation("XFEM540", 540)
T10.addLocation("XFEM810", 810)
TETRA10.addElrefe(T10)


# ------------------------------------------------------------
PENTA6 = MeshType(nbno=6, dim=3, code="PE6")

PE6 = Elrefe()
PE6.addLocation("NOEU", 6)
PE6.addLocation("NOEU_S", 6)
PE6.addLocation("FPG1", 1)
PE6.addLocation("FPG6", 6)
PE6.addLocation("FPG6B", 6)
PE6.addLocation("FPG6NOS", 12)
PE6.addLocation("FPG8", 8)
PE6.addLocation("FPG21", 21)
PE6.addLocation("FPG29", 29)
PE6.addLocation("XFEM72", 72)
PE6.addLocation("XFEM240", 240)
PE6.addLocation("XFEM480", 480)
PE6.addLocation("XFEM720", 720)
PE6.addLocation("XFEM960", 960)
PENTA6.addElrefe(PE6)


# ------------------------------------------------------------
PENTA15 = MeshType(nbno=15, dim=3, code="P15")

P15 = Elrefe()
P15.addLocation("NOEU", 15)
P15.addLocation("NOEU_S", 6)
P15.addLocation("FPG1", 1)
P15.addLocation("FPG6", 6)
P15.addLocation("FPG6B", 6)
P15.addLocation("FPG6NOS", 12)
P15.addLocation("FPG8", 8)
P15.addLocation("FPG21", 21)
P15.addLocation("FPG29", 29)
P15.addLocation("XFEM240", 240)
P15.addLocation("XFEM720", 720)
P15.addLocation("XFEM1440", 1440)
P15.addLocation("XFEM2160", 2160)
PENTA15.addElrefe(P15)

# ------------------------------------------------------------
PENTA18 = MeshType(nbno=18, dim=3, code="P18")

P18 = Elrefe()
P18.addLocation("NOEU", 18)
P18.addLocation("NOEU_S", 6)
P18.addLocation("FPG1", 1)
P18.addLocation("FPG6", 6)
P18.addLocation("FPG6B", 6)
P18.addLocation("FPG6NOS", 12)
P18.addLocation("FPG8", 8)
P18.addLocation("FPG21", 21)
P18.addLocation("FPG29", 29)
P18.addLocation("FPG95", 95)
PENTA18.addElrefe(P18)

# ------------------------------------------------------------
PYRAM5 = MeshType(nbno=5, dim=3, code="PY5")

PY5 = Elrefe()
PY5.addLocation("NOEU", 5)
PY5.addLocation("NOEU_S", 5)
PY5.addLocation("FPG1", 1)
PY5.addLocation("FPG1B", 1)
PY5.addLocation("FPG5", 5)
PY5.addLocation("FPG6", 6)
PY5.addLocation("FPG10", 10)
PY5.addLocation("FPG15", 15)
PY5.addLocation("FPG24", 24)
PY5.addLocation("FPG31", 31)
PY5.addLocation("XFEM180", 180)
PY5.addLocation("XFEM360", 360)
PY5.addLocation("XFEM540", 540)
PY5.addLocation("XFEM720", 720)
PYRAM5.addElrefe(PY5)


# ------------------------------------------------------------
PYRAM13 = MeshType(nbno=13, dim=3, code="P13")

P13 = Elrefe()
P13.addLocation("NOEU", 13)
P13.addLocation("NOEU_S", 5)
P13.addLocation("FPG1", 1)
P13.addLocation("FPG1B", 1)
P13.addLocation("FPG5", 5)
P13.addLocation("FPG6", 6)
P13.addLocation("FPG10", 10)
P13.addLocation("FPG15", 15)
P13.addLocation("FPG24", 24)
P13.addLocation("FPG31", 31)
P13.addLocation("XFEM180", 180)
P13.addLocation("XFEM540", 540)
P13.addLocation("XFEM1080", 1080)
P13.addLocation("XFEM1620", 1620)
PYRAM13.addElrefe(P13)


# ------------------------------------------------------------
HEXA8 = MeshType(nbno=8, dim=3, code="HE8")

HE8 = Elrefe()
HE8.addLocation("NOEU", 8)
HE8.addLocation("NOEU_S", 8)
HE8.addLocation("FPG1", 1)
HE8.addLocation("FPG8", 8)
HE8.addLocation("FPG8NOS", 16)
HE8.addLocation("FPG27", 27)
HE8.addLocation("FPG64", 64)
HE8.addLocation("XFEM480", 480)
HE8.addLocation("XFEM960", 960)
HE8.addLocation("XFEM1440", 1440)
HE8.addLocation("XFEM1920", 1920)
HEXA8.addElrefe(HE8)

POHOH8 = Elrefe()
POHOH8.addLocation("FPG8", 8)
HEXA8.addElrefe(POHOH8)

# ------------------------------------------------------------
HEXA20 = MeshType(nbno=20, dim=3, code="H20")

H20 = Elrefe()
H20.addLocation("NOEU", 20)
H20.addLocation("NOEU_S", 8)
H20.addLocation("FPG1", 1)
H20.addLocation("FPG8", 8)
H20.addLocation("FPG8NOS", 16)
H20.addLocation("FPG27", 27)
H20.addLocation("FPG64", 64)
H20.addLocation("XFEM480", 480)
H20.addLocation("XFEM1440", 1440)
H20.addLocation("XFEM2880", 2880)
H20.addLocation("XFEM4320", 4320)
HEXA20.addElrefe(H20)

POHOH20 = Elrefe()
POHOH20.addLocation("FPG27", 27)
HEXA20.addElrefe(POHOH20)


# ------------------------------------------------------------
HEXA27 = MeshType(nbno=27, dim=3, code="H27")

H27 = Elrefe()
H27.addLocation("NOEU", 27)
H27.addLocation("NOEU_S", 8)
H27.addLocation("FPG1", 1)
H27.addLocation("FPG8", 8)
H27.addLocation("FPG7", 7)
H27.addLocation("FPG27", 27)
H27.addLocation("FPG64", 64)
H27.addLocation("FPG125", 125)
H27.addLocation("FPG216", 216)
H27.addLocation("FPG343", 343)
H27.addLocation("FPG512", 512)
HEXA27.addElrefe(H27)

H10 = Elrefe()
H10.addLocation("FIS3", 3)
H10.addLocation("FIS4", 4)
H10.addLocation("FIS2", 2)
HEXA27.addElrefe(H10)

# ------------------------------------------------------------
TETRA15 = MeshType(nbno=15, dim=3, code="T15")

T15 = Elrefe()
T15.addLocation("NOEU", 15)
T15.addLocation("NOEU_S", 4)
T15.addLocation("FPG1", 1)
T15.addLocation("FPG4", 4)
T15.addLocation("FPG4NOS", 8)
T15.addLocation("FPG5", 5)
T15.addLocation("FPG11", 11)
T15.addLocation("FPG15", 15)
T15.addLocation("FPG24", 24)
T15.addLocation("FPG35", 35)
T15.addLocation("FPG46", 46)
T15.addLocation("FPG59", 59)
T15.addLocation("FPG74", 74)
T15.addLocation("FPG94", 94)
T15.addLocation("FPG117", 117)
T15.addLocation("FPG144", 144)
T15.addLocation("FPG204", 204)
TETRA15.addElrefe(T15)

# ------------------------------------------------------------
PENTA21 = MeshType(nbno=21, dim=3, code="P21")

P21 = Elrefe()
P21.addLocation("NOEU", 21)
P21.addLocation("NOEU_S", 6)
P21.addLocation("FPG1", 1)
P21.addLocation("FPG6", 6)
P21.addLocation("FPG6B", 6)
P21.addLocation("FPG6NOS", 12)
P21.addLocation("FPG8", 8)
P21.addLocation("FPG21", 21)
P21.addLocation("FPG29", 29)
P21.addLocation("FPG95", 95)
PENTA21.addElrefe(P21)

# ------------------------------------------------------------
PYRAM19 = MeshType(nbno=19, dim=3, code="P19")

P19 = Elrefe()
P19.addLocation("NOEU", 19)
P19.addLocation("NOEU_S", 5)
P19.addLocation("FPG1", 1)
P19.addLocation("FPG1B", 1)
P19.addLocation("FPG5", 5)
P19.addLocation("FPG6", 6)
P19.addLocation("FPG10", 10)
P19.addLocation("FPG15", 15)
P19.addLocation("FPG24", 24)
P19.addLocation("FPG31", 31)
P19.addLocation("FPG47", 47)
P19.addLocation("XFEM180", 180)
P19.addLocation("XFEM540", 540)
P19.addLocation("XFEM1080", 1080)
P19.addLocation("XFEM1620", 1620)
PYRAM19.addElrefe(P19)

# ------------------------------------------------------------
TR3QU4 = MeshType(nbno=7, dim=2, code="TQ7")


# ------------------------------------------------------------
QU4TR3 = MeshType(nbno=7, dim=2, code="QT7")


# ------------------------------------------------------------
TR6TR3 = MeshType(nbno=9, dim=2, code="TT1")


# ------------------------------------------------------------
TR3TR6 = MeshType(nbno=9, dim=2, code="TT2")


# ------------------------------------------------------------
TR6QU4 = MeshType(nbno=10, dim=2, code="TQ1")


# ------------------------------------------------------------
QU4TR6 = MeshType(nbno=10, dim=2, code="QT1")


# ------------------------------------------------------------
TR6QU8 = MeshType(nbno=14, dim=2, code="TQ2")


# ------------------------------------------------------------
QU8TR6 = MeshType(nbno=14, dim=2, code="QT2")


# ------------------------------------------------------------
TR6QU9 = MeshType(nbno=15, dim=2, code="TQ3")


# ------------------------------------------------------------
QU9TR6 = MeshType(nbno=15, dim=2, code="QT3")


# ------------------------------------------------------------
QU8TR3 = MeshType(nbno=11, dim=2, code="QT4")


# ------------------------------------------------------------
TR3QU8 = MeshType(nbno=11, dim=2, code="TQ4")


# ------------------------------------------------------------
QU8QU4 = MeshType(nbno=12, dim=2, code="QQ1")


# ------------------------------------------------------------
QU4QU8 = MeshType(nbno=12, dim=2, code="QQ2")


# ------------------------------------------------------------
QU8QU9 = MeshType(nbno=17, dim=2, code="QQ3")


# ------------------------------------------------------------
QU9QU8 = MeshType(nbno=17, dim=2, code="QQ4")


# ------------------------------------------------------------
QU9QU4 = MeshType(nbno=13, dim=2, code="QQ5")


# ------------------------------------------------------------
QU4QU9 = MeshType(nbno=13, dim=2, code="QQ6")


# ------------------------------------------------------------
QU9TR3 = MeshType(nbno=12, dim=2, code="QT5")


# ------------------------------------------------------------
TR3QU9 = MeshType(nbno=12, dim=2, code="TQ5")


# ------------------------------------------------------------
SEG32 = MeshType(nbno=5, dim=1, code="S32")


# ------------------------------------------------------------
SEG23 = MeshType(nbno=5, dim=1, code="S23")


# ------------------------------------------------------------
QU4QU4 = MeshType(nbno=8, dim=2, code="QQ7")


# ------------------------------------------------------------
TR3TR3 = MeshType(nbno=6, dim=2, code="TT3")


# ------------------------------------------------------------
HE8HE8 = MeshType(nbno=16, dim=3, code="HH1")


# ------------------------------------------------------------
PE6PE6 = MeshType(nbno=12, dim=3, code="PP1")


# ------------------------------------------------------------
TE4TE4 = MeshType(nbno=8, dim=3, code="TT4")


# ------------------------------------------------------------
QU8QU8 = MeshType(nbno=16, dim=2, code="QQ8")


# ------------------------------------------------------------
TR6TR6 = MeshType(nbno=12, dim=2, code="TT5")


# ------------------------------------------------------------
SE2TR3 = MeshType(nbno=5, dim=2, code="ST1")


# ------------------------------------------------------------
SE2TR6 = MeshType(nbno=8, dim=2, code="ST2")


# ------------------------------------------------------------
SE2QU4 = MeshType(nbno=6, dim=2, code="SQ1")


# ------------------------------------------------------------
SE2QU8 = MeshType(nbno=10, dim=2, code="SQ2")


# ------------------------------------------------------------
SE2QU9 = MeshType(nbno=11, dim=2, code="SQ3")


# ------------------------------------------------------------
SE3TR3 = MeshType(nbno=6, dim=2, code="ST3")


# ------------------------------------------------------------
SE3TR6 = MeshType(nbno=9, dim=2, code="ST4")


# ------------------------------------------------------------
SE3QU4 = MeshType(nbno=7, dim=2, code="SQ4")


# ------------------------------------------------------------
SE3QU8 = MeshType(nbno=11, dim=2, code="SQ5")


# ------------------------------------------------------------
SE3QU9 = MeshType(nbno=12, dim=2, code="SQ6")


# ------------------------------------------------------------
HEXA9 = MeshType(nbno=9, dim=3, code="HE9")
HE9 = Elrefe()
HE9.addLocation("NOEU", 9)
HE9.addLocation("NOEU_S", 8)
HE9.addLocation("LOB5", 5)
HE9.addLocation("LOB7", 7)
HE9.addLocation("FPG1", 1)
HE9.addLocation("FPG8", 8)
HEXA9.addElrefe(HE9)


# ------------------------------------------------------------
PENTA7 = MeshType(nbno=7, dim=3, code="PE7")
PE7 = Elrefe()
PE7.addLocation("NOEU", 7)
PE7.addLocation("NOEU_S", 6)
PE7.addLocation("LOB5", 5)
PE7.addLocation("LOB7", 7)
PE7.addLocation("FPG1", 1)
PENTA7.addElrefe(PE7)

# ------------------------------------------------------------
TR3SE2 = MeshType(nbno=5, dim=1, code="TS1")

# ------------------------------------------------------------
TR3SE3 = MeshType(nbno=6, dim=1, code="TS2")

# ------------------------------------------------------------
TR6SE2 = MeshType(nbno=8, dim=1, code="TS3")

# ------------------------------------------------------------
TR6SE3 = MeshType(nbno=9, dim=1, code="TS4")

# ------------------------------------------------------------
QU4SE2 = MeshType(nbno=6, dim=1, code="QS1")

# ------------------------------------------------------------
QU4SE3 = MeshType(nbno=7, dim=1, code="QS2")

# ------------------------------------------------------------
QU8SE2 = MeshType(nbno=10, dim=1, code="QS3")

# ------------------------------------------------------------
QU8SE3 = MeshType(nbno=11, dim=1, code="QS4")

# ------------------------------------------------------------
QU9SE2 = MeshType(nbno=11, dim=1, code="QS5")

# ------------------------------------------------------------
QU9SE3 = MeshType(nbno=12, dim=1, code="QS6")

# ------------------------------------------------------------
QU4TR7 = MeshType(nbno=11, dim=2, code="QT8")

# ------------------------------------------------------------
QU8TR7 = MeshType(nbno=15, dim=2, code="QT9")

# ------------------------------------------------------------
QU9TR7 = MeshType(nbno=16, dim=2, code="QT0")

# ------------------------------------------------------------
TR3TR7 = MeshType(nbno=10, dim=2, code="TT6")

# ------------------------------------------------------------
TR6TR7 = MeshType(nbno=13, dim=2, code="TT7")

# ------------------------------------------------------------
TRIA10 = MeshType(nbno=10, dim=2, code="TR1")


TR1 = Elrefe()
TR1.addLocation("NOEU", 10)
TR1.addLocation("NOEU_S", 3)
TR1.addLocation("FPG1", 1)
TR1.addLocation("FPG3", 3)
TR1.addLocation("FPG4", 4)
TR1.addLocation("FPG6", 6)
TR1.addLocation("FPG7", 7)
TR1.addLocation("FPG12", 12)
TR1.addLocation("FPG13", 13)
TR1.addLocation("FPG16", 16)
TR1.addLocation("FPG19", 19)
TR1.addLocation("FPG25", 25)
TR1.addLocation("FPG28", 28)
TR1.addLocation("FPG33", 33)
TR1.addLocation("FPG37", 37)
TR1.addLocation("FPG42", 42)
TRIA10.addElrefe(TR1)

# ------------------------------------------------------------
QUAD12 = MeshType(nbno=12, dim=2, code="Q12")

Q12 = Elrefe()
Q12.addLocation("NOEU", 12)
Q12.addLocation("NOEU_S", 4)
Q12.addLocation("FPG1", 1)
Q12.addLocation("FPG4", 4)
Q12.addLocation("FPG5", 5)
Q12.addLocation("FPG9", 9)
Q12.addLocation("FPG9COQ", 9)
Q12.addLocation("FPG16", 16)
Q12.addLocation("FPG25", 25)
Q12.addLocation("FPG36", 36)
Q12.addLocation("FPG49", 49)
Q12.addLocation("FPG64", 64)
QUAD12.addElrefe(Q12)

# ------------------------------------------------------------
TETRA20 = MeshType(nbno=20, dim=3, code="T20")

T20 = Elrefe()
T20.addLocation("NOEU", 20)
T20.addLocation("NOEU_S", 4)
T20.addLocation("FPG1", 1)
T20.addLocation("FPG4", 4)
T20.addLocation("FPG4NOS", 8)
T20.addLocation("FPG5", 5)
T20.addLocation("FPG11", 11)
T20.addLocation("FPG15", 15)
T20.addLocation("FPG24", 24)
T20.addLocation("FPG35", 35)
T20.addLocation("FPG46", 46)
T20.addLocation("FPG59", 59)
T20.addLocation("FPG74", 74)
T20.addLocation("FPG94", 94)
T20.addLocation("FPG117", 117)
T20.addLocation("FPG144", 144)
T20.addLocation("FPG204", 204)
TETRA20.addElrefe(T20)

# ------------------------------------------------------------
SEG44 = MeshType(nbno=8, dim=1, code="S44")

SEG43 = MeshType(nbno=7, dim=1, code="S43")

# For Nitsche 3D

TE4TR7 = MeshType(nbno=11, dim=2, code="T11")

T10TR7 = MeshType(nbno=17, dim=2, code="T17")

T20TR7 = MeshType(nbno=27, dim=2, code="T27")

HE8TR7 = MeshType(nbno=15, dim=2, code="H15")

H20TR7 = MeshType(nbno=27, dim=2, code="H26")

H27TR7 = MeshType(nbno=34, dim=2, code="H34")


# For Nitsche 2D

# ------------------------------------------------------------
Q12SE3 = MeshType(nbno=15, dim=1, code="QS7")

# ------------------------------------------------------------
T10SE3 = MeshType(nbno=13, dim=1, code="TS7")


# store all MeshType objects
ELREFS = objects_from_context(globals(), Elrefe)
MESHTYPES = objects_from_context(globals(), MeshType)
