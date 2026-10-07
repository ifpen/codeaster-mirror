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

# ELEMENTARY TREATMENT OF 2D COUPLING FEM-FEM
# AUGMENTED LAGRANGIAN FORMULATION

# Convention:
# EN1: Slave nodes - FEM
# EN2: Master nodes - FEM
# ----------------
# Modes locaux :
# ----------------


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        # Slave nodes - FEM
        ("EN1", ("DX", "DY", "DZ", "LAGX", "LAGY", "LAGZ")),
        # Master nodes - FEM
        ("EN2", ("DX", "DY", "DZ")),
    ),
)


MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)

MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)


# ------------------------------------------------------------
class CL_Q4Q4(Element):
    """
    THE CMQ4Q4 CLASS ELEMENT : SEG2/SEG2 (3D Edge / 3D edge )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD4 SLAVE  ELEMENT : 1-2-3-4 (DX,DY,DZ)
        QUAD4 MASTER ELEMENT : 5-6-7-8 (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QUAD44
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7, 8)))
    calculs = (
        OP.CHAR_MECA_CPL(
            te=523,
            para_in=(
                (SP.PDEPLMR, DDL_MECA),
                (SP.PDEPLPR, DDL_MECA),
                (SP.PGEOMER, LC.EGEOM3D),
                (SP.PPAIRR, LC.CPAIRR),
            ),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_CPL(
            te=523,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PPAIRR, LC.CPAIRR)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.RIGI_ELAS_CPL(
            te=523,
            para_in=((SP.PGEOMER, LC.EGEOM3D), (SP.PPAIRR, LC.CPAIRR)),
            para_out=((SP.PMATUUR, MMATUUR),),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
                (SP.PNEUT_F, LC.CNTINIF),
                (SP.PNEUT_R, LC.CNTINIR),
                (SP.PNEUT_I, LC.CNTINII),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM3D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENTINIF),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENTINIR),
                (OP.TOU_INI_ELNO.PNEUT_I, LC.ENTINII),
            ),
        ),
    )


# ------------------------------------------------------------
class CL_T3T3(CL_Q4Q4):
    """
    THE COT3T3 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA3/TRIA3 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        TRIA3 MASTER  ELEMENT : 4-5-6 (DX,DY,DZ)
        TRIA3 SLAVE ELEMENT : 1-2-3 (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.TRIA33
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6)))


# ------------------------------------------------------------
class CL_Q4T3(CL_Q4Q4):
    """
    THE COQ4T3 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD4/TRIA3 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD4 SLAVE  ELEMENT : 1-2-3-4 (DX,DY,DZ)
        TRIA3 MASTER ELEMENT : 5-6-7   (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU4TR3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7)))


# ------------------------------------------------------------
class CL_T3Q4(CL_Q4Q4):
    """
    THE COT3Q4 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA3/QUAD4 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        TRIA3 SLAVE  ELEMENT : 1-2-3     (DX,DY,DZ)
        QUAD4 MASTER ELEMENT : 4-5-6-7   (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.TR3QU4
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6, 7)))


# ------------------------------------------------------------
class CL_T6T3(CL_Q4Q4):
    """
    THE COT6T3 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA6/TRIA3 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        TRIA6 SLAVE  ELEMENT : 1-2-3-4-5-6     (DX,DY,DZ)
        TRIA3 MASTER ELEMENT : 7-8-9           (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.TR6TR3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6))), SetOfNodes("EN2", (7, 8, 9))


# ------------------------------------------------------------
class CL_T3T6(CL_Q4Q4):
    """
    THE COT3T6 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA3/TRIA6 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        TRIA3 SLAVE  ELEMENT : 1-2-3          (DX,DY,DZ)
        TRIA6 MASTER ELEMENT : 4-5-6-7-8-9    (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.TR3TR6
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6, 7, 8, 9)))


# ------------------------------------------------------------
class CL_T6Q4(CL_Q4Q4):
    """
    THE COT6Q4 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA6/QUAD4 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        TRIA6 SLAVE  ELEMENT : 1-2-3-4-5-6       (DX,DY,DZ)
        QUAD4 MASTER ELEMENT : 7-8-9-10          (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.TR6QU4
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)), SetOfNodes("EN2", (7, 8, 9, 10)))


# ------------------------------------------------------------
class CL_Q4T6(CL_Q4Q4):
    """
    THE COQ4T6 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD4/TRIA6 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD4 SLAVE  ELEMENT : 1-2-3-4       (DX,DY,DZ)
        TRIA6 MASTER ELEMENT : 5-6-7-8-9-10  (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU4TR6
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7, 8, 9, 10)))


# ------------------------------------------------------------
class CL_T6Q8(CL_Q4Q4):
    """
    THE COT6Q8 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA6/QUAD8 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        TRIA6 SLAVE  ELEMENT : 1-2-3-4-5-6           (DX,DY,DZ)
        QUAD8 MASTER ELEMENT : 7-8-9-10-11-12-13-14  (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.TR6QU8
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12, 13, 14)),
    )


# ------------------------------------------------------------
class CL_Q8T6(CL_Q4Q4):
    """
    THE COQ8T6 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD8/TRIA6 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD8 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8   (DX,DY,DZ)
        TRIA6 MASTER ELEMENT : 9-10-11-12-13-14  (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU8TR6
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14)),
    )


# ------------------------------------------------------------
class CL_T6Q9(CL_Q4Q4):
    """
    THE COT6Q9 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA6/QUAD9 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        TRIA6 SLAVE  ELEMENT : 1-2-3-4-5-6              (DX,DY,DZ)
        QUAD9 MASTER ELEMENT : 7-8-9-10-11-12-13-14-15  (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.TR6QU9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)),
        SetOfNodes("EN2", (7, 8, 9, 10, 11, 12, 13, 14, 15)),
    )


# ------------------------------------------------------------
class CL_Q9T6(CL_Q4Q4):
    """
    THE COQ9T6 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD9/TRIA6 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9  (DX,DY,DZ)
        TRIA6 MASTER ELEMENT : 10-11-12-13-14-15  (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU9TR6
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)),
        SetOfNodes("EN2", (10, 11, 12, 13, 14, 15)),
    )


# ------------------------------------------------------------
class CL_Q8T3(CL_Q4Q4):
    """
    THE COQ8T3 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD8/TRIA3 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD8 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8  (DX,DY,DZ)
        TRIA3 MASTER ELEMENT : 9-10-11          (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU8TR3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)), SetOfNodes("EN2", (9, 10, 11)))


# ------------------------------------------------------------
class CL_T3Q8(CL_Q4Q4):
    """
    THE COT3Q8 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA3/QUAD8 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        TRIA3 SLAVE  ELEMENT : 1-2-3                (DX,DY,DZ)
        QUAD8 MASTER ELEMENT : 4-5-6-7-8-9-10-11    (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.TR3QU8
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6, 7, 8, 9, 10, 11)))


# ------------------------------------------------------------
class CL_Q8Q4(CL_Q4Q4):
    """
    THE COQ8Q4 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD8/QUAD4 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD8 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8   (DX,DY,DZ)
        QUAD4 MASTER ELEMENT : 9-10-11-12        (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU8QU4
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)), SetOfNodes("EN2", (9, 10, 11, 12)))


# ------------------------------------------------------------
class CL_Q4Q8(CL_Q4Q4):
    """
    THE COQ4Q8 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD4/QUAD8 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD4 SLAVE  ELEMENT : 1-2-3-4              (DX,DY,DZ)
        QUAD8 MASTER ELEMENT : 5-6-7-8-9-10-11-12   (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU4QU8
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7, 8, 9, 10, 11, 12)))


# ------------------------------------------------------------
class CL_Q8Q9(CL_Q4Q4):
    """
    THE COQ8Q9 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD8/QUAD9 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD8 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8              (DX,DY,DZ)
        QUAD9 MASTER ELEMENT : 9-10-11-12-13-14-15-16-17    (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU8QU9
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16, 17)),
    )


# ------------------------------------------------------------
class CL_Q9Q8(CL_Q4Q4):
    """
    THE COQ9Q8 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD9/QUAD8 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9              (DX,DY,DZ)
        QUAD8 MASTER ELEMENT : 10-11-12-13-14-15-16-17        (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU9QU8
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)),
        SetOfNodes("EN2", (10, 11, 12, 13, 14, 15, 16, 17)),
    )


# ------------------------------------------------------------
class CL_Q9Q4(CL_Q4Q4):
    """
    THE COQ9Q4 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD9/QUAD4 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9      (DX,DY,DZ)
        QUAD4 MASTER ELEMENT : 10-11-12-13            (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU9QU4
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)), SetOfNodes("EN2", (10, 11, 12, 13)))


# ------------------------------------------------------------
class CL_Q4Q9(CL_Q4Q4):
    """
    THE COQ4Q9 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD4/QUAD9 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD4 SLAVE  ELEMENT : 1-2-3-4                    (DX,DY,DZ)
        QUAD9 MASTER ELEMENT : 5-6-7-8-9-10-11-12-13      (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU4QU9
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4)), SetOfNodes("EN2", (5, 6, 7, 8, 9, 10, 11, 12, 13)))


# ------------------------------------------------------------
class CL_Q9T3(CL_Q4Q4):
    """
    THE COQ9T3 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD9/TRIA3 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9   (DX,DY,DZ)
        TRIA3 MASTER ELEMENT : 10-11-12         (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QU9TR3
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)), SetOfNodes("EN2", (10, 11, 12)))


# ------------------------------------------------------------
class CL_T3Q9(CL_Q4Q4):
    """
    THE COT3Q9 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA3/QUAD9 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        TRIA3 SLAVE  ELEMENT : 1-2-3                     (DX,DY,DZ)
        QUAD9 MASTER ELEMENT : 4-5-6-7-8-9-10-11-12      (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.TR3QU9
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6, 7, 8, 9, 10, 11, 12)))


# ------------------------------------------------------------
class CL_Q8Q8(CL_Q4Q4):
    """
    THE COQ8Q8 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : TRIA3/QUAD9 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD8 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8                (DX,DY,DZ)
        QUAD8 MASTER ELEMENT : 9-10-11-12-13-14-14-15-16      (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QUAD88
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8)),
        SetOfNodes("EN2", (9, 10, 11, 12, 13, 14, 15, 16)),
    )


# ------------------------------------------------------------
class CL_Q9Q9(CL_Q4Q4):
    """
    THE COQ9Q9 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD9/QUAD9 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9                  (DX,DY,DZ)
        QUAD9 MASTER ELEMENT : 10-11-12-13-14-14-15-16-17-18      (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.QUAD99
    nodes = (
        SetOfNodes("EN1", (1, 2, 3, 4, 5, 6, 7, 8, 9)),
        SetOfNodes("EN2", (10, 11, 12, 13, 14, 15, 16, 17, 18)),
    )


# ------------------------------------------------------------
class CL_T6T6(CL_Q4Q4):
    """
    THE COQ9Q9 DERIVED FROM  CMQ4Q4 CLASS ELEMENT  : QUAD9/QUAD9 (3D Face / 3D Face )
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 3D : elementary treatments
    Local Numerotation :
        QUAD9 SLAVE  ELEMENT : 1-2-3-4-5-6-7-8-9                  (DX,DY,DZ)
        QUAD9 MASTER ELEMENT : 10-11-12-13-14-14-15-16-17-18      (DX,DY,DZ)
    Input parameters :
        PACCE_M - ACCELERATION at T-
        PVITE_M - VELOCITY at T-
        PDEPL_M - DISPL. at T-
        PVITE_P - VELOCITY at T+
        PDEPL_P - DISPL. at T+
        PGEOMER - CURRENT GEOMETRY
        PCONFR - FRICTIONAL CONTACT PARAMETERS
    Output parameters :
        PMATUNS : NON SYMMETRIC MATRIX (te=352)
        PMMATUR : SYMMETRIC MATRIX (te=352)
        PMMATUR : VECTOR OF CONTACT LOAD (te=352)
    """

    meshType = MT.TRIA66
    nodes = (SetOfNodes("EN1", (1, 2, 3, 4, 5, 6)), SetOfNodes("EN2", (7, 8, 9, 10, 11, 12)))


# ------------------------------------------------------------
class CL_POI3D(CL_Q4Q4):
    """
    CL_S3S3 DERIVED FROM THE CL_S2S2 CLASS ELEMENT : SEG2/SEG2
    LIAISON_ELEM / LAGRANGIAN / SEGMENT-TO-SEGMENT
        Coupling FEM/FEM Element in 2D : elementary treatments
    """

    meshType = MT.POI1
    nodes = (SetOfNodes("EN1", (1,)), SetOfNodes("EN2", ()))
    calculs = (
        OP.CHAR_MECA_CPL(
            te=524,
            para_in=((SP.PDEPLMR, DDL_MECA), (SP.PDEPLPR, DDL_MECA), (SP.PGEOMER, LC.EGEOM3D)),
            para_out=((SP.PVECTUR, MVECTUR),),
        ),
        OP.RIGI_CPL(te=524, para_in=((SP.PGEOMER, LC.EGEOM3D),), para_out=((SP.PMATUUR, MMATUUR),)),
        OP.RIGI_ELAS_CPL(
            te=524, para_in=((SP.PGEOMER, LC.EGEOM3D),), para_out=((SP.PMATUUR, MMATUUR),)
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=(
                (OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D),
                (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER),
                (SP.PNEUT_F, LC.CNTINIF),
                (SP.PNEUT_R, LC.CNTINIR),
                (SP.PNEUT_I, LC.CNTINII),
            ),
        ),
        OP.TOU_INI_ELNO(
            te=99,
            para_out=(
                (OP.TOU_INI_ELNO.PGEOM_R, LC.EGEOM3D),
                (OP.TOU_INI_ELNO.PNEUT_F, LC.ENTINIF),
                (OP.TOU_INI_ELNO.PNEUT_R, LC.ENTINIR),
                (OP.TOU_INI_ELNO.PNEUT_I, LC.ENTINII),
            ),
        ),
    )
