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


DDL_MECA = LocatedComponents(
    phys=PHY.DEPL_R,
    type="ELNO",
    diff=True,
    components=(
        # Shell nodes
        ("EN1", ("DX", "DY", "DZ", "DRX", "DRY", "DRZ")),
        # Volume nodes
        ("EN2", ("DX", "DY", "DZ")),
    ),
)

CCACOQU = LocatedComponents(
    phys=PHY.CACOQU_R, type="ELEM", components=("EP", "ALPHA", "BETA", "CTOR", "EXCENT", "INERTIE")
)
NGEOMER = LocatedComponents(phys=PHY.GEOM_R, type="ELNO", components=("X", "Y", "Z"))
MVECTUR = ArrayOfComponents(phys=PHY.VDEP_R, locatedComponents=DDL_MECA)
MMATUUR = ArrayOfComponents(phys=PHY.MDEP_R, locatedComponents=DDL_MECA)
MMATUNS = ArrayOfComponents(phys=PHY.MDNS_R, locatedComponents=DDL_MECA)


class RACS2T3(Element):
    """
    THE RACS2T3 CLASS ELEMENT : SEG2/TRIA3 ( EDGE / FACE )
    """

    meshType = MT.SE2TR3
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4, 5)))
    calculs = (
        OP.LIAI_CO_3D(
            te=231,
            para_in=((SP.PGEOMER, NGEOMER), (SP.PCACOQU, CCACOQU)),
            para_out=((SP.PMATUNS, MMATUNS),),
        ),
        OP.TOU_INI_ELEM(
            te=99,
            para_out=((OP.TOU_INI_ELEM.PGEOM_R, LC.CGEOM3D), (OP.TOU_INI_ELEM.PTEMP_R, LC.CTEMPER)),
        ),
        OP.TOU_INI_ELNO(te=99, para_out=((OP.TOU_INI_ELNO.PGEOM_R, NGEOMER),)),
    )


# ELEMENTS DE RACCORD COQUE/3D


class RACS2Q4(RACS2T3):
    """
    THE RACSQ4 CLASS ELEMENT : SEG2/QUA4 (2D EDGE / 3D FACE )
    """

    meshType = MT.SE2QU4
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4, 5, 6)))


class RACS2T6(RACS2T3):
    """
    THE RACSQ4 CLASS ELEMENT : SEG2/QUA4 (2D EDGE / 3D FACE )
    """

    meshType = MT.SE2TR6
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4, 5, 6, 7, 8)))


class RACS2Q8(RACS2T3):
    """
    THE RACSQ4 CLASS ELEMENT : SEG2/QUA4 (2D EDGE / 3D FACE )
    """

    meshType = MT.SE2QU8
    nodes = (SetOfNodes("EN1", (1, 2)), SetOfNodes("EN2", (3, 4, 5, 6, 7, 8, 9, 10)))


class RACS3T3(RACS2T3):
    """
    THE RACSQ4 CLASS ELEMENT : SEG3/TRIA6 (2D EDGE / 3D FACE )
    """

    meshType = MT.SE3TR3
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6)))


class RACS3T6(RACS2T3):
    """
    THE RACSQ4 CLASS ELEMENT : SEG3/TRIA6 (2D EDGE / 3D FACE )
    """

    meshType = MT.SE3TR6
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6, 7, 8, 9)))


class RACS3Q4(RACS2T3):
    """
    THE RACSQ4 CLASS ELEMENT : SEG3/QUA8 (2D EDGE / 3D FACE )
    """

    meshType = MT.SE3QU4
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6, 7)))


class RACS3Q8(RACS2T3):
    """
    THE RACSQ4 CLASS ELEMENT : SEG3/QUA8 (2D EDGE / 3D FACE )
    """

    meshType = MT.SE3QU8
    nodes = (SetOfNodes("EN1", (1, 2, 3)), SetOfNodes("EN2", (4, 5, 6, 7, 8, 9, 10, 11)))
