! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine cmhho(mesh_in, mesh_out, nb_list_elem, list_elem)
!
    use crea_maillage_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
!
    integer(kind=8), intent(in) :: nb_list_elem, list_elem(nb_list_elem)
    character(len=8), intent(in) :: mesh_in, mesh_out
!
! ----------------------------------------------------------------------
!         TRANSFORMATION DES MAILLES POUR HHO
! ----------------------------------------------------------------------
! IN        mesh_in   K8  NOM DU MAILLAGE INITIAL
! IN/JXOUT  mesh_out  K8  NOM DU MAILLAGE TRANSFORME
! IN        NBMA    I  NOMBRE DE MAILLES A TRAITER
! IN        LIMA    I  NUMERO ET TYPE DES MAILLES A TRAITER
! ----------------------------------------------------------------------
!
    type(Mmesh) :: mesh_hho
    character(len=8) :: conv_type(2)

!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Create new mesh
!
    call mesh_hho%init(mesh_in)
!
! - Add conversion
!
    conv_type = ["SEG2", "SEG3"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["SEG4", "SEG3"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TRIA3", "TRIA7"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TRIA6", "TRIA7"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TRIA10", "TRIA7 "]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["QUAD4", "QUAD9"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["QUAD8", "QUAD9"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["QUAD12", "QUAD9 "]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TETRA4 ", "TETRA15"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TETRA10", "TETRA15"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TETRA20", "TETRA15"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["HEXA8 ", "HEXA27"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["HEXA20", "HEXA27"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PYRAM5 ", "PYRAM19"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PYRAM13", "PYRAM19"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PENTA6 ", "PENTA21"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PENTA15", "PENTA21"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PENTA18", "PENTA21"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
!
! - Convert cells
!
    call mesh_hho%convert_cells(nb_list_elem, list_elem, ASTER_TRUE)
!
! - Check conformity
!
    call mesh_hho%check_conformity("F")
!
! - Fix cells
!
    call mesh_hho%fix(ASTER_FALSE, ASTER_TRUE, ASTER_TRUE, ASTER_FALSE, ASTER_TRUE, 1.0d-12)
!
! - Copy mesh
!
    call mesh_hho%copy_mesh(mesh_out)
!
! - Cleaning
!
    call mesh_hho%clean()
!
    call jedema()
end subroutine
