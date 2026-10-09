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

subroutine cmbqbq(main, maout, degree, info)
    use crea_maillage_module
!
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cargeo.h"
#include "asterfort/dismoi.h"
#include "asterfort/infoma.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8), intent(in) :: degree, info
    character(len=8), intent(in) :: main, maout
! ----------------------------------------------------------------------
!         TRANSFORMATION DES MAILLES
!-----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! IN        MAIN   K8  NOM DU MAILLAGE INITIAL
! IN/JXOUT  MAOUT  K8  NOM DU MAILLAGE TRANSFORME
! IN        NBMA    I  NOMBRE DE MAILLES A TRAITER
! IN        LIMA    I  NUMERO DES MAILLES A TRAITER
! ----------------------------------------------------------------------
!
!
    type(Mmesh) :: mesh_conv
    character(len=8) :: conv_type(2)
    integer(kind=8) :: nbma, ima, nbno
    integer(kind=8), pointer :: listCells(:) => null()
!
    call jemarq()
!
! - Create new mesh
!
    call mesh_conv%init(main, convert_max=logical(degree .ne. 1, kind=1))
!
! - Add conversion
!
    if (degree == 1) then
        conv_type = ["SEG3", "SEG2"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["SEG4", "SEG2"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA6", "TRIA3"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA7", "TRIA3"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA10", "TRIA3 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD8", "QUAD4"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD9", "QUAD4"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD12", "QUAD4 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TETRA10", "TETRA4 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TETRA20", "TETRA4 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["PENTA15", "PENTA6 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["PENTA18", "PENTA6 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["PYRAM13", "PYRAM5 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["HEXA20", "HEXA8 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["HEXA27", "HEXA8 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    elseif (degree == 2) then
        conv_type = ["SEG2", "SEG3"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["SEG4", "SEG3"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA3", "TRIA6"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA7", "TRIA6"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA10", "TRIA6 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD4", "QUAD8"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD9", "QUAD8"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD12", "QUAD8 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TETRA4 ", "TETRA10"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TETRA20", "TETRA10"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["PENTA6 ", "PENTA15"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["PENTA18", "PENTA15"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["PYRAM5 ", "PYRAM13"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["HEXA8 ", "HEXA20"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["HEXA27", "HEXA20"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    elseif (degree == 3) then
        conv_type = ["SEG2", "SEG3"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["SEG4", "SEG3"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA3", "TRIA6"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA6", "TRIA6"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA10", "TRIA6 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD4", "QUAD9"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD8", "QUAD9"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD12", "QUAD9 "]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TETRA4 ", "TETRA10"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TETRA20", "TETRA10"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["PENTA6 ", "PENTA18"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["PENTA15", "PENTA18"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["PYRAM5 ", "PYRAM13"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["HEXA8 ", "HEXA27"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["HEXA20", "HEXA27"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    elseif (degree == 4) then
        conv_type = ["SEG2", "SEG4"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["SEG3", "SEG4"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA3 ", "TRIA10"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA6 ", "TRIA10"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TRIA7 ", "TRIA10"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD4 ", "QUAD12"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD8 ", "QUAD12"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["QUAD9 ", "QUAD12"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TETRA4 ", "TETRA20"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
        conv_type = ["TETRA10", "TETRA20"]
        call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Convert cells
!
    call dismoi('NB_MA_MAILLA', main, 'MAILLAGE', repi=nbma)
    call dismoi('NB_NO_MAILLA', main, 'MAILLAGE', repi=nbno)
    call wkvect("&&CREAMA.MA", 'V V I', nbma, vi=listCells)
    do ima = 1, nbma
        listCells(ima) = ima
    end do
!
    call mesh_conv%convert_cells(nbma, listCells, ASTER_TRUE)
!
! - Check conformity
!
    call mesh_conv%check_conformity("A")
!
! - Copy mesh
!
    call mesh_conv%copy_mesh(maout)
!
! - Cleaning
!
    call mesh_conv%clean()
    call jedetr("&&CREAMA.MA")
!
! - Update parameters for modified mesh (bounding box and dimensions)
    call cargeo(maout)

! - Verbose
    call infoma(maout, info)
!
    call jedema()
!
end subroutine
