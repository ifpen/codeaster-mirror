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

! Total number of cell types
#define MT_NTYMAX  94
! Total number of real (physical) cell types
#define MT_NPHMAX  28
! Maximum number of nodes of all cells
#define MT_NNOMAX   27
#define MT_NNOMAX2D 12
#define MT_NNOMAX3D 27
! Maximum number of families for integration schemes
#define MT_NBFAMX  20
! Maximum number of Gauss points
#define MT_NBPGMX  1000
! Index of TYPE_ELEM (geometric)
#define MT_POI1    1
#define MT_SEG2    2
#define MT_SEG3    4
#define MT_SEG4    6
#define MT_TRIA3   7
#define MT_TRIA6   9
#define MT_TRIA7   11
#define MT_TRIA10  87
#define MT_QUAD4   12
#define MT_QUAD8   14
#define MT_QUAD9   16
#define MT_QUAD12  88
#define MT_TETRA4  18
#define MT_TETRA10 19
#define MT_TETRA15 28
#define MT_TETRA20 89
#define MT_PENTA6  20
#define MT_PENTA15 21
#define MT_PENTA18 22
#define MT_PENTA21 29
#define MT_PYRAM5  23
#define MT_PYRAM13 24
#define MT_PYRAM19 30
#define MT_HEXA8   25
#define MT_HEXA20  26
#define MT_HEXA27  27
#define MT_HEXA9   70
#define MT_PENTA7  71
! Maximum number of nodes on linear cells
#define MT_MAX_NBNODE_LINE  8
