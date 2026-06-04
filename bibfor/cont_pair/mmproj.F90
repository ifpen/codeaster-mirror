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

subroutine mmproj(cellCode, cellNbNode, cellDime, cellCoor, poinCoor, &
                  projtIterMaxi, projToleMaxi, toleou, dirapp, dir, &
                  ksi1, ksi2, tau1, tau2, iproj, &
                  projError)
!
!
    implicit none
#include "asterf_types.h"
#include "asterfort/mmnewd.h"
#include "asterfort/projOrthoNewton.h"
#include "asterfort/mmtole.h"
    character(len=8) :: cellCode
    integer(kind=8) :: cellDime
    integer(kind=8) :: cellNbNode
    real(kind=8) :: cellCoor(27)
    real(kind=8) :: poinCoor(3)
    aster_logical :: dirapp
    real(kind=8) :: dir(3)
    real(kind=8) :: ksi1, ksi2
    real(kind=8) :: tau1(3), tau2(3)
    real(kind=8) :: toleou
    integer(kind=8) :: iproj
    integer(kind=8) :: projError
    integer(kind=8) :: projtIterMaxi
    real(kind=8) :: projToleMaxi
!
! ----------------------------------------------------------------------
!
! ROUTINE APPARIEMENT
!
! FAIRE LA PROJECTION D'UN POINT SUR UNE MAILLE DONNEE
!
! ----------------------------------------------------------------------
!
!
! IN  ALIAS  : TYPE DE MAILLE
! IN  NNO    : NOMBRE DE NOEUD SUR LA MAILLE
! IN  NDIM   : DIMENSION DE LA MAILLE (2 OU 3)
! IN  COORMA : COORDONNEES DES NOEUDS DE LA MAILLE
! IN  COORPT : COORDONNEES DU NOEUD A PROJETER SUR LA MAILLE
! IN  ITEMAX : NOMBRE MAXI D'ITERATIONS DE NEWTON POUR LA PROJECTION
! IN  EPSMAX : RESIDU POUR CONVERGENCE DE NEWTON POUR LA PROJECTION
! IN  TOLEOU : TOLERANCE POUR LE PROJETE HORS MAILLE
! IN  DIRAPP : VAUT .TRUE. SI APPARIEMENT DANS UNE DIRECTION DE
!              RECHERCHE DONNEE (PAR DIR)
! IN  DIR    : DIRECTION D'APPARIEMENT
! OUT KSI1   : PREMIERE COORDONNEE PARAMETRIQUE DU POINT PROJETE
! OUT KSI2   : SECONDE COORDONNEE PARAMETRIQUE DU POINT PROJETE
! OUT TAU1   : PREMIER VECTEUR TANGENT EN (KSI1,KSI2)
! OUT TAU2   : SECOND VECTEUR TANGENT EN (KSI1,KSI2)
! OUT IPROJ  : VAUT 0 SI POINT PROJETE DANS L'ELEMENT
!                   1 SI POINT PROJETE DANS LA ZONE DEFINIE PAR TOLEOU
!                   2 SI POINT PROJETE EN DEHORS (EXCLUS)
! OUT NIVERR : RETOURNE UN CODE ERREUR
!                0  TOUT VA BIEN
!                1  ECHEC DE NEWTON
!
! ----------------------------------------------------------------------
!
    projError = 0
!
! --- ALGO DE NEWTON POUR LA PROJECTION SUIVANT UNE DIRECTION DONNEE
!
    if (dirapp) then
        call mmnewd(cellCode, cellNbNode, cellDime, cellCoor, poinCoor, &
                    projtIterMaxi, projToleMaxi, dir, ksi1, ksi2, &
                    tau1, tau2, projError)
    else
        call projOrthoNewton(cellCode, cellNbNode, cellDime, cellCoor, poinCoor, &
                             projtIterMaxi, projToleMaxi, &
                             ksi1, ksi2, &
                             tau1, tau2, &
                             projError)
    end if
    if (projError .gt. 0) then
        goto 999
    end if
!
! --- AJUSTEMENT PROJECTION HORS ZONE
!
    call mmtole(cellCode, cellNbNode, cellDime, cellCoor, toleou, &
                ksi1, ksi2, tau1, tau2, iproj)
!
999 continue
!
end subroutine
