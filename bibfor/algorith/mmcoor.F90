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
subroutine mmcoor(cellCode, cellNbNode, cellCoor, &
                  ksi1, ksi2, coorpt)
!
    implicit none
!
#include "asterfort/mmnonf.h"
!
    character(len=8), intent(in) :: cellCode
    integer(kind=8), intent(in)::  cellNbNode
    real(kind=8), intent(in) :: ksi1, ksi2
    real(kind=8), intent(in) :: cellCoor(27)
    real(kind=8), intent(out) :: coorpt(3)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE CONTACT (TOUTES METHODES - UTILITAIRE)
!
! CALCUL DES COORDONNEES D'UN POINT SUR UNE MAILLE A PARTIR
! DE SES COORDONNEES PARAMETRIQUES
!
! --------------------------------------------------------------------------------------------------
!
! IN  ALIAS  : TYPE DE MAILLE
! IN  NNO    : NOMBRE DE NOEUD SUR LA MAILLE
! IN  COORMA : COORDONNEES DES NOEUDS DE LA MAILLE
! IN  KSI1   : COORDONNEE PARAMETRIQUE KSI DU PROJETE
! IN  KSI2   : COORDONNEE PARAMETRIQUE ETA DU PROJETE
! OUT COORPT : COORDONNEES DU POINT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iDime, iNode
    real(kind=8) :: ff(9)
!
! --------------------------------------------------------------------------------------------------
!
    coorpt = 0.d0

! - FONCTIONS DE FORME
    call mmnonf(cellCode, ksi1, ksi2, ff)

! - COORDONNEES DU POINT
    do iDime = 1, 3
        do iNode = 1, cellNbNode
            coorpt(iDime) = ff(iNode)*cellCoor(3*(iNode-1)+iDime)+coorpt(iDime)
        end do
    end do
!
end subroutine
