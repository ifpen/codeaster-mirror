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
subroutine mmdonf(cellCode, ksi1, ksi2, dff)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrfdf.h"
!
    character(len=8), intent(in) :: cellCode
    real(kind=8), intent(in) :: ksi1, ksi2
    real(kind=8), intent(out) :: dff(2, 9)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE CONTACT (TOUTES METHODES - UTILITAIRE)
!
! CALCUL DES DERIVEES PREMIERES DES FONCTIONS DE FORME EN UN POINT
! DE L'ELEMENT DE REFERENCE
!
! --------------------------------------------------------------------------------------------------
!
! IN  ALIAS  : NOM D'ALIAS DE L'ELEMENT
! IN  NNO    : NOMBRE DE NOEUD DE L'ELEMENT
! IN  NDIM   : DIMENSION DE LA MAILLE (2 OU 3)
! IN  KSI1   : POINT DE CONTACT SUIVANT KSI1 DES
!               FONCTIONS DE FORME ET LEURS DERIVEES
! IN  KSI2   : POINT DE CONTACT SUIVANT KSI2 DES
!               FONCTIONS DE FORME ET LEURS DERIVEES
! OUT DFF    : DERIVEES PREMIERES DES FONCTIONS DE FORME EN XI YI
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: ksi(2), d2f(3, 9)
!
! --------------------------------------------------------------------------------------------------
!
    dff = 0.d0
    d2f = 0.d0
    ksi(1) = ksi1
    ksi(2) = ksi2
    call elrfdf(cellCode, ksi, d2f)
    dff(1:2, :) = d2f(1:2, :)
!
end subroutine
