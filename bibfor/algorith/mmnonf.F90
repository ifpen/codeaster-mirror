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
subroutine mmnonf(cellCode, ksi1, ksi2, ff)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrfvf.h"
!
    character(len=8), intent(in) :: cellCode
    real(kind=8), intent(in) :: ksi1, ksi2
    real(kind=8), intent(out) :: ff(9)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE CONTACT (TOUTES METHODES - UTILITAIRE)
!
! CALCUL DES FONCTIONS DE FORME EN UN POINT DE L'ELEMENT DE REFERENCE
!
! --------------------------------------------------------------------------------------------------
!
! IN  ALIAS  : NOM D'ALIAS DE L'ELEMENT
! IN  KSI1   : POINT DE CONTACT SUIVANT KSI1 DES
!               FONCTIONS DE FORME ET LEURS DERIVEES
! IN  KSI2   : POINT DE CONTACT SUIVANT KSI2 DES
!               FONCTIONS DE FORME ET LEURS DERIVEES
! OUT FF     : FONCTIONS DE FORMES EN XI,YI
!
! --------------------------------------------------------------------------------------------------
!
    ff = 0.d0
    call elrfvf(cellCode, [ksi1, ksi2], ff)
!
end subroutine
