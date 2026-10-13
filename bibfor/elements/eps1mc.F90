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
subroutine eps1mc(nno, ndim, nbsig, npg, ipoids, &
                  ivf, idfde, xyz, depl, nharm, &
                  eps1)
!.======================================================================
    implicit none
!
!      EPS1MC   -- CALCUL DES  DEFORMATIONS AUX POINTS D'INTEGRATION
!                  POUR LES ELEMENTS ISOPARAMETRIQUES
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NNO            IN     I        NOMBRE DE NOEUDS DE L'ELEMENT
!    NDIM           IN     I        DIMENSION DE L'ELEMENT (2 OU 3)
!    NBSIG          IN     I        NOMBRE DE CONTRAINTES ASSOCIE
!                                   A L'ELEMENT
!    NPG            IN     I        NOMBRE DE POINTS D'INTEGRATION
!                                   DE L'ELEMENT
!    IVF            IN     I        POINTEUR FONCTIONS DE FORME
!    IPOIDS         IN     I        POINTEUR POIDS D'INTEGRATION
!    IDFDE          IN     I        PT DERIVEES DES FONCTIONS DE FORME
!    XYZ(1)         IN     R        COORDONNEES DES CONNECTIVITES
!    DEPL(1)        IN     R        VECTEUR DES DEPLACEMENTS SUR
!                                   L'ELEMENT
!    NHARM          IN     R        NUMERO D'HARMONIQUE
!    EPS1(1)        OUT    R        DEFORMATIONS DU PREMIER ORDRE
!                                   AUX POINTS D'INTEGRATION
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "asterfort/bmatmc.h"
#include "FE_module.h"
    real(kind=8) :: xyz(1), depl(1), eps1(1)
    real(kind=8) :: nharm
! -----  VARIABLES LOCALES
    real(kind=8) :: b(6*MAX_BV_CG), jacgau
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     -----------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idfde, igau, ipoids, ivf, j, nbinco
    integer(kind=8) :: nbsig, ndim, nno, npg
    real(kind=8) :: s, undemi, zero
!-----------------------------------------------------------------------
    zero = 0.0d0
    undemi = 0.5d0
    nbinco = ndim*nno
!
    do i = 1, nbsig*npg
        eps1(i) = zero
    end do
!
! --- CALCUL DES DEFORMATIONS AUX POINTS D'INTEGRATION
! ---  BOUCLE SUR LES POINTS D'INTEGRATION
!      -----------------------------------
    do igau = 1, npg
!
!  --      CALCUL DE LA MATRICE B RELIANT LES DEFORMATIONS DU
!  --      PREMIER ORDRE AUX DEPLACEMENTS AU POINT D'INTEGRATION
!  --      COURANT : (EPS_1) = (B)*(UN)
!          ----------------------------
        call bmatmc(igau, nbsig, xyz, ipoids, ivf, &
                    idfde, nno, nharm, jacgau, b)
!
! ---      CALCUL DU VECTEUR DES COMPOSANTES DU TENSEUR DES
! ---      DEFORMATIONS AU POINT D'INTEGRATION COURANT
!          -------------------------------------------
        do i = 1, nbsig
!
            s = zero
!
            do j = 1, nbinco
                s = s+depl(j)*b((j-1)*nbsig+i)
            end do
!
            eps1(nbsig*(igau-1)+i) = s
        end do
!
        do i = 4, nbsig
            eps1(nbsig*(igau-1)+i) = undemi*eps1(nbsig*(igau-1)+i)
        end do
!
    end do
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
