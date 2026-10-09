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
subroutine te0398(option, nomte)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "MeshTypes_type.h"
!
    character(len=16) :: nomte, option
! ----------------------------------------------------------------------
! FONCTION REALISEE:  CALCUL DU GRADIENT AUX NOEUDS D'UN CHAMP SCALAIRE
!                      AUX NOEUDS A 9 COMPOSANTES
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!.......................................................................
!
!
    real(kind=8) :: dfdx(MT_NNOMAX3D), dfdy(MT_NNOMAX3D), dfdz(MT_NNOMAX3D)
    real(kind=8) :: jac, gradx, grady, gradz
!
    integer(kind=8) :: ndim, nno, nnos, npg
    integer(kind=8) :: ipoids, ivf, idfde, jgano, igeom, ineut, igr
    integer(kind=8) :: i, kp, ino
!
!
!
! DEB ------------------------------------------------------------------
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PNEUTER', 'L', ineut)
    call jevech('PGNEUTR', 'E', igr)
!
!   C'EST SUREMENT MIEUX D'INVERSER LES BOUCLES !!
!
!     BOUCLE SUR LES 9 CHAMPS SCALAIRES D'ENTREE
    do i = 1, 9
!
!       BOUCLE SUR LES POINTS DE GAUSS
        do kp = 1, npg
            call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                        jac, dfdx, dfdy, dfdz)
            gradx = 0.0d0
            grady = 0.0d0
            gradz = 0.0d0
            do ino = 1, nno
                gradx = gradx+zr(ineut-1+9*(ino-1)+i)*dfdx(ino)
                grady = grady+zr(ineut-1+9*(ino-1)+i)*dfdy(ino)
                gradz = gradz+zr(ineut-1+9*(ino-1)+i)*dfdz(ino)
            end do
            zr(igr-1+27*(kp-1)+3*(i-1)+1) = gradx
            zr(igr-1+27*(kp-1)+3*(i-1)+2) = grady
            zr(igr-1+27*(kp-1)+3*(i-1)+3) = gradz
        end do
!
    end do
!
!
! FIN ------------------------------------------------------------------
end subroutine
