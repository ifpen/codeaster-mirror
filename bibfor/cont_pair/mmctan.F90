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

subroutine mmctan(numema, alias, nno, ndim, coorma, &
                  coorno, epsmax, tau1, tau2, err_appa)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/infdbg.h"
#include "asterfort/reereg.h"
#include "asterfort/mmdonf.h"
#include "asterfort/mmtang.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: numema
    character(len=8) :: alias
    integer(kind=8) :: ndim, nno, err_appa
    real(kind=8) :: epsmax, coorno(3), coorma(27)
    real(kind=8) :: tau1(3), tau2(3)
!
! ----------------------------------------------------------------------
!
! ROUTINE APPARIEMENT (UTILITAIRE)
!
! CALCUL DES TANGENTES EN UN NOEUD D'UNE MAILLE
!
! ----------------------------------------------------------------------
!
!
! IN  NOMMAI : NOM DE LA MAILLE
! IN  ALIAS  : TYPE DE LA MAILLE
! IN  NNO    : NOMBRE DE NOEUDS DE LA MAILLE
! IN  NDIM   : DIMENSION DE LA MAILLE
! IN  COORMA : CORDONNNES DE LA MAILLE
! IN  COORNO : COORODNNEES DU NOEUD
! IN  ITEMAX : NOMBRE MAXI D'ITERATIONS DE NEWTON POUR LA PROJECTION
! IN  EPSMAX : RESIDU POUR CONVERGENCE DE NEWTON POUR LA PROJECTION
! OUT TAU1   : PREMIERE TANGENTE (NON NORMALISEE)
! OUT TAU2   : SECONDE TANGENTE (NON NORMALISEE)
!
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: niverr
    real(kind=8) :: ksi(3)
    real(kind=8) :: dff(2, 9)
!
! ----------------------------------------------------------------------
!
    call infdbg('APPARIEMENT', ifm, niv)
!
! --- INITIALISATIONS
!
    niverr = 0
!
! --- CALCUL DES VECTEURS TANGENTS DE LA MAILLE EN CE NOEUD
!
    call reereg('C', alias, nno, coorma, coorno, ndim, ksi, niverr, epsmax, 3)
!
! --- GESTION DES ERREURS LORS DU NEWTON LOCAL POUR LA PROJECTION
!
    if (niverr .eq. 1) then
        err_appa = 1
        if (niv .ge. 2) then
            call utmess('I', 'APPARIEMENT_13', si=numema, nr=3, valr=coorno)
        end if
    end if
!
    call mmdonf(alias, ksi(1), ksi(2), dff)
!
    tau1 = 0.d0
    tau2 = 0.d0
    call mmtang(ndim, nno, coorma, dff, tau1, tau2)
!
end subroutine
