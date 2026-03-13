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
subroutine te0138(option, nomte)
!    - FONCTION REALISEE:  CALCUL DES VECTEURS RESIDUS
!                          OPTION : 'RESI_THER_FLUXNL'
!                          ELEMENTS DE FACE 2D
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!   -------------------------------------------------------------------
!     ASTER INFORMATIONS:
!       04/04/02 (OB): CORRECTION BUG, OPTION NON PORTEE EN LUMPE
!       + MODIFS FORMELLES: IMPLICIT NONE, LAXI, IDENTATION...
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
! PARAMETRES D'APPEL
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/connec.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/foderi.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/vff2dn.h"
!
    character(len=16) :: option, nomte
!
!
! VARIABLES LOCALES
    character(len=8) :: coef, alias8
    real(kind=8) :: poids, r, nx, ny, alpha, alpha_b, rbid, tpg, tpg_b, coorse(18)
    real(kind=8) :: vectt(9)
    integer(kind=8) :: nno, nnos, jgano, ndim, kp, npg, ipoids, ivf, idfde, igeom, i, j
    integer(kind=8) :: l, li, iflux, iveres, nse, c(6, 9), ise, nnop2, itempi, itemps, btempi
    integer(kind=8) :: ibid
    aster_logical :: laxi
    character(len=8) :: elrefe
    real(kind=8) :: theta, delta_t
    aster_logical :: l_stat
!
!====
! 1.1 PREALABLES: RECUPERATION ADRESSES FONCTIONS DE FORMES...
!====
    call elref1(elrefe)
!
    if (lteatt('LUMPE', 'OUI')) then
        call teattr('S', 'ALIAS8', alias8, ibid)
        if (alias8(6:8) .eq. 'SE3') elrefe = 'SE2'
    end if
!
    call elrefe_info(elrefe=elrefe, fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PTEMPEI', 'L', itempi)
    call jevech('PTEMPER', 'L', btempi)
    call jevech('PFLUXNL', 'L', iflux)
    call jevech('PRESIDU', 'E', iveres)

!   SET STEADY OR TRANSIENT COMPUTATION
    l_stat = .false.
    theta = zr(itemps+2)
    delta_t = zr(itemps+1)
    ! FIXME: find a better way to define l_stat. Ideally it should be theta = -1
    ! See issue 34998
    if ((theta .eq. 1.d0) .and. (delta_t .lt. 0.d0)) then
        l_stat = .true.
    end if
!
! INITS.
    coef = zk8(iflux)
    if (coef(1:7) .eq. '&FOZERO') goto 100
    call connec(nomte, nse, nnop2, c)
    do i = 1, nnop2
        vectt(i) = 0.d0
    end do
!
! --- CALCUL ISO-P2 : BOUCLE SUR LES SOUS-ELEMENTS -------
!
    do ise = 1, nse
!
        do i = 1, nno
            do j = 1, 2
                coorse(2*(i-1)+j) = zr(igeom-1+2*(c(ise, i)-1)+j)
            end do
        end do
        do kp = 1, npg
            call vff2dn(ndim, nno, kp, ipoids, idfde, &
                        coorse, nx, ny, poids)
            if (laxi) then
                r = 0.d0
                do i = 1, nno
                    l = (kp-1)*nno+i
                    r = r+coorse(2*(i-1)+1)*zr(ivf+l-1)
                end do
                poids = poids*r
            end if
!
            tpg = 0.d0
            tpg_b = 0.d0
            do i = 1, nno
                l = (kp-1)*nno+i
                tpg = tpg+zr(itempi-1+c(ise, i))*zr(ivf+l-1)
                if (l_stat) tpg_b = tpg_b+zr(btempi-1+c(ise, i))*zr(ivf+l-1)
            end do
            call foderi(coef, tpg, alpha, rbid)
            call foderi(coef, tpg_b, alpha_b, rbid)
!
            do i = 1, nno
                li = ivf+(kp-1)*nno+i-1
                vectt(c(ise, i)) = vectt(c(ise, i))-poids*alpha*zr(li)
                if (l_stat) then
                    vectt(c(ise, i)) = vectt(c(ise, i)) &
                                       +poids*alpha_b*zr(li)
                end if
            end do
        end do
    end do
!
    do i = 1, nnop2
        zr(iveres-1+i) = vectt(i)
    end do
100 continue
! FIN ------------------------------------------------------------------
end subroutine
