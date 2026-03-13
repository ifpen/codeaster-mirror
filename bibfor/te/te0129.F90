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

subroutine te0129(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/foderi.h"
#include "asterfort/jevech.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS RESIDUS
!                          OPTION : 'RESI_THER_FLUXNL'
!                          ELEMENTS DE FACE 3D
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
    real(kind=8) :: nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9), jac
    real(kind=8) :: tpg, tpg_b, alpha, alpha_b, rbid
    integer(kind=8) :: ndim, nno, npg1, ipoids, ivf, idfdx, idfdy
    integer(kind=8) :: igeom, iflux, itempi, btempi, itemps, ino, jno, iveres
    integer(kind=8) :: nnos, jgano
    character(len=8) :: coef
    real(kind=8) :: theta, delta_t
    aster_logical :: l_stat
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idec, j, jdec, kdec, kp, ldec
!
!-----------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx+1
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PTEMPEI', 'L', itempi)
    call jevech('PTEMPER', 'L', btempi)
    call jevech('PFLUXNL', 'L', iflux)
    call jevech('PRESIDU', 'E', iveres)
    coef = zk8(iflux)
    if (coef(1:7) .eq. '&FOZERO') goto 999

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
!    CALCUL DES PRODUITS VECTORIELS OMI   OMJ
!
    do ino = 1, nno
        i = igeom+3*(ino-1)-1
        do jno = 1, nno
            j = igeom+3*(jno-1)-1
            sx(ino, jno) = zr(i+2)*zr(j+3)-zr(i+3)*zr(j+2)
            sy(ino, jno) = zr(i+3)*zr(j+1)-zr(i+1)*zr(j+3)
            sz(ino, jno) = zr(i+1)*zr(j+2)-zr(i+2)*zr(j+1)
        end do
    end do
!
    do kp = 1, npg1
        kdec = (kp-1)*nno*ndim
        ldec = (kp-1)*nno
        nx = 0.0d0
        ny = 0.0d0
        nz = 0.0d0
!
!   CALCUL DE LA NORMALE AU POINT DE GAUSS KP
!
        do i = 1, nno
            idec = (i-1)*ndim
            do j = 1, nno
                jdec = (j-1)*ndim
                nx = nx+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sx(i, j)
                ny = ny+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sy(i, j)
                nz = nz+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sz(i, j)
            end do
        end do
!
!   CALCUL DU JACOBIEN AU POINT DE GAUSS KP
!
        jac = sqrt(nx*nx+ny*ny+nz*nz)
!
        tpg = 0.d0
        tpg_b = 0.d0
        do i = 1, nno
            tpg = tpg+zr(itempi+i-1)*zr(ivf+ldec+i-1)
            if (l_stat) tpg_b = tpg_b+zr(btempi+i-1)*zr(ivf+ldec+i-1)
        end do
        call foderi(coef, tpg, alpha, rbid)
        call foderi(coef, tpg_b, alpha_b, rbid)
!
! ----- ON RAJOUTE DANS LE RESIDU LE TERME (1-THETA)*ALPHAP QUI NE
! ----- FIGURE PAS DANS LE 2ND MEMBRE LINEAIRE
!
        do i = 1, nno
            zr(iveres+i-1) = zr(iveres+i-1)-zr(ipoids+kp-1)*jac*alpha*zr(ivf+ldec+i-1)
            if (l_stat) then
                zr(iveres+i-1) = zr(iveres+i-1) &
                                 +zr(ipoids+kp-1)*jac*alpha_b*zr(ivf+ldec+i-1)
            end if
        end do
    end do
999 continue
! FIN ------------------------------------------------------------------
end subroutine
