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
! aslint: disable=W1306
!
subroutine nirmtd(ndim, nno1, nno2, nno3, npg, &
                  iw, vff2, vff3, ivf1, idff1, &
                  vu, vg, vp, jvGeom, materPara, &
                  matr)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/bmatmc.h"
#include "asterfort/dmatmc.h"
#include "asterfort/nbsigm.h"
#include "blas/dscal.h"
#include "jeveux.h"
#include "FE_module.h"
#include "MeshTypes_type.h"
!
    integer(kind=8) :: ndim, nno1, nno2, nno3, npg, iw, idff1
    integer(kind=8) :: vu(3, MT_NNOMAX), vg(MT_NNOMAX), vp(MT_NNOMAX)
    integer(kind=8) :: ivf1, jvGeom
    real(kind=8) :: vff2(nno2, npg), vff3(nno3, npg)
    real(kind=8) :: matr(*)
    type(Material_Para), intent(inout) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
! CALCUL DE LA RIGIDITE MECANIQUE POUR LES ELEMENTS INCOMPRESSIBLES POUR LES GRANDES DEFORMATIONS
! 3D/D_PLAN/AXIS
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO1    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AUX DEPLACEMENTS
! IN  NNO2    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AU GONFLEMENT
! IN  NNO3    : NOMBRE DE NOEUDS DE L'ELEMENT LIES A LA PRESSION
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : POIDS DES POINTS DE GAUSS
! IN  VFF2    : VALEUR  DES FONCTIONS DE FORME LIES AU GONFLEMENT
! IN  VFF3    : VALEUR  DES FONCTIONS DE FORME LIES A LA PRESSION
! IN  IDFF1   : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  VU      : TABLEAU DES INDICES DES DDL DE DEPLACEMENTS
! IN  VG      : TABLEAU DES INDICES DES DDL DE GONFLEMENT
! IN  VP      : TABLEAU DES INDICES DES DDL DE PRESSION
! IN  IGEOM   : POINTEUR SUR LES COORDONEES DES NOEUDS
! IN  MATE    : MATERIAU CODE
! OUT MATR    : MATRICE DE RIGIDITE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer(kind=8) :: kpg
    integer(kind=8) :: ia, na, ra, sa, ib, nb, rb, sb, ja, jb
    integer(kind=8) :: os, kk
    integer(kind=8) :: vuiana, vgra, vpsa
    integer(kind=8) :: nbsig
    real(kind=8) :: w
    real(kind=8) :: dsidep(2*ndim, 2*ndim)
    real(kind=8) :: b(2*ndim, MAX_BV_CG), def(2*ndim, nno1, ndim), deftr(nno1, ndim)
    real(kind=8) :: ddev(2*ndim, 2*ndim), devd(2*ndim, 2*ndim)
    real(kind=8) :: dddev(2*ndim, 2*ndim)
    real(kind=8) :: iddid, devdi(2*ndim), iddev(2*ndim)
    real(kind=8) :: t1, notime
    real(kind=8) :: idev(6, 6), idev2(4, 4), kr(6), kd(6)
    blas_int :: b_incx, b_n
!
    data kr/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
    data idev2/2.d0, -1.d0, -1.d0, 0.d0,&
     &                 -1.d0, 2.d0, -1.d0, 0.d0,&
     &                 -1.d0, -1.d0, 2.d0, 0.d0,&
     &                  0.d0, 0.d0, 0.d0, 3.d0/
    data idev/2.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0,&
     &                 -1.d0, 2.d0, -1.d0, 0.d0, 0.d0, 0.d0,&
     &                 -1.d0, -1.d0, 2.d0, 0.d0, 0.d0, 0.d0,&
     &                  0.d0, 0.d0, 0.d0, 3.d0, 0.d0, 0.d0,&
     &                  0.d0, 0.d0, 0.d0, 0.d0, 3.d0, 0.d0,&
     &                  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 3.d0/
!
! --------------------------------------------------------------------------------------------------
!
    notime = r8vide()

! - NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
    nbsig = nbsigm()
    do ia = 1, 3
        kd(ia) = 1.d0
        kd(ia+3) = 2.d0/rac2
    end do

! - CALCUL POUR CHAQUE POINT DE GAUSS
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)
!
! - CALCUL DES ELEMENTS GEOMETRIQUES
! - CALCUL DE DFDI,F,EPS,R(EN AXI) ET POIDS
        call bmatmc(kpg, nbsig, zr(jvGeom), iw, ivf1, &
                    idff1, nno1, 0.d0, w, b)
!
        do ia = 1, 2*ndim
            do ja = 1, nno1
                do na = 1, ndim
                    def(ia, ja, na) = b(ia, (ja-1)*ndim+na)*kd(ia)
                end do
            end do
        end do
!
! - CALCUL DE TRACE(B)
        do na = 1, nno1
            do ia = 1, ndim
                deftr(na, ia) = def(1, na, ia)+def(2, na, ia)+def(3, na, ia)
            end do
        end do

! - CALCUL DE LA MATRICE DE HOOKE (LE MATERIAU POUVANT
! - ETRE ISOTROPE, ISOTROPE-TRANSVERSE OU ORTHOTROPE)
        call dmatmc(materPara, "+", notime, &
                    nbsig, dsidep)
!
        b_n = to_blas_int(2*ndim-3)
        b_incx = to_blas_int(1)
        call dscal(b_n, rac2, dsidep(4, 1), b_incx)
        b_n = to_blas_int(2*ndim-3)
        b_incx = to_blas_int(1)
        call dscal(b_n, rac2, dsidep(4, 2), b_incx)
        b_n = to_blas_int(2*ndim-3)
        b_incx = to_blas_int(1)
        call dscal(b_n, rac2, dsidep(4, 3), b_incx)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, rac2, dsidep(1, 4), b_incx)
        b_n = to_blas_int(2*ndim-3)
        b_incx = to_blas_int(1)
        call dscal(b_n, 2.d0, dsidep(4, 4), b_incx)
        if (ndim .eq. 3) then
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            call dscal(b_n, rac2, dsidep(1, 5), b_incx)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            call dscal(b_n, rac2, dsidep(1, 6), b_incx)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            call dscal(b_n, 2.d0, dsidep(4, 5), b_incx)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            call dscal(b_n, 2.d0, dsidep(4, 6), b_incx)
        end if
!
        devd(:, :) = 0.d0
        ddev(:, :) = 0.d0
        dddev(:, :) = 0.d0
        if (ndim .eq. 3) then
            devd(1:6, 1:6) = matmul(idev/3.d0, dsidep(1:6, 1:6))
            ddev(1:6, 1:6) = matmul(dsidep(1:6, 1:6), idev/3.d0)
            dddev(1:6, 1:6) = matmul(devd(1:6, 1:6), idev/3.d0)
        else
            devd(1:4, 1:4) = matmul(idev2/3.d0, dsidep(1:4, 1:4))
            ddev(1:4, 1:4) = matmul(dsidep(1:4, 1:4), idev2/3.d0)
            dddev(1:4, 1:4) = matmul(devd(1:4, 1:4), idev2/3.d0)
        end if
!
! - CALCUL DE D^DEV:ID ET ID:D^DEV ET ID:D:ID/9.D0
        iddid = 0.d0
        do ia = 1, 2*ndim
            devdi(ia) = devd(ia, 1)+devd(ia, 2)+devd(ia, 3)
            iddev(ia) = ddev(1, ia)+ddev(2, ia)+ddev(3, ia)
            do ja = 1, 3
                iddid = iddid+kr(ia)*dsidep(ia, ja)
            end do
        end do
        iddid = iddid/9.d0
!
! - CALCUL DE LA MATRICE DE RIGIDITE
! - TERME K:UX
        do na = 1, nno1
            do ia = 1, ndim
                vuiana = vu(ia, na)
                os = (vuiana-1)*vuiana/2
!
! - TERME K:UU      KUU(NDIM,NNO1,NDIM,NNO1)
                do nb = 1, nno1
                    do ib = 1, ndim
                        if (vu(ib, nb) .le. vuiana) then
                            kk = os+vu(ib, nb)
                            t1 = 0.d0
                            do ja = 1, 2*ndim
                                do jb = 1, 2*ndim
                                    t1 = t1+def(ja, na, ia)*dddev(ja, jb)*def(jb, nb, ib)
                                end do
                            end do
                            matr(kk) = matr(kk)+w*t1
                        end if
                    end do
                end do
!
! - TERME K:UG      KUG(NDIM,NNO1,NNO2)
                t1 = 0.d0
                do ja = 1, 2*ndim
                    t1 = t1+def(ja, na, ia)*devdi(ja)
                end do
                t1 = t1/3.d0
!
                do rb = 1, nno2
                    if (vg(rb) .lt. vuiana) then
                        kk = os+vg(rb)
                        matr(kk) = matr(kk)+w*t1*vff2(rb, kpg)
                    end if
                end do
!
! - TERME K:UP      KUP(NDIM,NNO1,NNO3)
                do sb = 1, nno3
                    if (vp(sb) .lt. vuiana) then
                        kk = os+vp(sb)
                        t1 = deftr(na, ia)*vff3(sb, kpg)
                        matr(kk) = matr(kk)+w*t1
                    end if
                end do
            end do
        end do
!
! - TERME K:GX
        do ra = 1, nno2
            vgra = vg(ra)
            os = (vgra-1)*vgra/2
!
! - TERME K:GU      KGU(NDIM,NNO2,NNO1)
            do nb = 1, nno1
                do ib = 1, ndim
                    if (vu(ib, nb) .lt. vgra) then
                        kk = os+vu(ib, nb)
                        t1 = 0.d0
                        do jb = 1, 2*ndim
                            t1 = t1+iddev(jb)*def(jb, nb, ib)
                        end do
                        matr(kk) = matr(kk)+w*t1*vff2(ra, kpg)/3.d0
                    end if
                end do
            end do
!
! - TERME K:GG      KGG(NNO2,NNO2)
            do rb = 1, nno2
                if (vg(rb) .le. vgra) then
                    kk = os+vg(rb)
                    t1 = vff2(ra, kpg)*iddid*vff2(rb, kpg)
                    matr(kk) = matr(kk)+w*t1
                end if
            end do
!
! - TERME K:GP      KGP(NNO2,NNO3)
            do sb = 1, nno3
                if (vp(sb) .lt. vgra) then
                    kk = os+vp(sb)
                    t1 = -vff2(ra, kpg)*vff3(sb, kpg)
                    matr(kk) = matr(kk)+w*t1
                end if
            end do
        end do
!
! - TERME K:PX
        do sa = 1, nno3
            vpsa = vp(sa)
            os = (vpsa-1)*vpsa/2
!
! - TERME K:PU      KPU(NDIM,NNO3,NNO1)
            do nb = 1, nno1
                do ib = 1, ndim
                    if (vu(ib, nb) .lt. vpsa) then
                        kk = os+vu(ib, nb)
                        t1 = vff3(sa, kpg)*deftr(nb, ib)
                        matr(kk) = matr(kk)+w*t1
                    end if
                end do
            end do
!
! - TERME K:PG      KPG(NNO3,NNO2)
            do rb = 1, nno2
                if (vg(rb) .lt. vpsa) then
                    kk = os+vg(rb)
                    t1 = -vff3(sa, kpg)*vff2(rb, kpg)
                    matr(kk) = matr(kk)+w*t1
                end if
            end do
!
! - TERME K:PP = 0.D0      KPP(NNO3,NNO3)
        end do
    end do
end subroutine
