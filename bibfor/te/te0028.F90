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

subroutine te0028(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use Behaviour_module, only: behaviourOption
!
    use c_interface_plaq_mitc_j
    use iso_c_binding

    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmdlog.h"
#include "asterfort/nmgpfi.h"
#include "asterfort/nmgrla.h"
#include "asterfort/nmplxd.h"
#include "asterfort/nmtstm.h"
#include "asterfort/rcangm.h"
#include "asterfort/tecach.h"
#include "asterfort/tgveri.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "FE_module.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/dxroep.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/utpvgl.h"
#include "asterfort/dxtpgl.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------
! Elementary computation
!
! Elements: PLAQ_MITC
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------
!
    character(len=4) :: fami
    integer(kind=8) :: ndim, nno, igeom, imate, i, j
    integer(kind=8), parameter :: size_fenicsx = 42*42, size_aster = 30*30
    real(c_double) :: cst(4), coor(27), kappa, cdofs_f(12)
    ! NOTICE: see the size of the arrays in the C file: c_interface_plaq_mitc_j
    real(c_double), dimension(size_fenicsx) :: A0, A1, A2, A3, A4
    real(c_double), dimension(size_aster) :: A_int

    integer(c_int) :: ncst, ncd, ne0, ne1, ne2, ne3, nwinit
    integer(c_int) :: entities0(1), entities1(1), entities2(1), entities3(1)
!
    integer(kind=8), parameter :: size_mat = 30
    integer(kind=8) :: reorder(size_mat)
    real(kind=8) :: signs(size_mat)
    real(c_double) :: w_0(size_mat)
    real(kind=8), dimension(size_mat, size_mat) :: bint, bf
    real(kind=8) :: e, nu, e_, nu_, epais, rho, tiny
    parameter(tiny=1.d-8)
    integer(kind=8) :: elas_id, igau, ipoids, npg1
    integer(kind=8) :: nnos, ivf, idfde
    character(len=16) :: elas_keyword
! ---------------------------------------------------------------------
!
! - Finite element informations
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg1, jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
! - Initializations
!
    ndim = 2
    cst = 0.d0
    coor = 0.d0
    w_0 = 0.d0
    cdofs_f = 0.d0
    A_int = 0.d0
    A0 = 0.d0
    A1 = 0.d0
    A2 = 0.d0
    A3 = 0.d0
    A4 = 0.d0
    bint = 0.d0
    bf = 0.d0
    ![w1, θ_y1, -θ_x1, w2, θ_y2, -θ_x2, w3, θ_y3, -θ_x3, w4, θ_y4, -θ_x4,
    ! θ_y5, -θ_x5, γ_r1, p1, θ_y6, -θ_x6, γ_r2, p2,θ_y7, -θ_x7,
    ! γ_r3, p3, θ_y8, -θ_x8,γ_r4, p4, θ_y9, -θ_x9]
    signs = (/ &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, 1.d0, &
            1.d0, -1.d0, &
            1.d0, 1.d0, &
            1.d0, -1.d0, &
            1.d0, 1.d0, &
            1.d0, -1.d0, &
            1.d0, 1.d0, &
            1.d0, -1.d0 &
            /)
    reorder = (/ &
              19, 2, 1, &
              21, 6, 5, &
              22, 8, 7, &
              20, 4, 3, &
              12, 11, &
              24, 28, &
              16, 15, &
              26, 30, &
              14, 13, &
              25, 29, &
              10, 9, &
              23, 27, &
              18, 17 &
              /)
!
! - Geometry
!
    call jevech('PGEOMER', 'L', igeom)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', imate)

! ----- Get elastic parameters (only isotropic elasticity)
! FIXME: for instance E and NU are supposed to be equal for all quadpoints
!
    call get_elas_id(zi(imate), elas_id, elas_keyword)
    call get_elas_para(fami, zi(imate), '+', 1, 1, &
                       elas_id, elas_keyword, &
                       e_=e, nu_=nu)

    do igau = 2, npg1
        call get_elas_id(zi(imate), elas_id, elas_keyword)
        call get_elas_para(fami, zi(imate), '+', igau, 1, &
                           elas_id, elas_keyword, &
                           e_=e_, nu_=nu_)
        ASSERT(abs(e-e_) .le. tiny .and. abs(nu-nu_) .le. tiny)
    end do

    call dxroep(rho, epais)

! - Fill material parameters vector
    kappa = 5.0/6.0
    cst(1) = e
    cst(2) = nu
    cst(3) = kappa
    cst(4) = epais

! - Fill material entities vector
! - 4 entities car 4 arêtes dans un carré
    entities0(1) = 0
    entities1(1) = 1
    entities2(1) = 2
    entities3(1) = 3

! Remplissage du vecteur de coordonnées (3 coordonnées par nœud)
    do i = 0, 8
        coor(3*i+1) = zr(igeom+3*i)
        coor(3*i+2) = zr(igeom+3*i+1)
        coor(3*i+3) = zr(igeom+3*i+2)
    end do

    ! Remplissage des degrés (N1 , N2, N4, N3)
    cdofs_f(1) = coor(1)
    cdofs_f(2) = coor(2)
    cdofs_f(3) = coor(3)
!
    cdofs_f(4) = coor(10)
    cdofs_f(5) = coor(11)
    cdofs_f(6) = coor(12)
!
    cdofs_f(7) = coor(4)
    cdofs_f(8) = coor(5)
    cdofs_f(9) = coor(6)
!
    cdofs_f(10) = coor(7)
    cdofs_f(11) = coor(8)
    cdofs_f(12) = coor(9)
!
    nwinit = size(w_0)
    ncd = size(cdofs_f)
    ne0 = 1
    ne1 = 1
    ne2 = 1
    ne3 = 1
    ncst = size(cst)
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!MATRICE DE RIGIDITÉ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Appel de la fonction C++ - Matrice de rigidité -> Part1
    call BP4_qu9_Fortran(w_0, nwinit, cdofs_f, ncd, entities0, ne0, cst, ncst, A0)
! Appel de la fonction C++ - Matrice de rigidité -> Part2
    call BP5_qu9_Fortran(w_0, nwinit, cdofs_f, ncd, entities0, ne0, cst, ncst, A1)
    call BP5_qu9_Fortran(w_0, nwinit, cdofs_f, ncd, entities1, ne1, cst, ncst, A2)
    call BP5_qu9_Fortran(w_0, nwinit, cdofs_f, ncd, entities2, ne2, cst, ncst, A3)
    call BP5_qu9_Fortran(w_0, nwinit, cdofs_f, ncd, entities3, ne3, cst, ncst, A4)
!
! Remplissage de la matrice intermédiaire
    do i = 1, size_aster
        A_int(i) = A0(i)+A1(i)+A2(i)+A3(i)+A4(i)
    end do
! Remplissage de la matrice K à partir de A_int
    do i = 1, size_mat
        do j = 1, size_mat
            bint(i, j) = A_int((j-1)*size_mat+i)
        end do
    end do
! Reorganisation de la matrice
    do i = 1, size_mat
        do j = 1, size_mat
            bf(i, j) = signs(i)*signs(j)*bint(reorder(i), reorder(j))
        end do
    end do
!
    call writeMatrix('PMATUUR', size_mat, size_mat, ASTER_TRUE, bf)

end subroutine
