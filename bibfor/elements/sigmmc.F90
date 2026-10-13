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
subroutine sigmmc(materPara, &
                  nno, ndim, nbsig, npg, &
                  jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                  nodeCoor, nodeDisp, &
                  time, nharm, &
                  sigm)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/bmatmc.h"
#include "asterfort/dbudef.h"
#include "asterfort/dmatmc.h"
#include "jeveux.h"
#include "FE_module.h"
#include "MeshTypes_type.h"
!
    type(Material_Para), intent(inout) :: materPara
    integer(kind=8), intent(in) :: nno, ndim, nbsig, npg
    integer(kind=8), intent(in) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    real(kind=8), intent(in) :: nodeCoor(ndim*nno), nodeDisp(ndim*nno)
    real(kind=8), intent(in) :: time
    real(kind=8), intent(in) :: nharm
    real(kind=8), intent(out) :: sigm(nbsig*npg)
!
! --------------------------------------------------------------------------------------------------
!
! Compute stress at Gauss points
!
! --------------------------------------------------------------------------------------------------
!
! IO  materPara        : parameters of material
! In  nno              : number of nodes of element
! In  ndim             : dimension of element (2 ou 3)
! In  nbsig            : number of components for stress tensors (4 or 6)
! In  npg              : number of Gauss points
! In  jvGaussWeight    : adresse to Gauss point weight
! In  jvBaseFunc       : adresse to shape functions
! In  jvDBaseFunc      : adresse to derivative of shape functions
! In  nodeCoor         : coordinates of nodes
! In  nodeDisp         : displacements at nodes
! In  time             : current time
! In  nharm            : Fourier mode
! Out sigm             : values of stress tensor at integration points
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: kpg, nbinco
    real(kind=8) :: b(6*MAX_BV_CG), d(36), jacgau
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)
    ASSERT(nno .le. MT_NNOMAX)
    nbinco = ndim*nno
    sigm = 0.d0

! - Loop on Gauss points
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! ----- Compute [B] matrix (displacements to strains)
        call bmatmc(kpg, nbsig, nodeCoor, &
                    jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                    nno, nharm, jacgau, b)

! ----- Compute elasticity matrix
        call dmatmc(materPara, '+', time, &
                    nbsig, d)

! ----- Compute SIEF_ELGA
        call dbudef(nodeDisp, b, d, nbsig, nbinco, &
                    sigm(1+nbsig*(kpg-1)))
    end do
!
end subroutine
