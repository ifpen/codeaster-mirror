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
subroutine sigvmc(materPara, &
                  nno, ndim, nbsig, npg, &
                  jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                  nodeCoor, nodeDisp, &
                  time, nharm, &
                  sigm)
!
    use BehaviourStrain_type
    use MaterialPara_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/sigmmc.h"
#include "asterfort/sigtmc.h"
#include "jeveux.h"
#include "FE_module.h"
#include "MeshTypes_type.h"
!
    type(Material_Para), intent(inout) :: materPara
    integer(kind=8), intent(in) :: nno, ndim, nbsig, npg
    integer(kind=8), intent(in) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    real(kind=8), intent(in) :: nodeCoor(ndim*nno), nodeDisp(ndim*nno)
    real(kind=8), intent(in) :: time, nharm
    real(kind=8), intent(out) :: sigm(nbsig*npg)
!
! --------------------------------------------------------------------------------------------------
!
! Compute mechanical stress (without effect of external state variables)
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
    real(kind=8) :: sigmVarc(162)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)
    sigm = 0.d0

! - Compute total stress
    call sigmmc(materPara, &
                nno, ndim, nbsig, npg, &
                jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                nodeCoor, nodeDisp, &
                time, nharm, &
                sigm)

! - Compute stresses from external state variables
    call sigtmc(materPara, time, &
                nbsig, npg, ndim, &
                VARC_STRAIN_ALL, sigmVarc)

! - CALCUL DES CONTRAINTES TOTALES AUX POINTS D'INTEGRATION
    sigm(1:nbsig*npg) = sigm(1:nbsig*npg)-sigmVarc(1:nbsig*npg)
!
end subroutine
