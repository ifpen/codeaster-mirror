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

subroutine te0561(option, nomte)
!
    use FE_algebra_module
    use HHO_basis_module
    use HHO_eval_module
    use HHO_geometry_module
    use HHO_init_module, only: hhoInfoInitFace, hhoUpdateCoordinatesFace
    use HHO_quadrature_module
    use HHO_size_module, only: hhoMecaFaceDofs
    use HHO_type
    use HHO_utils_module
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/tefrep.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/readVector.h"
!
    character(len=16), intent(in) :: option, nomte
!
!---------------------------------------------------------------------------------------------------
!
! HHO METHODS
!
! Tangent matrix associated with a follower pressure.
!
! The HHO displacement is first evaluated at the Lagrange geometrical
! nodes. The current surface geometry and its tangent vectors are then
! computed with the Lagrange interpolation.
!
! Consequently, the tangent matrix must differentiate the chain:
!
! HHO coefficients
! -> nodal Lagrange displacements
! -> current geometrical tangent vectors
! -> current oriented surface vector
! -> follower-pressure load vector
!
!---------------------------------------------------------------------------------------------------
!
    type(HHO_Data) :: hhoData
    type(HHO_Face) :: hhoFace
    type(HHO_Face) :: hhoFaceDef
    type(HHO_Quadrature) :: hhoQuadFace
    type(HHO_basis_face) :: hhoBasisFace
!
    real(kind=8) :: lhs(MSIZE_FACE_VEC, MSIZE_FACE_VEC)
    real(kind=8) :: PresQP(MSIZE_QP_FACE)
    real(kind=8) :: BSEval(MSIZE_FACE_SCAL)
    real(kind=8) :: BSEvalNodes(MSIZE_FACE_SCAL, 9)
    real(kind=8), dimension(MSIZE_FACE_VEC) :: depl_prev, depl_incr, depl_curr
    real(kind=8) :: depl_nodes(3, 9), dNGeo(3, 8)
    real(kind=8) :: tau1(3), tau2(3)
    real(kind=8) :: geomGrad(2, MSIZE_FACE_SCAL)
    real(kind=8) :: skewTau1(3, 3), skewTau2(3, 3), block_ij(3, 3)
    real(kind=8) :: coef
!
    integer(kind=8) :: fbs, fbs_scal, celldim
    integer(kind=8) :: ipg, ino
    integer(kind=8) :: ibasis, jbasis
    integer(kind=8) :: idim, jdim
    integer(kind=8) :: irow, jcol
    integer(kind=8) :: j_pres
!
!---------------------------------------------------------------------------------------------------
! Retrieve HHO data
!---------------------------------------------------------------------------------------------------
!
    call hhoInfoInitFace( &
        hhoFace, hhoData, hhoQuadFace=hhoQuadFace)
!
    ASSERT(hhoQuadFace%nbQuadPoints <= MSIZE_QP_FACE)
    ASSERT(.not. hhoFace%l_axis)
!
    celldim = hhoFace%ndim+1
!
    call hhoMecaFaceDofs(hhoFace, hhoData, fbs)
    fbs_scal = fbs/celldim
!
    lhs = 0.d0
    PresQP = 0.d0
    BSEvalNodes = 0.d0
!
!---------------------------------------------------------------------------------------------------
! Current HHO displacement coefficients
!---------------------------------------------------------------------------------------------------
!
    call readVector('PDEPLMR', fbs, depl_prev)
    call readVector('PDEPLPR', fbs, depl_incr)
!
    call dcopy_1(fbs, depl_prev, depl_curr)
    call daxpy_1(fbs, 1.d0, depl_incr, depl_curr)
!
    call hhoBasisFace%initialize(hhoFace)
!
!---------------------------------------------------------------------------------------------------
! Evaluate the HHO basis functions at the Lagrange geometrical nodes
!---------------------------------------------------------------------------------------------------
!
    do ino = 1, hhoFace%nbnodes
        call hhoBasisFace%BSEval(hhoFace%coorno(1:3, ino), &
                                 0, hhoData%face_degree(), BSEvalNodes(1:, ino))
    end do
!
!---------------------------------------------------------------------------------------------------
! Compute the current nodal displacement and update the Lagrange geometry
!---------------------------------------------------------------------------------------------------
!
    depl_nodes = 0.d0
!
    do ino = 1, hhoFace%nbnodes
        depl_nodes(1:3, ino) = hhoEvalVecFace(hhoBasisFace, hhoData%face_degree(), &
                                              hhoFace%coorno(1:3, ino), depl_curr)
    end do
!
    hhoFaceDef = hhoFace
    call hhoUpdateCoordinatesFace(depl_nodes, hhoFaceDef)
!
!---------------------------------------------------------------------------------------------------
! Option
!---------------------------------------------------------------------------------------------------
!
    if (option == 'RIGI_MECA_PRSU_R') then
!
!---------------------------------------------------------------------------------------------------
! Evaluate pressure on the initial face quadrature
!---------------------------------------------------------------------------------------------------
!
        call jevech('PPRESSR', 'L', j_pres)
        call hhoFuncRScalEvalQp(hhoFace, hhoQuadFace, zr(j_pres), PresQP)
!
        do ipg = 1, hhoQuadFace%nbQuadPoints
!
!------------------------------------------------------------------------------------------------
! Current covariant tangent vectors computed with the Lagrange geometry
!------------------------------------------------------------------------------------------------
!
            call hhoTangentFaceQP(hhoFaceDef, hhoQuadFace%points_param(1:2, ipg), &
                                  tau1, tau2)
!
!------------------------------------------------------------------------------------------------
! HHO test basis on the initial face
!------------------------------------------------------------------------------------------------
!
            call hhoBasisFace%BSEval(hhoQuadFace%points(1:3, ipg), 0, hhoData%face_degree(), &
                                     BSEval)
!
!------------------------------------------------------------------------------------------------
! Differentiate the actual Lagrange geometry with respect to each HHO mode
!------------------------------------------------------------------------------------------------
!
            call hhoGeomDerivBasis(hhoFace%typema, hhoQuadFace%points_param(1:2, ipg), &
                                   dNGeo)
! --- Renumber derivative since normal and tangent can change of sign.
            if (hhoFace%reorientation) then
                if (hhoFace%ndim == 1) then
                    dNGeo(1, 1:8) = -dNGeo(1, 1:8)
                else
                    dNGeo(3, 1:8) = dNGeo(1, 1:8)
                    dNGeo(1, 1:8) = dNGeo(2, 1:8)
                    dNGeo(2, 1:8) = dNGeo(3, 1:8)
                    dNGeo(3, 1:8) = 0.d0
                end if
            end if
!
            geomGrad = 0.d0
            do jbasis = 1, fbs_scal
                do ino = 1, hhoFace%nbnodes
                    geomGrad(1, jbasis) = geomGrad(1, jbasis)+ &
                                          dNGeo(1, ino)*BSEvalNodes(jbasis, ino)
!
                    if (hhoFace%ndim == 2) then
                        geomGrad(2, jbasis) = geomGrad(2, jbasis) &
                                              +dNGeo(2, ino)*BSEvalNodes(jbasis, ino)
                    end if
                end do
            end do
!
!------------------------------------------------------------------------------------------------
! Cross-product matrix associated with current tau
!
! skewTau * v = tau x v
!------------------------------------------------------------------------------------------------
!
            skewTau1 = 0.d0
            skewTau2 = 0.d0
!
            if (hhoFace%ndim == 2) then
                skewTau1(1, 2) = -tau1(3)
                skewTau1(1, 3) = tau1(2)
                skewTau1(2, 1) = tau1(3)
                skewTau1(2, 3) = -tau1(1)
                skewTau1(3, 1) = -tau1(2)
                skewTau1(3, 2) = tau1(1)
!
                skewTau2(1, 2) = -tau2(3)
                skewTau2(1, 3) = tau2(2)
                skewTau2(2, 1) = tau2(3)
                skewTau2(2, 3) = -tau2(1)
                skewTau2(3, 1) = -tau2(2)
                skewTau2(3, 2) = tau2(1)
            else
                skewTau1(1, 2) = 1.d0
                skewTau1(2, 1) = -1.d0
            end if
!
            coef = hhoQuadFace%weights_param(ipg)*PresQP(ipg)
!
!------------------------------------------------------------------------------------------------
! Assemble the pressure tangent
!
! Kp_ij = integral_F0 p*phi_i*[ -D_1j [tau2]_x + D_2j [tau1]_x ] dS0
!
!------------------------------------------------------------------------------------------------
!
            do ibasis = 1, fbs_scal
                do jbasis = 1, fbs_scal
                    if (hhoFace%ndim == 2) then
                        block_ij = coef*BSEval(ibasis)*( &
                                   -geomGrad(1, jbasis)*skewTau2+geomGrad(2, jbasis)*skewTau1)
                    else
                        block_ij(1:2, 1:2) = coef*BSEval(ibasis) &
                                             *geomGrad(1, jbasis)*skewTau1(1:2, 1:2)
                    end if
                    do idim = 1, celldim
                        irow = (idim-1)*fbs_scal+ibasis
                        do jdim = 1, celldim
                            jcol = (jdim-1)*fbs_scal+jbasis
                            lhs(irow, jcol) = lhs(irow, jcol)+block_ij(idim, jdim)
!
                        end do
                    end do
                end do
            end do
!
        end do
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
!---------------------------------------------------------------------------------------------------
! Save the non-symmetric matrix
!---------------------------------------------------------------------------------------------------
!
    call writeMatrix('PMATUNS', fbs, fbs, ASTER_FALSE, lhs)
!
end subroutine
