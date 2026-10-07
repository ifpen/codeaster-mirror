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
module coupling_quadrature_module
!
    use FE_quadrature_module
    use FE_topo_module
    use HHO_quadrature_module
    use HHO_type
    use HHO_utils_module
!
    implicit none
!
    private
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/apnorm.h"
#include "asterfort/assert.h"
#include "asterfort/coupling_type.h"
#include "asterfort/coupling_type.h"
#include "asterfort/jevech.h"
#include "asterfort/latrco.h"
#include "asterfort/lcptga.h"
#include "asterfort/mmdonf.h"
#include "asterfort/mmmjac.h"
#include "asterfort/mmnewd.h"
#include "asterfort/mmnonf.h"
#include "asterfort/projInsideCell.h"
#include "MeshTypes_type.h"
#include "asterfort/reerel.h"
!
! --------------------------------------------------------------------------------------------------
!
! Coupling - Quadrature
!
! --------------------------------------------------------------------------------------------------
!
!
#define MAX_NB_INTE  8
#define PROJ_TOLE  1.0d-10

!
    private :: cplGetQuad, cplProjFEM2HHO
    public :: cplGetQuadFEMHHO, cplGetQuadFEMFEM, cplGetQuadOrderFEM
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
!
    integer(kind=8) function cplGetQuadOrderFEM(typema)
!
        implicit none
!
        character(len=8), intent(in) :: typema
!
!===================================================================================================
!    FEM/HHO - Get quadrature
!===================================================================================================
!
        select case (typema)
        case ("SE2", "TR3", "QU4")
            cplGetQuadOrderFEM = 1
        case ("SE3", "TR6", "TR7", "QU8", "QU9")
            cplGetQuadOrderFEM = 2
        case ("SE4", "TR1", "Q12")
            cplGetQuadOrderFEM = 3
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplGetQuad(FEFaceSl, max_order, FEQuadSl, order_)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        integer(kind=8), intent(in) :: max_order
        type(FE_Quadrature), intent(out) :: FEQuadSl
        integer(kind=8), intent(out), optional :: order_
!
!===================================================================================================
!    FEM/HHO - Get quadrature
!===================================================================================================
!
        integer(kind=8) :: nbPoinInte, jpair, iPoinInte, ndim
        integer(kind=8) :: nbPtGauss, nb_qp, i, order
        integer(kind=8) :: iTria, iGauss, nbTria, nbGauss
        real(kind=8) :: triaCoorSlav(2, 3), jaco, coorac(3)
        real(kind=8) :: gausWeightSlav(12), gausCoorSlav(2, 12)
        real(kind=8) :: shape_func(MT_NNOMAX2D), shape_dfunc(2, MT_NNOMAX2D)
        character(len=8) :: elga_fami
        real(kind=8) :: poinInteSlav(2, MAX_NB_INTE)
!
        call jevech('PPAIRR', 'L', jpair)
        nbPoinInte = int(zr(jpair-1+OFFSET_NB_PTS_INTER))
        ASSERT(nbPoinInte .le. MAX_NB_INTE)
        poinInteSlav = 0.d0
        do iPoinInte = 1, nbPoinInte
            poinInteSlav(1, iPoinInte) = zr(jpair-1+OFFSET_COOR_PTS_X+iPoinInte-1)
            poinInteSlav(2, iPoinInte) = zr(jpair-1+OFFSET_COOR_PTS_Y+iPoinInte-1)
        end do
!
        ndim = FEFaceSl%ndim+1
!
! - Triangulation of convex polygon defined by intersection points
        if (ndim .eq. 3) then
            if (nbPoinInte == 3) then
                nbTria = 1
            else
                nbTria = nbPoinInte
            end if
            select case (max_order)
            case (1, 2)
                ! order 2 - by triangle
                elga_fami = 'FPG3'
                nbPtGauss = 3
                order = 2
            case (3)
                ! order 3 - by triangle
                elga_fami = 'FPG4'
                nbPtGauss = 4
                order = 2
            case (4, 5)
                ! order 5 - by triangle
                elga_fami = 'FPG7'
                nbPtGauss = 7
                order = 5
            case (6, 7)
                ! order 7 - by triangle
                elga_fami = 'FPG13'
                nbPtGauss = 13
                order = 7
            case default
                ASSERT(ASTER_FALSE)
            end select

        elseif (ndim .eq. 2) then
            nbTria = 1
!
            select case (max_order)
            case (1, 2, 3)
                ! order 3
                elga_fami = 'FPG2'
                nbPtGauss = 2
                order = 3
            case (4, 5)
                ! order 5
                elga_fami = 'FPG3'
                nbPtGauss = 3
                order = 5
            case (6, 7)
                ! order 7
                elga_fami = 'FPG4'
                nbPtGauss = 4
                order = 7
            case (8, 9)
                ! order 7
                elga_fami = 'FPG5'
                nbPtGauss = 5
                order = 9
            case default
                ASSERT(ASTER_FALSE)
            end select
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (present(order_)) then
            order_ = order
        end if
!
        ASSERT(nbTria*nbPtGauss <= 64)
        FEQuadSl%nbQuadPoints = nbTria*nbPtGauss
        nb_qp = 0
!
! - Loop on triangles
        do iTria = 1, nbTria
! ----- Coordinates of current triangle (slave)
            triaCoorSlav = 0.d0
            if (ndim .eq. 3) then
                call latrco(iTria, nbPoinInte, poinInteSlav, triaCoorSlav)
            elseif (ndim .eq. 2) then
                ASSERT(nbPoinInte .eq. 2)
                triaCoorSlav(1:2, 1:2) = poinInteSlav(1:2, 1:2)
            end if

! ----- Get integration points for slave element
            call lcptga(ndim, triaCoorSlav, elga_fami, &
                        nbGauss, gausCoorSlav, gausWeightSlav)

! ----- Loop on integration points
            do iGauss = 1, nbGauss
                nb_qp = nb_qp+1
                FEQuadSl%points_param(1:2, nb_qp) = gausCoorSlav(1:2, iGauss)
                FEQuadSl%weights_param(nb_qp) = gausWeightSlav(iGauss)
!
! ------------- Get shape functions and first derivative only (for perf)
                call mmnonf(FEFaceSl%typemas, &
                            gausCoorSlav(1, iGauss), gausCoorSlav(2, iGauss), &
                            shape_func)
                call mmdonf(FEFaceSl%typemas, &
                            gausCoorSlav(1, iGauss), gausCoorSlav(2, iGauss), &
                            shape_dfunc)
!
                coorac = 0.d0
                do i = 1, FEFaceSl%nbnodes
                    coorac(1:3) = coorac(1:3)+FEFaceSl%coorno(1:3, i)*shape_func(i)
                end do
!
! ------------- Compute jacobian
                call mmmjac(ASTER_FALSE, FEFaceSl%nbnodes, ndim, &
                            FEFaceSl%typemas, FEFaceSl%coorno, &
                            shape_func, shape_dfunc, jaco)
!
                FEQuadSl%points(1:3, nb_qp) = coorac
                FEQuadSl%weights(nb_qp) = abs(jaco)*gausWeightSlav(iGauss)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplProjFEM2HHO(FEFaceSl, hhoFaceMa, FEQuadSl, order, hhoQuadMa)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        integer(kind=8), intent(in) :: order
        type(HHO_Quadrature), intent(out) :: hhoQuadMa
!
!===================================================================================================
!    FEM/HHO - Get quadrature
!===================================================================================================
!
        integer(kind=8) :: ndim, ipg, iret
        real(kind=8) :: norm_slav(3), tau_slav(3, 2), coor_qp_sl(2)
        real(kind=8) :: coor_qp_sl_re(3), ksi_line(2)
        real(kind=8) :: tau1_mast(3), tau2_mast(3), coor_nodes(3, MT_NNOMAX2D)
        character(len=8) :: type_mast
!
        hhoQuadMa%order = order
        hhoQuadMa%nbQuadPoints = FEQuadSl%nbQuadPoints
        ndim = FEFaceSl%ndim+1
!
        coor_nodes(1:3, 1:9) = hhoFaceMa%coorno
!
        call CellNameL2S(hhoFaceMa%typema, type_mast)
!
        do ipg = 1, hhoQuadMa%nbQuadPoints
            hhoQuadMa%weights(ipg) = FEQuadSl%weights(ipg)
            coor_qp_sl = FEQuadSl%points_param(1:2, ipg)
!
! ------ Compute outward slave normal (pairing configuration)
!
            call apnorm(FEFaceSl%nbnodes, FEFaceSl%typemas, ndim, FEFaceSl%coorno, &
                        coor_qp_sl(1), coor_qp_sl(2), norm_slav, tau_slav(1:3, 1), tau_slav(1:3, 2))
!
! ----- Return in real slave space (pairing configuration)
!
            coor_qp_sl_re = FEQuadSl%points(1:3, ipg)
!
! ----- Projection on master element
            call mmnewd(type_mast, hhoFaceMa%nbnodes, ndim, coor_nodes, &
                        coor_qp_sl_re, 100, PROJ_TOLE, norm_slav, ksi_line(1), &
                        ksi_line(2), tau1_mast, tau2_mast, iret)
            ASSERT(iret == 0)
!
! ----- Check that projected node is inside cell
!
            call projInsideCell(1e-6, ndim, type_mast, ksi_line, iret)
            ASSERT(iret == 0)
!
            hhoQuadMa%points_param(1:2, ipg) = ksi_line
            call reerel(type_mast, hhoFaceMa%nbnodes, 3, hhoFaceMa%coorno, &
                        hhoQuadMa%points_param(1:3, ipg), &
                        hhoQuadMa%points(1:3, ipg))
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplProjFEM2FEM(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(FE_Quadrature), intent(out) :: FEQuadMa
!
!===================================================================================================
!    FEM/FEM - Get quadrature
!===================================================================================================
!
        integer(kind=8) :: ndim, ipg, iret
        real(kind=8) :: norm_slav(3), tau_slav(3, 2), coor_qp_sl(2)
        real(kind=8) :: coor_qp_sl_re(3), ksi_line(2)
        real(kind=8) :: tau1_mast(3), tau2_mast(3)
!
        FEQuadMa%nbQuadPoints = FEQuadSl%nbQuadPoints
        ndim = FEFaceSl%ndim+1
!
        do ipg = 1, FEQuadMa%nbQuadPoints
            FEQuadMa%weights(ipg) = FEQuadSl%weights(ipg)
            coor_qp_sl = FEQuadSl%points_param(1:2, ipg)
!
! ------ Compute outward slave normal (pairing configuration)
!
            call apnorm(FEFaceSl%nbnodes, FEFaceSl%typemas, ndim, FEFaceSl%coorno, &
                        coor_qp_sl(1), coor_qp_sl(2), norm_slav, tau_slav(1:3, 1), tau_slav(1:3, 2))
!
! ----- Return in real slave space (pairing configuration)
!
            coor_qp_sl_re = FEQuadSl%points(1:3, ipg)
!
! ----- Projection on master element
            call mmnewd(FEFaceMa%typemas, FEFaceMa%nbnodes, ndim, FEFaceMa%coorno, &
                        coor_qp_sl_re, 100, PROJ_TOLE, norm_slav, ksi_line(1), &
                        ksi_line(2), tau1_mast, tau2_mast, iret)
            ASSERT(iret == 0)
!
! ----- Check that projected node is inside cell
!
            call projInsideCell(1e-6, ndim, FEFaceMa%typemas, ksi_line, iret)
            ASSERT(iret == 0)
!
            FEQuadMa%points_param(1:2, ipg) = ksi_line
            call reerel(FEFaceMa%typemas, FEFaceMa%nbnodes, 3, FEFaceMa%coorno, &
                        FEQuadMa%points_param(1:3, ipg), &
                        FEQuadMa%points(1:3, ipg))
        end do
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplGetQuadFEMHHO(FEFaceSl, hhoFaceMa, hhoData, FEQuadSl, hhoQuadMa)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoData
        type(FE_Quadrature), intent(out) :: FEQuadSl
        type(HHO_Quadrature), intent(out) :: hhoQuadMa
!
!===================================================================================================
!    FEM/HHO - Get quadrature
!===================================================================================================
!
        integer(kind=8) :: max_order, order_used
!
!       We add +1 for quadrature by security
        max_order = 2*(max(hhoData%face_degree(), cplGetQuadOrderFEM(FEFaceSl%typemas)))+1
!
        call cplGetQuad(FEFaceSl, max_order, FEQuadSl, order_used)
        call cplProjFEM2HHO(FEFaceSl, hhoFaceMa, FEQuadSl, order_used, hhoQuadMa)
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplGetQuadFEMFEM(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa
        type(FE_Quadrature), intent(out) :: FEQuadSl, FEQuadMa
!
!===================================================================================================
!    FEM/HHO - Get quadrature
!===================================================================================================
!
        integer(kind=8) :: max_order
!
        max_order = 2*(max(cplGetQuadOrderFEM(FEFaceSl%typemas), &
                           cplGetQuadOrderFEM(FEFaceMa%typemas)))+1
!
        call cplGetQuad(FEFaceSl, max_order, FEQuadSl)
        call cplProjFEM2FEM(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
end module
