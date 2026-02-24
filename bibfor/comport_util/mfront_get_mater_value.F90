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
subroutine mfront_get_mater_value(extern_addr, BEHInteg, relaComp, &
                                  fami, kpg, ksp, jvMaterCode, &
                                  propVale, nbProp)
!
    use Behaviour_type
    implicit none
!
#include "asterc/mgis_get_props.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/mat_proto.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: extern_addr
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=16), intent(in) :: relaComp
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode, nbProp
    real(kind=8), intent(out) :: propVale(nbProp)
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour (MFront)
!
! Get material properties
!
! --------------------------------------------------------------------------------------------------
!
! In  extern_addr       : address of the MGISBehaviour object
! In  BEHInteg         : parameters for integration of behaviour
! In  relaComp         : RELATION comportment
! In  fami             : Gauss family for integration point rule
! In  jvMaterCode      : coded material address
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  nbProp           : number of material properties
! Out propVale         : values of material properties
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: propCode(MGIS_MAX_PROPS), iProp
    character(len=16) :: propName(MGIS_MAX_PROPS)
    real(kind=8) :: zSumCold
    integer(kind=8) :: nbPhase, metaType, nbcoef
    integer(kind=8), parameter :: nbPara = 3
    real(kind=8) :: paraVale(nbPara)
    character(len=16), parameter :: paraName(nbPara) = (/'X', 'Y', 'Z'/)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbProp <= MGIS_MAX_PROPS)

! - Coordinates of current Gauss point
    paraVale = BEHInteg%behavESVA%behavESVAGeom%coorElga(kpg, :)

! - Get name of properties from MGIS
    call mgis_get_props(extern_addr, propName)

! - Get parameters
    if (relaComp .eq. 'MFRONT' .or. relaComp .eq. 'CZM_MFRONT') then
        nbcoef = nbProp
        call mat_proto(BEHInteg, &
                       fami, kpg, ksp, '+', jvMaterCode, 'MFRONT', &
                       nbcoef, propVale)
        ASSERT(nbcoef == nbProp)
    else
! ----- Get the properties values (enter under 'relaComp' in DEFI_MATERIAU)
        propVale(1:nbProp) = r8nnem()
        if (BEHInteg%behavESVA%tabcod(ZFERRITE) .eq. 1) then
            metaType = META_STEEL
            nbPhase = PSTEEL_NB
            call metaGetPhase(fami, '+', kpg, ksp, metaType, &
                              nbPhase, zSumCold_=zSumCold)
            do iProp = 1, nbProp
                if (propName(iProp) (1:4) .eq. 'meta') then
                    call rcvalb(fami, 1, 1, '+', jvMaterCode, &
                                ' ', relaComp, 1, 'META', [zSumCold], &
                                1, propName(iProp), propVale(iProp), propCode(iProp), 2)
                else
                    call rcvalb(fami, kpg, ksp, '+', jvMaterCode, &
                                ' ', relaComp, 0, ' ', [0.d0], &
                                1, propName(iProp), propVale(iProp), propCode(iProp), 1)
                end if
            end do
        elseif (BEHInteg%behavESVA%tabcod(ZALPHPUR) .eq. 1) then
            metaType = META_ZIRC
            nbPhase = PZIRC_NB
            call utmess('F', 'COMPOR4_24')
        else
            if (BEHInteg%behavESVA%lGeomInESVA) then
                call rcvalb(fami, kpg, ksp, '+', jvMaterCode, &
                            ' ', relaComp, 0, ' ', [0.d0], &
                            nbProp, propName, propVale, propCode, 1)
            else
                call rcvalb(fami, kpg, ksp, '+', jvMaterCode, &
                            ' ', relaComp, nbPara, paraName, paraVale, &
                            nbProp, propName, propVale, propCode, 1)
            end if
        end if
    end if

end subroutine
