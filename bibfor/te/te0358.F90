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
subroutine te0358(option, nomte)
!
    use MetallurgyMeca_module, only: metaGetType
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/jevech.h"
#include "asterfort/meta_vpta_coef.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
! Option: CHAR_MECA_META_Z
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1, nbSigm = 6
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8) :: iNode, iSigm, lgpg, iret
    real(kind=8) :: sigmo
    character(len=16) :: metaPhasName, relaComp, valk(2)
    integer(kind=8) :: jvMater, jvMaterCode
    integer(kind=8) :: metaType, nbPhases
    real(kind=8) :: young, nu, deuxmu
    integer(kind=8) :: j_sigm
    integer(kind=8) :: elasID
    integer(kind=8) :: j_poids, j_vf, j_dfde, jvGeom
    integer(kind=8) :: nbNode, kpg, npg, jtab(7)
    integer(kind=8) :: j_vectu
    integer(kind=8) :: j_vari
    real(kind=8) :: sig(nbSigm), sigdv(nbSigm)
    real(kind=8) :: dfdx(MT_NNOMAX3D), dfdy(MT_NNOMAX3D), dfdz(MT_NNOMAX3D), poids
    real(kind=8) :: coef, trans
    real(kind=8) :: zSumColdCurr
    real(kind=8) :: phasePrev(META_MECA_NBPHASE_MAXI), phaseCurr(META_MECA_NBPHASE_MAXI), temp
    aster_logical :: l_temp
    character(len=16) :: elasKeyword
    character(len=16) :: metaRela, metaGlob
    character(len=16), pointer :: compor(:) => null()
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
!
! --------------------------------------------------------------------------------------------------
!

! - Element reference
    call elrefe_info(fami=fami, nno=nbNode, npg=npg, &
                     jpoids=j_poids, jvf=j_vf, jdfde=j_dfde)
    call tecach('OOO', 'PVARIPR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)
    relaComp = compor(RELA_NAME)
    metaRela = compor(META_RELA)
    metaGlob = compor(META_GLOB)

! - Get type of phases
    metaPhasName = compor(META_PHAS)
    call metaGetType(metaType, nbPhases)
    ASSERT(nbPhases .le. META_MECA_NBPHASE_MAXI)

    if ((metaType .ne. META_NONE) .and. (relaComp .ne. 'META_LEMA_ANI') .and. &
        metaPhasName .ne. "VIDE") then

        ASSERT(nbPhases .ne. 0)

! ----- Check type of phases
        valk(1) = metaPhasName
        if (metaPhasName .eq. 'ACIERMECA') then
            if (metaType .ne. META_STEEL) then
                valk(2) = 'ZIRC'
                call utmess('F', 'COMPOR3_8', nk=2, valk=valk)
            end if
        elseif (metaPhasName .eq. 'ZIRCMECA') then
            if (metaType .ne. META_ZIRC) then
                valk(2) = 'ACIER'
                call utmess('F', 'COMPOR3_8', nk=2, valk=valk)
            end if
        elseif (metaPhasName .eq. 'ACIER_REVENUMECA') then
            if (metaType .ne. META_STEEL_TEMPER) then
                valk(2) = "ACIER_REVENU"
                call utmess('F', 'COMPOR3_8', nk=2, valk=valk)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Stresses
        call jevech('PCONTMR', 'L', j_sigm)

! ----- Internal variables
        call jevech('PVARIPR', 'L', j_vari)

! ----- Material parameters
        call jevech('PMATERC', 'L', jvMater)
        jvMaterCode = zi(jvMater)

! ----- Output vector
        call jevech('PVECTUR', 'E', j_vectu)
        do iNode = 1, nbNode
            zr(j_vectu+3*iNode-3) = 0.d0
            zr(j_vectu+3*iNode-2) = 0.d0
            zr(j_vectu+3*iNode-1) = 0.d0
        end do
!
        do kpg = 1, npg
! --------- Derived of shape functions
            call dfdm3d(nbNode, kpg, j_poids, j_dfde, zr(jvGeom), &
                        poids, dfdx, dfdy, dfdz)

! --------- Get current temperature
            call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                        ksp, temp, iret)
            l_temp = iret .eq. 0

! --------- Get phases
            phasePrev = 0.d0
            phaseCurr = 0.d0
            call metaGetPhase(fami, '-', kpg, ksp, metaType, &
                              nbPhases, phasePrev)
            call metaGetPhase(fami, '+', kpg, ksp, metaType, &
                              nbPhases, phaseCurr, zSumCold_=zSumColdCurr)

! --------- Get elastic parameters
            call get_elas_id(jvMaterCode, elasID, elasKeyword)
            call get_elas_para(fami, jvMaterCode, '+', kpg, ksp, &
                               elasID, elasKeyword, &
                               e_=young, nu_=nu)
            ASSERT(elasID .eq. ELAS_ISOT)
            deuxmu = young/(1.d0+nu)

! --------- Compute coefficients for second member
            call meta_vpta_coef(metaRela, metaGlob, &
                                lgpg, fami, kpg, jvMaterCode, &
                                l_temp, temp, metaType, nbPhases, phasePrev, &
                                phaseCurr, zSumColdCurr, young, deuxmu, coef, &
                                trans)

! --------- Compute stresses
            sigmo = 0.d0
            do iNode = 1, 3
                sigmo = sigmo+zr(j_sigm+(kpg-1)*nbSigm+iNode-1)
            end do
            sigmo = sigmo/3.d0
            do iSigm = 1, nbSigm
                sigdv(iSigm) = zr(j_sigm+(kpg-1)*nbSigm+iSigm-1)-sigmo*kron(iSigm)
                sig(iSigm) = coef*(1.5d0*trans*sigdv(iSigm))
                sig(iSigm) = deuxmu*sig(iSigm)
            end do

! --------- Second member
            do iNode = 1, nbNode
                zr(j_vectu+3*(iNode-1)) = zr(j_vectu+3*(iNode-1))+ &
                                          poids*(sig(1)*dfdx(iNode)+ &
                                                 sig(4)*dfdy(iNode)+ &
                                                 sig(5)*dfdz(iNode))
                zr(j_vectu+3*(iNode-1)+1) = zr(j_vectu+3*(iNode-1)+1)+ &
                                            poids*(sig(2)*dfdy(iNode)+ &
                                                   sig(4)*dfdx(iNode)+ &
                                                   sig(6)*dfdz(iNode))
                zr(j_vectu+3*(iNode-1)+2) = zr(j_vectu+3*(iNode-1)+2)+ &
                                            poids*(sig(3)*dfdz(iNode)+ &
                                                   sig(5)*dfdx(iNode)+ &
                                                   sig(6)*dfdy(iNode))
            end do
        end do
    end if
end subroutine
