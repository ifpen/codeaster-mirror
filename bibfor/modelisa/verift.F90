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

subroutine verift(fami, kpg, ksp, poum, jvMaterCode, &
                  materi_, iret_, epsth_, epsth_anis_, epsth_meta_, &
                  temp_prev_, temp_curr_, temp_refe_)
!
    use MetallurgyMeca_module, only: metaGetType
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elasth_para.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp
    character(len=*), intent(in) :: poum
    integer(kind=8), intent(in) :: jvMaterCode
    character(len=8), optional, intent(in) :: materi_
    integer(kind=8), optional, intent(out) :: iret_
    real(kind=8), optional, intent(out) :: epsth_, epsth_anis_(3), epsth_meta_
    real(kind=8), optional, intent(out) :: temp_prev_, temp_curr_, temp_refe_
!
! --------------------------------------------------------------------------------------------------
!
! Compute thermic dilatation
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  jvMaterCode  : coded material address
! In  poum         : parameters evaluation
!                     '-' for previous temperature
!                     '+' for current temperature
!                     'T' for current and previous temperature => epsth is increment
! In  kpg          : current point gauss
! In  ksp          : current "sous-point" gauss
! In  materi       : name of material if multi-material Gauss point (PMF)
! In  iret         : 0 if temperature defined
!                    1 if not
! Out epsth_       : thermic dilatation
! Out epsth_anis_  : non-isotropic thermic dilatation
! Out epsth_meta_  : multiphasic thermic dilatation
! Out temp_prev    : previous temperature
! Out temp_curr    : current temperature
! Out temp_refe    : reference temperature
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elem_name, materi
    integer(kind=8) :: iret_temp_prev, iret_temp_curr, iret_temp, iret_temp_refe
    real(kind=8) :: temp_prev, temp_refe, temp_curr
    real(kind=8) :: epsth, epsth_anis(3), epsth_meta
    real(kind=8) :: alpha_p(2)
    real(kind=8) :: alpha_l_p, alpha_t_p, alpha_n_p
    real(kind=8) :: alpha_c(2)
    real(kind=8) :: alpha_l_c, alpha_t_c, alpha_n_c
    real(kind=8) :: zSumColdPrev, zHotPrev, zSumColdCurr, zHotCurr, epsth_meta_h, epsth_meta_c
    real(kind=8) :: z_h_r, deps_ch_tref
    integer(kind=8) :: iadzi, iazk24
    integer(kind=8) :: elasID, iret_cmp, jvCompor, metaType, nbPhase
    character(len=16) :: elasKeyword, relaComp
!
! --------------------------------------------------------------------------------------------------
!
    materi = ' '
    if (present(materi_)) then
        materi = materi_
    end if
!
    iret_temp = 0
    iret_temp_prev = 0
    iret_temp_curr = 0
    iret_temp_refe = 0
    epsth = 0.d0
    epsth_anis = 0.d0
    epsth_meta = 0.d0
    temp_refe = 0.d0
    temp_curr = 0.d0
    temp_prev = 0.d0

! - No temperature -> thermic strain is zero
    call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                ksp, temp_curr, iret_temp)
    if (iret_temp .ne. 0) then
        goto 999
    end if

! - Get reference temperature
    call rcvarc(' ', 'TEMP', 'REF', fami, kpg, &
                ksp, temp_refe, iret_temp_refe)
    if (iret_temp_refe .eq. 1) then
        call tecael(iadzi, iazk24)
        elem_name = zk24(iazk24-1+3) (1:8)
        call utmess('F', 'COMPOR5_8', sk=elem_name)
    end if

! - Get type of elasticity (Isotropic/Orthotropic/Transverse isotropic)
    call get_elas_id(jvMaterCode, elasID, elasKeyword)

! - Get temperatures
    if (poum .eq. 'T' .or. poum .eq. '-') then
        call rcvarc(' ', 'TEMP', '-', fami, kpg, &
                    ksp, temp_prev, iret_temp_prev)
    end if
    if (poum .eq. 'T' .or. poum .eq. '+') then
        call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                    ksp, temp_curr, iret_temp_curr)
    end if

! - Get elastic parameters for thermic dilatation
    if (poum .eq. 'T' .or. poum .eq. '-') then
        if (iret_temp_prev .eq. 0) then
            call get_elasth_para(fami, jvMaterCode, '-', kpg, ksp, &
                                 elasID, elasKeyword, materi_=materi, &
                                 alpha=alpha_p, &
                                 alpha_l=alpha_l_p, &
                                 alpha_t=alpha_t_p, &
                                 alpha_n=alpha_n_p)
        end if
    end if
    if (poum .eq. 'T' .or. poum .eq. '+') then
        if (iret_temp_curr .eq. 0) then
            call get_elasth_para(fami, jvMaterCode, '+', kpg, ksp, &
                                 elasID, elasKeyword, materi_=materi, &
                                 alpha=alpha_c, &
                                 alpha_l=alpha_l_c, &
                                 alpha_t=alpha_t_c, &
                                 alpha_n=alpha_n_c)
        end if
    end if

! - Check non-isotropic material
    if (elasID .ne. ELAS_ISOT .and. elasID .ne. ELAS_SHELL .and. &
        elasID .ne. ELAS_MEMBRANE) then
        if (.not. present(epsth_anis_)) then
            call tecael(iadzi, iazk24)
            elem_name = zk24(iazk24-1+3) (1:8)
            call utmess('F', 'COMPOR5_9', sk=elem_name)
        end if
    end if

! - Check metallurgical material
    if (elasKeyword .eq. 'ELAS_META') then
        zSumColdPrev = 0.d0
        zHotPrev = 0.d0
        zSumColdCurr = 0.d0
        zHotCurr = 0.d0

! ----- Check behavior name (used only in metallurgical case)
        call tecach('NNO', 'PCOMPOR', 'L', iret_cmp, iad=jvCompor)
        if (iret_cmp .eq. 0) then
            relaComp = zk16(jvCompor)
        else
            relaComp = 'Unknown'
        end if

        if (.not. present(epsth_meta_)) then
            call tecael(iadzi, iazk24)
            elem_name = zk24(iazk24-1+3) (1:8)
            call utmess('F', 'COMPOR5_10', sk=elem_name)
        end if
!
        call metaGetType(metaType, nbPhase)
        if (nbPhase .eq. 0) then
            ASSERT(relaComp .eq. 'META_LEMA_ANI')
        else
            if (poum .eq. 'T' .or. poum .eq. '-') then
                if (iret_temp_prev .eq. 0) then
                    call metaGetPhase(fami, '-', kpg, ksp, &
                                      metaType, nbPhase, &
                                      zSumCold_=zSumColdPrev, &
                                      zHot_=zHotPrev)
                end if
            end if
            if (poum .eq. 'T' .or. poum .eq. '+') then
                if (iret_temp_prev .eq. 0) then
                    call metaGetPhase(fami, '+', kpg, ksp, &
                                      metaType, nbPhase, &
                                      zSumCold_=zSumColdCurr, &
                                      zHot_=zHotCurr)
                end if
            end if
        end if

    end if

! - Compute thermic strain
    if (poum .eq. 'T') then
        if (iret_temp_prev+iret_temp_curr .eq. 0) then
            if (elasID .eq. ELAS_ISOT .or. elasID .eq. ELAS_SHELL .or. &
                elasID .eq. ELAS_MEMBRANE) then
                if (elasKeyword .eq. 'ELAS_META') then
                    if (relaComp .eq. 'META_LEMA_ANI') then
                        epsth_meta = 0.d0
                    else
                        epsth_meta_h = zHotCurr*alpha_c(1)*(temp_curr-temp_refe)- &
                                       zHotPrev*alpha_p(1)*(temp_prev-temp_refe)
                        epsth_meta_c = zSumColdCurr*alpha_c(2)*(temp_curr-temp_refe)- &
                                       zSumColdPrev*alpha_p(2)*(temp_prev-temp_refe)
                        call get_elasth_para(fami, jvMaterCode, '+', kpg, ksp, &
                                             elasID, elasKeyword, materi_=materi, &
                                             z_h_r_=z_h_r, deps_ch_tref_=deps_ch_tref)
                        epsth_meta_h = epsth_meta_h+(1-z_h_r)*deps_ch_tref*(zHotPrev-zHotCurr)
                        epsth_meta_c = epsth_meta_c+z_h_r*deps_ch_tref*(zSumColdCurr-zSumColdPrev)
                        epsth_meta = epsth_meta_h+epsth_meta_c
                    end if

                else
                    epsth = alpha_c(1)*(temp_curr-temp_refe)-alpha_p(1)*(temp_prev-temp_refe)
                end if
            elseif (elasID .eq. ELAS_ORTH) then
                epsth_anis(1) = alpha_l_c*(temp_curr-temp_refe)-alpha_l_p*(temp_prev-temp_refe)
                epsth_anis(2) = alpha_t_c*(temp_curr-temp_refe)-alpha_t_p*(temp_prev-temp_refe)
                epsth_anis(3) = alpha_n_c*(temp_curr-temp_refe)-alpha_n_p*(temp_prev-temp_refe)

            elseif (elasID .eq. ELAS_ISTR) then
                epsth_anis(1) = alpha_l_c*(temp_curr-temp_refe)-alpha_l_p*(temp_prev-temp_refe)
                epsth_anis(2) = alpha_n_c*(temp_curr-temp_refe)-alpha_n_p*(temp_prev-temp_refe)

            else
                ASSERT(ASTER_FALSE)
            end if

        end if

    else if (poum .eq. '-') then
        if (iret_temp_prev .eq. 0) then
            if (elasID .eq. ELAS_ISOT .or. elasID .eq. ELAS_SHELL .or. &
                elasID .eq. ELAS_MEMBRANE) then
                if (elasKeyword .eq. 'ELAS_META') then
                    if (relaComp .eq. 'META_LEMA_ANI') then
                        epsth_meta = 0.d0
                    else
                        epsth_meta_h = zHotPrev*alpha_p(1)*(temp_prev-temp_refe)
                        epsth_meta_c = zSumColdPrev*alpha_p(2)*(temp_prev-temp_refe)
                        call get_elasth_para(fami, jvMaterCode, '+', kpg, ksp, &
                                             elasID, elasKeyword, materi_=materi, &
                                             z_h_r_=z_h_r, deps_ch_tref_=deps_ch_tref)
                        epsth_meta_h = epsth_meta_h-zHotPrev*(1-z_h_r)*deps_ch_tref
                        epsth_meta_c = epsth_meta_c+zSumColdPrev*z_h_r*deps_ch_tref
                        epsth_meta = epsth_meta_h+epsth_meta_c
                    end if
                else
                    epsth = alpha_p(1)*(temp_prev-temp_refe)
                end if

            elseif (elasID .eq. ELAS_ORTH) then
                epsth_anis(1) = alpha_l_p*(temp_prev-temp_refe)
                epsth_anis(2) = alpha_t_p*(temp_prev-temp_refe)
                epsth_anis(3) = alpha_n_p*(temp_prev-temp_refe)

            elseif (elasID .eq. ELAS_ISTR) then
                epsth_anis(1) = alpha_l_p*(temp_prev-temp_refe)
                epsth_anis(2) = alpha_n_p*(temp_prev-temp_refe)

            else
                ASSERT(ASTER_FALSE)
            end if
        end if

    else if (poum .eq. '+') then
        if (iret_temp_curr .eq. 0) then
            if (elasID .eq. ELAS_ISOT .or. elasID .eq. ELAS_SHELL .or. &
                elasID .eq. ELAS_MEMBRANE) then
                if (elasKeyword .eq. 'ELAS_META') then
                    if (relaComp .eq. 'META_LEMA_ANI') then
                        epsth_meta = 0.d0
                    else
                        epsth_meta_h = zHotCurr*alpha_c(1)*(temp_curr-temp_refe)
                        epsth_meta_c = zSumColdCurr*alpha_c(2)*(temp_curr-temp_refe)
                        call get_elasth_para(fami, jvMaterCode, '+', kpg, ksp, &
                                             elasID, elasKeyword, materi_=materi, &
                                             z_h_r_=z_h_r, deps_ch_tref_=deps_ch_tref)
                        epsth_meta_h = epsth_meta_h-zHotCurr*(1-z_h_r)*deps_ch_tref
                        epsth_meta_c = epsth_meta_c+zSumColdCurr*z_h_r*deps_ch_tref
                        epsth_meta = epsth_meta_h+epsth_meta_c
                    end if
                else
                    epsth = alpha_c(1)*(temp_curr-temp_refe)
                end if

            elseif (elasID .eq. ELAS_ORTH) then
                epsth_anis(1) = alpha_l_c*(temp_curr-temp_refe)
                epsth_anis(2) = alpha_t_c*(temp_curr-temp_refe)
                epsth_anis(3) = alpha_n_c*(temp_curr-temp_refe)

            elseif (elasID .eq. ELAS_ISTR) then
                epsth_anis(1) = alpha_l_c*(temp_curr-temp_refe)
                epsth_anis(2) = alpha_n_c*(temp_curr-temp_refe)

            else
                ASSERT(ASTER_FALSE)
            end if
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
!
999 continue

! - Output temperature
    if (present(temp_refe_)) then
        temp_refe_ = temp_refe
    end if
    if (present(temp_prev_)) then
        temp_prev_ = temp_prev
    end if
    if (present(temp_curr_)) then
        temp_curr_ = temp_curr
    end if

! - Output strains
    if (present(epsth_meta_)) then
        epsth_meta_ = epsth_meta
    end if
    if (present(epsth_anis_)) then
        epsth_anis_(1:3) = epsth_anis(1:3)
    end if
    if (present(epsth_)) then
        epsth_ = epsth
    end if

! - Output error
    if (present(iret_)) then
        iret_ = 0
        if ((iret_temp_prev+iret_temp_curr) .ne. 0) then
            iret_ = 1
        end if
        if (iret_temp .ne. 0) then
            iret_ = 1
        end if
    end if
!
end subroutine
