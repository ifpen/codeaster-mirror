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
subroutine calcPrepDataMeca(modelZ, materFieldZ, matecoZ, caraElemZ, &
                            disp_prev, disp_cumu_inst, vari_prev, sigm_prev, &
                            timePrev, timeCurr, &
                            ds_constitutive, ds_material, ds_system, &
                            hval_incr, hval_algo, &
                            vediri, vefnod, &
                            vevarc_prev, vevarc_curr)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nmch1p.h"
#include "asterfort/nmch2p.h"
#include "asterfort/dismoi.h"
#include "asterfort/nmcha0.h"
#include "asterfort/nmchex.h"
#include "asterfort/mvnume.h"
#include "asterfort/chpver.h"
#include "asterfort/nmvcle.h"
#include "asterfort/jeexin.h"
#include "asterfort/vrcomp.h"
#include "asterfort/vrcom2.h"
#include "asterfort/utmess.h"
#include "asterfort/gcncon.h"
#include "asterfort/nmvcre.h"
#include "asterfort/sgcomp.h"
!
    character(len=*), intent(in) :: modelZ, materFieldZ, matecoZ, caraElemZ
    character(len=19), intent(in) :: disp_prev, disp_cumu_inst
    character(len=19), intent(in) :: vari_prev, sigm_prev
    real(kind=8), intent(in) :: timePrev, timeCurr
    type(NL_DS_Constitutive), intent(inout) :: ds_constitutive
    type(NL_DS_Material), intent(out) :: ds_material
    type(NL_DS_System), intent(out) :: ds_system
    character(len=19), intent(out) :: hval_incr(:), hval_algo(:)
    character(len=19), intent(out) :: vediri, vefnod
    character(len=19), intent(out) :: vevarc_prev, vevarc_curr
!
! --------------------------------------------------------------------------------------------------
!
! Command CALCUL
!
! Prepare data for mechanics
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  materField       : name of material characteristics (field)
! In  caraElem         : name of elementary characteristics (field)
! In  disp_prev        : displacement at beginning of current step
! In  disp_cumu_inst   : displacement increment from beginning of step
! In  vari_prev        : internal variables at beginning of step
! In  sigm_prev        : stress at beginning of step
! In  timePrev        : time at beginning of step
! In  timeCurr        : time at end of step
! IO  ds_constitutive  : datastructure for constitutive laws management
! Out ds_material      : datastructure for material parameters
! Out ds_system        : datastructure for non-linear system management
! Out hval_incr        : hat-variable for incremental values fields
! Out hval_algo        : hat-variable for algorithms fields
! Out vediri           : name of elementary for reaction (Lagrange) vector
! Out vefnod           : name of elementary for forces vector (FORC_NODA)
! Out vevarc_prev      : name of elementary for external state variables at beginning of step
! Out vevarc_curr      : name of elementary for external state variables at end of step
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=19) :: disp_curr, varcPrev, varcCurr, sigm_curr, vari_curr, modelLigrel
    character(len=19) :: merigi, veinte
    character(len=24), parameter :: varcRefe = '&&OP0026.COMREF'
    aster_logical :: lModiVari
!
! --------------------------------------------------------------------------------------------------
!

! - Create "hat-variables"
    call nmch1p(hval_incr)
    call nmch2p(hval_algo)

! - Get LIGREL
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)

! - Put displacements in "hat-variables"
    call nmcha0('VALINC', 'DEPMOI', disp_prev, hval_incr)
    call nmcha0('SOLALG', 'DEPDEL', disp_cumu_inst, hval_algo)
!
! - Compute current displacements
!
    if (disp_prev .ne. ' ') then
        call nmchex(hval_incr, 'VALINC', 'DEPPLU', disp_curr)
        call mvnume(disp_prev, disp_cumu_inst, disp_curr)
    end if
!
! - Check and put stress in "hat-variables"
!
    if (sigm_prev .ne. ' ') then
        call sgcomp(ds_constitutive%compor, sigm_prev, modelLigrel, iret)
        if (iret .eq. 1) then
            call utmess('F', 'CALCUL1_6')
        end if
        call chpver('F', sigm_prev, 'ELGA', 'SIEF_R', iret)
        call nmcha0('VALINC', 'SIGMOI', sigm_prev, hval_incr)
    end if
!
! - Check and put internal variables in "hat-variables"
!
    if (vari_prev .ne. ' ') then
        call chpver('F', vari_prev, 'ELGA', 'VARI_R', iret)
        call nmcha0('VALINC', 'VARMOI', vari_prev, hval_incr)
    end if

! - Get name of fields for external state variables
    call nmchex(hval_incr, 'VALINC', 'COMMOI', varcPrev)
    call nmchex(hval_incr, 'VALINC', 'COMPLU', varcCurr)

! - Create fields for external state variables
    call nmvcle(modelZ, materFieldZ, caraElemZ, timeCurr, varcCurr)
    call nmvcle(modelZ, materFieldZ, caraElemZ, timePrev, varcPrev)
    call nmvcre(modelZ, materFieldZ, caraElemZ, varcRefe)

! - Checking number of internal variables
    call jeexin(ds_constitutive%compor(1:19)//'.CESD', iret)
    if (iret .gt. 0 .and. vari_prev .ne. ' ') then
        call vrcomp(ds_constitutive%compor, vari_prev, modelLigrel, iret, lModiVari_=lModiVari)
        if (iret .eq. 1) then
            call utmess('F', 'CALCUL1_5')
        end if
        if (lModiVari) then
            call vrcom2(ds_constitutive%compor, vari_prev, modelLigrel, ASTER_FALSE)
        end if
    end if

! - Datastructures name (automatic génération)
    call gcncon('_', sigm_curr)
    call gcncon('_', vari_curr)
    call gcncon('_', merigi)
    call gcncon('_', veinte)
    call gcncon('_', vediri)
    call gcncon('_', vefnod)
    call gcncon('_', vevarc_prev)
    call gcncon('_', vevarc_curr)
    call gcncon('_', ds_constitutive%comp_error)

! - Changeing names of variables
    call nmcha0('VALINC', 'SIGPLU', sigm_curr, hval_incr)
    call nmcha0('VALINC', 'VARPLU', vari_curr, hval_incr)

! - Prepare datastructures
    ds_material%mater = materFieldZ
    ds_material%mateco = matecoZ
    ds_material%varc_refe = varcRefe
    ds_system%merigi = merigi
    ds_system%veinte = veinte
!
end subroutine
