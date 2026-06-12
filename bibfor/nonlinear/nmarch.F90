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
subroutine nmarch(numeInst, &
                  model, ds_material, caraElem, listFuncActi, &
                  ds_print, sddisc, sdcrit, &
                  ds_measure, sderro, sddyna, sdpilo, ds_energy, &
                  ds_inout, ds_errorindic, ds_algorom_, lStoreInitState_)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
    use HHO_postpro_module, only: hhoPostDeplMeca
    implicit none
!
#include "asterf_types.h"
#include "asterfort/diinst.h"
#include "asterfort/dinuar.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmarc0.h"
#include "asterfort/nmarce.h"
#include "asterfort/nmarpc.h"
#include "asterfort/nmcrpc.h"
#include "asterfort/nmfinp.h"
#include "asterfort/nmleeb.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
#include "asterfort/romAlgoNLTableSave.h"
#include "asterfort/rsagsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/storeSaveLast.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpg.h"
!
    integer(kind=8), intent(in) :: numeInst
    character(len=24), intent(in) :: model
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=24), intent(in) :: caraElem
    integer(kind=8), intent(in) :: listFuncActi(*)
    type(NL_DS_Print), intent(in) :: ds_print
    character(len=19), intent(in) :: sddisc, sdcrit
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=24), intent(in) :: sderro
    character(len=19), intent(in) :: sddyna, sdpilo
    type(NL_DS_Energy), intent(in) :: ds_energy
    type(NL_DS_InOut), intent(in) :: ds_inout
    type(NL_DS_ErrorIndic), intent(in) :: ds_errorindic
    type(ROM_DS_AlgoPara), optional, intent(in) :: ds_algorom_
    aster_logical, intent(in), optional :: lStoreInitState_
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Storing results
!
! --------------------------------------------------------------------------------------------------
!
! In  numeInst         : index of current time step
! In  model            : name of model
! In  ds_material      : datastructure for material parameters
! In  caraElem         : name of elementary characteristics (field)
! In  listFuncActi     : list of active functionnalities
! In  ds_print         : datastructure for printing parameters
! In  sddisc           : datastructure for time discretization
! In  sdcrit           : name of datastructure to save convergence parameters
! IO  ds_measure       : datastructure for measure and statistics management
! In  sderro           : name of datastructure for error management (events)
! In  sddyna           : name of datastructure for dynamic parameters
! In  sdpilo           : continuation ("PILOTAGE") parameters datastructure
! In  ds_energy        : datastructure for energy management
! In  ds_inout         : datastructure for input/output management
! In  ds_errorindic    : datastructure for error indicator
! In  ds_algorom       : datastructure for ROM parameters
! In  lStoreInitState  : flag to store initial state
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, numeStore
    real(kind=8) :: timeCurr
    character(len=8) :: result
    aster_logical :: lForceStore, lprint, l_hho, lStoreInitState, lastTimeStep
    character(len=19) :: k19bid, sdarch
    character(len=24) :: listLoadResu
    character(len=4) :: etcalc
    integer(kind=8) :: numeReuse
!
! --------------------------------------------------------------------------------------------------
!
    lStoreInitState = ASTER_FALSE
    if (present(lStoreInitState_)) then
        lStoreInitState = lStoreInitState_
    end if

! - Initializations
    result = ds_inout%result
    listLoadResu = ds_inout%listLoadResu
    l_hho = isfonc(listFuncActi, 'HHO')
    sdarch = sddisc(1:14)//'.ARCH'

! - Loop state
    call nmleeb(sderro, 'CALC', etcalc)

! - Last step => storing
    lForceStore = ASTER_FALSE
    call nmfinp(sddisc, numeInst, lastTimeStep)
    lForceStore = lastTimeStep

! - Storing
    if (etcalc .eq. 'CONV' .or. etcalc .eq. 'STOP') then
        lForceStore = ASTER_TRUE
    end if

! - Print timer
    call uttcpg('IMPR', 'INCR')

! - Get index for storing
    call dinuar(result, sddisc, numeInst, lForceStore, &
                numeStore, numeReuse, lStoreInitState)

! - Current time
    timeCurr = diinst(sddisc, numeInst)

! - Save energy parameters in output table
    if (isfonc(listFuncActi, 'ENERGIE')) then
        call nmarpc(ds_energy, numeReuse, timeCurr)
    else
        call nmcrpc(ds_inout, numeReuse, timeCurr)
    end if

! - Print or not ?
    lprint = ds_print%l_print
    call storeSaveLast(sdarch, numeStore, timeCurr)

! - Storing
    if (numeStore .ge. 0) then
! ----- Begin timer
        call nmtime(ds_measure, 'Launch', 'Store')

! ----- Print head
        if (lprint) then
            call utmess('I', 'ARCHIVAGE_5')
        end if

! ----- Increased result datastructure if necessary
        call rsexch(' ', result, 'DEPL', numeStore, k19bid, iret)
        if (iret .eq. 110) then
            call rsagsd(result, 0)
        end if

! ----- Storing parameters
        call nmarc0(result, model, ds_material, caraElem, listFuncActi, &
                    sdcrit, sddyna, ds_errorindic, &
                    sdpilo, listLoadResu, numeStore, timeCurr)

! ----- Storing fields
        call nmarce(ds_inout, result, sddisc, timeCurr, numeStore, &
                    lForceStore, ds_print)

! ----- If HHO, we compute a post_processing
        if (l_hho) then
            call hhoPostDeplMeca(model, result, numeStore)
        end if

! ----- Storing reduced parameters table (ROM)
        if (present(ds_algorom_)) then
            if (ds_algorom_%l_rom) then
                if (numeStore .gt. 0) then
                    call romAlgoNLTableSave(numeStore, timeCurr, ds_algorom_)
                end if
            end if
        end if

! ----- End timer
        call nmtime(ds_measure, 'Stop', 'Store')
        call nmrinc(ds_measure, 'Store')
    end if
!
end subroutine
