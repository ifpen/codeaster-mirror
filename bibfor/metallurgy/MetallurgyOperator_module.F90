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
module MetallurgyOperator_module
! ==================================================================================================
    use Metallurgy_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: metaGetParameters, metaDelParameters
    public :: metaCompPhases, metaCompOtherOptions
    public :: metaGetInitialState, metaPrepTRCWorkingField
    public :: metaCompInitialField, metaPrepInitialState, metaCompMetaElno
    public :: metaComp, metaTemper
    private :: metaCompTransfField
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcop.h"
#include "asterfort/calcul.h"
#include "asterfort/cesvar.h"
#include "asterfort/chpver.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlima.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeveut.h"
#include "asterfort/mecact.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/mtdorc.h"
#include "asterfort/rcadme.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rs_get_liststore.h"
#include "asterfort/rs_getnume.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rslesd.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! metaGetInitialState
!
! Get field for initial phases
!
! In  resultName       : name of result datastructure
! Out metaInit         : name of field for initial phases
! Out numeFieldInit    : storing index of initial field in result datastructure
!                        (0 if already exists)
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetInitialState(resultName, metaInit, numeFieldInit)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: resultName
        character(len=24), intent(out) :: metaInit
        integer(kind=8), intent(out) :: numeFieldInit
! ----- Local
        character(len=8) :: resultInit
        character(len=24) :: selectCrit, metaFieldUser
        character(len=16), parameter :: factorKeyword = "ETAT_INIT"
        integer(kind=8) :: nbRet, iret, nbtrou, numeStoreInit
        real(kind=8) :: timeInit, selectTole
!   ------------------------------------------------------------------------------------------------
!
        metaInit = "Unknown"
        numeFieldInit = 0
        call getvid(factorKeyword, 'META_INIT_ELNO', iocc=1, scal=metaFieldUser, nbret=nbRet)
        if (nbRet .gt. 0) then
! --------- Check type
            call chpver('F', metaFieldUser, 'CART', 'VAR2_R', iret)
            metaInit = "&&SMEVOL_ZINIT"
            call copisd('CHAMP_GD', 'V', metaFieldUser, metaInit)
        else
            call getvid(factorKeyword, 'EVOL_THER', iocc=1, scal=resultInit, nbret=nbRet)
            if (resultInit .ne. resultName) then
                call utmess('F', 'METALLURGY1_11')
            end if
            call getvis(factorKeyword, 'NUME_INIT', iocc=1, scal=numeStoreInit, nbret=nbRet)
            if (nbRet .eq. 0) then
                call getvr8(factorKeyword, 'INST_INIT', iocc=1, scal=timeInit, nbret=nbRet)
                call getvr8(factorKeyword, 'PRECISION', iocc=1, scal=selectTole, nbret=nbRet)
                call getvtx(factorKeyword, 'CRITERE', iocc=1, scal=selectCrit, nbret=nbRet)
                call rs_getnume(resultName, timeInit, selectCrit, selectTole, numeStoreInit, nbtrou)
                if (nbtrou .eq. 0 .or. nbtrou .gt. 1) then
                    call utmess('F', 'METALLURGY1_51', sr=timeInit)
                end if
            end if
            call rsexch('F', resultName, 'META_ELNO', numeStoreInit, metaInit, iret)
            if (iret .ne. 0) then
                call utmess('F', 'METALLURGY1_52')
            end if
            numeFieldInit = numeStoreInit
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaPrepTRCWorkingField
!
! Prepare field to manage TRC curves in elementary computation
!
! IO  metaParaOperator : datastructure for parameters of CALC_META operator
!
! --------------------------------------------------------------------------------------------------
    subroutine metaPrepTRCWorkingField(metaParaOperator)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(inout) :: metaParaOperator
! ----- Local
        integer(kind=8) :: iMaterVale, materValeLen
        character(len=8), pointer :: materVale(:) => null()
        integer(kind=8) :: nbhist, icodre
        character(len=8) :: materPara
        integer(kind=8), parameter :: nbPara = 2
        character(len=8), parameter :: paraName(nbPara) = (/'I1  ', 'I2  '/)
        integer(kind=8) :: iadtrc(nbPara), adrsJv(nbPara)
!   ------------------------------------------------------------------------------------------------
!

! ----- Get access to material field
        call jeveuo(metaParaOperator%materialField(1:8)//'.CHAMP_MAT .VALE', 'E', vk8=materVale)
        call jelira(metaParaOperator%materialField(1:8)//'.CHAMP_MAT .VALE', 'LONMAX', materValeLen)

! ----- Get number of histories in TRC diagram
        nbhist = 0
        metaParaOperator%hasTRC = ASTER_FALSE
        iadtrc = 0
        do iMaterVale = 1, materValeLen
            materPara = materVale(iMaterVale)
            if (materPara .ne. '        ') then
                call rcadme(materPara, 'META_ACIER', 'TRC', iadtrc, icodre, 0)
                if (icodre .eq. 0) metaParaOperator%hasTRC = ASTER_TRUE
                nbhist = max(nbhist, iadtrc(1))
            end if
        end do

! ----- Create working field for TRC
        if (metaParaOperator%hasTRC) then
            call wkvect('&&SMEVOL_FTRC', 'V V R', 9*nbhist, adrsJv(1))
            call wkvect('&&SMEVOL_TRC', 'V V R', 15*nbhist, adrsJv(2))
            call jeveut('&&SMEVOL_FTRC', "E", adrsJv(1))
            call jeveut('&&SMEVOL_TRC', "E", adrsJv(2))
            call mecact('V', metaParaOperator%TRCField, 'LIGREL', &
                        metaParaOperator%modelLigrel, 'ADRSJEVN', &
                        ncmp=nbPara, lnomcmp=paraName, vi=adrsJv)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaCompInitialField
!
! Compute initial field of metallurgy
!
! From the map describing the initial phases provided by the user, we create the complete map
! with the other necessary parameters.
!
! In  metaParaOperator : datastructure for parameters of CALC_META operator
! In  temp             : temperature
! In  metaInitUser     : map describing the initial phases provided by the user
! In  metaInit         : name of map for initial metallurgy computation
!
! --------------------------------------------------------------------------------------------------
    subroutine metaCompInitialField(metaParaOperator, temp, &
                                    metaInitUser, metaInit)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(in) :: metaParaOperator
        character(len=24), intent(in) :: temp, metaInitUser
        character(len=24), intent(in) :: metaInit
! ----- Local
        integer(kind=8), parameter :: nbout = 1, nbInMax = 4
        character(len=8) :: lpaout(nbout), lpain(nbInMax)
        character(len=24) :: lchin(nbInMax)
        character(len=1), parameter :: base = "V"
        character(len=16), parameter :: option = "META_INIT_ELNO"
!   ------------------------------------------------------------------------------------------------
!
        lpain(1) = 'PMATERC'
        lchin(1) = metaParaOperator%materialCoding
        lpain(2) = 'PCOMPME'
        lchin(2) = metaParaOperator%comporMeta
        lpain(3) = 'PTEMPER'
        lchin(3) = temp
        lpain(4) = 'PPHASII'
        lchin(4) = metaInitUser
        lpaout(1) = 'PPHASOUT'
        call calcul('S', option, metaParaOperator%modelLigrel, nbInMax, lchin, &
                    lpain, nbout, metaInit, lpaout, base, &
                    'OUI')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaPrepInitialState
!
! Prepare initial state
!
! In  metaParaOperator : datastructure for parameters of CALC_META operator
! In  numeStore        : current index of storing to proceed
! In  metaInitUser     : map describing the initial phases provided by the user
!
! --------------------------------------------------------------------------------------------------
    subroutine metaPrepInitialState(metaParaOperator, numeStore, metaInitUser)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(in) :: metaParaOperator
        integer(kind=8), intent(in) :: numeStore
        character(len=24), intent(in) :: metaInitUser
! ----- Local
        integer(kind=8) :: iret, jvPara
        character(len=24), parameter :: metaInit = '&&SMEVOL.PHAS_META1'
        character(len=24) :: temp, resultField
        real(kind=8) :: time
!   ------------------------------------------------------------------------------------------------
!

! ----- Prepare metallurgical field (dynamic field)
        call copisd('CHAM_ELEM_S', 'V', metaParaOperator%comporMeta, metaInit)

! ----- Get temperature field
        call rsexch('F', metaParaOperator%resultName, 'TEMP', numeStore, temp, iret)

! ----- Get current time
        call rsadpa(metaParaOperator%resultName, 'L', 1, 'INST', numeStore, 0, sjv=jvPara)
        time = zr(jvPara)

! ----- Compute initial field of metallurgy
        call metaCompInitialField(metaParaOperator, temp, metaInitUser, metaInit)

! ----- Save map for initial metallurgy computation in result datastructure
        call rsexch(' ', metaParaOperator%resultName, 'META_ELNO', numeStore, resultField, iret)
        call copisd('CHAMP_GD', 'G', metaInit, resultField)
        call rsnoch(metaParaOperator%resultName, 'META_ELNO', numeStore)
        call utmess('I', 'ARCHIVAGE_6', sk='META_ELNO', si=numeStore, sr=time)

! ----- Save Save COMPORMETA in result datastructure
        call rsexch(' ', metaParaOperator%resultName, 'COMPORMETA', numeStore, resultField, iret)
        call copisd('CHAMP_GD', 'G', metaParaOperator%comporMeta, resultField)
        call rsnoch(metaParaOperator%resultName, 'COMPORMETA', numeStore)
        call utmess('I', 'ARCHIVAGE_6', sk='COMPORMETA', si=numeStore, sr=time)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaCompMetaElno
!
! Compute option META_ELNO
!
! In  metaParaOperator : datastructure for parameters of CALC_META operator
! In  time_1           : time at start of time step
! In  time_2           : time at end of time step
! In  temp_1           : temperature at start of time step
! In  temp_2           : temperature at end of time step
! In  metaIn           : input field of metallurgy
! Out metaOut          : output field of metallurgy
! In  time_0           : previous time (to estimate speed of temperature)
! In  temp_0           : previous temperature (to estimate speed of temperature)
!
! --------------------------------------------------------------------------------------------------
    subroutine metaCompMetaElno(metaParaOperator, &
                                time_1, time_2, &
                                temp_1, temp_2, &
                                metaIn, metaOut, &
                                forTemper, &
                                metaPrev_, &
                                time_0_, temp_0_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(in) :: metaParaOperator
        real(kind=8), intent(in) :: time_1, time_2
        character(len=24), intent(in) :: temp_1, temp_2
        character(len=24), intent(in) :: metaIn
        character(len=24), intent(out) :: metaOut
        aster_logical, intent(in) :: forTemper
        character(len=24), optional, intent(in) :: metaPrev_
        real(kind=8), optional, intent(in) :: time_0_
        character(len=24), optional, intent(in) :: temp_0_
! ----- Local
        integer(kind=8), parameter :: nbOutMax = 2, nbInMax = 9
        character(len=8) :: lpaout(nbOutMax), lpain(nbInMax)
        character(len=24) :: lchin(nbInMax), lchout(nbOutMax)
        integer(kind=8) :: nbOut, nbIn
        character(len=1), parameter :: base = "V"
        character(len=16), parameter :: option = "META_ELNO"
        character(len=24), parameter :: chtime = '&&SMEVOL.CH_INST_R'
        integer(kind=8), parameter :: nbTimePara = 3
        character(len=8), parameter :: timeParaName(nbTimePara) = &
                                       (/'INST    ', 'DELTA01 ', 'DELTA12 '/)
        real(kind=8) :: timeParaVale(nbTimePara)
        real(kind=8) :: delta01, delta12
!   ------------------------------------------------------------------------------------------------
!
        metaOut = "&&SMEVOL.PHAS_META3"
        lpain = " "
        lpaout = " "
        lchin = " "
        lchout = " "

! ----- Create map for time
        delta01 = r8nnem()
        delta12 = time_2-time_1
        if (present(time_0_)) then
            delta01 = time_1-time_0_
        end if
        timeParaVale(1) = time_1
        timeParaVale(2) = delta01
        timeParaVale(3) = delta12
        call mecact('V', chtime, 'MODELE', metaParaOperator%model, 'INST_R  ', &
                    ncmp=nbTimePara, lnomcmp=timeParaName, vr=timeParaVale)

! ----- Input fields
        lpain(1) = 'PMATERC'
        lchin(1) = metaParaOperator%materialCoding
        lpain(2) = 'PTEMPER'
        lchin(2) = temp_1
        lpain(3) = 'PTEMPIR'
        lchin(3) = temp_2
        lpain(4) = 'PTIMMTR'
        lchin(4) = chtime
        lpain(5) = 'PPHASIN'
        lchin(5) = metaIn
        nbIn = 5

        if (forTemper) then
            nbIn = nbIn+1
            lpain(nbIn) = 'PCOMPMT'
            lchin(nbIn) = metaParaOperator%comporMetaTemper
            nbIn = nbIn+1
            lpain(nbIn) = 'PPHASEP'
            lchin(nbIn) = metaPrev_
        else
            nbIn = nbIn+1
            lpain(nbIn) = 'PCOMPME'
            lchin(nbIn) = metaParaOperator%comporMeta
            nbIn = nbIn+1
            lpain(nbIn) = 'PTEMPAR'
            lchin(nbIn) = temp_0_
            nbIn = nbIn+1
            lpain(nbIn) = 'PFTRC'
            lchin(nbIn) = metaParaOperator%TRCField
        end if
        ASSERT(nbIn .le. nbInMax)

! ----- Output field
        lpaout(1) = 'PPHASOUT'
        lchout(1) = metaOut
        nbOut = 1
        ASSERT(nbOut .le. nbOutMax)

! ----- Field with dynamical size (VARI_R) => allocate from compor map
        call copisd('CHAM_ELEM_S', 'V', metaParaOperator%comporMeta, metaOut)
        if (forTemper) then
            call copisd('CHAM_ELEM_S', 'V', metaParaOperator%comporMetaTemper, metaOut)
        end if

! ----- Compute
        call calcul('S', option, metaParaOperator%modelLigrel, &
                    nbIn, lchin, lpain, &
                    nbOut, lchout, lpaout, &
                    base, 'OUI')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaComp
!
! Compute - Main
!
! In  metaParaOperator : datastructure for parameters of CALC_META operator
! In  numphi
!
! --------------------------------------------------------------------------------------------------
    subroutine metaComp(metaParaOperator, numphi)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(in) :: metaParaOperator
        integer(kind=8), intent(in) :: numphi
! ----- Local
        aster_logical, parameter :: forTemper = ASTER_FALSE
        character(len=24) :: resultField
        integer(kind=8) :: iStore, jvPara, iret
        character(len=24) :: metaIn, metaOut
        integer(kind=8) :: numeStore_0, numeStore_1, numeStore_2
        character(len=24) :: temp_0, temp_1, temp_2
        real(kind=8) :: time_0, time_1, time_2
!   ------------------------------------------------------------------------------------------------
!

! ----- Main loop to compute
        do iStore = 1, metaParaOperator%nbStore-2

            numeStore_0 = metaParaOperator%listStore(iStore)
            if (numeStore_0 .lt. metaParaOperator%listStore(numphi)) cycle

! --------- Compute metallurgy from time 1 to time 2
            numeStore_1 = metaParaOperator%listStore(iStore+1)
            numeStore_2 = metaParaOperator%listStore(iStore+2)

! --------- Get fields
            call rsexch('F', metaParaOperator%resultName, 'TEMP', numeStore_0, temp_0, iret)
            call rsexch('F', metaParaOperator%resultName, 'TEMP', numeStore_1, temp_1, iret)
            call rsexch('F', metaParaOperator%resultName, 'META_ELNO', numeStore_1, metaIn, iret)
            call rsexch('F', metaParaOperator%resultName, 'TEMP', numeStore_2, temp_2, iret)

! --------- Get times
            call rsadpa(metaParaOperator%resultName, 'L', 1, 'INST', numeStore_0, 0, sjv=jvPara)
            time_0 = zr(jvPara)
            call rsadpa(metaParaOperator%resultName, 'L', 1, 'INST', numeStore_1, 0, sjv=jvPara)
            time_1 = zr(jvPara)
            call rsadpa(metaParaOperator%resultName, 'L', 1, 'INST', numeStore_2, 0, sjv=jvPara)
            time_2 = zr(jvPara)

! --------- Compute option META_ELNO
            call metaCompMetaElno(metaParaOperator, &
                                  time_1, time_2, &
                                  temp_1, temp_2, &
                                  metaIn, metaOut, &
                                  forTemper, &
                                  time_0_=time_0, temp_0_=temp_0)

! --------- Save metallurgy field
            call rsexch(' ', metaParaOperator%resultName, 'META_ELNO', numeStore_2, resultField, &
                        iret)
            call copisd('CHAMP_GD', 'G', metaOut, resultField)
            call rsnoch(metaParaOperator%resultName, 'META_ELNO', numeStore_2)
            call utmess('I', 'ARCHIVAGE_6', sk='META_ELNO', si=numeStore_2, sr=time_2)

! --------- Save behaviour map
            call rsexch(' ', metaParaOperator%resultName, 'COMPORMETA', numeStore_2, resultField, &
                        iret)
            call copisd('CHAMP_GD', 'G', metaParaOperator%comporMeta, resultField)
            call rsnoch(metaParaOperator%resultName, 'COMPORMETA', numeStore_2)
            call utmess('I', 'ARCHIVAGE_6', sk='COMPORMETA', si=numeStore_2, sr=time_2)
        end do

! ----- Cleaning
        call jedetr('&&SMEVOL_FTRC')
        call jedetr('&&SMEVOL_TRC')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaTemper
!
! Compute tempering - Main
!
! In  metaParaOperator : datastructure for parameters of CALC_META operator
!
! --------------------------------------------------------------------------------------------------
    subroutine metaTemper(metaParaOperator)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(in) :: metaParaOperator
! ----- Local
        aster_logical, parameter :: forTemper = ASTER_TRUE
        character(len=24), parameter :: metaOutInit = "&&SMEVOL.PHAS_META4"
        character(len=24) :: resultField
        integer(kind=8) :: iStore, jvPara, iret
        character(len=24) :: metaIn, metaOut, metaPrev
        integer(kind=8) :: numeStore_1, numeStore_2, numeStoreInit
        character(len=24) :: temp_1, temp_2
        real(kind=8) :: time_1, time_2, timeInit
!   ------------------------------------------------------------------------------------------------
!

! ----- Initializations on first field
        numeStoreInit = metaParaOperator%listStore(1)
        call copisd('CHAM_ELEM_S', 'V', metaParaOperator%comporMetaTemper, metaOutInit)
        call rsexch('F', metaParaOperator%resultName, 'META_ELNO', numeStoreInit, metaIn, iret)
        call rsadpa(metaParaOperator%resultName, 'L', 1, 'INST', numeStoreInit, 0, sjv=jvPara)
        timeInit = zr(jvPara)
        call metaCompTransfField(metaParaOperator, metaIn, metaOutInit)
        call copisd('CHAMP_GD', 'G', metaOutInit, metaIn)
        call utmess('I', 'ARCHIVAGE_6', sk='META_ELNO', si=numeStoreInit, sr=timeInit)

! ----- Main loop to compute
        do iStore = 1, metaParaOperator%nbStore-1

! --------- Compute tempering from time 1 to time 2
            numeStore_1 = metaParaOperator%listStore(iStore)
            numeStore_2 = metaParaOperator%listStore(iStore+1)

! --------- Get fields (temperature)
            call rsexch('F', metaParaOperator%resultName, 'TEMP', numeStore_1, temp_1, iret)
            call rsexch('F', metaParaOperator%resultName, 'TEMP', numeStore_2, temp_2, iret)

! --------- Get fields (previous metallurgy with tempering)
            call rsexch('F', metaParaOperator%resultName, 'META_ELNO', numeStore_1, metaPrev, iret)

! --------- Get fields (previous metallurgy without tempering)
            call rsexch('F', metaParaOperator%resultName, 'META_ELNO', numeStore_2, metaIn, iret)

! --------- Get times
            call rsadpa(metaParaOperator%resultName, 'L', 1, 'INST', numeStore_1, 0, sjv=jvPara)
            time_1 = zr(jvPara)
            call rsadpa(metaParaOperator%resultName, 'L', 1, 'INST', numeStore_2, 0, sjv=jvPara)
            time_2 = zr(jvPara)

! --------- Compute option META_ELNO
            call metaCompMetaElno(metaParaOperator, &
                                  time_1, time_2, &
                                  temp_1, temp_2, &
                                  metaIn, metaOut, &
                                  forTemper, &
                                  metaPrev_=metaPrev)

! --------- Save metallurgy field
            call rsexch(' ', metaParaOperator%resultName, 'META_ELNO', numeStore_2, resultField, &
                        iret)
            call copisd('CHAMP_GD', 'G', metaOut, resultField)
            call rsnoch(metaParaOperator%resultName, 'META_ELNO', numeStore_2)
            call utmess('I', 'ARCHIVAGE_6', sk='META_ELNO', si=numeStore_2, sr=time_2)

! --------- Save behaviour map
            call rsexch(' ', metaParaOperator%resultName, 'COMPORMETA', numeStore_2, resultField, &
                        iret)
            call copisd('CHAMP_GD', 'G', metaParaOperator%comporMetaTemper, resultField)
            call rsnoch(metaParaOperator%resultName, 'COMPORMETA', numeStore_2)
            call utmess('I', 'ARCHIVAGE_6', sk='COMPORMETA', si=numeStore_2, sr=time_2)

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaGetParameters
!
! Get parameters of operator
!
! Out metaParaOperator : datastructure for parameters of CALC_META operator
!
! --------------------------------------------------------------------------------------------------
    subroutine metaGetParameters(metaParaOperator)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(out) :: metaParaOperator
! ----- Local
        integer(kind=8) :: nbOption, nbStore, numeStore0
        character(len=8) :: resultName, model, materialField
        character(len=16) :: resultType
        character(len=24) :: materialCoding
        integer(kind=8) :: nbocc, nbRet
!   ------------------------------------------------------------------------------------------------
!

! ----- Get output result
        call getvid(' ', 'RESULTAT', scal=resultName, nbret=nbRet)
        call gettco(resultName, resultType)
        ASSERT(resultType .eq. 'EVOL_THER')
        metaParaOperator%resultName = resultName

! ----- Get all storing index in result
        call rs_get_liststore(resultName, nbStore)
        if (nbStore .ne. 0) then
            call wkvect(metaParaOperator%listStoreJv, "V V I", &
                        nbStore, vi=metaParaOperator%listStore)
            call rs_get_liststore(resultName, nbStore, metaParaOperator%listStore)
        end if
        if (nbStore .lt. 2) then
            call utmess('F', 'METALLURGY1_10')
        end if
        metaParaOperator%nbStore = nbStore

! ----- Get main parameters
        materialCoding = ' '
        model = ' '
        materialField = ' '
        numeStore0 = metaParaOperator%listStore(1)
        call rslesd(resultName, numeStore0, model, materialField)
        if (materialField .ne. ' ') then
            call rcmfmc(materialField, materialCoding, l_ther_=ASTER_TRUE)
        end if
        metaParaOperator%materialCoding = materialCoding
        metaParaOperator%model = model
        metaParaOperator%materialField = materialField

! ----- Get options to compute
        call getvtx(' ', 'OPTION', nbval=0, nbret=nbRet)
        nbOption = -nbRet
        call wkvect(metaParaOperator%listOptionsJv, &
                    "V V K16", nbOption, vk16=metaParaOperator%listOption)
        call getvtx(' ', 'OPTION', nbval=nbOption, vect=metaParaOperator%listOption, nbret=nbRet)
        metaParaOperator%nbOption = nbOption

! ----- Metallirgical ?
        call getfac('COMPORTEMENT', nbocc)
        metaParaOperator%hasMetaLaw = nbocc .ne. 0

! ----- Tempering ?
        call getfac('REVENU', nbocc)
        metaParaOperator%hasTemperLaw = nbocc .ne. 0
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaDelParameters
!
! Delete parameters of operator
!
! In  metaParaOperator : datastructure for parameters of CALC_META operator
!
! --------------------------------------------------------------------------------------------------
    subroutine metaDelParameters(metaParaOperator)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(in) :: metaParaOperator
!   ------------------------------------------------------------------------------------------------
!
        call jedetr(metaParaOperator%listOptionsJv)
        call jedetr(metaParaOperator%listStoreJv)
        call detrsd('CARTE', metaParaOperator%comporMeta)
        call detrsd('CARTE', metaParaOperator%comporMetaTemper)
        if (metaParaOperator%hasTRC) then
            call jedetr('&&SMEVOL_FTRC')
            call jedetr('&&SMEVOL_TRC')
            call detrsd('CARTE', metaParaOperator%TRCField)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaCompPhases
!
! Computes phases
!
! IO  metaParaOperator : datastructure for parameters of CALC_META operator
!
! --------------------------------------------------------------------------------------------------
    subroutine metaCompPhases(metaParaOperator)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(inout) :: metaParaOperator
! ----- Local
        character(len=24) :: metaInitUser
        integer(kind=8) :: numeFieldInit
        character(len=16), parameter :: factorKeyword = "COMPORTEMENT"
        character(len=16), parameter :: factorKeywordTemper = "REVENU"
        integer(kind=8) :: numphi, numeStore_1, numeStore_2, iStore
        aster_logical :: lEmptyMap
!   ------------------------------------------------------------------------------------------------
!
        call dismoi('NOM_LIGREL', metaParaOperator%model, 'MODELE', &
                    repk=metaParaOperator%modelLigrel)

        if (metaParaOperator%hasMetaLaw) then
! --------- Get initial state
            call metaGetInitialState(metaParaOperator%resultName, metaInitUser, numeFieldInit)

! --------- Construct map for behaviour
            call mtdorc(factorKeyword, metaParaOperator%model, metaParaOperator%comporMeta, &
                        lEmptyMap)
            ASSERT(.not. lEmptyMap)

! --------- Prepare dynamic field from behaviour
            call detrsd('CHAM_ELEM_S', metaParaOperator%comporMeta)
            call cesvar(' ', metaParaOperator%comporMeta, &
                        metaParaOperator%modelLigrel, metaParaOperator%comporMeta)

! --------- Prepare field to manage TRC curves in elementary computation
            call metaPrepTRCWorkingField(metaParaOperator)

! --------- Compute initial field for metallurgy if required and store it
            if (numeFieldInit .eq. 0) then
                numphi = 1
                numeStore_1 = metaParaOperator%listStore(1)
                call metaPrepInitialState(metaParaOperator, numeStore_1, metaInitUser)
                numeStore_2 = metaParaOperator%listStore(2)
                call metaPrepInitialState(metaParaOperator, numeStore_2, metaInitUser)
            else
                numphi = 0
                do iStore = 2, metaParaOperator%nbStore
                    if (metaParaOperator%listStore(iStore) .eq. numeFieldInit) then
                        numphi = iStore-1
                    end if
                end do
            end if

! --------- Apply metallurgical behaviour
            call metaComp(metaParaOperator, numphi)
        end if

! ----- Tempering
        if (metaParaOperator%hasTemperLaw) then
! --------- Construct map behaviour (tempering)
            call mtdorc(factorKeywordTemper, metaParaOperator%model, &
                        metaParaOperator%comporMetaTemper, &
                        lEmptyMap)
            ASSERT(.not. lEmptyMap)

! --------- Prepare dynamic field from behaviour
            call detrsd('CHAM_ELEM_S', metaParaOperator%comporMetaTemper)
            call cesvar(' ', metaParaOperator%comporMetaTemper, metaParaOperator%modelLigrel, &
                        metaParaOperator%comporMetaTemper)

! --------- Apply tempering
            call metaTemper(metaParaOperator)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaCompOtherOptions
!
! Computes other options
!
! In  metaParaOperator : datastructure for parameters of CALC_META operator
! In  option           : option to compute
!
! --------------------------------------------------------------------------------------------------
    subroutine metaCompOtherOptions(metaParaOperator, option)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(in) :: metaParaOperator
        character(len=16), intent(in) :: option
! ----- Local
        character(len=16), parameter :: factorKeyword = "COMPORTEMENT"
        character(len=16), parameter :: resultType = "EVOL_THER"
        integer(kind=8) :: nbComp, iComp, nbRet, iret
        character(len=16) :: metaType
!   ------------------------------------------------------------------------------------------------
!
        nbComp = 0
        call getfac(factorKeyword, nbComp)
        do iComp = 1, nbComp
            call getvtx(factorKeyword, 'RELATION', iocc=iComp, scal=metaType, nbret=nbRet)
            if (nbRet .ne. 0) then
                if (metaType(1:5) .ne. 'ACIER' .and. option .eq. 'DURT_ELNO') then
                    call utmess('F', 'META1_3', sk=metaType)
                end if
            end if
        end do
        call calcop(option, metaParaOperator%listOptionsJv, metaParaOperator%resultName, &
                    metaParaOperator%resultName, metaParaOperator%listStoreJv, &
                    metaParaOperator%nbStore, resultType, iret)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaCompTransfField
!
! Transfer field from standard to tempering phases
!
! In  metaParaOperator : datastructure for parameters of CALC_META operator
! In  metaIn           : phases without tempering
! In  metaOut          : phases with tempering
!
! --------------------------------------------------------------------------------------------------
    subroutine metaCompTransfField(metaParaOperator, metaIn, metaOut)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(META_ParaOperator), intent(in) :: metaParaOperator
        character(len=24), intent(in) :: metaIn, metaOut
! ----- Local
        integer(kind=8), parameter :: nbout = 1, nbInMax = 1
        character(len=8) :: lpaout(nbout), lpain(nbInMax)
        character(len=24) :: lchin(nbInMax)
        character(len=1), parameter :: base = "V"
        character(len=16), parameter :: option = "META_TRAN_ELNO"
!   ------------------------------------------------------------------------------------------------
!
        lpain(1) = 'PPHASIN'
        lchin(1) = metaIn
        lpaout(1) = 'PPHASOUT'
        call calcul('S', option, metaParaOperator%modelLigrel, nbInMax, lchin, &
                    lpain, nbout, metaOut, lpaout, base, &
                    'OUI')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module MetallurgyOperator_module
