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
subroutine getExternalStateVariable(relaComp, relaCompPy, &
                                    l_mfront_offi, l_mfront_proto, &
                                    adrsMGIS, variExteCode)
!
    use NonLin_Datastructure_type
    use Behaviour_module
    implicit none
!
#include "asterc/lcextevari.h"
#include "asterc/lcinfo.h"
#include "asterc/mgis_get_esvs.h"
#include "asterc/mgis_get_number_of_esvs.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/iscode.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: relaComp, relaCompPy
    aster_logical, intent(in) :: l_mfront_offi, l_mfront_proto
    character(len=16), intent(in) :: adrsMGIS
    integer(kind=8), intent(out) :: variExteCode(2)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get external states variables
!
! --------------------------------------------------------------------------------------------------
!
! In  relaComp         : behaviour (RELATION keyword)
! In  relaCompPY       : behaviour (RELATION keyword) - For Python
! In  l_mfront_proto   : .true. if MFront prototype
! In  l_mfront_offi    : .true. if MFront official
! In  adrsMGIS         : address (hexadecimal) for the MGIS Behaviour
! Out variExteCode     : coded integers for external state variable
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbVarc, iVarc, idummy1, idummy2, iExteType
    integer(kind=8), parameter :: nbExteType = VARC_EXTE_NBTYPE
    character(len=64) :: varcNameExte(VARC_EXTE_NBMAXI)
    character(len=8) :: varcName
    integer(kind=8) :: tabcod(60)
    character(len=8), parameter :: varcNameList(nbExteType) = (/ &
                                   'ELTSIZE1', 'COORGA  ', &
                                   'GRADVELO', 'HYGR    ', 'NEUT1   ', &
                                   'NEUT2   ', 'TEMP    ', 'DTX     ', &
                                   'DTY     ', 'DTZ     ', 'X       ', &
                                   'Y       ', 'Z       ', 'SECH    ', &
                                   'HYDR    ', 'CORR    ', 'IRRA    ', &
                                   'EPSAXX  ', 'EPSAYY  ', 'EPSAZZ  ', &
                                   'EPSAXY  ', 'EPSAXZ  ', 'EPSAYZ  ', &
                                   'FERRITE ', 'PERLITE ', 'BAINITE ', &
                                   'MARTENS ', 'ALPHPUR ', 'ALPHBET ', &
                                   'TIME    ', 'TEMPREFE'/)
    aster_logical, parameter :: l_allow_mfront(nbExteType) = (/.true., .false., &
                                                               .false., .true., .true., &
                                                               .true., .true., .true., &
                                                               .true., .true., .true., &
                                                               .true., .true., .true., &
                                                               .true., .true., .true., &
                                                               .true., .true., .true., &
                                                               .true., .true., .true., &
                                                               .true., .true., .true., &
                                                               .true., .true., .true., &
                                                               .true., .true./)
!
! --------------------------------------------------------------------------------------------------
!
    variExteCode = 0

! - Get names of external state variables
    nbVarc = 0
    varcNameExte = ' '
    if (l_mfront_proto .or. l_mfront_offi) then
        call mgis_get_number_of_esvs(adrsMGIS, nbVarc)
        ASSERT(nbVarc .le. VARC_EXTE_NBMAXI)
        call mgis_get_esvs(adrsMGIS, varcNameExte)
    else
        call lcinfo(relaCompPy, idummy1, idummy2, nbVarc)
        ASSERT(nbVarc .le. VARC_EXTE_NBMAXI)
        call lcextevari(relaCompPy, nbVarc, varcNameExte)
    end if

! - Print
    if (nbVarc .gt. 0) then
        call utmess('I', 'COMPOR4_21', si=nbVarc, sk=relaComp)
        do iVarc = 1, nbVarc
            call utmess('I', 'COMPOR4_22', si=iVarc, sk=varcNameExte(iVarc))
        end do
    end if

! - Coding
    tabcod = 0
    do iVarc = 1, nbVarc
        do iExteType = 1, nbExteType
            varcName = getAsterVariableName(varcNameExte(iVarc))
            if (varcName .eq. varcNameList(iExteType)) then
                tabcod(iExteType) = 1
                if (.not. l_allow_mfront(iExteType) .and. &
                    (l_mfront_proto .or. l_mfront_offi)) then
                    call utmess('I', 'COMPOR2_25', sk=varcNameExte(iVarc))
                    tabcod(iExteType) = 0
                end if
            end if
        end do
    end do
    call iscode(tabcod, variExteCode, 60)
!
end subroutine
