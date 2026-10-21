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
subroutine nmvcle(modelZ, materFieldZ, caraElemZ, time, varcZ)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/mecact.h"
#include "asterfort/vrcins.h"
!
    character(len=*), intent(in) :: modelZ, materFieldZ, caraElemZ
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: varcZ
!
! --------------------------------------------------------------------------------------------------
!
! External state variables
!
! Create fields
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  materField       : name of material characteristics (field)
! In  caraElem         : name of elementary characteristics (field)
! In  time             : time
! In  varc             : name of external state variable field
!
! --------------------------------------------------------------------------------------------------
!
    character(len=2) :: codret
    character(len=8) :: model, materField, caraElem
    character(len=14) :: varc
    character(len=24) :: varcAll, varcTime, modelLigrel
!
! --------------------------------------------------------------------------------------------------
!
    varc = varcZ
    model = modelZ
    caraElem = caraElemZ
    materField = materFieldZ
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)

! - Old object deleted
    call detrsd('VARI_COM', varc)

! - Name of field for all external state variables
    varcAll = varc//'.TOUT'

! - Construct field for all external state variables
    call vrcins(model, materField, caraElem, time, varcAll, codret)

! - Construct field for time
    varcTime = varc//'.INST'
    call mecact('V', varcTime, 'MODELE', modelLigrel, 'INST_R', &
                ncmp=1, nomcmp='INST', sr=time)
!
end subroutine
