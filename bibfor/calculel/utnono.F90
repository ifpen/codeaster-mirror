! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

subroutine utnono(mess, mesh, type, nomgrp, nomobj, &
                  iret, lColleNd_)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: mesh, nomobj
    character(len=24) :: nomgrp
    character(len=*) :: mess, type
    integer(kind=8) :: iret
    aster_logical, optional, intent(in) :: lColleNd_
!
!     RENVOIE LE PREMIER NOEUD OU MAILLE CONTENU DANS UN GROUP_NO
!             OU UN GROUP_MA
!
! IN  : MESS   : TYPE DE MESSAGE UE L'ON VEUT IMPRIMER
!                'F'  MESSAGE FATAL
!                'E'  MESSAGE ERREUR
!                ' '  LE CODE RETOUR EST GERE PAR LE DEVELOPPEUR
! IN  : NOMMA  : NOM DU MAILLAGE.
! IN  : TYPE   : TRAITEMENT MAILLE OU NOEUD
! IN  : NOMGRP : NOM D'UN GROUP_NO OU D'UN GROUP_MA
! OUT : NOMOBJ : NOM DU NOEUD OU DE LA MAILLE
! OUT : IRET   : CODE RETOUR
!                 0 --> OK
!                10 --> LE GROUPE N'EXISTE PAS OU EST VIDE
!                 1 -->  PLUSIEURS NOEUDS OU MAILLES DANS LE GROUPE
! ----------------------------------------------------------------------
    character(len=24) :: valk(2)
!
    integer(kind=8) :: iret1, nbno, iad, ier
    character(len=1) :: typm
    character(len=8) :: knbno
    character(len=16) :: typgrp
    aster_logical :: lColleNd
!
    call jemarq()
    iret = 10
    typm = mess(1:1)
    nomobj = ' '
    lcolleNd = .false.
    if (present(lColleNd_)) then
        lColleNd = lColleNd_
    end if
!
    if (type(1:5) .eq. 'NOEUD') then
        typgrp = '.GROUPENO       '
    else if (type(1:6) .eq. 'MAILLE') then
        typgrp = '.GROUPEMA       '
    else
        goto 999
    end if
!
    call jeexin(jexnom(mesh//typgrp, nomgrp), iret1)
    if (iret1 .gt. 0) then
        call jelira(jexnom(mesh//typgrp, nomgrp), 'LONUTI', nbno)
    else
        nbno = 0
    end if
    if (nbno .eq. 0) then
        if (typm .eq. ' ') goto 999
        if (type(1:5) .eq. 'NOEUD') then
            call utmess(typm, 'ELEMENTS_67', sk=nomgrp)
        else
            call utmess(typm, 'ELEMENTS_73', sk=nomgrp)
        end if
        goto 999
    end if
!
    iret = 0
    ASSERT(nbno .gt. 0)
    if (nbno .ne. 1) then
        iret = 1
        if (typm .ne. ' ') then
            call codent(nbno, 'D', knbno)
            if (type(1:5) .eq. 'NOEUD') then
                valk(1) = nomgrp
                valk(2) = knbno
                call utmess(typm, 'CALCULEL5_20', nk=2, valk=valk)
            else
                valk(1) = nomgrp
                valk(2) = knbno
                call utmess(typm, 'CALCULEL5_21', nk=2, valk=valk)
            end if
            goto 999
        end if
    end if
!
    call jeveuo(jexnom(mesh//typgrp, nomgrp), 'L', iad)
    nomobj = int_to_char8(zi(iad), lcolleNd, mesh, "NOEUD")
!
999 continue
    call jedema()
end subroutine
