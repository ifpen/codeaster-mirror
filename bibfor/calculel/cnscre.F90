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

subroutine cnscre(maz, nomgdz, ncmp, licmp, basez, &
                  cnsz, undf0_)
! A_UTIL
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeecra.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeundf.h"
#include "asterfort/verigd.h"
#include "asterfort/wkvect.h"
#include "asterfort/isParallelMesh.h"
    character(len=*) :: maz, nomgdz, cnsz, basez
    integer(kind=8) :: ncmp
    character(len=*) :: licmp(ncmp)
    aster_logical, optional, intent(in) :: undf0_
! ------------------------------------------------------------------
! BUT : CREER UN CHAM_NO_S
! ------------------------------------------------------------------
!     ARGUMENTS:
! MAZ     IN/JXIN  K8  : SD MAILLAGE DE CNSZ
! NOMGDZ  IN       K8  : NOM DE LA GRANDEUR DE CNSZ
! NCMP    IN       I   : NOMBRE DE CMPS VOULUES DANS CNSZ
! LICMP   IN       L_K8: NOMS DES CMPS VOULUES DANS CNSZ
! BASEZ   IN       K1  : BASE DE CREATION POUR CNSZ : G/V/L
! CNSZ    IN/JXOUT K19 : SD CHAM_NO_S A CREER
! undf0   IN           : flag to init field with zeros
!     ------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=1) :: base
    character(len=3) :: tsca
    character(len=8) :: ma, nomgd
    character(len=19) :: cns
    integer(kind=8) :: nbno, jcnsk, jcnsd
    integer(kind=8) :: jcnsc, k, jcnsl, jcnsv, iret
    aster_logical :: undf0, l_pmesh, l_empty
!     ------------------------------------------------------------------
!
    call jemarq()
    cns = cnsz
    base = basez
    nomgd = nomgdz
    ma = maz

    undf0 = ASTER_FALSE
    if (present(undf0_)) then
        undf0 = undf0_
    end if
!
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbno)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    l_pmesh = isParallelMesh(ma)
!
!
!     -- SI CNS EXISTE DEJA, ON LE DETRUIT :
    call detrsd('CHAM_NO_S', cns)
!
!------------------------------------------------------------------
!     1- QUELQUES VERIFS :
!     ------------------------
    call verigd(nomgd, licmp, ncmp, iret)
    ASSERT(iret .le. 0)
!
!------------------------------------------------------------------
!     2- CREATION DE CNS.CNSK:
!     ------------------------
    call wkvect(cns//'.CNSK', base//' V K8', 3, jcnsk)
    zk8(jcnsk-1+1) = ma
    zk8(jcnsk-1+2) = nomgd
    zk8(jcnsk-1+3) = 'NOEU'
!
!------------------------------------------------------------------
!     3- CREATION DE CNS.CNSD:
!     ------------------------
    call wkvect(cns//'.CNSD', base//' V I', 2, jcnsd)
    zi(jcnsd-1+1) = nbno
    zi(jcnsd-1+2) = ncmp
!
!------------------------------------------------------------------
!     4- CREATION DE CNS.CNSC:
!     ------------------------
    l_empty = .false.
    if (ncmp .eq. 0 .and. l_pmesh) then
        l_empty = .true.
        ncmp = 1
    end if
    call wkvect(cns//'.CNSC', base//' V K8', ncmp, jcnsc)
    if (.not. l_empty) then
        do k = 1, ncmp
            zk8(jcnsc-1+k) = licmp(k)
        end do
    end if
!
!------------------------------------------------------------------
!     5- CREATION DE CNS.CNSL:
!     ------------------------
    call wkvect(cns//'.CNSL', base//' V L', max(nbno*ncmp, 1), jcnsl)
    call jeecra(cns//'.CNSL', 'LONUTI', nbno*ncmp)
    call jeundf(cns//'.CNSL')
!
!------------------------------------------------------------------
!     6- CREATION DE CNS.CNSV:
!     ------------------------
    call wkvect(cns//'.CNSV', base//' V '//tsca, max(nbno*ncmp, 1), jcnsv)
    call jeecra(cns//'.CNSV', 'LONUTI', nbno*ncmp)
    call jeundf(cns//'.CNSV', undf0)
!
    call jedema()
end subroutine
