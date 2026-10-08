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
! aslint: disable=W1502
!
subroutine lrvema(nomail, mfich, nochmd)
!
    use as_med_module, only: as_med_open
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_mfdfin.h"
#include "asterfort/as_mfdncn.h"
#include "asterfort/as_mficlo.h"
#include "asterfort/as_mmhnme.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/get_field_names_from_medfile.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/ulisog.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#ifdef ASTER_HAVE_MED
#include "med_parameter.hf"
#endif
#include "MeshTypes_type.h"
!
    integer(kind=8) :: mfich
    character(len=8) :: nomail
    character(len=64) :: nochmd
!
! --------------------------------------------------------------------------------------------------
!
!   BUT : ROUTINE DE LIRE RESU / LIRE_CHAMP QUI VERIFIE LA COHERENCE
!       ENTRE LE MAILLAGE FOURNI ET LES DONNEES DU FICHIER MED
!
! --------------------------------------------------------------------------------------------------
!
! IN  :
!   NOMAIL  K8   NOM DU MAILLAGE ASTER
!   MFICH   I    NUMERO DU FICHIER MED
!   NOCHMD  K64  NOM D'UN CHAMP REPOSANT SUR LE MAILLAGE MED A VERIFIER
!
! --------------------------------------------------------------------------------------------------
!
#ifdef ASTER_HAVE_MED
    integer(kind=8) :: iret, nmatyp, ncmp
    integer(kind=8) :: nbma, jnbtyp, jmatyp, nbtym, nbtv, codret
    med_idt :: idfimd
    integer(kind=8) :: i, j, iaux, jnbty2, nchmed
    integer(kind=8) :: vali(2), lnomam
    integer(kind=8), parameter :: edlect = 0, edconn = 1, edmail = 0, ednoda = 0
    character(len=1) :: k1b
    character(len=8) :: saux08
    character(len=64) :: nomamd
    character(len=200) :: nofimd
    character(len=255) :: kfic
    aster_logical :: lfirst
    character(len=16), pointer :: cname(:) => null()
    character(len=16), pointer :: cunit(:) => null()
    character(len=80), pointer :: v_names(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    character(len=8), parameter :: nomast(MT_NTYMAX) = (/'POI1    ', 'SEG2    ', 'SEG22   ', &
                                                         'SEG3    ', 'SEG33   ', 'SEG4    ', &
                                                         'TRIA3   ', 'TRIA33  ', 'TRIA6   ', &
                                                         'TRIA66  ', 'TRIA7   ', 'QUAD4   ', &
                                                         'QUAD44  ', 'QUAD8   ', 'QUAD88  ', &
                                                         'QUAD9   ', 'QUAD99  ', 'TETRA4  ', &
                                                         'TETRA10 ', 'PENTA6  ', 'PENTA15 ', &
                                                         'PENTA18 ', 'PYRAM5  ', 'PYRAM13 ', &
                                                         'HEXA8   ', 'HEXA20  ', 'HEXA27  ', &
                                                         'TETRA15 ', 'PENTA21 ', 'PYRAM19 ', &
                                                         'TR3QU4  ', 'QU4TR3  ', 'TR6TR3  ', &
                                                         'TR3TR6  ', 'TR6QU4  ', 'QU4TR6  ', &
                                                         'TR6QU8  ', 'QU8TR6  ', 'TR6QU9  ', &
                                                         'QU9TR6  ', 'QU8TR3  ', 'TR3QU8  ', &
                                                         'QU8QU4  ', 'QU4QU8  ', 'QU8QU9  ', &
                                                         'QU9QU8  ', 'QU9QU4  ', 'QU4QU9  ', &
                                                         'QU9TR3  ', 'TR3QU9  ', 'SEG32   ', &
                                                         'SEG23   ', 'QU4QU4  ', 'TR3TR3  ', &
                                                         'HE8HE8  ', 'PE6PE6  ', 'TE4TE4  ', &
                                                         'QU8QU8  ', 'TR6TR6  ', 'SE2TR3  ', &
                                                         'SE2TR6  ', 'SE2QU4  ', 'SE2QU8  ', &
                                                         'SE2QU9  ', 'SE3TR3  ', 'SE3TR6  ', &
                                                         'SE3QU4  ', 'SE3QU8  ', 'SE3QU9  ', &
                                                         'HEXA9   ', 'PENTA7  ', 'TR3SE2  ', &
                                                         'TR3SE3  ', 'TR6SE2  ', 'TR6SE3  ', &
                                                         'QU4SE2  ', 'QU4SE3  ', 'QU8SE2  ', &
                                                         'QU8SE3  ', 'QU9SE2  ', 'QU9SE3  ', &
                                                         'QU4TR7  ', 'QU8TR7  ', 'QU9TR7  ', &
                                                         'TR3TR7  ', 'TR6TR7  ', 'TRIA10  ', &
                                                         'QUAD12  ', 'TETRA20 ', 'SEG44   ', &
                                                         'SEG43   ', 'TE4TR7  ', 'T10TR7  ', &
                                                         'T20TR7  '/)
    integer(kind=8), parameter :: nummed(MT_NTYMAX) = (/ &
                                  MED_POINT1, MED_SEG2, MED_UNDEF_GEOTYPE, &
                                  MED_SEG3, MED_UNDEF_GEOTYPE, MED_SEG4, &
                                  MED_TRIA3, MED_UNDEF_GEOTYPE, MED_TRIA6, &
                                  MED_UNDEF_GEOTYPE, MED_TRIA7, MED_QUAD4, &
                                  MED_UNDEF_GEOTYPE, MED_QUAD8, MED_UNDEF_GEOTYPE, &
                                  MED_QUAD9, MED_UNDEF_GEOTYPE, MED_TETRA4, &
                                  MED_TETRA10, MED_PENTA6, MED_PENTA15, &
                                  MED_PENTA18, MED_PYRA5, MED_PYRA13, &
                                  MED_HEXA8, MED_HEXA20, MED_HEXA27, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, MED_UNDEF_GEOTYPE, &
                                  MED_UNDEF_GEOTYPE/)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
!     ON VERIFIE QUE LE MAILLAGE FOURNI ET CELUI
!     CONTENU DANS LE FICHIER MED ONT
!     - MEME TYPE DE MAILLE
!     - MEME NOMBRE DE MAILLE PAR TYPE
!     =================================
!
    call ulisog(mfich, kfic, k1b)
    if (kfic(1:1) .eq. ' ') then
        call codent(mfich, 'G', saux08)
        nofimd = 'fort.'//saux08
    else
        nofimd = kfic(1:200)
    end if
!
    nomamd = ' '
    call as_med_open(idfimd, nofimd, edlect, iaux)
    if (iaux .ne. 0) then
        lnomam = lxlgut(saux08)
        call utmess('F', 'MED_78', sk=saux08(1:lnomam))
    end if
!
    call as_mfdncn(idfimd, nochmd, ncmp, codret)
    if (codret .ne. 0) then
        call get_field_names_from_medfile(idfimd, '&&LRVEMA.CHAMNOM')
        call utmess('F+', 'MED_57', sk=nochmd)
        call jelira('&&LRVEMA.CHAMNOM', 'LONMAX', nchmed)
        call jeveuo('&&LRVEMA.CHAMNOM', 'L', vk80=v_names)
        do iaux = 1, nchmed
            call utmess('F+', 'MED2_2', sk=v_names(iaux))
        end do
        call utmess('F', 'VIDE_1')
    end if
    AS_ALLOCATE(vk16=cname, size=ncmp)
    AS_ALLOCATE(vk16=cunit, size=ncmp)
!
    call as_mfdfin(idfimd, nochmd, nomamd, nbtv, cunit(1), &
                   cname(1), codret)
!
    call wkvect('&&LRVERIMO_NBETYP1', 'V V I', int(MT_NTYMAX, 8), jnbtyp)
    call wkvect('&&LRVERIMO_NBETYP2', 'V V I', int(MT_NTYMAX, 8), jnbty2)
    do i = 1, MT_NTYMAX
        zi(jnbtyp+i-1) = 0
        if (nummed(i) .ne. 0) then
            call as_mmhnme(idfimd, nomamd, edconn, edmail, nummed(i), &
                           ednoda, nmatyp, iret)
            zi(jnbtyp+i-1) = nmatyp
        end if
    end do
!
    call as_mficlo(idfimd, iret)
    call dismoi('NB_MA_MAILLA', nomail, 'MAILLAGE', repi=nbma)
    call wkvect('&&LRVERIMO_NBMA_TYP', 'V V I', nbma, jmatyp)
    do i = 1, nbma
        zi(jmatyp+i-1) = 0
    end do
!
    call jeveuo(nomail//'.TYPMAIL', 'L', vi=typmail)
    do i = 1, nbma
        zi(jmatyp+i-1) = nummed(typmail(i))
    end do
!
    do i = 1, MT_NTYMAX
        nbtym = 0
        zi(jnbty2+i-1) = nbtym
        if (nummed(i) .ne. 0) then
            do j = 1, nbma
                if (zi(jmatyp+j-1) .eq. nummed(i)) then
                    nbtym = nbtym+1
                end if
            end do
        end if
        zi(jnbty2+i-1) = nbtym
    end do
!
    lfirst = .true.
    do i = 1, MT_NTYMAX
        if (nummed(i) .ne. 0) then
            if (zi(jnbtyp+i-1) .ne. zi(jnbty2+i-1) .and. lfirst) then
                lfirst = .false.
                call utmess('F+', 'MED_54')
                if (zi(jnbtyp+i-1) .lt. zi(jnbty2+i-1)) then
                    vali(1) = zi(jnbtyp+i-1)
                    vali(2) = zi(jnbty2+i-1)
                    call utmess('F', 'MED_59', sk=nomast(i), ni=int(2, 8), vali=vali)
                else
                    vali(1) = zi(jnbtyp+i-1)
                    vali(2) = zi(jnbty2+i-1)
                    call utmess('F', 'MED_61', sk=nomast(i), ni=int(2, 8), vali=vali)
                end if
!
            end if
        end if
    end do
!
    call jedetr('&&LRVERIMO_NBETYP1')
    call jedetr('&&LRVERIMO_NBETYP2')
    call jedetr('&&LRVERIMO_NBMA_TYP')
    AS_DEALLOCATE(vk16=cname)
    AS_DEALLOCATE(vk16=cunit)
!
    call jedema()
!
#endif
end subroutine
