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
subroutine dismme(questi, nomobz, repi, repkz, ierd)
    implicit none
!     --     DISMOI(MATR_ELEM OU VECT_ELEM)
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismmo.h"
#include "asterfort/dismre.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi
    character(len=*) :: nomobz, repkz
    character(len=32) :: repk
    character(len=19) :: nomob
! ----------------------------------------------------------------------
!    IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE CONCEPT MATR_ELEM
!    OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
!
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=7) :: typmat, kmpic, zero
    integer(kind=8) :: iret, i, i1, ialire, nbresu, iexi, size
    character(len=8) :: mo, mesh, pheno
    character(len=19) :: partit
    character(len=24), pointer :: rerr(:) => null()
!
!
!
    call jemarq()
    repk = ' '
    repi = 0
    ierd = 0
!
    nomob = nomobz
    call jeveuo(nomob//'.RERR', 'L', vk24=rerr)
    call jelira(nomob//'.RERR', 'LONMAX', size)
    mo = rerr(1) (1:8)
!
    if (questi .eq. 'NOM_MODELE') then
        repk = mo
!
    else if (questi .eq. 'TYPE_MATRICE') then
        repk = 'SYMETRI'
        call jeexin(nomob//'.RELR', iret)
        if (iret .gt. 0) then
            call jelira(nomob//'.RELR', 'LONUTI', nbresu)
            if (nbresu .gt. 0) call jeveuo(nomob//'.RELR', 'L', ialire)
            do i = 1, nbresu
                call jeexin(zk24(ialire-1+i) (1:19)//'.NOLI', iexi)
                if (iexi .eq. 0) goto 1
                call dismre(questi, zk24(ialire-1+i), repi, typmat, i1)
                if ((i1 .eq. 0) .and. (typmat .eq. 'NON_SYM')) then
                    repk = 'NON_SYM'
                    goto 999
                end if
1               continue
            end do
        end if
!
    else if (questi .eq. 'ZERO') then
        repk = 'OUI'
        call jeexin(nomob//'.RELR', iret)
        if (iret .gt. 0) then
            call jelira(nomob//'.RELR', 'LONUTI', nbresu)
            if (nbresu .gt. 0) call jeveuo(nomob//'.RELR', 'L', ialire)
            do i = 1, nbresu
                call jeexin(zk24(ialire-1+i) (1:19)//'.NOLI', iexi)
                if (iexi .eq. 0) goto 4
                call dismre(questi, zk24(ialire-1+i), repi, zero, i1)
                if ((i1 .eq. 0) .and. (zero .eq. 'NON')) then
                    repk = 'NON'
                    goto 999
                end if
4               continue
            end do
        end if
!
    else if (questi .eq. 'PARTITION') then
        repk = ' '
        call jeexin(nomob//'.RELR', iret)
        if (iret .gt. 0) then
            call jelira(nomob//'.RELR', 'LONUTI', nbresu)
            if (nbresu .gt. 0) call jeveuo(nomob//'.RELR', 'L', ialire)
            do i = 1, nbresu
                call jeexin(zk24(ialire-1+i) (1:19)//'.NOLI', iexi)
                if (iexi .eq. 0) goto 2
                call dismre(questi, zk24(ialire-1+i), repi, partit, i1)
                if (partit .ne. ' ' .and. repk .eq. ' ') repk = partit
                if (partit .ne. ' ') then
                    ASSERT(repk .eq. partit)
                end if
2               continue
            end do
        end if
!
    else if (questi .eq. 'MPI_COMPLET') then
        repk = ' '
        call jeexin(nomob//'.RELR', iret)
        if (iret .gt. 0) then
            call jelira(nomob//'.RELR', 'LONUTI', nbresu)
            if (nbresu .gt. 0) call jeveuo(nomob//'.RELR', 'L', ialire)
            do i = 1, nbresu
                call jeexin(zk24(ialire-1+i) (1:19)//'.NOLI', iexi)
                if (iexi .eq. 0) goto 21
                call dismre(questi, zk24(ialire-1+i), repi, kmpic, i1)
                if (i .eq. 1) then
                    repk = kmpic
                else
                    ASSERT(repk .eq. kmpic)
                end if
21              continue
            end do
        end if
!
    else if (questi .eq. 'NOM_MAILLA') then
        repk = ' '
        call jeexin(nomob//'.RELR', iret)
        if (iret .gt. 0) then
            call jelira(nomob//'.RELR', 'LONUTI', nbresu)
            if (nbresu .gt. 0) call jeveuo(nomob//'.RELR', 'L', ialire)
            do i = 1, nbresu
                call jeexin(zk24(ialire-1+i) (1:19)//'.NOLI', iexi)
                if (iexi .eq. 0) goto 22
                call dismre(questi, zk24(ialire-1+i), repi, mesh, i1)
                if (mesh .ne. ' ' .and. repk .eq. ' ') then
                    repk = mesh
                end if
                if (mesh .ne. ' ') then
                    ASSERT(repk .eq. mesh)
                end if
22              continue
            end do
        end if
!
    else if (questi .eq. 'PHENOMENE') then
        repk = ' '
        call jeexin(nomob//'.RELR', iret)
        if (iret .gt. 0) then
            call jelira(nomob//'.RELR', 'LONUTI', nbresu)
            if (nbresu .gt. 0) call jeveuo(nomob//'.RELR', 'L', ialire)
            do i = 1, nbresu
                call jeexin(zk24(ialire-1+i) (1:19)//'.NOLI', iexi)
                if (iexi .eq. 0) goto 23
                call dismre(questi, zk24(ialire-1+i), repi, pheno, i1)
                if (pheno .ne. ' ' .and. repk .eq. ' ') then
                    repk = pheno
                end if
                if (pheno .ne. ' ') then
                    ASSERT(repk .eq. pheno)
                end if
23              continue
            end do
        end if
!
    else if (questi .eq. 'SUR_OPTION') then
        ASSERT(size == 3)
        repk = rerr(2) (1:16)
!
    else if (questi .eq. 'NB_SS_ACTI') then
        call dismmo(questi, mo, repi, repk, ierd)
        if (repi > 0) then
            if (size == 1) then
            !! Case VECT_ELEM
                call jeexin(nomob//'.RELC', iret)
                if (iret == 0) then
                    repi = 0
                end if
            else
            !! Case MATR_ELEM
                if (rerr(3) (1:3) == "NON") then
                    repi = 0
                end if
            end if
        end if
    else
        ierd = 1
    end if
!
999 continue
    repkz = repk
    call jedema()
end subroutine
