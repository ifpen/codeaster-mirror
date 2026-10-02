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
   subroutine global_numbering_info(numddl)
!
       implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
#include "asterfort/jeexin.h"
!
       character(len=14) :: numddl
!--------------------------------------
!       Print info on matrix for HPC
! --------------------------------------
!
       integer(kind=8) :: ifm, niv, nequ_owned, ncoef_owned, iddl, nequ, rang
       integer(kind=8) :: vali(3), iret
       real(kind=8) :: valr(2), rcoef, requ
       mpi_int :: mrank
       integer(kind=8), pointer :: v_posdd(:) => null()
       integer(kind=8), pointer :: v_nequ(:) => null()
       integer(kind=8), pointer :: v_smdi(:) => null()
!
       call infniv(ifm, niv)
       call jemarq()

       call jeexin(numddl//'NUME.SMDI', iret)
!
       if (niv .ge. 1 .and. iret .ne. 0) then
! Comptage en HPC
           call asmpi_info(rank=mrank)
           rang = to_aster_int(mrank)
!
! Compteurs locaux limites aux lignes proprietaires.
           nequ_owned = 0
           ncoef_owned = 0

           call jeveuo(numddl//'.NUME.PDDL', 'L', vi=v_posdd)
           call jeveuo(numddl//'.NUME.NEEQ', 'L', vi=v_nequ)
           call jeveuo(numddl//'.NUME.SMDI', 'L', vi=v_smdi)
           nequ = v_nequ(1)

           do iddl = 1, nequ
               if (v_posdd(iddl) .ne. rang) cycle

               nequ_owned = nequ_owned+1

! Dans cette representation, SMDI(iddl) pointe sur la diagonale :
! la difference avec la diagonale precedente donne le nombre
! de termes stockes pour cette ligne.
               if (iddl == 1) then
                   ncoef_owned = ncoef_owned+v_smdi(1)
               else
                   ncoef_owned = ncoef_owned+v_smdi(iddl)-v_smdi(iddl-1)
               end if
           end do

           vali(1:2) = [nequ_owned, ncoef_owned]
           call asmpi_comm_vect("MPI_SUM", "I", 2, vi=vali)
!
           vali(3) = 2*vali(2)-vali(1)
           rcoef = vali(2)
           requ = vali(1)
           valr(1) = (100.d0*(2.d0*rcoef-requ))/(requ*requ)
           call utmess('I', 'FACTOR_2', ni=3, vali=vali, sr=valr(1))
       end if
!
       call jedema()
!
   end subroutine global_numbering_info
