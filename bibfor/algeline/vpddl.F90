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

subroutine vpddl(raide, masse, neq, nblagr, nbcine, &
                 neqact, dlagr, dbloq, ier)
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/pteddl.h"
#include "asterfort/typddl.h"
#include "asterfort/utmess.h"
!
    character(len=19), intent(in) :: masse, raide
    integer(kind=8), intent(in) :: neq
    integer(kind=8), intent(out) :: nblagr, nbcine, neqact, dlagr(neq), dbloq(neq), ier
!
!     ------------------------------------------------------------------
!     RENSEIGNEMENTS SUR LES DDL : LAGRANGE, BLOQUE, EXCLUS.
!     CONSTRUCTION DE TABLEAUX D'ENTIERS REPERANT LA POSITION DE CES DDL
!     ------------------------------------------------------------------
! IN  RAIDEUR : K  : NOM DE LA MATRICE DE "RAIDEUR"
! IN  MASSE   : K  : NOM DE LA MATRICE DE "MASSE"
! IN  NEQ     : IS : NPMBRE DE DDL
! OUT NBLAGR  : IS : NOMBRE DE DDL DE LAGRANGE
! OUT NBCINE  : IS : NOMBRE DE DDL BLOQUES PAR AFFE_CHAR_CINE
! OUT NEQACT  : IS : NOMBRE DE DDL ACTIFS
! OUT DLAGR   : IS : POSITION DES DDL DE LAGRANGE
! OUT DBLOQ   : IS : POSITION DES DDL BLOQUES PAR AFFE_CHAR_CINE
!
!
!
    integer(kind=8) :: iercon, nbprno, ieq, nba, nbb, nbl, nbliai, ifm, niv
    integer(kind=8) :: typelagr
    integer(kind=8) :: vali(4)
    character(len=14) :: nume
    integer(kind=8), pointer :: ccid(:) => null()
    character(len=24), pointer :: refn(:) => null()
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
    call jemarq()
!
!     ---RECUPERATION DU NIVEAU D'IMPRESSION---
    call infniv(ifm, niv)
!     -----------------------------------------
!
!     --- CALCUL DU NOMBRE DE LAGRANGES ---
!     -------------------------------------
!
!       --- RECUPERATION DU NOM DE LA NUMEROTATION ASSOCIEE AUX MATRICES
    call dismoi('NOM_NUME_DDL', raide, 'MATR_ASSE', repk=nume)
    call jeveuo(nume//'.NUME.REFN', 'L', vk24=refn)
    if (refn(4) .eq. "SIMPLE_LAGRANGE") then
        typelagr = 1
    else
        typelagr = 2
    end if
!
!       --- RECUPERATION DES POSITIONS DES DDL LAGRANGE : DLAGR
    call pteddl('NUME_DDL', nume, 1, 'LAGR    ', neq, &
                list_equa=dlagr)
!
!       --- CALCUL DU NOMBRE DE 'LAGRANGE': NBLAGR
    nblagr = 0
    do ieq = 1, neq
        nblagr = nblagr+dlagr(ieq)
    end do
!
!       --- INVERSION : DLAGR = 0 SI LAGRANGE ET 1 SINON
    do ieq = 1, neq
        dlagr(ieq) = abs(dlagr(ieq)-1)
    end do
!
!     --- DETECTION DES DDL BLOQUES PAR AFFE_CHAR_CINE ---
!     ----------------------------------------------------
!
    call typddl('ACLA', nume, neq, dbloq, nba, &
                nbb, nbl, nbliai)
!   blindage pour ELIM_LAGR (nba sort à 0 dans ce cas)
    if (nblagr .eq. 0) then
        neqact = neq
!
    else
        neqact = neq-(typelagr+1)*nblagr/typelagr
    end if
!
!       --- MISE A JOUR DE DBLOQ QUI VAUT 0 POUR TOUS LES DDL BLOQUES
    call jeexin(masse//'.CCID', iercon)
    nbcine = 0
    if (iercon .ne. 0) then
        call jeveuo(masse//'.CCID', 'E', vi=ccid)
        do ieq = 1, neq
            dbloq(ieq) = dbloq(ieq)*abs(ccid(ieq)-1)
        end do
!
!       --- CALCUL DU NOMBRE DE DDL BLOQUE PAR CETTE METHODE : NCINE ---
        do ieq = 1, neq
            nbcine = nbcine+ccid(ieq)
        end do
    end if
    neqact = neqact-nbcine
!   vérification sur le nombre de noeuds actifs
    if (neqact .le. 0) then
        call utmess('F', 'ALGELINE3_63')
    end if
!
!     --- SI NUMEROTATION GENERALISEE : PAS DE DDLS BLOQUES ---
!     ---------------------------------------------------------
    call jenonu(jexnom(nume//'.NUME.LILI', '&SOUSSTR'), nbprno)
    if (nbprno .ne. 0) then
        do ieq = 1, neq
            dbloq(ieq) = 1
        end do
    end if
!
!    -----IMPRESSION DES DDL -----
!
    if (niv .ge. 1) then
        vali(1) = neq
        vali(2) = nblagr
        if (nbcine .eq. 0) then
            vali(3) = neqact
            call utmess('I', 'ALGELINE7_17', ni=3, vali=vali)
        else
            vali(3) = nbcine
            vali(4) = neqact
            call utmess('I', 'ALGELINE7_18', ni=4, vali=vali)
        end if
    end if
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!
    ier = 0
    call jedema()
!
end subroutine
