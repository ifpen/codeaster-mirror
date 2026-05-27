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
subroutine mcopco(mesh, newgeo, cellNume, ksi1, &
                  ksi2, geom)
!
    implicit none
!
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mmelty.h"
#include "asterfort/mmvalp.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: newgeo
    integer(kind=8), intent(in) :: cellNume
    real(kind=8), intent(in) :: ksi1, ksi2
    real(kind=8), intent(out) :: geom(3)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE CONTINUE - APPARIEMENT - UTILITAIRE)
!
! CALCUL DES COORDONNEES D'UN NOEUD SUR UNE MAILLE A PARTIR
! DE SES COORDONNEES PARAMETRIQUES
!
! --------------------------------------------------------------------------------------------------
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  NEWGEO : COORDONNEES DE TOUS LES NOEUDS
! IN  NUMMAI : NUMERO ABSOLU DE LA MAILLE DANS LE MAILLAGE
! IN  KSI1   : COORDONNEE PARAMETRIQUE KSI DU PROJETE
! IN  KSI2   : COORDONNEE PARAMETRIQUE ETA DU PROJETE
! OUT GEOM   : COORDONNEES DU NOEUD
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) ::  jdes
    integer(kind=8) :: cellNbNode, iNode, nodeNume(9)
    real(kind=8) :: valeCell(27)
    character(len=8) :: cellCode
    real(kind=8), pointer :: vale(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    geom = 0.d0

! - ACCES AUX CHAMPS
    call jeveuo(newgeo(1:19)//'.VALE', 'L', vr=vale)
    call jeveuo(jexnum(mesh//'.CONNEX', cellNume), 'L', jdes)

! - INFOS SUR LA MAILLE
    call mmelty(mesh, cellNume, cellCode, cellNbNode)

! - NUMEROS ABSOLUS DES NOEUDS DE LA MAILLE
    do iNode = 1, cellNbNode
        nodeNume(iNode) = zi(jdes+iNode-1)
    end do

! - COORDONNEES DES NOEUDS DE LA MAILLE
    do iNode = 1, cellNbNode
        valeCell(3*(iNode-1)+1) = vale(1+3*(nodeNume(iNode)-1))
        valeCell(3*(iNode-1)+2) = vale(1+3*(nodeNume(iNode)-1)+1)
        valeCell(3*(iNode-1)+3) = vale(1+3*(nodeNume(iNode)-1)+2)
    end do

! - CALCUL DES COORDONNEES
    call mmvalp(cellCode, cellNbNode, ksi1, ksi2, valeCell, geom)
!
    call jedema()
end subroutine
