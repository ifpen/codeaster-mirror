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

subroutine ssrige(macrElem)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/assmam.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/merime.h"
#include "asterfort/numddl.h"
#include "asterfort/promor.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/sdmpic.h"
#include "asterfort/smosli.h"
#include "asterfort/ssriu1.h"
#include "asterfort/ssriu2.h"
!
    character(len=8), intent(in) :: macrElem
! ----------------------------------------------------------------------
!     BUT: TRAITER LE MOT CLEF "RIGI_MECA" DE L'OPERATEUR MACR_ELEM_STAT
!          CALCULER LA MATRICE DE RIGIDITE CONDENSEE DU MACR_ELEM_STAT.
!
!     IN: NOMU   : NOM DU MACR_ELEM_STAT
!
!     OUT: LES OBJETS SUIVANTS DU MACR_ELEM_STAT SONT CALCULES:
!          .PHI_IE ET .KP_EE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: renumRCMK = "RCMK"
    integer(kind=8) :: nbLoad
    integer(kind=8), parameter :: nbMatrElem = 1
    character(len=24), pointer :: listMatrElem(:) => null()
    real(kind=8) :: rtbloc
    character(len=1), parameter :: base = 'V'
    character(len=19), parameter :: matrElem = '&&MATEL'
    character(len=19), parameter :: comporMult = ' '
    character(len=8) :: model, caraElem, mate
    character(len=14) :: numeDof
    character(len=19) :: matrAsse
    character(len=24) :: mateco
    real(kind=8) :: time
    integer(kind=8), parameter :: modeFourier = 0
    integer(kind=8), pointer :: desm(:) => null()
    real(kind=8), pointer :: varm(:) => null()
    character(len=8), pointer :: refm(:) => null()
    character(len=24), pointer :: listLoadK24(:) => null()
    character(len=24) :: modelLigrel
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

    numeDof = macrElem
    matrAsse = macrElem//'.RIGIMECA'
!
    call jeveuo(macrElem//'.DESM', 'L', vi=desm)
    nbLoad = desm(6)
    ASSERT(nbLoad .le. 1)
!
    call jeveuo(macrElem//'.REFM', 'E', vk8=refm)
    model = refm(1)
    caraElem = refm(4)
    mate = refm(3)
    if (mate .eq. ' ') then
        mateco = ' '
    else
        call rcmfmc(mate, mateco, l_ther_=ASTER_FALSE)
    end if
!
    call jeveuo(macrElem//'.VARM', 'L', vr=varm)
    time = varm(2)

    if (nbLoad .ne. 0) then
        ASSERT(nbLoad .eq. 1)
        AS_ALLOCATE(vk24=listLoadK24, size=nbLoad)
        listLoadK24(1) = refm(10)
    end if

! - Elementary matrices for rigidity
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call merime(model, nbLoad, listLoadK24, mate, mateco, caraElem, &
                time, comporMult, matrElem, modeFourier, base, modelLigrel)

! - Get list of LIGREL from elementary matrices
    AS_ALLOCATE(vk24=listMatrElem, size=nbMatrElem)
    listMatrElem(1) = matrElem

! - Create numbering from elementary matrix
    call numddl(numeDof, renumRCMK, 'GG', nbMatrElem, listMatrElem, &
                lMacrElemZ_=ASTER_TRUE)
    AS_DEALLOCATE(vk24=listMatrElem)
!
!   -- ON MET LES DDLS INTERNES AVANT LES EXTERNES
!      AVANT DE CONSTRUIRE LE PROFIL :
    call ssriu1(macrElem)
    call promor(numeDof, 'G')
    rtbloc = varm(1)
    call smosli(numeDof//'.SMOS', numeDof//'.SLCS', 'G', rtbloc)
!
!   -- ASSEMBLAGE:
    call assmam('G', matrAsse, 1, matrElem, [1.d0], numeDof, 'ZERO', 1)

!   -- IL FAUT COMPLETER LA MATRICE SI LES CALCULS SONT DISTRIBUES:
    call sdmpic('MATR_ASSE', matrAsse)
!
    call ssriu2(macrElem)
!
!   -- MISE A JOUR DE .REFM(5) ET REFM(6)
    refm(5) = macrElem
    refm(6) = 'OUI_RIGI'
!
!
    call jedetr(macrElem//'      .NEWN')
    call jedetr(macrElem//'      .OLDN')
    call jedetr(numeDof//'     .ADNE')
    call jedetr(numeDof//'     .ADLI')
    call jedetr(matrAsse(1:19)//'.LILI')
    call detrsd('MATR_ELEM', matrElem)
    AS_DEALLOCATE(vk24=listLoadK24)
!
    call jedema()
end subroutine
