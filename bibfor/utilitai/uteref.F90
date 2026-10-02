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

subroutine uteref(chanom, typech, tyelas, nomte, lfichUniq, &
                  nomfpg, nnos, nno, nbpg, ndim, &
                  refcoo, gscoo, wg, nochmd, codret)
!
!-----------------------------------------------------------------------
!     UTILITAIRE - ELEMENT DE REFERENCE
!     --           -          ---
!-----------------------------------------------------------------------
!     ENTREES :
!       CHANOM : NOM ASTER DU CHAMP
!       TYPECH : TYPE DU CHAMP ('ELGA')
!       TYELAS : TYPE D'ELEMENT ASTER
!       NOMTE  : NOM DE L'ELEMENT FINI A EXAMINER
!     SORTIES :
!       NOMFPG : NOM DE LA FAMILLE DES POINTS DE GAUSS
!       NNOS   : NOMBRE DE NOEUDS SOMMETS
!       NNO    : NOMBRE DE NOEUDS TOTAL
!       NBPG   : NOMBRE DE POINTS DE GAUSS
!       NDIM   : DIMENSION DE L'ELEMENT
!       REFCOO : COORDONNEES DES NNO NOEUDS
!       GSCOO  : COORDONNEES DES POINTS DE GAUSS, SI CHAMP ELGA
!       WG     : POIDS DES POINTS DE GAUSS, SI CHAMP ELGA
!       CODRET : CODE DE RETOUR
!                0 : PAS DE PB
!                1 : LE CHAMP N'EST PAS DEFINI SUR CE TYPE D'ELEMENT
!     REMARQUE :
!     ON DOIT RETOURNER LES COORDONNEES SOUS LA FORME :
!     . ELEMENT 1D : X1 X2 ... ... XN
!     . ELEMENT 2D : X1 Y1 X2 Y2 ... ... XN YN
!     . ELEMENT 3D : X1 Y1 Z1 X2 Y2 Z2 ... ... XN YN ZN
!     C'EST CE QUE MED APPELLE LE MODE ENTRELACE
!     ON DOIT RETOURNER LES POIDS SOUS LA FORME :
!     WG1 WG2 ... ... WGN
!-----------------------------------------------------------------------
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/indik8.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/elraca.h"
#include "asterfort/elraga.h"
#include "asterfort/elref2.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/uteref2.h"
#include "MeshTypes_type.h"
!
    integer(kind=8) :: tyelas
    integer(kind=8) :: nnos, nno, nbpg, ndim
!
    real(kind=8) :: refcoo(*), gscoo(*), wg(*)
!
    character(len=8) :: typech
    character(len=16) :: nomte
    character(len=16) :: nomfpg
    character(len=19) :: chanom
    character(len=64) :: nochmd
    aster_logical :: lfichUniq
!
    integer(kind=8) :: codret
!
!
! 0.3. ==> VARIABLES LOCALES
!
    integer(kind=8) :: imolo
    integer(kind=8) :: itype, nb1, dimrefe
    integer(kind=8) :: ifm, nivinf, igrel, ntypef
    integer(kind=8) :: ierd, jliel, nbgrel
    integer(kind=8) :: iaux, vimolo(1)
!
    character(len=4) :: tych
    character(len=8) :: elrefe, famil, nomtypmail
    character(len=16) :: nomsym
    character(len=64) :: valk(3)
    character(len=19) :: ligrel, resu
!
    integer(kind=8), pointer :: celd(:) => null()
    character(len=24), pointer :: celk(:) => null()
!-----------------------------------------------------------------------
!     1- PREALABLES
!     ---------------
    call jemarq()
    codret = 0
!
    call infniv(ifm, nivinf)
    if (nivinf .gt. 1) then
        call utmess('I', 'UTILITAI5_39')
        write (ifm, 10) tyelas, nomte
10      format('ELEMENT FINI NUMERO', i6, ', DE NOM : ', a16)
    end if
    ASSERT(typech .eq. 'ELGA')
!
!     2- DETERMINATION DE ELREFE :
!     -----------------------------
!
    if (codret .eq. 0) then
!
        call dismoi('TYPE_CHAMP', chanom, 'CHAMP', repk=tych)
        ASSERT(tych .eq. typech)
        call jeveuo(chanom//'.CELK', 'L', vk24=celk)
        call jeveuo(chanom//'.CELD', 'L', vi=celd)
        ligrel = celk(1) (1:19)
        ASSERT(celk(3) (1:19) .eq. typech)
        nbgrel = celd(2)
!
    end if
!
!     3- DETERMINATION DU MODE_LOCAL (IMOLO) ASSOCIE AU CHAMP POUR
!        LE TYPE D'ELEMENT TYELAS :
!     -------------------------------------------------------------
!
    if (codret .eq. 0) then
        imolo = 0
        do igrel = 1, nbgrel
            call jeveuo(jexnum(ligrel//'.LIEL', igrel), 'L', jliel)
            call jelira(jexnum(ligrel//'.LIEL', igrel), 'LONMAX', nb1)
            itype = zi(jliel-1+nb1)
            if (itype .eq. tyelas) then
                imolo = celd(celd(4+igrel)+2)
                if (imolo .eq. 0) then
                    codret = 1
                    if (nivinf .gt. 1) then
                        write (ifm, *)&
         &          '==> LE CHAMP N''EST PAS DEFINI SUR CE TYPE D''ELEMENT'
                    end if
                end if
                goto 32
            end if
        end do
!
        if (.not. lfichUniq) then
            ASSERT(.false.)
        end if
!
32      continue

        if (lfichUniq) then
            vimolo(1) = imolo
            call asmpi_comm_vect('MPI_MAX', 'I', 1, vi=vimolo)
            imolo = vimolo(1)
        end if
!
    end if

    resu = chanom(1:8)
    call jeexin(resu//'.DESC', ierd)
    if (ierd .eq. 0) then
        nomsym = nochmd(1:16)
    else
        nomsym = chanom(1:16)
    end if

!
!     4- DETERMINATION DE LA FAMILLE DE POINTS DE GAUSS :
!     -------------------------------------------------------------
!
    if (codret .eq. 0) then
        call jenonu(jexnom('&CATA.TE.NOMTE', nomte), ntypef)
        call uteref2(ntypef, imolo, ASTER_TRUE, &
                     elrefe, famil, nomtypmail, nnos, nno, nbpg, ndim, dimrefe, &
                     refcoo, gscoo, wg, codret)

        nomfpg(1:8) = elrefe
        nomfpg(9:16) = famil

        if (codret .eq. 1) then
            call utmess("F", "MED2_11")
        else if (codret .eq. 2) then
            valk(1) = nomsym
            valk(2) = famil
            valk(3) = nomtypmail
            call utmess('F', 'MED2_5', nk=3, valk=valk)
        else if (codret .eq. 3) then
            call utmess('F', 'MED2_6', valk=nomsym)
        else
            ASSERT(codret .eq. 0)
        end if

    end if
!
!     6- IMPRESSION EVENTUELLE SUR LE FICHIER DE MESSAGES
!     ----------------------------------------------------------------
!
    if (codret .eq. 0) then
!
        if (nivinf .gt. 1) then
!
            write (ifm, 801) 'FAMILLE DE POINTS DE GAUSS', nomfpg
            write (ifm, 802) 'NOMBRE DE SOMMETS        ', nnos
            write (ifm, 802) 'NOMBRE DE NOEUDS         ', nno
            write (ifm, 802) 'NOMBRE DE POINTS DE GAUSS', nbpg
!
!     6.1. DIMENSION 1
!
            if (ndim .eq. 1) then
!                            123456789012345
                write (ifm, 601) 'NOEUDS         '
                do iaux = 1, nno
                    write (ifm, 711) iaux, refcoo(iaux)
                end do
                write (ifm, 721)
                write (ifm, 601) 'POINTS DE GAUSS'
                do iaux = 1, nbpg
                    write (ifm, 711) iaux, gscoo(iaux)
                end do
                write (ifm, 721)
!
!     6.2. DIMENSION 2
!
            else if (ndim .eq. 2) then
                write (ifm, 602) 'NOEUDS         '
                do iaux = 1, nno
                    write (ifm, 712) iaux, refcoo(ndim*(iaux-1)+1), &
                        refcoo(ndim*(iaux-1)+2)
                end do
                write (ifm, 722)
                write (ifm, 602) 'POINTS DE GAUSS'
                do iaux = 1, nbpg
                    write (ifm, 712) iaux, gscoo(ndim*(iaux-1)+1), &
                        gscoo(ndim*(iaux-1)+2)
                end do
                write (ifm, 722)
!
!     6.3. DIMENSION 3
!
            else
                write (ifm, 603) 'NOEUDS         '
                do iaux = 1, nno
                    write (ifm, 713) iaux, refcoo(ndim*(iaux-1)+1), &
                        refcoo(ndim*(iaux-1)+2), refcoo(ndim*(iaux-1)+3)
                end do
                write (ifm, 723)
                write (ifm, 603) 'POINTS DE GAUSS'
                do iaux = 1, nbpg
                    write (ifm, 713) iaux, gscoo(ndim*(iaux-1)+1), &
                        gscoo(ndim*(iaux-1)+2), gscoo(ndim*(iaux-1)+3)
                end do
                write (ifm, 723)
            end if
!
            write (ifm, 604)
            do iaux = 1, nbpg
                write (ifm, 711) iaux, wg(iaux)
            end do
            write (ifm, 721)
!
        end if
!
    end if
!
601 format(&
     &/, 28('*'),&
     &/, '*      COORDONNEES DES     *',&
     &/, '*      ', a15, '     *',&
     &/, 28('*'),&
     &/, '*  NUMERO  *       X       *',&
     &/, 28('*'))
602 format(&
     &/, 44('*'),&
     &/, '*       COORDONNEES DES ', a15, '    *',&
     &/, 44('*'),&
     &/, '*  NUMERO  *       X       *       Y       *',&
     &/, 44('*'))
603 format(&
     &/, 60('*'),&
     &/, '*            COORDONNEES DES ', a15,&
     &'               *',&
     &/, 60('*'),&
     &/, '*  NUMERO  *       X       *       Y       *',&
     &'       Z       *',&
     &/, 60('*'))
604 format(&
     &/, 28('*'),&
     &/, '*      POINTS DE GAUSS     *',&
     &/, '*  NUMERO  *     POIDS     *',&
     &/, 28('*'))
711 format('* ', i5, '    *', 1pg12.5, '    *')
712 format('* ', i5, 2('    *', 1pg12.5), '    *')
713 format('* ', i5, 3('    *', 1pg12.5), '    *')
721 format(28('*'))
722 format(44('*'))
723 format(60('*'))
801 format(a, ' : ', a)
802 format(a, ' : ', i4)
!
    call jedema()
end subroutine
