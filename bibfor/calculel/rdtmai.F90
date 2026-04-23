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

subroutine rdtmai(noma, nomare, base, corrn, corrm, bascor)
!
    use crea_maillage_module
!
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cargeo.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/juveca.h"
#include "asterfort/reliem.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=8) :: noma, nomare
    character(len=*) :: corrn, corrm
    character(len=1) :: base, bascor
!
!
! ======================================================================
!     BUT: REDUIRE UN MAILLAGE SUR UNE LISTE DE MAILLES
!
!  NOMA : IN  : MAILLAGE A REDUIRE
!  NOMARE : OUT : MAILLAGE REDUIT
!  BASE   : IN  : 'G' OU 'V' : BASE POUR LA CREATION DE NOMARE
!  CORRN  : IN/JXOUT : SI != ' ' : NOM DE L'OBJET QUI CONTIENDRA
!           LA CORRESPONDANCE INO_RE -> INO
!  CORRM  : IN/JXOUT : SI != ' ' : NOM DE L'OBJET QUI CONTIENDRA
!           LA CORRESPONDANCE IMA_RE -> IMA
!  BASCOR : IN  : 'G' OU 'V' : BASE POUR LA CREATION DE CORRN ET CORRM
! ======================================================================
!
    integer(kind=8) :: nbmaou, nbnoin, jnuma, jwk1, jconx2, ima, numa
    integer(kind=8) :: nbno, ino, nuno, iret, nbgma, jgma, nbgno
    integer(kind=8) :: jcorrm, n1, nbmain, jwk2, nbnoou, jmaor
    character(len=8) :: typmcl(2), ttgr
    character(len=16) :: motcle(2)
    character(len=24) :: grpname, grpma
    aster_logical :: lpmesh
    integer(kind=8), pointer :: vconnex(:) => null()
    integer(kind=8), pointer :: num_noeu_in(:) => null()
    type(Mmesh) :: mesh_res
!
    call jemarq()
!
    ASSERT(noma .ne. nomare)
    ASSERT(base .eq. 'V' .or. base .eq. 'G')
!
!
! -1- PRELIMINAIRES
!     ============
!
    lpmesh = isParallelMesh(noma)
!
! - Create new mesh
!
    call mesh_res%init(noma, 1, ASTER_FALSE)
!
!
! --- CALCUL DE LA LISTE DES MAILLES SUR LESQUELLES IL FAUT REDUIRE :
    motcle(1) = 'GROUP_MA'
    motcle(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
    call reliem(' ', noma, 'NU_MAILLE', 'RESTREINT', 1, &
                2, motcle, typmcl, '&&RDTMAI.NUM_MAIL_IN', nbmaou)
    if (nbmaou > 0) then
        call jeveuo('&&RDTMAI.NUM_MAIL_IN', 'L', jnuma)
    end if
!
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnoin)
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbmain)
!
! --- CREATION DE TABLEAUX DE TRAVAIL:
!     ZI(JWK1) :
!     - DIMENSIONNE AU NOMBRE DE NOEUDS DU MAILLAGE IN
!     - CORRESPONDANCE : NUMEROS DES NOEUDS MAILLAGE IN => MAILLAGE OUT
!     - EX: ZI(JWK1+INO1-1)=INO2
!         -> SI INO2!=0:LE NOEUD INO1 DU MAILLAGE IN CORRESPOND AU NOEUD
!                       INO2 DU MAILLAGE OUT.
!         -> SI INO2=0: LE NOEUD INO1 DU MAILLAGE IN N'EST PAS PRESENT
!                       DANS LE MAILLAGE OUT.
!
    call wkvect('&&RDTMAI_WORK_1', 'V V I', nbnoin, jwk1)
!
!     ZI(JWK2) : (L'INVERSE DE ZI(JWK1))
!     - DIMENSIONNE AU NOMBRE DE NOEUDS DU MAILLAGE IN
!     - CORRESPONDANCE : NUMEROS DES NOEUDS MAILLAGE OUT => MAILLAGE IN
!     - EX: ZI(JWK1+INO1-1)=INO2
!        -> LE NOEUD INO1 DU MAILLAGE OUT CORRESPOND AU NOEUD
!           INO2 DU MAILLAGE IN.
    call wkvect('&&RDTMAI_WORK_2', 'V V I', nbnoin, jwk2)
!
!
! ---  REMPLISSAGE DES TABLEAUX DE TRAVAIL
    call jeveuo(noma//'.CONNEX', 'L', vi=vconnex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
    nbnoou = 0
    do ima = 1, nbmaou
        numa = zi(jnuma+ima-1)
        nbno = zi(jconx2+numa)-zi(jconx2+numa-1)
        do ino = 1, nbno
            nuno = vconnex(zi(jconx2+numa-1)+ino-1)
            if (zi(jwk1+nuno-1) .eq. 0) then
                nbnoou = nbnoou+1
                zi(jwk1+nuno-1) = nbnoou
                zi(jwk2+nbnoou-1) = nuno
            end if
        end do
    end do
!
!   -- il faut ajouter les noeuds demandes par l'utlisateur :
    call reliem(' ', noma, 'NU_NOEUD', 'RESTREINT', 1, &
                1, ['GROUP_NO'], ['GROUP_NO'], '&&RDTMAI.NUM_NOEU_IN', n1)
    if (n1 .gt. 0) then
        call jeveuo('&&RDTMAI.NUM_NOEU_IN', 'L', vi=num_noeu_in)
        do ino = 1, n1
            nuno = num_noeu_in(ino)
            if (zi(jwk1+nuno-1) .eq. 0) then
                nbnoou = nbnoou+1
                zi(jwk1+nuno-1) = nbnoou
                zi(jwk2+nbnoou-1) = nuno
            end if
        end do
    end if
!
! - Trier les noeuds pour les avoir dans le même ordre que le maillage initial
!
    nbnoou = 0
    do ino = 1, nbnoin
        if (zi(jwk1+ino-1) > 0) then
            nbnoou = nbnoou+1
            zi(jwk1+ino-1) = nbnoou
            zi(jwk2+nbnoou-1) = ino
        else
            ! remove this node
            mesh_res%nodes(ino)%keep = ASTER_FALSE
        end if
    end do
!
! - Remove cells
!
    do ima = 1, nbmain
        mesh_res%cells(ima)%keep = ASTER_FALSE
    end do
!
    do ima = 1, nbmaou
        mesh_res%cells(zi(jnuma+ima-1))%keep = ASTER_TRUE
    end do
!
! - Update numbering
!
    call mesh_res%update()
!
! - Copy mesh
!
    call mesh_res%copy_mesh(nomare)
!
! - Remove groups if needed
!
    call getvtx('RESTREINT', 'TOUT_GROUP_MA', iocc=1, scal=ttgr, nbret=iret)
    if (ttgr .eq. 'NON') then
!       'TOUT_GROUP_MA'='NON'
        grpma = nomare//'.GROUPEMA       '
        grpname = nomare//'.PTRNOMMAI      '
        call jedetr(grpma)
        call jedetr(grpname)
        call getvtx('RESTREINT', 'GROUP_MA', iocc=1, nbval=0, nbret=nbgma)
        nbgma = -nbgma
        if (nbgma > 0) then
            call wkvect('&&RDTMAI_GRMA_FOURNIS', 'V V K24', nbgma, jgma)
            call getvtx('RESTREINT', 'GROUP_MA', iocc=1, nbval=nbgma, vect=zk24(jgma), &
                        nbret=iret)
            call mesh_res%copy_group_ma(grpma, grpname, nbgma, zk24(jgma))
            call jedetr('&&RDTMAI_GRMA_FOURNIS')
        end if
    end if
!
    call getvtx('RESTREINT', 'TOUT_GROUP_NO', iocc=1, scal=ttgr, nbret=iret)
    if (ttgr .eq. 'NON') then
        !       'TOUT_GROUP_NO'='NON'
        grpma = nomare//'.GROUPENO       '
        grpname = nomare//'.PTRNOMNOE      '
        call jedetr(grpma)
        call jedetr(grpname)
        call getvtx('RESTREINT', 'GROUP_NO', iocc=1, nbval=0, nbret=nbgno)
        nbgno = -nbgno
        if (nbgno > 0) then
            call wkvect('&&RDTMAI_GRNO_FOURNIS', 'V V K24', nbgno, jgma)
            call getvtx('RESTREINT', 'GROUP_NO', iocc=1, nbval=nbgno, vect=zk24(jgma), &
                        nbret=iret)
            call mesh_res%copy_group_no(grpma, grpname, nbgno, zk24(jgma))
            call jedetr('&&RDTMAI_GRNO_FOURNIS')
        end if
    end if
!
! - Cleaning
!
    call mesh_res%clean()
!
    call cargeo(nomare)
!
!
    if (base .eq. 'G') then
        call wkvect(nomare//'.MAOR', 'G V K8', 1, jmaor)
        zk8(jmaor) = noma
    end if
!     -- SI L'ON SOUHAITE RECUPERER LES TABLEAUX DE CORRESPONDANCE :
    if (corrn .ne. ' ') then
        if (nbnoou > 0) then
            call juveca('&&RDTMAI_WORK_2', nbnoou)
            call jedupo('&&RDTMAI_WORK_2', bascor, corrn, .false._1)
        end if
    end if
    if (corrm .ne. ' ') then
        if (nbmaou > 0) then
            call wkvect(corrm, bascor//' V I', nbmaou, jcorrm)
            do ima = 1, nbmaou
                zi(jcorrm-1+ima) = zi(jnuma+ima-1)
            end do
        end if
    end if
!
    call jedetr('&&RDTMAI_WORK_1')
    call jedetr('&&RDTMAI_WORK_2')
!
    call jedema()
!
end subroutine
