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

subroutine global_numbering_communicate(nume_equa)
    implicit none
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_recv_i.h"
#include "asterc/asmpi_send_i.h"
#include "asterc/asmpi_sendrecv_i.h"
#include "asterf_config.h"
#include "asterf_debug.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/create_graph_comm.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/isdeco.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/jeecra.h"
#include "asterfort/jecrec.h"
#include "asterfort/nbec.h"
#include "asterfort/parallel_ligrel_list.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=19) :: nume_equa

#ifdef ASTER_HAVE_MPI
!
    integer(kind=8) :: rang, nbproc, jrefn, iaux, nddll, i_join
    integer(kind=8) :: nb_comm, icmp, ico2, nbcmp
    integer(kind=8) :: idprn1, idprn2, ili, ili2, neql
    integer(kind=8) :: nec, numpro, jjoine, jjoinr, nbnoee, jaux, numno1, numno2, iec
    integer(kind=8) :: ncmpmx, iad, jcpnec, jencod, jenvoi1, lgenve1, lgenvr1, poscom
    integer(kind=8) :: nbddll, jnequ, nddl, jenco2, jcpne2, numpr2
    integer(kind=8) :: nbnoer, jrecep1, curpos, nbLigrTot, iLigrT
    integer(kind=8) :: jmlogl, nuno, nb_ddl_envoi, nbddl, domj_i
    integer(kind=8) :: nddlg, jenvoi2, jrecep2, numnoe, gd, iddl
    integer(kind=8) :: ifm, niv, vali(5), ino, nno, nb_node, nlag, exipr
    integer(kind=8) :: lgenve2, lgenvr2, jnujoi1, jnujoi2, iret, iret1, iret2, nlili, nlili2
    mpi_int :: mrank, mnbproc, mpicou, tag4, numpr4, n4e, n4r
    integer(kind=8), pointer :: v_noex(:) => null()
    integer(kind=8), pointer :: v_nugll(:) => null()
    integer(kind=8), pointer :: v_posdd(:) => null()
    integer(kind=8), pointer :: v_deeq(:) => null()
    integer(kind=8), pointer :: v_deeg(:) => null()
    integer(kind=8), pointer :: v_nuls(:) => null()
    integer(kind=8), pointer :: v_comm(:) => null()
    integer(kind=8), pointer :: v_tag(:) => null()
    integer(kind=8), pointer :: v_dom(:) => null()
    integer(kind=8), pointer :: v_gco(:) => null()
    integer(kind=8), pointer :: v_refp(:) => null()
    integer(kind=4), pointer :: v_pgid(:) => null()
    integer(kind=8), pointer :: v_lilt(:) => null()
    aster_logical :: lligrel_cp
!
    character(len=8) :: mesh, k8bid, nomgdr
    character(len=19) :: nomlig, comm_name, tag_name, joints, meshj
    character(len=24) :: nonulg, domj, recv, send, name, gcom, pgid, lilt
    character(len=32) :: nojoie, nojoir
    character(len=24), pointer :: tco(:) => null()
!
!----------------------------------------------------------------------
!
!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS PRNO DES S.D. LIGREL
!     REPERTORIEES DANS LE CHAMP LILI DE NUME_DDL ET A LEURS ADRESSES
!     ZZPRNO(ILI,NUNOEL,1) = NUMERO DE L'EQUATION ASSOCIEES AU 1ER DDL
!                            DU NOEUD NUNOEL DANS LA NUMEROTATION LOCALE
!                            AU LIGREL ILI DE .LILI
!     ZZPRNO(ILI,NUNOEL,2) = NOMBRE DE DDL PORTES PAR LE NOEUD NUNOEL
!     ZZPRNO(ILI,NUNOEL,2+1) = 1ER CODE
!     ZZPRNO(ILI,NUNOEL,2+NEC) = NEC IEME CODE
!
#define zzprno(ili, nunoel, l) zi(idprn1 - 1 + zi(idprn2 + ili - 1) + (nunoel - 1)*(nec + 2) + l-1)
!
    call jemarq()
!
    call infniv(ifm, niv)
!
    call asmpi_comm('GET', mpicou)
    call asmpi_info(rank=mrank, size=mnbproc)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(mnbproc)
    DEBUG_MPI('crnggc', rang, nbproc)

!   RECUPERATION DU NOM DU MAILLAGE DANS LE BUT D'OBTENIR LE JOINT
    call jeveuo(nume_equa//'.REFN', 'L', jrefn)
    mesh = zk24(jrefn)
    nomgdr = zk24(jrefn+1)
    call parallel_ligrel_list(nume_equa, 'G')
!
!   SI LE MAILLAGE N'EST PAS PARALLELE, ON SORT
    if (.not. isParallelMesh(mesh)) goto 999
    call jeveuo(mesh//'.NOEX', 'L', vi=v_noex)
!
    call dismoi('NB_EQUA', nume_equa, 'NUME_EQUA', repi=neql)
!
    call jeveuo(nume_equa//'.NULG', 'E', vi=v_nugll)
    call jeveuo(nume_equa//'.PDDL', 'E', vi=v_posdd)
    call jeexin(nume_equa//'.REFP', iret)
    if (iret .ne. 0) then
        call jeveuo(nume_equa//'.REFP', 'L', vi=v_refp)
    end if
!
! -- Création du graphe de comm
    meshj = mesh//'.JOIN'
    name = nume_equa(1:14)//'.JOIN.DOMJ'
    call gnomsd(nume_equa, name, 10, 14)
    joints = name(1:19)
    zk24(jrefn-1+5) = joints
    domj = joints//'.DOMJ'
    send = joints//'.SEND'
    recv = joints//'.RECV'
    gcom = joints//'.GCOM'
    pgid = joints//'.PGID'
    lilt = nume_equa//'.LILT'
    call jeveuo(lilt, 'L', vi=v_lilt)
    call jelira(lilt, 'LONMAX', ival=nbLigrTot)
    comm_name = '&CRNUGG.COMM'
    tag_name = '&CRNUGG.TAG'
    call create_graph_comm(mesh, 'MAILLAGE_P', nb_comm, comm_name, tag_name)
!
    if (nb_comm > 0) then
        call jedupo(meshj//'.DOMJ', 'G', domj, ASTER_FALSE)
        call jedupo(meshj//'.GCOM', 'G', gcom, ASTER_FALSE)
        call jedupo(meshj//'.PGID', 'G', pgid, ASTER_FALSE)
        call jecrec(send, 'G V I', 'NU', 'DISPERSE', 'VARIABLE', nb_comm)
        call jecrec(recv, 'G V I', 'NU', 'DISPERSE', 'VARIABLE', nb_comm)
        call jeveuo(meshj//'.DOMJ', 'L', vi=v_dom)
        call jeveuo(comm_name, 'L', vi=v_comm)
        call jeveuo(tag_name, 'L', vi=v_tag)
        call jeveuo(gcom, 'L', vi=v_gco)
        call jeveuo(pgid, 'L', vi4=v_pgid)
        mpicou = to_mpi_int(v_gco(1))
    end if

!     !!!! IL PEUT ETRE INTERESSANT DE STOCKER CES INFOS
!     !!!! EN CAS DE CONSTRUCTION MULTIPLE DE NUMEDDL
!
!   RECHERCHE DES ADRESSES DU .PRNO DE .NUME
    call jeveuo(nume_equa//'.PRNO', 'E', idprn1)
    call jelira(nume_equa//'.PRNO', 'NMAXOC', nlili, k8bid)
    call jelira(nume_equa//'.LILI', 'NOMMAX', nlili2, k8bid)
    if (nlili .eq. 1) then
        ! si la longueur de la collection '.PRNO' est 1, on ne peut pas utiliser
        ! l'attribut '.LONCUM'. On crée un objet temporaire contenant le
        ! décalage pour le seul de la collection qui vaut donc 1
        call wkvect('&&CRNLGC.IDPRN', 'V V I', 1, idprn2)
        zi(idprn2) = 1
    else
        call jeveuo(jexatr(nume_equa//'.PRNO', 'LONCUM'), 'L', idprn2)
    end if
    call jeveuo(nume_equa//'.DEEQ', 'L', vi=v_deeq)

!   !!! VERIFIER QU'IL N'Y A PAS DE MACRO-ELTS
!
    call dismoi('NUM_GD_SI', nume_equa, 'NUME_EQUA', repi=gd)
    nec = nbec(gd)
    call wkvect('&&CRNUGG.NEC', 'V V I', nec, jencod)
    call wkvect('&&CRNUGG.NEC2', 'V V I', nec, jenco2)

    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgdr), 'L', iad)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgdr), 'LONMAX', ncmpmx, k8bid)
    call wkvect('&&CRNUGG.CMP', 'V V I', ncmpmx, jcpnec)
    call wkvect('&&CRNUGG.CMP2', 'V V I', ncmpmx, jcpne2)
!
!   Il faut maintenant communiquer les numeros partages
    do i_join = 1, nb_comm
        domj_i = v_comm(i_join)
        numpro = v_dom(domj_i)
        numpr2 = v_pgid(numpro+1)
        if (numpro .ne. -1) then
            nojoie = jexnum(meshj//'.SEND', domj_i)
            call jelira(nojoie, 'LONMAX', nbnoee, k8bid)
            nojoir = jexnum(meshj//'.RECV', domj_i)
            call jeveuo(nojoir, 'L', jjoinr)
            call jelira(nojoir, 'LONMAX', nbnoer, k8bid)
            nbnoee = nbnoee/2
            nbnoer = nbnoer/2
!
!       DES DEUX COTES LES NOEUDS NE SONT PAS DANS LE MEME ORDRE ?
            tag4 = to_mpi_int(v_tag(i_join))
            numpr4 = to_mpi_int(numpr2)
            lgenve1 = nbnoee*(1+nec)+1
            lgenvr1 = nbnoer*(1+nec)+1
            call wkvect('&&CRNUGG.NOEUD_NEC_E1', 'V V I', lgenvr1, jenvoi1)
            call wkvect('&&CRNUGG.NOEUD_NEC_R1', 'V V I', lgenve1, jrecep1)
!
            lgenve2 = nbnoee*(1+nec)+1
            lgenvr2 = nbnoer*(1+nec)+1
!
!       On commence par envoyer, le but final est de recevoir les numeros de ddl
!       On boucle donc sur les noeuds a recevoir
            nb_ddl_envoi = 0
            do jaux = 1, nbnoer
                poscom = (jaux-1)*(1+nec)+1
!           Numero local (pour le proc courant)
                numno1 = zi(jjoinr+2*(jaux-1))
!           Numero local (pour le proc d'en face)
                numno2 = zi(jjoinr+2*jaux-1)
                zi(jenvoi1+poscom) = numno2
                do iec = 1, nec
                    zi(jenvoi1+poscom+iec) = zzprno(1, numno1, 2+iec)
                end do
!           On ne souhaite recevoir les num de ddl globaux que du processeur
!           qui les possedent (cela explique le if( v_noex(numno1).eq.numpr2 ) )
!           Malgré tout, on envoie tout au proc d'en face (numeros de noeuds et
!           entiers codes) car il s'attend a avoir une liste dimensionnee a la
!           taille du raccord (cela explique le remplissage de zi(jenvoi1+...)
                if (v_noex(numno1) .eq. numpro) then
                    nb_ddl_envoi = nb_ddl_envoi+zzprno(1, numno1, 2)
                end if
            end do
            zi(jenvoi1) = nb_ddl_envoi
            n4e = to_mpi_int(lgenvr1)
            n4r = to_mpi_int(lgenve1)
            call asmpi_sendrecv_i(zi(jenvoi1), n4e, numpr4, tag4, &
                                  zi(jrecep1), n4r, numpr4, tag4, mpicou)

            call wkvect('&&CRNUGG.NUM_DDL_GLOB_E', 'V V I', zi(jrecep1)+1, jenvoi2)
            call wkvect('&&CRNUGG.NUM_DDL_GLOB_R', 'V V I', zi(jenvoi1)+1, jrecep2)

            nbddl = 0
            if (zi(jrecep1) > 0) then
!           On dimensionne le raccord a la taille du nombre de ddl attendus
!           par le proc d'en face
                call jeecra(jexnum(send, domj_i), 'LONMAX', zi(jrecep1))
                call jeveuo(jexnum(send, domj_i), 'E', jnujoi1)
!
                do jaux = 1, nbnoee
                    poscom = (jaux-1)*(1+nec)+1
                    numno1 = zi(jrecep1+poscom)
!               Si on ne possede pas le noeud : on envoie rien (cf.
!               if( v_noex(numno1).ne.numpr2 ) then plus bas
                    if (v_noex(numno1) .ne. rang) then
                        cycle
                    end if
!
                    nddl = zzprno(1, numno1, 1)
                    nbcmp = zzprno(1, numno1, 2)
                    if (nbcmp .eq. 0) then
                        nddlg = -1
                    else
                        nddlg = v_nugll(nddl)
                    end if
!
!               Recherche des composantes demandees
                    do iec = 1, nec
                        zi(jencod+iec-1) = zzprno(1, numno1, 2+iec)
                        zi(jenco2+iec-1) = zi(jrecep1+poscom+iec)
                    end do
                    call isdeco(zi(jencod), zi(jcpnec), ncmpmx)
                    call isdeco(zi(jenco2), zi(jcpne2), ncmpmx)
                    ico2 = 0
                    do icmp = 1, ncmpmx
!                   Comme on est sur un noeud possede par le proc courant,
!                   si le proc d'en face attend un ddl (ie composante presente),
!                   il faut necessairement que la composante soit aussi presente
!                   sur ce proc (d'ou ASSERT(zi(jcpnec+icmp-1) .ne. 0))
                        if (zi(jcpne2+icmp-1) .eq. 1) then
                            ASSERT(zi(jcpnec+icmp-1) .ne. 0)
                            ASSERT(nddlg .ne. -1)
                            zi(jenvoi2+nbddl+1) = nddlg+ico2
                            zi(jnujoi1+nbddl) = nddl+ico2
                            nbddl = nbddl+1
                            ico2 = ico2+1
                        end if
                    end do
                end do
                ASSERT(zi(jrecep1) .eq. nbddl)
                zi(jenvoi2) = nbddl
            end if
!
            n4e = to_mpi_int(nbddl+1)
            n4r = to_mpi_int(nb_ddl_envoi+1)
            call asmpi_sendrecv_i(zi(jenvoi2), n4e, numpr4, tag4, &
                                  zi(jrecep2), n4r, numpr4, tag4, mpicou)

            if (zi(jrecep2) > 0) then
                call jeecra(jexnum(recv, domj_i), 'LONMAX', nb_ddl_envoi)
                call jeveuo(jexnum(recv, domj_i), 'E', jnujoi2)
!
                curpos = 1
                do jaux = 1, nbnoer
                    numno1 = zi(jjoinr+2*(jaux-1))
!               On ne regarde que les noeuds possedes par le proc d'en face
!               puisque c'est ce qu'il a envoye
                    if (v_noex(numno1) .ne. numpro) then
                        cycle
                    end if
                    nddll = zzprno(1, numno1, 1)
                    nbcmp = zzprno(1, numno1, 2)
                    do icmp = 0, nbcmp-1
                        v_nugll(nddll+icmp) = zi(jrecep2+curpos)
                        v_posdd(nddll+icmp) = numpro
                        zi(jnujoi2+curpos-1) = nddll+icmp
                        curpos = curpos+1
                    end do
                end do
                ASSERT(curpos-1 .eq. nb_ddl_envoi)
                ASSERT(curpos-1 .eq. zi(jrecep2))
            end if
!
            call jedetr('&&CRNUGG.NUM_DDL_GLOB_E')
            call jedetr('&&CRNUGG.NUM_DDL_GLOB_R')
!
            call jedetr('&&CRNUGG.NOEUD_NEC_E1')
            call jedetr('&&CRNUGG.NOEUD_NEC_R1')
        end if
    end do
    call jedetr(comm_name)
    call jedetr(tag_name)

    do iLigrT = 2, nbLigrTot
        if (v_lilt(iLigrT) .eq. -1) then
            cycle
        else
            ili = v_lilt(iLigrT)
        end if
        call jenuno(jexnum(nume_equa//'.LILI', ili), nomlig)
        lligrel_cp = .false.
        call jeexin(nomlig//'._TCO', iret)
        if (iret .ne. 0) then
            call jeveuo(nomlig//'._TCO', "L", vk24=tco)
            lligrel_cp = (tco(1) .eq. 'LIGREL_CP')
        end if
        if (lligrel_cp) then
            ili2 = v_refp(ili)
        else
            ili2 = ili
        end if
        call jeexin(jexnum(nume_equa//'.PRNO', ili), exipr)
        call create_graph_comm(nomlig, 'LIGREL', nb_comm, comm_name, tag_name)
        if (nb_comm > 0) then
            call jeveuo(comm_name, 'L', vi=v_comm)
            call jeveuo(tag_name, 'L', vi=v_tag)
            call dismoi('JOINTS', nomlig, 'LIGREL', repk=joints, arret='F')
            domj = joints//'.DOMJ'
            send = joints//'.SEND'
            recv = joints//'.RECV'
            pgid = joints//'.PGID'
            gcom = joints//'.GCOM'
            call jeveuo(pgid, 'L', vi4=v_pgid)
            call jeveuo(domj, 'L', vi=v_dom)
            call jeveuo(gcom, 'L', vi=v_gco)
            mpicou = to_mpi_int(v_gco(1))
            do i_join = 1, nb_comm
                domj_i = v_comm(i_join)
                numpro = v_dom(domj_i)
                numpr2 = v_pgid(numpro+1)
                numpr4 = to_mpi_int(numpr2)
                tag4 = to_mpi_int(v_tag(i_join))
                nojoie = jexnum(send, domj_i)

                call jeexin(nojoie, iret1)
                if (iret1 .ne. 0) then
                    if (lligrel_cp) then
                        call jeveuo(nojoie, 'L', jjoine)
                        call jelira(nojoie, 'LONMAX', nbnoee, k8bid)
                        call wkvect('&&CRNUGL.NUM_DDL_GLOB_E', 'V V I', 2*nbnoee, jnujoi1)
                        do jaux = 1, nbnoee
                            numnoe = zi(jjoine+jaux-1)
                            nddll = zzprno(1, numnoe, 1)
                            zi(jnujoi1+2*(jaux-1)) = v_nugll(nddll)
                            zi(jnujoi1+2*(jaux-1)+1) = zzprno(1, numnoe, 2)
                        end do
                        n4e = to_mpi_int(2*nbnoee)
                    else
                        call jeveuo(nojoie, 'L', jjoine)
                        call jelira(nojoie, 'LONMAX', nbnoee, k8bid)
                        call wkvect('&&CRNUGL.NUM_DDL_GLOB_E', 'V V I', nbnoee, jnujoi1)
                        do jaux = 1, nbnoee
                            numnoe = -zi(jjoine+jaux-1)
                            ASSERT(zzprno(ili2, numnoe, 2) .eq. 1)
                            nddll = zzprno(ili2, numnoe, 1)
                            zi(jnujoi1+jaux-1) = v_nugll(nddll)
                        end do
                        n4e = to_mpi_int(nbnoee)
                    end if
                end if

                nojoir = jexnum(recv, domj_i)
                call jeexin(nojoir, iret2)
                if (iret2 .ne. 0) then
                    call jeveuo(nojoir, 'L', jjoinr)
                    call jelira(nojoir, 'LONMAX', nbnoer, k8bid)
                    if (lligrel_cp) then
                        call wkvect('&&CRNUGL.NUM_DDL_GLOB_R', 'V V I', 2*nbnoer, jnujoi2)
                        n4r = to_mpi_int(2*nbnoer)
                    else
                        call wkvect('&&CRNUGL.NUM_DDL_GLOB_R', 'V V I', nbnoer, jnujoi2)
                        n4r = to_mpi_int(nbnoer)
                    end if
                end if
                if (rang .lt. numpro) then
                    if (iret1 .ne. 0) then
                        call asmpi_send_i(zi(jnujoi1), n4e, numpr4, tag4, mpicou)
                    end if
                    if (iret2 .ne. 0) then
                        call asmpi_recv_i(zi(jnujoi2), n4r, numpr4, tag4, mpicou)
                    end if
                else if (rang .gt. numpro) then
                    if (iret2 .ne. 0) then
                        call asmpi_recv_i(zi(jnujoi2), n4r, numpr4, tag4, mpicou)
                    end if
                    if (iret1 .ne. 0) then
                        call asmpi_send_i(zi(jnujoi1), n4e, numpr4, tag4, mpicou)
                    end if
                end if

                if (iret2 .ne. 0 .and. exipr .ne. 0) then
                    if (lligrel_cp) then
                        do jaux = 1, nbnoer
                            numnoe = -zi(jjoinr+jaux-1)
                            nbddl = zi(jnujoi2+2*(jaux-1)+1)
                            nddll = zzprno(ili2, numnoe, 1)
                            do iddl = 0, nbddl-1
                                v_nugll(nddll+iddl) = zi(jnujoi2+2*(jaux-1))+iddl
                                v_posdd(nddll+iddl) = numpro
                            end do
                        end do
                    else
                        do jaux = 1, nbnoer
                            numnoe = -zi(jjoinr+jaux-1)
                            ASSERT(zzprno(ili2, numnoe, 2) .eq. 1)
                            nddll = zzprno(ili2, numnoe, 1)
                            v_nugll(nddll) = zi(jnujoi2+jaux-1)
                            v_posdd(nddll) = numpro
                        end do
                    end if
                end if
                call jedetr('&&CRNUGL.NUM_DDL_GLOB_E')
                call jedetr('&&CRNUGL.NUM_DDL_GLOB_R')
            end do
        end if
        call jedetr(comm_name)
        call jedetr(tag_name)
    end do
!
    call jedetc('V', '&&CRNUGG', 1)
    call jedetr('&&CRNLGC.IDPRN')
!
! --- Vérification de la numérotation
!
!   NOMBRE DE DDL LOCAUX
    call jeveuo(nume_equa//'.NEQU', 'L', jnequ)
    nbddll = zi(jnequ)
    do iaux = 1, nbddll
        nuno = v_deeq((iaux-1)*2+1)
        if (nuno .le. 0) cycle
        if (nuno .ne. 0) then
            if (v_posdd(iaux) .ne. v_noex(nuno)) then
                call utmess('F', 'ASSEMBLA_9')
            end if
        else
            ASSERT(v_posdd(iaux) .ne. -1)
        end if
        ASSERT(v_nugll(iaux) .ne. -1)
    end do
!
! --- Affichage inconnue systeme
    if (niv .ge. 1) then
        call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_node)
        nno = 0
        do ino = 1, nb_node
            if (v_noex(ino) == rang) then
                if (zi(idprn1-1+(ino-1)*(2+nec)+2) .gt. 0) nno = nno+1
            end if
        end do
        call asmpi_comm_vect('MPI_SUM', 'I', sci=nno)
!
        call jeexin(nume_equa//'.MDLA', iret)
        if (iret .ne. 0) then
            call jelira(nume_equa//'.MDLA', 'LONMAX', nlag, k8bid)
            nlag = nlag/3
        else
            nlag = 0
        end if
        call asmpi_comm_vect('MPI_SUM', 'I', sci=nlag)
!
        vali(1) = zi(jnequ+1)
        vali(2) = zi(jnequ+1)-nlag
        vali(3) = nno
        vali(4) = nlag
        vali(5) = nlag/2
        call utmess('I', 'FACTOR_1', ni=5, vali=vali)
    end if
!
! --- Pour debuggage en hpc
    if (ASTER_FALSE) then
        nonulg = mesh//'.NUNOLG'
        call jeveuo(nonulg, 'L', jmlogl)
        do iaux = 0, nbddll-1
            nuno = v_deeq(iaux*2+1)
            if (nuno .ne. 0) nuno = zi(jmlogl+nuno-1)+1
! numero ddl local, numéro noeud local, numéro noeud global, num composante du noeud,
!            num ddl global, num ddl seq, num proc proprio
            write (130+rang, *) iaux, v_deeq(iaux*2+1), v_deeg(iaux*2+1), v_deeq(iaux*2+2), &
                v_nugll(iaux+1), v_nuls(iaux+1), v_posdd(1+iaux)

            write (190+rang, *) iaux, v_deeg(iaux*2+1), v_nuls(iaux+1), v_nugll(iaux+1)
        end do
        flush (130+rang)
        flush (190+rang)
    end if
!
999 continue
!
    call jedema()
#else
    character(len=19) :: kbid
    kbid = nume_equa
#endif
!
end subroutine
