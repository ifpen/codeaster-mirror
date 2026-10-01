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
function nbddlMaxMa(nume_ddlz, matr_assez, nbmat, v_name_mat) result(maxDDLMa)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jaexin.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/nbno.h"
#include "asterfort/parti0.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: nume_ddlz, matr_assez
    integer(kind=8), intent(in)  :: nbmat
    character(len=*), intent(in) :: v_name_mat(nbmat)
    integer(kind=8)              :: maxDDLMa
!     ------------------------------------------------------------------
!     RETOURNE LE NOMBRE DE DDL MAXIMUM POUR UN ELEMENT
!     ------------------------------------------------------------------
! IN  K14 nume_ddlz     : NOM DE LA SD_NUME_DDL
! IN  K19 matr_assez    : NOM DE LA MATR_ASSE
! IN  I   nbmat         : nombre de matr_elem dans la MATR_ASSE
!     ------------------------------------------------------------------
    character(len=3)  :: answer
    character(len=8)  :: mesh, model, model_elem, nogdco, nogdsi, partition
    character(len=14) :: nume_ddl
    character(len=19) :: nomligrel, matr_elem, resu_elem, matr_asse
    integer(kind=8) :: iconx1, iconx2, iel, iret, nnoe
    integer(kind=8) :: igrel, numa, ino, n1, n12, nddl1, rang, jrefa, jdesc
    integer(kind=8) :: nddlt, nel, nec, mode, nugd, imat, nbssa, iamail
    integer(kind=8) :: nb_resu_elem, iresu, ilima, ilinu, nbproc, ilinu_ref
    aster_logical :: l_dgrel, l_distme, l_matd, lligrel_cp
    mpi_int :: mrank, msize
!
    integer(kind=8), pointer :: v_adne(:) => null()
    integer(kind=8), pointer :: v_adli(:) => null()
    integer(kind=8), pointer :: v_numsd(:) => null()
    integer(kind=8), pointer :: v_prn1(:) => null()
    integer(kind=8), pointer :: v_prn2(:) => null()
    integer(kind=8), pointer :: sssa(:) => null()
    integer(kind=8), pointer :: v_refp(:) => null()
    integer(kind=8), pointer :: v_crco(:) => null()
    character(len=24), pointer :: v_relr(:) => null()
    character(len=24), pointer :: v_nomlig(:) => null()
    character(len=24), pointer :: v_prtk(:) => null()
    character(len=24), pointer :: tco(:) => null()
!
!-----------------------------------------------------------------------
!     FONCTIONS LOCALES D'ACCES AUX DIFFERENTS CHAMPS DES
!     S.D. MANIPULEES DANS LE SOUS PROGRAMME
!-----------------------------------------------------------------------
!---- FONCTION D ACCES AU CHAMP CONNEX DE LA S.D. MAILLA DE TYPE
!     MAILLAGE
!     ZZCONX(IMAIL,J) = NUMERO DANS LA NUMEROTATION DU MAILLAGE
!         DU NOEUD J DE LA MAILLE IMAIL
#define zzconx(imail,j) zi(iconx1-1+zi(iconx2+imail-1)+j-1)
!
!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS LIEL DES S.D. LIGREL
!     REPERTORIEES DANS LE REPERTOIRE TEMPORAIRE .MATAS.LILI
!     ZZLIEL(ILI,IGREL,J) =
!      SI LA JIEME MAILLE DU LIEL IGREL DU LIGREL ILI EST:
!          -UNE MAILLE DU MAILLAGE : SON NUMERO DANS LE MAILLAGE
!          -UNE MAILLE TARDIVE : -POINTEUR DANS LE CHAMP .NEMA
!
#define zzliel(ili,igrel,j) zi(v_adli(1+3*(ili-1)+1)-1+ \
    zi(v_adli(1+3*(ili-1)+2)+igrel-1)+j-1)
!
!---- NBRE DE GROUPES D'ELEMENTS (DE LIEL) DU LIGREL ILI
!
#define zzngel(ili) v_adli(1+3*(ili-1))
!
!---- NBRE D ELEMENTS DU LIEL IGREL DU LIGREL ILI DU REPERTOIRE TEMP.
!     .MATAS.LILI(DIM DU VECTEUR D'ENTIERS .LILI(ILI).LIEL(IGREL) )
!
#define zznelg(ili,igrel) zi(v_adli(1+3*(ili-1)+2)+igrel)- \
    zi(v_adli(1+3*(ili-1)+2)+igrel-1)-1
!
!---- NBRE D ELEMENTS SUPPLEMENTAIRE (.NEMA) DU LIGREL ILI DU
!     REPERTOIRE TEMPORAIRE .MATAS.LILI
!
!
!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS NEMA DES S.D. LIGREL
!     REPERTORIEES DANS LE REPERTOIRE TEMPO. .MATAS.LILI
!
#define zznema(ili,iel,j) zi(v_adne(1+3*(ili-1)+1)-1+ \
    zi(v_adne(1+3*(ili-1)+2)+iel-1)+j-1)
!
!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS PRNO DES S.D. LIGREL
!     REPERTORIEES DANS nume_ddl.LILI DE LA S.D. NUME_DDL ET A LEURS ADRESSES
!     ZZPRNO(ILI,NUNOEL,1) = NUMERO DE L'EQUATION ASSOCIEES AU 1ER DDL
!                            DU NOEUD NUNOEL DANS LA NUMEROTATION LOCALE
!                            AU LIGREL ILI DE .LILI
!     ZZPRNO(ILI,NUNOEL,2) = NOMBRE DE DDL PORTES PAR LE NOEUD NUNOEL
!     ZZPRNO(ILI,NUNOEL,2+1) = 1ER CODE
!     ZZPRNO(ILI,NUNOEL,2+NEC) = NEC IEME CODE
!
#define zzprno(ili,nunoel,l) v_prn1(v_prn2(ili)+(nunoel-1)*(nec+2)+l-1)
!----------------------------------------------------------------------
!
    call jemarq()
!
    matr_asse = matr_assez
    nume_ddl = nume_ddlz
!
! --- Verification nume_ddl
!
    call jeexin(matr_asse//'.REFA', iret)
    if (iret .gt. 0) then
        call jeveuo(matr_asse//'.REFA', 'L', jrefa)
        ASSERT(zk24(jrefa-1+2) (1:14) .eq. nume_ddl)
    end if
!
    call dismoi('NOM_MODELE', nume_ddl, 'NUME_DDL', repk=model)
    call dismoi('NOM_MAILLA', nume_ddl, 'NUME_DDL', repk=mesh)
    call dismoi('NOM_GD', nume_ddl, 'NUME_DDL', repk=nogdco)
    call dismoi('NOM_GD_SI', nogdco, 'GRANDEUR', repk=nogdsi)
    call dismoi('NUM_GD_SI', nogdsi, 'GRANDEUR', repi=nugd)
!
    nec = nbec(nugd)
!
    call jeexin(model//'.MODELE    .SSSA', iret)
    if (iret .gt. 0) then
        call jeveuo(model//'.MODELE    .SSSA', 'L', vi=sssa)
    end if
!
    call jeexin(mesh//'.CONNEX', iret)
    if (iret .gt. 0) then
        call jeveuo(mesh//'.CONNEX', 'L', iconx1)
        call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', iconx2)
    else
        iconx1 = 0
        iconx2 = 0
    end if
!
! --- calcul de l_matd et jnueq
!
    call dismoi('MATR_DISTRIBUEE', nume_ddl, 'NUME_DDL', repk=answer)
    l_matd = (answer .eq. 'OUI')
!
    call jelira(nume_ddl//'.NUME.REFN', 'LONMAX', n1)
    ASSERT(n1 .eq. 5)
!
! --- Acces au NUME_EQUA
!
    call jeexin(nume_ddl//'.NUME.REFP', iret)
    if (iret .ne. 0) then
        call jeveuo(nume_ddl//'.NUME.REFP', 'L', vi=v_refp)
    end if
    if (l_matd) then
        call jeveuo(nume_ddl//'.NUML.PRNO', 'L', vi=v_prn1)
        call jeveuo(jexatr(nume_ddl//'.NUML.PRNO', 'LONCUM'), 'L', vi=v_prn2)
    else
        call jeveuo(nume_ddl//'.NUME.PRNO', 'L', vi=v_prn1)
        call jeveuo(jexatr(nume_ddl//'.NUME.PRNO', 'LONCUM'), 'L', vi=v_prn2)
    end if
!
! --- objet matr_asse.lili
!
    call jeveuo(matr_asse//'.ADNE', 'L', vi=v_adne)
    call jeveuo(matr_asse//'.ADLI', 'L', vi=v_adli)
!
! --- calcul de l_distme, partition, l_dgrel, v_numsd
!
    call parti0(nbmat, v_name_mat, partition)
!
    call exisd("PARTITION", partition, iret)
    if (iret .ne. 0) then
        l_distme = ASTER_TRUE
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
        call jeveuo(partition//'.PRTK', 'L', vk24=v_prtk)
        l_dgrel = (v_prtk(1) .eq. 'SOUS_DOMAINE') .or. (v_prtk(1) .eq. 'GROUP_ELEM')
        if (.not. l_dgrel) then
            call jeveuo(partition//'.NUPR', 'L', vi=v_numsd)
        end if
    else
        rang = 0
        nbproc = 1
        l_distme = ASTER_FALSE
        l_dgrel = ASTER_FALSE
    end if
!
    if (l_matd) then
        ASSERT(l_distme)
    end if
!
! --- Calcul du nombre maximum de ddl par maille
!
!
! --- boucle sur les matr_elem
!
    maxDDLMa = 0
!
    do imat = 1, nbmat
        matr_elem = v_name_mat(imat)
        call dismoi('NOM_MODELE', matr_elem, 'MATR_ELEM', repk=model_elem)
!
        if (model_elem .ne. model) then
            call utmess('F', 'ASSEMBLA_5')
        end if
!
! --- traitement des elements finis classiques
!
        call jeexin(matr_elem//'.RELR', iret)
        if (iret .ne. 0) then
!
            call jelira(matr_elem//'.RELR', 'LONUTI', nb_resu_elem)
            if (nb_resu_elem .gt. 0) then
                call jeveuo(matr_elem//'.RELR', 'L', vk24=v_relr)
!
                do iresu = 1, nb_resu_elem
                    resu_elem = v_relr(iresu) (1:19)
!
                    call jeexin(resu_elem//'.DESC', iret)
                    if (iret .ne. 0) then
!
                        call dismoi('MPI_COMPLET', resu_elem, 'RESUELEM', repk=answer)
                        if (answer .eq. 'NON') then
                            ASSERT(l_distme)
                        end if
!
                        call jeveuo(resu_elem//'.NOLI', 'L', vk24=v_nomlig)
                        nomligrel = v_nomlig(1) (1:19)

                        call jenonu(jexnom(matr_asse//'.LILI', nomligrel), ilima)
                        call jenonu(jexnom(nume_ddl//'.NUME.LILI', nomligrel), ilinu)
                        ilinu_ref = ilinu
                        lligrel_cp = .false.
                        call jeexin(nomligrel//'._TCO', iret)
                        if (iret .ne. 0) then
                            call jeveuo(nomligrel//'._TCO', "L", vk24=tco)
                            lligrel_cp = (tco(1) .eq. 'LIGREL_CP')
                        end if
                        if (lligrel_cp) then
                            call asmpi_info(size=msize)
                            if (msize > 1) then
                                call jeveuo(jexnum(nume_ddl//'.NUME.CRCO', ilinu), 'L', vi=v_crco)
                            end if
                            ilinu_ref = v_refp(ilinu)
                        end if
                        call jeexin(nomligrel//'.NTCM', iret)

                        call jeveuo(resu_elem//'.DESC', 'L', jdesc)
!
! --- boucle sur les grels du ligrel
!
                        do igrel = 1, zzngel(ilima)
                            if (l_dgrel .and. mod(igrel, nbproc) .ne. rang) goto 100
!
! --- il se peut que le grel igrel soit vide
                            call jaexin(jexnum(resu_elem//'.RESL', igrel), iret)
                            if (iret .eq. 0) goto 100
!
                            mode = zi(jdesc+igrel+1)
                            if (mode .gt. 0) then
                                nnoe = nbno(mode)
!
! --- nombre d'elements du grel igrel du ligrel nomligrel/ilima
!
                                nel = zznelg(ilima, igrel)
!
! --- boucle sur les elements du grel
!
                                do iel = 1, nel
                                    nddlt = 0
                                    numa = zzliel(ilima, igrel, iel)
!
! --- SI LES CALCULS ONT ETE DISTRIBUES
!
                                    if (l_distme .and. .not. l_dgrel) then
!       SI ON EST DANS UN CALCUL DISTRIBUE, ON SE POSE
!       LA QUESTION DE L'APPARTENANCE DE LA MAILLE NUMA AUX
!       DONNEES ATTRIBUEES AU PROC SI MAILLE PHYSIQUE: CHAQUE PROC
!       NE TRAITE QUE CELLES ASSOCIEES AUX SD QUI LUI SONT ATTRIBUES
!       SI MAILLE TARDIVE: ELLES SONT TRAITEES PAR LE PROC 0
                                        if (numa .gt. 0) then
                                            if (v_numsd(numa) .ne. rang) goto 200
                                        else
                                            if (rang .ne. 0) goto 200
                                        end if
                                    end if
!
                                    if (numa .gt. 0) then
! --- MAILLE DU MAILLAGE
                                        if (lligrel_cp) then
                                            ASSERT(.false.)
                                        end if
                                        do ino = 1, nnoe
                                            n1 = zzconx(numa, ino)
                                            nddl1 = zzprno(1, n1, 2)
                                            nddlt = nddlt+nddl1
                                        end do
                                    else
! --- MAILLE TARDIVE
                                        numa = -numa
                                        do ino = 1, nnoe
! --- N1 : INDICE DU NOEUDS DS LE .NEMA DU LIGREL DE CHARGE GLOBAL OU LOCAL
                                            n1 = zznema(ilima, numa, ino)
                                            if (lligrel_cp .and. n1 .lt. 0) then
                                                n12 = -v_crco(-n1)
                                            else
                                                n12 = n1
                                            end if
                                            if (n1 .lt. 0) then
! --- NOEUD TARDIF
                                                n1 = -n12
                                                nddl1 = zzprno(ilinu_ref, n1, 2)
                                            else
! --- NOEUD PHYSIQUE
                                                nddl1 = zzprno(1, n12, 2)
                                            end if
                                            nddlt = nddlt+nddl1
                                        end do
                                    end if
!
200                                 continue
                                    maxDDLMa = max(maxDDLMa, nddlt)
                                end do
                            end if
100                         continue
                        end do
                    end if
                end do
            end if
        end if
!
! --- traitement des sous-structures statiques
!
        call dismoi('NB_SS_ACTI', matr_elem, 'MATR_ELEM', repi=nbssa)
!
        if (nbssa > 0) then
!
! --- boucle sur les macro-elements
!
            do iel = 1, nbssa
                nddlt = 0
                if (sssa(iel) .ne. 0) then
                    call jeveuo(jexnum(mesh//'.SUPMAIL', iel), 'L', iamail)
                    call jelira(jexnum(mesh//'.SUPMAIL', iel), 'LONMAX', nnoe)

                    do ino = 1, nnoe
                        n1 = zi(iamail-1+ino)
                        ASSERT(n1 .ne. 0)
                        nddl1 = zzprno(1, n1, 2)
                        nddlt = nddlt+nddl1
                    end do
                end if
                maxDDLMa = max(maxDDLMa, nddlt)
            end do
        end if
    end do
!
    call jedema()
!
end function
