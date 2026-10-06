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
! aslint: disable=W1504
!
subroutine ircame(ifi, nochmd, chanom, typech, ligrel, &
                  nbcmp, nomcmp, etiqcp, partie, numpt, &
                  instan, numord, adsk, adsd, adsc, &
                  adsv, adsl, nbenec, lienec, sdcarm, &
                  carael, field_type, nbCmpDyna, lfichUniq, codret)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/asmpi_allgather_i.h"
#include "asterc/asmpi_allgatherv_char16.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/ircam1.h"
#include "asterfort/ircmpr.h"
#include "asterfort/irelst.h"
#include "asterfort/iremed_filtre.h"
#include "asterfort/irmail.h"
#include "asterfort/irmpga.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/lrmtyp.h"
#include "asterfort/mdexch.h"
#include "asterfort/mdexma.h"
#include "asterfort/mdnoma.h"
#include "asterfort/ulisog.h"
#include "asterfort/utlicm.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "MeshTypes_type.h"
!
    character(len=8) :: typech, sdcarm, carael
    character(len=19) :: chanom, ligrel
    character(len=64) :: nochmd
    character(len=*) :: nomcmp(*), partie, etiqcp
    integer(kind=8) :: nbcmp, numpt, numord, ifi
    integer(kind=8) :: adsk, adsd, adsc, adsv, adsl
    integer(kind=8) :: nbenec
    integer(kind=8) :: lienec(*)
    integer(kind=8) :: typent, tygeom
    real(kind=8) :: instan
    character(len=16), intent(in) :: field_type
    integer(kind=8), intent(inout) :: nbCmpDyna
    aster_logical :: lfichUniq
    integer(kind=8) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE D'UN CHAMP - FORMAT MED
!
! --------------------------------------------------------------------------------------------------
!
!     ENTREES :
!       IFI    : UNITE LOGIQUE D'IMPRESSION DU CHAMP
!       NOCHMD : NOM MED DU CHAMP A ECRIRE
!       PARTIE: IMPRESSION DE LA PARTIE IMAGINAIRE OU REELLE POUR
!               UN CHAMP COMPLEXE
!       CHANOM : NOM ASTER DU CHAM A ECRIRE
!       TYPECH : TYPE DU CHAMP ('NOEU', 'ELNO', 'ELGA')
!       LIGREL : LIGREL ASSOCIE AU CHAMP
!       NBCMP  : NOMBRE DE COMPOSANTES A ECRIRE. S'IL EST NUL, ON
!                PREND TOUT
!       NOMCMP : NOMS DES COMPOSANTES A ECRIRE
!       NUMPT  : NUMERO DE PAS DE TEMPS
!       INSTAN : VALEUR DE L'INSTANT A ARCHIVER
!       NUMORD : NUMERO D'ORDRE DU CHAMP
!       ADSK, D, ... : ADRESSES DES TABLEAUX DES CHAMPS SIMPLIFIES
!       NBVATO : NOMBRE DE VALEURS TOTALES
!       NBENEC : NOMBRE D'ENTITES A ECRIRE (O, SI TOUTES)
!       LIENEC : LISTE DES ENTITES A ECRIRE SI EXTRAIT
!       SDCARM : CARA_ELEM (UTILE POUR LES SOUS-POINTS)
!       CARAEL : NOM CARA ELEM
!       LFICHUNIQ : LOGICAL FICHIER UNIQUE
! In  field_type       : type of field (symbolic name in result datastructure)
!     SORTIES:
!       CODRET : CODE DE RETOUR (0 : PAS DE PB, NON NUL SI PB)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6), parameter :: nompro = 'IRCAME'
    integer(kind=8), parameter :: ednoeu = 3, edmail = 0, ednoma = 4, typnoe = 0
    character(len=1)   :: saux01
    character(len=8)   :: saux08, nomaas, nomtyp(MT_NTYMAX), nom_sd_fu
    character(len=16)  :: formar
    character(len=24)  :: ntlcmp, ntncmp, ntucmp, ntproa, nmcmfi, ncaimi, ncaimk
    character(len=24)  :: indcmp, nosdpr
    character(len=32)  :: nomgd
    character(len=64)  :: nomamd, dschmd
    character(len=200) :: nofimd
    character(len=255) :: kfic
    med_idt :: ifimed
    integer(kind=8) :: nbtyp, ifm, nivinf, lnomam, jindir
    integer(kind=8) :: ncmpve, nvalec, nbprof, nbvato, ncmprf
    integer(kind=8) :: nbimpr, jnocm1, jnocm2, nbcmp2, icmp1, icmp2
    integer(kind=8) :: adcaii, adcaik, ncmpvl, jnocm3
    integer(kind=8) :: iaux, jaux, nrimpr, jtest, nb_cmp_tot, jnocmp, numcmp, cmpt
    integer(kind=8) :: existc, nbcmfi, nbval, nbcmpmax, vnbcmp(1), rang, nbproc, iproc
    aster_logical :: lgaux, existm
    integer(kind=8) :: iCmp
    character(len=8), pointer :: cmpUserName(:) => null()
    character(len=8), pointer :: cmpCataName(:) => null()
    character(len=8), pointer :: v_ma(:) => null()
    character(len=16), pointer :: v_nomcmp(:) => null()
    character(len=16), pointer :: v_nomcm2(:) => null()
    real(kind=8) :: start_time, end_time
    mpi_int :: mrank, mnbproc, world, one4, taille
    mpi_int, pointer :: v_count(:) => null()
    mpi_int, pointer :: v_displ(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: typgeo(MT_NTYMAX)
    integer(kind=8) :: nnotyp(MT_NTYMAX)
    integer(kind=8) :: modnum(MT_NTYMAX)
    integer(kind=8) :: numnoa(MT_NTYMAX, MT_NNOMAX), nuanom(MT_NTYMAX, MT_NNOMAX)
    integer(kind=8) :: renumd(MT_NTYMAX)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! 1. PREALABLES
!   1.1. ==> RECUPERATION DU NIVEAU D'IMPRESSION
    call infniv(ifm, nivinf)
!   1.2. ==> NOMS DES TABLEAUX DE TRAVAIL
!             12   345678   9012345678901234
    ntlcmp = '&&'//nompro//'.LISTE_N0MCMP   '
    ntncmp = '&&'//nompro//'.NOMCMP         '
    ntucmp = '&&'//nompro//'.UNITECMP       '
    ntproa = '&&'//nompro//'.PROFIL_ASTER   '
    nmcmfi = '&&'//nompro//'.NOMCMP_FICHIER '
    ncaimi = '&&'//nompro//'.CARAC_NOMBRES__'
    ncaimk = '&&'//nompro//'.CARAC_CHAINES__'
!   1.3. ==> NOM DU FICHIER MED
    call ulisog(ifi, kfic, saux01)
    if (kfic(1:1) .eq. ' ') then
        call codent(ifi, 'G', saux08)
        nofimd = 'fort.'//saux08
    else
        nofimd = kfic(1:200)
    end if
!
    if (nivinf .gt. 1) then
        call cpu_time(start_time)
        write (ifm, *) '<', nompro, '> DEBUT ECRITURE DU CHAMP MED ', typech(1:4), ' :'
        write (ifm, *) '<', nompro, '> NOM DU FICHIER MED : ', nofimd
    end if
!   1.4. ==> LES NOMBRES CARACTERISTIQUES
    nbvato = zi(adsd)
    ncmprf = zi(adsd+1)
!
! 2. LE MAILLAGE
!   2.1. ==> NOM DU MAILLAGE ASTER
    nomaas = zk8(adsk-1+1)
!   Generate name of mesh for MED
    call mdnoma(nomamd, lnomam, nomaas, codret)
!   Creation des vecteurs jeveux pour les filtres med
    nom_sd_fu = ' '
!   2.3. ==> CE MAILLAGE EST-IL DEJA PRESENT DANS LE FICHIER ?
    iaux = 0
    ifimed = 0
    call mdexma(nofimd, ifimed, nomamd, iaux, existm, jaux, codret)
!   2.4. ==> SI LE MAILLAGE EST ABSENT, ON L'ECRIT
    if (.not. existm) then
        saux08 = 'MED     '
        lgaux = .false.
        formar = ' '
        if (lfichUniq) then
            if (.not. isParallelMesh(nomaas)) then
                call utmess('F', 'MED3_5')
            end if
            nom_sd_fu = '&&IRMHD2'
            call jedetc('V', nom_sd_fu, 1)
            call iremed_filtre(nomaas, nom_sd_fu, 'V', ASTER_TRUE)
        end if
        call irmail(saux08, ifi, iaux, nomaas, lgaux, ligrel, nivinf, formar, lfichUniq, &
                    nom_sd_fu)
    else
        if (lfichUniq) then
            if (.not. isParallelMesh(nomaas)) then
                call utmess('F', 'MED3_5')
            end if
            nom_sd_fu = '&&IRMHD2'
            call jeveuo(nom_sd_fu//'.NOMA', "L", vk8=v_ma)
            ASSERT(v_ma(1) == nomaas)
            ! call jeexin(nom_sd_fu//'.NOMA', codret)
            if (codret .eq. 0) then
                !call iremed_filtre(nomaas, nom_sd_fu, 'V', ASTER_TRUE)
            end if
        end if
    end if
!
! 3. PREPARATION DU CHAMP A ECRIRE
!   3.1. ==> NUMEROS, NOMS ET UNITES DES COMPOSANTES A ECRIRE
    if (nbCmp > 0) then
        AS_ALLOCATE(vk8=cmpUserName, size=nbcmp)
!
        do iCmp = 1, nbCmp
            cmpUserName(iCmp) = nomcmp(iCmp)
        end do
    end if
!
    if (ncmprf > 0) then
        AS_ALLOCATE(vk8=cmpCataName, size=ncmprf)
        do iCmp = 1, ncmprf
            cmpCataName(iCmp) = zk8(adsc-1+iCmp)
        end do
    else
        ! avoid leaks
        if (nbCmp > 0) AS_DEALLOCATE(vk8=cmpUserName)
        goto 998
    end if
    call utlicm(zk8(adsk+1), &
                nbcmp, cmpUserName, &
                ncmprf, cmpCataName, &
                ncmpve, ntlcmp, &
                ntncmp, ntucmp)
    AS_DEALLOCATE(vk8=cmpUserName)
    AS_DEALLOCATE(vk8=cmpCataName)
    if (ncmpve .gt. 200) then
        call utmess('A', 'MED_99', sk=nochmd)
        goto 999
    end if
    indcmp = '&&IRCAME.CMPINDIR'

    call dismoi('NOM_GD', chanom, 'CHAMP', repk=nomgd)
    dschmd = trim(nomgd)//'-'//trim(field_type)

    if (lfichUniq) then
        call jeveuo(ntncmp, 'L', jnocm1)
        call jeveuo(ntucmp, 'L', jnocm2)
        call jecreo('&&IRCAME.CMPLOC', 'V N K16')
        call jeecra('&&IRCAME.CMPLOC', 'NOMMAX', ncmpve)
        cmpt = 0
        do icmp1 = 1, ncmpve
            call jecroc(jexnom('&&IRCAME.CMPLOC', zk16(jnocm1+icmp1-1)))
        end do
        one4 = to_mpi_int(1)
        if (nomgd .eq. 'VARI_R') then
            vnbcmp(1) = ncmpve
            call asmpi_comm_vect('MPI_MAX', 'I', 1, vi=vnbcmp)
            nbcmpmax = vnbcmp(1)
        else
            call dismoi('NB_CMP_MAX', nomgd, 'GRANDEUR', repi=nbcmpmax)
        end if
        call asmpi_comm('GET', world)
        call asmpi_info(rank=mrank, size=mnbproc)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(mnbproc)
        call wkvect('&&IRCAME.TEST', 'V V I', nbproc, jtest)
        call wkvect('&&IRCAME.COUNT', 'V V S', nbproc, vi4=v_count)
        call wkvect('&&IRCAME.DISPL', 'V V S', nbproc+1, vi4=v_displ)
        call wkvect('&&IRCAME.NOMCMP2', 'V V K16', ncmpve, jnocmp)
        call wkvect('&&IRCAME.NOMCMP3', 'V V K16', ncmpve, jnocm3)
        call asmpi_allgather_i([ncmpve], one4, zi(jtest), one4, world)
        nb_cmp_tot = 0
        do iproc = 0, nbproc-1
            nb_cmp_tot = nb_cmp_tot+zi(jtest+iproc)
            v_count(iproc+1) = to_mpi_int(zi(jtest+iproc))
            v_displ(iproc+2) = to_mpi_int(nb_cmp_tot)
        end do
        do icmp1 = 1, ncmpve
            zk16(jnocmp+icmp1-1) = zk16(jnocm1+icmp1-1)
            zk16(jnocm3+icmp1-1) = zk16(jnocm2+icmp1-1)
        end do
        call wkvect('&&IRCAME.NOMFAG', 'V V K16', nb_cmp_tot, vk16=v_nomcmp)
        call wkvect('&&IRCAME.NOMFA2', 'V V K16', nb_cmp_tot, vk16=v_nomcm2)
        call jedetr(ntncmp)
        call wkvect(ntncmp, 'V V K16', nb_cmp_tot, jnocm1)
        call jedetr(ntucmp)
        call wkvect(ntucmp, 'V V K16', nb_cmp_tot, jnocm2)
        call wkvect(indcmp, 'V V I', nb_cmp_tot, jindir)
        taille = to_mpi_int(ncmpve)
        call asmpi_allgatherv_char16(zk16(jnocmp), taille, v_nomcmp, v_count, v_displ, world)
        call asmpi_allgatherv_char16(zk16(jnocm3), taille, v_nomcm2, v_count, v_displ, world)
        call jecreo('&&IRCAME.PTRNOM', 'V N K16')
        call jeecra('&&IRCAME.PTRNOM', 'NOMMAX', nb_cmp_tot)
        cmpt = 0
        do icmp1 = 1, nb_cmp_tot
            call jenonu(jexnom('&&IRCAME.PTRNOM', v_nomcmp(icmp1)), numcmp)
            if (numcmp .eq. 0) then
                call jecroc(jexnom('&&IRCAME.PTRNOM', v_nomcmp(icmp1)))
                cmpt = cmpt+1
                zk16(jnocm1+cmpt-1) = v_nomcmp(icmp1)
                zk16(jnocm2+cmpt-1) = v_nomcm2(icmp1)
                call jenonu(jexnom('&&IRCAME.CMPLOC', v_nomcmp(icmp1)), numcmp)
                if (numcmp .ne. 0) then
                    zi(jindir+icmp1-1) = cmpt
                end if
            end if
        end do
        call jedetr('&&IRCAME.TEST')
        call jedetr('&&IRCAME.COUNT')
        call jedetr('&&IRCAME.DISPL')
        call jedetr('&&IRCAME.NOMCMP2')
        call jedetr('&&IRCAME.NOMCMP3')
        call jedetr('&&IRCAME.NOMFAG')
        call jedetr('&&IRCAME.NOMFA2')
        call jedetr('&&IRCAME.PTRNOM')
        call jedetr('&&IRCAME.CMPLOC')
        ncmpvl = ncmpve
        ncmpve = cmpt
    else
        call wkvect(indcmp, 'V V I', ncmpve, jindir)
        do icmp1 = 1, ncmpve
            zi(jindir+icmp1-1) = icmp1
        end do
        ncmpvl = ncmpve
    end if
!   ON REMPLACE LES NOMS DES COMPOSANTES
    if (etiqcp .ne. ' ') then
        call jeveuo(ntncmp, 'L', jnocm1)
        call jeveuo(etiqcp, 'L', jnocm2)
        call jelira(etiqcp, 'LONMAX', nbcmp2)
        nbcmp2 = nbcmp2/2
        do icmp1 = 1, ncmpve
            do icmp2 = 1, nbcmp2
                if (zk16(jnocm1+icmp1-1) .eq. zk16(jnocm2+2*(icmp2-1))) then
                    zk16(jnocm1+icmp1-1) = zk16(jnocm2+2*icmp2-1)
                    goto 10
                end if
            end do
10          continue
        end do
    end if
!   3.2. ==> RECUPERATION DES NB/NOMS/NBNO/NBITEM DES TYPES DE MAILLESDANS CATALOGUE
!            RECUPERATION DES TYPES GEOMETRIE CORRESPONDANT POUR MED
!            VERIF COHERENCE AVEC LE CATALOGUE
    call lrmtyp(nbtyp, nomtyp, nnotyp, typgeo, renumd, &
                modnum, nuanom, numnoa)
    nosdpr = '&&IRCAME.PROFIL_SIZE'
!   3.3. ==> DEFINITIONS DES IMPRESSIONS ET CREATION DES PROFILS EVENTUELS
    call ircmpr(nofimd, typech, nbimpr, ncaimi, ncaimk, &
                ncmprf, ncmpvl, ntlcmp, nbvato, nbenec, &
                lienec, adsd, adsl, nomaas, ligrel, &
                typgeo, nomtyp, ntproa, chanom, sdcarm, &
                field_type, nom_sd_fu, nosdpr)

! - Get number of components
    if (field_type .eq. 'VARI_ELGA') then
        if (nbCmpDyna .eq. 0) then
            nbCmpDyna = ncmpve
        else
            if (nbCmpDyna .gt. ncmpve) then
                call utmess('F', 'MED_31', sk=field_type)
            end if
        end if
    end if
!
    if (nbimpr .gt. 0) then
        call jeveuo(ncaimi, 'E', adcaii)
        call jeveuo(ncaimk, 'E', adcaik)
!       3.4. ==> CARACTERISATION DES SUPPORTS QUAND CE NE SONT PAS DES NOEUDS
        if (typech(1:4) .eq. 'ELGA' .or. typech(1:4) .eq. 'ELEM') then
            if (sdcarm .ne. ' ' .and. typech(1:4) .eq. 'ELGA') then
!
                if (nom_sd_fu .ne. ' ') call utmess('F', 'MED3_3')
!
                call irelst(nofimd, chanom, nochmd, typech, nomaas, &
                            nomamd, nbimpr, zi(adcaii), zk80(adcaik), sdcarm, carael)
            end if
            call irmpga(nofimd, chanom, nochmd, typech, nomtyp, &
                        nbimpr, zi(adcaii), zk80(adcaik), modnum, nuanom, &
                        sdcarm, carael, lfichUniq, field_type, codret)
        end if
!
!       Reperage du champ : existe-t-il deja ?
!       on doit parcourir toutes les impressions possibles pour ce champ
        existc = 0
        do nrimpr = 1, nbimpr
            if (codret .eq. 0) then
                tygeom = zi(adcaii+10*nrimpr-2)
                if (tygeom .eq. typnoe) then
                    typent = ednoeu
                else
                    if (typech .eq. 'ELNO') then
                        typent = ednoma
                    else
                        typent = edmail
                    end if
                end if
                nvalec = zi(adcaii+10*nrimpr-4)
                call jedetr(nmcmfi)
                ifimed = 0
                call mdexch(nofimd, ifimed, nochmd, numpt, numord, &
                            ncmpve, ntncmp, nvalec, typent, tygeom, &
                            jaux, nbcmfi, nmcmfi, nbval, nbprof, &
                            codret)
                existc = max(existc, jaux)
            end if
        end do
!
!       ECRITURE SI C'EST POSSIBLE
        if (existc .le. 2) then
            call ircam1(nofimd, nochmd, existc, ncmprf, numpt, &
                        instan, numord, adsd, adsv, adsl, &
                        adsk, partie, indcmp, ncmpve, ntlcmp, &
                        ntncmp, ntucmp, ntproa, nbimpr, zi(adcaii), &
                        zk80(adcaik), typech, nomamd, nomtyp, modnum, &
                        numnoa, lfichUniq, nom_sd_fu, dschmd, nosdpr, &
                        codret)
        else
            call utmess('F', 'MED2_4', sk=nochmd, sr=instan)
        end if
    else
        call utmess('A', 'MED_82', sk=nochmd)
    end if
    call jedetr(nosdpr)
!
999 continue
    if (nivinf .gt. 1) then
        call cpu_time(end_time)
        write (ifm, *) '<', nompro, '> FIN ECRITURE DU CHAMP MED EN ', &
            end_time-start_time, "sec."
        write (ifm, *) ' '
    end if
!
    call jedetr('&&'//nompro//'.LISTE_N0MCMP   ')
    call jedetr('&&'//nompro//'.NOMCMP         ')
    call jedetr('&&'//nompro//'.UNITECMP       ')
    call jedetr('&&'//nompro//'.PROFIL_ASTER   ')
    call jedetr('&&'//nompro//'.NOMCMP_FICHIER ')
    call jedetr('&&'//nompro//'.CARAC_NOMBRES__')
    call jedetr('&&'//nompro//'.CARAC_CHAINES__')
    call jedetr(indcmp)
!
998 continue
    call jedema()
!
end subroutine
