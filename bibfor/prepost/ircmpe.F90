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
subroutine ircmpe(nofimd, ncmpve, numcmp, exicmp, nbvato, &
                  nbmaec, limaec, adsd, adsl, nbimpr, &
                  ncaimi, ncaimk, tyefma, typmai, typgeo, &
                  nomtyp, typech, profas, promed, prorec, &
                  nroimp, chanom, sdcarm, field_type, nosdfu, &
                  nosdpr)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "asterc/asmpi_allgather_i.h"
#include "asterc/asmpi_allgatherv_char16.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/celfpg.h"
#include "asterfort/cesexi.h"
#include "asterfort/infniv.h"
#include "asterfort/ircael.h"
#include "asterfort/ircmpf.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/wkvect.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: nbvato, ncmpve, numcmp(ncmpve), nbmaec, typmai(*), adsd
    integer(kind=8) :: limaec(*), nbimpr, typgeo(*), profas(nbvato), tyefma(*), start
    integer(kind=8) :: nroimp(nbvato), promed(nbvato), prorec(nbvato), adsl, nbtot, nbloc
    character(len=*) :: nofimd, nosdpr
    character(len=8) :: nomtyp(*), typech, sdcarm, nosdfu
    character(len=19) :: chanom
    character(len=24) :: ncaimi, ncaimk
    aster_logical :: exicmp(nbvato)
    character(len=16), intent(in) :: field_type
!
! --------------------------------------------------------------------------------------------------
!
!     ECRITURE D'UN CHAMP - MAILLES ET PROFIL SUR LES ELEMENTS
!
! --------------------------------------------------------------------------------------------------
!
!   ENTREES :
!       NOFIMD      : NOM DU FICHIER MED
!       NCMPVE      : NOMBRE DE COMPOSANTES VALIDES EN ECRITURE
!       NUMCMP      : NUMEROS DES COMPOSANTES VALIDES
!       EXICMP      : EXISTENCE DES COMPOSANTES PAR MAILLES
!       NBVATO      : NOMBRE DE VALEURS TOTALES
!       NBMAEC      : NOMBRE D'ENTITES A ECRIRE (O, SI TOUTES)
!       LIMAEC      : LISTE DES ENTITES A ECRIRE SI EXTRAIT
!       ADSD,ADSL   : ADRESSES DES TABLEAUX DES CHAMPS SIMPLIFIES
!       TYPMAI      : TYPE ASTER POUR CHAQUE MAILLE
!       TYEFMA      : NRO D'ELEMENT FINI OU DE MAILLE ASSOCIE A CHAQUE MAILLE
!       TYPGEO      : TYPE GEOMETRIQUE DE MAILLE ASSOCIEE AU TYPE ASTER
!       NOMTYP      : NOM DES TYPES DE MAILLES ASTER
!       PROREC      : PROFIL RECIPROQUE. AUXILIAIRE.
!       field_type  : type of field (symbolic name in result datastructure)
!
!   SORTIES :
!       NBIMPR : NOMBRE D'IMPRESSIONS
!       NCAIMI : STRUCTURE ASSOCIEE AU TABLEAU CAIMPI
!               CAIMPI : ENTIERS POUR CHAQUE IMPRESSION
!                   CAIMPI(1,I) = TYPE D'EF / MAILLE ASTER (0, SI NOEUD)
!                   CAIMPI(2,I) = NOMBRE DE POINTS (GAUSS OU NOEUDS)
!                   CAIMPI(3,I) = NOMBRE DE SOUS-POINTS
!                   CAIMPI(4,I) = NOMBRE DE COUCHES
!                   CAIMPI(5,I) = NOMBRE DE SECTEURS
!                   CAIMPI(6,I) = NOMBRE DE FIBTRES
!                   CAIMPI(7,I) = NOMBRE DE MAILLES A ECRIRE
!                   CAIMPI(8,I) = TYPE DE MAILLES ASTER (0, SI NOEUD)
!                   CAIMPI(9,I) = TYPE GEOMETRIQUE AU SENS MED
!                   CAIMPI(10,I) = NOMBRE TOTAL DE MAILLES IDENTIQUES
!       NCAIMK : STRUCTURE ASSOCIEE AU TABLEAU CAIMPK
!               CAIMPK : CARACTERES POUR CHAQUE IMPRESSION
!                   CAIMPK(1,I) = NOM DE LA LOCALISATION ASSOCIEE
!                   CAIMPK(2,I) = NOM DU PROFIL AU SENS MED
!                   CAIMPK(3,I) = NOM DE L'ELEMENT DE STRUCTURE
!       PROFAS : PROFIL ASTER.
!                C'EST LA LISTE DES NUMEROS ASTER DES ELEMENTS POUR LESQUELS LE CHAMP EST DEFINI
!       PROMED : PROFIL MED.
!                C'EST LA LISTE DES NUMEROS MED DES ELEMENTS POUR LESQUELS LE CHAMP EST DEFINI
!       NROIMP : NUMERO DE L'IMPRESSION ASSOCIEE A CHAQUE MAILLE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=80), parameter :: ednopf = ' ', ednoga = ' '
    integer(kind=8), parameter :: nmaxfi = 10
    integer(kind=8) :: nugrfi(nmaxfi), nugrf2(nmaxfi)
    integer(kind=8) :: nmaty0(MT_NTYMAX), adraux(MT_NTYMAX)
    integer(kind=8) :: ifm, niv, ibid, iret, i_fpg, jaux, kaux, ima, jcesd, laux
    integer(kind=8) :: jcesc, jcesl, jcesv, nrefma
    integer(kind=8) :: nrcmp, nrpg, nrsp, nbpg, nbsp, nb_fpg, typmas, nbimp0, nrimpr
    integer(kind=8) :: nbtcou, nbqcou, nbsec, nbfib, jnbimpr, jnbimpr3
    integer(kind=8) :: adcaii, adcaik, nbgrf, nbmaect, nbimprl(1), nb_cmp_tot
    integer(kind=8) :: nbgrf2, nbtcou2, nbqcou2, nbsec2, nbfib2, ima2, nbimprt
    integer(kind=8) :: igrfi, imafib, rang, nbproc, jma, jnbma, ityp, nbmal, iaux
    integer(kind=8), pointer :: v_nb_ma(:) => null()
    integer(kind=8), pointer :: v_nbimpr2(:) => null()
    integer(kind=8), pointer :: v_indir(:) => null()
    character(len=16) :: nomfpg
    character(len=64) :: noprof
    aster_logical :: exicar, grfidt, elga_sp, okgrcq, oktuy
    character(len=16), pointer :: fpg_name(:) => null()
    character(len=16), pointer :: nofpgma(:) => null()
    character(len=16), pointer :: v_nofpg(:) => null()
    character(len=16), pointer :: v_nofpg_l(:) => null()
    mpi_int :: mrank, msize, mpicou, taille
    aster_logical :: lficUniq, lnbmal, lnbmaec, cara_ele_used
    real(kind=8) :: start_time, end_time
    mpi_int, pointer :: v_count(:) => null()
    mpi_int, pointer :: v_displ(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    !
    nmaty0(:) = 0
    !
    if (niv .ge. 2) then
        call cpu_time(start_time)
        write (ifm, 805) 'DEBUT DE IRCMPE'
    end if
    !
    ! ON REMPLIT UN PREMIER TABLEAU PAR MAILLE :
    !   * VRAI DES QU'UNE DES COMPOSANTES DU CHAMP EST PRESENTE SUR LA MAILLE
    !   * FAUX SINON
    ! REMARQUE : ON EXAMINE LES NCMPVE COMPOSANTES QUI SONT DEMANDEES,
    !    MAIS IL FAUT BIEN TENIR COMPTE DE LA NUMEROTATION DE REFERENCE
    laux = adsd+1
    cifpg: do i_fpg = 1, nbvato
        laux = laux+4
        nbpg = zi(laux)
        nbsp = zi(laux+1)
        do nrcmp = 1, ncmpve
            kaux = numcmp(nrcmp)
            do nrpg = 1, nbpg
                do nrsp = 1, nbsp
                    call cesexi('C', adsd, adsl, i_fpg, nrpg, &
                                nrsp, kaux, jaux)
                    if (jaux .gt. 0) then
                        exicmp(i_fpg) = .true.
                        cycle cifpg
                    end if
                end do
            end do
        end do
    end do cifpg
    !
    ! PROFAS : LISTE DES MAILLES POUR LESQUELS ON AURA IMPRESSION
    !    UNE MAILLE EST PRESENTE SI ET SEULEMENT SI AU MOINS UNE COMPOSANTE
    !    Y EST DEFINIE ET SI ELLE FAIT PARTIE DU FILTRAGE DEMANDE
    !
    nb_fpg = 0
    lficUniq = .false._1
    if (nosdfu .ne. ' ') then
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
        call jeveuo(nosdfu//'.MAIL', 'L', jma)
        lficUniq = .true._1
    end if
    if (lficUniq) then
        nbmaect = nbmaec
        call asmpi_comm_vect('MPI_SUM', 'I', nbval=1, sci=nbmaect)
        lnbmaec = nbmaect .eq. 0
    else
        lnbmaec = nbmaec .eq. 0
    end if
    !
    ! SANS FILTRAGE : C'EST LA LISTE DES MAILLES QUI POSSEDENT UNE COMPOSANTE VALIDE
    if (lnbmaec) then
        do i_fpg = 1, nbvato
            if (exicmp(i_fpg)) then
                if (lficUniq) then
                    if (zi(jma+i_fpg-1) .gt. 0) then
                        nb_fpg = nb_fpg+1
                        profas(nb_fpg) = i_fpg
                    end if
                else
                    nb_fpg = nb_fpg+1
                    profas(nb_fpg) = i_fpg
                end if
            end if
        end do
        !
        ! AVEC FILTRAGE : C'EST LA LISTE DES MAILLES REQUISES ET AVEC UNE COMPOSANTE VALIDE
    else
        do jaux = 1, nbmaec
            i_fpg = limaec(jaux)
            if (exicmp(i_fpg)) then
                if (lficUniq) then
                    if (zi(jma+i_fpg-1) .gt. 0) then
                        nb_fpg = nb_fpg+1
                        profas(nb_fpg) = i_fpg
                    end if
                else
                    nb_fpg = nb_fpg+1
                    profas(nb_fpg) = i_fpg
                end if
            end if
        end do
    end if
    !
    ! CARACTERISATIONS DES IMPRESSIONS
    !   ON TRIE SELON DEUX CRITERES :
    !       1.  LE NOMBRE DE SOUS-POINTS
    !       2.  LE TYPE D'ELEMENT FINI POUR UN CHAMP ELGA, OU LE TYPE DE LA
    !           MAILLE, POUR UN AUTRE TYPE DE CHAMP. LE TABLEAU TYEFMA VAUT DONC
    !           EFMAI OU TYPMAI A L'APPEL , SELON LE TYPE DE CHAMP.
    !
    ! TABLEAU DES CARACTERISATIONS ENTIERES DES IMPRESSIONS ALLOCATION INITIALE
    nbimp0 = 20
    i_fpg = 10*nbimp0
    call wkvect(ncaimi, 'V V I', i_fpg, adcaii)
    !
    ! PARCOURS DES MAILLES QUI PASSENT LE FILTRE
    nbimpr = 0
    ! SI ON EST SUR UN CHAMP ELGA, LE TRI DOIT SE FAIRE SUR LES FAMILLES DE POINTS DE GAUSS
    if (typech(1:4) .eq. 'ELGA') then
        call celfpg(chanom, '&&IRCMPE.NOFPGMA', ibid)
        if (nb_fpg .ne. 0) then
            AS_ALLOCATE(vk16=fpg_name, size=nb_fpg)
        end if
        call jeveuo('&&IRCMPE.NOFPGMA', 'L', vk16=nofpgma)
    end if
    !
    call jeexin(sdcarm//'.CANBSP    .CESV', iret)
    exicar = .false.
    if (iret .ne. 0 .and. typech(1:4) .eq. 'ELGA') then
        call jeveuo(sdcarm//'.CANBSP    .CESD', 'L', jcesd)
        call jeveuo(sdcarm//'.CANBSP    .CESC', 'L', jcesc)
        call jeveuo(sdcarm//'.CANBSP    .CESL', 'L', jcesl)
        call jeveuo(sdcarm//'.CANBSP    .CESV', 'L', jcesv)
        exicar = .true.
    end if
    !
    cara_ele_used = .false.
    do i_fpg = 1, nb_fpg
        ima = profas(i_fpg)
        nrefma = tyefma(ima)
        !
        laux = adsd+4*ima+1
        nbpg = zi(laux)
        ! Les types de mailles ont été modifiés avant pour les cas particuliers
        ! NE TIENT PAS COMPTE DES CAS PARTICULIER DE IRCMPR (PENTA21, HEXA9, ..)
        nbsp = zi(laux+1)
        if (typech(1:4) .eq. 'ELNO') then
            ! For HEXA9 (COQUE_SOLIDE element)
            if (nbpg .eq. 9) then
                if (typmai(ima) .eq. MT_HEXA8) then
                    nbpg = 8
                end if
            end if
            ! For PENTA21
            if (nbpg .eq. 21) then
                if (typmai(ima) .eq. MT_PENTA18) then
                    nbpg = 18
                end if
            end if
            ! For PYRAM19
            if (nbpg .eq. 19) then
                if (typmai(ima) .eq. MT_PYRAM13) then
                    nbpg = 13
                end if
            end if
            ! For TETRA15
            if (nbpg .eq. 15) then
                if (typmai(ima) .eq. MT_TETRA10) then
                    nbpg = 10
                end if
            end if
            ! For TETRA20
            if (nbpg .eq. 20) then
                if (typmai(ima) .eq. MT_TETRA4) then
                    nbpg = 4
                end if
            end if
            ! For QUAD12
            if (nbpg .eq. 12) then
                if (typmai(ima) .eq. MT_QUAD4) then
                    nbpg = 4
                end if
            end if
            ! For TRIA10
            if (nbpg .eq. 10) then
                if (typmai(ima) .eq. MT_TRIA3) then
                    nbpg = 3
                end if
            end if
        end if

        nomfpg = 'a fac'
        elga_sp = .false.
        imafib = 0
        if (typech(1:4) .eq. 'ELGA') then
            nomfpg = nofpgma(ima)
            if ((nbsp .ge. 1) .and. exicar) then
                call ircael(jcesd, jcesl, jcesv, jcesc, ima, &
                            nbqcou, nbtcou, nbsec, nbfib, nbgrf, nugrfi)
                if (nbsp .gt. 1 .or. nbqcou .eq. 1 .or. nbfib .eq. 1) then
                    elga_sp = .true.
                    if (nbfib .ne. 0) imafib = ima
                end if
            end if
        end if
        !
        ! RECHERCHE D'UNE IMPRESSION SEMBLABLE
        do jaux = 1, nbimpr
            if (typech(1:4) .eq. 'ELGA') then
                ! Pour les ELGA, tri sur les familles de points de gauss
                if (.not. elga_sp) then
                    ! si on n'a pas de cara_elem, le cas est simple
                    ! on compare le nom de la famille de pg et nbsp
                    if (fpg_name(jaux) .eq. nomfpg .and. &
                        zi(adcaii+10*(jaux-1)+2) .eq. nbsp) then
                        nrimpr = jaux
                        goto 423
                    end if
                else
                    if (fpg_name(jaux) .eq. nomfpg) then
                        ! Les différents cas : PMF, Coques, Tuyaux, Grille
                        if (nbfib .ne. 0) then
                            ! Pour les PMF, on compare aussi les groupes de fibres
                            ima2 = zi(adcaii+10*(jaux-1)+5)
                            call ircael(jcesd, jcesl, jcesv, jcesc, ima2, &
                                        nbqcou2, nbtcou2, nbsec2, nbfib2, nbgrf2, nugrf2)
                            if ((nbfib2 .eq. nbfib) .and. (nbgrf2 .eq. nbgrf)) then
                                grfidt = .true.
                                do igrfi = 1, nbgrf2
                                    if (nugrf2(igrfi) .ne. nugrfi(igrfi)) grfidt = .false.
                                end do
                                if (grfidt) then
                                    nrimpr = jaux
                                    goto 423
                                end if
                            end if
                        else
                            ! Coque ou grille
                            okgrcq = (nbqcou .ge. 1) .and. (nbsp .ge. 1)
                            ! tuyau
                            oktuy = (nbtcou .ge. 1) .and. (nbsec .ge. 1)
                            !
                            kaux = adcaii+10*(jaux-1)
                            if (oktuy .and. (zi(kaux+3) .eq. nbtcou) .and. &
                                (zi(kaux+4) .eq. nbsec)) then
                                ! Tuyaux : même nb de couche et de secteur ==> même nbsp
                                nrimpr = jaux
                                goto 423
                            else if (okgrcq .and. (zi(kaux+3) .eq. nbqcou) .and. &
                                     (zi(kaux+2) .eq. nbsp) .and. &
                                     (zi(kaux+1) .eq. nbpg)) then
                                ! Coques ou Grilles : même nb de couche, de sous-point, de pt gauss
                                !   Coques  nbsp = 3*nbqcou
                                !   Grilles nbqcou=nbsp=1
                                nrimpr = jaux
                                goto 423
                            else
                                ! autre cas
                                ! on compare le nom de la famille de pg et nbsp
                                if (zi(adcaii+10*(jaux-1)+2) .eq. nbsp) then
                                    nrimpr = jaux
                                    goto 423
                                end if
                            end if
                        end if
                    end if
                end if
            else
                if (zi(adcaii+10*(jaux-1)) .eq. nrefma .and. &
                    zi(adcaii+10*(jaux-1)+2) .eq. nbsp) then
                    nrimpr = jaux
                    goto 423
                end if
            end if
        end do
        !
        ! ON CREE UNE NOUVELLE IMPRESSION SI ON DEPASSE LA LONGUEUR RESERVEE, ON DOUBLE
        if (nbimpr .eq. nbimp0) then
            nbimp0 = 2*nbimp0
            call juveca(ncaimi, 10*nbimp0)
            call jeveuo(ncaimi, 'E', adcaii)
        end if
        !
        nbimpr = nbimpr+1
        jaux = adcaii+10*(nbimpr-1)
        ! CAIMPI(1,I) = TYPE D'EF / MAILLE ASTER (0, SI NOEUD)
        zi(jaux) = nrefma
        ! CAIMPI(2,I) = NOMBRE DE POINTS (DE GAUSS OU NOEUDS)
        zi(jaux+1) = nbpg
        ! CAIMPI(3,I) = NOMBRE DE SOUS-POINTS
        zi(jaux+2) = nbsp
        ! CAIMPI(4,I) = NOMBRE DE COUCHES  : Tuyaux ou Coques-Grilles
        ! CAIMPI(5,I) = NOMBRE DE SECTEURS : Tuyaux
        if (elga_sp) then
            if (imafib .eq. 0) then
                ! Tuyaux
                if ((nbsec .ge. 1) .and. (nbtcou .ge. 1)) then
                    zi(jaux+3) = nbtcou
                    zi(jaux+4) = nbsec
                    if (nbsp .ne. (2*nbsec+1)*(2*nbtcou+1)) then
                        call utmess('F', 'MED2_12', sk=field_type)
                    end if
                    ! Coques
                else if ((nbqcou .ge. 1) .and. (nbsp .eq. 3*nbqcou)) then
                    zi(jaux+3) = nbqcou
                    zi(jaux+4) = 0
                    ! Grilles
                else if ((nbqcou .eq. 1) .and. (nbsp .eq. 1)) then
                    zi(jaux+3) = nbqcou
                    zi(jaux+4) = 0
                else
                    ! Ce n'est pas un élément à sous-point
                    zi(jaux+3) = 0
                    zi(jaux+4) = 0
                end if
            else
                zi(jaux+3) = 0
                zi(jaux+4) = 0
                if (nbsp .ne. nbfib) then
                    call utmess('F', 'MED2_12', sk=field_type)
                end if
            end if
            cara_ele_used = .true.
        else
            zi(jaux+3) = 0
            zi(jaux+4) = 0
        end if
        ! CAIMPI(6,I) = NUMERO DE LA MAILLE 'EXEMPLE' POUR LES PMF
        zi(jaux+5) = imafib
        ! CAIMPI(7,I) = NOMBRE DE MAILLES A ECRIRE
        zi(jaux+6) = 0
        ! CAIMPI(8,I) = TYPE DE MAILLES ASTER (0, SI NOEUD)
        zi(jaux+7) = typmai(ima)
        ! CAIMPI(9,I) = TYPE GEOMETRIQUE AU SENS MED
        zi(jaux+8) = typgeo(typmai(ima))
        !
        if (typech(1:4) .eq. 'ELGA') then
            fpg_name(nbimpr) = nomfpg
        end if
        nrimpr = nbimpr
        !
        ! MEMORISATION DE L'IMPRESSION DE LA MAILLE, CUMUL DU NOMBRE DE MAILLES POUR L'IMPRESSION
423     continue
        !
        nroimp(ima) = nrimpr
        jaux = adcaii+10*(nrimpr-1)+6
        zi(jaux) = zi(jaux)+1
    end do

    if (exicar .and. .not. cara_ele_used) then
        call utmess('A', 'MED2_13', sk=field_type)
    end if
    !
    ! Le but du bloc suivant est d'avoir les memes impressions a realiser
    ! sur tous les procs (quite à ce que certaines soient vides)
    ! Il faut donc communiquer pour savoir les impressions de chaque procs
    if (lficUniq) then
        call asmpi_comm('GET', mpicou)
        nbimprt = nbimpr
        call asmpi_comm_vect('MPI_SUM', 'I', 1, 0, sci=nbimprt)
        nbimprl(1) = nbimpr
        call wkvect('&&IRCMPE.NBIMPR', 'V V I', nbproc+1, jnbimpr)
        taille = 1
        call asmpi_allgather_i(nbimprl, taille, zi(jnbimpr+1), taille, mpicou)
        call wkvect('&&IRCMPE.COUNT', 'V V S', nbproc, vi4=v_count)
        call wkvect('&&IRCMPE.DISPL', 'V V S', nbproc+1, vi4=v_displ)
        nb_cmp_tot = 0
        do iaux = 2, nbproc+1
            nb_cmp_tot = nb_cmp_tot+zi(jnbimpr+iaux-1)
            v_count(iaux-1) = to_mpi_int(zi(jnbimpr+iaux-1))
            v_displ(iaux) = to_mpi_int(nb_cmp_tot)
            zi(jnbimpr+iaux-1) = zi(jnbimpr+iaux-1)+zi(jnbimpr+iaux-2)
        end do
        zi(jnbimpr) = 0
        call wkvect('&&IRCMPE.NBIMPR2', 'V V I', 10*nbimprt, vi=v_nbimpr2)
        do iaux = 1, nbimpr
            do jaux = 0, 9
                v_nbimpr2(10*(zi(jnbimpr+rang)+iaux-1)+jaux+1) = zi(adcaii+10*(iaux-1)+jaux)
            end do
        end do
        call jedetr('&&IRCMPE.NBIMPR')
        call wkvect('&&IRCMPE.NBIMPR', 'V V I', nbimprt, jnbimpr)
        call wkvect('&&IRCMPE.NBIMPR3', 'V V I', 10*nbimprt, jnbimpr3)
        call wkvect('&&IRCMPE.NBIMPR4', 'V V K16', 10*nbimprt, vk16=v_nofpg_l)
        call asmpi_comm_vect('MPI_SUM', 'I', 10*nbimprt, 0, vi=v_nbimpr2)
        taille = to_mpi_int(nbimpr)
        call wkvect('&&IRCMPE.NOFPG', 'V V K16', nb_cmp_tot, vk16=v_nofpg)
        if (typech .eq. 'ELGA') then
            call asmpi_allgatherv_char16(fpg_name, taille, v_nofpg, v_count, v_displ, mpicou)
        end if
        call jedetr('&&IRCMPE.COUNT')
        call jedetr('&&IRCMPE.DISPL')
        nbimprl(1) = 0
        do iaux = 1, nbimprt
            if (zi(jnbimpr+iaux-1) .eq. -1) cycle
            nomfpg = v_nofpg(iaux)
            nrefma = v_nbimpr2(10*(iaux-1)+1)
            nbsp = v_nbimpr2(10*(iaux-1)+3)
            do jaux = iaux+1, nbimprt
                if (typech .eq. 'ELGA') then
                    if (nomfpg .eq. v_nofpg(jaux) .and. &
                        nbsp .eq. v_nbimpr2(10*(jaux-1)+3)) then
                        zi(jnbimpr+jaux-1) = -1
                        cycle
                    end if
                else
                    if (nrefma .eq. v_nbimpr2(10*(jaux-1)+1) .and. &
                        nbsp .eq. v_nbimpr2(10*(jaux-1)+3)) then
                        zi(jnbimpr+jaux-1) = -1
                        cycle
                    end if
                end if
            end do
            do jaux = 0, 9
                zi(jnbimpr3+10*nbimprl(1)+jaux) = v_nbimpr2(10*(iaux-1)+jaux+1)
            end do
            zi(jnbimpr3+10*nbimprl(1)+6) = 0
            v_nofpg_l(nbimprl(1)+1) = nomfpg
            nbimprl(1) = nbimprl(1)+1
        end do
        if (nbimpr .ne. 0) then
            call wkvect('&&IRCMPE.INDIR', 'V V I', nbimpr, vi=v_indir)
        end if
        do iaux = 1, nbimprl(1)
            nomfpg = v_nofpg_l(iaux)
            nrefma = zi(jnbimpr3+10*(iaux-1))
            nbsp = zi(jnbimpr3+10*(iaux-1)+2)
            nbpg = zi(jnbimpr3+10*(iaux-1)+1)
            do jaux = 1, nbimpr
                if (typech .eq. 'ELGA') then
                    if (nomfpg .eq. fpg_name(jaux) .and. &
                        nbsp .eq. zi(adcaii+10*(jaux-1)+2)) then
                        zi(jnbimpr3+10*(iaux-1)+6) = zi(adcaii+10*(jaux-1)+6)
                        ASSERT(nbimpr .ne. 0)
                        v_indir(jaux) = iaux
                        cycle
                    end if
                else
                    if (nrefma .eq. zi(adcaii+10*(jaux-1)) .and. &
                        nbsp .eq. zi(adcaii+10*(jaux-1)+2)) then
                        zi(jnbimpr3+10*(iaux-1)+6) = zi(adcaii+10*(jaux-1)+6)
                        ASSERT(nbimpr .ne. 0)
                        v_indir(jaux) = iaux
                        cycle
                    end if
                end if
            end do
        end do
        call jedetr(ncaimi)
        nbimpr = nbimprl(1)
        call wkvect(ncaimi, 'V V I', 10*nbimpr, adcaii)
        do iaux = 1, 10*nbimprl(1)
            zi(adcaii+iaux-1) = zi(jnbimpr3+iaux-1)
        end do
        call jedetr('&&IRCMPE.NBIMPR')
        call jedetr('&&IRCMPE.NBIMPR2')
        call jedetr('&&IRCMPE.NBIMPR3')
        call jedetr('&&IRCMPE.NBIMPR4')
        call jedetr('&&IRCMPE.NOFPG')
    else
        nbimprt = nbimpr
        if (nbimpr .ne. 0) then
            call wkvect('&&IRCMPE.INDIR', 'V V I', nbimpr, vi=v_indir)
            do jaux = 1, nbimpr
                v_indir(jaux) = jaux
            end do
        end if
    end if

    !
    if (typech(1:4) .eq. 'ELGA') then
        AS_DEALLOCATE(vk16=fpg_name)
        call jedetr('&&IRCMPE.NOFPGMA')
    end if
    if (nbimprt .eq. 0) then
        goto 999
    end if
    !
    ! CONVERSION DU PROFIL EN NUMEROTATION MED
    !   PROMED : ON STOCKE LES VALEURS DES NUMEROS DES MAILLES AU SENS MED PAR TYPE DE MAILLES.
    !            IL FAUT REORDONNER LE TABLEAU PROFAS PAR IMPRESSION SUCCESSIVE :
    !               LE TABLEAU EST ORGANISE EN SOUS-TABLEAU CORRESPONDANT A CHAQUE IMPRESSION.
    !               ON REPERE CHAQUE DEBUT DE SOUS-TABLEAU AVEC ADRAUX.
    !
    !   PROREC : C'EST LA LISTE RECIPROQUE.
    !            POUR LA MAILLE NUMERO IAUX EN NUMEROTATION ASTER, ON A SA POSITION DANS LE
    !            TABLEAU DES VALEURS S'IL FAIT PARTIE DE LA LISTE SINON C'EST 0.
    do i_fpg = 1, nb_fpg
        ima = profas(i_fpg)
        prorec(ima) = i_fpg
    end do
    !
    ! ADRESSES DANS LE TABLEAU PROFAS
    !   ADRAUX(IAUX) = ADRESSE DE LA FIN DE LA ZONE DE L'IMPRESSION PRECEDENTE, IAUX-1
    adraux(1) = 0
    ASSERT(nbimpr .le. MT_NTYMAX)
    do i_fpg = 2, nbimpr
        adraux(i_fpg) = adraux(i_fpg-1)+zi(adcaii+10*(i_fpg-2)+6)
    end do
    !
    ! DECOMPTE DU NOMBRE DE MAILLES PAR TYPE DE MAILLES ASTER
    !   NMATY0(IAUX) =  NUMERO MED DE LA MAILLE COURANTE, DANS LA CATEGORIE ASTER IAUX.
    !                   A LA FIN, NMATY0(IAUX) VAUT LE NOMBRE DE MAILLES PAR TYPE DE MAILLES ASTER,
    !                   POUR TOUTES LES MAILLES DU MAILLAGE
    !   ADRAUX(JAUX) =  ADRESSE DANS LES TABLEAUX PROMED ET PROFAS DE LA MAILLE COURANTE ASSOCIEE
    !                   A L'IMPRESSION NUMERO JAUX
    if (.not. lficUniq) then
        do ima = 1, nbvato
            typmas = typmai(ima)
            nmaty0(typmas) = nmaty0(typmas)+1
            if (prorec(ima) .ne. 0) then
                jaux = v_indir(nroimp(ima))
                adraux(jaux) = adraux(jaux)+1
                ASSERT(adraux(jaux) .le. nbvato)
                promed(adraux(jaux)) = nmaty0(typmas)
                profas(adraux(jaux)) = ima
            end if
        end do
    else
        do ima = 1, nbvato
            typmas = typmai(ima)
            nmaty0(typmas) = nmaty0(typmas)+1
            if (zi(jma+ima-1) .gt. 0) then
                if (prorec(ima) .ne. 0) then
                    jaux = v_indir(nroimp(ima))
                    adraux(jaux) = adraux(jaux)+1
                    ASSERT(adraux(jaux) .le. nbvato)
                    promed(adraux(jaux)) = nmaty0(typmas)
                    profas(adraux(jaux)) = ima
                end if
            end if
        end do
    end if
    call jedetr('&&IRCMPE.INDIR')
    !
    ! MEMORISATION DANS LES CARACTERISTIQUES DE L'IMPRESSION
    !
    ! NOMBRE DE MAILLES DU MEME TYPE
    do i_fpg = 1, nbimpr
        jaux = adcaii+10*(i_fpg-1)
        typmas = zi(jaux+7)
        ! CAIMPI(10,I) = NOMBRE DE MAILLES IDENTIQUES
        zi(jaux+9) = nmaty0(typmas)
    end do
    !
    !CARACTERISTIQUES CARACTERES
    i_fpg = 3*nbimpr
    call wkvect(ncaimk, 'V V K80', i_fpg, adcaik)
    do i_fpg = 1, nbimpr
        jaux = adcaik+2*(i_fpg-1)
        ! CAIMPK(1,I) = NOM DE LA LOCALISATION ASSOCIEE
        zk80(jaux) = ednoga
        ! CAIMPK(2,I) = NOM DU PROFIL AU SENS MED
        zk80(jaux+1) = ednopf
        ! CAIMPK(3,I) = NOM DE L'ELEMENT DE STRUCTURE
        zk80(jaux+2) = ednopf
    end do
    !
    ! STOCKAGE DES EVENTUELS PROFILS DANS LE FICHIER MED
    kaux = 1
    if (lficUniq) then
        call jeveuo(nosdfu//'.MATY', 'L', jnbma)
    else
        jnbma = 0
    end if
    !
    call wkvect(nosdpr, "V V I", 3*nbimpr, vi=v_nb_ma)
    do i_fpg = 1, nbimpr
        jaux = adcaii+10*(i_fpg-1)
        ! SI LE NOMBRE DE MAILLES A ECRIRE EST >0
        if (zi(jaux+6) .gt. 0 .or. lficUniq) then
            if (lficUniq) then
                ityp = zi(jaux+7)
                nbmal = zi(jnbma+3*(ityp-1)+1)
                iaux = 0
                if (nbmal .ne. zi(jaux+6)) iaux = 1
                call asmpi_comm_vect('MPI_MAX', 'I', 1, 0, sci=iaux)
                lnbmal = .false.
                if (iaux .eq. 1) lnbmal = .true._1
            else
                lnbmal = (zi(jaux+6) .ne. zi(jaux+9))
                ityp = 0
            end if
            ! SI LE NOMBRE DE MAILLES A ECRIRE EST DIFFERENT
            ! DU NOMBRE TOTAL DE MAILLES DE MEME TYPE:
            if (lnbmal) then
                if (lficUniq) then
                    call ircmpf(nofimd, zi(jaux+6), profas(kaux), noprof, nosdfu, &
                                1, ityp, nbtot, nbloc, start)
                else
                    call ircmpf(nofimd, zi(jaux+6), promed(kaux), noprof, nosdfu, &
                                1, ityp, nbtot, nbloc, start)
                end if
                ! CAIMPK(2,I) = NOM DU PROFIL AU SENS MED
                zk80(adcaik+3*(i_fpg-1)+1) = noprof
                v_nb_ma(3*i_fpg-2) = nbtot
                v_nb_ma(3*i_fpg-1) = nbloc
                v_nb_ma(3*i_fpg) = start
            end if
            ! KAUX := POINTEUR PERMETTANT DE SE PLACER DANS PROMED POUR LA PROCHAINE IMPRESSION
            kaux = kaux+zi(jaux+6)
        end if
    end do
    !
    if (niv .ge. 2) then
        if (typech(1:4) .eq. 'ELGA') then
            write (ifm, 801)
        else
            write (ifm, 804)
        end if
        do i_fpg = 1, nbimpr
            jaux = adcaii+10*(i_fpg-1)
            if (zi(jaux+6) .gt. 0) then
                write (ifm, 802) nomtyp(zi(jaux+7)), zi(jaux+6), zi(jaux+1), zi(jaux+2)
            end if
        end do
        write (ifm, 803)
        write (ifm, 805) 'FIN DE IRCMPE'
        call cpu_time(end_time)
        write (ifm, *) "=========== EN ", end_time-start_time, "sec"
    end if
801 format(4x, 65('*'), /, 4x, '*  TYPE DE *', 22x, 'NOMBRE DE', 21x, '*',&
     &/, 4x, '*  MAILLE  *  VALEURS   * POINT(S) DE GAUSS *',&
     &     '   SOUS_POINT(S)   *', /, 4x, 65('*'))
802 format(4x, '* ', a8, ' *', i11, ' *', i15, '    *', i15, '    *')
803 format(4x, 65('*'))
804 format(4x, 65('*'), /, 4x, '*  TYPE DE *', 22x, 'NOMBRE DE', 21x, '*',&
     &/, 4x, '*  MAILLE  *  VALEURS   *      POINTS       *',&
     &     '   SOUS_POINT(S)   *', /, 4x, 65('*'))
805 format(/, 4x, 10('='), a, 10('='),/)
!
999 continue
!
end subroutine
