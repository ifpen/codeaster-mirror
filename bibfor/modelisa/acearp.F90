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

subroutine acearp(nbocc, infdonn, infcarte, grplmax, zjdlm)
!
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM
!
!     AFFECTATION DES CARACTERISTIQUES POUR LES ELEMENTS DISCRET PAR RAIDEUR REPARTIE
!
! --------------------------------------------------------------------------------------------------
!
    use cara_elem_parameter_module
    use cara_elem_info_type
    use cara_elem_carte_type
!
    implicit none
    integer(kind=8)         :: nbocc, zjdlm(*)
    type(cara_elem_info)    :: infdonn
    type(cara_elem_carte)   :: infcarte(*)
    character(len=24)       :: grplmax(*)
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/affdis.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/r8inir.h"
#include "asterfort/rairep.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!&<
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: nbcar, nbval, nrd, GroupeMaxOccur
    parameter(nbcar=100, nbval=12, nrd=2)
    integer(kind=8) :: jdc(3), jdv(3), iunite, ifm, iretour
    integer(kind=8) :: jdcinf, jdvinf
    integer(kind=8) :: ii, idecal, in, inbn, ino, inoe, ioc, irep
    integer(kind=8) :: irgno, isym, itbmp, itbno, iv
    integer(kind=8) :: jd, jj, jn
    integer(kind=8) :: ll, ldgm, ldnm, lokm, lorep, nbnma
    integer(kind=8) :: nbparno, nb_noeu_disc, nc, ncar, ncmp, ndim, ng, ngp, nma, nrep, nval
    integer(kind=8) :: dimcar, nummail, NoeudMaxMaille, MailleMaxOccur
    integer(kind=8) :: vali(2)
    integer(kind=8) :: jdme
! --------------------------------------------------------------------------------------------------
    real(kind=8)      :: val(nbval), eta, vale(nbval), rirot(3)
    character(len=1)  :: kma(3)
    character(len=8)  :: nomnoe, nommai, nomu, car(nbcar), lamass, noma
    character(len=16) :: rep, repdis(nrd)
    character(len=19) :: cart(3), cartdi
    character(len=24) :: nogp
! --------------------------------------------------------------------------------------------------
    aster_logical :: transl, trarot, okunite
!
    data repdis/'GLOBAL          ', 'LOCAL           '/
    data kma/'K', 'M', 'A'/
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    if (nbocc .le. 0) goto 999
!   Récupère les informations du concept
    nomu = infdonn%nomu
    noma = infdonn%maillage
    ndim = infdonn%dimmod
    jdme = infdonn%jmodmail
!   Pour les discrets c'est obligatoirement du 2D ou 3D
    ASSERT((ndim .eq. 2) .or. (ndim .eq. 3))
!   Si c'est un maillage partionné ==> PLOUF
    if ( infdonn%IsParaMesh ) then
        call utmess('F', 'AFFECARAELEM_99', sk='RIGI_PARASOL')
    endif
!
    GroupeMaxOccur = infdonn%GroupeMaxOccur
    MailleMaxOccur = infdonn%MailleMaxOccur
    NoeudMaxMaille = infdonn%NoeudMaxMaille
!
!   Nombre maximum de noeud concernés
!       (nombre maximum de maille par occurrence) * (nombre maximum de noeud pour les mailles)
    nbparno = MailleMaxOccur*NoeudMaxMaille
!
    call wkvect('&&TMPTABNO', 'V V K8', nbparno, itbno)
    call wkvect('&&TMPRIGNO', 'V V R', 6*nbparno, irgno)
    call wkvect('&&TMPTABMP', 'V V K8', nbparno, itbmp)
!
!   Les cartes sont déjà construites : ace_crea_carte
    cartdi = infcarte(ACE_CAR_DINFO)%nom_carte
    jdcinf = infcarte(ACE_CAR_DINFO)%adr_cmp
    jdvinf = infcarte(ACE_CAR_DINFO)%adr_val
    dimcar = infcarte(ACE_CAR_DINFO)%nbr_cmp
!
    cart(1) = infcarte(ACE_CAR_DISCK)%nom_carte
    jdc(1)  = infcarte(ACE_CAR_DISCK)%adr_cmp
    jdv(1)  = infcarte(ACE_CAR_DISCK)%adr_val
!
    cart(2) = infcarte(ACE_CAR_DISCM)%nom_carte
    jdc(2)  = infcarte(ACE_CAR_DISCM)%adr_cmp
    jdv(2)  = infcarte(ACE_CAR_DISCM)%adr_val
!
    cart(3) = infcarte(ACE_CAR_DISCA)%nom_carte
    jdc(3)  = infcarte(ACE_CAR_DISCA)%adr_cmp
    jdv(3)  = infcarte(ACE_CAR_DISCA)%adr_val
!
    ifm = infdonn%ivr(4)
!   Boucle sur les occurrences de rigi_parasol
    do ioc = 1, nbocc
        eta = 0.0d0
!       Par défaut on est dans le repère global, matrices symétriques
        irep = 1; isym = 1; rep = repdis(1)
!
        call getvtx('RIGI_PARASOL', 'GROUP_MA', iocc=ioc, nbval=GroupeMaxOccur, &
                    vect=grplmax, nbret=ng)
        call getvtx('RIGI_PARASOL', 'CARA', iocc=ioc, nbval=nbcar, vect=car, nbret=ncar)
        call getvr8('RIGI_PARASOL', 'VALE', iocc=ioc, nbval=nbval, vect=val, nbret=nval)
        call getvtx('RIGI_PARASOL', 'REPERE', iocc=ioc, scal=rep, nbret=nrep)
        call getvtx('RIGI_PARASOL', 'GROUP_MA_POI1', iocc=ioc, scal=nogp, nbret=ngp)
        if (ngp .eq. 0) then
            call getvtx('RIGI_PARASOL', 'GROUP_MA_SEG2', iocc=ioc, scal=nogp, nbret=ngp)
        end if
        ASSERT(ngp .ne. 0)
        ASSERT(ncar .ge. 1)
!
        if (nrep .ne. 0) then
            do ii = 1, nrd
                if (rep .eq. repdis(ii)) irep = ii
            end do
        end if
!       Unité pour imprimer les valeur des discrets
        iunite = -1
        call getvis('RIGI_PARASOL', 'UNITE', iocc=ioc, scal=iunite, nbret=iretour)
        okunite = (iunite .gt. 0) .and. (ndim .eq. 3)
        if (okunite) then
            write (iunite, 100) rep, ioc
        end if
!       GROUP_MA = toutes les mailles de tous les groupes de mailles
        if (ng .le. 0) goto 30
        idecal = 0
        do nc = 1, ncar
            ! Discrets seulement en translation
            transl = (car(nc)(1:7) .eq. 'K_T_D_N') .or. (car(nc)(1:7) .eq. 'K_T_D_L') .or. &
                     (car(nc)(1:7) .eq. 'A_T_D_N') .or. (car(nc)(1:7) .eq. 'A_T_D_L')
            ! Discrets en translation et rotation
            trarot = (car(nc)(1:8) .eq. 'K_TR_D_N') .or. (car(nc)(1:8) .eq. 'K_TR_D_L') .or. &
                     (car(nc)(1:8) .eq. 'A_TR_D_N') .or. (car(nc)(1:8) .eq. 'A_TR_D_L')
            !
            if (transl .eqv. trarot) then
                call utmess('F', 'MODELISA_17', sk=car(nc))
            end if
            ! Si 2 caractéristiques
            if (nc .eq. 2) then
                ! A et K : pas (K et K) ou (A et A)
                if (car(1)(1:1) .eq. car(2)(1:1)) then
                    call utmess('F', 'MODELISA_16')
                endif
                ! La modélisation doit être [A|K]_LaMême
                if (transl) lokm = 7
                if (trarot) lokm = 8
                if (car(1)(2:lokm) .ne. car(2)(2:lokm)) then
                    call utmess('F', 'MODELISA_16')
                endif
            end if
            !
            vale(:) = 0.0
            if (transl) then
                lamass = 'M'//car(nc)(2:7)
                ! En 3D          1  2  3     1  2  3
                !   ncar=1  :   vx vy vz
                !   ncar=2  :   vx vy vz    vx vy vz
                ! En 2D
                !   ncar=1  :   vx vy
                !   ncar=2  :   vx vy       vx vy
                if ( ndim .eq. 3 )then
                    if (idecal+3 .gt. nval) then
                        call utmess('F', 'DISCRETS_21')
                    end if
                    do jj = 1, 3
                        vale(jj) = val(jj+idecal)
                    end do
                    idecal = idecal+3
                else
                    if (idecal+2 .gt. nval) then
                        call utmess('F', 'DISCRETS_21')
                    end if
                    do jj = 1, 2
                        vale(jj) = val(jj+idecal)
                    enddo
                    idecal = idecal+2
                endif
                call rairep(noma, ioc, car(nc), vale, ng, &
                            grplmax, zjdlm, nbparno, zk8(itbno), zr(irgno), rirot, ndim, jdme)
            else if (trarot) then
                lamass = 'M'//car(nc)(2:8)
                ! En 3D          1  2  3  4  5  6     1  2  3  4  5  6
                !   ncar=1  :   vx vy vz rx ry rz
                !   ncar=2  :   vx vy vz rx ry rz    vx vy vz rx ry rz
                ! En 2D
                !   ncar=1  :   vx vy          rz
                !   ncar=2  :   vx vy          rz    vx vy          rz
                if ( ndim .eq. 3 )then
                    if (idecal+6 .gt. nval) then
                        call utmess('F', 'DISCRETS_21')
                    end if
                    do jj = 1, 6
                        vale(jj) = val(jj+idecal)
                    end do
                    idecal = idecal+6
                else
                    if (idecal+3 .gt. nval) then
                        call utmess('F', 'DISCRETS_21')
                    end if
                    vale(1) = val(1+idecal)
                    vale(2) = val(2+idecal)
                    vale(6) = val(3+idecal)
                    idecal = idecal+3
                endif
                call rairep(noma, ioc, car(nc), vale, ng, &
                            grplmax, zjdlm, nbparno, zk8(itbno), zr(irgno), rirot, ndim, jdme)
            else
                ASSERT(.false.)
            end if
            zk8(itbmp:itbmp+nbparno-1) = ' '
!
            nb_noeu_disc = 0
            lokm = 0
            if (transl) lokm = 7
            if (trarot) lokm = 8
            if (car(nc)(lokm:lokm) .eq. 'N') nb_noeu_disc = 1
            if (car(nc)(lokm:lokm) .eq. 'L') nb_noeu_disc = 2
            ASSERT((nb_noeu_disc .gt. 0) .and. (lokm .gt. 0))
!
!           Le groupe de maille de DISCRET (POI1 ou SEG2)
            call jelira(jexnom(noma//'.GROUPEMA', nogp), 'LONMAX', nma)
            call jeveuo(jexnom(noma//'.GROUPEMA', nogp), 'L', ldgm)
!
            if (nma .ne. nbparno) then
                vali(1) = nbparno
                vali(2) = nma
                call utmess('F', 'MODELISA2_10', sk=nogp, ni=2, vali=vali)
            end if
!           Boucle sur les mailles du groupe 'nogp'
            do in = 0, nma-1
!               Le numéro de la maille
                nummail = zi(ldgm+in)
!               Vérification que la maille DISCRET fait partie du maillage
                if (nummail .le. 0) then
                    nommai = '????'
                    call utmess('F', 'AFFECARAELEM_25', si=ioc, nk=2, valk=[nogp,nommai])
                endif
!               Nom de la maille du discret
                nommai = int_to_char8(nummail)
!               Vérification que la maille DISCRET fait partie du modèle
                if ( zjdlm(nummail) .eq. 0 ) then
                    call utmess('F', 'AFFECARAELEM_25', si=ioc, nk=2, valk=[nogp,nommai])
                endif
                zjdlm(nummail) = -abs(zjdlm(nummail))
!               Récupère le nombre de noeud de la maille
                call jelira(jexnum(noma//'.CONNEX', nummail), 'LONMAX', nbnma)
                call jeveuo(jexnum(noma//'.CONNEX', nummail), 'L', ldnm)
!               Si la maille n'a pas le bon nombre de noeud
                if (nbnma .ne. nb_noeu_disc) then
                    call utmess('F', 'MODELISA_20', sk=nommai)
                end if
!               Boucle sur le nb de noeud de la maille
                do inbn = 1, nbnma
                    inoe = zi(ldnm+inbn-1)
!                   Nom du noeud
                    nomnoe = int_to_char8(inoe)
!                   Vérification que le noeud fait partie de la surface (sortie de rairep)
                    do ino = 1, nbparno
                        if (zk8(itbno+ino-1) .eq. nomnoe) then
                            zk8(itbmp+ino-1) = nommai
                            goto 22
                        end if
                    end do
                end do
!               Si on passe ici aucun des noeuds du discret appartient à la surface.
!               Ce n'est pas normal ==> F
                write (ifm, *) 'GROUP_MA :', (' '//grplmax(ii), ii=1, ng)
                call utmess('F', 'MODELISA_21', sk=nomnoe)
22              continue
            end do
!
!           Vérification que les discrets sont fixés au radier
            if (nc .eq. 1) then
                do ino = 1, nbparno
                    if (zk8(itbmp+ino-1) .eq. ' ') then
                        nomnoe = int_to_char8(ino)
                        call utmess('F', 'MODELISA2_8', sk=nomnoe)
                    end if
                end do
            end if
!
!           Préparation des impressions dans le fichier message
            lorep = 5
            if (irep .eq. 1) lorep = 6
            if (okunite) then
                if (transl) then
                    write(iunite, 105) car(nc)(1:lokm)
                else
                    write(iunite, 106) car(nc)(1:lokm), rirot(1), rirot(2), rirot(3)
                end if
            end if
!
            if (okunite) then
                do ii = 1, nbparno
                    iv = 1
                    jd = itbmp+ii-1
                    jn = itbno+ii-1
                    if (nb_noeu_disc .eq. 1) then
                        if (transl) then
                            write(iunite, 110) 'NOEUD', zk8(jn), car(nc)(1:lokm), &
                                (zr(irgno+6*(ii-1)+jj), jj=0, 2), repdis(irep)(1:lorep)
                        else
                            write(iunite, 111) 'NOEUD', zk8(jn), car(nc)(1:lokm), &
                                (zr(irgno+6*(ii-1)+jj), jj=0, 5), repdis(irep)(1:lorep)
                        end if
                    else
                        if (transl) then
                            write(iunite, 110) 'MAILLE', zk8(jd), car(nc)(1:lokm), &
                                (zr(irgno+6*(ii-1)+jj), jj=0, 2), repdis(irep)(1:lorep)
                        else
                            write(iunite, 111) 'MAILLE', zk8(jd), car(nc)(1:lokm), &
                                (zr(irgno+6*(ii-1)+jj), jj=0, 5), repdis(irep)(1:lorep)
                        end if
                    end if
                end do
            end if
!
            do ii = 1, nbparno
                iv = 1
                jd = itbmp+ii-1
                jn = itbno+ii-1
!
!               Affectation des valeurs réparties
                call affdis(ndim, irep, eta, car(nc), zr(irgno+6*ii-6), &
                            jdc, jdv, infdonn%ivr, iv, kma, &
                            ncmp, ll, jdcinf, jdvinf, isym)
                call nocart(cartdi, 3, dimcar, mode='NOM', nma=1, limano=[zk8(jd)])
                call nocart(cart(ll), 3, ncmp, mode='NOM', nma=1, limano=[zk8(jd)])
!               affectation de matrice masse nulle
                iv = 1
                call r8inir(nbval, 0.0d0, vale, 1)
                call affdis(ndim, irep, eta, lamass, vale, &
                            jdc, jdv, infdonn%ivr, iv, kma, &
                            ncmp, ll, jdcinf, jdvinf, isym)
                call nocart(cartdi, 3, dimcar, mode='NOM', nma=1, limano=[zk8(jd)])
                call nocart(cart(ll), 3, ncmp, mode='NOM', nma=1, limano=[zk8(jd)])
            end do
        end do
        if (idecal .ne. nval) then
            call utmess('F', 'DISCRETS_21')
        end if
30      continue
    end do
!
    call jedetr('&&TMPTABNO')
    call jedetr('&&TMPRIGNO')
    call jedetr('&&TMPTABMP')
!
999 continue
    call jedema()
!
100 format(/, ' <DISCRET> MATRICES AFFECTEES AUX ELEMENTS DISCRET ', &
            '(REPERE ', a6, '), OCCURRENCE ', i4)
105 format(/, ' PAS DE REPARTITION EN ROTATION POUR DES ', a,/)
106 format(/, ' RAIDEURS DE ROTATION A REPARTIR POUR DES ', a, / &
            , '  RX: ', 1pe12.5, ' RY: ', 1pe12.5, ' RZ: ', 1pe12.5,/)
110 format(' _F(', a, '=''', a8, ''', CARA=''', a, ''',', /, &
           '    VALE=(', 3(1x, 1pe12.5, ','), '),', /, &
           '    REPERE=''', a, '''),')
111 format(' _F(', a, '=''', a8, ''', CARA=''', a, ''',', /, &
           '    VALE=(', 3(1x, 1pe12.5, ','), /, &
           '          ', 3(1x, 1pe12.5, ','), '),', /, &
           '    REPERE=''', a, '''),')
! 210 format(a,2x,a8,2x,a,3(2x,1pe12.5),2x,a)
! 211 format(a,2x,a8,2x,a,3(2x,1pe12.5),2(1x,1pe12.5),2x,a)
!&>
end subroutine
