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

subroutine te0041(option, nomte)
    implicit none
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!
!        CALCUL DES MATRICES DE RAIDEUR, MASSE, AMORTISSEMENT
!                   POUR  LES ELEMENTS DISCRETS
!
! --------------------------------------------------------------------------------------------------
!
!     option : nom de l'option a calculer
!        pour discrets : symétriques  et non-symetriques
!           RIGI_MECA   MASS_MECA  MASS_MECA_DIAG  AMOR_MECA
!        pour discrets : symétriques
!           RIGI_MECA_HYST   RIGI_MECA_TANG   RIGI_FLUI_STRU
!           M_GAMMA          MASS_FLUI_STRU   MASS_MECA_EXPLI
!
!     nomte  : nom du type d'élément
!           MECA_DIS_T_N      MECA_DIS_T_L
!           MECA_DIS_TR_N     MECA_DIS_TR_L
!           MECA_2D_DIS_T_N   MECA_2D_DIS_T_L
!           MECA_2D_DIS_TR_N  MECA_2D_DIS_TR_L
!
! --------------------------------------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/dikpkt.h"
#include "asterfort/infdis.h"
#include "asterfort/infted.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/matrot.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvala.h"
#include "asterfort/tecach.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2plg.h"
#include "asterfort/utmess.h"
#include "asterfort/utpplg.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpngl.h"
#include "asterfort/utpslg.h"
#include "asterfort/vecma.h"
#include "asterfort/Behaviour_type.h"
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8)         ::  nddlm, nl1, nl2, infodi, ntermx
    parameter(nddlm=12, nl1=nddlm*(nddlm+1)/2, nl2=nddlm*nddlm, ntermx=144)
!
    real(kind=8)    :: mata1(nl1), mata2(nl1), mata3(nl2), mata4(nl2)
!
    integer(kind=8)         :: ibid, itype, irep, nbterm, nno, nc, ndim, nddl, i, j, iret
    integer(kind=8)         :: jdr, jdm, lorien, jdc, iacce, ivect, jma, jVNonLin, jvarip
    real(kind=8)    :: pgl(3, 3), matv1(nl1), matp(nddlm, nddlm)
    real(kind=8)    :: eta, r8bid, xrota, enfoncement_max, riginor_choc, amornor_choc, amortan_choc
    real(kind=8)    :: tempo(ntermx)
    complex(kind=8) :: hyst, dcmplx
    aster_logical :: lNonLinear, lDisChoc, lDisChocEndo, IsCoulomb, IsCoin2D
!
    character(len=8)    :: k8bid
    character(len=24)   :: valk(2)
!
    integer(kind=8)             :: icodre(5), icodre2(2), icodre3(1), tecro2
    real(kind=8)        :: valres(5), valres2(2), valres3(1)
    character(len=16)   :: nomres(5), nomres2(2)
    character(len=16), pointer :: compor(:) => null()
    real(kind=8)        :: kp, kt1, kt2, indicChoc
! --------------------------------------------------------------------------------------------------
    aster_logical       :: assemble_amor
! --------------------------------------------------------------------------------------------------
!   Ce sont bien des éléments discrets
    ASSERT(lteatt('DIM_TOPO_MODELI', '-1'))
!
!   On vérifie que les caractéristiques ont été affectées
!   Le code du discret
    call infdis('CODE', ibid, r8bid, nomte)
!   Le code stoke dans la carte
    call infdis('TYDI', infodi, r8bid, k8bid)
    if (infodi .ne. ibid) then
        call utmess('F+', 'DISCRETS_25', sk=nomte)
        call infdis('DUMP', ibid, r8bid, 'F+')
    end if
!
    assemble_amor = ASTER_FALSE
    lDisChoc = ASTER_FALSE
    lDisChocEndo = ASTER_FALSE
    if (option .eq. 'AMOR_MECA') then
        call tecach('ONO', 'PNONLIN', 'L', iret, iad=jVNonLin)
        lNonLinear = ASTER_FALSE
        if (iret .eq. 0) then
            call jevech('PNONLIN', 'L', jVNonLin)
            lNonLinear = zi(jVNonLin) .eq. 1
        end if
        if (lNonLinear) then
! --------- Nonlinear cases (=> only for DIS_CHOC/CHOC_ENDO_PENA)
            call jevech('PCOMPOR', 'L', vk16=compor)
            lDisChoc = compor(RELA_NAME) .eq. 'DIS_CHOC'
            lDisChocEndo = compor(RELA_NAME) .eq. 'CHOC_ENDO_PENA'
            assemble_amor = lDisChoc .or. lDisChocEndo
        else
! --------- Linear case (=> for all cases)
            assemble_amor = ASTER_TRUE
        end if
    end if
!
    if (option .eq. 'RIGI_MECA') then
        call infdis('SYMK', infodi, r8bid, k8bid)
    else if (option .eq. 'MASS_MECA') then
        call infdis('SYMM', infodi, r8bid, k8bid)
    else if (option .eq. 'MASS_MECA_DIAG') then
        call infdis('SYMM', infodi, r8bid, k8bid)
    else if (option .eq. 'AMOR_MECA') then
        call infdis('SYMA', infodi, r8bid, k8bid)
    else
!       Pour les autres options toutes les matrices doivent être symétriques
        call infdis('SKMA', infodi, r8bid, k8bid)
        if (infodi .ne. 3) then
            valk(1) = option
            valk(2) = nomte
            call utmess('F', 'DISCRETS_32', nk=2, valk=valk)
        end if
        ! Elles sont toutes symétriques
        infodi = 1
    end if
!
!   Informations sur les discrets :
!       nbterm   = nombre de coefficients dans K
!       nno      = nombre de noeuds
!       nc       = nombre de composante par noeud
!       ndim     = dimension de l'élément
!       itype    = type de l'élément
    call infted(nomte, infodi, nbterm, nno, nc, ndim, itype)
!   Nombre de DDL par noeuds
    nddl = nno*nc
!
    if (infodi .eq. 1 .and. option .eq. 'RIGI_MECA_HYST') then
        call jevech('PRIGIEL', 'L', jdr)
        call jevech('PMATUUC', 'E', jdm)
        call infdis('ETAK', ibid, eta, k8bid)
        hyst = dcmplx(1.0, eta)
        do i = 1, nbterm
            zc(jdm+i-1) = zr(jdr+i-1)*hyst
        end do
        goto 999
    end if
!
!   Matrice de passage global vers local
    call jevech('PCAORIE', 'L', lorien)
    call matrot(zr(lorien), pgl)
    xrota = abs(zr(lorien))+abs(zr(lorien+1))+abs(zr(lorien+2))
!
!   Matrices symétriques
    if (infodi .eq. 1) then
        matv1 = 0.0
        mata1 = 0.0
        mata2 = 0.0
!
        if ((option .eq. 'RIGI_MECA') .or. &
            (option .eq. 'RIGI_MECA_TANG') .or. &
            (option .eq. 'RIGI_FLUI_STRU')) then
!           Discret de type raideur
            call infdis('DISK', infodi, r8bid, k8bid)
            if (infodi .eq. 0) then
                call utmess('A+', 'DISCRETS_27', sk=nomte)
                call infdis('DUMP', ibid, r8bid, 'A+')
            end if
            call jevech('PCADISK', 'L', jdc)
            call infdis('REPK', irep, r8bid, k8bid)
            call jevech('PMATUUR', 'E', jdm)
        else if ((option .eq. 'MASS_MECA') .or. &
                 (option .eq. 'MASS_MECA_DIAG') .or. &
                 (option .eq. 'MASS_MECA_EXPLI') .or. &
                 (option .eq. 'MASS_FLUI_STRU')) then
!           Discret de type masse
            call infdis('DISM', infodi, r8bid, k8bid)
            if (infodi .eq. 0) then
                call utmess('A+', 'DISCRETS_26', sk=nomte)
                call infdis('DUMP', ibid, r8bid, 'A+')
            end if
            call jevech('PCADISM', 'L', jdc)
            call infdis('REPM', irep, r8bid, k8bid)
            call jevech('PMATUUR', 'E', jdm)
        else if (option .eq. 'AMOR_MECA') then
!           Discret de type amortissement
            call infdis('DISA', infodi, r8bid, k8bid)
            if (infodi .eq. 0) then
                call utmess('A+', 'DISCRETS_28', sk=nomte)
                call infdis('DUMP', ibid, r8bid, 'A+')
            end if
            call jevech('PCADISA', 'L', jdc)
            call infdis('REPA', irep, r8bid, k8bid)
            call jevech('PMATUUR', 'E', jdm)
!           Traitement du cas de assemble_amor
            if (assemble_amor) then
                if (ndim .ne. 3) goto 666
                call tecach('NNN', 'PMATERC', 'L', iret, iad=jma)
                if ((jma .eq. 0) .or. (iret .ne. 0)) goto 666
                if (lDisChoc) then
                    tecro2 = 0
                    ! Récupération des paramètres matériau DIS_CONTACT
                    nomres(1) = 'RIGI_NOR'
                    nomres(2) = 'AMOR_NOR'
                    nomres(3) = 'AMOR_TAN'
                    nomres(4) = 'COULOMB'
                    nomres(5) = 'CONTACT'
                    valres = 0.0
                    call rcvala(zi(jma), ' ', 'DIS_CONTACT', 0, ' ', &
                                [0.0d0], 5, nomres, valres, icodre, 0)
                    riginor_choc = valres(1)
                    amornor_choc = valres(2)
                    amortan_choc = valres(3)
                    ! --- Vérification du frottement de Coulomb
                    IsCoulomb = ASTER_FALSE
                    if (icodre(4) .eq. 0) then
                        if (valres(4) .gt. r8prem()) then
                            IsCoulomb = ASTER_TRUE
                        end if
                    end if
                    ! --- Vérification du type de contact (1D ou COIN_2D)
                    IsCoin2D = ASTER_FALSE
                    if (icodre(5) .eq. 0) then
                        if (nint(valres(5)) .ne. 0) then
                            IsCoin2D = ASTER_TRUE
                        end if
                    end if

                    ! Récupération de la matrice tangente
                    if ((lDisChoc) .and. (IsCoulomb) .and. (.not. IsCoin2D)) then
                        ! Cas de DIS_CHOC avec matrice tangente non symétrique
                        call tecach('ONO', 'PRIGINS', 'L', iret, iad=jdr)
                        if (jdr .eq. 0) goto 666
                        call utpngl(nno, nc, pgl, zr(jdr), matv1)
                    else
                        ! Cas d'une matrice tangente symétrique
                        call tecach('ONO', 'PRIGIEL', 'L', iret, iad=jdr)
                        if (jdr .eq. 0) goto 666
                        call utpsgl(nno, nc, pgl, zr(jdr), matv1)
                    end if

                    ! Récupération des raideurs élastiques en parallèle
                    call dikpkt(zi(jma), 'DIS_CONTACT', kp, kt1, kt2)

                else if (lDisChocEndo) then
                    ! Récupération des paramètres matériau DIS_CHOC_ENDO
                    call jevech('PVARIPG', 'L', jvarip)
                    enfoncement_max = zr(jvarip)
                    nomres2(1) = 'RIGI_NOR'; nomres2(2) = 'AMOR_NOR'
                    call rcvala(zi(jma), ' ', 'DIS_CHOC_ENDO', 1, 'DX', &
                                [-enfoncement_max], 2, nomres2, valres2, icodre2, &
                                1)
                    icodre(1) = icodre2(1)
                    icodre(2) = icodre2(2)
                    riginor_choc = valres2(1)
                    amornor_choc = valres2(2)

                    !   Type d'amortissement inclus ou exclus
                    tecro2 = 0
                    call rcvala(zi(jma), ' ', 'DIS_CHOC_ENDO', 0, ' ', &
                                [0.0d0], 1, ['CRIT_AMOR'], valres3, icodre3, &
                                1)
                    tecro2 = nint(valres3(1))

                    ! Récupération de la matrice tangente (symétrique)
                    call tecach('ONO', 'PRIGIEL', 'L', iret, iad=jdr)
                    if (jdr .eq. 0) goto 666
                    call utpsgl(nno, nc, pgl, zr(jdr), matv1)

                    ! Récupération des raideurs élastiques en parallèle
                    call dikpkt(zi(jma), 'DIS_CHOC_ENDO', kp, kt1, kt2)

                end if

                ! Prise en compte de l'amortissement de choc
                if (icodre(1) .eq. 0) then
                    if (abs(riginor_choc) > r8prem()) then
                        ! Définition du facteur "indicateur de choc" (1 si choc, 0 sinon)
                        ! (contribution élastique à retrancher à la matrice tangente)
                        indicChoc = (matv1(1)-kp)/riginor_choc
                        if ((icodre(2) .eq. 0) .and. (tecro2 .ne. 1)) then
                            mata1(1) = indicChoc*amornor_choc
                        end if
                        if (lDisChoc) then
                            if (icodre(3) .eq. 0) then
                                mata1(3) = indicChoc*amortan_choc
                            end if
                            mata1(6) = mata1(3)
                        end if
                    end if
                end if
                if (nno .eq. 2 .and. nc .eq. 3) then
                    mata1(7) = -mata1(1)
                    mata1(10) = mata1(1)
                    mata1(15) = mata1(3)
                    mata1(21) = mata1(3)
                    mata1(12) = -mata1(3)
                    mata1(18) = -mata1(3)
                else if (nno .eq. 2 .and. nc .eq. 6) then
                    mata1(22) = -mata1(1)
                    mata1(28) = mata1(1)
                    mata1(36) = mata1(3)
                    mata1(45) = mata1(3)
                    mata1(30) = -mata1(3)
                    mata1(39) = -mata1(3)
                end if
                call utpslg(nno, nc, pgl, mata1, mata2)
666             continue
            end if
        else if (option .eq. 'M_GAMMA') then
!           Discret de type masse
            call infdis('DISM', infodi, r8bid, k8bid)
            if (infodi .eq. 0) then
                call utmess('A+', 'DISCRETS_26', sk=nomte)
                call infdis('DUMP', ibid, r8bid, 'A+')
            end if
            call jevech('PCADISM', 'L', jdc)
            call infdis('REPM', irep, r8bid, k8bid)
            call jevech('PACCELR', 'L', iacce)
            call jevech('PVECTUR', 'E', ivect)
        else
!           Option de calcul invalide
            ASSERT(.false.)
        end if
!
        if (irep .eq. 1) then
!           Repère global ==> pas de rotation
            if (option .eq. 'M_GAMMA') then
                call vecma(zr(jdc), nbterm, matp, nddl)
                call pmavec('ZERO', nddl, matp, zr(iacce), zr(ivect))
            else
                do i = 1, nbterm
                    zr(jdm+i-1) = zr(jdc+i-1)+mata2(i)
                end do
            end if
        else if (irep .eq. 2) then
!           Local ==> Global
            if (xrota <= r8prem()) then
!               Angles quasi nuls  ===>  pas de rotation
                if (option .eq. 'M_GAMMA') then
                    call vecma(zr(jdc), nbterm, matp, nddl)
                    call pmavec('ZERO', nddl, matp, zr(iacce), zr(ivect))
                else
                    do i = 1, nbterm
                        zr(jdm+i-1) = zr(jdc+i-1)+mata2(i)
                    end do
                end if
            else
!               Angles non nuls  ===>  rotation
                if (option .eq. 'M_GAMMA') then
                    if (ndim .eq. 3) then
                        call utpslg(nno, nc, pgl, zr(jdc), matv1)
                    else if (ndim .eq. 2) then
                        call ut2mlg(nno, nc, pgl, zr(jdc), matv1)
                    end if
                    call vecma(matv1, nbterm, matp, nddl)
                    call pmavec('ZERO', nddl, matp, zr(iacce), zr(ivect))
                else
                    if (ndim .eq. 3) then
                        call utpslg(nno, nc, pgl, zr(jdc), zr(jdm))
                        if (assemble_amor) then
                            do i = 1, nbterm
                                zr(jdm+i-1) = zr(jdm+i-1)+mata2(i)
                            end do
                        end if
                    else if (ndim .eq. 2) then
                        call ut2mlg(nno, nc, pgl, zr(jdc), zr(jdm))
                    end if
                end if
            end if
        end if

!   Matrices non-symétriques
    else
        matv1 = 0.0
        mata3 = 0.0
        mata4 = 0.0
        tempo = 0.0
!
        if (option .eq. 'RIGI_MECA') then
!           Discret de type raideur
            call infdis('DISK', infodi, r8bid, k8bid)
            if (infodi .eq. 0) then
                call utmess('A+', 'DISCRETS_27', sk=nomte)
                call infdis('DUMP', ibid, r8bid, 'A+')
            end if
            call jevech('PCADISK', 'L', jdc)
            call infdis('REPK', irep, r8bid, k8bid)
            call jevech('PMATUNS', 'E', jdm)
        else if ((option .eq. 'MASS_MECA') .or. &
                 (option .eq. 'MASS_MECA_DIAG')) then
!           Discret de type masse
            call infdis('DISM', infodi, r8bid, k8bid)
            if (infodi .eq. 0) then
                call utmess('A+', 'DISCRETS_26', sk=nomte)
                call infdis('DUMP', ibid, r8bid, 'A+')
            end if
            call jevech('PCADISM', 'L', jdc)
            call infdis('REPM', irep, r8bid, k8bid)
            call jevech('PMATUNS', 'E', jdm)
        else if (option .eq. 'AMOR_MECA') then
!           Discret de type amortissement
            call infdis('DISA', infodi, r8bid, k8bid)
            if (infodi .eq. 0) then
                call utmess('A+', 'DISCRETS_28', sk=nomte)
                call infdis('DUMP', ibid, r8bid, 'A+')
            end if
            call jevech('PCADISA', 'L', jdc)
            call infdis('REPA', irep, r8bid, k8bid)
            call jevech('PMATUNS', 'E', jdm)
            if (assemble_amor) then
                if (ndim .ne. 3) goto 777
                call tecach('NNN', 'PMATERC', 'L', iret, iad=jma)
                if ((jma .eq. 0) .or. (iret .ne. 0)) goto 777
                ! Récupération des paramètres matériau DIS_CONTACT
                nomres(1) = 'RIGI_NOR'
                nomres(2) = 'AMOR_NOR'
                nomres(3) = 'AMOR_TAN'
                nomres(4) = 'COULOMB'
                nomres(5) = 'CONTACT'
                valres = 0.0
                call rcvala(zi(jma), ' ', 'DIS_CONTACT', 0, ' ', &
                            [0.0d0], 5, nomres, valres, icodre, 0)
                ! --- Vérification du frottement de Coulomb
                IsCoulomb = ASTER_FALSE
                if (icodre(4) .eq. 0) then
                    if (valres(4) .gt. r8prem()) then
                        IsCoulomb = ASTER_TRUE
                    end if
                end if
                ! --- Vérification du type de contact (1D ou COIN_2D)
                IsCoin2D = ASTER_FALSE
                if (icodre(5) .eq. 0) then
                    if (nint(valres(5)) .ne. 0) then
                        IsCoin2D = ASTER_TRUE
                    end if
                end if

                ! Récupération de la matrice tangente
                if ((lDisChoc) .and. (IsCoulomb) .and. (.not. IsCoin2D)) then
                    ! Cas de DIS_CHOC avec matrice tangente non symétrique
                    call tecach('ONO', 'PRIGINS', 'L', iret, iad=jdr)
                    if (jdr .eq. 0) goto 777
                    call utpngl(nno, nc, pgl, zr(jdr), matv1)
                else
                    ! Cas d'une matrice tangente symétrique
                    call tecach('ONO', 'PRIGIEL', 'L', iret, iad=jdr)
                    if (jdr .eq. 0) goto 777
                    call utpsgl(nno, nc, pgl, zr(jdr), matv1)
                end if

                if (icodre(1) .eq. 0) then
                    if (abs(valres(1)) > r8prem()) then
                        if (icodre(2) .eq. 0) then
                            mata3(1) = matv1(1)*valres(2)/valres(1)
                        end if
                        if (icodre(3) .eq. 0) then
                            mata3(3) = matv1(1)*valres(3)/valres(1)
                        end if
                    end if
                end if
                if (nno .eq. 2 .and. nc .eq. 3) then
                    mata3(19) = -mata3(1)
                    mata3(22) = mata3(1)
                    mata3(29) = mata3(3)
                    mata3(36) = mata3(3)
                    mata3(26) = -mata3(3)
                    mata3(33) = -mata3(3)
                    mata3(8) = mata3(3)
                    mata3(15) = mata3(3)
                    mata3(4) = mata3(19)
                    mata3(11) = mata3(26)
                    mata3(18) = mata3(33)
                    mata3(3) = 0.d0
                else if (nno .eq. 2 .and. nc .eq. 6) then
                    mata3(73) = -mata3(1)
                    mata3(79) = mata3(1)
                    mata3(92) = mata3(3)
                    mata3(105) = mata3(3)
                    mata3(86) = -mata3(3)
                    mata3(99) = -mata3(3)
                    mata3(7) = mata3(73)
                    mata3(20) = mata3(86)
                    mata3(33) = mata3(99)
                    mata3(14) = mata3(3)
                    mata3(27) = mata3(3)
                    mata3(3) = 0.d0
                end if
                call utpplg(nno, nc, pgl, mata3, mata4)
777             continue
            end if
        else
!           Option de calcul invalide
            ASSERT(.false.)
        end if
!
        if (irep .eq. 1) then
!           Repère global ==> pas de rotation
            do i = 1, nbterm
                tempo(i) = zr(jdc+i-1)+mata4(i)
            end do
        else if (irep .eq. 2) then
!           Local ==> global
            if (xrota <= r8prem()) then
!               Angles quasi nuls  ===>  pas de rotation
                do i = 1, nbterm
                    tempo(i) = zr(jdc+i-1)+mata4(i)
                end do
            else
!               Angles non nuls  ===>  rotation
                if (ndim .eq. 3) then
                    call utpplg(nno, nc, pgl, zr(jdc), tempo)
                    if (assemble_amor) then
                        do i = 1, nbterm
                            tempo(i) = tempo(i)+mata4(i)
                        end do
                    end if
                else if (ndim .eq. 2) then
                    call ut2plg(nno, nc, pgl, zr(jdc), tempo)
                end if
            end if
        end if
!
        do i = 1, nddl
            do j = 1, nddl
                zr(jdm+(i-1)*nddl+j-1) = tempo((j-1)*nddl+i)
            end do
        end do
    end if
!
999 continue
!
end subroutine
