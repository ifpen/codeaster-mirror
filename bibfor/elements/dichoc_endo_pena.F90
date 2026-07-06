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
subroutine dichoc_endo_pena(for_discret, iret)
!
! -------------------------------------------------------------------------------
!
!     RELATION DE COMPORTEMENT "CHOC_ENDO_PENA" : COMPORTEMENT DISCRET CHOC NON-LINEAIRE
!
! --------------------------------------------------------------------------------------------------
!
! IN    for_discret : voir l'appel
! OUT   iret        : code retour
!
! --------------------------------------------------------------------------------------------------
!
    use te0047_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/diraidklv.h"
#include "asterfort/diklvraid.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvala.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/Behaviour_type.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! -------------------------------------------------------------------------------
!
    integer(kind=8) :: jdc, irep, imat, ivarim, ii, ivitp, idepen, iviten, neq, igeom, ivarip
    integer(kind=8) :: iretlc, ifono, imatsym
    integer(kind=8) :: icontp, iadzi, iazk24
!
    real(kind=8) :: dvl(for_discret%nno*for_discret%nc), dpe(for_discret%nno*for_discret%nc)
    real(kind=8) :: dve(for_discret%nno*for_discret%nc)
    real(kind=8) :: klv(for_discret%nbt), force(3), fl(for_discret%nno*for_discret%nc), raide(6)
    real(kind=8) :: r8bid
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: k8bid
    character(len=24) :: messak(6)
!
    aster_logical :: resi
! -------------------------------------------------------------------------------
    integer(kind=8), parameter :: nbres = 3
    real(kind=8) :: valres(nbres)
    integer(kind=8) :: codres(nbres)
    character(len=8) :: nomres(nbres)
    integer(kind=8) :: nbpar
    real(kind=8) :: valpar
    character(len=8) :: nompar
    integer(kind=8) :: tecro2
!
! -------------------------------------------------------------------------------
!   Variables internes
    integer(kind=8), parameter :: nbvari = 3
    real(kind=8) :: varmo(nbvari), varpl(nbvari)
! -------------------------------------------------------------------------------
    real(kind=8) :: xl(6), xd(3), rignor, amornor_in, deplace, ld, seuil, seuil2
    real(kind=8) :: enfoncement_resi, enfoncement_resi_moins, enfoncement
    real(kind=8) :: ftry, vitesse
    real(kind=8) :: indic_charge, indic_charge_moins, enfoncement_max
    blas_int :: b_incx, b_incy, b_n
! -------------------------------------------------------------------------------
!
    iret = 0
    resi = for_discret%lVect .or. for_discret%lSigm .or. for_discret%lVari
!   Seulement en 3D, sur un segment, avec seulement de la translation
    if ((for_discret%nomte(1:12) .ne. 'MECA_DIS_T_L') .or. (for_discret%ndim .ne. 3) .or. &
        (for_discret%nno .ne. 2) .or. (for_discret%nc .ne. 3)) then
        call jevech('PCOMPOR', 'L', vk16=compor)
        messak(1) = for_discret%nomte
        messak(2) = for_discret%option
        messak(3) = compor(INCRELAS)
        messak(4) = compor(RELA_NAME)
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_22', nk=5, valk=messak)
    end if
!   Nombre de degré de liberté
    neq = for_discret%nno*for_discret%nc
!   Paramètres en entrée
    call jevech('PCADISK', 'L', jdc)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imat)
!
    call infdis('REPK', irep, r8bid, k8bid)
!   irep = 1 = matrice en repère global ==> passer en local
    if (irep .eq. 1) then
        if (for_discret%ndim .eq. 3) then
            call utpsgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
        else if (for_discret%ndim .eq. 2) then
            call ut2mgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
        end if
    else
        b_n = to_blas_int(for_discret%nbt)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
    end if
!   Récupération des termes diagonaux : raide = klv(i,i)
    call diraidklv(for_discret%nomte, raide, klv)
!
!   Champ de vitesse
    call tecach('ONO', 'PVITPLU', 'L', iretlc, iad=ivitp)
    if (iretlc .eq. 0) then
        if (for_discret%ndim .eq. 3) then
            call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(ivitp), dvl)
        else if (for_discret%ndim .eq. 2) then
            call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(ivitp), dvl)
        end if
    else
        dvl(:) = 0.0d0
    end if
!   Champ de déplacement d'entrainement
    call tecach('ONO', 'PDEPENT', 'L', iretlc, iad=idepen)
    if (iretlc .eq. 0) then
        if (for_discret%ndim .eq. 3) then
            call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(idepen), dpe)
        else if (for_discret%ndim .eq. 2) then
            call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(idepen), dpe)
        end if
    else
        dpe(:) = 0.0d0
    end if
!   Champ de vitesse d'entrainement
    call tecach('ONO', 'PVITENT', 'L', iretlc, iad=iviten)
    if (iretlc .eq. 0) then
        if (for_discret%ndim .eq. 3) then
            call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(iviten), dve)
        else if (for_discret%ndim .eq. 2) then
            call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(iviten), dve)
        end if
    else
        dve(:) = 0.d0
    end if
!
!   Variables internes
    call jevech('PVARIMR', 'L', ivarim)
    do ii = 1, nbvari
        varmo(ii) = zr(ivarim+ii-1)
    end do
!
! -------------------------------------------------------------------------------
!   Relation de comportement de choc
!
!   Coordonnees du discret dans le repère local
    xl(:) = 0.0
    if (for_discret%ndim .eq. 3) then
        call utpvgl(for_discret%nno, 3, for_discret%pgl, zr(igeom), xl)
    else if (for_discret%ndim .eq. 2) then
        call ut2vgl(for_discret%nno, 2, for_discret%pgl, zr(igeom), xl)
    end if
!
!   Caractéristiques du matériau
    nbpar = 0; nompar = ' '; valpar = 0.d0
    valres(:) = 0.0; nomres(:) = ' '
    nomres(1) = 'DIST_1'; nomres(2) = 'DIST_2'
    call rcvala(zi(imat), ' ', 'DIS_CHOC_ENDO', nbpar, nompar, &
                [valpar], 2, nomres, valres, codres, &
                0, nan='NON')
!
!   calcul du jeu
    if (for_discret%nno .eq. 2) then
!       longueur du discret
        xd(1:3) = xl(1+for_discret%ndim:2*for_discret%ndim)-xl(1:for_discret%ndim)
        ld = xd(1)-valres(1)-valres(2)
    else
        ld = valres(3)-valres(1)
    end if
!
!   recuperation variables internes t(-)
    enfoncement_max = varmo(1)
    enfoncement_resi_moins = varmo(2)
    indic_charge_moins = varmo(3)
!
!   calcul de l'enfoncement
    if (resi) then
        if (for_discret%nno .eq. 1) then
            deplace = for_discret%ulm(1)+for_discret%dul(1)+dpe(1)
            vitesse = dvl(1)+dve(1)
        else
            deplace = (for_discret%ulm(1+for_discret%nc)-for_discret%ulm(1)+ &
                       for_discret%dul(1+for_discret%nc)-for_discret%dul(1)+ &
                       dpe(1+for_discret%nc)-dpe(1))
            vitesse = (dvl(1+for_discret%nc)-dvl(1))+ &
                      (dve(1+for_discret%nc)-dve(1))
        end if
    else
        if (for_discret%nno .eq. 1) then
            deplace = for_discret%ulm(1)+dpe(1)
            vitesse = dvl(1)+dve(1)
        else
            deplace = (for_discret%ulm(1+for_discret%nc)-for_discret%ulm(1))+ &
                      (dpe(1+for_discret%nc)-dpe(1))
            vitesse = (dvl(1+for_discret%nc)-dvl(1))+ &
                      (dve(1+for_discret%nc)-dve(1))
        end if
    end if
    enfoncement = deplace+ld
    if (enfoncement < enfoncement_max) then
        enfoncement_max = enfoncement
    end if
!
!   Type d'amortissement inclus ou exclus
    tecro2 = 0
    call rcvala(zi(imat), ' ', 'DIS_CHOC_ENDO', 0, ' ', &
                [0.0d0], 1, ['CRIT_AMOR'], valres, codres, &
                1)
    tecro2 = nint(valres(1))
!   récupération de l'effort enveloppe, de la raideur et de l'amortissement
    nomres(1) = 'FX'; nomres(2) = 'RIGI_NOR'; nomres(3) = 'AMOR_NOR'
    call rcvala(zi(imat), ' ', 'DIS_CHOC_ENDO', 1, 'DX', &
                [-enfoncement], nbres, nomres, valres, codres, &
                1)
    seuil = valres(1)
!
    call rcvala(zi(imat), ' ', 'DIS_CHOC_ENDO', 1, 'DX', &
                [-enfoncement_max], nbres, nomres, valres, codres, &
                1)
    seuil2 = valres(1)
    rignor = valres(2)
    amornor_in = 0.0
    if (tecro2 .eq. 1) then
        amornor_in = valres(3)
    end if
!
!   indic_charge [0, 1, 2] : [pas de contact, contact élastique, sur le seuil]
    force(:) = 0.0
    if (enfoncement-enfoncement_resi_moins > 0) then
        force(1) = 0.0
        indic_charge = 0.0
        enfoncement_resi = enfoncement_resi_moins
    else
! correction de l'enfoncement residuel (vitesse = 0)
        ftry = rignor*(enfoncement-enfoncement_resi_moins)
        if (abs(ftry) .gt. seuil2) then
            enfoncement_resi = enfoncement+seuil2/rignor
        else
            enfoncement_resi = enfoncement_resi_moins
        end if
! limitation de f au seuil (en prenant en compte l'amortissement)
        ftry = rignor*(enfoncement-enfoncement_resi)+amornor_in*vitesse
        if (abs(ftry) .ge. seuil) then
            force(1) = -seuil
            indic_charge = 2.0
            raide(1) = rignor
        else if (ftry .gt. 0.d0) then
            force(1) = 0.0
            indic_charge = 0.0
        else
            force(1) = ftry
            indic_charge = 1.0
            raide(1) = rignor
        end if
    end if
! Actualisation de la matrice de raideur
    call diklvraid(for_discret%nomte, klv, raide)
    !
! stockage des contraintes
    if (for_discret%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        if (for_discret%nno .eq. 1) then
            zr(icontp-1+1) = force(1)
            zr(icontp-1+2) = force(2)
            if (for_discret%ndim .eq. 3) then
                zr(icontp-1+3) = force(3)
            end if
        else if (for_discret%nno .eq. 2) then
            zr(icontp-1+1) = force(1)
            zr(icontp-1+1+for_discret%nc) = force(1)
            zr(icontp-1+2) = force(2)
            zr(icontp-1+2+for_discret%nc) = force(2)
            if (for_discret%ndim .eq. 3) then
                zr(icontp-1+3) = force(3)
                zr(icontp-1+3+for_discret%nc) = force(3)
            end if
        end if
    end if
! stockage des efforts
    if (for_discret%lVect) then
        fl(:) = 0.d0
        if (for_discret%nno .eq. 1) then
            fl(1) = force(1)
            fl(2) = force(2)
            if (for_discret%ndim .eq. 3) then
                fl(3) = force(3)
            end if
        else if (for_discret%nno .eq. 2) then
            fl(1) = -force(1)
            fl(1+for_discret%nc) = force(1)
            fl(2) = -force(2)
            fl(2+for_discret%nc) = force(2)
            if (for_discret%ndim .eq. 3) then
                fl(3) = -force(3)
                fl(3+for_discret%nc) = force(3)
            end if
        end if
        call jevech('PVECTUR', 'E', ifono)
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        end if
    end if
! stockage variables internes
    if (for_discret%lVari) then
        varpl(1) = enfoncement_max
        varpl(2) = enfoncement_resi
        varpl(3) = indic_charge
        call jevech('PVARIPR', 'E', ivarip)
        if (for_discret%nno .eq. 1) then
            do ii = 1, nbvari
                zr(ivarip+ii-1) = varpl(ii)
            end do
        else
            do ii = 1, nbvari
                zr(ivarip+ii-1) = varpl(ii)
                zr(ivarip+ii-1+nbvari) = varpl(ii)
            end do
        end if
    end if
! stockage matrice tangente
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imatsym)
        if (for_discret%ndim .eq. 3) then
            call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imatsym))
        else if (for_discret%ndim .eq. 2) then
            call ut2mlg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imatsym))
        end if
    end if
    !
end subroutine
