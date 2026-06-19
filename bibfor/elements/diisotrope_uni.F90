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

subroutine diisotrope_uni(for_discret, iret, ipi, jmat, ivarim, icontm, klv, raide)
! --- Extraction de diisotrope à partir de la ligne 172 ---
!     Cette procédure est indépendante : toutes les variables nécessaires
!     sont passées en arguments et les constantes locales sont redéclarées.
    use te0047_type
    implicit none

#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/diraidklv.h"
#include "asterfort/diklvraid.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvala.h"
#include "asterfort/rk5adp.h"
#include "asterfort/tecael.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "asterfort/disc_isotr.h"
#include "blas/dcopy.h"

    type(te0047_dscr), intent(in)    :: for_discret
    integer(kind=8), intent(inout) :: iret
    integer(kind=8), intent(in)    :: ipi, jmat, ivarim, icontm
    real(kind=8), intent(inout) :: klv(78), raide(6)

! --- paramètres / tableaux locaux ---
    integer(kind=8), parameter :: nbpara = 2, nbfct = 1*5
    integer(kind=8), parameter :: nbequa = 14
    real(kind=8), parameter :: precis = 1.0d-08
    integer(kind=8), parameter :: lmat = 9, lfct = 10
    integer(kind=8) :: ldcfct(nbfct)
    real(kind=8) :: ldcpar(nbpara)
    character(len=8) :: ldccar(1)
    real(kind=8) :: y0(nbequa), dy0(nbequa), resu(nbequa*2), ynorme(nbequa)
    real(kind=8) :: fl(12), klc(144)

    integer(kind=8) :: idf, iloi, nbvale, jprol, jvale, tecro
    integer(kind=8) :: jtp, jtm, icarcr, nbdecp
    integer(kind=8) :: imat, icontp, ivarip, ifono
    integer(kind=8) :: ii, kk, neq
    real(kind=8) :: temps0, temps1, dtemps, errmax

    integer(kind=8) :: codret(1)
    real(kind=8) :: valres(1)

    neq = for_discret%nno*for_discret%nc

    ldcfct(:) = -1
    iloi = -1
    idf = zi(ipi)+zi(ipi+1)

    do kk = 1, zi(ipi+2)
        if ('FX' .eq. zk16(zi(ipi+3)+idf+kk-1)) then
            ! lmat=9, lfct=10
            ldcfct(1) = ipi+lmat-1+lfct*(kk-1)
            iloi = 1
        else if ('FTAN' .eq. zk16(zi(ipi+3)+idf+kk-1)) then
            ldcfct(1) = ipi+lmat-1+lfct*(kk-1)
            iloi = 2
        end if
    end do
    ASSERT(ldcfct(1) .ne. -1)
    ASSERT(iloi .ne. -1)

! ==================== code extrait (début ex-ligne 172) ====================

! Nombre de point de la fonction
    nbvale = zi(ldcfct(1))

! Adresse des informations sur le type de fonction
    jprol = zi(ldcfct(1)+1)
! Adresse des valeurs de la fonction
    jvale = zi(ldcfct(1)+2)
! Type d'écrouissage
    tecro = 0
    if (iloi == 2) then
        call rcvala(jmat, ' ', 'DIS_ECRO_TRAC', 0, ' ', &
 &                  [0.0d0], 1, ['ECROUISSAGE'], valres, codret, &
 &                  1)
        tecro = nint(valres(1))
    end if
! Pour l'intégration
    ldcfct(1) = nbvale
    ldcfct(2) = jprol
    ldcfct(3) = jvale
    ldcfct(4) = iloi
    ldcfct(5) = tecro
! les incréments de déplacement sont nuls
!     ==> récupération de la matrice tangente précédente, si possible
!     ==> si pas possible, pente initiale de la courbe
    if (for_discret%lMatrPred) then
        if (iloi == 1) then
! La tangente est donnée par la pente initiale
            raide(1) = zr(jvale+nbvale+1)/zr(jvale+1)
! Tangente précédente
            if (abs(zr(ivarim-1+15)) .gt. r8miem()) then
                raide(1) = zr(ivarim-1+15)
            end if
        else
! La tangente est donnée par la pente initiale
            raide(2) = zr(jvale+nbvale+1)/zr(jvale+1)
            raide(3) = raide(2)
! Tangente précédente
            if (abs(zr(ivarim-1+16)) .gt. r8miem()) then
                raide(2) = zr(ivarim-1+16)
            end if
            if (abs(zr(ivarim-1+17)) .gt. r8miem()) then
                raide(3) = zr(ivarim-1+17)
            end if
        end if
        resu(1) = zr(icontm)
        resu(2) = zr(icontm+1)
        resu(3) = zr(icontm+2)
        goto 800
    end if
!     !
! loi de comportement non-linéaire : récupération du temps + et - , calcul de dt
    call jevech('PINSTPR', 'L', jtp)
    call jevech('PINSTMR', 'L', jtm)
    temps0 = zr(jtm)
    temps1 = zr(jtp)
    dtemps = temps1-temps0
! contrôle de rk5 : découpage successif, erreur maximale
    call jevech('PCARCRI', 'L', icarcr)
! nombre d'itérations maxi (ITER_INTE_MAXI=-20 par défaut)
    nbdecp = abs(nint(zr(icarcr)))

! tolérance de convergence (RESI_INTE)
    errmax = zr(icarcr+2)
! équations du système :
!                  1 2 3       4 5 6  7      8         9 10 11  12 13 14
!     y0   : force(x,y,z) depl(x,y,z) Dissip pcum deplp(x,y,z)  X(x,y,z)   15 16 17
!     vari : force(x,y,z) depl(x,y,z) Dissip pcum deplp(x,y,z) dX(x,y,z)  tangente(x,y,z)
    y0(:) = 0.0d0
    dy0(:) = 0.0d0
    if (iloi == 1) then
! Déplacements précédent + incréments
        if (for_discret%nno == 1) then
            y0(4) = for_discret%ulm(1)
            dy0(4) = for_discret%dul(1)/dtemps
        else
            y0(4) = (for_discret%ulm(1+for_discret%nc)-for_discret%ulm(1))
            dy0(4) = (for_discret%dul(1+for_discret%nc)-for_discret%dul(1))/dtemps
        end if
! Efforts précédents
        y0(1) = zr(icontm)
    else
! Déplacements précédent + incréments
        if (for_discret%nno == 1) then
            y0(5) = for_discret%ulm(2)
            dy0(5) = for_discret%dul(2)/dtemps
            y0(6) = for_discret%ulm(3)
            dy0(6) = for_discret%dul(3)/dtemps
        else
            y0(5) = (for_discret%ulm(2+for_discret%nc)-for_discret%ulm(2))
            dy0(5) = (for_discret%dul(2+for_discret%nc)-for_discret%dul(2))/dtemps
            y0(6) = (for_discret%ulm(3+for_discret%nc)-for_discret%ulm(3))
            dy0(6) = (for_discret%dul(3+for_discret%nc)-for_discret%dul(3))/dtemps
        end if
! Efforts précédents
        y0(2) = zr(icontm+1)
        y0(3) = zr(icontm+2)
    end if
! Récupération des variables internes

    y0(7:14) = zr(ivarim-1+7:ivarim-1+14)
! Le seuil élastique et le déplacement correspondant
    ldcpar(1) = zr(jvale+nbvale+1)
    ldcpar(2) = zr(jvale+1)

! Norme pour le critère d'erreur
    ynorme(1:3) = ldcpar(1)/10.0d0
    ynorme(4:6) = ldcpar(2)/10.0d0
    ynorme(7) = ldcpar(1)*ldcpar(2)/100.0d0
    ynorme(8:11) = ldcpar(2)/10.0d0
    ynorme(12:14) = ldcpar(1)/100.0d0

    call rk5adp(nbequa, ldcpar, ldcfct, ldccar, temps0, &
 &              dtemps, nbdecp, errmax, y0, dy0, &
 &              disc_isotr, resu, iret, ynorme)

! resu(1:nbeq)            : variables intégrées
! resu(nbeq+1:2*nbeq)     : d(resu)/d(t) a t+dt
    if (iret .ne. 0) goto 999
    !
! calcul de la tangente au comportement
    if (iloi == 1) then
        raide(1) = ldcpar(1)/ldcpar(2)
        if (abs(resu(nbequa+4)) .gt. precis) then
            raide(1) = resu(nbequa+1)/resu(nbequa+4)
        end if
    else
        raide(2) = ldcpar(1)/ldcpar(2)
        raide(3) = raide(2)
        if (abs(resu(nbequa+5)) .gt. precis) then
            raide(2) = resu(nbequa+2)/resu(nbequa+5)
        end if
        if (abs(resu(nbequa+6)) .gt. precis) then
            raide(3) = resu(nbequa+3)/resu(nbequa+6)
        end if
    end if
! actualisation de la matrice quasi-tangente
800 continue
    !
! Actualisation des termes diagonaux
    call diklvraid(for_discret%nomte, klv, raide)
! Actualisation de la matrice quasi-tangente
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        if (for_discret%ndim .eq. 3) then
            call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        else if (for_discret%ndim .eq. 2) then
            call ut2mlg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        end if
    end if
    !
    if (for_discret%lVect .or. for_discret%lSigm) then
! demi-matrice klv transformée en matrice pleine klc
        call vecma(klv, for_discret%nbt, klc, neq)

! calcul de fl = klc.dul (incrément d'effort)
        call pmavec('ZERO', neq, klc, for_discret%dul, fl)
    end if
! Calcul des efforts généralisés
    if (for_discret%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, for_discret%nc
                zr(icontp-1+ii) = fl(ii)+zr(icontm-1+ii)
            end do
            if (iloi == 1) then
                zr(icontp-1+1) = resu(1)
            else
                zr(icontp-1+2) = resu(2)
                zr(icontp-1+3) = resu(3)
            end if
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                zr(icontp-1+ii) = -fl(ii)+zr(icontm-1+ii)
                zr(icontp-1+ii+for_discret%nc) = fl(ii+for_discret%nc)+&
                &                                zr(icontm-1+ii+for_discret%nc)
            end do
            if (iloi == 1) then
                zr(icontp-1+1) = resu(1)
                zr(icontp-1+1+for_discret%nc) = resu(1)
            else
                zr(icontp-1+2) = resu(2)
                zr(icontp-1+3) = resu(3)
                zr(icontp-1+2+for_discret%nc) = resu(2)
                zr(icontp-1+3+for_discret%nc) = resu(3)
            end if
        end if
    end if
! Calcul des forces nodales
    if (for_discret%lVect) then
        call jevech('PVECTUR', 'E', ifono)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, for_discret%nc
                fl(ii) = fl(ii)+zr(icontm-1+ii)
            end do
            if (iloi == 1) then
                fl(1) = resu(1)
            else
                fl(2) = resu(2)
                fl(3) = resu(3)
            end if
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                fl(ii) = fl(ii)-zr(icontm-1+ii)
                fl(ii+for_discret%nc) = fl(ii+for_discret%nc)+zr(icontm-1+ii+for_discret%nc)
            end do
            if (iloi == 1) then
                fl(1) = -resu(1)
                fl(1+for_discret%nc) = resu(1)
            else
                fl(2) = -resu(2)
                fl(3) = -resu(3)
                fl(2+for_discret%nc) = resu(2)
                fl(3+for_discret%nc) = resu(3)
            end if
        end if
! forces nodales aux noeuds 1 et 2 (repère global)
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        end if
    end if
! mise à jour des variables internes
    if (for_discret%lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        zr(ivarip-1+1:ivarip-1+14) = resu(1:14)
        zr(ivarip-1+15) = raide(1)
        zr(ivarip-1+16) = raide(2)
        zr(ivarip-1+17) = raide(3)
        if (for_discret%nno .eq. 2) then
            zr(ivarip-1+20:ivarip-1+36) = zr(ivarip-1+1:ivarip-1+17)
        end if
    end if
999 continue
end subroutine diisotrope_uni
