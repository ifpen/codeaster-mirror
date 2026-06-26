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

subroutine te0147(option, nomte)
!
!------------------------------------------------------------------------------------------
! FONCTION REALISEE:
!
!      CALCUL DU TAUX DE RESTITUTION D'ENERGIE ELEMENTAIRE
!      BORDS ELEMENTS ISOPARAMETRIQUES 2D/3D AVEC CHARGEMENT DE BORD
!      PRESSION-CISAILLEMENT ET FORCE REPARTIE
!
!      ELEMENTS DE BORDS ISOPARAMETRIQUES 2D/3D
!
!      OPTION : 'CALC_G'          (LOCAL,CHARGES REELLES)
!               'CALC_G_F'        (LOCAL,CHARGES FONCTIONS)
!               'CALC_K_G'        (LOCAL,CHARGES REELLES)
!               'CALC_K_G_F'      (LOCAL,CHARGES FONCTIONS)
!               'CALC_KJ_G'       (LOCAL, CHARGES REELES)
!               'CALC_KJ_G_F'     (LOCAL, CHARGES FONCTIONS)
!
! ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!------------------------------------------------------------------------------------------
!
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/coor_cyl.h"
#include "asterfort/deffk.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/lteatt.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvad2.h"
#include "asterc/r8prem.h"
#include "asterfort/utmess.h"
#include "asterfort/tecach.h"
#include "asterfort/thetapdg.h"
!
! =====================================================================
!                       DECLARATION DES VARIABLES
! =====================================================================
!
    integer(kind=8)           :: i, j, kp, k, ind, iret
    integer(kind=8)           :: ndim, nno, nnos, npg, compt
    integer(kind=8)           :: ipoids, ivf, idfde, jgano
    integer(kind=8)           :: ithet, igthet, igeom, idepl
    integer(kind=8)           :: ipref, itemps, iforf, ipres, iforc
    integer(kind=8)           :: icode, imate, jlsn, jlst, ibalo, ideg, ilag
    integer(kind=8)           :: nbpara, reeldim, icodre(3), iiscplan
    real(kind=8)      :: xno1, xno2, yno1, yno2, d1, d2
    real(kind=8)      :: epsi, valpar(4), coor(18)
    real(kind=8)      :: a1(3), a2(3), a3(3), i1(3), i2(3), depl(3)
    real(kind=8)      :: dford1(3), dford2(3), dfor(3), coorg(3), forcg(3)
    real(kind=8)      :: vf, dfde, dxde, dyde, dsde, poids, dsde2
    real(kind=8)      :: th1, th2, dth1d1, dth2d2, th0, t(3), tcla
    real(kind=8)      :: gradth0, gradt(3), dfxde, dfyde, fx, fy, dtdm(3, 4)
    real(kind=8)      :: press, presg(2), prod, divt, forc, thet, absno
    real(kind=8)      :: a1norm, a3norm, i2norm, dfdx(9), dfdy(9)
    real(kind=8)      :: presno, cisano, p(3, 3), invp(3, 3), devres(3)
    real(kind=8)      :: coeff_K1K2, coeff_K3, e, nu, mu, lsng, lstg
    real(kind=8)      :: ka, phig, prsc, rg, valres(3), tcla1, tcla2, tcla3
    real(kind=8)      :: u1g(3), u2g(3), u3g(3), prod1, prod2, fkpo(3, 3)
    character(len=4)  :: fami
    character(len=8)  :: nompar(4), discr
    character(len=16) :: option, nomte, nomres(3)
!
    aster_logical :: axi, fonc, l_not_zero, iscplan
!
    real(kind=8), pointer :: presn(:) => null()
    real(kind=8), pointer :: forcn(:) => null()
    real(kind=8), pointer :: ffp(:) => null()
    real(kind=8), pointer :: dfdi(:) => null()
!
! =====================================================================
!                       INITIALISATION PARAMETRES
!                       Cas 2D/3D
!                       Option calc_G/calc_K_G/calc_KJ_G
! =====================================================================
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
!-- Initialisation des paramètres
!
!-- Dimention du modèle et pas de l'element de bord
    reeldim = ndim+1
    nbpara = reeldim+1
    epsi = r8prem()
    axi = ASTER_FALSE
    tcla = 0.d0
    tcla1 = 0.d0
    tcla2 = 0.d0
    tcla3 = 0.d0
    coeff_K1K2 = 0.d0
    coeff_K3 = 0.d0
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
!
!-- Allocation dynamique : vecteurs 2D/3D
    if (reeldim .eq. 3) then
        AS_ALLOCATE(vr=presn, size=nno)
        AS_ALLOCATE(vr=dfdi, size=reeldim*nno)
    else
!------ Pression/cisaillement
        AS_ALLOCATE(vr=presn, size=reeldim*nno)
    end if
    AS_ALLOCATE(vr=forcn, size=reeldim*nno)
    AS_ALLOCATE(vr=ffp, size=nno)
!
!-- Cas 2D / Element de bord 1D
    if (reeldim .eq. 2) then
        if (lteatt('AXIS', 'OUI')) axi = ASTER_TRUE
    end if
!
! =====================================================================
!                       CALCUL DE THETA
! =====================================================================
!
!-- Recuperation des parametres pour créer theta
    call jevech('PTHETAR', 'L', ithet)  ! champ d'entree
!
!-- Recuperation du champ de sortie
    call jevech('PGTHETA', 'E', igthet) ! champ de sortie
!
!   Recuperation du degré si LEGENDRE ou abscisse si LAGRANGE
    if (reeldim .eq. 3) then
!
        call tecach('ONO', 'PDEG', 'L', iret, iad=ideg)
!
        if (iret .eq. 0) then
            discr = "LEGENDRE"
        else
            discr = "LINEAIRE"
            call jevech('PLAG', 'L', ilag)
        end if
!
    else
        discr = "2D"
    end if
!
!-- Valeur de theta0(r) nulle sur element : pas de calcul
    compt = 0
    do i = 1, nno
        thet = zr(ithet-1+6*(i-1)+1)
        if (thet .lt. epsi) compt = compt+1
    end do
    if (compt .eq. nno) goto 999
!-- Si LAGRANGE, test sur les abscisses curvilignes
    compt = 0
    if (discr == "LINEAIRE") then
        do i = 1, nno
            absno = zr(ithet-1+6*(i-1)+5)
!---------- Cas fond fermé
            if (zr(ilag) .gt. zr(ilag+1)) then
                if (absno .ge. zr(ilag) .and. absno .le. zr(ithet-1+6*(i-1)+6)) then
                    compt = compt+1
                elseif (absno .ge. zr(ilag+1) .and. absno .le. zr(ilag+2)) then
                    compt = compt+1
                end if
            else
                if (absno .ge. zr(ilag) .and. absno .lt. zr(ilag+2)) compt = compt+1
            end if
        end do
        if (compt .eq. 0) goto 999
    end if
!
! =====================================================================
!                  RECUPERATION DES CHAMPS LOCAUX
! =====================================================================
!

    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLAR', 'L', idepl)
    call jevech('PMATERC', 'L', imate)
    if (option .eq. 'CALC_K_G' .or. option .eq. 'CALC_K_G_F') then
        call jevech('PBASLOR', 'L', ibalo)
        call jevech('PLSN', 'L', jlsn)
        call jevech('PLST', 'L', jlst)
    end if
!
    if (reeldim .eq. 2) then
        call tecach('ONO', 'ISCPLAN', 'L', iret, iad=iiscplan)

        if (zi(iiscplan) == 1) then
            iscplan = .true.
        else
            iscplan = .false.
        end if
    end if
!
    if (option .eq. 'CALC_G_F' .or. option .eq. 'CALC_K_G_F' .or. option .eq. 'CALC_KJ_G_F') then
        fonc = ASTER_TRUE
        call jevech('PPRESSF', 'L', ipref)
        call jevech('PINSTR', 'L', itemps)
        if (reeldim .eq. 3) then
            call jevech('PFF2D3D', 'L', iforf)
            nompar(1) = 'X'
            nompar(2) = 'Y'
            nompar(3) = 'Z'
            nompar(4) = 'INST'
            valpar(4) = zr(itemps)
        else
            call jevech('PFF1D2D', 'L', iforf)
            nompar(1) = 'X'
            nompar(2) = 'Y'
            nompar(3) = 'INST'
            valpar(3) = zr(itemps)
        end if
    else
        fonc = ASTER_FALSE
        call jevech('PPRESSR', 'L', ipres)
        if (reeldim .eq. 3) then
            call jevech('PFR2D3D', 'L', iforc)

        else
            call jevech('PFR1D2D', 'L', iforc)
        end if
    end if
!
! =====================================================================
!      SI CHARGE FONCTION RECUPERATION DES VALEURS AUX NOEUDS
! =====================================================================
!
    if (fonc) then
!------ Dimension réel du modèle et pas de l'element de bord
        do i = 1, nno
!
            do j = 1, reeldim
                valpar(j) = zr(igeom+reeldim*(i-1)+j-1)
            end do
!
            if (reeldim .eq. 3) then
                call fointe('FM', zk8(ipref), nbpara, nompar, valpar, &
                            presn(i), icode)
            else
                do j = 1, reeldim
                    call fointe('FM', zk8(ipref+j-1), nbpara, nompar, valpar, &
                                presn(2*(i-1)+j), icode)
                end do
            end if
!
            do j = 1, reeldim
                call fointe('FM', zk8(iforf+j-1), nbpara, nompar, valpar, &
                            forcn(reeldim*(i-1)+j), icode)
            end do
        end do
    end if
!
! =====================================================================
!     CALCUL DU REPERE LOCAL ORTHONORME DANS LE CAS 3D
! =====================================================================
!
    if (reeldim .eq. 3) then
!
        do j = 1, reeldim
            a1(j) = zr(igeom+3*(2-1)+j-1)-zr(igeom+3*(1-1)+j-1)
            a2(j) = zr(igeom+3*(3-1)+j-1)-zr(igeom+3*(1-1)+j-1)
        end do
!
!------ Determination du troisieme vecteur de la base
        a3(1) = a1(2)*a2(3)-a1(3)*a2(2)
        a3(2) = a1(3)*a2(1)-a1(1)*a2(3)
        a3(3) = a1(1)*a2(2)-a1(2)*a2(1)
!
!------ Calcul du repere orthonorme
        i2(1) = a3(2)*a1(3)-a3(3)*a1(2)
        i2(2) = a3(3)*a1(1)-a3(1)*a1(3)
        i2(3) = a3(1)*a1(2)-a3(2)*a1(1)
!
        a1norm = sqrt(a1(1)*a1(1)+a1(2)*a1(2)+a1(3)*a1(3))
        i2norm = sqrt(i2(1)*i2(1)+i2(2)*i2(2)+i2(3)*i2(3))
        a3norm = sqrt(a3(1)*a3(1)+a3(2)*a3(2)+a3(3)*a3(3))
!
        do i = 1, reeldim
            i1(i) = a1(i)/a1norm
            i2(i) = i2(i)/i2norm
            a3(i) = a3(i)/a3norm
        end do
!
        do i = 1, nno
            coor(2*i-1) = 0.d0
            coor(2*i) = 0.d0
            do j = 1, reeldim
                coor(2*i-1) = coor(2*i-1)+(zr(igeom+3*(i-1)+j-1)- &
                                           zr(igeom+j-1))*i1(j)
                coor(2*i) = coor(2*i)+(zr(igeom+3*(i-1)+j-1)- &
                                       zr(igeom+j-1))*i2(j)
            end do
        end do
    end if
!
! =====================================================================
!           BOUCLE PRINCIPALE SUR LES POINTS DE GAUSS
! =====================================================================
!
    do kp = 1, npg
!
        k = (kp-1)*nno
!
        ! ===========================================
        !              INITIALISATION
        ! ===========================================
!
        depl(:) = 0.d0
        dford1(:) = 0.d0
        dford2(:) = 0.d0
        dfor(:) = 0.d0
        coorg(:) = 0.d0
        forcg(:) = 0.d0
        dxde = 0.d0
        dyde = 0.d0
        th1 = 0.d0
        th2 = 0.d0
        dth1d1 = 0.d0
        dth2d2 = 0.d0
        th0 = 0.d0
        t(:) = 0.d0
        gradth0 = 0.d0
        gradt(:) = 0.d0
        dfxde = 0.d0
        dfyde = 0.d0
        fx = 0.d0
        fy = 0.d0
        press = 0.d0
        presg(:) = 0.d0
        p(:, :) = 0.d0
        invp(:, :) = 0.d0
        prsc = 0.d0
        u1g(:) = 0.d0
        u2g(:) = 0.d0
        u3g(:) = 0.d0
        lsng = 0.d0
        lstg = 0.d0
!
!------ Coordonnées des points de Gauss
        do i = 1, nno
            do j = 1, reeldim
                coorg(j) = coorg(j)+zr(ivf+k+i-1)* &
                           zr(igeom+reeldim*(i-1)+j-1)
            end do
        end do
!
        ! =================================================
        ! RECUPERATION DES FF, DERIVEES FF, DEPL, POIDS
        ! =================================================
!
        if (reeldim .eq. 3) then
            call dfdm2d(nno, kp, ipoids, idfde, coor, &
                        poids, dfdx, dfdy)
            do i = 1, nno
                dfdi(i) = dfdx(i)
                dfdi(i+nno) = dfdy(i)
                dfdi(i+2*nno) = 0.d0
                do j = 1, reeldim
                    depl(j) = depl(j)+zr(ivf+k+i-1)*zr(idepl+reeldim*(i-1)+j-1)
                end do
            end do
        else
            do i = 1, nno
!-------------- Voir pour utiliser dfdmd1 et mutualiser la routine
!-------------- 2D : pression/cisaillement donc fointe différent entre les deux
                vf = zr(ivf+k+i-1)
                dfde = zr(idfde+k+i-1)
                dxde = dxde+dfde*zr(igeom+2*(i-1))
                dyde = dyde+dfde*zr(igeom+2*(i-1)+1)
!-------------- Paramètres de theta et gradtheta au pdg
                th0 = th0+vf*zr(ithet-1+6*(i-1)+1)
                gradth0 = gradth0+dfde*zr(ithet-1+6*(i-1)+1)
                do j = 1, reeldim
                    t(j) = t(j)+vf*zr(ithet-1+6*(i-1)+j+1)
                    gradt(j) = gradt(j)+dfde*zr(ithet-1+6*(i-1)+j+1)
                    depl(j) = depl(j)+vf*zr(idepl+reeldim*(i-1)+j-1)
                end do
            end do
!
!---------- Jacobien dans le cas 2D (element 1D)
            dsde = sqrt(dxde*dxde+dyde*dyde)
            dsde2 = dsde*dsde
!
            if (axi) then
                if (coorg(1) .lt. r8prem()) then
                    call utmess('F', 'RUPTURE0_56')
                else
                    poids = zr(ipoids+kp-1)*dsde*coorg(1)
                end if
            else
                poids = zr(ipoids+kp-1)*dsde
            end if
!
        end if
!
        ! ===========================================
        !    CALCUL DU CHARGEMENT ET DE SON GRADIENT
        ! ===========================================
!
        if (fonc) then
!
            do j = 1, reeldim
                valpar(j) = coorg(j)
            end do
!
            ! ==============
            !  FONCTION 3D
            ! ==============
            if (reeldim .eq. 3) then
                call fointe('FM', zk8(ipref), nbpara, nompar, valpar, &
                            press, icode)
                do j = 1, reeldim
                    call fointe('FM', zk8(iforf+j-1), nbpara, nompar, valpar, &
                                forcg(j), icode)
                end do
!
!-------------- Valeurs de la derivée du chargement aux pdg
                do i = 1, nno
                    do j = 1, 3
                        dford1(j) = dford1(j)+(forcn(3*(i-1)+j)-presn(i)*a3(j))*dfdx(i)
                        dford2(j) = dford2(j)+(forcn(3*(i-1)+j)-presn(i)*a3(j))*dfdy(i)
                    end do
                end do
!
                ! ==============
                !  FONCTION 2D
                ! ==============
            else
                do j = 1, reeldim
                    call fointe('FM', zk8(ipref+j-1), nbpara, nompar, valpar, &
                                presg(j), icode)
                    call fointe('FM', zk8(iforf+j-1), nbpara, nompar, valpar, &
                                forcg(j), icode)
                end do
!
!-------------- Valeurs de la derivée du chargement aux pdg
                do i = 1, nno
                    dfde = zr(idfde+k+i-1)
                    presno = presn(2*(i-1)+1)
                    cisano = presn(2*(i-1)+2)
                    dfxde = dfxde+dfde*(forcn(2*(i-1)+1)- &
                                        (dyde*presno-dxde*cisano)/dsde)
                    dfyde = dfyde+dfde*(forcn(2*(i-1)+2)+ &
                                        (dxde*presno+dyde*cisano)/dsde)
                end do
            end if
        else
            ! ==============
            !     CAS 3D
            ! ==============
            if (reeldim .eq. 3) then
                do i = 1, nno
                    press = press+zr(ipres+i-1)*zr(ivf+k+i-1)
                    do j = 1, reeldim
                        forcg(j) = forcg(j)+zr(iforc+3*(i-1)+j-1)*zr(ivf+k+i-1)
                    end do
                end do
                ! ==============
                !    CAS 2D
                ! ==============
            else
                do i = 1, nno
                    do j = 1, reeldim
                        presg(j) = presg(j)+zr(ipres+2*(i-1)+j-1)*zr(ivf+k+i-1)
                        forcg(j) = forcg(j)+zr(iforc+2*(i-1)+j-1)*zr(ivf+k+i-1)
                    end do
                end do
            end if
!
        end if
!
!       ! ===========================================
        !  CONSTRUCTION DE THETA , GRADTHETA AU PDG
        ! ===========================================
!
        if (reeldim .eq. 3) then
!
!---------- Calcul de DTDM pour LEGENDRE ou LAGRANGE
            call thetapdg(reeldim, nno, discr, &
                          zr(ivf+k-1+1:ivf+k-1+nno), &
                          dfdi, ideg, ilag, ithet, dtdm)
!
!---------- Calcul de theta dans la base I1,I2
            th1 = dot_product(dtdm(1:3, 4), i1)
            th2 = dot_product(dtdm(1:3, 4), i2)
            dth1d1 = dot_product(dtdm(1:3, 1), i1)
            dth2d2 = dot_product(dtdm(1:3, 2), i2)
!
!---------- Calcul de la divergence
            divt = dth1d1+dth2d2
        else
            th1 = (th0*t(1)*dxde)/dsde2
            th2 = (th0*t(2)*dyde)/dsde2
            dth1d1 = ((gradth0*t(1)+th0*gradt(1))*dxde)/dsde2
            dth2d2 = ((gradth0*t(2)+th0*gradt(2))*dyde)/dsde2
!
!---------- Calcul de la divergence
            divt = dth1d1+dth2d2
        end if
!
        ! ===========================================
        !      RECUPERATION DES DONNEES MATERIAU
        ! ===========================================
!
        call rcvad2(fami, kp, 1, '+', zi(imate), 'ELAS', &
                    3, nomres, valres, devres, icodre)
!
        if ((icodre(1) .ne. 0) .or. (icodre(2) .ne. 0)) then
            call utmess('F', 'RUPTURE1_25')
        end if
!
        e = valres(1)
        nu = valres(2)
        mu = e/(2.d0*(1.d0+nu))
!
        if (reeldim .eq. 3 .or. lteatt('AXIS', 'OUI') .or. (.not. iscplan)) then
            ka = 3.d0-4.d0*nu
            coeff_K1K2 = e/(1.d0-nu*nu)
            coeff_K3 = 2.d0*mu
        else
!----------- Contrainte plane
            ka = (3.d0-nu)/(1.d0+nu)
            coeff_K1K2 = e
        end if
!
        ! ===========================================
        !         CALCUL DES SIFS ; OPTION K
        ! ===========================================
!
        if (option .eq. 'CALC_K_G' .or. option .eq. 'CALC_K_G_F') then
!
            ! ===========================================
            !      CALCUL DES COORDONNEE CYLINDRIQUE
            ! ===========================================
!
            do i = 1, nno
                lsng = lsng+zr(jlsn-1+i)*zr(ivf+k+i-1)
                lstg = lstg+zr(jlst-1+i)*zr(ivf+k+i-1)
                ffp(i) = zr(ivf-1+nno*(kp-1)+i)
            end do
!
            if (reeldim == 2) then
!-------------- On construit le vecteur normal à partir du vecteur direction
!-------------- pour s'assurer d'avoir un repère direct (cas 2d uniquement)
                do i = 1, nno
                    zr(ibalo-1+6*(i-1)+5) = -zr(ibalo-1+6*(i-1)+4)
                    zr(ibalo-1+6*(i-1)+6) = zr(ibalo-1+6*(i-1)+3)
                end do
            end if
!
            call coor_cyl(reeldim, nno, zr(ibalo), zr(igeom), ffp, &
                          p, invp, rg, phig, l_not_zero)
!
            if (.not. l_not_zero) then
                call utmess('F', 'RUPTURE1_75')
            end if
!
!---------- On utilise les level sets pour déterminer le signe de phi
!---------- On détermine si on est sur levre sup ou sur levre inf
            if ((abs(lsng) .lt. 1.0d-8) .and. (lstg .lt. 0.0d0)) then
!
                if (reeldim .eq. 2) then
                    xno1 = zr(igeom)
                    yno1 = zr(igeom+1)
                    xno2 = zr(igeom+2)
                    yno2 = zr(igeom+3)
                    d1 = ((xno1-zr(ibalo-1+1))*(xno1-zr(ibalo-1+1)))+ &
                         ((yno1-zr(ibalo-1+2))*(yno1-zr(ibalo-1+2)))
                    d2 = ((xno2-zr(ibalo-1+1))*(xno2-zr(ibalo-1+1)))+ &
                         ((yno2-zr(ibalo-1+2))*(yno2-zr(ibalo-1+2)))
                    if (d2 .gt. d1) then
                        phig = -1.0d0*abs(phig)
                    else
                        phig = abs(phig)
                    end if
                else
!------------------ Produit scalaire entre vecteur normale au fond de fissures
!------------------ et la normale de l'element
                    do i = 1, 3
                        prsc = prsc+p(i, 2)*a3(i)
                    end do
                    if (prsc .gt. r8prem()) then
                        phig = -1.0d0*abs(phig)
                    else
                        phig = abs(phig)
                    end if
                end if
!
            end if
!
            ! ===========================================
            !      CALCUL DES CHAMPS SINGULIERS
            ! ===========================================
!
            call deffk(ka, mu, rg, phig, reeldim, fkpo)
!
            do i = 1, reeldim
                do ind = 1, reeldim
                    u1g(i) = u1g(i)+p(i, ind)*fkpo(1, ind)
                    u2g(i) = u2g(i)+p(i, ind)*fkpo(2, ind)
                    if (reeldim .eq. 3) then
                        u3g(i) = u3g(i)+p(i, ind)*fkpo(3, ind)
                    end if
                end do
            end do
        end if
!
        ! ===========================================
        !   CALCUL DU TAUX DE RESTITUTION G ET K*
        ! ===========================================
!
        divt = dth1d1+dth2d2
!
        if (axi) divt = divt+((th0*t(1))/coorg(1))
!
        if (reeldim .eq. 3) then
            do j = 1, 3
                dfor(j) = dfor(j)+dford1(j)*th1+dford2(j)*th2
            end do
            do j = 1, 3
                forc = forcg(j)-press*a3(j)
                tcla = tcla+poids*(forc*divt+dfor(j))*depl(j)
                if (option .eq. 'CALC_K_G' .or. option .eq. 'CALC_K_G_F') then
                    tcla1 = tcla1+0.5d0*poids*(forc*divt+dfor(j))*u1g(j)
                    tcla2 = tcla2+0.5d0*poids*(forc*divt+dfor(j))*u2g(j)
                    tcla3 = tcla3+0.5d0*poids*(forc*divt+dfor(j))*u3g(j)
                end if
            end do
        else
            fx = forcg(1)-(dyde*presg(1)-dxde*presg(2))/dsde
            fy = forcg(2)+(dxde*presg(1)+dyde*presg(2))/dsde
            prod = (divt*fx+dfxde*(th1+th2))*depl(1)+(divt*fy+dfyde*(th1+th2))*depl(2)
            tcla = tcla+prod*poids
!
            if (option .eq. 'CALC_K_G' .or. option .eq. 'CALC_K_G_F') then
                prod1 = (divt*fx+dfxde*(th1+th2))*u1g(1)+(divt*fy+dfyde*(th1+th2))*u1g(2)
                prod2 = (divt*fx+dfxde*(th1+th2))*u2g(1)+(divt*fy+dfyde*(th1+th2))*u2g(2)
                tcla1 = tcla1+0.5d0*poids*prod1
                tcla2 = tcla2+0.5d0*poids*prod2
            end if
        end if
!
    end do

    ! ===========================================================
    !         CALCUL DE K1 A PARTIR DE LA FORMULE D'IRWIN ET DE G
    !         OPTION G
    ! ===========================================================
!
    if (option == 'CALC_G' .or. option == 'CALC_G_F') then
        tcla1 = tcla
    end if

!
!-- Exit sur valeur de theta nulle sur element
999 continue
!
!-- Assemblage final
!
    zr(igthet) = tcla
    zr(igthet+1) = tcla1*coeff_K1K2
    zr(igthet+2) = tcla2*coeff_K1K2
    zr(igthet+3) = tcla3*coeff_K3

    zr(igthet+4) = tcla1*sqrt(coeff_K1K2)
    zr(igthet+5) = tcla2*sqrt(coeff_K1K2)
    zr(igthet+6) = tcla3*sqrt(coeff_K3)
!
    AS_DEALLOCATE(vr=forcn)
    AS_DEALLOCATE(vr=presn)
    AS_DEALLOCATE(vr=ffp)
    if (reeldim .eq. 3) then
        AS_DEALLOCATE(vr=dfdi)
    end if
!
end subroutine
