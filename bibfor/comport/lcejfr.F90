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
! aslint: disable=W0413,W1306
!
subroutine lcejfr(BEHinteg, fami, kpg, ksp, ndim, &
                  mate, option, epsm, deps, sigma, &
                  dsidep, vim, vip, typmod, instam, &
                  instap)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/matinv.h"
#include "asterfort/pmavec.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    type(Behaviour_Integ), intent(in) :: BEHinteg
    integer(kind=8) :: mate, ndim, kpg, ksp
    real(kind=8) :: epsm(6), deps(6)
    real(kind=8) :: sigma(6), dsidep(6, 6)
    real(kind=8) :: vim(*), vip(*), instam, instap
    character(len=8) :: typmod(*)
    character(len=16) :: option
    character(len=*) :: fami
!
!-----------------------------------------------------------------------
!     LOI DE COMPORTEMENT  MOHR-COULOMB
!     POUR LES ELEMENTS DE JOINT ET JOINT_HYME 2D ET 3D.
!
! DEUX TYPES DE CALCUL SONT POSSIBLES:
!
! 1) MODELISATIONS *_JOINT_HYME, AVEC PARAMETRE HYDRO PRESENTS
!    CALCUL COUPLE HYDRO-MECANIQUE SUR UN MAILLAGE QUADRATIQUE
!
! 2) MODELISATIONS *_JOINT, PAS DE PARAMETRES HYDRO POSSIBLE
!    CALCUL MECA (AVEC ENVENTUELLEMENT PRES_FLUIDE)
!    SUR UN MAILLAGE LINEAIRE OU QUADRATIQUE
!
! IN : EPSM SAUT INSTANT MOINS ET GRAD PRESSION ET PRES FLUIDE SI HYME
! IN : DEPS INCR DE SAUT ET INCR GRAD PRESS ET INC PRES FL SI HYME
! IN : MATE, OPTION, VIM, COOROT,INSTAM, INSTAP
! OUT : SIGMA , DSIDEP , VIP
!-----------------------------------------------------------------------
    integer(kind=8) :: nbpa
    parameter(nbpa=10)
    integer(kind=8) :: cod(nbpa)
    integer(kind=8) :: i, j, n, ifplas, kronec, ifouv
    real(kind=8) :: kn, kt, kappa, mu, adhe, a(ndim), plasti(ndim), dplas(ndim)
    real(kind=8) :: lambda, dlam, val(nbpa), criter
    real(kind=8) :: abstau, tau(ndim), coefd, coefhd, r8bid
    real(kind=8) :: inst, valpar(ndim+1), rhof, visf, amin, presfl, presg
    real(kind=8) :: gp(ndim-1), gploc(ndim), gpglo(ndim), fhloc(ndim)
    real(kind=8) :: fhglo(ndim), doset, oset, sciage
    real(kind=8) :: coorot(ndim+ndim*ndim), invrot(ndim, ndim), rigart
    real(kind=8) :: sigma_eff
    character(len=8) :: nompar(ndim+1)
    character(len=16) :: nom(nbpa)
    character(len=1) :: poum
    aster_logical :: resi, rigi, elas, ifpahm, ifhyme
    blas_int :: b_incx, b_incy, b_n
!
! OPTION CALCUL DU RESIDU OU CALCUL DE LA MATRICE TANGENTE
! CALCUL DE CONTRAINTE (RESIDU)
    resi = option(1:9) .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
! CALCUL DE LA MATRICE TANGEANTE (RIGIDITE)
    rigi = option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RIGI_MECA'
! CALCUL DE LA MATRICE ELASTIQUE A LA PLACE DE LA MATRICE TANGENTE
    elas = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'RIGI_MECA_ELAS'
!
! INDICATEUR AVEC/SANS HYDRO
    if (typmod(2) .eq. 'EJ_HYME') ifhyme = .true.
    if (typmod(2) .eq. 'ELEMJOIN') ifhyme = .false.
!
! #####################################
! RECUPERATION DES PARAMETRES PHYSIQUES
! #####################################
    nom(1) = 'K_N'
    nom(2) = 'MU'
    nom(3) = 'ADHESION'
    nom(4) = 'K_T'
    nom(5) = 'PENA_TANG'
    nom(6) = 'PRES_FLUIDE'
    nom(7) = 'RHO_FLUIDE'
    nom(8) = 'VISC_FLUIDE'
    nom(9) = 'OUV_MIN'
    nom(10) = 'SCIAGE'
!
    if (option .eq. 'RIGI_MECA_TANG') then
        poum = '-'
    else
        poum = '+'
    end if
!
! INSTANT DE CALCUL T- OU T+
    inst = instam
    if (resi) inst = instap
!
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_FROT', 0, ' ', [0.d0], &
                3, nom, val, cod, 2)
! DEFINITION DE PARAMETRES PHYSIQUE:
!     PENTE ELASTIQUE NORMALE
    kn = val(1)
!     COEFFICIENT DE FROTTEMENT
    mu = val(2)
!     ADHESION
    adhe = val(3)
!
! PENTE ELASTIQUE TANGENTIELLE
! (SI ELLE N'EST PAS DEFINI ALORS K_T=K_N)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_FROT', 0, ' ', [0.d0], &
                1, nom(4), val(4), cod(4), 0)
    if (cod(4) .eq. 0) then
        kt = val(4)
    else
        kt = kn
    end if
! PARAMETRE PENA_TANG
! (SI IL N'EST PAS DEFINI ALORS KAPPA=(K_N+K_T)*1E-6)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_FROT', 0, ' ', [0.d0], &
                1, nom(5), val(5), cod(5), 0)
    if (cod(5) .eq. 0) then
        kappa = val(5)
    else
        kappa = (kn+kt)*1.d-6
    end if
!
!
! DEFINITION DES PARAMETRES POUR LA RECUPERATION DES FONCTIONS
    coorot = 0.d0
    do i = 1, ndim
        coorot(i) = BEHinteg%behavESVA%behavESVAGeom%coorElga(kpg, i)
    end do
    do i = 1, ndim*ndim
        coorot(ndim+i) = BEHinteg%behavESVA%behavESVAOther%rotpg(i)
    end do
!
    nompar(1) = 'INST'
    nompar(2) = 'X'
    nompar(3) = 'Y'
    valpar(1) = inst
    valpar(2) = coorot(1)
    valpar(3) = coorot(2)
    if (ndim .eq. 3) then
        nompar(4) = 'Z'
        valpar(4) = coorot(3)
    end if
!
! RECUPERATION DE LA PRESS FLUIDE IMPOSEE (FCT DE L'ESPACE ET DU TEMPS)
!
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_FROT', ndim+1, nompar, valpar, &
                1, nom(6), val(6), cod(6), 0)
    if (cod(6) .eq. 0) then
        presfl = val(6)
    else
        presfl = 0.d0
    end if
!
! RECUPERATION DE LA TAILLE DE SCIE = SCIAGE (FONCTION DE L'ESPACE ET DU TEMPS)
!
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_FROT', ndim+1, nompar, valpar, &
                1, nom(10), val(10), cod(10), 0)
!
    if (cod(10) .eq. 0) then
        sciage = val(10)
    else
        sciage = 0.d0
    end if
!
! RECUPERATION DE LA MASSE VOL ET DE LA VISCO (MODELISATION JOINT HM)
!--------------------------------------------------------------------
!
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_FROT', 0, ' ', [0.d0], &
                1, nom(7), val(7), cod(7), 0)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_FROT', 0, ' ', [0.d0], &
                1, nom(8), val(8), cod(8), 0)
    call rcvalb(fami, kpg, ksp, poum, mate, &
                ' ', 'JOINT_MECA_FROT', 0, ' ', [0.d0], &
                1, nom(9), val(9), cod(9), 0)
!
    if (cod(7) .eq. 0) rhof = val(7)
    if (cod(8) .eq. 0) visf = val(8)
    if (cod(9) .eq. 0) amin = val(9)
!
! INDICATEUR SI LES PARAMETRES HYDRO SONT RENSEIGNES
    ifpahm = (cod(7) .eq. 0) .and. (cod(8) .eq. 0) .and. (cod(9) .eq. 0)
!
! VERIFICATION DE LA PRESENCE/ABSENCE DES PARAMETRES
! EN FONCTION DE LA MODELISATION MECA PUR OU HYDRO MECA
!
    if (ifhyme) then
!       POUR LE CALCUL HYDRO => PAS DE SCIAGE
        if (cod(10) .eq. 0) then
            call utmess('F', 'ALGORITH17_14')
        end if
!       POUR LE CALCUL HYDRO => PAS DE PRES_FLUIDE
        if ((cod(6) .eq. 0)) then
            call utmess('F', 'ALGORITH17_14')
        end if
!       POUR LE CALCUL HYDRO => PRESENCE DE PARA HYDRO
        if (.not. ifpahm) then
            call utmess('F', 'ALGORITH17_15')
        end if
    else
!       POUR LE CALCUL MECA => PAS DE PARAMETRE HYDRO
        if (ifpahm) then
            call utmess('F', 'ALGORITH17_16')
        end if
    end if
!
!
! #####################################
! INITIALISATION DE VARIABLES
! #####################################
!
! CALCUL DU SAUT EN T+ OU T- EN FOCTION DE L'OPTION DE CALCUL
!     A=AM
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, epsm, b_incx, a, b_incy)
!     A=A+DA
    if (resi) then
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deps, b_incx, a, &
                   b_incy)
    end if
!
! DANS LE CAS DU SCIAGE
! INITIALISATION DU POINT D'EQUILIBRE POUR LA LDC (OFFSET)
!-----------------------
! L'EPASSEUR SCIEE EST DIMINUEE DE L'OUVERTURE INITALE DE JOINT
    doset = 0.d0
    sciage = sciage-max(0., epsm(1))
    if (sciage .gt. 0.d0) doset = doset-sciage
    oset = vim(10)+doset
!
! LA LDC EST DEFINIE PAR RAPPORT A NOUVEAU POINT D'EQUILIBRE
    a(1) = a(1)-oset
!
!
! GRADIENT DE PRESSION ET PRESSION EN T- OU T+
    if (ifhyme) then
!
        do n = 1, ndim-1
            gp(n) = epsm(ndim+n)
            if (resi) gp(n) = gp(n)+deps(ndim+n)
        end do
!
        presg = epsm(2*ndim)
        if (resi) presg = presg+deps(2*ndim)
!
    end if
!
! INITIALISATION DE VARIABLE INTERNES
!     LAMBDA = accumutation de deplacement tang en valeur absolue
!     PLASTI = vecteur de deplacement plastique tangentiel
    lambda = vim(1)
    do i = 2, ndim
        plasti(i) = vim(i+1)
    end do
!     IDICATEUR DE PLASTIFICATION A L'INSTANT ACTUEL
    if (elas) then
        ifplas = 0
    else
        ifplas = nint(vim(2))
    end if
!     INDICATEUR D'OUVERTURE
    ifouv = nint(vim(5))
!
! #####################################
! CALCUL DE LA CONTRAINTE
! #####################################
!
!     INITIALISATION DE LA CONTRAINTE A ZERO
    call r8inir(6, 0.d0, sigma, 1)
!
! CALCUL DE LA CONTRAINTE HYDRO : DEBIT (LOI CUBIQUE)
    if (ifhyme) then
        do n = 1, ndim-1
            sigma(ndim+n) = -rhof*gp(n)*(max(amin, a(1)+amin))**3/(12*visf)
        end do
    end if
!
!     CONTRAINTE NORMALE DE CONTACT PENALISE +
!     INDICATEUR D'OUVERTURE COMPLETE
!     (CONTRAINTE TANGENTIEL EST MISE A ZERO)
    if (kn*a(1) .lt. (adhe/mu)) then
        ifouv = 0
        sigma_eff = kn*a(1)
    else
        ifouv = 1
        sigma_eff = adhe/mu
    end if
!
!     CONTRAINTE TANGENTIELLE
    do i = 2, ndim
        sigma(i) = kt*(a(i)-plasti(i))
    end do
!
!     DIRECTION DE GLISSEMENT = SIGMA TANG SANS INCREMENT PLASTIQUE
    do i = 2, ndim
        tau(i) = sigma(i)
    end do
!
!     MODULE DE SIGMA TANGENTE SANS INCREMENT PLASTIQUE
!     NB:SI ABSTAU==0 ON EST TOUJOURS DANS LE REGIME ELASTIQUE
    abstau = 0.d0
    do i = 2, ndim
        abstau = abstau+tau(i)**2
    end do
    abstau = sqrt(abstau)
!
! ###########################################################
! SI ON CALCULE LA RIGIDITE SEULEMENT, ON Y SAUTE DIRECTEMENT
! ###########################################################
    if (.not. resi) goto 5000
!
!     CRITERE DE PLASTICITE  NB: SIGMA(1)<0 EN COMPRESSION
    criter = abstau+mu*sigma_eff-kappa*lambda-adhe
!     VERIFICATION DE CRITERE DE PLASTICITE
    if (criter .le. 0.d0) then
!     PAS DE PLASTICITE
        ifplas = 0
        do i = 2, ndim
            dplas(i) = 0.d0
        end do
        dlam = 0.d0
    else
!     AVEC LA PLASTICITE
        ifplas = 1
        do i = 2, ndim
            dplas(i) = criter/(kt+kappa)*tau(i)/abstau
            sigma(i) = sigma(i)-kt*dplas(i)
        end do
        dlam = criter/(kt+kappa)
    end if
!
!     PRISE EN COMPTE DE LA PRESSION DE FLUIDE EVENTUELLE
!     PRESFL : IMPOSEE, PRESG : CALCULEE
!
    if (ifhyme) then
        sigma(1) = sigma_eff-presg
    else
        sigma(1) = sigma_eff-presfl
    end if
!
!
! ACTUALISATION DES VARIABLES INTERNES
!   V1 :  LE DEPLACEMENT PLASTIQUE CUMULE (SANS ORIENTATION) LAMBDA:
!            LAMBDA NE PEUT QU'AUGMENTER
!   V2 : INDICATEUR DE PLASTIFICATION (0 : NON, 1 : OUI)
!   V3-V4 :  VECTEUR DE DEPLACEMENT TANG PAR RAPPORT AU POINT DE DEPART
!                        (INDIQUE LA POSITION D'EQUILIBRE ACTUELLE)
!   V5    : INDICATEUR D'OUVERTURE COMPLETE
!                         (CONTRAINTE TANGENTIEL EST MIS A ZERO)
!   V6    : MODULE DE LA CONTRAINTE TANGENTE
!   V7 A V9 : VALEURS DU SAUT
!   V10 : EPAISSEUR DU JOINT
!   V11 : CONTRAINTE MECANIQUE NORMALE (SANS PRESSION DE FLUIDE)
!   V12 A V14 : COMPOSANTES DU GRAD DE PRESSION DANS LE REPERE GLOBAL
!   V15 A V17 : COMPOSANTES DU FLUX HYDRO DANS LE REPERE GLOBAL
!   V18 : PRESSION DE FLUIDE INTERPOLEE
!
    vip(1) = lambda+dlam
    vip(2) = ifplas
    do i = 2, ndim
        vip(i+1) = vim(i+1)+dplas(i)
    end do
    vip(5) = max(nint(vim(5)), ifouv)
!
    vip(6) = 0.d0
    do i = 2, ndim
        vip(6) = vip(6)+sigma(i)**2
    end do
    vip(6) = sqrt(vip(6))
!
    vip(7) = a(1)+oset
    do i = 2, ndim
        vip(i+6) = a(i)
    end do
!
!     CALCUL DU NOUVEAU POINT D'EQUILIBRE V10 EN CAS DE SCIAGE
!     LE SCIAGE FAIT DIMINUER L'EPAISSEUR DU JOINT
!     => OSET EST DECROISSANT (SCIAGE)
    vip(10) = oset
!
!     VISUALISATION DES FLUX, DES GRAD DE PRESSION ET DE LA PRESSION
!     DANS LE REPERE GLOBAL
    if (ifhyme) then
!
        gploc(1) = 0.d0
        gploc(2) = gp(1)
        if (ndim .eq. 3) then
            gploc(3) = gp(2)
        end if
!
        fhloc(1) = 0.d0
        fhloc(2) = sigma(ndim+1)
        if (ndim .eq. 3) then
            fhloc(3) = sigma(2*ndim-1)
        end if
!
        call matinv('S', ndim, coorot(ndim+1), invrot, r8bid)
        call pmavec('ZERO', ndim, invrot, gploc, gpglo)
        call pmavec('ZERO', ndim, invrot, fhloc, fhglo)
!
!       CONTRAINTE MECANIQUE NORMALE SANS PRESSION DE FLUIDE CALCULEE
!       ON ANNULE SON INFLUENCE
        vip(11) = sigma(1)+presg
!
        vip(12) = gpglo(1)
        vip(13) = gpglo(2)
        vip(15) = fhglo(1)
        vip(16) = fhglo(2)
!
        if (ndim .eq. 3) then
            vip(14) = gpglo(3)
            vip(17) = fhglo(3)
        else
            vip(14) = 0.d0
            vip(17) = 0.d0
        end if
!
!       PRESSION DE FLUIDE CALCULEE AUX NOEUDS (DDL) ET INTERPOL AU PG
        vip(18) = presg
!
    else
!
!       CONTRAINTE MECANIQUE NORMALE SANS PRESSION DE FLUIDE IMPOSEE
!       ON ANNULE SON INFLUENCE
        vip(11) = sigma(1)+presfl
!
        vip(12) = 0.d0
        vip(13) = 0.d0
        vip(14) = 0.d0
        vip(15) = 0.d0
        vip(16) = 0.d0
        vip(17) = 0.d0
!
!       PRESSION DE FLUIDE IMPOSEE AU PG :
        vip(18) = presfl
!
    end if
!
!
! #####################################
! CALCUL DE LA MATRICE TANGENTE
! #####################################
!
5000 continue
    if (.not. rigi) goto 999
!
!     INITIALISATION DE DSIDEP
    call r8inir(6*6, 0.d0, dsidep, 1)
!
! CALCUL DE LA MATRICE TANGENTE HYDRO
!------------------------------------
!
    if (ifhyme) then
!
!       TERME : DW/DGP  (POUR KTAN P P)
        do n = 1, ndim-1
!
            dsidep(ndim+n, ndim+n) = -rhof*(max(amin, a(1)+amin))**3/(12*visf)
!
        end do
!
!       TERME : DW/DDELTA_N  (POUR KTAN P U)
        do n = 1, ndim-1
!
            if (a(1) .lt. 0.d0) then
                dsidep(ndim+n, 1) = 0.d0
            else
                dsidep(ndim+n, 1) = -3*rhof*gp(n)*(a(1)+amin)**2/(12*visf)
            end if
!
        end do
!
    end if
!
! DSIGMA_N/DDELTA_N
    if (ifouv .eq. 0) then
        dsidep(1, 1) = kn
    else
        dsidep(1, 1) = 0.d0
    end if
! DSIGMA_N/DDELTA_T
    do i = 2, ndim
        dsidep(1, i) = 0.d0
    end do
! DSIGMA_T/DDELTA_N
    do i = 2, ndim
        if ((ifouv .eq. 0) .and. (ifplas .eq. 1)) then
            dsidep(i, 1) = -tau(i)*mu*kn*kt/abstau/(kt+kappa)
        else
            dsidep(i, 1) = 0.d0
        end if
    end do
! DSIGMA_T/DDELTA_T
    if (ifplas .eq. 1) then
        coefhd = -(kappa*lambda+adhe-mu*sigma_eff)*kt**2/abstau**3/(kt+kappa)
        coefd = kappa*kt/(kt+kappa)-coefhd*abstau**2
        do j = 2, ndim
            do i = j, ndim
                if (i .eq. j) then
                    kronec = 1
                else
                    kronec = 0
                end if
                dsidep(j, i) = coefhd*tau(j)*tau(i)+coefd*kronec
                dsidep(i, j) = coefhd*tau(i)*tau(j)+coefd*kronec
            end do
        end do
    else
        do i = 2, ndim
            dsidep(i, i) = kt
        end do
    end if
!
!
!CC MATRICE TANGENTE DE CONTACT OUVERT
!
! DANS LE CAS OU L'ELEMENT EST TOTALEMENT CASSE ON INTRODUIT UNE
! RIGIDITE ARTIFICIELLE DANS LA MATRICE TANGENTE POUR ASSURER
! LA CONVERGENCE
!
    rigart = 1.d-8
!
!     POUR LE JOINT OUVERT LA PARTIE NORMALE EST CORRIGEE
    if (ifouv .eq. 1) then
!       COMPLETEMENT CASSE NORMALE
        dsidep(1, 1) = kn*rigart
    end if
!     POUR LE JOINT SANS ECROUISSAGE LA PARTIE TANGENTIELLE EST DECALEE
    if (kappa .eq. 0.d0) then
!     COMPLETEMENT CASSE TANGENTE
        do i = 2, ndim
            dsidep(i, i) = dsidep(i, i)+kt*rigart
        end do
    end if
!
!
999 continue
!
end subroutine
