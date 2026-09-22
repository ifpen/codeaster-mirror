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
! aslint: disable=W0104, W1306
!
module vmis_isot_nl_module

    use scalar_newton_module, only: &
        newton_state, &
        utnewt

    use tenseur_dime_module, only: &
        proten, &
        kron, &
        deviator, &
        voigt, &
        identity

    implicit none
    private
    public:: CONSTITUTIVE_LAW, Init, InitGradVari, InitViscoPlasticity, Integrate, PathFollowing

#include "asterf_types.h"
#include "asterc/r8gaem.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/pinorm.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/zerop2.h"

! --------------------------------------------------------------------------------------------------

    ! Material characteristics

    type MATERIAL
        real(kind=8) :: lambda, deuxmu, troismu_yng, troisk
        real(kind=8) :: r0, rh, r1, g1, r2, g2, rk, p0, gk
        real(kind=8) :: eps_luders, sig_luders
        real(kind=8) :: c = 0.d0
        real(kind=8) :: r = 0.d0
        real(kind=8) :: v0 = 0.d0
        real(kind=8) :: q = 2.d0
    end type MATERIAL

    ! VMIS_ISOT_NL class
    type CONSTITUTIVE_LAW
        integer(kind=8)       :: exception = 0
        aster_logical :: elas, rigi, vari, pred, pilo
        aster_logical :: grvi = ASTER_FALSE
        aster_logical :: visc = ASTER_FALSE
        aster_logical :: luders = ASTER_FALSE
        aster_logical :: uniax = ASTER_FALSE
        integer(kind=8)       :: ndimsi, itemax
        real(kind=8)  :: cvuser
        real(kind=8)  :: phi = 0.d0
        real(kind=8)  :: telq
        real(kind=8)  :: kam
        type(MATERIAL):: mat
    end type CONSTITUTIVE_LAW

contains

! =====================================================================
!  OBJECT CREATION AND INITIALISATION
! =====================================================================

    function Init(ndimsi, option, fami, kpg, ksp, imate, itemax, precvg) &
        result(self)

        implicit none

        integer(kind=8), intent(in)          :: kpg, ksp, imate, itemax, ndimsi
        real(kind=8), intent(in)    :: precvg
        character(len=16), intent(in):: option
        character(len=*), intent(in) :: fami
        type(CONSTITUTIVE_LAW)      :: self
! --------------------------------------------------------------------------------------------------
! ndimsi    symmetric tensor dimension (2*ndim)
! option    computation option
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! itemax    max number of iterations for the solver
! precvg    required accuracy (with respect to stress level))
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter   :: nbel = 2, nbec = 10
! --------------------------------------------------------------------------------------------------
        integer(kind=8)             :: iok(nbel+nbec)
        real(kind=8)        :: valel(nbel), valec(nbec)
        real(kind=8)        :: r8nan
        character(len=16)   :: nomel(nbel), nomec(nbec)
! --------------------------------------------------------------------------------------------------
        data nomel/'E', 'NU'/
        data nomec/'R0', 'RH', 'R1', 'GAMMA_1', 'R2', 'GAMMA_2', 'RK', 'P0', 'GAMMA_M', &
            'EPSP_LUDERS'/
! --------------------------------------------------------------------------------------------------

        ! Variables non initialisees
        r8nan = r8nnem()
        self%telq = r8nan
        self%kam = r8nan

        ! Parametres generaux
        self%ndimsi = ndimsi
        self%itemax = itemax
        self%cvuser = precvg

        ! 1D ou non
        self%uniax = ndimsi .eq. 1

        ! Options de calcul
        self%elas = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'FULL_MECA_ELAS'
        self%rigi = option .eq. 'RIGI_MECA_TANG' .or. option .eq. 'RIGI_MECA_ELAS' &
                    .or. option .eq. 'FULL_MECA' .or. option .eq. 'FULL_MECA_ELAS'
        self%vari = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'FULL_MECA' &
                    .or. option .eq. 'RAPH_MECA'
        self%pred = option .eq. 'RIGI_MECA_ELAS' .or. option .eq. 'RIGI_MECA_TANG'
        self%pilo = option .eq. 'PILO_PRED_ELAS'

        ! Elasticity material parameters
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ELAS', 0, ' ', [0.d0], nbel, nomel, valel, &
                    iok, 2)
        self%mat%lambda = valel(1)*valel(2)/((1+valel(2))*(1-2*valel(2)))
        self%mat%deuxmu = valel(1)/(1+valel(2))
        self%mat%troisk = valel(1)/(1.d0-2.d0*valel(2))
        self%mat%troismu_yng = merge(valel(1), 1.5d0*self%mat%deuxmu, self%uniax)

        ! Hardening material parameters (with default values)
        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'ECRO_NL', 0, ' ', [0.d0], nbec, nomec, &
                    valec, iok, 0)

        if (iok(1) .ne. 0) call utmess('F', 'COMPOR1_52')
        self%mat%r0 = valec(1)
        self%mat%rh = merge(valec(2), 0.d0, iok(2) .eq. 0)
        self%mat%r1 = merge(valec(3), 0.d0, iok(3) .eq. 0)
        self%mat%g1 = merge(valec(4), 0.d0, iok(4) .eq. 0)
        self%mat%r2 = merge(valec(5), 0.d0, iok(5) .eq. 0)
        self%mat%g2 = merge(valec(6), 0.d0, iok(6) .eq. 0)
        self%mat%rk = merge(valec(7), 0.d0, iok(7) .eq. 0)
        self%mat%p0 = merge(valec(8), 0.d0, iok(8) .eq. 0)
        self%mat%gk = merge(valec(9), 1.d0, iok(9) .eq. 0)
        self%mat%eps_luders = merge(valec(10), -1.d0, iok(10) .eq. 0)

        ! Contrainte de reference pour les calculs de precision / residus
        self%mat%sig_luders = f_ecro_loca(self, max(self%mat%eps_luders, 0.d0))

    end function Init

! =====================================================================
!  COMPLEMENTARY INITIALISATION FOR GRAD_VARI
! =====================================================================

    subroutine InitGradVari(self, fami, kpg, ksp, imate, lag, apg)

        implicit none

        integer(kind=8), intent(in)          :: kpg, ksp, imate
        real(kind=8), intent(in)     :: lag, apg
        character(len=*), intent(in) :: fami
        type(CONSTITUTIVE_LAW), intent(inout)      :: self
! --------------------------------------------------------------------------------------------------
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! lag       Lagrangian value
! apg       nonlocal hardening variable
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter:: nb = 2
! --------------------------------------------------------------------------------------------------
        integer(kind=8)             :: iok(nb)
        real(kind=8)        :: vale(nb)
        character(len=16)   :: nom(nb)
! --------------------------------------------------------------------------------------------------
        data nom/'C_GRAD_VARI', 'PENA_LAGR'/
! --------------------------------------------------------------------------------------------------

        self%grvi = ASTER_TRUE

        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'NON_LOCAL', 0, ' ', [0.d0], nb, nom, vale, &
                    iok, 2)
        self%mat%c = vale(1)
        self%mat%r = vale(2)
        self%phi = lag+self%mat%r*apg

    end subroutine InitGradVari

! =====================================================================
!  COMPLEMENTARY INITIALISATION FOR VISCOPLASTICITY
! =====================================================================

    subroutine InitViscoPlasticity(self, fami, kpg, ksp, imate, deltat)

        implicit none

        integer(kind=8), intent(in)          :: kpg, ksp, imate
        real(kind=8), intent(in)             :: deltat
        character(len=*), intent(in)         :: fami
        type(CONSTITUTIVE_LAW), intent(inout):: self
! --------------------------------------------------------------------------------------------------
! fami      Gauss point set
! kpg       Gauss point number
! ksp       Layer number (for structure elements)
! imate     material pointer
! deltat    time increment (instap - instam)
! --------------------------------------------------------------------------------------------------
        integer(kind=8), parameter:: nb = 2
! --------------------------------------------------------------------------------------------------
        integer(kind=8)             :: iok(nb)
        real(kind=8)        :: vale(nb)
        character(len=16)   :: nom(nb)
! --------------------------------------------------------------------------------------------------
        data nom/'N', 'K'/
! --------------------------------------------------------------------------------------------------

        self%visc = ASTER_TRUE

        call rcvalb(fami, kpg, ksp, '+', imate, ' ', 'NORTON', 0, ' ', [0.d0], nb, nom, vale, &
                    iok, 2)
        self%mat%q = 1.d0/vale(1)
        self%mat%v0 = vale(2)/deltat**self%mat%q

    end subroutine InitViscoPlasticity

! =====================================================================
!  INTEGRATION OF THE CONSTITUTIVE LAW (MAIN ROUTINE)
! =====================================================================

    subroutine Integrate(self, eps, vim, sig, vip, deps_sig, vinl, dphi_sig, deps_vi, dphi_vi)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)         :: eps(:), vim(:)
        real(kind=8), intent(out)         :: sig(:), vip(:), deps_sig(:, :)
        real(kind=8), intent(out), optional:: dphi_sig(:), deps_vi(:), dphi_vi, vinl
! --------------------------------------------------------------------------------------------------
! eps       strain at the end of the time step
! vim       internal variables at the beginning of the time step
! sig       stress at the end of the time step
! vip       internal variables at the end of the time step
! deps_sig  derivee dsig / deps
! vinl      nonlocal internal variable (=vi(1)) (grad_vari)
! dphi_sig  derivee dsig / dphi   (grad_vari)
! deps_vi   derivee dka  / deps  (grad_vari)
! dphi_vi   derivee dka  / dphi  (grad_vari)
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: state
        real(kind=8)    :: kam, ka, epm(self%ndimsi), ep(self%ndimsi), rac2(self%ndimsi)
        real(kind=8)    :: vdum1(self%ndimsi), vdum2(self%ndimsi), rdum
        type(MATERIAL)  :: realMat
! --------------------------------------------------------------------------------------------------

        ASSERT(present(vinl) .eqv. self%grvi)
        ASSERT(present(dphi_sig) .eqv. self%grvi)
        ASSERT(present(deps_vi) .eqv. self%grvi)
        ASSERT(present(dphi_sig) .eqv. self%grvi)

! unpack internal variables
        rac2 = voigt(self%ndimsi)
        kam = vim(1)
        epm = vim(3:2+self%ndimsi)*rac2

! plasticity behaviour integration

        ! Mode plateau de Luders ou non
        self%luders = kam .lt. self%mat%eps_luders
        do
            ! Si plateau Luders -> premier passage en plasticite parfaite
            if (self%luders) then
                realMat = self%mat
                self%mat%r0 = realMat%sig_luders
                self%mat%rh = 0.d0
                self%mat%r1 = 0.d0
                self%mat%g1 = 0.d0
                self%mat%r2 = 0.d0
                self%mat%g2 = 0.d0
                self%mat%rk = 0.d0
                self%mat%p0 = 0.d0
                self%mat%gk = 1.d0
            end if

            if (self%grvi) then
                call ComputePlasticity(self, eps, kam, epm, state, ka, ep, sig, deps_sig, &
                                       dphi_sig, deps_vi, dphi_vi)
            else if (self%uniax) then
                    call ComputePlasticity1D(self, eps(1), kam, epm(1), state, ka, ep(1), sig(1), &
                                         deps_sig(1, 1))
                else
                call ComputePlasticity(self, eps, kam, epm, state, ka, ep, sig, deps_sig, &
                                       vdum1, vdum2, rdum)
            end if
            if (self%exception .ne. 0) goto 999

            ! Juste pour etre sur
            if (.not. self%luders) then
                ASSERT(ka .ge. self%mat%eps_luders)
            end if

            if (.not. self%luders .or. ka .le. self%mat%eps_luders) exit

            ! La solution n'est pas sur le plateau de Luders: 2e passage en ecrouissage standard
            self%mat = realMat
            self%luders = ASTER_FALSE
        end do

! pack internal variables
        if (self%grvi) vinl = ka

        if (self%vari) then
            vip(1) = ka
            vip(2) = state
            vip(3:8) = 0
            vip(3:2+self%ndimsi) = ep/rac2

            ! Complément des composantes plastiques par incompressibilité
            if (self%uniax) vip(4:5) = -vip(3)/2
        end if

999     continue
    end subroutine Integrate

! =====================================================================
!  PLASTICITY COMPUTATION AND TANGENT OPERATORS
! =====================================================================

    subroutine ComputePlasticity(self, eps, kam, epm, state, ka, ep, &
                                 t, deps_t, dphi_t, deps_ka, dphi_ka)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)             :: eps(:)
        real(kind=8), intent(in)             :: kam, epm(:)
        integer(kind=8), intent(out)                 :: state
        real(kind=8), intent(out)            :: ka, ep(:)
        real(kind=8), intent(out)            :: t(:), deps_t(:, :)
        real(kind=8), intent(out)            :: dphi_t(:), deps_ka(:), dphi_ka
! --------------------------------------------------------------------------------------------------
! eps       deformation a la fin du pas de temps
! ka        variable d'ecrouissage kappa (in=debut pas de temps, out=fin)
! state     etat pendant le pas (0=elastique, 1=plastique, 2=singulier) (in=debut, out=fin)
! ep        deformation plastique (in=debut, out=fin)
! t         contrainte en fin de pas de temps
! deps_t    derivee dt / deps
! dphi_t    derivee dt / dphi   (grad_vari)
! deps_ka   derivee dka / deps  (grad_vari)
! dphi_ka   derivee dka / dphi  (grad_vari)
! --------------------------------------------------------------------------------------------------
        aster_logical   :: visc
        integer(kind=8)         :: ite, itev
        real(kind=8)    :: kr(size(eps)), presig
        real(kind=8)    :: mve, dkas, rks, rvs, mvs, dka
        real(kind=8)    :: tel(size(eps)), telh, telq, teld(size(eps)), dep(size(eps))
        real(kind=8)    :: pdev(size(eps), size(eps))
        real(kind=8)    :: equ, dka_equ, dv_equ, rk, rv, dkax, res, rvx, mh
        real(kind=8)    :: n(size(eps)), deps_telq(size(eps)), deps_n(size(eps), size(eps))
        real(kind=8)    :: dtelq_ka, dka_mv
        type(newton_state):: mem
! --------------------------------------------------------------------------------------------------

!   Initialisation
        kr = kron(self%ndimsi)

!   Contrainte elastique
        tel = self%mat%lambda*sum(eps(1:3)-epm(1:3))*kr+self%mat%deuxmu*(eps-epm)
        telh = sum(tel(1:3))/3.d0
        teld = tel-telh*kr
        telq = sqrt(1.5d0*dot_product(teld, teld))

!   Copie des informations dans self pour utilisation des fonctions du module
        self%kam = kam
        self%telq = telq

!   Bornes pour kappa et valeurs correspondantes pour la fonction MV
        dkas = telq/self%mat%troismu_yng
        mve = f_m_hat(self, kam)
        rks = f_ecro(self, kam+dkas)
        rvs = f_visco(self, dkas)
        mvs = -(rks+rvs)

!   Seuil de convergence absolu
        presig = self%mat%sig_luders*self%cvuser

! ======================================================================
!               INTEGRATION DE LA LOI DE COMPORTEMENT
! ======================================================================

! --------------------------------------------------------------------------------------------------
!  REGIME ELASTIQUE
! --------------------------------------------------------------------------------------------------

!   Tir elastique
        if (mve .le. presig) then
            state = 0
            dka = 0.d0
            dep = 0
            goto 800
        end if

! --------------------------------------------------------------------------------------------------
!  REGIME ECOULEMENT SINGULIER
! --------------------------------------------------------------------------------------------------

        if (mvs .le. 0) goto 200

        dka = dkas
        dkax = dkas-(rks+rvs)/(self%mat%r+self%mat%rh)
        visc = ASTER_FALSE

        do ite = 1, self%itemax

            ! Evaluation du terme plastique et du terme visqueux
            if (visc) then
                dka = inv_visco(self, rv)
            else
                rv = f_visco(self, dka)
            end if
            rk = f_ecro(self, kam+dka)

            ! Test de convergence
            res = rk+rv
            if (abs(res) .le. presig) exit

            ! Borne max discriminante plas / visc identifiee
            if (res .gt. 0 .and. rk .le. presig) then
                visc = ASTER_TRUE
                rvx = rv
                itev = ite-1
            end if

            ! Resolution de la plasticite seule (sans viscosite)
            if (.not. visc) then
                equ = rk
                dka_equ = dka_ecro(self, kam+dka)
                dka = utnewt(dka, equ, dka_equ, ite, mem, xmin=dkas, xmax=dkax)

                ! Resolution de la viscoplasticite (avec changement de variable)
            else
                equ = rk+rv
                dv_equ = dka_ecro(self, kam+dka)*dv_inv_visco(self, rv)+1
                rv = utnewt(rv, equ, dv_equ, ite-itev, mem, xmin=rvs, xmax=rvx)
            end if

        end do
        if (ite .gt. self%itemax) then
            self%exception = 1
            goto 999
        end if

        state = 2
        dep = teld/self%mat%deuxmu
        goto 800

! --------------------------------------------------------------------------------------------------
!  REGIME ECOULEMENT REGULIER
! --------------------------------------------------------------------------------------------------

200     continue
        dka = 0.d0
        visc = ASTER_FALSE

        do ite = 1, self%itemax

            ! Evaluation du terme plastique et du terme visqueux
            if (visc) then
                dka = inv_visco(self, rv)
            else
                rv = f_visco(self, dka)
            end if
            mh = f_m_hat(self, kam+dka)

            ! Test de convergence
            res = -mh+rv
            if (abs(res) .le. presig) exit

            ! Borne max discriminante plas / visc identifiee
            if (res .gt. 0 .and. -mh .le. presig) then
                visc = ASTER_TRUE
                rvx = rv
                itev = ite-1
            end if

            ! Resolution de la plasticite seule (sans viscosite)
            if (.not. visc) then
                equ = -mh
                dka_equ = -dka_m_hat(self, kam+dka)
                dka = utnewt(dka, equ, dka_equ, ite, mem, xmin=0.d0, xmax=dkas)

                ! Resolution de la viscoplasticite (avec changement de variable)
            else
                equ = -mh+rv
                dv_equ = -dka_m_hat(self, kam+dka)*dv_inv_visco(self, rv)+1
                rv = utnewt(rv, equ, dv_equ, ite-itev, mem, xmin=0.d0, xmax=rvx)
            end if

        end do
        if (ite .gt. self%itemax) then
            self%exception = 1
            goto 999
        end if

        state = 1
        dep = 1.5d0*dka*teld/telq

800     continue
        ka = kam+dka
        ep = epm+dep
        t = tel-self%mat%deuxmu*dep

! ======================================================================
!                           MATRICES TANGENTES
! ======================================================================

        if (.not. self%rigi) goto 999

        dphi_t = 0
        deps_ka = 0
        dphi_ka = 0

        ! Etat pour choisir l'operateur tangent en phase de prediction
        if (self%pred) then

            ! Regimes frontieres mve=0 avec telq=0 si singulier ou telq>0 sinon
            if (mve .lt. -presig) then
                state = 0
            else if (telq .le. presig) then
                state = 2
            else
                state = 1
            end if

        end if

        ! Regime elastique
        if (state .eq. 0 .or. self%elas) then
            deps_t = self%mat%lambda*proten(kr, kr)+self%mat%deuxmu*identity(self%ndimsi)

            ! Regime plastique regulier
        else if (state .eq. 1) then

            ! Projecteur deviatorique
            pdev = identity(self%ndimsi)-proten(kr, kr)/3.d0

            ! Quantites liees a la contrainte elastique
            n = teld/telq
            deps_telq = self%mat%troismu_yng*n
            deps_n = (self%mat%deuxmu*pdev-proten(n, deps_telq))/telq

            ! Variations de kappa
            dka_mv = dka_m_hat(self, kam+dka)-dka_visco(self, dka)
            dtelq_ka = -dtelq_m_hat(self, kam+dka)/dka_mv
            deps_ka = dtelq_ka*deps_telq
            dphi_ka = -dphi_m_hat(self, kam+dka)/dka_mv

            ! Operateurs tangents
            deps_t = self%mat%lambda*proten(kr, kr)+self%mat%deuxmu*identity(self%ndimsi) &
                     -self%mat%troismu_yng*dka*deps_n &
                     -self%mat%troismu_yng*proten(n, deps_ka)

            dphi_t = -self%mat%troismu_yng*n*dphi_ka

            ! Regime singulier (dphi_t et deps_ka sont nulles)
        else if (state .eq. 2) then
            deps_t = self%mat%troisk/3.d0*proten(kr, kr)
            dphi_ka = -dphi_ecro(self, kam+dka)/(dka_ecro(self, kam+dka)+dka_visco(self, dka))

        else
            ASSERT(ASTER_FALSE)
        end if

999     continue

    end subroutine ComputePlasticity

! =====================================================================
!  1D PLASTICITY COMPUTATION AND TANGENT OPERATORS
! =====================================================================

    subroutine ComputePlasticity1D(self, eps, kam, epm, state, ka, ep, &
                                   t, deps_t)

        implicit none

        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in)             :: eps
        real(kind=8), intent(in)             :: kam, epm
        integer(kind=8), intent(out)         :: state
        real(kind=8), intent(out)            :: ka, ep
        real(kind=8), intent(out)            :: t, deps_t
! --------------------------------------------------------------------------------------------------
! eps       deformation a la fin du pas de temps
! ka        variable d'ecrouissage kappa (in=debut pas de temps, out=fin)
! state     etat pendant le pas (0=elastique, 1=plastique) (in=debut, out=fin)
! ep        deformation plastique (in=debut, out=fin)
! t         contrainte en fin de pas de temps
! deps_t    derivee dt / deps
! --------------------------------------------------------------------------------------------------
        aster_logical   :: visc
        integer(kind=8) :: ite, itev
        real(kind=8)    :: presig
        real(kind=8)    :: mve, dkas, rks, rvs, mvs, dka
        real(kind=8)    :: tel, telq, dep, sg
        real(kind=8)    :: equ, dka_equ, dv_equ, rv, res, rvx, mh
        real(kind=8)    :: deps_telq, dtelq_ka, deps_ka, dka_mv
        type(newton_state):: mem
! --------------------------------------------------------------------------------------------------

!   Contrainte elastique
        tel = self%mat%troismu_yng*(eps-epm)
        telq = abs(tel)
        sg = sign(1.d0, tel)

!   Copie des informations dans self pour utilisation des fonctions du module
        self%kam = kam
        self%telq = telq

!   Bornes pour kappa et valeurs correspondantes pour la fonction MV
        dkas = telq/self%mat%troismu_yng
        mve = f_m_hat(self, kam)
        rks = f_ecro(self, kam+dkas)
        rvs = f_visco(self, dkas)
        mvs = -(rks+rvs)

        ! Pas de grad_vari en 1D
        ASSERT(mvs .le. 0)

!   Seuil de convergence absolu
        presig = self%mat%sig_luders*self%cvuser

! ======================================================================
!               INTEGRATION DE LA LOI DE COMPORTEMENT
! ======================================================================

! --------------------------------------------------------------------------------------------------
!  REGIME ELASTIQUE
! --------------------------------------------------------------------------------------------------

!   Tir elastique
        if (mve .le. presig) then
            state = 0
            dka = 0.d0
            dep = 0
            goto 800
        end if

! --------------------------------------------------------------------------------------------------
!  REGIME ECOULEMENT REGULIER
! --------------------------------------------------------------------------------------------------

        dka = 0.d0
        visc = ASTER_FALSE

        do ite = 1, self%itemax

            ! Evaluation du terme plastique et du terme visqueux
            if (visc) then
                dka = inv_visco(self, rv)
            else
                rv = f_visco(self, dka)
            end if
            mh = f_m_hat(self, kam+dka)

            ! Test de convergence
            res = -mh+rv
            if (abs(res) .le. presig) exit

            ! Borne max discriminante plas / visc identifiee
            if (res .gt. 0 .and. -mh .le. presig) then
                visc = ASTER_TRUE
                rvx = rv
                itev = ite-1
            end if

            ! Resolution de la plasticite seule (sans viscosite)
            if (.not. visc) then
                equ = -mh
                dka_equ = -dka_m_hat(self, kam+dka)
                dka = utnewt(dka, equ, dka_equ, ite, mem, xmin=0.d0, xmax=dkas)

                ! Resolution de la viscoplasticite (avec changement de variable)
            else
                equ = -mh+rv
                dv_equ = -dka_m_hat(self, kam+dka)*dv_inv_visco(self, rv)+1
                rv = utnewt(rv, equ, dv_equ, ite-itev, mem, xmin=0.d0, xmax=rvx)
            end if

        end do
        if (ite .gt. self%itemax) then
            self%exception = 1
            goto 999
        end if

        state = 1
        dep = dka*sg

800     continue
        ka = kam+dka
        ep = epm+dep
        t = tel-self%mat%troismu_yng*dep

! ======================================================================
!                           MATRICES TANGENTES
! ======================================================================

        if (.not. self%rigi) goto 999

        ! Etat pour choisir l'operateur tangent en phase de prediction
        if (self%pred) then

            ! Regimes frontieres mve=0 avec telq=0 si singulier ou telq>0 sinon
            if (mve .lt. -presig) then
                state = 0
            else
                state = 1
            end if

        end if

        ! Regime elastique
        if (state .eq. 0 .or. self%elas) then
            deps_t = self%mat%troismu_yng

            ! Regime plastique
        else
            ! Quantites liees a la contrainte elastique
            deps_telq = sg*self%mat%troismu_yng

            ! Variations de kappa
            dka_mv = dka_m_hat(self, kam+dka)-dka_visco(self, dka)
            dtelq_ka = -dtelq_m_hat(self, kam+dka)/dka_mv
            deps_ka = dtelq_ka*deps_telq

            ! Operateurs tangents
            deps_t = self%mat%troismu_yng*(1.d0-sg*deps_ka)
        end if

999     continue

    end subroutine ComputePlasticity1D

! =====================================================================
!  PATH FOLLOWING
! =====================================================================

    subroutine PathFollowing(self, dtau, vim, eps0, eps1, copilo)

        implicit none
        type(CONSTITUTIVE_LAW), intent(inout):: self
        real(kind=8), intent(in) :: dtau, vim(:), eps0(:), eps1(:)
        real(kind=8), intent(out) :: copilo(:)
! ---------------------------------------------------------------------
! dtau          target path-following increment
! vim           internal variables at t-
! eps0          constant strain
! eps1          path-following strain
! copilo        path-following coefficients
! ---------------------------------------------------------------------
        real(kind=8), parameter:: rac2_3 = sqrt(2.d0/3.d0)
! ---------------------------------------------------------------------
        real(kind=8)    :: kr(self%ndimsi), rac2(self%ndimsi)
        real(kind=8)    :: kam, ka, dka, epm(self%ndimsi)
        real(kind=8):: s0(self%ndimsi), s1(self%ndimsi), rk, gk
! ---------------------------------------------------------------------

        ! Initialisation
        kr = kron(self%ndimsi)
        rac2 = voigt(self%ndimsi)

        ! Normalisation: 5% plasticity increment corresponds to dtau=1
        dka = dtau*5.d-2

        ! unpack internal variables
        kam = vim(1)
        epm = vim(3:2+self%ndimsi)*rac2

        ! elastic stresses: telq = norm(s0+eta*s1)
        s0 = self%mat%troismu_yng*(eps0-epm)
        s1 = self%mat%troismu_yng*eps1
        if (self%ndimsi .ne. 1) then
            s0 = rac2_3*deviator(s0)
            s1 = rac2_3*deviator(s1)
        end if

        ! target threshold
        ka = kam+dka
        rk = f_ecro_loca(self, max(self%mat%eps_luders, ka))
        gk = self%mat%troismu_yng*dka+rk

        ! Path-following coefficients through norm2(eta*s1+s0) = gk
        call pinorm(s1, s0, gk, dtau, copilo)

    end subroutine PathFollowing

! ----------------------------------------------------------------------------------------
!  Liste des fonctions intermediaires et leurs derivees
! ----------------------------------------------------------------------------------------

! --> Derivees: dka_*, dteh_*, dteq_*, dphi_*

    function f_visco(self, dka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::dka
        res = self%mat%v0*dka**self%mat%q
    end function f_visco

    function dka_visco(self, dka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::dka

        if (self%visc .and. self%pred) then
            ! infini -> on renvoie une valeur tres grande
            ! (ok car utilise uniquement pour l'operateur tangent)
            res = r8gaem()
        else
            res = self%mat%v0*self%mat%q*dka**(self%mat%q-1)
        end if

    end function dka_visco

    function inv_visco(self, v) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::v
        res = (v/self%mat%v0)**(1/self%mat%q)
    end function inv_visco

    function dv_inv_visco(self, v) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::v
        res = 1/(self%mat%q*self%mat%v0)*(v/self%mat%v0)**(1/self%mat%q-1)
    end function dv_inv_visco

    function f_ecro_loca(self, ka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::ka
        res = self%mat%r0+self%mat%rh*ka+ &
              self%mat%r1*(1-exp(-self%mat%g1*ka))+ &
              self%mat%r2*(1-exp(-self%mat%g2*ka))+ &
              self%mat%rk*(ka+self%mat%p0)**self%mat%gk
    end function f_ecro_loca

    function dka_ecro_loca(self, ka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::ka

        res = self%mat%rh+ &
              self%mat%r1*self%mat%g1*exp(-self%mat%g1*ka)+ &
              self%mat%r2*self%mat%g2*exp(-self%mat%g2*ka)+ &
              self%mat%rk*self%mat%gk*(ka+self%mat%p0)**(self%mat%gk-1)

    end function dka_ecro_loca

    function f_ecro(self, ka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::ka
        res = f_ecro_loca(self, ka)+self%mat%r*ka-self%phi
    end function f_ecro

    function dka_ecro(self, ka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::ka

        res = dka_ecro_loca(self, ka)+self%mat%r

    end function dka_ecro

    function dphi_ecro(self, ka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::ka
        res = -1.d0
    end function dphi_ecro

    function f_m_hat(self, ka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::ka
        res = self%telq-self%mat%troismu_yng*(ka-self%kam)-f_ecro(self, ka)
    end function f_m_hat

    function dka_m_hat(self, ka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::ka
        res = -self%mat%troismu_yng-dka_ecro(self, ka)
    end function dka_m_hat

    function dtelq_m_hat(self, ka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::ka
        res = 1.d0
    end function dtelq_m_hat

    function dphi_m_hat(self, ka) result(res)
        implicit none
        type(CONSTITUTIVE_LAW), intent(in):: self
        real(kind=8):: res
        real(kind=8), intent(in)::ka
        res = -dphi_ecro(self, ka)
    end function dphi_m_hat

end module vmis_isot_nl_module
