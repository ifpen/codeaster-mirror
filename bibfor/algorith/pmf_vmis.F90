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

subroutine pmf_vmis(for_pmf, nf, nbvalc, &
                    compor, crit, defam, defap, varim, &
                    varimp, contm, defm, ddefp, modf, &
                    sigf, varip, codret)
!
! --------------------------------------------------------------------------------------------------
!
!               COMPORTEMENT DES ÉLÉMENTS DE POUTRE MULTI-FIBRES
!
! --------------------------------------------------------------------------------------------------
!
!   IN
!       compor  : information sur le comportement du groupe de fibres
!       crit    : critères de convergence locaux
!       nf      : nombre de fibres du groupe
!       nbvalc  : nombre de variable internes
!       defam   : déformations anélastiques a l'instant précédent
!       defap   : déformations anélastiques a l'instant du calcul
!       varim   : variables internes moins
!       varimp  : variables internes itération précédente (pour DE BORST)
!       contm   : contraintes moins par fibre
!       defm    : déformation à l'instant du calcul précédent
!       ddefp   : incrément de déformation
!
!   OUT
!       modf    : module tangent des fibres
!       sigf    : contrainte a l'instant actuel des fibres
!       varip   : variables internes a l'instant actuel
!       codret  : code retour (0 c'est ok)
!
! --------------------------------------------------------------------------------------------------
!
    use pmfcom_type
!
    implicit none

#include "MultiFiber_type.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/nm1dci.h"
#include "asterfort/nm1dis.h"
#include "asterfort/paeldt.h"
#include "asterfort/rcexistvarc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/vmci1d.h"
!
! --------------------------------------------------------------------------------------------------
!
    type(pmfcom_user), intent(in) :: for_pmf
    integer(kind=8)      :: nf, nbvalc, codret
    real(kind=8) :: contm(nf), defm(nf), ddefp(nf), modf(nf), sigf(nf)
    real(kind=8) :: varimp(nbvalc*nf), varip(nbvalc*nf), varim(nbvalc*nf)
    real(kind=8) :: crit(*), defap(*), defam(*)
!
    character(len=24) :: compor(*)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbval = 2
    integer(kind=8)             :: icodre(nbval)
    real(kind=8)        :: valres(nbval)
!
    integer(kind=8)      :: ksp, fib, ivari, nbvari_grfibre
    real(kind=8) :: ep, em, depsth
    real(kind=8) :: depsm, nu
!
    character(len=4)  :: fami
    character(len=8)  :: materi
    character(len=16) :: rela_comp, algo
!
    aster_logical     :: istemp
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)             :: kpg
    integer(kind=8)             :: debsp
    integer(kind=8)             :: icdmat
    real(kind=8)        :: instam
    real(kind=8)        :: instap
    real(kind=8)        :: epsm
    character(len=16)   :: option
!
! --------------------------------------------------------------------------------------------------
    kpg = for_pmf%kpg
    icdmat = for_pmf%icdmat
    option = for_pmf%option
    debsp = for_pmf%debsp
    instam = for_pmf%instam
    instap = for_pmf%instap
    epsm = for_pmf%epsm
! --------------------------------------------------------------------------------------------------
    codret = 0
    fami = 'RIGI'
    materi = compor(MULTI_FIBER_MATER) (1:8)
    rela_comp = compor(MULTI_FIBER_RELA) (1:16)
    algo = compor(MULTI_FIBER_ALGO) (1:16)
    read (compor(MULTI_FIBER_NBVARI), '(I24)') nbvari_grfibre
!   Vérification du nombre de fibre
    ASSERT(nbvari_grfibre .le. nbvalc)
!
! --------------------------------------------------------------------------------------------------
    if (rela_comp .eq. 'VMIS_CINE_GC') then
        ! Température ou pas ?
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call rcvalb(fami, 1, 1, '+', icdmat, materi, 'ELAS', &
                        0, '', [0.d0], 1, ['E'], valres, icodre, 1)
            ep = valres(1)
            em = ep
            ! Boucle sur chaque fibre
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                depsm = ddefp(fib)
                call vmci1d('RIGI', kpg, ksp, icdmat, em, &
                            ep, contm(fib), depsm, varim(ivari), option, &
                            materi, sigf(fib), varip(ivari), modf(fib))
            end do
        else
            ! Boucle sur chaque fibre avec température
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth)
                depsm = ddefp(fib)-depsth
                call vmci1d('RIGI', kpg, ksp, icdmat, em, &
                            ep, contm(fib), depsm, varim(ivari), option, &
                            materi, sigf(fib), varip(ivari), modf(fib))
            end do
        end if
!
! --------------------------------------------------------------------------------------------------
    else if (rela_comp .eq. 'VMIS_CINE_LINE') then
        ! Température ou pas ?
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call rcvalb(fami, 1, 1, '+', icdmat, materi, 'ELAS', &
                        0, '', [0.d0], 1, ['E'], valres, icodre, 1)
            ep = valres(1)
            em = ep
            ! Boucle sur chaque fibre
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                depsm = ddefp(fib)
                call nm1dci('RIGI', kpg, ksp, icdmat, em, &
                            ep, contm(fib), depsm, varim(ivari), option, &
                            materi, sigf(fib), varip(ivari), modf(fib))
            end do
        else
            ! Boucle sur chaque fibre
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth)
                depsm = ddefp(fib)-depsth
                call nm1dci('RIGI', kpg, ksp, icdmat, em, &
                            ep, contm(fib), depsm, varim(ivari), option, &
                            materi, sigf(fib), varip(ivari), modf(fib))
            end do
        end if
!
! --------------------------------------------------------------------------------------------------
    else if ((rela_comp .eq. 'VMIS_ISOT_LINE') .or. &
             (rela_comp .eq. 'VMIS_ISOT_TRAC')) then
        ! Température ou pas ?
        istemp = rcexistvarc('TEMP')
        if (.not. istemp) then
            call rcvalb(fami, 1, 1, '+', icdmat, materi, 'ELAS', &
                        0, '', [0.d0], 1, ['E'], valres, icodre, 1)
            ep = valres(1)
            em = ep
            ! Boucle sur chaque fibre
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                depsm = ddefp(fib)
                call nm1dis('RIGI', kpg, ksp, icdmat, em, &
                            ep, contm(fib), depsm, varim(ivari), option, &
                            rela_comp, materi, sigf(fib), varip(ivari), modf(fib))
            end do
        else
            ! Boucle sur chaque fibre
            do fib = 1, nf
                ivari = nbvalc*(fib-1)+1
                ksp = debsp-1+fib
                call paeldt(kpg, ksp, fami, 'T', icdmat, materi, em, ep, nu, depsth)
                depsm = ddefp(fib)-depsth
                call nm1dis('RIGI', kpg, ksp, icdmat, em, &
                            ep, contm(fib), depsm, varim(ivari), option, &
                            rela_comp, materi, sigf(fib), varip(ivari), modf(fib))
            end do
        end if
    else
        call utmess('F', 'ELEMENTS2_39', sk=rela_comp)
    end if
!
end subroutine
