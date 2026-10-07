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
module HHO_precalc_module
!
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/chpchd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/inical.h"
#include "asterfort/megeom.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic
!
! Specific routines for precomputation
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: hhoPreCalcOp, hhoPreCalc, hhoAddInputField, hhoPreCalcBase
!
! --------------------------------------------------------------------------------------------------
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoPreCalc(model)
!
        implicit none
!
        character(len=8), intent(in) :: model
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic
!
! Precomputation of the operators and basis
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
! --- Precompute basis
!
        call hhoPreCalcBase(model)
!
! --- Precompute Operators
!
        call hhoPreCalcOp(model)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoPreCalcOp(model)
!
        implicit none
!
        character(len=8), intent(in) :: model
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic
!
! Precomputation of the operators
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), parameter :: nbin = 2
        integer(kind=8), parameter :: nbout = 2
        character(len=8)  :: lpain(nbin), lpaout(nbout)
        character(len=19) :: lchin(nbin), lchout(nbout)
        character(len=19) :: ligrel_model
        character(len=16) :: option
        character(len=1)  :: base
        character(len=24) :: chgeom
!
! --------------------------------------------------------------------------------------------------
!
        base = 'G'
        option = 'HHO_PRECALC_OP'
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel_model)
!
! --- Init fields
!
        call inical(nbin, lpain, lchin, nbout, lpaout, lchout)
!
! --- Geometry field
!
        call megeom(model, chgeom)
!
! --- Input fields
!
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom(1:19)
        lpain(2) = 'PCHHOBS'
        lchin(2) = model(1:8)//'.HHO.BASE'
!
! --- Output fields
!
        lpaout(1) = 'PCHHOGT'
        lchout(1) = model(1:8)//'.HHO.GRAD'
        lpaout(2) = 'PCHHOST'
        lchout(2) = model(1:8)//'.HHO.STAB'
!
! - Compute
!
        call calcul('S', option, ligrel_model, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, base, &
                    'OUI')
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoPreCalcBase(model)
!
        implicit none
!
        character(len=8), intent(in) :: model
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic
!
! Precomputation of the operators
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), parameter :: nbin = 1
        integer(kind=8), parameter :: nbout = 1
        character(len=8)  :: lpain(nbin), lpaout(nbout)
        character(len=19) :: lchin(nbin), lchout(nbout)
        character(len=19) :: ligrel_model, chbase, chelno
        character(len=16) :: option
        character(len=1)  :: base
        character(len=24) :: chgeom
!
! --------------------------------------------------------------------------------------------------
!
        base = 'G'
        option = 'HHO_PRECALC_BS'
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel_model)
        chelno = model(1:8)//'.HHO.TMP'
        chbase = model(1:8)//'.HHO.BASE'
!
! --- Init fields
!
        call inical(nbin, lpain, lchin, nbout, lpaout, lchout)
!
! --- Geometry field
!
        call megeom(model, chgeom)
!
! --- Input fields
!
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom(1:19)
!
! --- Output fields
!
        lpaout(1) = 'PCHHOBO'
        lchout(1) = chelno
!
! - Compute
!
        call calcul('S', option, ligrel_model, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, base, 'OUI')
!
! - ELNO -> NOEU (average)
!
        call chpchd(chelno, "NOEU", ' ', 'NON', base, chbase)
        call detrsd("CHAM_ELEM", chelno)

!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoAddInputField(model, mxchin, lchin, lpain, nchin)
!
        implicit none
!
        integer(kind=8), intent(in) :: mxchin
        character(len=*), intent(in) :: model
        character(len=*), intent(inout) :: lpain(mxchin)
        character(len=*), intent(inout) :: lchin(mxchin)
        integer(kind=8), intent(inout) :: nchin
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic
!
! Add input field for CALCUL.F90
!
! --------------------------------------------------------------------------------------------------
!
! In  model  : name of model
! In  mxchin : maximum number of input fields
! IO  lpain  : list of parameters
! IO  lchin  : list of fields
! IO  nbin   : number of input fields
!
! --------------------------------------------------------------------------------------------------
        character(len=16) :: repk
!
        call dismoi('EXI_HHO', model(1:8)//".MODELE", 'LIGREL', repk=repk)

        if (repk == "OUI") then
!
! --- Output fields
!
            ASSERT(nchin+3 .le. mxchin)
            lpain(nchin+1) = 'PCHHOGT'
            lchin(nchin+1) = model(1:8)//'.HHO.GRAD'
            lpain(nchin+2) = 'PCHHOST'
            lchin(nchin+2) = model(1:8)//'.HHO.STAB'
            lpain(nchin+3) = 'PCHHOBS'
            lchin(nchin+3) = model(1:8)//'.HHO.BASE'
            nchin = nchin+3
        end if
!
    end subroutine
!
end module
