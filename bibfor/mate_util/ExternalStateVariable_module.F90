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
! ==================================================================================================
!
! Module for External State Variables
!
! ==================================================================================================
!
module ExternalStateVariable_module
! ==================================================================================================
    use ExternalStateVariable_type
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: creaCata, creaJvObjects, creaMaps, debugMaps
    public :: freeCata, fillJvObjects, fillMaps, shrinkMaps
    public :: getAccessToDescriptiveMap, getParametersOnCell
    public :: existInList, extendList
    private :: sameVariables
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/afva01.h"
#include "asterfort/alcart.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cesexi.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/ExternalStateVariable_type.h"
#include "asterfort/getelem.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/lteatt.h"
#include "asterfort/nocart.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/xvarc_temp.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! creaCata
!
! Create catalog of all external state variables
!
! Ptr exteVariCata     : catalog of external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine creaCata(exteVariCata)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(EXTE_VARI_CATA), pointer :: exteVariCata(:)
! ----- Locals
        integer(kind=8) :: iExteStatVari, nbCmp, iCmp
        character(len=8) :: exteVariName, physQuantity
!   ------------------------------------------------------------------------------------------------
!
! ----- Allocate catalog of external state variables
        allocate (exteVariCata(EXTEVARI_NB_MAXI))

! ----- Definition
        do iExteStatVari = 1, EXTEVARI_NB_MAXI
! --------- Main properties
            exteVariName = listExteStatVari(iExteStatVari)
            physQuantity = listPhysQuantity(iExteStatVari)
            nbCmp = listNbCmp(iExteStatVari)
            ASSERT(nbCmp .le. EXTEVARI_NBCMP_MAXI)
            exteVariCata(iExteStatVari)%name = exteVariName
            exteVariCata(iExteStatVari)%physQuantity = physQuantity
            exteVariCata(iExteStatVari)%nbCmp = nbCmp

! --------- List of components
            allocate (exteVariCata(iExteStatVari)%listCmp(nbCmp))
            if (exteVariName .eq. 'TEMP') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "TEMP"
                exteVariCata(iExteStatVari)%listCmp(2)%physQuantityCmp = "TEMP_MIL"
                exteVariCata(iExteStatVari)%listCmp(3)%physQuantityCmp = "TEMP_INF"
                exteVariCata(iExteStatVari)%listCmp(4)%physQuantityCmp = "TEMP_SUP"
                exteVariCata(iExteStatVari)%listCmp(5)%physQuantityCmp = "DTX"
                exteVariCata(iExteStatVari)%listCmp(6)%physQuantityCmp = "DTY"
                exteVariCata(iExteStatVari)%listCmp(7)%physQuantityCmp = "DTZ"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "TEMP"
                exteVariCata(iExteStatVari)%listCmp(2)%nameCmp = "TEMP_MIL"
                exteVariCata(iExteStatVari)%listCmp(3)%nameCmp = "TEMP_INF"
                exteVariCata(iExteStatVari)%listCmp(4)%nameCmp = "TEMP_SUP"
                exteVariCata(iExteStatVari)%listCmp(5)%nameCmp = "DTX"
                exteVariCata(iExteStatVari)%listCmp(6)%nameCmp = "DTY"
                exteVariCata(iExteStatVari)%listCmp(7)%nameCmp = "DTZ"
                ASSERT(nbCmp .eq. 7)
                exteVariCata(iExteStatVari)%fieldType = 'TEMP'

            elseif (exteVariName .eq. 'NEUT1') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "X1"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "NEUT1"
                ASSERT(nbCmp .eq. 1)
                exteVariCata(iExteStatVari)%fieldType = 'NEUT'

            elseif (exteVariName .eq. 'NEUT2') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "X1"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "NEUT2"
                ASSERT(nbCmp .eq. 1)
                exteVariCata(iExteStatVari)%fieldType = 'NEUT'

            elseif (exteVariName .eq. 'NEUT3') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "X1"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "NEUT3"
                ASSERT(nbCmp .eq. 1)
                exteVariCata(iExteStatVari)%fieldType = 'NEUT'

            elseif (exteVariName .eq. 'GEOM') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "X"
                exteVariCata(iExteStatVari)%listCmp(2)%physQuantityCmp = "Y"
                exteVariCata(iExteStatVari)%listCmp(3)%physQuantityCmp = "Z"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "X"
                exteVariCata(iExteStatVari)%listCmp(2)%nameCmp = "Y"
                exteVariCata(iExteStatVari)%listCmp(3)%nameCmp = "Z"
                ASSERT(nbCmp .eq. 3)
                exteVariCata(iExteStatVari)%fieldType = 'GEOM'

            elseif (exteVariName .eq. 'CORR') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "CORR"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "CORR"
                ASSERT(nbCmp .eq. 1)
                exteVariCata(iExteStatVari)%fieldType = 'CORR'

            elseif (exteVariName .eq. 'IRRA') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "IRRA"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "IRRA"
                ASSERT(nbCmp .eq. 1)
                exteVariCata(iExteStatVari)%fieldType = 'IRRA'

            elseif (exteVariName .eq. 'HYDR') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "HYDR"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "HYDR"
                ASSERT(nbCmp .eq. 1)
                exteVariCata(iExteStatVari)%fieldType = 'HYDR_ELNO'

            elseif (exteVariName .eq. 'SECH') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "SECH"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "SECH"
                ASSERT(nbCmp .eq. 1)
                exteVariCata(iExteStatVari)%fieldType = 'SECH'

            elseif (exteVariName .eq. 'PTOT') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "PTOT"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "PTOT"
                ASSERT(nbCmp .eq. 1)
                exteVariCata(iExteStatVari)%fieldType = 'PTOT'

            elseif (exteVariName .eq. 'EPSA') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "EPXX"
                exteVariCata(iExteStatVari)%listCmp(2)%physQuantityCmp = "EPYY"
                exteVariCata(iExteStatVari)%listCmp(3)%physQuantityCmp = "EPZZ"
                exteVariCata(iExteStatVari)%listCmp(4)%physQuantityCmp = "EPXY"
                exteVariCata(iExteStatVari)%listCmp(5)%physQuantityCmp = "EPXZ"
                exteVariCata(iExteStatVari)%listCmp(6)%physQuantityCmp = "EPYZ"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "EPSAXX"
                exteVariCata(iExteStatVari)%listCmp(2)%nameCmp = "EPSAYY"
                exteVariCata(iExteStatVari)%listCmp(3)%nameCmp = "EPSAZZ"
                exteVariCata(iExteStatVari)%listCmp(4)%nameCmp = "EPSAXY"
                exteVariCata(iExteStatVari)%listCmp(5)%nameCmp = "EPSAXZ"
                exteVariCata(iExteStatVari)%listCmp(6)%nameCmp = "EPSAYZ"
                ASSERT(nbCmp .eq. 6)
                exteVariCata(iExteStatVari)%fieldType = 'EPSA'

            elseif (exteVariName .eq. 'M_ACIER') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "V1"
                exteVariCata(iExteStatVari)%listCmp(2)%physQuantityCmp = "V2"
                exteVariCata(iExteStatVari)%listCmp(3)%physQuantityCmp = "V3"
                exteVariCata(iExteStatVari)%listCmp(4)%physQuantityCmp = "V4"
                exteVariCata(iExteStatVari)%listCmp(5)%physQuantityCmp = "V5"
                exteVariCata(iExteStatVari)%listCmp(6)%physQuantityCmp = "V6"
                exteVariCata(iExteStatVari)%listCmp(7)%physQuantityCmp = "V7"
                exteVariCata(iExteStatVari)%listCmp(8)%physQuantityCmp = "V8"
                exteVariCata(iExteStatVari)%listCmp(9)%physQuantityCmp = "V9"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "PFERRITE"
                exteVariCata(iExteStatVari)%listCmp(2)%nameCmp = "PPERLITE"
                exteVariCata(iExteStatVari)%listCmp(3)%nameCmp = "PBAINITE"
                exteVariCata(iExteStatVari)%listCmp(4)%nameCmp = "PMARTENS"
                exteVariCata(iExteStatVari)%listCmp(5)%nameCmp = "PAUSTENI"
                exteVariCata(iExteStatVari)%listCmp(6)%nameCmp = "PCOLDSUM"
                exteVariCata(iExteStatVari)%listCmp(7)%nameCmp = "TAUSTE"
                exteVariCata(iExteStatVari)%listCmp(8)%nameCmp = "TRANSF"
                exteVariCata(iExteStatVari)%listCmp(9)%nameCmp = "TACIER"
                ASSERT(nbCmp .eq. 9)
                exteVariCata(iExteStatVari)%fieldType = 'META_ELNO'

            elseif (exteVariName .eq. 'M_ZIRC') then
                exteVariCata(iExteStatVari)%listCmp(1)%physQuantityCmp = "V1"
                exteVariCata(iExteStatVari)%listCmp(2)%physQuantityCmp = "V2"
                exteVariCata(iExteStatVari)%listCmp(3)%physQuantityCmp = "V3"
                exteVariCata(iExteStatVari)%listCmp(4)%physQuantityCmp = "V4"
                exteVariCata(iExteStatVari)%listCmp(5)%physQuantityCmp = "V5"
                exteVariCata(iExteStatVari)%listCmp(1)%nameCmp = "ALPHPUR"
                exteVariCata(iExteStatVari)%listCmp(2)%nameCmp = "ALPHBETA"
                exteVariCata(iExteStatVari)%listCmp(3)%nameCmp = "BETA"
                exteVariCata(iExteStatVari)%listCmp(4)%nameCmp = "TZIRC"
                exteVariCata(iExteStatVari)%listCmp(5)%nameCmp = "TEMPS"
                ASSERT(nbCmp .eq. 5)
                exteVariCata(iExteStatVari)%fieldType = 'META_ELNO'
            else
                ASSERT(ASTER_FALSE)
            end if
        end do

! ----- Debug if required
        if (EXTEVARI_DBG_READ .eq. 1) then
            write (6, *) 'Description de toutes les variables de commande'
            write (6, *) '==============================================='
            write (6, *) 'Nombre total de variables de commande disponibles:', EXTEVARI_NB_MAXI
            do iExteStatVari = 1, EXTEVARI_NB_MAXI
                write (6, *) 'Variable de commande :', iExteStatVari
                write (6, *) '> Nom           :', exteVariCata(iExteStatVari)%name
                write (6, *) '> GRANDEUR      :', exteVariCata(iExteStatVari)%physQuantity
                write (6, *) '> Type du champ :', exteVariCata(iExteStatVari)%fieldType
                write (6, *) '> NB_CMP        :', exteVariCata(iExteStatVari)%nbCmp
                do iCmp = 1, exteVariCata(iExteStatVari)%nbCmp
                    write (6, *) '> Nombre de composantes :', iCmp
                    write (6, *) '>> CMP_GD   :', &
                        exteVariCata(iExteStatVari)%listCmp(iCmp)%physQuantityCmp
                    write (6, *) '>> CMP_VARC :', &
                        exteVariCata(iExteStatVari)%listCmp(iCmp)%nameCmp
                end do
            end do
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! creaJvObjects
!
! Create JEVEUX objects in material field
!
! In  jvBase           : JEVEUX base where to create objects
! In  mateField        : name of material field (CHAM_MATER)
! In  exteVariAffe     : datastructure for assigned external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine creaJvObjects(jvBase, mateField, exteVariAffe)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: jvBase
        character(len=8), intent(in) :: mateField
        type(EXTE_VARI_AFFE), intent(in) :: exteVariAffe
! ----- Locals
        integer(kind=8) :: jvDummy, nbCmpTotal
        character(len=24) :: jvNameCVRCNom, jvNameCVRCVarc, jvNameCVRCGd, jvNameCVRCCmp
!   ------------------------------------------------------------------------------------------------
!
        nbCmpTotal = exteVariAffe%nbCmpTotal
        jvNameCVRCNom = mateField//'.CVRCNOM'
        jvNameCVRCVarc = mateField//'.CVRCVARC'
        jvNameCVRCGd = mateField//'.CVRCGD'
        jvNameCVRCCmp = mateField//'.CVRCCMP'
        call wkvect(jvNameCVRCNom, jvBase//' V K8', nbCmpTotal, jvDummy)
        call wkvect(jvNameCVRCVarc, jvBase//' V K8', nbCmpTotal, jvDummy)
        call wkvect(jvNameCVRCGd, jvBase//' V K8', nbCmpTotal, jvDummy)
        call wkvect(jvNameCVRCCmp, jvBase//' V K8', nbCmpTotal, jvDummy)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! creaMaps
!
! Create maps objects in material field
!
! Two maps for all external state variable (TEMP, SECH, EPSA, etc..)
!    First maps: VALE_REFE
!    Second maps: other parameters
!
! In  jvBase           : JEVEUX base where to create objects
! In  mesh             : name of mesh
! In  mateField        : name of material field (CHAM_MATER)
! Ptr exteVariCata     : catalog of external state variables
! In  exteVariAffe     : datastructure for assigned external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine creaMaps(jvBase, mesh, mateField, exteVariCata, exteVariAffe)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: jvBase
        character(len=8), intent(in) :: mesh, mateField
        type(EXTE_VARI_CATA), pointer :: exteVariCata(:)
        type(EXTE_VARI_AFFE), intent(in) :: exteVariAffe
! ----- Locals
        integer(kind=8) :: iAffe, nbAffe, cataIndx, iret, iCmp
        character(len=8) :: exteVariName, cmpName
        character(len=19) :: exteVariMap1, exteVariMap2
        character(len=8), pointer :: mapNcmp(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        nbAffe = exteVariAffe%nbAffe
        ASSERT(nbAffe .ne. 0)
        do iAffe = 1, nbAffe
            cataIndx = exteVariAffe%exteVariList(iAffe)%cataIndx
            exteVariName = exteVariCata(cataIndx)%name

! --------- Create first map: for each external state variable, VALE_REFE
            exteVariMap1 = mateField//'.'//exteVariName//'.1'
            call exisd('CARTE', exteVariMap1, iret)
            if (iret .eq. 0) then
                call alcart(jvBase, exteVariMap1, mesh, 'NEUT_R')
            end if
            call jeveuo(exteVariMap1//'.NCMP', 'E', vk8=mapNcmp)
            cmpName = 'X'
            do iCmp = 1, EXTEVARI_NBCMP_MAXI
                call codent(iCmp, 'G', cmpName(2:8))
                mapNcmp(iCmp) = cmpName
            end do

! --------- Create second map: for each external state variable, other parameters
            exteVariMap2 = mateField//'.'//exteVariName//'.2'
            call exisd('CARTE', exteVariMap2, iret)
            if (iret .eq. 0) then
                call alcart(jvBase, exteVariMap2, mesh, 'NEUT_K16')
            end if
            call jeveuo(exteVariMap2//'.NCMP', 'E', vk8=mapNcmp)
            cmpName = 'Z'
            do iCmp = 1, EXTEVARI_MAP2_SIZE
                call codent(iCmp, 'G', cmpName(2:8))
                mapNcmp(iCmp) = cmpName
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! freeCata
!
! Free catalog of all external state variables
!
! Ptr exteVariCata     : catalog of external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine freeCata(exteVariCata)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(EXTE_VARI_CATA), pointer :: exteVariCata(:)
! ----- Locals
        integer(kind=8) :: iExteStatVari
!   ------------------------------------------------------------------------------------------------
!
        do iExteStatVari = 1, EXTEVARI_NB_MAXI
            deallocate (exteVariCata(iExteStatVari)%listCmp)
        end do
        deallocate (exteVariCata)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! fillJvObjects
!
! Fill JEVEUX objects in material field
!
! In  mateField        : name of material field (CHAM_MATER)
! Ptr exteVariCata     : catalog of external state variables
! In  exteVariAffe     : datastructure for assigned external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine fillJvObjects(mateField, exteVariCata, exteVariAffe)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: mateField
        type(EXTE_VARI_CATA), pointer :: exteVariCata(:)
        type(EXTE_VARI_AFFE), intent(in) :: exteVariAffe
! ----- Locals
        character(len=24) :: jvNameCVRCNom, jvNameCVRCVarc, jvNameCVRCGd, jvNameCVRCCmp
        character(len=8), pointer :: CVRCVarc(:) => null(), CVRCNom(:) => null()
        character(len=8), pointer :: CVRCGd(:) => null(), CVRCCmp(:) => null()
        integer(kind=8) :: cataIndx, indxCmp
        integer(kind=8) :: iCmp, iAffe
        integer(kind=8) :: nbCmp, nbCmpTotal, nbAffe
        character(len=8) :: exteVariName, physQuantity
        aster_logical, pointer :: v_active(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        nbAffe = exteVariAffe%nbAffe
        nbCmpTotal = exteVariAffe%nbCmpTotal
        ASSERT(nbAffe .ne. 0)

! ----- Access to objects
        jvNameCVRCNom = mateField//'.CVRCNOM'
        jvNameCVRCVarc = mateField//'.CVRCVARC'
        jvNameCVRCGd = mateField//'.CVRCGD'
        jvNameCVRCCmp = mateField//'.CVRCCMP'
        call jeveuo(jvNameCVRCNom, 'E', vk8=CVRCNom)
        call jeveuo(jvNameCVRCVarc, 'E', vk8=CVRCVarc)
        call jeveuo(jvNameCVRCGd, 'E', vk8=CVRCGd)
        call jeveuo(jvNameCVRCCmp, 'E', vk8=CVRCCmp)

! ----- Set values
        AS_ALLOCATE(vl=v_active, size=EXTEVARI_NB_MAXI)
        v_active = ASTER_FALSE
        indxCmp = 0
        do iAffe = 1, nbAffe
            exteVariName = exteVariAffe%exteVariList(iAffe)%exteVariName
            cataIndx = exteVariAffe%exteVariList(iAffe)%cataIndx
            physQuantity = exteVariCata(cataIndx)%physQuantity
            nbCmp = exteVariCata(cataIndx)%nbCmp
            if (.not. v_active(cataIndx)) then
                do iCmp = 1, nbCmp
                    CVRCNom(iCmp+indxCmp) = &
                        exteVariCata(cataIndx)%listCmp(iCmp)%nameCmp
                    CVRCCmp(iCmp+indxCmp) = &
                        exteVariCata(cataIndx)%listCmp(iCmp)%physQuantityCmp
                    CVRCVarc(iCmp+indxCmp) = exteVariName
                    CVRCGd(iCmp+indxCmp) = physQuantity
                end do
                indxCmp = indxCmp+nbCmp
                v_active(cataIndx) = ASTER_TRUE
            end if
        end do
        ASSERT(indxCmp .eq. nbCmpTotal)
        AS_DEALLOCATE(vl=v_active)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! fillMaps
!
! Fill maps in material field
!
! In  mateField        : name of material field (CHAM_MATER)
! In  mesh             : name of mesh
! In  model            : name of model
! Ptr exteVariCata     : catalog of external state variables
! In  exteVariAffe     : datastructure for assigned external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine fillMaps(mateField, mesh, model, &
                        exteVariCata, exteVariAffe)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: mateField, mesh, model
        type(EXTE_VARI_CATA), pointer :: exteVariCata(:)
        type(EXTE_VARI_AFFE), intent(in) :: exteVariAffe
! ----- Locals
        character(len=16), parameter :: factorKeyword = "AFFE_VARC"
        character(len=19) :: exteVariMap1, exteVariMap2
        integer(kind=8) :: cataIndx
        integer(kind=8) :: iCmp, iAffe
        integer(kind=8) :: nbCmp, nbAffe, nbCellAffe, iCellAffe
        integer(kind=8) :: jvCell
        character(len=8) :: exteVariName, affeType, answer
        real(kind=8), pointer :: map1Valv(:) => null()
        character(len=16), pointer :: map2Valv(:) => null()
        real(kind=8) :: valeRefe
        character(len=8) :: funcResult, dsUser
        character(len=19) :: ligrel
        character(len=16) :: fieldType, funcExtrLeft, funcExtrRight
        character(len=24), parameter :: listCell = '&&AFVARC.LIST_CELL'
        aster_logical :: onAllCells, hasTHM
        character(len=16) :: elemTypeName
        integer(kind=8) :: elemTypeNume, cellNume
        integer(kind=8), pointer :: cellAffectedByModel(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        nbAffe = exteVariAffe%nbAffe
        do iAffe = 1, nbAffe
! --------- Parameters from user
            exteVariName = exteVariAffe%exteVariList(iAffe)%exteVariName
            valeRefe = exteVariAffe%exteVariList(iAffe)%valeRefe
            affeType = exteVariAffe%exteVariList(iAffe)%affeType
            fieldType = exteVariAffe%exteVariList(iAffe)%fieldType
            dsUser = exteVariAffe%exteVariList(iAffe)%dsUser
            funcResult = exteVariAffe%exteVariList(iAffe)%funcResult
            funcExtrLeft = exteVariAffe%exteVariList(iAffe)%funcExtrLeft
            funcExtrRight = exteVariAffe%exteVariList(iAffe)%funcExtrRight

! --------- Parameters from catalog
            cataIndx = exteVariAffe%exteVariList(iAffe)%cataIndx
            nbCmp = exteVariCata(cataIndx)%nbCmp

! --------- Set value on map 1
            ASSERT(nbCmp .le. EXTEVARI_NBCMP_MAXI)
            exteVariMap1 = mateField//'.'//exteVariName//'.1'
            call jeveuo(exteVariMap1//'.VALV', 'E', vr=map1Valv)
            do iCmp = 1, nbCmp
                map1Valv(iCmp) = valeRefe
            end do

! --------- Set values on map 2
            exteVariMap2 = mateField//'.'//exteVariName//'.2'
            call jeveuo(exteVariMap2//'.VALV', 'E', vk16=map2Valv)
            map2Valv(1) = exteVariName
            if (affeType .eq. "CHAMP") then
                map2Valv(2) = "CHAMP"
                map2Valv(3) = dsUser
                map2Valv(4) = " "
                map2Valv(5) = " "
                map2Valv(6) = " "
                map2Valv(7) = " "
            else if (affeType .eq. "EVOL") then
                map2Valv(2) = "EVOL"
                map2Valv(3) = dsUser
                map2Valv(4) = fieldType
                map2Valv(5) = funcExtrLeft
                map2Valv(6) = funcExtrRight
                map2Valv(7) = funcResult
            else
                ASSERT(ASTER_FALSE)
            end if

! --------- Topological assignment
            nbCellAffe = 0
            onAllCells = ASTER_FALSE
            call getelem(mesh, factorKeyword, iAffe, ' ', listCell, nbCellAffe, &
                         model=model, onAllCells_=onAllCells)
            if (nbCellAffe .gt. 0) then
                if (onAllCells) then
                    call nocart(exteVariMap1, 1, nbCmp)
                    call nocart(exteVariMap2, 1, EXTEVARI_MAP2_SIZE)
                else
                    call jeveuo(listCell, 'L', jvCell)
                    call nocart(exteVariMap1, 3, nbCmp, &
                                mode='NUM', nma=nbCellAffe, limanu=zi(jvCell))
                    call nocart(exteVariMap2, 3, EXTEVARI_MAP2_SIZE, &
                                mode='NUM', nma=nbCellAffe, limanu=zi(jvCell))
                end if

! --------- Some checks
                if (exteVariName .eq. "TEMP") then
! ------------- For XFEM: change TEMP to TEMP_ELGA
                    call xvarc_temp(affeType, dsUser, funcExtrLeft, funcExtrRight, &
                                    funcResult, nbAffe, exteVariMap2)

! ------------- For THM: no temperature
                    if (model .ne. " ") then
                        call dismoi('EXI_THM', model, 'MODELE', repk=answer)
                        hasTHM = answer .eq. "OUI"
                        if (hasTHM) then
                            call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel)
                            call jeveuo(ligrel//'.TYFE', 'L', vi=cellAffectedByModel)
                            if (onAllCells) then
                                call utmess('F', 'VARC1_1')
                            else
                                if (nbCellAffe .ne. 0) then
                                    do iCellAffe = 1, nbCellAffe
                                        cellNume = zi(jvCell-1+iCellAffe)
                                        elemTypeNume = cellAffectedByModel(cellNume)
                                        call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), &
                                                    elemTypeName)
                                        if (lteatt('TYPMOD2', 'THM', typel=elemTypeName)) then
                                            call utmess('F', 'VARC1_1')
                                        end if
                                    end do
                                end if

                            end if
                        end if
                    end if
                end if
                call jedetr(listCell)
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! shrinkMaps
!
! Shrink number of components to save memory
!
! In  mateField        : name of material field (CHAM_MATER)
! In  exteVariAffe     : datastructure for assigned external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine shrinkMaps(mateField, exteVariAffe)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: mateField
        type(EXTE_VARI_AFFE), intent(in) :: exteVariAffe
! ----- Locals
        character(len=24) :: jvNameCVRCNom, jvNameCVRCVarc, jvNameCVRCGd, jvNameCVRCCmp
        character(len=8), pointer :: CVRCVarc(:) => null(), CVRCNom(:) => null()
        character(len=8), pointer :: CVRCGd(:) => null(), CVRCCmp(:) => null()
        integer(kind=8) :: nbAffe, nbCmp, nbCmpTotal, nbCmpToDelete, physNbCmpMaxi
        integer(kind=8) :: iCmpNew, iCmp
        integer(kind=8), pointer :: cmpToDelete(:) => null()
        character(len=8) :: exteVariName, affeType, dsUser, nameCmp
        character(len=16) :: fieldType
        character(len=19) :: exteVariMap2
        character(len=16), pointer :: map2Vale(:) => null()
        integer(kind=8), pointer :: map2Desc(:) => null()
        aster_logical :: lCmpToDelete, l_other
!   ------------------------------------------------------------------------------------------------
!
        nbAffe = exteVariAffe%nbAffe
        nbCmpTotal = exteVariAffe%nbCmpTotal
        ASSERT(nbAffe .ne. 0)

! ----- Access to objects
        jvNameCVRCNom = mateField//'.CVRCNOM'
        jvNameCVRCVarc = mateField//'.CVRCVARC'
        jvNameCVRCGd = mateField//'.CVRCGD'
        jvNameCVRCCmp = mateField//'.CVRCCMP'
        call jeveuo(jvNameCVRCNom, 'E', vk8=CVRCNom)
        call jeveuo(jvNameCVRCVarc, 'E', vk8=CVRCVarc)
        call jeveuo(jvNameCVRCGd, 'E', vk8=CVRCGd)
        call jeveuo(jvNameCVRCCmp, 'E', vk8=CVRCCmp)

! ----- Detect components in TEMP that are not TEMP (TEMp_MIL, TEMP_SUP, etc.)
        AS_ALLOCATE(vi=cmpToDelete, size=nbCmpTotal)
        nbCmpToDelete = 0
        do iCmp = 1, nbCmpTotal
            exteVariName = CVRCVarc(iCmp)
            nameCmp = CVRCNom(iCmp)
            if (exteVariName .eq. 'TEMP') then
                if (nameCmp .ne. 'TEMP') then
                    nbCmpToDelete = nbCmpToDelete+1
                    cmpToDelete(iCmp) = 1
                end if
            end if
        end do

        if (nbCmpToDelete .ne. 0) then
            exteVariName = "TEMP"

! --------- Access to map
            exteVariMap2 = mateField//'.'//exteVariName//'.2'
            call jeveuo(exteVariMap2//'.DESC', 'E', vi=map2Desc)
            call jeveuo(exteVariMap2//'.VALE', 'E', vk16=map2Vale)
            nbCmp = map2Desc(3)

! --------- Get number of components for this physical quantity
            call jelira(jexnom('&CATA.GD.NOMCMP', 'NEUT_K16'), 'LONMAX', physNbCmpMaxi)

! --------- Detect when TEMP and LAGR are not used
            lCmpToDelete = ASTER_TRUE
            do iCmp = 1, nbCmp
                exteVariName = map2Vale(physNbCmpMaxi*(iCmp-1)+1) (1:8)
                ASSERT(exteVariName .eq. 'TEMP')
                affeType = map2Vale(physNbCmpMaxi*(iCmp-1)+2) (1:8)
                dsUser = map2Vale(physNbCmpMaxi*(iCmp-1)+3) (1:8)
                fieldType = map2Vale(physNbCmpMaxi*(iCmp-1)+4)
                call afva01(affeType, dsUser, fieldType, l_other)
                if (l_other) then
                    lCmpToDelete = ASTER_FALSE
                    exit
                end if
            end do

! --------- ON SUPPRIME CE QUI NE SERT A RIEN
            if (lCmpToDelete) then
                iCmpNew = 0
                do iCmp = 1, nbCmpTotal
                    if (cmpToDelete(iCmp) .eq. 0) then
                        iCmpNew = iCmpNew+1
                        CVRCNom(iCmpNew) = CVRCNom(iCmp)
                        CVRCVarc(iCmpNew) = CVRCVarc(iCmp)
                        CVRCGd(iCmpNew) = CVRCGd(iCmp)
                        CVRCCmp(iCmpNew) = CVRCCmp(iCmp)
                    end if
                end do
                call juveca(jvNameCVRCNom, iCmpNew)
                call juveca(jvNameCVRCVarc, iCmpNew)
                call juveca(jvNameCVRCGd, iCmpNew)
                call juveca(jvNameCVRCCmp, iCmpNew)
                call jeveuo(jvNameCVRCNom, 'E', vk8=CVRCNom)
                call jeveuo(jvNameCVRCVarc, 'E', vk8=CVRCVarc)
                call jeveuo(jvNameCVRCGd, 'E', vk8=CVRCGd)
                call jeveuo(jvNameCVRCCmp, 'E', vk8=CVRCCmp)
                ASSERT(iCmpNew .eq. nbCmpTotal-nbCmpToDelete)
            end if
        end if
        AS_DEALLOCATE(vi=cmpToDelete)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getAccessToDescriptiveMap
!
! Get access to descriptive map for current external state variable
!
! --------------------------------------------------------------------------------------------------
    subroutine getAccessToDescriptiveMap(mateField, exteVariName, map2S, &
                                         map2SJvCESD, map2SJvCESL, map2SCESV, map2SCESK)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: mateField, exteVariName
        character(len=19), intent(in) :: map2S
        integer(kind=8), intent(out) :: map2SJvCESD, map2SJvCESL
        character(len=16), pointer :: map2SCESV(:)
        character(len=8), pointer :: map2SCESK(:)
! ----- Locals
        character(len=19) :: exteVariMap2
        integer(kind=8) :: iret
!   ------------------------------------------------------------------------------------------------
!
        exteVariMap2 = mateField//'.'//exteVariName//'.2'
        call carces(exteVariMap2, 'ELEM', ' ', 'V', map2S, 'A', iret)
        ASSERT(iret .eq. 0)
        call jeveuo(map2S//'.CESD', 'L', map2SJvCESD)
        call jeveuo(map2S//'.CESL', 'L', map2SJvCESL)
        call jeveuo(map2S//'.CESV', 'L', vk16=map2SCESV)
        call jeveuo(map2S//'.CESK', 'L', vk8=map2SCESK)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getParametersOnCell
!
! Get parameters on cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getParametersOnCell(map2SJvCESD, map2SJvCESL, map2SCESV, numeCell, &
                                   exteVariName, exteVariDesc, noValue)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: map2SJvCESD, map2SJvCESL
        character(len=16), pointer :: map2SCESV(:)
        integer(kind=8), intent(in) :: numeCell
        character(len=8), intent(in) :: exteVariName
        type(EXTE_VARI_DESC), intent(out) :: exteVariDesc
        aster_logical, intent(out) :: noValue
! ----- Locals
        integer(kind=8) :: jvAdrs
!   ------------------------------------------------------------------------------------------------
!
        call cesexi('C', map2SJvCESD, map2SJvCESL, numeCell, 1, 1, 1, jvAdrs)

        if (jvAdrs .gt. 0) then
            noValue = ASTER_FALSE
            exteVariDesc%exteVariName = exteVariName
            exteVariDesc%affeType = map2SCESV(jvAdrs-1+2) (1:8)
            exteVariDesc%dsUser = map2SCESV(jvAdrs-1+3) (1:8)
            exteVariDesc%fieldType = map2SCESV(jvAdrs-1+4)
            exteVariDesc%funcExtrLeft = map2SCESV(jvAdrs-1+5) (1:8)
            exteVariDesc%funcExtrRight = map2SCESV(jvAdrs-1+6) (1:8)
            exteVariDesc%funcResult = map2SCESV(jvAdrs-1+7) (1:8)
        else
            noValue = ASTER_TRUE
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! sameVariables
!
! Test to check is two external state variables are the same (except VALE_REFE)
!
! In  exteVariDesc1    : descriptor of first external state variable
! In  exteVariDesc2    : descriptor of second external state variable
! Out same             : .TRUE. if the same
!
! --------------------------------------------------------------------------------------------------
    subroutine sameVariables(exteVariDesc1, exteVariDesc2, same)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(EXTE_VARI_DESC), intent(in) :: exteVariDesc1, exteVariDesc2
        aster_logical, intent(out) :: same
!   ------------------------------------------------------------------------------------------------
!
        same = ASTER_TRUE
        if ((exteVariDesc1%exteVariName .ne. exteVariDesc2%exteVariName) .or. &
            (exteVariDesc1%affeType .ne. exteVariDesc2%affeType) .or. &
            (exteVariDesc1%fieldType .ne. exteVariDesc2%fieldType) .or. &
            (exteVariDesc1%dsUser .ne. exteVariDesc2%dsUser) .or. &
            (exteVariDesc1%funcExtrLeft .ne. exteVariDesc2%funcExtrLeft) .or. &
            (exteVariDesc1%funcExtrRight .ne. exteVariDesc2%funcExtrRight) .or. &
            (exteVariDesc1%funcResult .ne. exteVariDesc2%funcResult)) then
            same = ASTER_FALSE
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! existInList
!
! Test if external state variables already exists in list
!
! In  exteVariDesc     : descriptor of external state variable
! Ptr exteVariList     : list of external state variable
! In  listSize         : size of list
! Out exists           : .TRUE. if exists
!
! --------------------------------------------------------------------------------------------------
    subroutine existInList(exteVariDesc, exteVariList, listSize, exists)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(EXTE_VARI_DESC), intent(in) :: exteVariDesc
        type(EXTE_VARI_DESC), pointer :: exteVariList(:)
        integer(kind=8), intent(in) :: listSize
        aster_logical, intent(out) :: exists
! ----- Locals
        integer(kind=8) :: iVari
        aster_logical :: same
!   ------------------------------------------------------------------------------------------------
!
        exists = ASTER_FALSE
        do iVari = 1, listSize
            call sameVariables(exteVariDesc, exteVariList(iVari), same)
            if (same) then
                exists = ASTER_TRUE
                exit
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! extendList
!
! Extend list of external state variables
!
! Ptr exteVariList     : list of external state variablef
! In  listSize         : size of list
! Out newListSize      : new size of list
!
! --------------------------------------------------------------------------------------------------
    subroutine extendList(exteVariList, listSize, newListSize)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(EXTE_VARI_DESC), pointer :: exteVariList(:)
        integer(kind=8), intent(in) :: listSize
        integer(kind=8), intent(out) :: newListSize
! ----- Locals
        type(EXTE_VARI_DESC), pointer :: exteVariCopy(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        newListSize = 2*listSize
        allocate (exteVariCopy(listSize))
        exteVariCopy(:) = exteVariList(:)
        deallocate (exteVariList)
        allocate (exteVariList(newListSize))
        exteVariList(:listSize) = exteVariCopy(:listSize)
        deallocate (exteVariCopy)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! debugMaps
!
! Print maps for debug
!
! In  mateField        : name of material field (CHAM_MATER)
!
! --------------------------------------------------------------------------------------------------
    subroutine debugMaps(mateField)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: mateField
! ----- Locals
        character(len=24) :: jvNameCVRCNom, jvNameCVRCVarc, jvNameCVRCGd, jvNameCVRCCmp
        character(len=8), pointer :: CVRCVarc(:) => null(), CVRCNom(:) => null()
        character(len=8), pointer :: CVRCGd(:) => null(), CVRCCmp(:) => null()
        integer(kind=8) :: nbCmpTotal, iCmpTotal
!   ------------------------------------------------------------------------------------------------
!

! ----- Access to objects
        jvNameCVRCNom = mateField//'.CVRCNOM'
        jvNameCVRCVarc = mateField//'.CVRCVARC'
        jvNameCVRCGd = mateField//'.CVRCGD'
        jvNameCVRCCmp = mateField//'.CVRCCMP'
        call jeveuo(jvNameCVRCNom, 'L', vk8=CVRCNom)
        call jeveuo(jvNameCVRCVarc, 'L', vk8=CVRCVarc)
        call jeveuo(jvNameCVRCGd, 'L', vk8=CVRCGd)
        call jeveuo(jvNameCVRCCmp, 'L', vk8=CVRCCmp)

! ----- Main parameters
        call jelira(jvNameCVRCVarc, 'LONMAX', nbCmpTotal)

! ----- Print
        write (6, *) 'Contenu des maps pour les variables de commande'
        write (6, *) '==============================================='
        write (6, *) 'Nombre total de composantes:', nbCmpTotal
        do iCmpTotal = 1, nbCmpTotal
            write (6, *) 'Variable de commande :', iCmpTotal
            write (6, *) '> Nom                :', CVRCVarc(iCmpTotal)
            write (6, *) '> Physical quantity  :', CVRCGd(iCmpTotal)
            write (6, *) '> Name of component  :', CVRCNom(iCmpTotal)
            write (6, *) '> Type of component  :', CVRCCmp(iCmpTotal)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module ExternalStateVariable_module
