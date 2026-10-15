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
! Module for kinematic loads
!
! ==================================================================================================
!
module LoadKinematic_module
! ==================================================================================================
    use KineListRela_module
    use KineLoadUtility_module
    use KineListRela_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: kineLoadGlueMeshMeca, kineLoadGlueMeshTher
    public :: kineLoadLinkProj
    private :: kineLoadMeshProjVoVo, kineLoadMeshProjShVo, kineLoadMeshProjShSh
    private :: kineLoadMeshProjVoSh, kineLoadMeshLinkVoSh
    private :: kineLoadMeshLinkVoVo, kineLoadMeshLinkShVo, kineLoadMeshLinkShSh
    private :: kineLoadGlueMeshMecaPara, kineLoadGlueMeshMecaLine
    private :: kineLoadGlueMeshTherPara, kineLoadGlueMeshTherLine
    private :: kineLoadLinkProjPara
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/base3n.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/elrfvf.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/pj2dco.h"
#include "asterfort/pj3dco.h"
#include "asterfort/pj3dcoSupInf.h"
#include "asterfort/pj4dco.h"
#include "asterfort/utmess.h"
#include "MeshTypes_type.h"

contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! kineLoadGlueMeshMeca
!
! Main subroutine to glue mesh for mechanic (LIAISON_MAIL)
!
! In  loadName         : name of load
! In  model            : model
! In  valeType         : affected value type (real, complex or function)
! In  lVerbose         : flag to advanced debug print
! Out listLineRela     : name of datastructure for list of linear relation
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGlueMeshMeca(loadNameZ, modelZ, valeTypeZ, lVerbose, listLineRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: loadNameZ, modelZ, valeTypeZ
        aster_logical, intent(in) :: lVerbose
        character(len=19), intent(out) :: listLineRela
! - Local
        character(len=16), parameter :: factorKeyword = 'LIAISON_MAIL'
        character(len=8), parameter :: meshDebugJv = '&&CALIRC'
        character(len=24) :: cellMastJv, geomMastJv
        character(len=24) :: nodeSlavJv, geomSlavJv
        character(len=8) :: mesh, chnorm, slavTransf(3), mastTransf(3)
        character(len=16) :: linkType
        character(len=4) :: valeType
        character(len=16), parameter :: corres = '&&CALIRC.CORRES'
        character(len=16), parameter :: corre1 = '&&CALIRC.CORRE1'
        character(len=16), parameter :: corre2 = '&&CALIRC.CORRE2'
        character(len=16), parameter :: corre3 = '&&CALIRC.LISV1'
        integer(kind=8) :: geomDime, meshNbNode
        integer(kind=8) :: iOcc, nocc
        aster_logical :: useNormal, lElimMult, lApplyRota, lMastTransf, lSlavTransf, useDisp
        real(kind=8) :: rotaMatr(3, 3)
        aster_logical :: lDistMaxi
        real(kind=8) :: distMaxi, distAlarm, thickness
        integer(kind=8) :: nbCellMast, nbNodeSlav, nbCellSlav
        character(len=8), pointer :: dofName(:) => null()
        integer(kind=8), pointer :: cellMast(:) => null()
        integer(kind=8), pointer :: nodeSlav(:) => null(), cellSlav(:) => null()
        integer(kind=8), pointer :: nodeElim(:) => null(), nodeSkinToBody(:) => null()
        real(kind=8), pointer :: normSlav(:) => null()
!   NOMBRE MAX DE TERMES D'UNE RELATION LINEAIRE EN 3D = 2*27 + 3 = 57 and 165 for COQUE_MASSIF
        integer(kind=8), parameter :: nbTermMaxi = 165
        type(KINE_LIST_RELA) :: kineListRela
!   ------------------------------------------------------------------------------------------------
!
        valeType = valeTypeZ
        ASSERT(valeType .eq. 'REEL')
        call dismoi('NOM_MAILLA', modelZ, 'MODELE', repk=mesh)
        call dismoi('DIM_GEOM', modelZ, 'MODELE', repi=geomDime)
        if (.not. (geomDime .eq. 2 .or. geomDime .eq. 3)) then
            call utmess('F', 'CHARGES7_6')
        end if
        call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=meshNbNode)

! - Name of datastructure for list of linear relations
        listLineRela = '&&CALIRC.RLLISTE'

! - Working vectors
        AS_ALLOCATE(vi=nodeElim, size=meshNbNode)

! - Create object for list of linear relations
        call kineListRelaCreate('Implicit', nbTermMaxi, listLineRela, kineListRela)

        if (lVerbose) then
            call copisd('MAILLAGE', 'V', mesh, meshDebugJv)
        end if

        call getfac(factorKeyword, nocc)
        do iOcc = 1, nocc
! ----- Get main parameters from user
            call kineLoadGlueMeshMecaPara(factorKeyword, iOcc, geomDime, &
                                          linkType, useNormal, lElimMult, &
                                          lMastTransf, mastTransf, &
                                          lSlavTransf, slavTransf, &
                                          lDistMaxi, distMaxi, distAlarm, &
                                          chnorm, thickness, &
                                          useDisp, dofName)

! ----- Get list of master cells
            call kineLoadGetMasterCells(modelZ, mesh, &
                                        factorKeyword, iOcc, &
                                        cellMastJv, nbCellMast, cellMast)

! ----- Get list of slave nodes
            if (useNormal) then
                call kineLoadGetSlaveNodesFromSkin(modelZ, mesh, meshNbNode, geomDime, &
                                                   factorKeyword, iOcc, &
                                                   nodeSlavJv, nbNodeSlav, nodeSlav, &
                                                   nbCellSlav, cellSlav, &
                                                   nodeSkinToBody)
            else
                call kineLoadGetSlaveNodes(mesh, &
                                           factorKeyword, iOcc, &
                                           nodeSlavJv, nbNodeSlav, nodeSlav)
            end if

! ----- Compute normals at slave nodes
            if (useNormal) then
                call kineLoadNormFromSkinCells(mesh, geomDime, &
                                               nbCellSlav, cellSlav, &
                                               nbNodeSlav, nodeSlav, &
                                               normSlav)
            end if

! ----- Manage multiple nodes
            if (lElimMult) then
                call kineLoadElimMult(iOcc, nodeSlavJv, nbNodeSlav, nodeSlav, nodeElim)
            end if

! ----- Read geometric transformation
            call kineLoadReadTransf(mesh, geomDime, &
                                    factorKeyword, iOcc, &
                                    nodeSlavJv, nbNodeSlav, geomSlavJv, &
                                    lApplyRota, rotaMatr)
            if (lApplyRota) then
                ASSERT(linkType .eq. 'MASSIF')
            end if
            geomMastJv = ' '
            if (linkType .eq. 'COQUE') then
                call kineLoadApplyTransf(loadNameZ, mesh, meshNbNode, &
                                         factorKeyword, iOcc, &
                                         lMastTransf, mastTransf, &
                                         lSlavTransf, slavTransf, &
                                         nbNodeSlav, nodeSlav, &
                                         lVerbose, meshDebugJv, &
                                         geomSlavJv, geomMastJv)
            end if

! ----- Compute projection of mesh (PROJ_CHAMP)
            if (linkType .eq. 'MASSIF') then
                call kineLoadMeshProjVoVo(modelZ, geomDime, meshNbNode, &
                                          lDistMaxi, distMaxi, distAlarm, &
                                          nbNodeSlav, geomSlavJv, nodeSlav, &
                                          nbCellMast, cellMast, &
                                          corres)

            else if (linkType .eq. 'COQUE_MASSIF') then
                call kineLoadMeshProjShVo(modelZ, meshNbNode, &
                                          lDistMaxi, distMaxi, distAlarm, &
                                          chnorm, thickness, &
                                          nbNodeSlav, geomSlavJv, nodeSlav, &
                                          nbCellMast, cellMast, &
                                          corres, corre1, corre2, corre3)

            elseif (linkType .eq. 'COQUE') then
                call kineLoadMeshProjShSh(modelZ, meshNbNode, &
                                          lDistMaxi, distMaxi, distAlarm, &
                                          nbNodeSlav, geomSlavJv, nodeSlav, &
                                          nbCellMast, geomMastJv, cellMast, &
                                          corres)

            elseif (linkType .eq. 'MASSIF_COQUE') then
                call kineLoadMeshProjVoSh(modelZ, meshNbNode, &
                                          lDistMaxi, distMaxi, distAlarm, &
                                          nbNodeSlav, geomSlavJv, nodeSlav, &
                                          nbCellMast, geomMastJv, cellMast, &
                                          corres)

            else
                ASSERT(ASTER_FALSE)
            end if

! ----- Apply linear relations
            call kineLoadGlueMeshMecaLine(linkType, mesh, geomDime, &
                                          useDisp, dofName, &
                                          lApplyRota, rotaMatr, &
                                          useNormal, normSlav, nodeSkinToBody, &
                                          corres, corre1, corre2, corre3, &
                                          kineListRela)

! ----- Clean
            call detrsd('CORRESP_2_MAILLA', corres)
            call detrsd('CORRESP_2_MAILLA', corre1)
            call detrsd('CORRESP_2_MAILLA', corre2)
            call jedetr(geomMastJv)
            call jedetr(geomSlavJv)
            call jedetr(cellMastJv)
            call jedetr(nodeSlavJv)
            call jedetr(corre3)
            AS_DEALLOCATE(vi=nodeSkinToBody)
            AS_DEALLOCATE(vk8=dofName)
        end do

! - Clean
        call jedetr(meshDebugJv)
        AS_DEALLOCATE(vi=nodeElim)
        call kineListRelaDelete(kineListRela)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGlueMeshTher
!
! Main subroutine to glue mesh for thermic (LIAISON_MAIL)
!
! In  model            : model
! In  valeType         : affected value type (real, complex or function)
! Out listLineRela     : name of datastructure for list of linear relation
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGlueMeshTher(modelZ, valeTypeZ, listLineRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: modelZ, valeTypeZ
        character(len=19), intent(out) :: listLineRela
! - Local
        character(len=16), parameter :: factorKeyword = 'LIAISON_MAIL'
        character(len=24) :: cellMastJv
        character(len=24) :: nodeSlavJv, geomSlavJv
        character(len=8) :: mesh
        character(len=4) :: valeType
        character(len=16), parameter :: corres = '&&CALIRC.CORRES'
        integer(kind=8) :: geomDime, meshNbNode
        integer(kind=8) :: iOcc, nocc
        aster_logical :: lElimMult, lDistMaxi
        real(kind=8) :: distMaxi, distAlarm
        integer(kind=8) :: nbCellMast, nbNodeSlav
        integer(kind=8), pointer :: cellMast(:) => null(), nodeSlav(:) => null()
        integer(kind=8), pointer :: nodeElim(:) => null()
!   NOMBRE MAX DE TERMES D'UNE RELATION LINEAIRE = 2*27 + 3 = 57
        integer(kind=8), parameter :: nbTermMaxi = 2*MT_NNOMAX+3
        type(KINE_LIST_RELA) :: kineListRela
!   ------------------------------------------------------------------------------------------------
!
        valeType = valeTypeZ
        ASSERT(valeType .eq. 'REEL')
        call dismoi('NOM_MAILLA', modelZ, 'MODELE', repk=mesh)
        call dismoi('DIM_GEOM', modelZ, 'MODELE', repi=geomDime)
        if (.not. (geomDime .eq. 2 .or. geomDime .eq. 3)) then
            call utmess('F', 'CHARGES7_6')
        end if
        call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=meshNbNode)

! - Name of datastructure for list of linear relations
        listLineRela = '&&CALIRC.RLLISTE'

! - Working vectors
        AS_ALLOCATE(vi=nodeElim, size=meshNbNode)

! - Create object for list of linear relations
        call kineListRelaCreate('Implicit', nbTermMaxi, listLineRela, kineListRela)

        call getfac(factorKeyword, nocc)
        do iOcc = 1, nocc
! ----- Get main parameters from user
            call kineLoadGlueMeshTherPara(factorKeyword, iOcc, &
                                          lElimMult, lDistMaxi, distMaxi, distAlarm)

! ----- Get list of master cells
            call kineLoadGetMasterCells(modelZ, mesh, &
                                        factorKeyword, iOcc, &
                                        cellMastJv, nbCellMast, cellMast)

! ----- Get list of slave nodes
            call kineLoadGetSlaveNodes(mesh, &
                                       factorKeyword, iOcc, &
                                       nodeSlavJv, nbNodeSlav, nodeSlav)

! ----- Manage multiple nodes
            if (lElimMult) then
                call kineLoadElimMult(iOcc, nodeSlavJv, nbNodeSlav, nodeSlav, nodeElim)
            end if

! ----- Apply geometric transformation
            call kineLoadReadTransf(mesh, geomDime, &
                                    factorKeyword, iOcc, &
                                    nodeSlavJv, nbNodeSlav, geomSlavJv)

! ----- Compute projection of mesh (PROJ_CHAMP)
            call kineLoadMeshProjVoVo(modelZ, geomDime, meshNbNode, &
                                      lDistMaxi, distMaxi, distAlarm, &
                                      nbNodeSlav, geomSlavJv, nodeSlav, &
                                      nbCellMast, cellMast, &
                                      corres)

! ----- Apply linear relations
            call kineLoadGlueMeshTherLine(mesh, corres, kineListRela)

! ----- Clean
            call detrsd('CORRESP_2_MAILLA', corres)
            call jedetr(geomSlavJv)
            call jedetr(cellMastJv)
            call jedetr(nodeSlavJv)
        end do

! - Clean
        AS_DEALLOCATE(vi=nodeElim)
        call kineListRelaDelete(kineListRela)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGlueMeshMecaPara
!
! Get parameters to glue meshes (mechanic)
!
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! In  geomDime         : geometric dimension (2 or 3)
! Out linkType         : type of link (MASSIF, COQUE, COQUE_MASSIF, MASSIF_COQUE)
! Out useNormal        : flag to link with normal to slave side
! Out lElimMult        : flag to eliminate multiple linear relations
! Out lMastTransf      : apply geometric transformation to master side for COQUE
! Out mastTransf       : geometric transformation to master side for COQUE
! Out lSlavTransf      : apply geometric transformation to slave side for COQUE
! Out slavTransf       : geometric transformation to slave side for COQUE
! Out lDistMaxi        : flag to set maximum value of distance for projection
! Out distMaxi         : maximum value of distance for projection
! Out distAlarm        : value to emit an alarm of distance for projection
! Out chnorm           : field of normals for COQUE_MASSIF
! Out thickness        : thickness of shells for COQUE_MASSIF
! Out useDisp          : flag to use only "pure" displacement dof DX/DY/DZ/DRX/DRY/DRZ
! Ptr dofName          : name of dof to link
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGlueMeshMecaPara(factorKeywordZ, iOcc, geomDime, &
                                        linkType, useNormal, lElimMult, &
                                        lMastTransf, mastTransf, &
                                        lSlavTransf, slavTransf, &
                                        lDistMaxi, distMaxi, distAlarm, &
                                        chnorm, thickness, &
                                        useDisp, dofName)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: factorKeywordZ
        integer(kind=8), intent(in) :: iOcc, geomDime
        character(len=16), intent(out) :: linkType
        aster_logical, intent(out) :: useNormal, lElimMult
        aster_logical, intent(out) :: lMastTransf, lSlavTransf
        character(len=8), intent(out) :: mastTransf(3), slavTransf(3)
        aster_logical, intent(out) :: lDistMaxi
        real(kind=8), intent(out) :: distMaxi, distAlarm
        character(len=8), intent(out) :: chnorm
        real(kind=8), intent(out) :: thickness
        aster_logical, intent(out) :: useDisp
        character(len=8), pointer :: dofName(:)
! - Local
        character(len=8) :: answer
        integer(kind=8) :: nbRet, nbDof, iDof
!   ------------------------------------------------------------------------------------------------
!
        call getvtx(factorKeywordZ, 'TYPE_RACCORD', iocc=iocc, scal=linkType, nbret=nbRet)
        if (geomDime .eq. 2) then
            ASSERT(linkType .eq. 'MASSIF')
        end if

        useDisp = ASTER_FALSE
        useNormal = ASTER_FALSE
        call getvtx(factorKeywordZ, 'DDL', iocc=iocc, nbval=0, nbret=nbDof)
        nbDof = -nbDof
        if (nbDof .eq. 0) then
            useDisp = ASTER_TRUE
            useNormal = ASTER_FALSE
        else
            AS_ALLOCATE(vk8=dofName, size=nbDof)
            call getvtx(factorKeywordZ, 'DDL', iocc=iocc, nbval=nbDof, vect=dofName)
            if (nbDof .eq. 1) then
                if (dofName(1) .eq. 'DNOR') then
                    useNormal = ASTER_TRUE
                end if
            else
                do iDof = 1, nbDof
                    if (dofName(iDof) .eq. 'DNOR') then
                        call utmess('F', 'CHARGES7_2')
                    end if
                end do
            end if
        end if

        call getvtx(factorKeywordZ, 'ELIM_MULT', iocc=iocc, scal=answer, nbret=nbRet)
        lElimMult = answer .eq. 'OUI'
        distMaxi = 0.d0
        distAlarm = 0.d0
        call getvr8(factorKeywordZ, 'DISTANCE_MAX', iocc=iocc, scal=distMaxi, nbret=nbRet)
        lDistMaxi = nbRet .eq. 1
        call getvr8(factorKeywordZ, 'DISTANCE_ALARME', iocc=iocc, scal=distAlarm, nbret=nbRet)
        if (nbRet .eq. 0) then
            distAlarm = -1.d0
        end if

! - Get parameters for linkType = 'COQUE'
        lMastTransf = ASTER_FALSE
        mastTransf = ' '
        lSlavTransf = ASTER_FALSE
        slavTransf = ' '
        if (linkType .eq. 'COQUE') then
            call getvid(factorKeywordZ, 'TRANSF_GEOM_MAIT', iocc=iocc, nbval=3, &
                        vect=mastTransf, nbret=nbRet)
            lMastTransf = nbret .gt. 0
            call getvid(factorKeywordZ, 'TRANSF_GEOM_ESCL', iocc=iocc, nbval=3, &
                        vect=slavTransf, nbret=nbRet)
            lSlavTransf = nbret .gt. 0
        end if

! - Get parameters for linkType = 'COQUE_MASSIF'
        chnorm = ' '
        thickness = 0.d0
        if (linkType .eq. 'COQUE_MASSIF') then
            call getvid(factorKeywordZ, 'CHAM_NORMALE', iocc=iocc, scal=chnorm)
            call getvr8(factorKeywordZ, 'EPAIS', iocc=iocc, scal=thickness)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGlueMeshTherPara
!
! Get parameters to glue meshes (thermic)
!
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! Out lElimMult        : flag to eliminate multiple linear relations
! Out lDistMaxi        : flag to set maximum value of distance for projection
! Out distMaxi         : maximum value of distance for projection
! Out distAlarm        : value to emit an alarm of distance for projection
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGlueMeshTherPara(factorKeywordZ, iOcc, &
                                        lElimMult, &
                                        lDistMaxi, distMaxi, distAlarm)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
        aster_logical, intent(out) ::  lElimMult
        aster_logical, intent(out) :: lDistMaxi
        real(kind=8), intent(out) :: distMaxi, distAlarm
! - Local
        character(len=8) :: answer
        integer(kind=8) :: nbRet
!   ------------------------------------------------------------------------------------------------
!
        call getvtx(factorKeywordZ, 'ELIM_MULT', iocc=iocc, scal=answer, nbret=nbRet)
        lElimMult = answer .eq. 'OUI'
        distMaxi = 0.d0
        distAlarm = 0.d0
        call getvr8(factorKeywordZ, 'DISTANCE_MAX', iocc=iocc, scal=distMaxi, nbret=nbRet)
        lDistMaxi = nbRet .eq. 1
        call getvr8(factorKeywordZ, 'DISTANCE_ALARME', iocc=iocc, scal=distAlarm, nbret=nbRet)
        if (nbRet .eq. 0) then
            distAlarm = -1.d0
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadMeshProjVoVo
!
! Calculates the correspondence slave/master for 'MASSIF'
!
! In  model            : model
! In  geomDime         : geometric dimension (2 or 3)
! In  meshNbNode       : total number of nodes of mesh
! In  lDistMaxi        : flag to set maximum value of distance for projection
! In  distMaxi         : maximum value of distance for projection
! In  distAlarm        : value to emit an alarm of distance for projection
! In  nbNodeSlav       : number of slave nodes
! In  geomSlavJv       : name of JEVEUX object for coordinates of slave nodes after transformation
! Ptr nodeSlav         : pointer to list of slave nodes
! In  nbCellMast       : number of master cells
! Ptr cellMast         : pointer to list of master cells
! In  corres           : name of object for correspondence slave/master
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadMeshProjVoVo(modelZ, geomDime, meshNbNode, &
                                    lDistMaxi, distMaxi, distAlarm, &
                                    nbNodeSlav, geomSlavJv, nodeSlav, &
                                    nbCellMast, cellMast, &
                                    corres)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: modelZ
        integer(kind=8), intent(in) :: geomDime, meshNbNode
        aster_logical, intent(inout) :: lDistMaxi
        real(kind=8), intent(in) :: distMaxi, distAlarm
        integer(kind=8), intent(in) :: nbNodeSlav, nbCellMast
        character(len=24), intent(in) :: geomSlavJv
        integer(kind=8), pointer :: nodeSlav(:), cellMast(:)
        character(len=16), intent(in) :: corres
! - Local
        character(len=8) :: model
        integer(kind=8) :: nbLink
!   ------------------------------------------------------------------------------------------------
!
        model = modelZ
        if (geomDime .eq. 2) then
            call pj2dco('PARTIE', model, model, nbCellMast, cellMast, &
                        nbNodeSlav, nodeSlav, ' ', geomSlavJv, corres, &
                        lDistMaxi, distMaxi, distAlarm)
        else if (geomDime .eq. 3) then
            call pj3dco('PARTIE', model, model, nbCellMast, cellMast, &
                        nbNodeSlav, nodeSlav, ' ', geomSlavJv, corres, &
                        lDistMaxi, distMaxi, distAlarm)
        else
            ASSERT(ASTER_FALSE)
        end if
        call jelira(corres//'.PJEF_NB', 'LONMAX', nbLink)
        ASSERT(nbLink .eq. meshNbNode)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadMeshProjShSh
!
! Calculates the correspondence slave/master for 'COQUE'
!
! In  model            : model
! In  meshNbNode       : total number of nodes of mesh
! In  lDistMaxi        : flag to set maximum value of distance for projection
! In  distMaxi         : maximum value of distance for projection
! In  distAlarm        : value to emit an alarm of distance for projection
! In  nbNodeSlav       : number of slave nodes
! In  geomSlavJv       : name of JEVEUX object for coordinates of slave nodes after transformation
! Ptr nodeSlav         : pointer to list of slave nodes
! In  nbCellMast       : number of master cells
! In  geomMastJv       : name of JEVEUX object for coordinates of master nodes after transformation
! Ptr cellMast         : pointer to list of master cells
! In  corres           : name of JEVEUX object for correspondence slave/master
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadMeshProjShSh(modelZ, meshNbNode, &
                                    lDistMaxi, distMaxi, distAlarm, &
                                    nbNodeSlav, geomSlavJv, nodeSlav, &
                                    nbCellMast, geomMastJv, cellMast, &
                                    corres)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: modelZ
        integer(kind=8), intent(in) :: meshNbNode
        aster_logical, intent(inout) :: lDistMaxi
        real(kind=8), intent(in) :: distMaxi, distAlarm
        integer(kind=8), intent(in) :: nbNodeSlav, nbCellMast
        character(len=24), intent(in) :: geomSlavJv, geomMastJv
        integer(kind=8), pointer :: nodeSlav(:), cellMast(:)
        character(len=16), intent(in) :: corres
! - Local
        character(len=8) :: model
        integer(kind=8) :: nbLink
!   ------------------------------------------------------------------------------------------------
!
        model = modelZ
        call pj4dco('PARTIE', model, model, nbCellMast, cellMast, &
                    nbNodeSlav, nodeSlav, geomMastJv, geomSlavJv, corres, &
                    lDistMaxi, distMaxi, distAlarm)
        call jelira(corres//'.PJEF_NB', 'LONMAX', nbLink)
        ASSERT(nbLink .eq. meshNbNode)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadMeshProjShVo
!
! Calculates the correspondence slave/master for 'COQUE_MASSIF'
!
! In  model            : model
! In  meshNbNode       : total number of nodes of mesh
! In  lDistMaxi        : flag to set maximum value of distance for projection
! In  distMaxi         : maximum value of distance for projection
! In  distAlarm        : value to emit an alarm of distance for projection
! In  chnorm           : field of normals for COQUE_MASSIF
! In  thickness        : thickness of shells for COQUE_MASSIF
! In  nbNodeSlav       : number of slave nodes
! In  geomSlavJv       : name of JEVEUX object for coordinates of slave nodes after transformation
! Ptr nodeSlav         : pointer to list of slave nodes
! In  nbCellMast       : number of master cells
! Ptr cellMast         : pointer to list of master cells
! In  corres           : name of JEVEUX object for correspondence slave/master
! In  corre1           : name of JEVEUX object for correspondence slave/master (superior face)
! In  corre2           : name of JEVEUX object for correspondence slave/master (inferior face)
! In  corre3           : name of JEVEUX object for normal thickness
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadMeshProjShVo(modelZ, meshNbNode, &
                                    lDistMaxi, distMaxi, distAlarm, &
                                    chnormZ, thickness, &
                                    nbNodeSlav, geomSlavJv, nodeSlav, &
                                    nbCellMast, cellMast, &
                                    corres, corre1, corre2, corre3)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: modelZ
        integer(kind=8), intent(in):: meshNbNode
        aster_logical, intent(inout) :: lDistMaxi
        real(kind=8), intent(in) :: distMaxi, distAlarm
        character(len=*), intent(in) :: chnormZ
        real(kind=8), intent(in) :: thickness
        integer(kind=8), intent(in) :: nbNodeSlav, nbCellMast
        character(len=24), intent(in) :: geomSlavJv
        integer(kind=8), pointer :: nodeSlav(:), cellMast(:)
        character(len=16), intent(in) :: corres, corre1, corre2, corre3
! - Local
        character(len=8) :: model
        integer(kind=8) :: nbLink
!   ------------------------------------------------------------------------------------------------
!
        model = modelZ

! - Compute geometric projection (3D <=> Shell) on medium plane
        call pj3dco('PARTIE', &
                    model, model, &
                    nbCellMast, cellMast, &
                    nbNodeSlav, nodeSlav, &
                    ' ', geomSlavJv, corres, &
                    lDistMaxi, distMaxi, distAlarm)

! - Compute geometric projection (3D <=> Shell) on superior and inferior planes
        call pj3dcoSupInf(model, meshNbNode, &
                          nbCellMast, cellMast, &
                          nbNodeSlav, nodeSlav, geomSlavJv, &
                          corre1, corre2, corre3, &
                          chnormZ, thickness)

        call jelira(corres//'.PJEF_NB', 'LONMAX', nbLink)
        ASSERT(nbLink .eq. meshNbNode)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadMeshProjVoSh
!
! Calculates the correspondence slave/master for 'MASSIF_COQUE'
!
! In  model            : model
! In  meshNbNode       : total number of nodes of mesh
! In  lDistMaxi        : flag to set maximum value of distance for projection
! In  distMaxi         : maximum value of distance for projection
! In  distAlarm        : value to emit an alarm of distance for projection
! In  nbNodeSlav       : number of slave nodes
! In  geomSlavJv       : name of JEVEUX object for coordinates of slave nodes after transformation
! Ptr nodeSlav         : pointer to list of slave nodes
! In  nbCellMast       : number of master cells
! In  geomMastJv       : name of JEVEUX object for coordinates of master nodes after transformation
! Ptr cellMast         : pointer to list of master cells
! In  corres           : name of JEVEUX object for correspondence slave/master
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadMeshProjVoSh(modelZ, meshNbNode, &
                                    lDistMaxi, distMaxi, distAlarm, &
                                    nbNodeSlav, geomSlavJv, nodeSlav, &
                                    nbCellMast, geomMastJv, cellMast, &
                                    corres)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: modelZ
        integer(kind=8), intent(in) :: meshNbNode
        aster_logical, intent(inout) :: lDistMaxi
        real(kind=8), intent(in) :: distMaxi, distAlarm
        integer(kind=8), intent(in) :: nbNodeSlav, nbCellMast
        character(len=24), intent(in) :: geomSlavJv, geomMastJv
        integer(kind=8), pointer :: nodeSlav(:), cellMast(:)
        character(len=16), intent(in) :: corres
! - Local
        character(len=8) :: model
        integer(kind=8) :: nbLink
!   ------------------------------------------------------------------------------------------------
!
        model = modelZ
        call pj4dco('PARTIE', model, model, nbCellMast, cellMast, &
                    nbNodeSlav, nodeSlav, geomMastJv, geomSlavJv, corres, &
                    lDistMaxi, distMaxi, distAlarm)
        call jelira(corres//'.PJEF_NB', 'LONMAX', nbLink)
        ASSERT(nbLink .eq. meshNbNode)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGlueMeshMecaLine
!
! Create list of linear relations to glue meshes (mechanic)
!
! In  linkType         : type of link (MASSIF, COQUE, COQUE_MASSIF, MASSIF_COQUE)
! In  mesh             : mesh
! In  geomDime         : geometric dimension (2 or 3)
! In  useDisp          : flag to use only "pure" displacement dof DX/DY/DZ/DRX/DRY/DRZ
! Ptr dofName          : name of dof to link
! In  lApplyRota       : flag when rotation is applied
! In  rotaMatr         : rotation matrix
! In  useNormal        : flag to link with normal to slave side
! Ptr normSlav         : pointer to list of slave normals
! Ptr nodeSkinToBody   : pointer between body and skin nodes
! In  corres           : name of JEVEUX object for correspondence slave/master
! In  corre1           : name of JEVEUX object for correspondence slave/master (superior face)
! In  corre2           : name of JEVEUX object for correspondence slave/master (inferior face)
! In  corre3           : name of JEVEUX object for normal thickness
! IO  kineListRela     : object for list of linear relations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGlueMeshMecaLine(linkType, meshZ, geomDime, &
                                        useDisp, dofName, &
                                        lApplyRota, rotaMatr, &
                                        useNormal, normSlav, nodeSkinToBody, &
                                        corres, corre1, corre2, corre3, &
                                        kineListRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=16), intent(in) :: linkType
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in) :: geomDime
        aster_logical, intent(in) :: useDisp
        character(len=8), pointer :: dofName(:)
        aster_logical, intent(in) :: lApplyRota
        real(kind=8), intent(in) :: rotaMatr(3, 3)
        aster_logical, intent(in) :: useNormal
        real(kind=8), pointer :: normSlav(:)
        integer(kind=8), pointer :: nodeSkinToBody(:)
        character(len=16), intent(in) :: corres, corre1, corre2, corre3
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
! - Local
        integer(kind=8), parameter :: nbPairMaxi = 3, nbDofMaxi = 4
        character(len=8) :: mesh
        integer(kind=8) :: nbDofSlav, nbDofMast, nbEqua, iDof, nbDof
        character(len=8), dimension(nbPairMaxi, nbDofMaxi) :: dofLinkName

!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ

! - Create list of equations
        if (linkType .eq. 'MASSIF') then
            if (useNormal) then
                nbDofSlav = 1
                nbDofMast = 1
                nbEqua = 1
                dofLinkName(1, 1) = 'DEPL'
                dofLinkName(1, 2) = 'DEPL'
            elseif (useDisp) then
                if (lApplyRota) then
                    nbDofMast = 1
                    nbDofSlav = geomDime
                    nbEqua = geomDime
                    dofLinkName(1, 1) = 'DX' ! Mast
                    dofLinkName(1, 2) = 'DX' ! Slav
                    dofLinkName(1, 3) = 'DY' ! Slav
                    dofLinkName(2, 1) = 'DY' ! Mast
                    dofLinkName(2, 2) = 'DX' ! Slav
                    dofLinkName(2, 3) = 'DY' ! Slav
                    if (geomDime .eq. 3) then
                        dofLinkName(1, 4) = 'DZ' ! Slav
                        dofLinkName(2, 4) = 'DZ' ! Slav
                        dofLinkName(3, 1) = 'DZ' ! Mast
                        dofLinkName(3, 2) = 'DX' ! Slav
                        dofLinkName(3, 3) = 'DY' ! Slav
                        dofLinkName(3, 4) = 'DZ' ! Slav
                    end if
                else
                    nbDofMast = 1
                    nbDofSlav = 1
                    nbEqua = geomDime
                    dofLinkName(1, 1) = 'DX'
                    dofLinkName(1, 2) = 'DX'
                    dofLinkName(2, 1) = 'DY'
                    dofLinkName(2, 2) = 'DY'
                    if (geomDime .eq. 3) then
                        dofLinkName(3, 1) = 'DZ'
                        dofLinkName(3, 2) = 'DZ'
                    end if
                end if
            else
                nbDof = size(dofName)
                nbDofMast = 1
                nbDofSlav = 1
                nbEqua = nbDof
                do iDof = 1, nbDof
                    dofLinkName(iDof, 1) = dofName(iDof)
                    dofLinkName(iDof, 2) = dofName(iDof)
                end do
            end if
        end if

! - Create list of linear relations
        if (linkType .eq. 'MASSIF') then
            call kineLoadMeshLinkVoVo(mesh, geomDime, &
                                      lApplyRota, rotaMatr, &
                                      useNormal, normSlav, nodeSkinToBody, &
                                      corres, &
                                      nbEqua, nbDofSlav, nbDofMast, dofLinkName, &
                                      kineListRela)

        elseif (linkType .eq. 'MASSIF_COQUE') then
            call kineLoadMeshLinkVoSh(mesh, geomDime, &
                                      corres, kineListRela)

        elseif (linkType .eq. 'COQUE_MASSIF') then
            call kineLoadMeshLinkShVo(mesh, geomDime, &
                                      corres, corre1, corre2, corre3, &
                                      kineListRela)

        elseif (linkType .eq. 'COQUE') then
            call kineLoadMeshLinkShSh(mesh, geomDime, &
                                      corres, useDisp, dofName, kineListRela)

        else
            ASSERT(ASTER_FALSE)

        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGlueMeshTherLine
!
! Create list of linear relations to glue meshes (thermic)
!
! In  mesh             : mesh
! In  corres           : name of JEVEUX object for correspondence slave/master
! IO  kineListRela     : object for list of linear relations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGlueMeshTherLine(meshZ, corres, kineListRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ
        character(len=16), intent(in) :: corres
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
! - Local
        character(len=80), parameter :: title = 'LIAISON_MAILLE'
        character(len=8) :: mesh
        integer(kind=8) :: shifNodeMast, iLink, iNodeMast
        integer(kind=8) :: nbLink
        integer(kind=8) :: jvPjefNb, jvPjefNu, jvPjefCf
        integer(kind=8) :: nodeSlavNume, nodeMastNume, nbNodeMast, nbTerm
        character(len=8) :: nodeSlavName, nodeMastName
        real(kind=8) :: nodeMastCoef
        integer(kind=8), pointer :: pjefNb(:) => null(), pjefNu(:) => null()
        real(kind=8), pointer :: pjefCf(:) => null()
        aster_logical :: lcolle
        integer(kind=8) :: ier
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        lcolle = .false.
        call jeexin(meshZ(1:8)//'.NOMNOE', ier)
        if (ier .ne. 0) then
            lcolle = .true.
        end if

! - Acces to link parameters
        call jelira(corres//'.PJEF_NB', 'LONMAX', nbLink)
        call jeveuo(corres//'.PJEF_NB', 'L', jvPjefNb)
        call jeveuo(corres//'.PJEF_NU', 'L', jvPjefNu)
        call jeveuo(corres//'.PJEF_CF', 'L', jvPjefCf)
        call jeveuo(corres//'.PJEF_NB', 'L', vi=pjefNb)
        call jeveuo(corres//'.PJEF_NU', 'L', vi=pjefNu)
        call jeveuo(corres//'.PJEF_CF', 'L', vr=pjefCf)

        shifNodeMast = 0
        do iLink = 1, nbLink
            nbNodeMast = pjefNb(iLink)

            if (nbNodeMast .eq. 0) cycle

! ----- Init list of relations
            call kineListRelaInit(kineListRela)

! ----- Current slave node
            nodeSlavNume = iLink
            nodeSlavName = int_to_char8(nodeSlavNume, lcolle, meshZ, "NOEUD")

! ----- Coefficient for slave node
            kineListRela%nodeName(1) = nodeSlavName
            kineListRela%coefMultReal(1) = -1.d0

! ----- DOF to link
            kineListRela%dofName(1:nbNodeMast+1) = 'TEMP'

! ----- Coefficients for master nodes
            do iNodeMast = 1, nbNodeMast
                nodeMastNume = pjefNu(shifNodeMast+iNodeMast)
                nodeMastCoef = pjefCf(shifNodeMast+iNodeMast)
                nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
!           SI LA RELATION EST UNE TAUTOLOGIE, ON NE L'ECRIT PAS :
                if (nodeMastNume .eq. nodeSlavNume) then
                    if (abs(kineListRela%coefMultReal(iNodeMast+1)-1.0d0) .lt. 1.0d-02) then
                        call utmess('A', 'CHARGES7_1')
                        goto 290
                    end if
                end if
                kineListRela%nodeName(iNodeMast+1) = nodeMastName
                kineListRela%coefMultReal(iNodeMast+1) = nodeMastCoef
            end do

! ----- Affect linear relation
            nbTerm = nbNodeMast+1
            call kineListRelaSave(title, nbTerm, kineListRela, epsiDebg_=ASTER_TRUE)

290         continue
            shifNodeMast = shifNodeMast+nbNodeMast

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGlueMeshSetLCS
!
! Compute locale coordinate system
!
! In  geomDime         : geometric dimension (2 or 3)
! In  nbNodeMast       : number of master nodes
! In  iLink            : current link
! Ptr normSlav         : pointer to list of slave normals
! Ptr nodeSkinToBody   : pointer between body and skin nodes
! In  lApplyRota       : flag when rotation is applied
! In  rotaMatr         : rotation matrix
! IO  kineListRela     : object for list of linear relations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGlueMeshSetLCS(geomDime, nbNodeMast, &
                                      iLink, &
                                      normSlav, nodeSkinToBody, &
                                      lApplyRota, rotaMatr, &
                                      kineListRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: geomDime, nbNodeMast
        integer(kind=8), intent(in) :: iLink
        real(kind=8), pointer :: normSlav(:)
        integer(kind=8), pointer :: nodeSkinToBody(:)
        aster_logical, intent(in) :: lApplyRota
        real(kind=8), intent(in) :: rotaMatr(3, 3)
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
! - Local
        integer(kind=8) :: iNodeMast, iGeomDime, jGeomDime
        real(kind=8) :: normal(3)
!   ------------------------------------------------------------------------------------------------
!
        kineListRela%LCSType(1:nbNodeMast+1) = geomDime
        if (lApplyRota) then
            normal = 0.d0
            do iGeomDime = 1, geomDime
                do jGeomDime = 1, geomDime
                    normal(iGeomDime) = normal(iGeomDime)+ &
                                        rotaMatr(jGeomDime, iGeomDime)* &
                                        normSlav(geomDime*(nodeSkinToBody(iLink)-1)+jGeomDime)
                end do
            end do
            do iGeomDime = 1, geomDime
                kineListRela%LCSVale(iGeomDime) = &
                    normSlav(geomDime*(nodeSkinToBody(iLink)-1)+iGeomDime)
            end do
            do iNodeMast = 2, nbNodeMast+1
                do iGeomDime = 1, geomDime
                    kineListRela%LCSVale(3*(iNodeMast-1)+iGeomDime) = -normal(iGeomDime)
                end do
            end do
        else
            do iNodeMast = 1, nbNodeMast+1
                do iGeomDime = 1, geomDime
                    kineListRela%LCSVale(3*(iNodeMast-1)+iGeomDime) = &
                        normSlav(geomDime*(nodeSkinToBody(iLink)-1)+iGeomDime)
                end do
            end do
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadMeshLinkVoVo
!
! Create linear relations for 'MASSIF'
!
! Each equation is a relation between one slave node and several master nodes:
!    CoefSlav_1 * NodeSlav_Dof1 + CoefSlav_2 * NodeSlav_Dof2  + ... =
!      (i=1, nbNodeMast) {CoefMast_1 * NodeMast^i_Dof1 + CoefMast_2 * NodeMast^i_Dof2  + ... }
!
!      nbEqua: number of equations
!      nbDofSlav: number of dof for slave node for each equation
!      nbDofMast: number of dof for master nodes for each equation
!      dofName: list of dof for slave and master nodes
!
! In  mesh             : mesh
! In  geomDime         : geometric dimension (2 or 3)
! In  lApplyRota       : flag when rotation is applied
! In  rotaMatr         : rotation matrix
! In  useNormal        : flag to link with normal to slave side
! Ptr normSlav         : pointer to list of slave normals
! Ptr nodeSkinToBody   : pointer between body and skin nodes
! In  corres           : name of JEVEUX object for correspondence slave/master
! In  nbEqua           : number of equations
! In  nbDofSlav        : number of dof for slave node for each equation
! In  nbDofMast        : number of dof for master nodes for each equation
! Ptr dofLinkName      : pointer to dof for slave and master nodes for each equation
! IO  kineListRela     : object for list of linear relations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadMeshLinkVoVo(meshZ, geomDime, &
                                    lApplyRota, rotaMatr, &
                                    useNormal, normSlav, nodeSkinToBody, &
                                    corres, &
                                    nbEqua, nbDofSlav, nbDofMast, dofLinkName, &
                                    kineListRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in) :: geomDime
        aster_logical, intent(in) :: lApplyRota
        real(kind=8), intent(in) :: rotaMatr(3, 3)
        aster_logical, intent(in) :: useNormal
        real(kind=8), pointer :: normSlav(:)
        integer(kind=8), pointer :: nodeSkinToBody(:)
        character(len=16), intent(in) :: corres
        integer(kind=8), intent(in) :: nbEqua, nbDofSlav, nbDofMast
        character(len=8), dimension(:, :) :: dofLinkName
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
! - Local
        character(len=80), parameter :: title = 'LIAISON_MAIL'
        integer(kind=8) :: iNodeMast, iLink, iTerm
        integer(kind=8) :: shifNodeMast, nbLink, nodeSlavNume, nodeMastNume
        character(len=8) :: nodeSlavName, nodeMastName
        character(len=8) :: mesh, dofSlavName, dofMastName
        real(kind=8) :: nodeMastCoef
        integer(kind=8) :: nbNodeMast, nbTerm
        integer(kind=8) :: iDofMast, iDofSlav, iEqua
        integer(kind=8), pointer :: pjefNu(:) => null(), pjefNb(:) => null()
        real(kind=8), pointer :: pjefCf(:) => null()
        aster_logical :: lcolle
        integer(kind=8) :: ier
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        lcolle = .false.
        call jeexin(meshZ(1:8)//'.NOMNOE', ier)
        if (ier .ne. 0) then
            lcolle = .true.
        end if

! - Acces to link parameters
        call jelira(corres//'.PJEF_NB', 'LONMAX', nbLink)
        call jeveuo(corres//'.PJEF_NB', 'L', vi=pjefNb)
        call jeveuo(corres//'.PJEF_NU', 'L', vi=pjefNu)
        call jeveuo(corres//'.PJEF_CF', 'L', vr=pjefCf)

! - Create linear relations
        shifNodeMast = 0
        do iLink = 1, nbLink
            nbNodeMast = pjefNb(iLink)
            if (nbNodeMast .eq. 0) cycle

! ----- Current slave node
            nodeSlavNume = iLink
            nodeSlavName = int_to_char8(nodeSlavNume, lcolle, meshZ, "NOEUD")

! ----- Set local coordinate system
            if (useNormal) then
                call kineLoadGlueMeshSetLCS(geomDime, nbNodeMast, &
                                            iLink, &
                                            normSlav, nodeSkinToBody, &
                                            lApplyRota, rotaMatr, &
                                            kineListRela)
            end if

! ----- Loop on equations
            do iEqua = 1, nbEqua
                do iDofMast = 1, nbDofMast
! ------------- Number of terms
                    nbTerm = (nbNodeMast+nbDofSlav)
                    ASSERT(nbTerm .le. kineListRela%nbTermMaxi)

! ------------- Init list of relations
                    call kineListRelaInit(kineListRela)

! ------------- Coefficient for slave node
                    iTerm = 0
                    do iDofSlav = 1, nbDofSlav
                        iTerm = iTerm+1
                        dofSlavName = dofLinkName(iEqua, 1+iDofSlav)
                        kineListRela%nodeName(iTerm) = nodeSlavName
                        if (lApplyRota) then
                            kineListRela%coefMultReal(iTerm) = -rotaMatr(iDofSlav, iEqua)
                        else
                            kineListRela%coefMultReal(iTerm) = -1.d0
                        end if
                        kineListRela%dofName(iTerm) = dofSlavName
                    end do

! ------------- Coefficients for master nodes
                    dofMastName = dofLinkName(iEqua, 1)
                    do iNodeMast = 1, nbNodeMast
                        nodeMastNume = pjefNu(shifNodeMast+iNodeMast)
                        nodeMastCoef = pjefCf(shifNodeMast+iNodeMast)
                        nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
                        if (nodeMastNume .eq. nodeSlavNume) then
                            if (abs(nodeMastCoef-1.0d0) .lt. 1.0d-02) then
                                call utmess('A', 'CHARGES7_1')
                                goto 294
                            end if
                        end if
                        iTerm = iTerm+1
                        kineListRela%nodeName(iTerm) = nodeMastName
                        kineListRela%coefMultReal(iTerm) = nodeMastCoef
                        kineListRela%dofName(iTerm) = dofMastName
                    end do
                    ASSERT(iTerm .eq. nbTerm)

! ------------- Affect linear relation
                    call kineListRelaSave(title, nbTerm, kineListRela, epsiDebg_=ASTER_TRUE)

294                 continue
                end do
            end do

! ----- Next master nodes
            shifNodeMast = shifNodeMast+nbNodeMast
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGlueMeshVoSH
!
! Create linear relations for 'MASSIF_COQUE'
!
! In  mesh             : mesh
! In  geomDime         : geometric dimension (2 or 3)
! In  corres           : name of JEVEUX object for correspondence slave/master
! IO  kineListRela     : object for list of linear relations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadMeshLinkVoSh(meshZ, geomDime, &
                                    corres, kineListRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in) :: geomDime
        character(len=16), intent(in) :: corres
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
! - Local
        character(len=80), parameter :: title = 'LIAISON_MAIL-MASSIF_COQUE'
        character(len=8), parameter :: dofDispName(3) = (/'DX', 'DY', 'DZ'/)
        character(len=8), parameter :: dofRotaName(2, 3) = &
                                       reshape((/'DRY', 'DRZ', &
                                                 'DRZ', 'DRX', &
                                                 'DRX', 'DRY'/), (/2, 3/))
        integer(kind=8), parameter :: dofRotaIndx(2, 3) = &
                                      reshape((/3, 2, &
                                                1, 3, &
                                                2, 1/), (/2, 3/))
        aster_logical :: lCoque3d, lcolle
        integer(kind=8) :: iNodeMast, iNodeMastDisp, iNodeMastRota
        integer(kind=8) :: iLink, iDime, iTerm, shiftNodeMast, ier
        integer(kind=8) :: nbLink, nodeSlavNume, nodeMastNume, cellMastNume
        character(len=8) :: nodeSlavName, nodeMastName, dofName
        character(len=8) :: mesh, elrefa
        real(kind=8) :: nodeMastCoef
        integer(kind=8) :: nbNodeMast, nbNodeMastDisp, nbNodeMastRota, nbTerm
        real(kind=8), pointer :: meshVale(:) => null()
        real(kind=8) :: a(3), n2(3), an2(3), xr3(2), ff(8)
        integer(kind=8), pointer :: pjefNu(:) => null(), pjefNb(:) => null(), pjefM1(:) => null()
        real(kind=8), pointer :: pjefCf(:) => null(), pjefCo(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ

        lcolle = .false.
        call jeexin(meshZ(1:8)//'.NOMNOE', ier)
        if (ier .ne. 0) then
            lcolle = .true.
        end if

! - Acces to link parameters
        call jelira(corres//'.PJEF_NB', 'LONMAX', nbLink)
        call jeveuo(corres//'.PJEF_M1', 'L', vi=pjefM1)
        call jeveuo(corres//'.PJEF_NB', 'L', vi=pjefNb)
        call jeveuo(corres//'.PJEF_NU', 'L', vi=pjefNu)
        call jeveuo(corres//'.PJEF_CF', 'L', vr=pjefCf)
        call jeveuo(corres//'.PJEF_CO', 'L', vr=pjefCo)

! - Access to mesh
        call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=meshVale)
        ASSERT(geomDime .eq. 3)

! - Create linear relations
        shiftNodeMast = 0
        do iLink = 1, nbLink
! ----- Master cell
            cellMastNume = pjefM1(iLink)
            if (cellMastNume .eq. 0) cycle

! ----- Current slave node
            nodeSlavNume = iLink
            nodeSlavName = int_to_char8(nodeSlavNume, lcolle, meshZ, "NOEUD")

! ----- Number of master nodes
            nbNodeMast = pjefNb(iLink)
            ASSERT(nbNodeMast .ge. 3 .and. nbNodeMast .le. 9)
            lCoque3d = (nbNodeMast .eq. 7 .or. nbNodeMast .eq. 9)

! ----- Number of master nodes with displacement and rotation dof (careful with COQUE_3D !)
            nbNodeMastDisp = nbNodeMast
            if (lCoque3d) then
                nbNodeMastDisp = nbNodeMast-1
            end if
            nbNodeMastRota = nbNodeMast

! ----- Get shape functions
            if (lCoque3d) then
                xr3(1) = pjefCo(3*(nodeSlavNume-1)+1)
                xr3(2) = pjefCo(3*(nodeSlavNume-1)+2)
                if (nbNodeMastDisp .eq. 6) then
                    elrefa = 'TR6'
                else
                    elrefa = 'QU8'
                end if
                call elrfvf(elrefa, xr3, ff, nbNodeMastDisp)
            else
                do iNodeMast = 1, nbNodeMast
                    ff(iNodeMast) = pjefCf(shiftNodeMast+iNodeMast)
                end do
            end if

! ----- Compute vector between shell node and volumic
            do iDime = 1, 3
                a(iDime) = 0.d0
                do iNodeMast = 1, nbNodeMast
                    nodeMastNume = pjefNu(shiftNodeMast+iNodeMast)
                    nodeMastCoef = pjefCf(shiftNodeMast+iNodeMast)
                    a(iDime) = a(iDime)+nodeMastCoef*meshVale(3*(nodeMastNume-1)+iDime)
                end do
                n2(iDime) = meshVale(3*(nodeSlavNume-1)+iDime)
                an2(iDime) = n2(iDime)-a(iDime)
            end do

! ----- Total number of terms for each link
            nbTerm = 1+2*nbNodeMastRota+nbNodeMastDisp
            ASSERT(nbTerm .le. kineListRela%nbTermMaxi)

! ----- Set local coordinate system
            do iDime = 1, 3
                dofName = dofDispName(iDime)

! --------- Init list of relations
                call kineListRelaInit(kineListRela)

! --------- Slave node (displacements only)
                iTerm = 1
                kineListRela%dofName(iTerm) = dofName
                kineListRela%nodeName(iTerm) = nodeSlavName
                kineListRela%coefMultReal(iTerm) = -1.d0
                iTerm = iTerm+1

! --------- Loop on master nodes (displacements)
                do iNodeMastDisp = 1, nbNodeMastDisp
                    nodeMastNume = pjefNu(shiftNodeMast+iNodeMastDisp)
                    nodeMastCoef = ff(iNodeMastDisp)
                    nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
                    if (nodeMastNume .eq. nodeSlavNume) then
                        if (abs(nodeMastCoef-1.0d0) .lt. 1.0d-02) then
                            call utmess('A', 'CHARGES7_1')
                            goto 295
                        end if
                    end if
                    kineListRela%nodeName(iTerm) = nodeMastName
                    kineListRela%dofName(iTerm) = dofName
                    kineListRela%coefMultReal(iTerm) = nodeMastCoef
                    iTerm = iTerm+1
                end do

! --------- Loop on master nodes (rotations)
                do iNodeMastRota = 1, nbNodeMastRota
                    nodeMastNume = pjefNu(shiftNodeMast+iNodeMastRota)
                    nodeMastCoef = pjefCf(shiftNodeMast+iNodeMastRota)
                    nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
                    kineListRela%nodeName(iTerm-1+1) = nodeMastName
                    kineListRela%nodeName(iTerm-1+2) = nodeMastName
                    kineListRela%dofName(iTerm-1+1) = dofRotaName(1, iDime)
                    kineListRela%dofName(iTerm-1+2) = dofRotaName(2, iDime)
                    kineListRela%coefMultReal(iTerm-1+1) = &
                        +nodeMastCoef*an2(dofRotaIndx(1, iDime))
                    kineListRela%coefMultReal(iTerm-1+2) = &
                        -nodeMastCoef*an2(dofRotaIndx(2, iDime))
                    iTerm = iTerm+2
                end do

                ASSERT(iTerm .eq. nbTerm+1)

! --------- Affect linear relation
                call kineListRelaSave(title, nbTerm, kineListRela)

295             continue

            end do
! ----- Next master node
            shiftNodeMast = shiftNodeMast+nbNodeMast

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadMeshLinkShVo
!
! Create linear relations for 'COQUE_MASSIF'
!
! In  mesh             : mesh
! In  geomDime         : geometric dimension (2 or 3)
! In  corres           : name of JEVEUX object for correspondence slave/master
! In  corres           : name of JEVEUX object for correspondence slave/master
! In  corre1           : name of JEVEUX object for correspondence slave/master (superior face)
! In  corre2           : name of JEVEUX object for correspondence slave/master (inferior face)
! In  corre3           : name of JEVEUX object for normal thickness
! IO  kineListRela     : object for list of linear relations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadMeshLinkShVo(meshZ, geomDime, &
                                    corres, corre1, corre2, corre3, &
                                    kineListRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in) :: geomDime
        character(len=16), intent(in) :: corres, corre1, corre2, corre3
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
! - Local
        character(len=80), parameter :: title = 'LIAISON_MAIL-COQUE_MASSIF'
        character(len=8), parameter :: dofDispName(3) = (/'DX', 'DY', 'DZ'/)
        integer(kind=8) :: iNodeMast, iNodeMastSup, iNodeMastInf, iTerm, iDime
        integer(kind=8) :: shiftNodeMast, shiftNodeMastSup, shiftNodeMastInf, cellMastNume
        integer(kind=8) :: iLink, nbLink, nodeSlavNume, nodeMastNume, ier
        character(len=8) :: nodeSlavName, nodeMastName, dofName
        character(len=8) :: mesh
        real(kind=8) :: nodeMastCoef, thickness, normalVect(3), mat33(3, 3)
        integer(kind=8) :: nbNodeMast, nbNodeMastSup, nbNodeMastInf, nbTerm
        real(kind=8), pointer :: meshVale(:) => null()
        integer(kind=8), pointer :: pjefNu(:) => null(), pjefNb(:) => null(), pjefM1(:) => null()
        real(kind=8), pointer :: pjefCf(:) => null(), pjefCo(:) => null()
        integer(kind=8), pointer :: pjef1Nu(:) => null(), pjef1Nb(:) => null()
        real(kind=8), pointer :: pjef1Cf(:) => null()
        integer(kind=8), pointer :: pjef2Nu(:) => null(), pjef2Nb(:) => null()
        real(kind=8), pointer :: pjef2Cf(:) => null()
        real(kind=8), pointer :: thickNormal(:) => null()
        aster_logical :: lcolle
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ

        lcolle = .false.
        call jeexin(meshZ(1:8)//'.NOMNOE', ier)
        if (ier .ne. 0) then
            lcolle = .true.
        end if

! - Acces to link parameters
        call jelira(corres//'.PJEF_NB', 'LONMAX', nbLink)
        call jeveuo(corres//'.PJEF_M1', 'L', vi=pjefM1)
        call jeveuo(corres//'.PJEF_NB', 'L', vi=pjefNb)
        call jeveuo(corres//'.PJEF_NU', 'L', vi=pjefNu)
        call jeveuo(corres//'.PJEF_CF', 'L', vr=pjefCf)
        call jeveuo(corres//'.PJEF_CO', 'L', vr=pjefCo)
        call jeveuo(corre1//'.PJEF_NB', 'L', vi=pjef1Nb)
        call jeveuo(corre1//'.PJEF_NU', 'L', vi=pjef1Nu)
        call jeveuo(corre1//'.PJEF_CF', 'L', vr=pjef1Cf)
        call jeveuo(corre2//'.PJEF_NB', 'L', vi=pjef2Nb)
        call jeveuo(corre2//'.PJEF_NU', 'L', vi=pjef2Nu)
        call jeveuo(corre2//'.PJEF_CF', 'L', vr=pjef2Cf)
        call jeveuo(corre3, 'L', vr=thickNormal)

! - Access to mesh
        call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=meshVale)
        ASSERT(geomDime .eq. 3)

! - Create linear relations
        shiftNodeMast = 0
        shiftNodeMastSup = 0
        shiftNodeMastInf = 0
        do iLink = 1, nbLink
! ----- Master cell
            cellMastNume = pjefM1(iLink)
            if (cellMastNume .eq. 0) cycle

! ----- Current slave node
            nodeSlavNume = iLink
            nodeSlavName = int_to_char8(nodeSlavNume, lcolle, meshZ, "NOEUD")

! ----- Number of master nodes
            nbNodeMast = pjefNb(iLink)

! ----- Link for DISPLACEMENTS
            nbTerm = 1+nbNodeMast
            ASSERT(nbTerm .le. kineListRela%nbTermMaxi)

            do iDime = 1, 3
                dofName = dofDispName(iDime)

! --------- Init list of relations
                call kineListRelaInit(kineListRela)

! --------- Slave node (displacements)
                iTerm = 1
                kineListRela%nodeName(iTerm) = nodeSlavName
                kineListRela%coefMultReal(iTerm) = -1.d0
                kineListRela%dofName(iTerm) = dofName

! --------- Master nodes (displacements)
                do iNodeMast = 1, nbNodeMast
                    nodeMastNume = pjefNu(shiftNodeMast+iNodeMast)
                    nodeMastCoef = pjefCf(shiftNodeMast+iNodeMast)
                    nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
                    if (nodeMastNume .eq. nodeSlavNume) then
                        if (abs(nodeMastCoef-1.0d0) .lt. 1.0d-02) then
                            call utmess('A', 'CHARGES7_1')
                            goto 291
                        end if
                    end if
                    iTerm = iTerm+1
                    kineListRela%nodeName(iTerm) = nodeMastName
                    kineListRela%dofName(iTerm) = dofName
                    kineListRela%coefMultReal(iTerm) = nodeMastCoef
                end do

! --------- Affect linear relation
                call kineListRelaSave(title, nbTerm, kineListRela, epsiDebg_=ASTER_TRUE)

291             continue
            end do

! ----- Nodes linked to superior and inferior faces
            nbNodeMastSup = pjef1Nb(nodeSlavNume)
            nbNodeMastInf = pjef2Nb(nodeSlavNume)

! ----- Link for ROTATIONS
            nbTerm = 3*(1+nbNodeMastInf+nbNodeMastSup)
            ASSERT(nbTerm .le. kineListRela%nbTermMaxi)

! ----- Init list of relations
            call kineListRelaInit(kineListRela)

! ----- Thickness of shell
            normalVect = thickNormal(3*(iLink-1)+1:3*(iLink-1)+3)
            thickness = sqrt(normalVect(1)**2+normalVect(2)**2+normalVect(3)**2)
            call base3n(normalVect, mat33)

! ----- Slave node
            iTerm = 1
            kineListRela%nodeName(iTerm) = nodeSlavName
            kineListRela%coefMultReal(iTerm) = -1.d0*mat33(1, 2)
            kineListRela%dofName(iTerm) = 'DRX'
            iTerm = iTerm+1
            kineListRela%nodeName(iTerm) = nodeSlavName
            kineListRela%coefMultReal(iTerm) = -1.d0*mat33(2, 2)
            kineListRela%dofName(iTerm) = 'DRY'
            iTerm = iTerm+1
            kineListRela%nodeName(iTerm) = nodeSlavName
            kineListRela%coefMultReal(iTerm) = -1.d0*mat33(3, 2)
            kineListRela%dofName(iTerm) = 'DRZ'

! ----- Master nodes on superior face (translations)
            do iNodeMastSup = 1, nbNodeMastSup
                nodeMastNume = pjef1Nu(shiftNodeMastSup+iNodeMastSup)
                nodeMastCoef = pjef1Cf(shiftNodeMastSup+iNodeMastSup)
                nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DX'
                kineListRela%coefMultReal(iTerm) = -(nodeMastCoef*mat33(1, 3))/thickness
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DY'
                kineListRela%coefMultReal(iTerm) = -(nodeMastCoef*mat33(2, 3))/thickness
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DZ'
                kineListRela%coefMultReal(iTerm) = -(nodeMastCoef*mat33(3, 3))/thickness
            end do

! ----- Master nodes on inferior face (translations)
            do iNodeMastInf = 1, nbNodeMastInf
                nodeMastNume = pjef2Nu(shiftNodeMastInf+iNodeMastInf)
                nodeMastCoef = pjef2Cf(shiftNodeMastInf+iNodeMastInf)
                nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DX'
                kineListRela%coefMultReal(iTerm) = +(nodeMastCoef*mat33(1, 3))/thickness
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DY'
                kineListRela%coefMultReal(iTerm) = +(nodeMastCoef*mat33(2, 3))/thickness
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DZ'
                kineListRela%coefMultReal(iTerm) = +(nodeMastCoef*mat33(3, 3))/thickness
            end do
            ASSERT(nbTerm .eq. iTerm)

! ----- Affect linear relation
            call kineListRelaSave(title, nbTerm, kineListRela, epsiDebg_=ASTER_FALSE)

! ----- Init list of relations
            call kineListRelaInit(kineListRela)

! ----- Slave node
            iTerm = 1
            kineListRela%nodeName(iTerm) = nodeSlavName
            kineListRela%coefMultReal(iTerm) = -1.d0*mat33(1, 3)
            kineListRela%dofName(iTerm) = 'DRX'
            iTerm = iTerm+1
            kineListRela%nodeName(iTerm) = nodeSlavName
            kineListRela%coefMultReal(iTerm) = -1.d0*mat33(2, 3)
            kineListRela%dofName(iTerm) = 'DRY'
            iTerm = iTerm+1
            kineListRela%nodeName(iTerm) = nodeSlavName
            kineListRela%coefMultReal(iTerm) = -1.d0*mat33(3, 3)
            kineListRela%dofName(iTerm) = 'DRZ'

! ----- Master nodes on superior face (rotations)
            do iNodeMastSup = 1, nbNodeMastSup
                nodeMastNume = pjef1Nu(shiftNodeMastSup+iNodeMastSup)
                nodeMastCoef = pjef1Cf(shiftNodeMastSup+iNodeMastSup)
                nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DX'
                kineListRela%coefMultReal(iTerm) = +(nodeMastCoef*mat33(1, 2))/thickness
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DY'
                kineListRela%coefMultReal(iTerm) = +(nodeMastCoef*mat33(2, 2))/thickness
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DZ'
                kineListRela%coefMultReal(iTerm) = +(nodeMastCoef*mat33(3, 2))/thickness
            end do

! ----- Master nodes on inferior face (translations)
            do iNodeMastInf = 1, nbNodeMastInf
                nodeMastNume = pjef2Nu(shiftNodeMastInf+iNodeMastInf)
                nodeMastCoef = pjef2Cf(shiftNodeMastInf+iNodeMastInf)
                nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DX'
                kineListRela%coefMultReal(iTerm) = -(nodeMastCoef*mat33(1, 2))/thickness
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DY'
                kineListRela%coefMultReal(iTerm) = -(nodeMastCoef*mat33(2, 2))/thickness
                iTerm = iTerm+1
                kineListRela%nodeName(iTerm) = nodeMastName
                kineListRela%dofName(iTerm) = 'DZ'
                kineListRela%coefMultReal(iTerm) = -(nodeMastCoef*mat33(3, 2))/thickness
            end do
            ASSERT(nbTerm .eq. iTerm)

! ----- Affect linear relation
            call kineListRelaSave(title, nbTerm, kineListRela, epsiDebg_=ASTER_FALSE)

! ----- Next master nodes
            shiftNodeMast = shiftNodeMast+nbNodeMast
            shiftNodeMastSup = shiftNodeMastSup+nbNodeMastSup
            shiftNodeMastInf = shiftNodeMastInf+nbNodeMastInf
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadMeshLinkShSh
!
! Create linear relations for 'COQUE'
!
! In  mesh             : mesh
! In  geomDime         : geometric dimension (2 or 3)
! In  corres           : name of JEVEUX object for correspondence slave/master
! IO  kineListRela     : object for list of linear relations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadMeshLinkShSh(meshZ, geomDime, &
                                    corres, useDisp, dofName, kineListRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in)                 :: geomDime
        character(len=16), intent(in)       :: corres
        aster_logical, intent(in)           :: useDisp
        character(len=8), pointer           :: dofName(:)
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
! - Local
        character(len=80), parameter :: title = 'LIAISON_MAIL-COQUE'
        character(len=8), parameter :: dofDispName(3) = (/'DX', 'DY', 'DZ'/)
        character(len=8), parameter :: dofRotaName(3) = (/'DRX', 'DRY', 'DRZ'/)
        integer(kind=8) :: iNodeMast, iTerm, iDime, iDof, nbDof
        integer(kind=8) :: shiftNodeMast
        integer(kind=8) :: iLink, nbLink, nodeSlavNume, nodeMastNume, cellMastNume
        character(len=8) :: nodeSlavName, nodeMastName, LocDofName
        character(len=8) :: mesh
        real(kind=8) :: nodeMastCoef
        integer(kind=8) :: nbNodeMast, nbTerm
        real(kind=8), pointer :: meshVale(:) => null()
        integer(kind=8), pointer :: pjefNu(:) => null(), pjefNb(:) => null(), pjefM1(:) => null()
        real(kind=8), pointer :: pjefCf(:) => null(), pjefCo(:) => null()
        aster_logical :: TakeDof, lcolle
        integer(kind=8) :: ier
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        lcolle = .false.
        call jeexin(meshZ(1:8)//'.NOMNOE', ier)
        if (ier .ne. 0) then
            lcolle = .true.
        end if

! - Access to link parameters
        call jelira(corres//'.PJEF_NB', 'LONMAX', nbLink)
        call jeveuo(corres//'.PJEF_M1', 'L', vi=pjefM1)
        call jeveuo(corres//'.PJEF_NB', 'L', vi=pjefNb)
        call jeveuo(corres//'.PJEF_NU', 'L', vi=pjefNu)
        call jeveuo(corres//'.PJEF_CF', 'L', vr=pjefCf)
        call jeveuo(corres//'.PJEF_CO', 'L', vr=pjefCo)

! - Access to mesh
        call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=meshVale)
        ASSERT(geomDime .eq. 3)

! - Create linear relations
        shiftNodeMast = 0
        do iLink = 1, nbLink
! ----- Master cell
            cellMastNume = pjefM1(iLink)
            if (cellMastNume .eq. 0) cycle

! ----- Current slave node
            nodeSlavNume = iLink
            nodeSlavName = int_to_char8(nodeSlavNume, lcolle, meshZ, "NOEUD")

! ----- Number of master nodes
            nbNodeMast = pjefNb(iLink)
            ASSERT(nbNodeMast .ge. 3 .and. nbNodeMast .le. 9)

! ----- Link for DISPLACEMENTS
            nbTerm = 1+nbNodeMast
            ASSERT(nbTerm .le. kineListRela%nbTermMaxi)

! ----- Number of Dof to consider
            nbDof = size(dofName)
            ! Not necessary because the rule is in the catalog, but ...
            !   useDisp : True, if DDL is not used ==> (not useDisp) DDL is used
            if (.not. useDisp) then
                ASSERT((nbDof .ge. 1) .and. (nbDof .le. 6))
            end if

            do iDime = 1, 3
                LocDofName = dofDispName(iDime)
! --------- Is this Dof is taken into account
                if (.not. useDisp) then
                    TakeDof = ASTER_FALSE
                    do iDof = 1, nbDof
                        if (LocDofName .eq. dofName(iDof)) then
                            TakeDof = ASTER_TRUE
                            exit
                        end if
                    end do
                    if (.not. TakeDof) then
                        cycle
                    end if
                end if

! --------- Init list of relations
                call kineListRelaInit(kineListRela)

! --------- Slave node (displacements)
                iTerm = 1
                kineListRela%nodeName(iTerm) = nodeSlavName
                kineListRela%coefMultReal(iTerm) = -1.d0
                kineListRela%dofName(iTerm) = LocDofName

! --------- Master nodes (displacements)
                do iNodeMast = 1, nbNodeMast
                    nodeMastNume = pjefNu(shiftNodeMast+iNodeMast)
                    nodeMastCoef = pjefCf(shiftNodeMast+iNodeMast)
                    nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
                    if (nodeMastNume .eq. nodeSlavNume) then
                        if (abs(nodeMastCoef-1.0d0) .lt. 1.0d-02) then
                            call utmess('A', 'CHARGES7_1')
                            goto 292
                        end if
                    end if
                    iTerm = iTerm+1
                    kineListRela%nodeName(iTerm) = nodeMastName
                    kineListRela%dofName(iTerm) = LocDofName
                    kineListRela%coefMultReal(iTerm) = nodeMastCoef
                end do

! --------- Affect linear relation
                call kineListRelaSave(title, nbTerm, kineListRela, epsiDebg_=ASTER_TRUE)

292             continue

            end do

! ----- Link for ROTATIONS
            nbTerm = 1+nbNodeMast
            ASSERT(nbTerm .le. kineListRela%nbTermMaxi)

            do iDime = 1, 3
                LocDofName = dofRotaName(iDime)
! --------- Is this Dof is taken into account
                if (.not. useDisp) then
                    TakeDof = ASTER_FALSE
                    do iDof = 1, nbDof
                        if (LocDofName .eq. dofName(iDof)) then
                            TakeDof = ASTER_TRUE
                            exit
                        end if
                    end do
                    if (.not. TakeDof) then
                        cycle
                    end if
                end if

! --------- Init list of relations
                call kineListRelaInit(kineListRela)

! --------- Slave node (displacements)
                iTerm = 1
                kineListRela%nodeName(iTerm) = nodeSlavName
                kineListRela%coefMultReal(iTerm) = -1.d0
                kineListRela%dofName(iTerm) = LocDofName

! --------- Master nodes (displacements)
                do iNodeMast = 1, nbNodeMast
                    nodeMastNume = pjefNu(shiftNodeMast+iNodeMast)
                    nodeMastCoef = pjefCf(shiftNodeMast+iNodeMast)
                    nodeMastName = int_to_char8(nodeMastNume, lcolle, meshZ, "NOEUD")
                    if (nodeMastNume .eq. nodeSlavNume) then
                        if (abs(nodeMastCoef-1.0d0) .lt. 1.0d-02) then
                            call utmess('A', 'CHARGES7_1')
                            goto 293
                        end if
                    end if
                    iTerm = iTerm+1
                    kineListRela%nodeName(iTerm) = nodeMastName
                    kineListRela%dofName(iTerm) = LocDofName
                    kineListRela%coefMultReal(iTerm) = nodeMastCoef
                end do

! --------- Affect linear relation
                call kineListRelaSave(title, nbTerm, kineListRela, epsiDebg_=ASTER_TRUE)

293             continue

            end do

! ----- Next master nodes
            shiftNodeMast = shiftNodeMast+nbNodeMast
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadLinkProj
!
! Main subroutine to glue mesh from projection (LIAISON_PROJ)
!
! In  model            : model
! In  valeType         : affected value type (real, complex or function)
! Out listLineRela     : name of datastructure for list of linear relation
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadLinkProj(modelZ, valeTypeZ, listLineRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: modelZ, valeTypeZ
        character(len=19), intent(out) :: listLineRela
! - Local
!   NOMBRE MAX DE TERMES D'UNE RELATION LINEAIRE EN 3D = 1 + 3*27 (max maille: 27 noeuds)
        integer(kind=8), parameter :: nbTermMaxi = 1+3*MT_NNOMAX
        character(len=80), parameter :: title = 'LIAISON_PROJ'
        character(len=16), parameter :: factorKeyword = 'LIAISON_PROJ'
        character(len=8), parameter :: physQuanName = "DEPL_R"
        character(len=19) :: modelLigrel
        character(len=8) :: mesh
        character(len=16) :: meshLink
        character(len=4) :: valeType
        integer(kind=8) :: geomDime
        integer(kind=8) :: iOcc, nocc, nbDof
        character(len=8), pointer :: dofName(:) => null()
        integer(kind=8) :: iNode1, iNode2, iDof
        integer(kind=8) :: nbNode1, nbNode2
        integer(kind=8) :: cell1Nume, cell1TypeNume, node1Nume, node2Nume
        integer(kind=8) :: idecal, iTerm, nbTerm
        character(len=8) :: node1Name, node2Name, dofLocaName
        aster_logical :: hasExcent, oneDofDoesntExist
        aster_logical, pointer :: dofExist(:) => null()
        integer(kind=8), pointer :: cell1List(:) => null(), cellType(:) => null()
        integer(kind=8), pointer :: node1List(:) => null(), node2List(:) => null()
        real(kind=8), pointer :: nodeCoef(:) => null(), nodeCoor(:) => null()
        real(kind=8) :: coeffi, xyzom(3), coefZero
        type(KINE_LIST_RELA) :: kineListRela
        integer(kind=8) :: physQuanNbCmp, jvPhysQuanCmpName, nbec, jvPrnm
        character(len=1) :: dofUnknown
        integer(kind=8), parameter :: nbDofRota = 3
        character(len=8), pointer :: dofRotaName(:) => null()
        aster_logical :: lcolle
        integer(kind=8) :: ier

!   ------------------------------------------------------------------------------------------------
!
        valeType = valeTypeZ
        ASSERT(valeType .eq. 'REEL')

! - List of rotation DOF
        AS_ALLOCATE(vk8=dofRotaName, size=nbDofRota)
        dofRotaName(1) = "DRX"
        dofRotaName(2) = "DRY"
        dofRotaName(3) = "DRZ"

! - About model
        call dismoi('NOM_MAILLA', modelZ, 'MODELE', repk=mesh)
        call dismoi('DIM_GEOM', modelZ, 'MODELE', repi=geomDime)
        call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)
        call jeveuo(modelLigrel//'.PRNM', 'L', jvPrnm)

        lcolle = .false.
        call jeexin(mesh(1:8)//'.NOMNOE', ier)
        if (ier .ne. 0) then
            lcolle = .true.
        end if

! - About mesh
        call jeveuo(mesh//'.TYPMAIL', 'L', vi=cellType)
        call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=nodeCoor)
        if (.not. (geomDime .eq. 2 .or. geomDime .eq. 3)) then
            call utmess('F', 'CHARGES7_6')
        end if

! - Name of datastructure for list of linear relations
        listLineRela = '&&CALIRC.RLLISTE'

! - Get information about physical quantity
        call jeveuo(jexnom('&CATA.GD.NOMCMP', physQuanName), 'L', jvPhysQuanCmpName)
        call jelira(jexnom('&CATA.GD.NOMCMP', physQuanName), 'LONMAX', physQuanNbCmp)
        call dismoi('NB_EC', physQuanName, 'GRANDEUR', repi=nbec)

! - Create object for list of linear relations
        call kineListRelaCreate('Implicit', nbTermMaxi, listLineRela, kineListRela)
        coefZero = kineListRela%coefMultTole

        call getfac(factorKeyword, nocc)
        do iOcc = 1, nocc

! ----- Get main parameters from user
            call kineLoadLinkProjPara(factorKeyword, iOcc, mesh, &
                                      meshLink, hasExcent, &
                                      nbDof, dofName, dofUnknown, dofExist)

! ----- Acces to mesh link
            call jelira(meshLink//'.PJEF_NB', 'LONMAX', nbNode2)
            call jeveuo(meshLink//'.PJEF_NB', 'L', vi=node2List)
            call jeveuo(meshLink//'.PJEF_M1', 'L', vi=cell1List)
            call jeveuo(meshLink//'.PJEF_NU', 'L', vi=node1List)
            call jeveuo(meshLink//'.PJEF_CF', 'L', vr=nodeCoef)

! ----- Loop on slave nodes
            idecal = 0
            do iNode2 = 1, nbNode2

! --------- Current slave node
                node2Nume = iNode2
                node2Name = int_to_char8(node2Nume, lcolle, mesh, "NOEUD")
                cell1Nume = cell1List(node2Nume)
                if (cell1Nume .eq. 0) cycle
                if (hasExcent) then
                    cell1TypeNume = cellType(cell1Nume)
                    if (cell1TypeNume .ne. MT_TRIA3 .and. cell1TypeNume .ne. MT_QUAD4) then
                        call utmess('F', 'CHARGES7_12')
                    end if
                    if (geomDime .ne. 3) then
                        call utmess('F', 'CHARGES7_15')
                    end if
                end if

! --------- Check components on current slave node
                call kineLoadCheckCmpOnNode(jvPrnm, node2Nume, &
                                            physQuanNbCmp, jvPhysQuanCmpName, &
                                            nbDof, nbec, dofName, &
                                            dofExist, oneDofDoesntExist)

! --------- Relation for current slave node
                do iDof = 1, nbDof
! ------------- Init list of relations
                    call kineListRelaInit(kineListRela)
                    iTerm = 0
                    dofLocaName = dofName(iDof)

                    if (dofExist(iDof)) then
                        iTerm = iTerm+1
                        kineListRela%nodeName(iTerm) = node2Name
                        kineListRela%coefMultReal(iTerm) = -1.d0
                        kineListRela%dofName(iTerm) = dofLocaName
                    else
                        if (dofUnknown .ne. " ") then
                            call utmess(dofUnknown, "CHARGES7_14")
                        end if
                        cycle
                    end if

! ------------- None of the DOF exists on slave node
                    if (iTerm .eq. 0) then
                        goto 130
                    end if

! ------------- Loop on master nodes
                    nbNode1 = node2List(node2Nume)
                    do iNode1 = 1, nbNode1

! ----------------- Current master node
                        node1Nume = node1List(idecal+iNode1)
                        node1Name = int_to_char8(node1Nume, lcolle, mesh, "NOEUD")
                        coeffi = nodeCoef(idecal+iNode1)
                        if (node1Nume .eq. node2Nume) then
                            call utmess("A", "CHARGES7_13")
                            goto 130
                        end if

! ----------------- Check components on current master node
                        if (hasExcent) then

                            if (dofExist(iDof)) then
                                iTerm = iTerm+1
                                kineListRela%nodeName(iTerm) = node1Name
                                kineListRela%coefMultReal(iTerm) = coeffi
                                kineListRela%dofName(iTerm) = dofLocaName
                            else
                                if (dofUnknown .ne. " ") then
                                    call utmess(dofUnknown, "CHARGES7_14")
                                end if
                                cycle
                            end if

! --------------------- Est-ce que les ddl de rotation existent sur le noeud d'en face ?
                            call kineLoadCheckCmpOnNode(jvPrnm, node1Nume, &
                                                        physQuanNbCmp, jvPhysQuanCmpName, &
                                                        nbDofRota, nbec, dofRotaName, &
                                                        dofExist, oneDofDoesntExist)
                            if (oneDofDoesntExist) then
                                call utmess('F', 'CHARGES7_16')
                            end if

! --------------------- U1 = UI2 + DRI2^O2M ==> -UI2 + U1*al(i) - DRI2^O2M*al(i)
                            xyzom(1) = nodeCoor(3*(node2Nume-1)+1)-nodeCoor(3*(node1Nume-1)+1)
                            xyzom(2) = nodeCoor(3*(node2Nume-1)+2)-nodeCoor(3*(node1Nume-1)+2)
                            xyzom(3) = nodeCoor(3*(node2Nume-1)+3)-nodeCoor(3*(node1Nume-1)+3)

! --------------------- Add relations
                            if (dofLocaName .eq. 'DX') then
                                call kineLoadApplyEccentricity(3, node1Name, "DRY", &
                                                               coeffi, coefZero, xyzom, &
                                                               iTerm, kineListRela)
                                call kineLoadApplyEccentricity(2, node1Name, "DRZ", &
                                                               -coeffi, coefZero, xyzom, &
                                                               iTerm, kineListRela)
                            elseif (dofLocaName .eq. 'DY') then
                                call kineLoadApplyEccentricity(3, node1Name, "DRX", &
                                                               -coeffi, coefZero, xyzom, &
                                                               iTerm, kineListRela)
                                call kineLoadApplyEccentricity(1, node1Name, "DRZ", &
                                                               coeffi, coefZero, xyzom, &
                                                               iTerm, kineListRela)
                            elseif (dofLocaName .eq. 'DZ') then
                                call kineLoadApplyEccentricity(2, node1Name, "DRX", &
                                                               coeffi, coefZero, xyzom, &
                                                               iTerm, kineListRela)
                                call kineLoadApplyEccentricity(1, node1Name, "DRY", &
                                                               -coeffi, coefZero, xyzom, &
                                                               iTerm, kineListRela)
                            end if
                        else
                            if (abs(coeffi) .gt. coefZero) then
                                call kineLoadCheckCmpOnNode(jvPrnm, node1Nume, &
                                                            physQuanNbCmp, jvPhysQuanCmpName, &
                                                            nbDof, nbec, dofName, &
                                                            dofExist, oneDofDoesntExist)
                                if (dofExist(iDof)) then
                                    iTerm = iTerm+1
                                    kineListRela%nodeName(iTerm) = node1Name
                                    kineListRela%coefMultReal(iTerm) = coeffi
                                    kineListRela%dofName(iTerm) = dofName(iDof)
                                else
                                    if (dofUnknown .ne. " ") then
                                        call utmess(dofUnknown, "CHARGES7_14")
                                    end if
                                    cycle
                                end if
                            end if
                        end if
                    end do

! ------------- Affect linear relation
                    nbTerm = iTerm
                    call kineListRelaSave(title, nbTerm, kineListRela, epsiDebg_=ASTER_TRUE)

                end do
130             continue
                idecal = idecal+nbNode1
            end do

! ----- Clean
            AS_DEALLOCATE(vk8=dofName)
            AS_DEALLOCATE(vl=dofExist)
        end do

! - Clean
        call kineListRelaDelete(kineListRela)
        AS_DEALLOCATE(vk8=dofRotaName)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadLinkProjPara
!
! Get parameters for LIAISON_PROJ
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadLinkProjPara(factorKeywordZ, iOcc, mesh, &
                                    meshLink, hasExcent, &
                                    nbDof, dofName, dofUnknown, dofExist)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
        character(len=8), intent(in) :: mesh
        character(len=16), intent(out) :: meshLink
        integer(kind=8), intent(out) :: nbDof
        character(len=8), pointer :: dofName(:)
        aster_logical, pointer :: dofExist(:)
        aster_logical, intent(out) :: hasExcent
        character(len=1), intent(out) :: dofUnknown
! - Local
        integer(kind=8) :: iret
        character(len=16) :: answer
        character(len=8) :: mesh1, mesh2
        character(len=24), pointer :: meshLinkPjxx(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        call getvtx(factorKeywordZ, 'DDL', iocc=iOcc, nbval=0, nbret=nbDof)
        nbDof = -nbDof
        ASSERT(nbDof .ge. 1)
        AS_ALLOCATE(size=nbDof, vk8=dofName)
        AS_ALLOCATE(size=nbDof, vl=dofExist)
        call getvtx(factorKeywordZ, 'DDL', iocc=iocc, nbval=nbDof, vect=dofName)
        call getvid(factorKeywordZ, 'MATR_PROJECTION', iocc=iocc, scal=meshLink)
        call getvtx(factorKeywordZ, 'TYPE', iocc=iocc, scal=answer, nbret=iret)
        hasExcent = ASTER_FALSE
        if (iret .ne. 0) then
            hasExcent = answer .eq. 'EXCENTREMENT'
        end if
        !call getvtx(factorKeywordZ, 'DDL_EXIST', iocc=iocc, scal=answer, nbret=iret)
        dofUnknown = "F"
        ! if (iret .ne. 0) then
        !     if (answer .eq. "ERREUR") dofUnknown = "F"
        !     if (answer .eq. "ALARME") dofUnknown = "A"
        !     if (answer .eq. "IGNORE") dofUnknown = " "
        ! endif

        call jeveuo(meshLink//'.PJXX_K1', 'L', vk24=meshLinkPjxx)
        mesh1 = meshLinkPjxx(1) (1:8)
        mesh2 = meshLinkPjxx(2) (1:8)
        if ((mesh .ne. mesh1) .or. (mesh .ne. mesh2)) then
            call utmess('F', 'CHARGES7_11')
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module LoadKinematic_module
