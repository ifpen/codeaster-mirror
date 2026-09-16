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
subroutine pairWrap(method, &
                    mesh, newgeo, mastConxInvName, &
                    mastNeighName, slavNeighName, &
                    pairTole, areaTole, distRatio, verbosity, &
                    nbCellMast, listCellMast, &
                    nbCellSlav, listCellSlav, &
                    nbNodeMast, listNodeMast, &
                    nbPairZone, baseName)
!
    use MeshPairing_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/aplcpgn.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedetr.h"
#include "asterfort/mesh_pairing_type.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: method
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: newgeo, mastConxInvName
    character(len=24), intent(in) :: mastNeighName, slavNeighName
    real(kind=8), intent(in) :: pairTole, areaTole, distRatio
    integer(kind=8), intent(in) :: verbosity
    integer(kind=8), intent(in) :: nbCellMast, listCellMast(nbCellMast)
    integer(kind=8), intent(in) :: nbCellSlav, listCellSlav(nbCellSlav)
    integer(kind=8), intent(in) :: nbNodeMast, listNodeMast(nbNodeMast)
    integer(kind=8), intent(out) :: nbPairZone
    character(len=8), intent(in) :: baseName
!
! --------------------------------------------------------------------------------------------------
!
! Pairing segment to segment
!
! Building a list of paired cells by PANG method
!
! --------------------------------------------------------------------------------------------------
!
! In  method           : method of pairing
! In  mesh             : mesh
! In  newgeo           : updated coordinates of nodes
! In  mastConxInvName  : name of object for inverse connectivity of master cells on current zone
! In  mastNeighName    : name of object for neighbours of master cells
! In  slavNeighName    : name of object for neighbours of slave cells
! In  pairTole         : tolerance for projection (all operations ! )
! In  distRatio        : tolerance from DIST_RATIO
! In  verbosity        : level of verbosity
! In  nbCellSlav       : number of slave cells
! In  listCellSlav     : list of slave cells
! In  nbCellMast       : number of master cells
! In  listCellMast     : list of master cells
! In  nbNodeMast       : number of master nodes
! In  listNodeMast     : list of master nodes
! Out nbPairZone       : number of paired cells
! In  baseName         : JEVEUX base name for output objects
! In  zonePair         : name of datastructure for list of paired cells
! In  zoneNbPoinInte   : name of datastructure for number of integration points on slave cell
! In  zonePoinInte     : name of datastructure for coordinates of integration points on slave cell
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: pairDime, jvData, spaceDime
    type(MESH_PAIRING) :: meshPairing
    character(len=24) :: zonePair, zoneNbPoinInte, zonePoinInte
!
! --------------------------------------------------------------------------------------------------
!
    nbPairZone = 0

! - Protection
    if (nbCellSlav .eq. 0 .or. nbCellMast .eq. 0) then
        call utmess('F', 'MESH4_1')
    end if

! - Flags for debug
    meshPairing%debug = verbosity .ge. 2

! - Create main object for pairing
! - 5 is an average of pair by face - Avoid to allocate to large object
    pairDime = 5*nbCellSlav
    call pairAllocate(pairDime, meshPairing)

! - Tolerances
    meshPairing%distRatio = distRatio
    meshPairing%pairTole = pairTole
    meshPairing%areaTole = areaTole

! - Name of output objects
    zonePair = baseName(1:8)//".LISTPAIRS"
    zoneNbPoinInte = baseName(1:8)//".NBPOIN"
    zonePoinInte = baseName(1:8)//".INTERSLPTS"
    call jedetr(zonePair)
    call jedetr(zoneNbPoinInte)
    call jedetr(zonePoinInte)

! - Get space dimension
    call dismoi('DIM_GEOM', mesh, 'MAILLAGE', repi=spaceDime)
    spaceDime = spaceDime
    if (spaceDime .eq. 2) then
        spaceDime = 2
    else
        spaceDime = 3
    end if
    meshPairing%spaceDime = spaceDime

! - Pairing (fast version)
    if (method == PAIR_FAST) then
        call fastPair(mesh, newgeo, mastConxInvName, &
                      mastNeighName, slavNeighName, &
                      nbCellSlav, nbCellMast, nbNodeMast, &
                      listCellSlav, listCellMast, listNodeMast, &
                      meshPairing)
    elseif (method == PAIR_OLD) then
        call aplcpgn(mesh, newgeo, &
                     mastConxInvName, mastNeighName, slavNeighName, &
                     pairTole, distRatio, &
                     nbCellMast, listCellMast, &
                     nbCellSlav, listCellSlav, &
                     listNodeMast, nbNodeMast, &
                     meshPairing)
    elseif (method == PAIR_ROBUST) then
        call robustPair(mesh, newgeo, &
                        nbCellSlav, nbCellMast, &
                        listCellSlav, listCellMast, &
                        meshPairing)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Save results for this zone
    nbPairZone = meshPairing%nbPair
    if (nbPairZone > 0) then
        call wkvect(zonePair, 'G V I', 2*nbPairZone, jvData)
        zi(jvData-1+1:jvData-1+2*nbPairZone) = &
            meshPairing%Pair(1:2*nbPairZone)
        call wkvect(zoneNbPoinInte, 'G V I', nbPairZone, jvData)
        zi(jvData-1+1:jvData-1+nbPairZone) = &
            meshPairing%nbPoinInte(1:nbPairZone)
        call wkvect(zonePoinInte, 'G V R', 2*MAX_NB_INTE*nbPairZone, jvData)
        zr(jvData-1+1:jvData-1+2*MAX_NB_INTE*nbPairZone) = &
            meshPairing%poinInte(1:2*MAX_NB_INTE*nbPairZone)
    end if

! - Clean memory
    call pairDeallocate(meshPairing)
!
end subroutine
