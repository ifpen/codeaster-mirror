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
subroutine pj4dco(typeSelect, &
                  entity1, entity2, &
                  nbCellSelect1, listCellSelect1, &
                  nbNodeSelect2, listNodeSelect2, &
                  geom1, geom2, corrMesh, &
                  l_dmax, dmax, dala, &
                  listInterc_, nbInterc_)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/inslri.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/pj2dtr.h"
#include "asterfort/pj3dfb.h"
#include "asterfort/pj4dap.h"
#include "asterfort/pjxxut.h"
#include "asterfort/pjloin.h"
#include "asterfort/utimsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*), intent(in) :: typeSelect
    character(len=8), intent(in) :: entity1, entity2
    integer(kind=8), intent(in) :: nbCellSelect1, listCellSelect1(*)
    integer(kind=8), intent(in) ::  nbNodeSelect2, listNodeSelect2(*)
    character(len=*), intent(in) :: geom1, geom2
    character(len=16), intent(in)  :: corrMesh
    aster_logical, intent(inout) :: l_dmax
    real(kind=8), intent(in) :: dmax, dala
    character(len=16), optional, intent(in)  :: listInterc_
    integer(kind=8), optional, intent(in)  :: nbInterc_
!
! --------------------------------------------------------------------------------------------------
!
! Create corresp_2_mailla datastructure
!
! Projection of nodes from entity1 to entity2 when entity1 is 2.5D
!
! --------------------------------------------------------------------------------------------------
!
! In  typeSelect       : type of selection (all or restricted) TOUT or PARTIE
! In  entity1          : name of first entity (model or mesh)
! In  entity2          : name of second entity (model or mesh)
! In  nbCellSelect1    : number of cells to select in first entity if typeSelect = 'PARTIE'
! In  listCellSelect1  : list of cells to select in first entity if typeSelect = 'PARTIE'
! In  nbNodeSelect2    : number of nodes to select in second entity if typeSelect = 'PARTIE'
! In  listNodeSelect2  : list of nodes to select in second entity if typeSelect = 'PARTIE'
! In  geom1            : name of object with coordinates of nodes for entity 1
! In  geom2            : name of object with coordinates of nodes for entity 2
! In  coorMesh         : name of corresp_2_mailla datastructure
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: dbg = ASTER_FALSE
    integer(kind=8), parameter :: spacedim = 3
    character(len=8) :: mesh1, mesh2, nodeName2
    character(len=14), parameter :: boite = '&&PJ4DCO.BOITE'
    character(len=16), parameter :: corrMeshTemp = '&&PJ4DCO.CORRESP'
    character(len=16) :: listInterc
    integer(kind=8) :: nbCellType
    integer(kind=8) :: cellListType(MT_NTYMAX)
    character(len=8) :: cellListCode(MT_NTYMAX)
    integer(kind=8) :: ifm, niv, nbNode1, nbNode2, nbCell1, nbCell2, k
    integer(kind=8) :: nbTria, iatr3
    integer(kind=8) :: iCell1, iNode2
    integer(kind=8) :: iacoo1, iacoo2, ino
    integer(kind=8) :: iabtco, jxxk1, iaconu, iacocf, iacotr
    integer(kind=8) :: ilcnx1
    integer(kind=8) :: iaconb, cellTypeNume, idecal, cellLink1, nbtrou, nbInterc
    aster_logical :: loin, loin2
    real(kind=8) :: dmin, cobary(3)
    integer(kind=8), parameter :: nbmax = 5
    integer(kind=8) :: tino2m(nbmax), nbnod, nbnodm
    real(kind=8) :: tdmin2(nbmax)
    integer(kind=8), pointer :: bt3dlc(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: bt3ddi(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8), pointer :: bt3dnb(:) => null()
    real(kind=8), pointer :: bt3dvr(:) => null()
    integer(kind=8), pointer :: listCell1(:) => null(), listNode2(:) => null()
    integer(kind=8), pointer :: lino_loin(:) => null()
    integer(kind=8), pointer :: vinterc(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)

! - Prepare list of entities
    call pjxxut('2D', typeSelect, &
                entity1, entity2, &
                nbCellSelect1, listCellSelect1, &
                nbNodeSelect2, listNodeSelect2, &
                mesh1, mesh2, &
                nbCellType, cellListType, cellListCode)
    call jeveuo('&&PJXXCO.LIMA1', 'L', vi=listCell1)
    call jeveuo('&&PJXXCO.LINO2', 'L', vi=listNode2)

! - Access to meshes
    call dismoi('NB_NO_MAILLA', mesh1, 'MAILLAGE', repi=nbNode1)
    call dismoi('NB_NO_MAILLA', mesh2, 'MAILLAGE', repi=nbNode2)
    call dismoi('NB_MA_MAILLA', mesh1, 'MAILLAGE', repi=nbCell1)
    call dismoi('NB_MA_MAILLA', mesh2, 'MAILLAGE', repi=nbCell2)
    call jeveuo(mesh1//'.TYPMAIL', 'L', vi=typmail)

!   -- l'objet lino_loin contiendra la liste des noeuds projetes un peu loin
    AS_ALLOCATE(vi=lino_loin, size=nbNode2)

! - Cut cells in TRIA3: number of TRIA3
    nbTria = 0
    do iCell1 = 1, nbCell1
        if (listCell1(iCell1) .ne. 0) then
            cellTypeNume = typmail(iCell1)
            if (cellTypeNume .eq. cellListType(1)) then
                nbTria = nbTria+1
            else if (cellTypeNume .eq. cellListType(2)) then
                nbTria = nbTria+1
            else if (cellTypeNume .eq. cellListType(3)) then
                nbTria = nbTria+1
            else if (cellTypeNume .eq. cellListType(4)) then
                nbTria = nbTria+2
            else if (cellTypeNume .eq. cellListType(5)) then
                nbTria = nbTria+2
            else if (cellTypeNume .eq. cellListType(6)) then
                nbTria = nbTria+2
            else
                call utmess('F', 'PROJECTION4_1')
            end if
        end if
    end do
    if (nbTria .eq. 0) then
        call utmess('F', 'PROJECTION4_55')
    end if

! - Create object to cut cells in TRIA3
!         long(v)=1+4*ntr3
!         v(1) : ntr3(=nombre de tria3)
!         v(1+4(i-1)+1) : numero du 1er  noeud du ieme tria3
!         v(1+4(i-1)+2) : numero du 2eme noeud du ieme tria3
!         v(1+4(i-1)+3) : numero du 3eme noeud du ieme tria3
!         v(1+4(i-1)+4) : numero de la maille mere du ieme tria3
!
    call wkvect('&&PJXXCO.TRIA3', 'V V I', 1+4*nbTria, iatr3)
    zi(iatr3-1+1) = nbTria
    call jeveuo(mesh1//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(mesh1//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
    nbTria = 0
    do iCell1 = 1, nbCell1
        if (listCell1(iCell1) .ne. 0) then
            cellTypeNume = typmail(iCell1)
            if ((cellTypeNume .eq. cellListType(1)) .or. &
                (cellTypeNume .eq. cellListType(2)) .or. &
                (cellTypeNume .eq. cellListType(3))) then
                nbTria = nbTria+1
                zi(iatr3+(nbTria-1)*4+4) = iCell1
                zi(iatr3+(nbTria-1)*4+1) = connex(1+zi(ilcnx1-1+iCell1)-2+1)
                zi(iatr3+(nbTria-1)*4+2) = connex(1+zi(ilcnx1-1+iCell1)-2+2)
                zi(iatr3+(nbTria-1)*4+3) = connex(1+zi(ilcnx1-1+iCell1)-2+3)
            else if ((cellTypeNume .eq. cellListType(4)) .or. &
                     (cellTypeNume .eq. cellListType(5)) .or. &
                     (cellTypeNume .eq. cellListType(6))) then
                nbTria = nbTria+1
                zi(iatr3+(nbTria-1)*4+4) = iCell1
                zi(iatr3+(nbTria-1)*4+1) = connex(1+zi(ilcnx1-1+iCell1)-2+1)
                zi(iatr3+(nbTria-1)*4+2) = connex(1+zi(ilcnx1-1+iCell1)-2+2)
                zi(iatr3+(nbTria-1)*4+3) = connex(1+zi(ilcnx1-1+iCell1)-2+3)
                nbTria = nbTria+1
                zi(iatr3+(nbTria-1)*4+4) = iCell1
                zi(iatr3+(nbTria-1)*4+1) = connex(1+zi(ilcnx1-1+iCell1)-2+1)
                zi(iatr3+(nbTria-1)*4+2) = connex(1+zi(ilcnx1-1+iCell1)-2+3)
                zi(iatr3+(nbTria-1)*4+3) = connex(1+zi(ilcnx1-1+iCell1)-2+4)
            else
                call utmess('F', 'PROJECTION4_1')
            end if
        end if
    end do

! - Get coordinates of nodes
    if (geom1 .eq. ' ') then
        call jeveuo(mesh1//'.COORDO    .VALE', 'L', iacoo1)
    else
        call jeveuo(geom1, 'L', iacoo1)
    end if
    if (geom2 .eq. ' ') then
        call jeveuo(mesh2//'.COORDO    .VALE', 'L', iacoo2)
    else
        call jeveuo(geom2, 'L', iacoo2)
    end if

! - Create boxes
    call pj3dfb(boite, '&&PJXXCO.TRIA3', zr(iacoo1), zr(iacoo2))
    call jeveuo(boite//'.BT3DDI', 'L', vi=bt3ddi)
    call jeveuo(boite//'.BT3DVR', 'L', vr=bt3dvr)
    call jeveuo(boite//'.BT3DNB', 'L', vi=bt3dnb)
    call jeveuo(boite//'.BT3DLC', 'L', vi=bt3dlc)
    call jeveuo(boite//'.BT3DCO', 'L', iabtco)

! Boxes for 3D case
!     boite_3d (k14) ::= record
!      .bt3ddi   : ojb s v i  long=3
!      .bt3dvr   : ojb s v r  long=9
!      .bt3dnb   : ojb s v i  long=nx*ny*nz
!      .bt3dlc   : ojb s v i  long=1+nx*ny*nz
!      .bt3dco   : ojb s v i  long=*
!      .bt3ddi(1) : nx=nombre de boites dans la direction x
!      .bt3ddi(2) : ny=nombre de boites dans la direction y
!      .bt3ddi(3) : nz=nombre de boites dans la direction z
!      .bt3dvr(1) : xmin     .bt3dvr(2) : xmax
!      .bt3dvr(3) : ymin     .bt3dvr(4) : ymax
!      .bt3dvr(5) : zmin     .bt3dvr(6) : zmax
!      .bt3dvr(7) : dx = (xmax-xmin)/nbx
!      .bt3dvr(8) : dy = (ymax-ymin)/nby
!      .bt3dvr(9) : dz = (zmax-zmin)/nbz
!      .bt3dnb    : longueurs des boites
!      .bt3dnb(1) : nombre de tetr4 contenus dans la boite(1,1,1)
!      .bt3dnb(2) : nombre de tetr4 contenus dans la boite(2,1,1)
!      .bt3dnb(3) : ...
!      .bt3dnb(nx*ny*nz) : nombre de tetr4 contenus dans la boite(nx,ny,nz)
!      .bt3dlc    : longueurs cumulees de .bt3dco
!      .bt3dlc(1) : 0
!      .bt3dlc(2) : bt3dlc(1)+nbtr3(boite(1,1))
!      .bt3dlc(3) : bt3dlc(2)+nbtr3(boite(2,1))
!      .bt3dlc(4) : ...
!      .bt3dco    : contenu des boites
!       soit   nbtr3 =nbtr3(boite(p,q,r)=bt3dnb((r-1)*ny*nx+(q-1)*nx+p)
!              debtr3=bt3dlc((r-1)*ny*nx+(q-1)*nx+p)
!        do k=1,nbtr3
!          tr3=.bt3dco(debtr3+k)
!        done
!        TR3 EST LE NUMERO DU KIEME TETR4 DE LA BOITE (P,Q,R)

! - Create temporary datastructure
    call wkvect(corrMeshTemp//'.PJXX_K1', 'V V K24', 5, jxxk1)
    zk24(jxxk1-1+1) = mesh1
    zk24(jxxk1-1+2) = mesh2
    zk24(jxxk1-1+3) = 'COLLOCATION'
    call wkvect(corrMeshTemp//'.PJEF_NB', 'V V I', nbNode2, iaconb)
    call wkvect(corrMeshTemp//'.PJEF_NU', 'V V I', 4*nbNode2, iaconu)
    call wkvect(corrMeshTemp//'.PJEF_CF', 'V V R', 3*nbNode2, iacocf)
    call wkvect(corrMeshTemp//'.PJEF_TR', 'V V I', nbNode2, iacotr)

! - Look for cell cellLink1 in entity1 for each node from entity 2
    idecal = 0
    loin2 = .false.
    nbnod = 0
    nbnodm = 0
    do iNode2 = 1, nbNode2
        if (listNode2(iNode2) .eq. 0) cycle
        call pj4dap(iNode2, zr(iacoo2), zr(iacoo1), zi(iatr3), &
                    cobary, cellLink1, nbtrou, bt3ddi, bt3dvr, &
                    bt3dnb, bt3dlc, zi(iabtco), &
                    l_dmax, dmax, dala, loin, dmin)
        if (loin) then
!           on regarde si le noeud est deja projete par une autre
!           occurrence de VIS_A_VIS
            if (present(nbInterc_)) then
                if (nbInterc_ .ne. 0) then
                    ASSERT(present(listInterc_))
                    call jeveuo(listInterc_, 'L', vi=vinterc)
                    do ino = 1, nbInterc_
                        if (iNode2 .eq. vinterc(ino)) then
                            loin = .false.
                            l_dmax = .true.
                            nbtrou = 0
                            nodeName2 = int_to_char8(iNode2)
                            call utmess('A', 'CALCULEL5_47', si=vinterc(nbInterc_+1), sk=nodeName2)
                            exit
                        end if
                    end do
                end if
            end if
        end if
        if (loin) then
            loin2 = .true.
            nbnodm = nbnodm+1
            lino_loin(nbnodm) = iNode2
        end if
        call inslri(nbmax, nbnod, tdmin2, tino2m, dmin, iNode2)
        if (l_dmax .and. (nbtrou .eq. 0)) then
            zi(iaconb-1+iNode2) = 3
            zi(iacotr-1+iNode2) = 0
            cycle
        end if
        if (nbtrou .eq. 0) then
            nodeName2 = int_to_char8(iNode2)
            call utmess('F', 'PROJECTION4_56', sk=nodeName2)
        end if
        zi(iaconb-1+iNode2) = 3
        zi(iacotr-1+iNode2) = cellLink1
        do k = 1, 3
            zi(iaconu-1+idecal+k) = zi(iatr3+4*(cellLink1-1)+k)
            zr(iacocf-1+idecal+k) = cobary(k)
        end do
        idecal = idecal+zi(iaconb-1+iNode2)
    end do

! - For alarm, see ticket 16186

! - Transform corrMeshTemp in corrMesh (real cells)
    if (present(listInterc_)) then
        ASSERT(present(nbInterc_))
        listInterc = listInterc_
        nbInterc = nbInterc_
        call pj2dtr(corrMeshTemp, corrMesh, &
                    cellListType, cellListCode, &
                    zr(iacoo1), zr(iacoo2), &
                    spacedim, dala, &
                    listInterc, nbInterc)
    else
        listInterc = ' '
        nbInterc = 0
        call pj2dtr(corrMeshTemp, corrMesh, &
                    cellListType, cellListCode, &
                    zr(iacoo1), zr(iacoo2), &
                    spacedim, dala)
    end if

! - Debug
    if (dbg) then
        call utimsd(ifm, 2, ASTER_FALSE, ASTER_TRUE, '&&PJ4DCO', 1, ' ')
        call utimsd(ifm, 2, ASTER_FALSE, ASTER_TRUE, corrMesh, 1, ' ')
    end if

! - Clean
    call detrsd('CORRESP_2_MAILLA', corrMeshTemp)
    call jedetr(boite//'.BT3DDI')
    call jedetr(boite//'.BT3DVR')
    call jedetr(boite//'.BT3DNB')
    call jedetr(boite//'.BT3DLC')
    call jedetr(boite//'.BT3DCO')
    call jedetr('&&PJXXCO.LIMA1')
    call jedetr('&&PJXXCO.LIMA2')
    call jedetr('&&PJXXCO.LINO1')
    call jedetr('&&PJXXCO.LINO2')
    call jedetr('&&PJXXCO.TRIA3')
    AS_DEALLOCATE(vi=lino_loin)
!
    call jedema()
end subroutine
