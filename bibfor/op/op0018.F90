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
subroutine op0018()
!
    use HHO_precalc_module, only: hhoPreCalc
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/adalig.h"
#include "asterfort/ajlipa.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/asmpi_any.h"
#include "asterfort/assert.h"
#include "asterfort/cetucr.h"
#include "asterfort/cormgi.h"
#include "asterfort/crevge.h"
#include "asterfort/dismoi.h"
#include "asterfort/fetcrf.h"
#include "asterfort/fetskp.h"
#include "asterfort/getelem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/gnoms3.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/initel.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/modelCheck.h"
#include "asterfort/modelGetFEType.h"
#include "asterfort/modelPrint.h"
#include "asterfort/ssafmo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
! COMMAND:  AFFE_MODELE
!
! --------------------------------------------------------------------------------------------------
!
! REMARQUES ET RESTRICTIONS D UTILISATION:
!       LES SEULES VERIFICATIONS FAITES ( FAUX=EXIT ), PORTENT SUR:
!       - L AFFECTATION D ELEMENTS FINIS A TOUTES LES MAILLES DEMANDEES
!       - L AFFECTATION D ELEMENTS FINIS SUR UNE MAILLE AU MOINS
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: dim_topo_curr, dim_topo_init
    integer(kind=8) :: ifm, niv
    character(len=8) :: mesh, model, partsd
    character(len=8) :: name_elem, z_quasi_zero, methode
    character(len=16) :: k16dummy, name_type_geom, repk, valk(2)
    character(len=16) :: phenom, phenomRead, modeli, list_modelisa(10), keywordfact, modeli_in
    character(len=19) :: ligrel
    character(len=24) :: kdis
    character(len=32) :: phemod
    character(len=24), parameter :: list_elem = '&&OP0018.LIST_ELEM'
    integer(kind=8), pointer :: p_list_elem(:) => null()
    integer(kind=8) :: nb_elem
    aster_logical :: l_elem, l_grandeur_cara, lParallelMesh
    aster_logical :: l_calc_rigi, l_need_neigh, has_hho
    aster_logical :: lCheckJacobian, lCheckFSINorms, lCheckPlaneity, lCheckFSIFormulation
    integer(kind=8) :: ielem, iaffe
    integer(kind=8) :: vali(4), ico, idx_modelisa
    integer(kind=8), pointer :: p_cata_dim(:) => null()
    integer(kind=8), pointer :: p_cata_model(:) => null()
    character(len=24) :: mesh_type_geom
    integer(kind=8), pointer :: p_mesh_type_geom(:) => null()
    integer(kind=8), pointer :: p_wk_mail1(:) => null()
    integer(kind=8), pointer :: p_wk_mail2(:) => null()
    character(len=24) :: model_liel
    integer(kind=8), pointer :: p_model_liel(:) => null()
    character(len=24) :: model_maille
    integer(kind=8), pointer :: p_model_maille(:) => null()
    character(len=8), pointer :: p_model_lgrf(:) => null()
    integer(kind=8), pointer :: p_model_nbno(:) => null()
    integer(kind=8) :: lont_liel, nb_grel, nb_elem_affe, nb_mesh_elem, nb_elem_affe_min
    integer(kind=8) :: nb_elem_naffe, nbproc, nbpart
    integer(kind=8) :: nb_affe, nb_affe_ss, nbocc, n1
    integer(kind=8) :: long_grel, nb_modelisa, nume_type_poi1, nume_grel
    integer(kind=8) :: nume_elem, idx_in_grel, nume_type_model, nume_type_geom
    mpi_int :: mrank, msize
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)

! - Initializations
    l_elem = .false.

! - Get command parameters
    call getres(model, k16dummy, k16dummy)
    ligrel = model//'.MODELE'

! - Get mesh
    call getvid(' ', 'MAILLAGE', scal=mesh)
    call dismoi('PARALLEL_MESH', mesh, 'MAILLAGE', repk=repk)
    lParallelMesh = repk .eq. 'OUI'

! - Check jacobians
    call getvtx(' ', 'VERI_JACOBIEN', scal=repk)
    lCheckJacobian = repk .eq. 'OUI'

! - Check FSI normals
    call getvtx(' ', 'VERI_NORM_IFS', scal=repk)
    lCheckFSINorms = repk .eq. 'OUI'

! - Check FSI Formulation
    call getvtx(' ', 'VERI_FORMULATION_IFS', scal=repk)
    lCheckFSIFormulation = repk .eq. 'OUI'

! - Check planeity
    call getvtx(' ', 'VERI_PLAN', scal=repk)
    lCheckPlaneity = repk .eq. 'OUI'

! - Grandeurs caracteristiques
    keywordfact = 'GRANDEUR_CARA'
    call getfac(keywordfact, nbocc)
    l_grandeur_cara = nbocc .gt. 0

! - AFFE_SOUS_STRUC
    keywordfact = 'AFFE_SOUS_STRUC'
    call getfac(keywordfact, nb_affe_ss)

! - AFFE
    keywordfact = 'AFFE'
    call getfac(keywordfact, nb_affe)

! - Access to catalog
    call jeveuo('&CATA.TM.TMDIM', 'L', vi=p_cata_dim)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), nume_type_poi1)

! - Common definition for model SD
    call wkvect(ligrel//'.LGRF', 'G V K8', 4, vk8=p_model_lgrf)
    call wkvect(ligrel//'.NBNO', 'G V I', 1, vi=p_model_nbno)
    p_model_lgrf(1) = mesh
    p_model_nbno(1) = 0

! - Get phenomenon
    if (nb_affe .gt. 0) then
        call getvtx('AFFE', 'PHENOMENE', iocc=1, scal=phenom)
    else if (nb_affe_ss .gt. 0) then
        call getvtx('AFFE_SOUS_STRUC', 'PHENOMENE', iocc=1, scal=phenom)
    end if
    call jeecra(ligrel//'.LGRF', 'DOCU', cval=phenom(1:4))
!
    if (nb_affe .ne. 0) then
        keywordfact = 'AFFE'
! ----- Access to mesh objects
        mesh_type_geom = mesh//'.TYPMAIL'
        call jelira(mesh_type_geom, 'LONMAX', nb_mesh_elem)
        call jeveuo(mesh_type_geom, 'L', vi=p_mesh_type_geom)
! ----- Name of objects for model
        model_maille = ligrel//'.TYFE'
        model_liel = ligrel//'.LIEL'
! ----- Create main objects for model
        call wkvect(model_maille, 'G V I', nb_mesh_elem, vi=p_model_maille)
! ----- Working objects
        AS_ALLOCATE(vi=p_wk_mail1, size=nb_mesh_elem)
        AS_ALLOCATE(vi=p_wk_mail2, size=nb_mesh_elem)
! ----- Loop on AFFE keyword
        do iaffe = 1, nb_affe
            nb_elem = 0
            dim_topo_init = -99
            p_wk_mail2(1:nb_mesh_elem) = 0
! --------- Get Only ONE phénomene !
            call getvtx(keywordfact, 'PHENOMENE', iocc=iaffe, scal=phenomRead)
            if (phenomRead .ne. phenom) then
                call utmess('F', 'MODELE1_11')
            end if
! --------- Get modelisation
            call getvtx(keywordfact, 'MODELISATION', iocc=iaffe, nbval=10, vect=list_modelisa, &
                        nbret=nb_modelisa)
            ASSERT(nb_modelisa .eq. 1)
! --------- Get current finite element
            modeli_in = list_modelisa(1)
            call modelGetFEType(iaffe, phenom, modeli_in, idx_modelisa, modeli)
            call jeveuo(jexnum('&CATA.'//phenom, idx_modelisa), 'L', vi=p_cata_model)
            phemod = phenom//modeli
! --------- Get elements
            call jedetr(list_elem)
            call getelem(mesh, keywordfact, iaffe, ' ', list_elem, nb_elem)
            ASSERT(lParallelMesh .or. nb_elem .gt. 0)
! --------- Check dimensions
            call dismoi('DIM_TOPO', phemod, 'PHEN_MODE', repi=dim_topo_curr)
            if (dim_topo_init .eq. -99) then
                dim_topo_init = dim_topo_curr
            else
                if (dim_topo_init .ne. dim_topo_curr) then
                    call utmess('F', 'MODELE1_1')
                end if
            end if
! --------- Check 2D modelisations
            if (modeli(1:4) .eq. 'AXIS' .or. modeli(1:4) .eq. 'PLAN' &
                .or. modeli(2:6) .eq. '_PLAN') then
                call dismoi('Z_QUASI_ZERO', mesh, 'MAILLAGE', repk=z_quasi_zero)
                if (z_quasi_zero .ne. 'OUI') then
                    call utmess('A', 'MODELE1_3')
                end if
            end if
! --------- Affect on elements
            if (nb_elem .ne. 0) then
                l_elem = .true.
                call jeveuo(list_elem, 'L', vi=p_list_elem)
                do ielem = 1, nb_elem
                    nume_elem = p_list_elem(ielem)
                    nume_type_geom = p_mesh_type_geom(nume_elem)
                    if (p_cata_model(nume_type_geom) .gt. 0) then
                        p_model_maille(nume_elem) = p_cata_model(nume_type_geom)
                    end if
                    p_wk_mail1(nume_elem) = 1
                    if (p_cata_dim(nume_type_geom) .eq. dim_topo_init) then
                        p_wk_mail2(nume_elem) = 1
                    end if
                end do
            end if
! --------- ON VERIFIE QU'A CHAQUE OCCURENCE DE AFFE, LES MAILLES
! --------- "PRINCIPALES" (QUI N'ETAIENT PAS DEJA AFFECTEES) ONT BIEN ETE AFFECTEES
! --------- PAR DES ELEMENTS
! --------- (PB DES MODELISATIONS A "TROUS") :
            ico = 0
            do nume_elem = 1, nb_mesh_elem
                if ((p_wk_mail2(nume_elem) .eq. 1) .and. (p_model_maille(nume_elem) .eq. 0)) then
                    ico = ico+1
                end if
            end do
            if (ico .gt. 0) then
                vali(1) = iaffe
                vali(2) = ico
                vali(3) = dim_topo_init
                call utmess('A', 'MODELE1_70', ni=3, vali=vali)
            end if

! --------- Check if user elements have been affected
            nb_elem_naffe = 0
            do ielem = 1, nb_mesh_elem
                nume_elem = ielem
                if (p_wk_mail1(nume_elem) .eq. 1) then
                    if (p_model_maille(nume_elem) .eq. 0) then
                        nb_elem_naffe = nb_elem_naffe+1
                        name_elem = int_to_char8(nume_elem)
                        nume_type_geom = p_mesh_type_geom(nume_elem)
                        call jenuno(jexnum('&CATA.TM.NOMTM', nume_type_geom), name_type_geom)
                        if (niv .eq. 2) then
                            valk(1) = name_elem
                            valk(2) = name_type_geom
                            call utmess('I', 'MODELE1_2', nk=2, valk=valk)
                        end if
                    end if
                end if
            end do
        end do
!
! ----- Count number of GREL
!
        nb_grel = 0
!
! ----- Count number of GREL - Elements
!
        nb_elem_affe = 0
        nume_type_model = 0
        do ielem = 1, nb_mesh_elem
            nume_elem = ielem
            if (p_model_maille(nume_elem) .ne. 0) then
                nb_elem_affe = nb_elem_affe+1
                if (p_model_maille(nume_elem) .ne. nume_type_model) then
                    nume_type_model = p_model_maille(nume_elem)
                    nb_grel = nb_grel+1
                end if
            end if
        end do
!
! ----- Printing informations
!
        if (l_elem) then
            vali(1) = nb_mesh_elem
            vali(2) = nb_elem_affe+nb_elem_naffe
            vali(3) = nb_elem_affe
            call utmess('I', 'MODELE1_4', sk=mesh, ni=3, vali=vali)
        end if
        if (lParallelMesh) then
            nb_elem_affe_min = nb_elem_affe
            call asmpi_comm_vect('MPI_MIN', 'I', sci=nb_elem_affe_min)
            if (nb_elem_affe_min .eq. 0) then
                call utmess('F', 'MODELE1_6', sk=mesh)
            end if
        else
            if (nb_elem_affe .eq. 0) then
                call utmess('F', 'MODELE1_6', sk=mesh)
            end if
        end if
!
! ----- Create LIEL
!
        lont_liel = nb_grel+nb_elem_affe
        call jecrec(model_liel, 'G V I', 'NU', 'CONTIG', 'VARIABLE', nb_grel)
        call jeecra(model_liel, 'LONT', lont_liel)
        call jeveuo(model_liel, 'E', vi=p_model_liel)
!
! ----- Store GREL in LIEL - Elements
!
        nume_type_model = 0
        nume_grel = 0
        long_grel = 0
        idx_in_grel = 0
        do nume_elem = 1, nb_mesh_elem
            if (p_model_maille(nume_elem) .ne. 0) then
!
! ------------- Create new GREL
!
                if (p_model_maille(nume_elem) .ne. nume_type_model .and. nume_type_model &
                    .ne. 0) then
                    nume_grel = nume_grel+1
                    long_grel = long_grel+1
                    idx_in_grel = idx_in_grel+1
                    p_model_liel(idx_in_grel) = nume_type_model
                    call jecroc(jexnum(model_liel, nume_grel))
                    call jeecra(jexnum(model_liel, nume_grel), 'LONMAX', long_grel)
                    long_grel = 0
                end if
!
! -------------- Add element in GREL
!
                long_grel = long_grel+1
                idx_in_grel = idx_in_grel+1
                p_model_liel(idx_in_grel) = nume_elem
                nume_type_model = p_model_maille(nume_elem)
            end if
!
! --------- Last element
!
            if (nume_elem .eq. nb_mesh_elem .and. long_grel .ne. 0) then
                nume_grel = nume_grel+1
                long_grel = long_grel+1
                idx_in_grel = idx_in_grel+1
                p_model_liel(idx_in_grel) = nume_type_model
                call jecroc(jexnum(model_liel, nume_grel))
                call jeecra(jexnum(model_liel, nume_grel), 'LONMAX', long_grel)
            end if
        end do
!
        AS_DEALLOCATE(vi=p_wk_mail1)
        AS_DEALLOCATE(vi=p_wk_mail2)
    end if
!
! - AFFE_SOUS_STRUCT
!
    if (nb_affe_ss .gt. 0) then
        call ssafmo(model)
    end if
!
! - Init elements for this LIGREL
!
    call initel(ligrel, l_calc_rigi)
    if ((.not. l_calc_rigi) .and. (nb_affe .ne. 0)) then
        call utmess('A', 'MODELE1_64', sk=model)
    end if
!
! - Print model information
!
    if (nb_affe .gt. 0) then
        call modelPrint(model)
    end if
!
! - Automatic GREL size adaptation
!
!   -- Danger : on ne peut pas appeler les differentes routines
!      dans n'importe quel ordre :
!        * cormgi doit etre appelee apres adalig
!        * fetcrf doit etre appelee apres initel
    call asmpi_info(rank=mrank, size=msize)
    nbproc = to_aster_int(msize)
    call getvtx('DISTRIBUTION', 'METHODE', iocc=1, scal=kdis, nbret=n1)
    ASSERT(n1 .eq. 1)
    if (nbproc .eq. 1) then
        kdis = 'CENTRALISE'
    end if
    if (lParallelMesh .and. kdis .ne. 'CENTRALISE') then
        call utmess('F', 'MODELE1_99', nk=2, valk=[kdis, mesh])
    end if
    if (kdis .ne. 'CENTRALISE') then
!   name of partition
        partsd = "PART"
        call gnoms3(partsd, 5, 8, ".PRTK")
        p_model_lgrf(2) = partsd
    end if
    if (kdis .eq. 'SOUS_DOMAINE') then
        call getvtx('DISTRIBUTION', 'PARTITIONNEUR', iocc=1, scal=methode, nbret=n1)
        ASSERT(n1 .eq. 1)
        call getvis('DISTRIBUTION', 'NB_SOUS_DOMAINE', iocc=1, scal=nbpart, nbret=n1)
        if (n1 .eq. 0) then
            nbpart = nbproc
        end if
        call fetskp(model, methode, nbpart)
        call fetcrf(model, nbpart)
    end if
    if (kdis .eq. 'SOUS_DOMAINE') then
        call dismoi('PARTITION', ligrel, 'LIGREL', repk=partsd)
        call adalig(ligrel, partsd)
    else
        call adalig(ligrel)
    end if
!
! - Set element/(IGREL,IM) object
!
    call cormgi('G', ligrel)
!
! - Creation de la partition :
!
    call ajlipa(model, 'G', kdis)
!
! - Check model
!
    call modelCheck(model, lCheckJacobian, lCheckFSINorms, lCheckPlaneity, lCheckFSIFormulation)
!
! - Create grandeurs caracteristiques
!
    if (l_grandeur_cara) then
        keywordfact = 'GRANDEUR_CARA'
        call cetucr(keywordfact, model)
    end if
!
! - Need neighbours ?
!
    call dismoi('BESOIN_VOISIN', ligrel, 'LIGREL', repk=repk)
    l_need_neigh = repk .eq. 'OUI'
!
! - Create SD_VOISINAGE if necessary
!
    if (l_need_neigh) then
!       Voisins + methodes MAIL_XXXX ne fonctionne pas
        if (kdis(1:5) .eq. 'MAIL_') call utmess('F', 'MODELE1_13')
        call crevge(ligrel, 'G')
    end if
!
! - Precomputation for HHO
!
    call dismoi('EXI_HHO', ligrel, 'LIGREL', repk=repk)
    has_hho = repk == "OUI"
    has_hho = asmpi_any(has_hho, ASTER_TRUE)
    if (has_hho) then
        call hhoPreCalc(model)
    end if
!
!
    call jedema()
end subroutine
