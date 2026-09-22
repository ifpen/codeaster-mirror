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
subroutine exfonc(listFuncActi, ds_algopara, solver, ds_contact, &
                  ds_constitutive, sddyna, nlDynaDamping, &
                  mater, model)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/dismoi.h"
#include "asterfort/exi_thms.h"
#include "asterfort/exisd.h"
#include "asterfort/getvtx.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ndynlo.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: listFuncActi(*)
    character(len=19), intent(in) :: solver, sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    type(NL_DS_Contact), intent(in) :: ds_contact
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=24), intent(in) :: mater, model
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Initializations
!
! Check compatibility of some functionnalities
!
! --------------------------------------------------------------------------------------------------
!
! In  listFuncActi     : list of active functionnalities
! In  ds_algopara      : datastructure for algorithm parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  solver           : datastructure for solver parameters
! In  ds_contact       : datastructure for contact management
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! In  mater            : name of material
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: reac_incr, reac_iter, iexi
    aster_logical :: l_cont, lallv, l_cont_cont, l_cont_disc, lpena, leltc, l_cont_lac, l_iden_rela
    aster_logical :: l_pilo, l_line_search, lmacr, l_unil, l_diri_undead
    aster_logical :: l_vibr_mode, l_buckling, lexpl, l_xfem, lmodim, l_mult_front
    aster_logical :: l_cont_gcp, lpetsc, lamg, limpex, l_matr_distr, lgcpc
    aster_logical :: londe, l_dyna, l_grot_gdep, l_newt_krylov, l_mumps, l_rom
    aster_logical :: l_energy, lproj, lmatdi, lldsp, lResiCompRela, lResiRefeRela
    aster_logical :: l_unil_pena, l_cont_acti, l_hho, l_undead, lDampModal, lthms, limpl
    aster_logical :: l_state_init, l_reuse
    aster_logical :: check_frot, check_rupt, check_endo
    character(len=24) :: typilo, metres, char24
    character(len=16) :: reli_meth, matrix_pred, matrix_corr
    character(len=3) :: mfdet
    character(len=8) :: partit
    character(len=24), pointer :: slvk(:) => null()
    integer(kind=8), pointer :: slvi(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Active functionnalites
    l_xfem = isfonc(listFuncActi, 'XFEM')
    l_cont_cont = isfonc(listFuncActi, 'CONT_CONTINU')
    l_cont_disc = isfonc(listFuncActi, 'CONT_DISCRET')
    l_cont = isfonc(listFuncActi, 'CONTACT')
    l_cont_lac = isfonc(listFuncActi, 'CONT_LAC')
    l_unil = isfonc(listFuncActi, 'LIAISON_UNILATER')
    l_pilo = isfonc(listFuncActi, 'PILOTAGE')
    l_line_search = isfonc(listFuncActi, 'RECH_LINE')
    lmacr = isfonc(listFuncActi, 'MACR_ELEM_STAT')
    l_vibr_mode = isfonc(listFuncActi, 'MODE_VIBR')
    l_buckling = isfonc(listFuncActi, 'CRIT_STAB')
    londe = ndynlo(sddyna, 'ONDE_PLANE')
    l_dyna = ndynlo(sddyna, 'DYNAMIQUE')
    limpl = ndynlo(sddyna, 'IMPLICITE')
    lexpl = isfonc(listFuncActi, 'EXPLICITE')
    l_grot_gdep = isfonc(listFuncActi, 'GD_ROTA')
    lDampModal = nlDynaDamping%lDampModal
    limpex = isfonc(listFuncActi, 'IMPLEX')
    l_newt_krylov = isfonc(listFuncActi, 'NEWTON_KRYLOV')
    l_rom = isfonc(listFuncActi, 'ROM')
    l_energy = isfonc(listFuncActi, 'ENERGIE')
    lproj = isfonc(listFuncActi, 'PROJ_MODAL')
    lmatdi = isfonc(listFuncActi, 'MATR_DISTRIBUEE')
    leltc = isfonc(listFuncActi, 'ELT_CONTACT')
    lResiCompRela = isfonc(listFuncActi, 'RESI_COMP')
    lResiRefeRela = isfonc(listFuncActi, 'RESI_REFE')
    lgcpc = isfonc(listFuncActi, 'GCPC')
    lpetsc = isfonc(listFuncActi, 'PETSC')
    lldsp = isfonc(listFuncActi, 'LDLT_SP')
    l_mumps = isfonc(listFuncActi, 'MUMPS')
    l_mult_front = isfonc(listFuncActi, 'MULT_FRONT')
    l_diri_undead = isfonc(listFuncActi, 'DIRI_UNDEAD')
    l_matr_distr = isfonc(listFuncActi, 'MATR_DISTRIBUEE')
    l_hho = isfonc(listFuncActi, 'HHO')
    l_undead = isfonc(listFuncActi, 'NEUM_UNDEAD') .or. isfonc(listFuncActi, 'DIRI_UNDEAD')
    l_state_init = isfonc(listFuncActi, 'ETAT_INIT')
    l_reuse = isfonc(listFuncActi, 'REUSE')
!
! - Get algorithm parameters
!
    reac_iter = ds_algopara%reac_iter
    reac_incr = ds_algopara%reac_incr
    matrix_pred = ds_algopara%matrix_pred
    matrix_corr = ds_algopara%matrix_corr
    reli_meth = ds_algopara%line_search%method
    check_frot = ds_constitutive%lJoiFrot
    check_rupt = ds_constitutive%lJoiRupt
    check_endo = ds_constitutive%lJoiEndo

!
! - Get solver parameters
!
    call jeveuo(solver//'.SLVK', 'E', vk24=slvk)
    call jeveuo(solver//'.SLVI', 'E', vi=slvi)
    metres = slvk(1)
    lamg = ((slvk(2) .eq. 'ML') .or. (slvk(2) .eq. 'BOOMER'))
!
! - Contact (DISCRETE)
!
    if (l_cont_disc) then
        lmodim = cfdisl(ds_contact%sdcont_defi, 'MODI_MATR_GLOB')
        lallv = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')
        lpena = cfdisl(ds_contact%sdcont_defi, 'CONT_PENA')
        l_cont_gcp = cfdisl(ds_contact%sdcont_defi, 'CONT_GCP')
        l_cont_acti = cfdisl(ds_contact%sdcont_defi, 'CONT_ACTI')
        if (l_pilo) then
            call utmess('F', 'MECANONLINE_43')
        end if
        if (l_line_search .and. (.not. lallv)) then
            call utmess('A', 'MECANONLINE3_89')
        end if
        if (lgcpc .or. lpetsc) then
            if (.not. (lallv .or. lpena .or. l_cont_gcp)) then
                call utmess('F', 'MECANONLINE3_90', sk=metres)
            end if
            if (l_cont_gcp .and. .not. lldsp) then
                call utmess('F', 'MECANONLINE3_88')
            end if
        end if
        if (reac_incr .eq. 0) then
            if (lmodim) then
                call utmess('F', 'CONTACT_88')
            end if
        end if
        if ((l_vibr_mode .or. l_buckling) .and. lmodim) then
            call utmess('F', 'MECANONLINE5_14')
        end if
    end if
!
! - Contact (CONTINUE)
!
    if (l_cont_cont) then
        if (l_pilo) then
            call utmess('F', 'MECANONLINE3_92')
        end if
        if (l_line_search) then
            call utmess('F', 'MECANONLINE3_91')
        end if
        if (lamg) then
            call utmess('F', 'MECANONLINE3_97', sk=slvk(2))
        end if
        if (lpetsc .and. lmatdi) then
            call utmess('F', 'MECANONLINE3_98')
        end if
        if (lDampModal) then
            call utmess('F', 'MECANONLINE3_93')
        end if
    end if

! - Joint elements
    if (matrix_corr .eq. 'ELASTIQUE') then
        if (check_frot .or. check_rupt .or. check_endo) then
            call utmess('F', 'MECANONLINE3_3')
        end if
    end if
!
! - Contact (LAC)
!
    if (l_cont_lac) then
        l_iden_rela = ds_contact%l_iden_rela
        if (l_iden_rela .and. l_mult_front) then
            call utmess('F', 'MECANONLINE3_99')
        end if
        if (l_matr_distr) then
            call utmess('F', 'CONTACT2_19')
        end if
        if ((lpetsc .or. lgcpc) .and. .not. lldsp) then
            call utmess('F', 'MECANONLINE3_87')
        end if
        if (lResiRefeRela) then
            call utmess('F', 'CONTACT2_21')
        end if
    end if
!
! - Contact: excluion CONTACT+DISTRIBUTION/MODEL AUTRE QUE CENTRALISE (SDNV105C en // issue25915)
!
    if (l_cont) then
        if (limpl) then

            call dismoi('PARTITION', model, 'MODELE', repk=partit)
            call exisd('PARTITION', partit, iexi)
            if (iexi .ne. 0) then
                call utmess('F', 'CONTACT3_46')
            end if
        end if
    end if
!
! - Unilateral link
!
    if (l_unil) then
        l_unil_pena = cfdisl(ds_contact%sdcont_defi, 'UNIL_PENA')
        if (l_unil_pena) then
            lmodim = .true.
            if (reac_incr .eq. 0) then
                if (lmodim) then
                    call utmess('F', 'CONTACT_88')
                end if
            end if
        end if
        if (l_pilo) then
            call utmess('F', 'MECANONLINE3_94')
        end if
        if (l_line_search) then
            call utmess('A', 'MECANONLINE3_95')
        end if
        if (lgcpc .or. lpetsc) then
            call utmess('F', 'MECANONLINE3_96', sk=slvk(1))
        end if
    end if
!
! - Dirichlet undead loads
!
    if (l_diri_undead) then
        if (l_pilo) then
            call utmess('F', 'MECANONLINE5_42')
        end if
        if (l_line_search) then
            call utmess('F', 'MECANONLINE5_39')
        end if
        if (l_dyna) then
            call utmess('F', 'MECANONLINE5_40')
        end if
        if (reac_iter .ne. 1) then
            call utmess('F', 'MECANONLINE5_41')
        end if
    end if
!
! - Post-treatment (buckling, ...)
!
    if (l_vibr_mode .or. l_buckling) then
        if (lgcpc .or. lpetsc) then
            call utmess('F', 'FACTOR_52', sk=slvk(1))
        end if
        if (leltc) then
            call utmess('F', 'MECANONLINE5_3')
        end if
    end if
!
! - Explicit solver
!
    if (lexpl) then
        if (l_cont) then
            call utmess('F', 'MECANONLINE5_22')
        end if
        if (l_unil) then
            call utmess('F', 'MECANONLINE5_23')
        end if
        if (l_grot_gdep) then
            call utmess('A', 'MECANONLINE5_24')
        end if
    end if
!
! - Dynamic
!
    if (l_dyna) then
        if (lResiCompRela) then
            call utmess('F', 'MECANONLINE5_53')
        end if
        if (l_pilo) then
            call utmess('F', 'MECANONLINE5_25')
        end if
        if (l_xfem) then
            call utmess('F', 'MECANONLINE5_28')
        end if
        if (limpex) then
            call utmess('F', 'MECANONLINE5_33')
        end if
        char24 = ''
        lthms = exi_thms(model, .true._1, char24, 0)
        if (lthms) then
            call utmess('F', 'MECANONLINE5_16')
        end if
    end if
!
! - Continuation methods (PILOTAGE)
!
    if (l_pilo) then
        call getvtx('PILOTAGE', 'TYPE', iocc=1, scal=typilo)
        if (l_line_search) then
            if (typilo .eq. 'DDL_IMPO') then
                call utmess('F', 'MECANONLINE5_34')
            end if
        end if
        if ((matrix_pred .eq. 'DEPL_CALCULE') .or. (matrix_pred .eq. 'EXTRAPOLE')) then
            call utmess('F', 'MECANONLINE5_36')
        end if
        call dismoi('VARC_F_INST', mater, 'CHAM_MATER', repk=mfdet)
        if (mfdet .eq. 'OUI') then
            call utmess('F', 'CALCULEL2_58', nk=1, valk=mater(1:8))
        end if
    end if
    if (l_line_search) then
        if ((reli_meth .eq. 'PILOTAGE') .and. (.not. l_pilo)) then
            call utmess('F', 'MECANONLINE5_35')
        end if
    end if
!
! - NEWTON_KRYLOV
!
    if (l_newt_krylov) then
        if (l_pilo) then
            call utmess('F', 'MECANONLINE5_48')
        end if
        if ((.not. lgcpc) .and. (.not. lpetsc)) then
            call utmess('F', 'MECANONLINE5_51')
        end if
    end if
!
! - ROM
!
    if (l_rom) then
        if (l_pilo) then
            call utmess('F', 'ROM5_69')
        end if
        if (l_line_search) then
            call utmess('F', 'ROM5_34')
        end if
        if (l_dyna) then
            call utmess('F', 'ROM5_70')
        end if
        if (l_cont) then
            call utmess('F', 'ROM5_71')
        end if
    end if
!
! - Energy
!
    if (l_energy) then
        if (lproj) then
            call utmess('F', 'MECANONLINE5_6')
        end if
        if (lmatdi) then
            call utmess('F', 'MECANONLINE5_8')
        end if
        if (leltc) then
            call utmess('F', 'MECANONLINE5_15')
        end if
    end if
!
! --- SI ON A BESOIN DE FACTORISER SIMULTANEMENT DEUX MATRICES AVEC LE SOLVEUR MUMPS ON LUI
!     SIGNALE AFIN QU'IL OPTIMISE AU MIEUX LA MEMOIRE POUR CHACUNES D'ELLES.
!     CE N'EST VRAIMENT UTILE QUE SI SOLVEUR/GESTION_MEMOIRE='AUTO'.
!
    if (l_mumps) then
        if (l_vibr_mode .or. l_buckling) then
            ASSERT(slvi(6) .ge. 0)
            slvi(6) = 2
        end if
    end if
!
! - HHO
!
    if (l_hho) then
        if (l_cont) then
            call utmess('F', 'MECANONLINE5_60')
        end if
        if (l_unil) then
            call utmess('F', 'MECANONLINE5_61')
        end if
        if (l_pilo) then
            call getvtx('PILOTAGE', 'TYPE', iocc=1, scal=typilo)
            if (typilo .ne. 'DEFORMATION' .and. typilo .ne. "PRED_ELAS") then
                call utmess('F', 'MECANONLINE5_62')
            end if
        end if
        if (l_rom) then
            call utmess('F', 'MECANONLINE5_65')
        end if
        if (l_xfem) then
            call utmess('F', 'MECANONLINE5_66')
        end if
        if (lmacr) then
            call utmess('F', 'MECANONLINE5_67')
        end if
        if (l_vibr_mode .or. l_buckling) then
            call utmess('F', 'MECANONLINE5_68')
        end if
        if (lResiCompRela) then
            call utmess('F', 'MECANONLINE5_73')
        end if
    end if
!
    call jedema()
!
end subroutine
