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
! Module for management of computation of post-processing
!
! ==================================================================================================
!
module postComp_module
! ==================================================================================================
    use FED_module
    use listLoad_module
    use listLoad_type
    use mesh_module
    use postComp_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: setResuPara, initPara, getPara, getInputFields, getNumbering
    public :: compForcNoda, compFextNoda, compMGamma
    public :: getListLoadsUser
    public :: getLigrel, checkOption, initCompRest, deleteCompRest
    public :: getCellSiroElem
    private :: getScalType, getCompLigrel, getMassMatrix, getListLoads
    private :: getListLoadsResu, setListLoads
    private :: initParaNoda
! ==================================================================================================
    private
#include "asterc/isnnem.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asasve.h"
#include "asterfort/ascomb.h"
#include "asterfort/ascova.h"
#include "asterfort/assert.h"
#include "asterfort/calcop.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/dylach.h"
#include "asterfort/exlima.h"
#include "asterfort/getelem.h"
#include "asterfort/gettco.h"
#include "asterfort/gnomsd.h"
#include "asterfort/ischar.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jerazo.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lislec.h"
#include "asterfort/lisnnb.h"
#include "asterfort/mcmult.h"
#include "asterfort/memam2.h"
#include "asterfort/mrmult.h"
#include "asterfort/mtdscr.h"
#include "asterfort/nmdoch.h"
#include "asterfort/ntdoch.h"
#include "asterfort/numecn.h"
#include "asterfort/oriem1.h"
#include "asterfort/pteddl.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rs_get_caraelem.h"
#include "asterfort/rs_get_mate.h"
#include "asterfort/rs_get_model.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmamo.h"
#include "asterfort/utmess.h"
#include "asterfort/vecgme.h"
#include "asterfort/vechme.h"
#include "asterfort/vefnme_cplx.h"
#include "asterfort/vefpme.h"
#include "asterfort/verif_bord.h"
#include "asterfort/vrcins.h"
#include "asterfort/vtcreb.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/zaxpy.h"
#include "jeveux.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! initPara
!
! Init parameters to compute option
!
! In  option            : option to compute
! IO  postComp          : general type for post-processing
!
! --------------------------------------------------------------------------------------------------
    subroutine initPara(optionZ, postComp)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: optionZ
        type(POST_COMP), intent(inout) :: postComp
! ----- Local
        character(len=8) :: answer, modelRefe, caraElemRefe
        character(len=16) :: option
        integer(kind=8) :: codret, numeStore0
        aster_logical :: lRDM
!   ------------------------------------------------------------------------------------------------
!
        option = optionZ(1:16)

! ----- Get first storing index
        numeStore0 = postComp%postCompResu%listStore(1)

! ----- Get reference model
        call rs_get_model(postComp%postCompResu%resultIn, numeStore0, modelRefe, codret)
        if (codret .eq. 4) then
            call utmess('A', 'POSTCOMP1_2')
        end if
        if (codret .lt. 0) then
            call utmess('A', 'POSTCOMP1_3', sk=option)
        end if

! ----- Detect structural elements
        call dismoi('EXI_RDM', modelRefe, 'MODELE', repk=answer)
        lRDM = answer(1:3) .eq. 'OUI'

! ----- Get reference CARA_ELEM
        call rs_get_caraelem(postComp%postCompResu%resultIn, numeStore0, caraElemRefe, codret)
        if (codret .eq. 4) then
            call utmess('A', 'POSTCOMP1_6')
        end if
        if (codret .lt. 0) then
            if (lRDM) then
                call utmess('F', 'POSTCOMP1_14', sk=option)
            end if
        end if

! ----- Get parameters about STRX_ELGA
        call dismoi('EXI_STRX', modelRefe, 'MODELE', repk=answer)
        postComp%postCompPara%lElemStrx = (answer .eq. 'OUI')
        call dismoi('EXI_STR2', modelRefe, 'MODELE', repk=answer)
        postComp%postCompPara%lElemStrxComp = (answer .eq. 'OUI')

! ----- Reduction of computation on GROUP_MA ?
        call getCompLigrel(modelRefe, &
                           postComp%postCompPara%lReduComp, postComp%postCompPara%calcLigrel)

! ----- Get scalar type of field (complex or real)
        call getScalType(postComp%postCompResu, numeStore0, postComp%scalType)
        postComp%lCplx = postComp%scalType .eq. "C"

! ----- Get descriptor for mass matrix
        call getMassMatrix(postComp%postCompResu, &
                           postComp%postCompPara%numeDofRefe, postComp%postCompPara%jvMassMatr)

! ----- Specific parameters for "nodal" post-treatment
        if (postComp%lPostNoda) then
            call initParaNoda(option, modelRefe, numeStore0, postComp)
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initParaNoda
!
! Init parameters to compute option - Nodal options
!
! In  option            : option to compute
! In  modelRefe         : reference model
! In  numeStore0        : index for first storage index
! IO  postComp          : general type for post-processing
!
! --------------------------------------------------------------------------------------------------
    subroutine initParaNoda(optionZ, modelRefeZ, numeStore0, postComp)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: optionZ, modelRefeZ
        integer(kind=8), intent(in) :: numeStore0
        type(POST_COMP), intent(inout) :: postComp
! ----- Local
        character(len=8) :: modelRefe
        character(len=16) :: option
        character(len=24) :: listLoadUser, listLoadResu
        aster_logical :: lLoadsFromUser, lLoadsFromResu, lLoadsHarmo
        complex(kind=8), parameter :: cun = dcmplx(1.D0, 0.D0)
        complex(kind=8), parameter :: cmun = dcmplx(-1.D0, 0.D0)
!   ------------------------------------------------------------------------------------------------
!
        option = optionZ
        modelRefe = modelRefeZ

! ----- Do need loads ?
        postComp%lNeedLoads = option .eq. "REAC_NODA" .or. option .eq. "FEXT_NODA"
        if (postComp%lNeedLoads) then
            call getListLoadsUser(postComp%postCompResu, modelRefe, &
                                  lLoadsFromUser, listLoadUser, lLoadsHarmo)
            call getListLoadsResu(postComp%postCompResu, numeStore0, lLoadsFromResu, listLoadResu)
            call setListLoads(lLoadsFromUser, lLoadsFromResu, listLoadUser, listLoadResu, &
                              lLoadsHarmo, postComp%postCompPara)
            if (.not. postComp%postCompPara%lLoadsAreDefined) then
                call utmess("I", 'POSTCOMP1_15', sk=option)
            end if
        end if

! ----- No damping effect on force
        if ((option .eq. 'REAC_NODA') .and. &
            ((postComp%postCompResu%resultType .eq. 'DYNA_TRANS') .or. &
             (postComp%postCompResu%resultType .eq. 'DYNA_HARMO'))) then
            call utmess('I', 'POSTCOMP1_5')
        end if

! ----- Select options to compute
! ----- REAC_NODA = FORC_NODA - FEXT_NODA
        postComp%postCompNoda%lForcNoda = option .eq. "FORC_NODA" .or. option .eq. "REAC_NODA"
        postComp%postCompNoda%lReacNoda = option .eq. "REAC_NODA"
        postComp%postCompNoda%lFextNoda = option .eq. "REAC_NODA" .or. option .eq. "FEXT_NODA"
        postComp%postCompNoda%lMGamma = option .eq. "REAC_NODA" .or. option .eq. "M_GAMMA"

! ----- Select coefficients
        postComp%postCompNoda%coefFextR = 1.d0
        postComp%postCompNoda%coefFextC = cun
        postComp%postCompNoda%coefMGamR = 1.D0
        postComp%postCompNoda%coefMGamC = cun
        if (postComp%postCompNoda%lReacNoda) then
            postComp%postCompNoda%coefFextR = -1.d0
            postComp%postCompNoda%coefFextC = cmun
            postComp%postCompNoda%coefMGamR = 1.D0
            postComp%postCompNoda%coefMGamC = cun
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setResuPara
!
! Set parameters of results
!
! In  resultIn          : name of datastructure for input results
! In  resultOut         : name of datastructure for output results
! In  resultType        : type of results datastructure
! In  nbStore           : number of storing indexes
! Ptr listStore         : pointer to list of storing indexes
! In  listStoreJv       : JEVEUX name object for list of storing indexes
! IO  postCompResu      : calculation of a post-processing option - Parameters of results
!
! --------------------------------------------------------------------------------------------------
    subroutine setResuPara(resultInZ, resultOutZ, resultTypeZ, &
                           nbStore, listStore, listStoreJvZ, &
                           postCompResu)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: resultInZ, resultOutZ, resultTypeZ
        integer(kind=8), intent(in) :: nbStore
        integer(kind=8), pointer :: listStore(:)
        character(len=*), intent(in) :: listStoreJvZ
        type(POST_COMP_RESU), intent(inout) :: postCompResu
! ----- Locals
        integer(kind=8) :: jvPara
!   ------------------------------------------------------------------------------------------------
!
        postCompResu%resultIn = resultInZ
        postCompResu%resultOut = resultOutZ
        postCompResu%resultType = resultTypeZ
        postCompResu%listStoreJv = listStoreJvZ
        postCompResu%nbStore = nbStore
        postCompResu%listStore => listStore

! ----- Detect transient results
        postCompResu%lTransient = ASTER_FALSE
        if (postCompResu%resultType .eq. 'EVOL_ELAS' .or. &
            postCompResu%resultType .eq. 'EVOL_NOLI' .or. &
            postCompResu%resultType .eq. 'DYNA_TRANS') then
            postCompResu%lTransient = ASTER_TRUE
        end if

! ----- Detect thermal results
        postCompResu%lTher = ASTER_FALSE
        if (postCompResu%resultType .eq. 'EVOL_THER') then
            postCompResu%lTher = ASTER_TRUE
        end if

! ----- Detect mode results
        postCompResu%lMode = ASTER_FALSE
        if (postCompResu%resultType .eq. 'MODE_MECA') then
            call rsadpa(resultInZ, 'L', 1, 'TYPE_MODE', 1, 0, sjv=jvPara)
            postCompResu%lMode = ASTER_TRUE
            postCompResu%modeType = zk16(jvPara)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getScalType
!
! Get scalar type of field (complex or real)
!
! In  postCompResu      : calculation of a post-processing option - Parameters of results
! In  numeStore0        : first storage index
! Out scalType          : scalar type of field (complex or real)
!
! --------------------------------------------------------------------------------------------------
    subroutine getScalType(postCompResu, numeStore0, scalType)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(POST_COMP_RESU), intent(in) :: postCompResu
        integer(kind=8), intent(in) :: numeStore0
        character(len=1), intent(out) :: scalType
! ----- Locals
        character(len=24) :: disp
        integer(kind=8) :: iret
!   ------------------------------------------------------------------------------------------------
!
        scalType = "R"
        if (postCompResu%resultType .eq. 'DYNA_HARMO' .or. &
            postCompResu%resultType .eq. "MODE_MECA") then
            call rsexch(' ', postCompResu%resultIn, 'DEPL', numeStore0, disp, iret)
            if (iret .eq. 0) then
                call jelira(disp(1:19)//'.VALE', 'TYPE', cval=scalType)
            else
                call utmess('F', 'POSTCOMP1_4')
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getCompLigrel
!
! Reduce computation on part of mesh
!
! In  model             : model
! Out lReduComp         : flag when option is not compute on whole model
! Out calcLigrel        : ligrel to compute option
!
! --------------------------------------------------------------------------------------------------
    subroutine getCompLigrel(model, lReduComp, calcLigrel)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: model
        aster_logical, intent(out) :: lReduComp
        character(len=24), intent(out) :: calcLigrel
! ----- Locals
        character(len=24) :: modelLigrel
!   ------------------------------------------------------------------------------------------------
!
        calcLigrel = " "
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
        call exlima(' ', 0, 'V', model, calcLigrel)
        lReduComp = calcLigrel .ne. modelLigrel
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getInputFields
!
! Get input fields to compute option
!
! In  option            : option to compute
! In  numeStore         : current storing index
! In  postComp          : general type for post-processing
!
! --------------------------------------------------------------------------------------------------
    subroutine getInputFields(optionZ, numeStore, postComp)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: optionZ
        integer(kind=8), intent(in) :: numeStore
        type(POST_COMP), intent(inout) :: postComp
! ----- Locals
        integer(kind=8) :: iret, nbEqua
        character(len=2) :: codret
        character(len=16) :: option
!   ------------------------------------------------------------------------------------------------
!
        option = optionZ

! ----- Get displacements
        call rsexch(' ', postComp%postCompResu%resultIn, 'DEPL', &
                    numeStore, postComp%postCompFields%disp, iret)
        if (iret .ne. 0) then
            call utmess('F', 'POSTCOMP1_9', sk=option, si=numeStore)
        end if

! ----- Create zero vector
        call copisd('CHAMP_GD', 'V', postComp%postCompFields%disp, postComp%postCompFields%vectZero)
        call jelira(postComp%postCompFields%vectZero(1:19)//'.VALE', 'LONMAX', nbEqua)
        call jerazo(postComp%postCompFields%vectZero(1:19)//'.VALE', nbEqua, 1)
        postComp%postCompFields%nbEqua = nbEqua

! ----- Get speed (if not present => zero)
        call rsexch(' ', postComp%postCompResu%resultIn, 'VITE', &
                    numeStore, postComp%postCompFields%vite, iret)
        if (iret .eq. 0) then
            postComp%postCompFields%hasVite = ASTER_TRUE
        else
            postComp%postCompFields%vite = postComp%postCompFields%vectZero
        end if

! ----- Get acceleration (if not present => zero)
        call rsexch(' ', postComp%postCompResu%resultIn, 'ACCE', &
                    numeStore, postComp%postCompFields%acce, iret)
        if (iret .eq. 0) then
            postComp%postCompFields%hasAcce = ASTER_TRUE
        else
            postComp%postCompFields%acce = postComp%postCompFields%vectZero
        end if

! ----- Get stresses (SIEF_ELGA) - Compute and save if doesn't exist
        call rsexch(' ', postComp%postCompResu%resultIn, 'SIEF_ELGA', numeStore, &
                    postComp%postCompFields%sigm, iret)
        if (iret .ne. 0) then
            call rsexch(' ', postComp%postCompResu%resultOut, 'SIEF_ELGA', numeStore, &
                        postComp%postCompFields%sigm, iret)
            if (iret .ne. 0) then
! ------------- Compute and save
                if (postComp%postCompResu%resultType .eq. "EVOL_NOLI") then
                    call utmess("F", "CALCCHAMP1_4")
                else
                    call computeField(postComp%postCompResu, "SIEF_ELGA", option)
                    call rsexch(' ', postComp%postCompResu%resultOut, 'SIEF_ELGA', numeStore, &
                                postComp%postCompFields%sigm, iret)
                end if
            end if
        end if

! ----- Get STRX_ELGA - Compute if doesn't exist (and it's possible for this element !)
        if (postComp%postCompPara%lElemStrx) then
            call rsexch(' ', postComp%postCompResu%resultIn, 'STRX_ELGA', numeStore, &
                        postComp%postCompFields%strx, iret)
            if (iret .ne. 0 .and. postComp%postCompPara%lElemStrxComp) then
! ------------- Compute and save
                call computeField(postComp%postCompResu, "STRX_ELGA", option)
                call rsexch(' ', postComp%postCompResu%resultOut, 'STRX_ELGA', numeStore, &
                            postComp%postCompFields%strx, iret)
            end if
        end if

! ----- Prepare external state variables
        call vrcins(postComp%postCompPara%model, postComp%postCompPara%materField, &
                    postComp%postCompPara%caraElem, postComp%postCompPara%time, &
                    postComp%postCompPara%chvarc, codret)

! ----- Get field for non-linear behaviour
        call rsexch(' ', postComp%postCompResu%resultIn, 'COMPORTEMENT', numeStore, &
                    postComp%postCompPara%compor, iret)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! computeField
!
! Compute field
!
! In  postCompResu      : calculation of a post-processing option - Parameters of results
! In  option            : option to compute field
! In  optionOrig        : original option requiring this field
!
! --------------------------------------------------------------------------------------------------
    subroutine computeField(postCompResu, optionZ, optionOrigZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(POST_COMP_RESU), intent(in) :: postCompResu
        character(len=*), intent(in) :: optionZ, optionOrigZ
! ----- Locals
        integer(kind=8) :: cret
        character(len=16) :: option
!   ------------------------------------------------------------------------------------------------
!
        option = optionZ
        call calcop(option, ' ', &
                    postCompResu%resultIn, postCompResu%resultOut, &
                    postCompResu%listStoreJv, postCompResu%nbStore, &
                    postCompResu%resultType, cret, 'V')
        if (cret .ne. 0) then
            call utmess("F", "CALCCHAMP_52", nk=2, valk=[optionOrigZ(1:16), option])
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getPara
!
! Get parameters to compute option
!
! In  option            : option to compute
! In  numeStore         : current storing index
! IO  postComp          : general type for post-processing
!
! --------------------------------------------------------------------------------------------------
    subroutine getPara(optionZ, numeStore, postComp)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: optionZ
        integer(kind=8), intent(in) :: numeStore
        type(POST_COMP), intent(inout) :: postComp
! ----- Local
        character(len=16) :: option
        integer(kind=8) :: codret, jvPara
!   ------------------------------------------------------------------------------------------------
!
        option = optionZ

! ----- Get time
        if (postComp%postCompResu%lTransient) then
            call rsadpa(postComp%postCompResu%resultIn, 'L', 1, 'INST', &
                        numeStore, 0, sjv=jvPara, iStop=0)
            postComp%postCompPara%time = zr(jvPara)
            postComp%postCompPara%hasTime = postComp%postCompPara%time .ne. r8vide()
        end if

! ----- Get frequency
        if (postComp%postCompResu%resultType .eq. 'DYNA_HARMO') then
            call rsadpa(postComp%postCompResu%resultIn, 'L', 1, 'FREQ', &
                        numeStore, 0, sjv=jvPara, iStop=0)
            postComp%postCompPara%freq = zr(jvPara)
            postComp%postCompPara%hasFreq = postComp%postCompPara%freq .ne. r8vide()
        end if

! ----- Get (square) pulsation
        if (postComp%postCompResu%resultType .eq. 'MODE_MECA') then
            call rsadpa(postComp%postCompResu%resultIn, 'L', 1, 'OMEGA2', &
                        numeStore, 0, sjv=jvPara, iStop=0)
            postComp%postCompPara%omega = zr(jvPara)
            postComp%postCompPara%hasOmega = postComp%postCompPara%omega .ne. r8vide()
        end if

! ----- Get Fourier Harmonic
        if (postComp%postCompResu%resultType(1:8) .eq. 'FOURIER_') then
            call rsadpa(postComp%postCompResu%resultIn, 'L', 1, 'NUME_MODE', &
                        numeStore, 0, sjv=jvPara, iStop=0)
            postComp%postCompPara%nh = zi(jvPara)
            if (postComp%postCompPara%nh .eq. isnnem()) then
                call utmess('F', 'POSTCOMP1_18')
            end if
        end if

! ----- Get model
        call rs_get_model(postComp%postCompResu%resultIn, numeStore, &
                          postComp%postCompPara%model, codret)
        if (codret .eq. 4) then
            call utmess('I', 'POSTCOMP1_2')
        end if
        if (codret .lt. 0) then
            call utmess('F', 'POSTCOMP1_3', sk=option)
        end if

! ----- Get elem. characteristics
        call rs_get_caraelem(postComp%postCompResu%resultIn, numeStore, &
                             postComp%postCompPara%caraElem, codret)
        if (codret .eq. 4) then
            call utmess('I', 'POSTCOMP1_6')
        end if

! ----- Get material field
        if (postComp%postCompResu%resultType .ne. 'DYNA_HARMO') then
            call rs_get_mate(postComp%postCompResu%resultIn, numeStore, &
                             postComp%postCompPara%materField, codret)
            if (codret .eq. 4) then
                call utmess('I', 'POSTCOMP1_7')
            end if
            if (codret .lt. 0) then
                call utmess("A", 'POSTCOMP1_8', sk=option)
            end if
            if (postComp%postCompPara%materField .ne. " ") then
                call rcmfmc(postComp%postCompPara%materField, postComp%postCompPara%materCode, &
                            l_ther_=postComp%postCompResu%lTher)
            end if
        end if

! ----- Get loads
        if (postComp%lNeedLoads) then
            if (.not. postComp%postCompPara%lLoadsHarmo) then
                call getListLoads(postComp, numeStore)
                if (.not. postComp%postCompPara%lLoadsAreDefined) then
                    call utmess("I", 'POSTCOMP1_15', sk=option)
                end if
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMassMatrix
!
! Get mass matrix
!
! In  postCompResu      : calculation of a post-processing option - Parameters of results
! Out numeDof           : numbering from matrices
! Out jvMassMatr        : JEVEUX base to mass matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine getMassMatrix(postCompResu, numeDof, jvMassMatr)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(POST_COMP_RESU), intent(in) :: postCompResu
        character(len=24), intent(out) :: numeDof
        integer(kind=8), intent(out) :: jvMassMatr
! ----- Local
        integer(kind=8) :: iret
        character(len=16) :: typeSD
        character(len=19) :: massGene
        character(len=24) :: modalBase, massMatrix
        character(len=24), pointer :: refa(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        massMatrix = " "
        numeDof = " "
        if (postCompResu%resultType .eq. 'MODE_MECA' .or. &
            postCompResu%resultType .eq. 'DYNA_TRANS') then
            call dismoi('NUME_DDL', postCompResu%resultIn, 'RESU_DYNA', repk=numeDof, arret='C')
            call jeexin(postCompResu%resultIn//'           .REFD', iret)
            if (iret .ne. 0) then
                if (iret .ne. 0) then
                    call dismoi('REF_MASS_PREM', postCompResu%resultIn, 'RESU_DYNA', &
                                repk=massMatrix, arret='C')
                    if (massMatrix .ne. ' ') then
                        call mtdscr(massMatrix)
                    end if
                end if
                call dismoi('BASE_MODALE', postCompResu%resultIn, 'RESU_DYNA', repk=modalBase)
                call gettco(modalBase, typeSD)
                if (typeSD(1:16) .eq. 'MAILLAGE_SDASTER') then
                    call dismoi('REF_MASS_PREM', postCompResu%resultIn, 'RESU_DYNA', &
                                repk=massGene, arret='C')
                    if (massGene .ne. ' ') then
                        call gettco(massGene, typeSD)
                        if (typeSD .eq. 'MATR_ASSE_ELIM_R') then
                            call jeveuo(massGene//'.REFA', 'L', vk24=refa)
                            massMatrix = refa(20) (1:19)
                            call mtdscr(massMatrix)
                            call dismoi('NOM_NUME_DDL', massMatrix, 'MATR_ASSE', repk=numeDof)
                        end if
                    end if
                end if
            end if
        else if (postCompResu%resultType .eq. 'DYNA_HARMO') then
            call jeexin(postCompResu%resultIn//'           .REFD', iret)
            if (iret .ne. 0) then
                call dismoi('REF_MASS_PREM', postCompResu%resultIn, 'RESU_DYNA', &
                            repk=massMatrix, arret='C')
                call mtdscr(massMatrix)
            end if
        end if
        if (massMatrix .ne. " ") then
            call jeveuo(massMatrix(1:19)//'.&INT', 'L', jvMassMatr)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getNumbering
!
! Get numbering
!
! In  option            : option to compute
! In  postComp          : general type for post-processing
! Out numeDof           : numbering from matrices
!
! --------------------------------------------------------------------------------------------------
    subroutine getNumbering(optionZ, postComp, numeDof)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: optionZ
        type(POST_COMP), intent(in) :: postComp
        character(len=14), intent(out) :: numeDof
! ----- Local
        character(len=24) :: numeEqua
        character(len=16) :: option
!   ------------------------------------------------------------------------------------------------
!
        numeDof = " "
        option = optionZ

! ----- Get NUME_EQUA
        numeEqua = " "
        if (postComp%postCompResu%resultType .eq. 'MODE_MECA' .or. &
            postComp%postCompResu%resultType .eq. 'DYNA_TRANS') then
            numeEqua = postComp%postCompPara%numeDofRefe(1:14)//'.NUME'
        else
            call numecn(postComp%postCompPara%model, postComp%postCompFields%disp, numeEqua)
        end if
        if (numeEqua .eq. " ") then
            call utmess('F', "POSTCOMP1_10", sk=option)
        end if

! ----- Get NUME_DOF object
        numeDof = numeEqua(1:14)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compForcNoda
!
! Compute FORC_NODA
!
! In  postComp          : general type for post-processing
! In  numeDofIn         : numbering
! In  fieldOut          : JEVEUX name of output field
!
! --------------------------------------------------------------------------------------------------
    subroutine compForcNoda(postComp, numeDofIn, fieldOutZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(POST_COMP), intent(in) :: postComp
        character(len=14), intent(in) :: numeDofIn
        character(len=*), intent(in) :: fieldOutZ
! ----- Local
        character(len=16), parameter :: option = "FORC_NODA"
        character(len=14) :: numeDof
        integer(kind=8) :: iEqua, jfi, jfr, jfo, numeNbEqua
        character(len=24) :: fieldOut, numeEqua
        character(len=24) :: veFnod(2)
        character(len=24) :: vaFnodCR, vaFnodCI, vaFnodR
        real(kind=8), pointer :: vafnodCRVale(:) => null(), vafnodCIVale(:) => null()
        complex(kind=8), pointer :: fieldCVale(:) => null()
        real(kind=8), pointer :: vafnodRVale(:) => null()
        real(kind=8), pointer :: fieldRVale(:) => null()
        blas_int :: b_incx, b_incy, b_n
!   ------------------------------------------------------------------------------------------------
!
        fieldOut = fieldOutZ
        b_n = to_blas_int(postComp%postCompFields%nbEqua)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)

! ----- Access to output field
        if (postComp%lCplx) then
            call jeveuo(fieldOut(1:19)//'.VALE', 'E', vc=fieldCVale)
        else
            call jeveuo(fieldOut(1:19)//'.VALE', 'E', vr=fieldRVale)
        end if

! ----- Compute elementary vectors for FORC_NODA (real or complex)
        veFnod(1) = "&&VEFNOR"
        veFnod(2) = "&&VEFNOC"
        call vefnme_cplx(option, 'V', &
                         postComp%postCompPara%model, postComp%postCompPara%materCode, &
                         postComp%postCompPara%caraElem, postComp%postCompPara%compor, &
                         postComp%postCompPara%timePrev, postComp%postCompPara%time, &
                         postComp%postCompPara%nh, postComp%postCompPara%calcLigrel, &
                         postComp%postCompPara%chvarc, &
                         postComp%postCompFields%sigmPrev, postComp%postCompFields%sigm, &
                         postComp%postCompFields%strx, postComp%postCompFields%disp, &
                         veFnod)

! ----- Get numbering
        call dismoi('NB_EQUA', numeDofIn, 'NUME_DDL', repi=numeNbEqua)
        numeDof = numeDofIn
        if (numeNbEqua .ne. postComp%postCompFields%nbEqua) then
            call dismoi('NUME_EQUA', postComp%postCompFields%disp, 'CHAM_NO', repk=numeEqua)
            numeDof = numeEqua(1:14)
        end if

        call dismoi('NB_EQUA', numeDof, 'NUME_DDL', repi=numeNbEqua)
        if (numeNbEqua .ne. postComp%postCompFields%nbEqua) then
            call utmess("F", "CALCCHAMP1_5", &
                        ni=2, vali=[numeNbEqua, postComp%postCompFields%nbEqua])
        end if

! ----- Pre-assemblying
        vaFnodCR = " "
        vaFnodCI = " "
        vaFnodR = " "
        if (postComp%lCplx) then
            call vtcreb(veFnod(1), 'V', 'R', nume_ddlz=numeDof)
            call asasve(veFnod(1), numeDof, 'R', vaFnodCR)
            call vtcreb(veFnod(2), 'V', 'R', nume_ddlz=numeDof)
            call asasve(veFnod(2), numeDof, 'R', vaFnodCI)
        else
            call asasve(veFnod(1), numeDof, 'R', vaFnodR)
        end if

! ----- Assemblying
        if (postComp%lCplx) then
            call jeveuo(fieldOut(1:19)//'.VALE', 'E', vc=fieldCVale)
            call jeveuo(vaFnodCR, 'L', jfr)
            call jeveuo(zk24(jfr) (1:19)//'.VALE', 'L', vr=vafnodCRVale)
            call jeveuo(vaFnodCI, 'L', jfi)
            call jeveuo(zk24(jfi) (1:19)//'.VALE', 'L', vr=vafnodCIVale)
            do iEqua = 1, postComp%postCompFields%nbEqua
                fieldCVale(iEqua) = dcmplx(vafnodCRVale(iEqua), vafnodCIVale(iEqua))
            end do
        else
            call jeveuo(fieldOut(1:19)//'.VALE', 'E', vr=fieldRVale)
            call jeveuo(vaFnodR, 'L', jfo)
            call jeveuo(zk24(jfo) (1:19)//'.VALE', 'L', vr=vafnodRVale)
            call dcopy(b_n, vafnodRVale, b_incx, fieldRVale, b_incy)
        end if

! ----- Cleaning
        call detrsd('VECT_ELEM', veFnod(1))
        call detrsd('VECT_ELEM', veFnod(2))
        call jedetr(vaFnodCR(1:8))
        call jedetr(vaFnodCI(1:8))
        call jedetr(vaFnodR(1:8))
        call jedetr(veFnod(1) (1:8)//'           .REFE')
        call jedetr(veFnod(2) (1:8)//'           .REFE')
        call jedetr(veFnod(1) (1:8)//'           .DESC')
        call jedetr(veFnod(2) (1:8)//'           .DESC')
        call jedetr(veFnod(1) (1:8)//'           .VALE')
        call jedetr(veFnod(2) (1:8)//'           .VALE')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compFextNoda
!
! Compute FEXT_NODA
!
! In  postComp          : general type for post-processing
! In  optionUser        : option required by user
! In  numeDof           : numbering
! In  fieldOut          : JEVEUX name of output field
!
! --------------------------------------------------------------------------------------------------
    subroutine compFextNoda(postComp, optionUserZ, numeDof, fieldOutZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(POST_COMP), intent(in) :: postComp
        character(len=*), intent(in) :: optionUserZ
        character(len=14), intent(in) :: numeDof
        character(len=*), intent(in) :: fieldOutZ
! ----- Local
        character(len=16) :: optionUser
        character(len=24) :: fieldOut, loadInfoJv, loadNameJv, loadFuncJv
        character(len=24) :: vechmp, vachmp, cnchmpR
        character(len=24) :: vecgmp, vacgmp, cncgmp
        character(len=24) :: vefpip, vafpip, cnfpip
        character(len=8), parameter :: paraName = "FREQ"
        character(len=19), parameter :: vediri = '&&VEDIRI', veneum = '&&VENEUM'
        character(len=19), parameter :: vevoch = '&&VEVOCH', vassec = '&&VASSEC'
        character(len=19), parameter :: cnchmpC = '&&CCFNRN.CHARGE'
        real(kind=8) :: timePara(3)
        real(kind=8) :: coefFextR
        complex(kind=8) :: coefFextC
        blas_int :: b_incx, b_incy, b_n
        complex(kind=8), pointer :: fieldCVale(:) => null()
        real(kind=8), pointer :: fieldRVale(:) => null()
        complex(kind=8), pointer :: cnchmpCVale(:) => null()
        real(kind=8), pointer :: cnchmpRVale(:) => null()
        real(kind=8), pointer :: cncgmpVale(:) => null()
        real(kind=8), pointer :: cnfpipVale(:) => null()
        integer(kind=8) :: nbEqua
        character(len=1) :: stopCalc
!   ------------------------------------------------------------------------------------------------
!
        fieldOut = fieldOutZ
        optionUser = optionUserZ
        nbEqua = postComp%postCompFields%nbEqua
        b_n = to_blas_int(nbEqua)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        coefFextC = postComp%postCompNoda%coefFextC
        coefFextR = postComp%postCompNoda%coefFextR

! ----- Access to output field
        if (postComp%lCplx) then
            call jeveuo(fieldOut(1:19)//'.VALE', 'E', vc=fieldCVale)
        else
            call jeveuo(fieldOut(1:19)//'.VALE', 'E', vr=fieldRVale)
        end if

! ----- Alarm if no skin cells
        if (postComp%postCompPara%lReduComp) then
            call verif_bord(postComp%postCompPara%model, postComp%postCompPara%calcLigrel)
            stopCalc = 'C'
        else
            stopCalc = 'S'
        end if

! ----- Complex case
        if (postComp%lCplx) then
            if (postComp%postCompPara%lReduComp) then
                call utmess('F', 'POSTCOMP1_13', sk=optionUser)
            end if
! --------- Compute elementary vectors and pre-assemblying
            call dylach(postComp%postCompPara%model, postComp%postCompPara%materField, &
                        postComp%postCompPara%materCode, postComp%postCompPara%caraElem, &
                        postComp%postCompPara%listLoad, &
                        numeDof, vediri, veneum, vevoch, vassec)

! --------- Create output field
            call vtcreb(cnchmpC, 'V', 'C', nume_ddlz=numeDof)

! --------- Assemblying
            call ascomb(postComp%postCompPara%listLoad, veneum, 'C', &
                        paraName, postComp%postCompPara%freq, &
                        cnchmpC)
        end if

! ----- Real case
        vechmp = " "
        vecgmp = " "
        vefpip = " "
        vachmp = " "
        vacgmp = " "
        vafpip = " "
        cnchmpR = " "
        cncgmp = " "
        cnfpip = " "
        if (.not. postComp%lCplx) then
            loadNameJv = postComp%postCompPara%listLoad(1:19)//".LCHA"
            loadInfoJv = postComp%postCompPara%listLoad(1:19)//".INFC"
            loadFuncJv = postComp%postCompPara%listLoad(1:19)//".FCHA"
            timePara = 0.d0
            timePara(1) = postComp%postCompPara%time

! --------- Compute elementary vectors and assemblying (dead loads)
            call vechme(stopCalc, &
                        postComp%postCompPara%model, postComp%postCompPara%caraElem, &
                        postComp%postCompPara%materField, postComp%postCompPara%materCode, &
                        loadNameJv, loadInfoJv, &
                        timePara, &
                        vechmp, varcCurrZ_=postComp%postCompPara%chvarc, &
                        ligrelCalcZ_=postComp%postCompPara%calcLigrel, &
                        nharm_=postComp%postCompPara%nh)
            call asasve(vechmp, numeDof, 'R', vachmp)
            call ascova('D', vachmp, loadFuncJv, 'INST', postComp%postCompPara%time, &
                        'R', cnchmpR)

! --------- Compute elementary vectors and assemblying (undead loads)
            call vecgme('S', &
                        postComp%postCompPara%model, postComp%postCompPara%caraElem, &
                        postComp%postCompPara%materField, postComp%postCompPara%materCode, &
                        postComp%postCompPara%compor, &
                        loadNameJv, loadInfoJv, &
                        postComp%postCompPara%time, postComp%postCompPara%time, &
                        postComp%postCompFields%disp, postComp%postCompFields%vectZero, &
                        postComp%postCompFields%vite, postComp%postCompFields%acce, &
                        postComp%postCompFields%strx, &
                        vecgmp, &
                        ligrelCalcZ_=postComp%postCompPara%calcLigrel)
            call asasve(vecgmp, numeDof, 'R', vacgmp)
            call ascova('D', vacgmp, loadFuncJv, 'INST', postComp%postCompPara%time, &
                        'R', cncgmp)

! --------- Compute elementary vectors and assemblying ("PILOTAGE")
            if (postComp%postCompResu%resultType .eq. 'EVOL_NOLI') then
                call vefpme('S', &
                            postComp%postCompPara%model, postComp%postCompPara%caraElem, &
                            postComp%postCompPara%materField, postComp%postCompPara%materCode, &
                            loadNameJv, loadInfoJv, &
                            timePara, &
                            postComp%postCompFields%disp, postComp%postCompFields%vectZero, &
                            postComp%postCompFields%vectZero, &
                            vefpip, postComp%postCompPara%calcLigrel)
                call asasve(vefpip, numeDof, 'R', vafpip)
                call ascova('D', vafpip, loadFuncJv, 'INST', postComp%postCompPara%time, &
                            'R', cnfpip)
            end if
        end if

! ----- Access to load fields
        if (postComp%lCplx) then
            call jeveuo(cnchmpC(1:19)//'.VALE', 'L', vc=cnchmpCVale)
        else
            call jeveuo(cnchmpR(1:19)//'.VALE', 'L', vr=cnchmpRVale)
            call jeveuo(cncgmp(1:19)//'.VALE', 'L', vr=cncgmpVale)
            if (postComp%postCompResu%resultType .eq. 'EVOL_NOLI') then
                call jeveuo(cnfpip(1:19)//'.VALE', 'L', vr=cnfpipVale)
            end if
        end if

! ----- Combine
        if (postComp%lCplx) then
            call jeveuo(fieldOut(1:19)//'.VALE', 'E', vc=fieldCVale)
            call zaxpy(b_n, coefFextC, cnchmpCVale, b_incx, fieldCVale, b_incy)
        else
            call jeveuo(fieldOut(1:19)//'.VALE', 'E', vr=fieldRVale)
            call daxpy(b_n, coefFextR, cnchmpRVale, b_incx, fieldRVale, b_incy)
            call daxpy(b_n, coefFextR, cncgmpVale, b_incx, fieldRVale, b_incy)
            if (postComp%postCompResu%resultType .eq. 'EVOL_NOLI') then
                call daxpy(b_n, coefFextR*postComp%postCompPara%eta, cnfpipVale, b_incx, &
                           fieldRVale, b_incy)
            end if
        end if

! ----- Cleaning
        call detrsd('VECT_ELEM', veneum)
        call detrsd('VECT_ELEM', vechmp)
        call detrsd('VECT_ELEM', vecgmp)
        call detrsd('VECT_ELEM', vefpip)
        call jedetr(vachmp(1:8)//'.ASCOVA')
        call jedetr(vacgmp(1:8)//'.ASCOVA')
        call jedetr(vafpip(1:8)//'.ASCOVA')
        call detrsd('CHAMP_GD', cnchmpR)
        call detrsd('CHAMP_GD', cnchmpC)
        call detrsd('CHAMP_GD', cncgmp)
        call detrsd('CHAMP_GD', cnfpip)
        call detrsd('CHAMP_GD', cnchmpR(1:8)//'.ASCOVA')
        call detrsd('CHAMP_GD', cncgmp(1:8)//'.ASCOVA')
        call detrsd('CHAMP_GD', cnfpip(1:8)//'.ASCOVA')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compMGamma
!
! Compute M_GAMMA
!
! In  postComp          : general type for post-processing
! In  optionUser        : option required by user
! In  numeStore         : current storing index
! In  numeDof           : numbering
! In  fieldOut          : JEVEUX name of output field
!
! --------------------------------------------------------------------------------------------------
    subroutine compMGamma(postComp, optionUserZ, numeStore, numeDof, fieldOutZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(POST_COMP), intent(in) :: postComp
        character(len=*), intent(in) :: optionUserZ
        integer(kind=8), intent(in) :: numeStore
        character(len=14), intent(in) :: numeDof
        character(len=*), intent(in) :: fieldOutZ
! ----- Locals
        character(len=16) :: optionUser
        character(len=24) :: fieldOut
        character(len=16) :: modeType, resultType
        integer(kind=8) :: nbEqua, jvPara, jref, iComp, iEqua, numeSeleEqua
        real(kind=8), pointer :: dispVale(:) => null(), vectWorkR(:) => null()
        real(kind=8), pointer :: acceValeR(:) => null()
        complex(kind=8), pointer :: acceValeC(:) => null(), vectWorkC(:) => null()
        real(kind=8), pointer :: fieldRVale(:) => null()
        complex(kind=8), pointer :: fieldCVale(:) => null()
        real(kind=8) :: omega2
        real(kind=8) :: coefMGamR
        complex(kind=8) :: coefMGamC
        blas_int :: b_incx, b_incy, b_n
        integer(kind=8), pointer :: seleDof(:) => null()
        real(kind=8), pointer :: valeDof(:) => null()
        integer(kind=8), parameter :: nbComp = 3
        character(len=8), parameter :: compName(nbComp) = (/'DX', 'DY', 'DZ'/)
        real(kind=8) :: coefVale(nbComp)
        character(len=24) :: vemgam, vamgam
        real(kind=8), pointer :: cnmgamVale(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        fieldOut = fieldOutZ
        optionUser = optionUserZ
        resultType = postComp%postCompResu%resultType
        modeType = postComp%postCompResu%modeType
        nbEqua = postComp%postCompFields%nbEqua
        b_n = to_blas_int(nbEqua)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        coefMGamC = postComp%postCompNoda%coefMGamC
        coefMGamR = postComp%postCompNoda%coefMGamR

! ----- Access to output field
        if (postComp%lCplx) then
            call jeveuo(fieldOut(1:19)//'.VALE', 'E', vc=fieldCVale)
        else
            call jeveuo(fieldOut(1:19)//'.VALE', 'E', vr=fieldRVale)
        end if

! ----- Working vector
        if (postComp%lCplx) then
            AS_ALLOCATE(vc=vectWorkC, size=nbEqua)
        else
            AS_ALLOCATE(vr=vectWorkR, size=nbEqua)
        end if

        if (resultType .eq. 'MODE_MECA' .and. modeType(1:8) .eq. 'MODE_DYN' .and. &
            .not. postComp%lCplx) then

! --------- Get pulsation
            omega2 = postComp%postCompPara%omega
            if (.not. postComp%postCompPara%hasOmega) then
                call utmess('F', 'POSTCOMP1_19', sk=optionUser, si=numeStore)
            end if

! --------- Access to displacement
            call jeveuo(postComp%postCompFields%disp(1:19)//'.VALE', 'L', vr=dispVale)
            if (postComp%postCompPara%jvMassMatr .eq. 0) then
                call utmess('F', 'POSTCOMP1_16', sk=optionUser)
            end if

! --------- Product omega² * disp * Mass
            call mrmult('ZERO', postComp%postCompPara%jvMassMatr, dispVale, vectWorkR, 1, .true._1)
            call daxpy(b_n, -coefMGamR*omega2, vectWorkR, b_incx, fieldRVale, b_incy)

        elseif (resultType .eq. 'MODE_MECA' .and. modeType(1:8) .eq. 'MODE_STA' &
                .and. .not. postComp%lCplx) then
            call rsadpa(postComp%postCompResu%resultIn, 'L', 1, 'TYPE_DEFO', numeStore, &
                        0, sjv=jvPara)
            if (zk16(jvPara) (1:9) .eq. 'FORC_IMPO') then
                call rsadpa(postComp%postCompResu%resultIn, 'L', 1, 'NUME_DDL', numeStore, &
                            0, sjv=jvPara)
                numeSeleEqua = zi(jvPara)
                fieldRVale(numeSeleEqua) = fieldRVale(numeSeleEqua)-1.d0

            else if (zk16(jvPara) (1:9) .eq. 'ACCE_IMPO') then
                if (postComp%postCompPara%jvMassMatr .eq. 0) then
                    call utmess('F', 'POSTCOMP1_16', sk=optionUser)
                end if

! ------------- Get coefficients
                coefVale = 0.d0
                call rsadpa(postComp%postCompResu%resultIn, 'L', 1, 'COEF_X', numeStore, &
                            0, sjv=jvPara)
                coefVale(1) = zr(jvPara)
                call rsadpa(postComp%postCompResu%resultIn, 'L', 1, 'COEF_Y', numeStore, &
                            0, sjv=jvPara)
                coefVale(2) = zr(jvPara)
                call rsadpa(postComp%postCompResu%resultIn, 'L', 1, 'COEF_Z', numeStore, &
                            0, sjv=jvPara)
                coefVale(3) = zr(jvPara)

! ------------- Select DOF DX/DY/DZ in numbering
                AS_ALLOCATE(vi=seleDof, size=3*nbEqua)
                AS_ALLOCATE(vr=valeDof, size=nbEqua)
                call pteddl('NUME_DDL', numeDof, nbComp, compName, nbEqua, &
                            tabl_equa=seleDof)
                do iComp = 1, nbComp
                    do iEqua = 1, nbEqua
                        valeDof(iEqua) = valeDof(iEqua)+ &
                                         seleDof(nbEqua*(iComp-1)+iEqua)*coefVale(iComp)
                    end do
                end do

! ------------- Product M.Gamma
                call mrmult('ZERO', postComp%postCompPara%jvMassMatr, valeDof, vectWorkR, 1, &
                            .true._1)
                call daxpy(b_n, -1.d0*coefMGamR, vectWorkR, b_incx, fieldRVale, b_incy)
                AS_DEALLOCATE(vi=seleDof)
                AS_DEALLOCATE(vr=valeDof)

            end if
        else if (resultType .eq. 'DYNA_TRANS') then
            if (postComp%postCompFields%hasAcce) then
                call jeveuo(postComp%postCompFields%acce(1:19)//'.VALE', 'L', vr=acceValeR)
                if (postComp%postCompPara%jvMassMatr .eq. 0) then
                    call utmess('F', 'POSTCOMP1_16', sk=optionUser)
                end if
                call mrmult('ZERO', postComp%postCompPara%jvMassMatr, acceValeR, vectWorkR, 1, &
                            .true._1)
                call daxpy(b_n, coefMGamR, vectWorkR, b_incx, fieldRVale, b_incy)
            else
                call utmess('F', 'POSTCOMP1_17', sk=optionUser, si=numeStore)
            end if

        else if (resultType .eq. 'DYNA_HARMO' .or. &
                 (resultType .eq. 'MODE_MECA' .and. &
                  modeType(1:8) .eq. 'MODE_STA' .and. &
                  postComp%lCplx)) then
            if (postComp%postCompFields%hasAcce) then
                call jeveuo(postComp%postCompFields%acce(1:19)//'.VALE', 'L', vc=acceValeC)
                if (postComp%postCompPara%jvMassMatr .eq. 0) then
                    call utmess('F', 'POSTCOMP1_16', sk=optionUser)
                end if
                call mcmult('ZERO', postComp%postCompPara%jvMassMatr, acceValeC, vectWorkC, 1, &
                            .true._1)
                call zaxpy(b_n, coefMGamC, vectWorkC, b_incx, fieldCVale, b_incy)

            else
                call utmess('F', 'POSTCOMP1_17', sk=optionUser, si=numeStore)
            end if

        else if (resultType .eq. 'EVOL_NOLI') then
            vemgam = "&&CCFNRN.VEMGAM"
            if (postComp%postCompFields%hasAcce) then
                call memam2("M_GAMMA", &
                            postComp%postCompPara%model, postComp%postCompPara%materField, &
                            postComp%postCompPara%materCode, postComp%postCompPara%caraElem, &
                            postComp%postCompPara%compor, postComp%postCompPara%time, &
                            postComp%postCompFields%acce, vemgam, 'V', &
                            postComp%postCompPara%calcLigrel)
                call asasve(vemgam, numeDof, 'R', vamgam)
                call jeveuo(vamgam, 'L', jref)
                call jeveuo(zk24(jref) (1:19)//'.VALE', 'L', vr=cnmgamVale)
                call daxpy(b_n, coefMGamR, cnmgamVale, b_incx, fieldRVale, b_incy)
            end if
            call detrsd('VECT_ELEM', vemgam)
        end if
        if (postComp%lCplx) then
            AS_DEALLOCATE(vc=vectWorkC)
        else
            AS_DEALLOCATE(vr=vectWorkR)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getListLoads
!
! Get list of loads
!
! IO  postComp          : general type for post-processing
! In  numeStore         : current storing index
!
! --------------------------------------------------------------------------------------------------
    subroutine getListLoads(postComp, numeStore)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(POST_COMP), intent(inout) :: postComp
        integer(kind=8), intent(in) :: numeStore
! ----- Local
        character(len=24) :: listLoadResu
        aster_logical :: hasPiloLoads, lPilo1, lPilo2
        aster_logical :: lLoadsFromResu, lLoadsAreDefined
        integer(kind=8) :: jvPara
!   ------------------------------------------------------------------------------------------------
!
        if (.not. postComp%postCompPara%lLoadsFromUser) then
            call getListLoadsResu(postComp%postCompResu, numeStore, lLoadsFromResu, listLoadResu)
            lLoadsAreDefined = lLoadsFromResu
            postComp%postCompPara%lLoadsAreDefined = lLoadsAreDefined
            postComp%postCompPara%listLoad = listLoadResu
            postComp%postCompPara%lLoadsFromUser = ASTER_FALSE
        end if
        hasPiloLoads = ASTER_FALSE
        if (postComp%postCompPara%lLoadsAreDefined) then
            lPilo1 = ischar(postComp%postCompPara%listLoad, 'DIRI', 'PILO')
            lPilo2 = ischar(postComp%postCompPara%listLoad, 'NEUM', 'PILO')
            hasPiloLoads = lPilo1 .or. lPilo2
        end if
        postComp%postCompPara%hasPiloLoads = hasPiloLoads

! ----- Get ETA (PILOTAGE)
        if (postComp%postCompPara%hasPiloLoads) then
            call rsadpa(postComp%postCompResu%resultIn, 'L', 1, 'ETA_PILOTAGE', &
                        numeStore, 0, sjv=jvPara, iStop=0)
            postComp%postCompPara%eta = zr(jvPara)
            postComp%postCompPara%hasEta = postComp%postCompPara%eta .ne. r8vide()
        else
            postComp%postCompPara%eta = r8vide()
            postComp%postCompPara%hasEta = ASTER_FALSE
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getListLoadsUser
!
! Get list of loads from user (command file)
!
! In  postCompResu      : calculation of a post-processing option - Parameters of results
! In  model             : name of model
! Out lLoadsFromUser    : flag when loads are defined by user
! Out listLoadUser      : name of datastructure for list of loads from user
! Out lLoadsHarmo       : flag for loads in harmonic case
!
! --------------------------------------------------------------------------------------------------
    subroutine getListLoadsUser(postCompResu, modelZ, &
                                lLoadsFromUser, listLoadUser, lLoadsHarmo)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(POST_COMP_RESU), intent(in) :: postCompResu
        character(len=*), intent(in) :: modelZ
        aster_logical, intent(out) :: lLoadsFromUser
        character(len=24), intent(out) :: listLoadUser
        aster_logical, intent(out) :: lLoadsHarmo
! ----- Local
        character(len=16), parameter :: loadKeyword = 'EXCIT'
        character(len=1), parameter :: jvBase = "V"
        character(len=8) :: model
        integer(kind=8) :: nbLoadUser
        type(ListLoad_Prep) :: listLoadPrep
!   ------------------------------------------------------------------------------------------------
!
        listLoadUser = "&&CCFNRN.LOADUSER"
        lLoadsFromUser = ASTER_FALSE
        lLoadsHarmo = ASTER_FALSE
        model = modelZ
        nbLoadUser = 0

        if (postCompResu%resultType .eq. "DYNA_HARMO") then
! --------- Get loads from user (for dynamic)
            call lislec(loadKeyword, 'MECANIQUE', jvBase, listLoadUser(1:19))
            call lisnnb(listLoadUser(1:19), nbLoadUser)
            lLoadsHarmo = ASTER_TRUE
        else
! --------- Get loads/BC from user
            listLoadPrep%model = model
            listLoadPrep%lHasPilo = ASTER_FALSE
            listLoadPrep%funcIsCplx = ASTER_FALSE
            if (postCompResu%lTher) then
                call ntdoch(listLoadPrep, listLoadUser, jvBase)
            else
                call nmdoch(listLoadPrep, listLoadUser, jvBase)
            end if
            call getNbLoadsFromList(listLoadUser, nbLoadUser)
        end if
        lLoadsFromUser = nbLoadUser .ne. 0
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getListLoadsResu
!
! Get list of loads from result datastructure
!
! In  postCompResu      : calculation of a post-processing option - Parameters of results
! In  numeStore         : current storing index
! Out lLoadsFromResu    : flag when loads are defined in result datastructure
! Out listLoadResu      : name of datastructure for list of loads from results datastructure
!
! --------------------------------------------------------------------------------------------------
    subroutine getListLoadsResu(postCompResu, numeStore, lLoadsFromResu, listLoadResu)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(POST_COMP_RESU), intent(in) :: postCompResu
        integer(kind=8), intent(in) :: numeStore
        aster_logical, intent(out) :: lLoadsFromResu
        character(len=24), intent(out) :: listLoadResu
! ----- Local
        integer(kind=8) :: nbLoadResu
        integer(kind=8) :: jvPara
!   ------------------------------------------------------------------------------------------------
!
        listLoadResu = "&&CCFNRN.LOADRESU"
        lLoadsFromResu = ASTER_FALSE
        nbLoadResu = 0

        if (postCompResu%resultType .ne. "DYNA_HARMO") then
! --------- Get loads from result datastructure
            call rsadpa(postCompResu%resultIn, 'L', 1, 'EXCIT', numeStore, &
                        0, sjv=jvPara, istop=0)
            listLoadResu = zk24(jvPara)
        end if

! ----- Get number of loads
        call getNbLoadsFromList(listLoadResu, nbLoadResu)
        lLoadsFromResu = nbLoadResu .ne. 0
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setListLoads
!
! Set list of loads
!
! In  lLoadsFromUser    : flag when loads are defined by user
! In  listLoadUser      : name of datastructure for list of loads from user
! In  lLoadsFromResu    : flag when loads are defined in result datastructure
! In  listLoadResu      : name of datastructure for list of loads from results datastructure
! In  lLoadsHarmo       : flag for loads in harmonic case
! IO  postCompPara      : calculation of a post-processing option - Parameters
!
! --------------------------------------------------------------------------------------------------
    subroutine setListLoads(lLoadsFromUser, lLoadsFromResu, listLoadUser, listLoadResu, &
                            lLoadsHarmo, postCompPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: lLoadsFromUser, lLoadsFromResu, lLoadsHarmo
        character(len=24), intent(in) :: listLoadUser, listLoadResu
        type(POST_COMP_PARA), intent(inout) :: postCompPara
! ----- Local
        character(len=24) :: listLoad
        aster_logical :: lConsistent, lLoadsAreDefined
        aster_logical :: hasPiloLoads, lPilo1, lPilo2, fromResu, fromUser
!   ------------------------------------------------------------------------------------------------
!
        hasPiloLoads = ASTER_FALSE
        lLoadsAreDefined = ASTER_FALSE
        listLoad = " "

! ----- Check consistency
        fromUser = lLoadsFromUser
        fromResu = lLoadsFromResu
        if (fromUser .and. fromResu) then
            fromUser = ASTER_TRUE
            fromResu = ASTER_FALSE
            call checkConsistency(listLoadUser, listLoadResu, lConsistent)
            if (.not. lConsistent) then
                call utmess('A', 'CHARGES9_53')
            end if
        end if

! ----- Name of list of loads
        if (fromUser) then
            listLoad = listLoadUser
        elseif (fromResu) then
            listLoad = listLoadResu
        else
            listLoad = " "
        end if
        lLoadsAreDefined = fromResu .or. fromUser

! ----- Get "PILOTAGE" loads ?
        hasPiloLoads = ASTER_FALSE
        if (lLoadsAreDefined) then
            if (.not. lLoadsHarmo) then
                lPilo1 = ischar(listLoad, 'DIRI', 'PILO')
                lPilo2 = ischar(listLoad, 'NEUM', 'PILO')
                hasPiloLoads = lPilo1 .or. lPilo2
            end if
        end if

! ----- Set values
        postCompPara%lLoadsAreDefined = lLoadsAreDefined
        postCompPara%listLoad = listLoad
        postCompPara%hasPiloLoads = hasPiloLoads
        postCompPara%lLoadsFromUser = fromUser
        postCompPara%lLoadsHarmo = lLoadsHarmo
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLigrel
!
! Get ligrel to compute field
!
! In  option            : option to compute
! In  model             : model
! In  result            : result datastructure
! In  nbStore           : number of storing indexes
! In  jvBase            : jeveux base to create ligrel
! IO  postCompRest      : object for management of restriction
! Out ligrel            : finite elements descriptor
!
! --------------------------------------------------------------------------------------------------
    subroutine getLigrel(option, model, result, jvBaseZ, postCompRest, &
                         ligrel)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: option
        character(len=8), intent(in) :: model, result
        character(len=*), intent(in) :: jvBaseZ
        type(POST_COMP_REST), intent(inout) :: postCompRest
        character(len=24), intent(out) :: ligrel
! ----- Locals
        character(len=16), parameter :: factorKeyword = " "
        integer(kind=8), parameter :: iOccZero = 0
        integer(kind=8) :: nbLigr, iLigr
        character(len=8) :: jvBaseLigrel, jvBase
        character(len=24) :: noojb, modelLigrel
        character(len=8) :: mesh
        character(len=24), parameter :: listCellJv = "&&MEDOM2.LISTE_MAILLES"
        integer(kind=8) :: nbCell
        integer(kind=8), pointer :: listCell(:) => null()
        integer(kind=8) :: nbCellSiro
        integer(kind=8), pointer :: listCellSiro(:) => null()
        aster_logical :: partialMesh, createFED, modelExist
!   ------------------------------------------------------------------------------------------------
!
        nbLigr = postCompRest%nbLigr
        createFED = ASTER_FALSE
        jvBase = jvBaseZ

! ----- Model already used ?
        modelExist = ASTER_FALSE
        jvBaseLigrel = ' '
        do iLigr = 1, nbLigr
            if (model .eq. postCompRest%listModel(iLigr)) then
                modelExist = ASTER_TRUE
                jvBaseLigrel = postCompRest%listJvBase(iLigr)
            end if
        end do

! ----- Get ligrel
        if (modelExist .and. (jvBaseLigrel .eq. jvBase)) then
            ASSERT(nbLigr .gt. 0)
! --------- Already used => get ligrel
            ligrel = postCompRest%listFED(nbLigr)
        else
! --------- We need a model
            if (model .eq. " ") then
                call utmess('F', 'CALCCHAMP1_10')
            end if
            call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
            call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! --------- Original definition
            partialMesh = hasCellsDefinedFromCmd(factorKeyword, iOccZero)
            createFED = ASTER_FALSE

! --------- Generate name of ligrel
            if (partialMesh .or. option .eq. "SIRO_ELEM") then
                noojb = '12345678.LIGR000000.LIEL'
                call gnomsd(result, noojb, 14, 19)
                ligrel = noojb(1:19)
            else
                ligrel = modelLigrel
            end if

! --------- Get list of cells
            if (partialMesh) then
! ------------- Get list of cells from user
                call getelem(mesh, factorKeyword, iOccZero, 'F', listCellJv, nbCell, model=model)
                ASSERT(nbCell .ge. 1)
                call jeveuo(listCellJv, 'L', vi=listCell)
                createFED = ASTER_TRUE
            else
                createFED = ASTER_FALSE
            end if

! --------- Create FED (ligrel) or get it from model
            if (option .eq. "SIRO_ELEM") then
                if (.not. partialMesh) then
                    call utmamo(model, nbCell, listCellJv)
                    ASSERT(nbCell .ge. 1)
                    call jeveuo(listCellJv, 'L', vi=listCell)
                end if
                call getCellSiroElem(model, nbCell, listCell, &
                                     nbCellSiro, listCellSiro)
                call createFEDFromList(model, jvBase, ligrel, &
                                       nbCellSiro, listCellSiro)
                AS_DEALLOCATE(vi=listCellSiro)
                createFED = ASTER_TRUE
            else
                if (createFED) then
                    call createFEDFromList(model, jvBase, ligrel, &
                                           nbCell, listCell)
                end if
            end if
            call jedetr(listCellJv)

! --------- Add new ligrel in list
            if (createFED) then
                nbLigr = nbLigr+1
                postCompRest%nbLigr = nbLigr
                ASSERT(nbLigr .le. postCompRest%nbLigrMaxi)
                postCompRest%listFED(nbLigr) = ligrel
                postCompRest%listModel(nbLigr) = model
                postCompRest%listJvBase(nbLigr) = jvBase(1:1)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! checkOption
!
! Check option (error, alarm, etc.)
!
! In  option            : option to compute
! In  ligrel            : finite element descriptor (ligrel)
! In  resultType        : type of results datastructure
!
! --------------------------------------------------------------------------------------------------
    subroutine checkOption(optionZ, ligrelZ, resultType)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: optionZ, ligrelZ
        character(len=16), intent(in) :: resultType
! ----- Locals
        aster_logical :: lPipe, lShell
        character(len=8) :: answer
!   ------------------------------------------------------------------------------------------------
!
        if (optionZ(1:4) .eq. "EPSP") then
            call dismoi('EXI_TUYAU', ligrelZ, 'LIGREL', repk=answer)
            lPipe = answer .eq. "OUI"
            call dismoi('EXI_COQUE', ligrelZ, 'LIGREL', repk=answer)
            lShell = answer .eq. "OUI"
            if (lPipe .or. lShell) then
                call utmess('I', 'ELEMENTS3_13')
            end if
        end if
        if (optionZ .eq. 'SIEQ_ELGA') then
            if (resultType .eq. 'FOURIER_ELAS') then
                call utmess('F', 'CALCULEL6_83', sk=optionZ(1:16))
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getCellSiroElem
!
! Get cells only for SIRO_ELEM option
!
! In  model             : model
! In  nbCell            : number of cells
! Ptr listCell          : list of cells
! Out nbCellSiro        : number of cells for SIRO_ELEM (skin + 3D support)
! Ptr listCellSiro      : list of cells for SIRO_ELEM (skin + 3D support)
!
! --------------------------------------------------------------------------------------------------
    subroutine getCellSiroElem(modelZ, nbCell, listCell, &
                               nbCellSiro, listCellSiro)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: modelZ
        integer(kind=8), intent(in) :: nbCell
        integer(kind=8), pointer :: listCell(:)
        integer(kind=8), intent(out) :: nbCellSiro
        integer(kind=8), pointer :: listCellSiro(:)
! ----- Locals
        character(len=16), parameter :: option = "SIGM_ELNO"
        character(len=8), parameter :: paramIn = "PSIEFNOR"
        character(len=8) :: model, mesh
        integer(kind=8) :: cellDime, cellSkinDime
        integer(kind=8) :: iCellSkin, nbCellSkin
        character(len=24) :: modelLigrel
        integer(kind=8), pointer :: listCellSigm(:) => null()
        integer(kind=8), pointer :: listCellSkin(:) => null()
        integer(kind=8) :: nbCellSigm
        integer(kind=8), pointer :: listCellSupp(:) => null()
        integer(kind=8) :: nbCellSupp, cellNumeSupp, cellNumeSkin
        aster_logical :: lCell1d, lCell2d
!   ------------------------------------------------------------------------------------------------
!
        model = modelZ

! ----- Access to model FED
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! ----- Access to mesh
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! ----- Prepare dimensions
        call dismoi('DIM_GEOM', mesh, 'MAILLAGE', repi=cellDime)
        cellSkinDime = cellDime-1
        lCell1d = ASTER_FALSE
        lCell2d = ASTER_FALSE
        if (cellSkinDime .eq. 1) then
            lCell1d = ASTER_TRUE
        else if (cellSkinDime .eq. 2) then
            lCell2d = ASTER_TRUE
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Get only skin cells
        call getSkinCell(mesh, cellSkinDime, &
                         nbCell, listCell, &
                         nbCellSkin, listCellSkin)
        if (nbCellSkin .eq. 0) then
            call utmess('F', 'CALCULEL5_54')
        end if

! ----- Create list of cells with SIGM_ELNO option
        call getAllCellsWithOption(modelLigrel, option, paramIn, &
                                   nbCellSigm, listCellSigm)

! ----- Get 3D support cells
        nbCellSupp = nbCellSkin
        AS_ALLOCATE(vi=listCellSupp, size=nbCellSupp)
        call getSkinCellSupport(mesh, nbCellSkin, listCellSkin, lCell2d, lCell1d, &
                                listCellSupp, &
                                nbCellSigm, listCellSigm)

! ----- Create final list of cells (skin + 3D support)
        nbCellSiro = 2*nbCellSkin
        AS_ALLOCATE(vi=listCellSiro, size=nbCellSiro)
        do iCellSkin = 1, nbCellSkin
            cellNumeSkin = listCellSkin(iCellSkin)
            cellNumeSupp = listCellSupp(iCellSkin)
            if (cellNumeSupp .gt. 0) then
                if (cellDime .eq. 3) then
                    call oriem1(mesh, '3D', cellNumeSkin, cellNumeSupp)
                end if
                listCellSiro(iCellSkin) = cellNumeSkin
                listCellSiro(nbCellSkin+iCellSkin) = cellNumeSupp
                ! listCellSupp(iCellSkin) = cellNumeSupp
            end if
        end do

! ----- Cleaning
        AS_DEALLOCATE(vi=listCellSkin)
        AS_DEALLOCATE(vi=listCellSigm)
        AS_DEALLOCATE(vi=listCellSupp)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initCompRest
!
! Initialisation of datastructure for managment of restriction of computation (partial mesh)
!
! In  nbStore           : number of storing indexes
! IO  postCompRest      : object for management of restriction
!
! --------------------------------------------------------------------------------------------------
    subroutine initCompRest(nbStore, postCompRest)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nbStore
        type(POST_COMP_REST), intent(inout) :: postCompRest
!   ------------------------------------------------------------------------------------------------
!
        postCompRest%nbLigr = 0
        postCompRest%nbLigrMaxi = 2*nbStore
        AS_ALLOCATE(vk24=postCompRest%listFED, size=2*nbStore)
        AS_ALLOCATE(vk8=postCompRest%listModel, size=2*nbStore)
        AS_ALLOCATE(vk8=postCompRest%listJvBase, size=2*nbStore)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! deleteCompRest
!
! Delete of datastructure for managment of restriction of computation (partial mesh)
!
! IO  postCompRest      : object for management of restriction
!
! --------------------------------------------------------------------------------------------------
    subroutine deleteCompRest(postCompRest)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(POST_COMP_REST), intent(inout) :: postCompRest
!   ------------------------------------------------------------------------------------------------
!
        AS_DEALLOCATE(vk24=postCompRest%listFED)
        AS_DEALLOCATE(vk8=postCompRest%listModel)
        AS_DEALLOCATE(vk8=postCompRest%listJvBase)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module postComp_module
