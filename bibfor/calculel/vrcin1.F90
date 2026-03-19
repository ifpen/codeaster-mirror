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
subroutine vrcin1(model, mateField, caraElem, timeCurr, codret, nompar)
!
    use ExternalStateVariablePrep_module
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/cesces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesvar.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnsces.h"
#include "asterfort/codent.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/fointe.h"
#include "asterfort/indk80.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/manopg.h"
#include "asterfort/rsinch.h"
#include "asterfort/vrcin_elno.h"
#include "asterfort/wkvect.h"
#include "asterfort/xvrcin.h"
!
    character(len=8), intent(in) :: model, mateField, caraElem
    real(kind=8), intent(in) :: timeCurr
    character(len=2), intent(out) :: codret
    character(len=*), intent(in) :: nompar

! ======================================================================
!   BUT : FAIRE L'INTERPOLATION AU TEMPS INST DES DIFFERENTS CHAMPS
!         DE VARIABLES DE COMMANDE.
!         CES CHAMPS SONT PASSES AUX POINTS DE GAUSS MECANIQUE.
!         (INIT_VARC/PVARCPR)
!
!   IN :
!     MODELE (K8)  IN/JXIN : SD MODELE
!     CHMAT  (K8)  IN/JXIN : SD CHAM_MATER
!     CARELE  (K8)  IN/JXIN : SD CARA_ELEM
!     INST   (R)   IN      : VALEUR DE L'INSTANT
!     nompar (k8)  in      : nom du parametre parmi (PVARCPR, PVARCNO)
!                            servant a  allouer le cham_elem "chvarc"
!                            PVARCPR => LISTE_CH(i) = cham_elem_s/ELGA
!                            PVARCNO => LISTE_CH(i) = cham_elem_s/ELNO
!
!   OUT :
!       - CREATION DE CHMAT//'.LISTE_CH' V V K24  LONG=NBCHS
!          .LISTE_CH(I) : IEME CHAM_ELEM_S / EL** PARTICIPANT A
!                         LA CREATION DU CHAMP DE CVRC
!       - CREATION DE CHMAT//'.LISTE_SD' V V K16  LONG=7*NBCHS
!          .LISTE_SD(7*(I-1)+1) : /'EVOL' /'CHAMP' :
!               TYPE DE LA SD DONT EST ISSU .LISTE_CHS(I)
!          .LISTE_SD(7*(I-1)+2) : NOMSD
!               NOM DE L'EVOL (OU DU CHAMP) DONT EST ISSU .LISTE_CHS(I)
!          .LISTE_SD(7*(I-1)+3) : NOMSYM / ' '
!               SI 'EVOL' : NOM SYMBOLIQUE DU CHAMP DONT EST
!               ISSU .LISTE_CHS(I). SINON : ' '
!          .LISTE_SD(7*(I-1)+4) : VARC
!               VARC ASSOCIE A .LISTE_CHS(I).
!          .LISTE_SD(7*(I-1)+5) : PROLGA
!               (SI EVOL : TYPE DE PROLONGEMENT A GAUCHE)
!          .LISTE_SD(7*(I-1)+6) : PROLDR
!               (SI EVOL : TYPE DE PROLONGEMENT A DROITE)
!          .LISTE_SD(7*(I-1)+7) : FINST (OU ' ')
!               (SI EVOL : FONCTION DE TRANSFORMATION DU TEMPS)
!        CODRET (K2) : POUR CHAQUE RESULTAT, 'OK' SI ON A TROUVE,
!                                            'NO' SINON
! ----------------------------------------------------------------------

    real(kind=8) :: timeInter
    integer(kind=8) :: n1, ibid, nbCellMesh, jcesd1, jcesl1, iad, lonk80
    integer(kind=8) :: itrou, nbk80, iCellMesh, jlk80, iret, nbchs, jlissd, ichs
    integer(kind=8) :: jlisch, nval1
    aster_logical :: l_xfem, l_elga
    character(len=8) :: exteVariName, mesh, affeType
    character(len=8) :: funcExtrRight, funcExtrLeft, resultUser, funcResult

    character(len=8) :: physQuan, fieldUserDisc, dsUser
    character(len=16) :: fieldType
    character(len=19) :: exteVariMap2, chs, cesmod, celmod, modelFED, mnoga, dceli
    character(len=19) :: cns1, fieldUser, celtmp

    character(len=19), parameter :: ces1 = '&&VRCIN1.CES1'

    character(len=80) :: stringStore, stringStorePrev
    character(len=16), pointer :: cesv(:) => null()
    real(kind=8), parameter :: prec = 1.0d-10
    character(len=8), parameter :: crit = 'ABSOLU'
    save nval1

    character(len=24) :: jvNameCVRCVarc, jvNameCVRCGd
    character(len=8), pointer :: CVRCVarc(:) => null()
    character(len=8), pointer :: CVRCGd(:) => null()
    integer(kind=8) :: nbCmpTotal, iCmpTotal
! ----------------------------------------------------------------------

    call jemarq()
    call infmaj()
!
!   nom du parametre "nompar" servant a allouer le cham_elem "celmod" :
!   PVARCPR <-> ELGA (par defaut dans vrcins)
!   PVARCNO <-> ELNO
    ASSERT((nompar .eq. 'PVARCPR') .or. (nompar .eq. 'PVARCNO'))
    l_elga = .true.
    if (nompar .eq. 'PVARCNO') then
        l_elga = .false.
    end if
!
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCellMesh)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelFED)

    jvNameCVRCVarc = mateField//'.CVRCVARC'
    jvNameCVRCGd = mateField//'.CVRCGD'
    call jelira(jvNameCVRCVarc, 'LONMAX', nbCmpTotal)
    call jeveuo(jvNameCVRCVarc, 'L', vk8=CVRCVarc)
    call jelira(jvNameCVRCVarc, 'LONMAX', nbCmpTotal)

    call jeveuo(jvNameCVRCGd, 'L', vk8=CVRCGd)

    codret = 'OK'

!     1. CREATION DE CHMAT.LISTE_SD :
!     -------------------------------
    call jeexin(mateField//'.LISTE_SD', iret)
    if (iret .eq. 0) then
        exteVariName = ' '
        stringStorePrev = ' '
        nbk80 = 0
        lonk80 = 5
        call wkvect('&&VRCIN1.LK80', 'V V K80', lonk80, jlk80)
        do iCmpTotal = 1, nbCmpTotal
            if (CVRCVarc(iCmpTotal) .eq. exteVariName) cycle

! --------- Access to map of current external state variable
            exteVariName = CVRCVarc(iCmpTotal)
            exteVariMap2 = mateField//'.'//exteVariName//'.2'
            call carces(exteVariMap2, 'ELEM', ' ', 'V', ces1, 'A', iret)
            ASSERT(iret .eq. 0)
            call jeveuo(ces1//'.CESD', 'L', jcesd1)
            call jeveuo(ces1//'.CESL', 'L', jcesl1)
            call jeveuo(ces1//'.CESV', 'L', vk16=cesv)

            do iCellMesh = 1, nbCellMesh
                call cesexi('C', jcesd1, jcesl1, iCellMesh, 1, 1, 1, iad)
                if (iad .le. 0) cycle
                iad = iad-1

! ------------- Get parameter of external state variable on current cell
                affeType = cesv(iad+2) (1:8)
                dsUser = cesv(iad+3) (1:8)
                fieldType = cesv(iad+4)
                funcExtrLeft = cesv(iad+5) (1:8)
                funcExtrRight = cesv(iad+6) (1:8)
                funcResult = cesv(iad+7) (1:8)
                ASSERT((affeType .eq. 'EVOL') .or. (affeType .eq. 'CHAMP'))

! ------------- Create composite string to store parameters of external state variable
                stringStore = ' '
                stringStore(1:8) = affeType
                stringStore(9:16) = dsUser
                stringStore(17:32) = fieldType
                stringStore(33:40) = exteVariName
                stringStore(41:48) = funcExtrLeft
                stringStore(49:56) = funcExtrRight
                stringStore(57:64) = funcResult

                if (stringStore .ne. stringStorePrev) then
                    stringStorePrev = stringStore
                    itrou = indk80(zk80(jlk80), stringStore, 1, nbk80)
                    if (itrou .le. 0) then
                        nbk80 = nbk80+1
                        if (nbk80 .gt. lonk80) then
                            lonk80 = 2*lonk80
                            call juveca('&&VRCIN1.LK80', lonk80)
                            call jeveuo('&&VRCIN1.LK80', 'E', jlk80)
                        end if
                        zk80(jlk80-1+nbk80) = stringStore
                    end if
                end if
            end do
            call detrsd('CHAM_ELEM_S', ces1)
        end do

        nbchs = nbk80
        if (nbchs .eq. 0) then
            call jedetr('&&VRCIN1.LK80')
            goto 999
        end if

! ----- nbchs: number of different affectations of external state variable
        call wkvect(mateField//'.LISTE_SD', 'V V K16', 7*nbchs, jlissd)
        do ichs = 1, nbchs
            stringStore = zk80(jlk80-1+ichs)
            zk16(jlissd-1+7*(ichs-1)+1) = stringStore(1:8)
            zk16(jlissd-1+7*(ichs-1)+2) = stringStore(9:16)
            zk16(jlissd-1+7*(ichs-1)+3) = stringStore(17:32)
            zk16(jlissd-1+7*(ichs-1)+4) = stringStore(33:40)
            zk16(jlissd-1+7*(ichs-1)+5) = stringStore(41:48)
            zk16(jlissd-1+7*(ichs-1)+6) = stringStore(49:56)
            zk16(jlissd-1+7*(ichs-1)+7) = stringStore(57:64)
        end do
        call jedetr('&&VRCIN1.LK80')
    end if

!   2. CREATION DE CHMAT.LISTE_CH :
!   -------------------------------
    call jeveuo(mateField//'.LISTE_SD', 'L', jlissd)
    call jelira(mateField//'.LISTE_SD', 'LONMAX', n1)
    nbchs = n1/7
    ASSERT(n1 .eq. 7*nbchs)
    call jedetr(mateField//'.LISTE_CH')
    call wkvect(mateField//'.LISTE_CH', 'V V K24', nbchs, jlisch)
    chs = mateField//'.CHS000'

!   2.0.1  CREATION DE CESMOD :
!   ---------------------------
    cesmod = model//'.VRC.CESMOD'
!   --  cesmod n'est pas detruit pour gagner du temps
    call exisd('CHAM_ELEM_S', cesmod, iret)
    if (iret .eq. 0) then
        celmod = '&&VRCIN1.CELMOD'
        dceli = '&&VRCIN1.DCELI'
        call cesvar(caraElem, ' ', modelFED, dceli)
        call alchml(modelFED, 'INIT_VARC', nompar, 'V', celmod, &
                    iret, dceli)
        ASSERT(iret .eq. 0)
        call detrsd('CHAMP', dceli)
        call celces(celmod, 'V', cesmod)
        call jelira(celmod//'.CELV', 'LONMAX', nval1)
        call detrsd('CHAMP', celmod)
    end if

!   2.0.2  CREATION DE MNOGA :
!   ---------------------------
    mnoga = model//'.VRC.MNOGA'
    call exisd('CHAM_ELEM_S', mnoga, iret)
    if (iret .eq. 0) call manopg(model, modelFED, 'INIT_VARC', 'PVARCPR', mnoga)

    do ichs = 1, nbchs

! ----- Generate name of field for current state variable affectation
        call codent(ichs, 'D0', chs(13:15))
        zk24(jlisch-1+ichs) = chs

! ----- Get current state variable affectation
        affeType = zk16(jlissd-1+7*(ichs-1)+1) (1:8)
        exteVariName = zk16(jlissd-1+7*(ichs-1)+4) (1:8)
        itrou = indik8(CVRCVarc, exteVariName, 1, nbCmpTotal)
        ASSERT(itrou .gt. 0)
        physQuan = CVRCGd(itrou)

!       2.1 INTERPOLATION EN TEMPS => NOMCH
!       ------------------------------------
        resultUser = '        '
        if (affeType .eq. 'EVOL') then
!           -- SI TYSD='EVOL', ON INTERPOLE AU TEMPS INST
            resultUser = zk16(jlissd-1+7*(ichs-1)+2) (1:8)
            fieldType = zk16(jlissd-1+7*(ichs-1)+3)
            funcExtrLeft = zk16(jlissd-1+7*(ichs-1)+5) (1:8)
            funcExtrRight = zk16(jlissd-1+7*(ichs-1)+6) (1:8)
            funcResult = zk16(jlissd-1+7*(ichs-1)+7) (1:8)
            fieldUser = '&&VRCIN1.NOMCH'

! --------- Get time value for interpolate
            timeInter = timeCurr
            if (funcResult .ne. ' ') then
                call fointe('F', funcResult, 1, ['INST'], [timeCurr], &
                            timeInter, ibid)
            end if

! --------- Interpolate field for time value
            call rsinch(resultUser, fieldType, 'INST', timeInter, fieldUser, &
                        funcExtrRight, funcExtrLeft, 2, 'V', prec, crit, iret)
            ASSERT(iret .le. 12)
            if (iret .ge. 10) then
                codret = 'NO'
                goto 999
            end if

        else
            ASSERT(affeType .eq. 'CHAMP')
            fieldUser = zk16(jlissd-1+7*(ichs-1)+2)

        end if

! ----- Check conformity of field from user
        call checkField(mesh, model, physQuan, &
                        exteVariName, fieldUser, fieldUserDisc)

!       2.2.1 Cas particulier ou l'on ne souhaite pas ELGA mais ELNO
!       --------------------------------------
        if (.not. l_elga) then
            call vrcin_elno(fieldUser, cesmod, chs)
!           on passe a l'iteration suivante de la boucle sur LISTE_CH
            cycle
        end if

!       2.2.2 PASSAGE AUX POINTS DE GAUSS => CHS (cas general ELGA)
!       --------------------------------------
        if (fieldUserDisc .eq. 'CART') then
            call carces(fieldUser, 'ELGA', cesmod, 'V', chs, 'A', iret)
            ASSERT(iret .eq. 0)

        else if (fieldUserDisc .eq. 'NOEU') then
            cns1 = '&&VRCIN1.CNS1'
            call cnocns(fieldUser, 'V', cns1)
            call cnsces(cns1, 'ELGA', cesmod, mnoga, 'V', chs)
            call detrsd('CHAM_NO_S', cns1)

        else if ((fieldUserDisc .eq. 'ELNO') .or. (fieldUserDisc .eq. 'ELEM')) then
            call celces(fieldUser, 'V', ces1)
            call cesces(ces1, 'ELGA', cesmod, mnoga, ' ', 'V', chs)
            call detrsd('CHAM_ELEM_S', ces1)

        else if (fieldUserDisc .eq. 'ELGA') then
            celtmp = '&&VRCIN1.CELTMP'
            l_xfem = .false.
            call xvrcin(modelFED, fieldUser, resultUser, fieldType, celtmp, l_xfem)
            if (l_xfem) then
                call celces(celtmp, 'V', chs, l_copy_nan_=ASTER_FALSE)
                call detrsd('CHAM_ELEM', celtmp)
            else
                call celces(fieldUser, 'V', chs)
            end if

        else
            ASSERT(.false.)
        end if
        if (affeType .eq. 'EVOL') then
            call detrsd('CHAMP', fieldUser)
        end if
    end do

999 continue
    call jedema()
end subroutine
