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
subroutine vrcins(modelz, chmatz, carelz, inst, chvarc, &
                  codret, nompaz, basez)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/imprsd.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcin1.h"
#include "asterfort/vrcin2.h"
!
    character(len=2) :: codret
    character(len=19) :: chvarc
    character(len=*) :: chmatz, carelz, modelz
    character(len=*), intent(in), optional :: nompaz
    character(len=1), intent(in), optional :: basez
    real(kind=8) :: inst
! ======================================================================
!   BUT : FABRIQUER LE CHAMP DE VARIABLES DE COMMANDE CORRESPONDANT A
!         UN INSTANT DONNE.
!   ARGUMENTS :
!   MODELZ (K8)  IN/JXIN : SD MODELE
!   CHMATZ (K8)  IN/JXIN : SD CHAM_MATER
!   CARELZ (K8)  IN/JXIN : SD CARA_ELEM (SOUS-POINTS)
!   INST   (R)   IN      : VALEUR DE L'INSTANT
!   CHVARC (K19) IN/JXOUT: SD CHAM_ELEM/EL** CONTENANT LES VARC
!   CODRET (K2)  OUT : POUR CHAQUE RESULTAT, 'OK' SI ON A TROUVE,
!                                            'NO' SINON
!   nompaz (k8)  in (optional) : nom du parametre parmi (PVARCPR, PVARCNO)
!                                servant a  allouer le cham_elem "chvarc"
!                        PVARCPR => chvarc = cham_elem/ELGA (par defaut)
!                        PVARCNO => chvarc = cham_elem/ELNO
!
! ----------------------------------------------------------------------
!
!
    integer(kind=8) :: iret, ichs, nbchs, jcesd, jcesl
    integer(kind=8) :: nbcmp, kcmp, kcvrc
    integer(kind=8) :: nbma, ima, nbpt, nbsp, ipt, isp, iad, iad1
    integer(kind=8) :: jce1d, jce1l, nncp, n1, k
    real(kind=8) :: valeur, rundef
    character(len=19) :: chvars, modelFED, chs
    character(len=16) :: nomte
    aster_logical :: avrc, dbg
    integer(kind=8) :: ibid, nbcvrc, nute, jmaille, vali(2)
    character(len=8) :: model, chmat, carele, exteVariName1, varc2, nocmp1, nocmp2, nomail
    character(len=8) :: nompar
    character(len=8), pointer :: cvrcvarc(:) => null()
    character(len=24), pointer :: liste_ch(:) => null()
    character(len=8), pointer :: cesc(:) => null()
    character(len=16), pointer :: liste_sd(:) => null()
    character(len=8), pointer :: cvrccmp(:) => null()
    character(len=1) :: base
    integer(kind=8), pointer :: cesvi(:) => null()
    real(kind=8), pointer :: ce1v(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    aster_logical :: exival
! ----------------------------------------------------------------------
!
    call jemarq()
!
    if (present(basez)) then
        base = basez
    else
        base = 'V'
    end if
!
!   nom du parametre "nompar" servant a allouer le cham_elem "chvarc" :
!   PVARCPR <-> ELGA (par defaut)
!   PVARCNO <-> ELNO
    nompar = 'PVARCPR'
    if (present(nompaz)) then
        ASSERT((nompaz .eq. 'PVARCPR') .or. (nompaz .eq. 'PVARCNO'))
        nompar = nompaz
    end if
!
    chmat = chmatz
    carele = carelz
    model = modelz

    call detrsd('CHAM_ELEM', chvarc)
!
    call jeexin(chmat//'.CVRCVARC', iret)
!     AVRC : .TRUE. SI AFFE_MATERIAU/AFFE_VARC EST UTILISE
    avrc = (iret .gt. 0)
    if (.not. avrc) goto 999
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelFED)
    call jeveuo(modelFED//'.TYFE', 'L', jmaille)
!
!
!   1. interpolation en temps :
!      fabrication d'une liste de cham_elem_s / elga
!      contenant les vrc a l'instant inst
!      calcul de  chmat.liste_ch(:) et chmat.liste_sd(:)
!   -----------------------------------------------------
    call vrcin1(model, chmat, carele, inst, codret, &
                nompar)
!
    call jeexin(chmat//'.LISTE_SD', iret)
    ASSERT(iret .ne. 0)

!
!   2. allocation du champ_elem_s resultat (chvars)
!      calcul de chmat.cesvi
!      (cette etape est economisee d'un instant a l'autre)
!   -------------------------------------------------------------
    chvars = chmat//'.CHVARS'
    call jeexin(chmat//'.CESVI', iret)
    if (iret .eq. 0) call vrcin2(model, chmat, carele, chvars, nompar)
!
!
!   3. concatenation des champs de .liste_ch  dans chvars :
!   -----------------------------------------------------
    call jeveuo(chmat//'.LISTE_CH', 'L', vk24=liste_ch)
    call jelira(chmat//'.LISTE_CH', 'LONMAX', nbchs)
    call jeveuo(chmat//'.LISTE_SD', 'L', vk16=liste_sd)
    call jeveuo(chmat//'.CVRCVARC', 'L', vk8=cvrcvarc)
    call jeveuo(chmat//'.CVRCCMP', 'L', vk8=cvrccmp)
    call jelira(chmat//'.CVRCCMP', 'LONMAX', nbcvrc)
!
    call jeveuo(chvars//'.CESD', 'L', jce1d)
    call jeveuo(chvars//'.CESL', 'E', jce1l)
    call jeveuo(chmat//'.CESVI', 'L', vi=cesvi)
!
!   -- il faut remettre cesv a nan:
    rundef = r8nnem()
    call jeveuo(chvars//'.CESV', 'E', vr=ce1v)
    call jelira(chvars//'.CESV', 'LONMAX', n1)
    do k = 1, n1
        ce1v(k) = rundef
        zl(jce1l-1+k) = .true.
    end do
!

!
    do ichs = 1, nbchs
        chs = liste_ch(ichs) (1:19)
        exteVariName1 = liste_sd(7*(ichs-1)+4) (1:8)
        call jeveuo(chs//'.CESD', 'L', jcesd)
        call jeveuo(chs//'.CESL', 'L', jcesl)
        call jeveuo(chs//'.CESV', 'L', vr=cesv)
        call jeveuo(chs//'.CESC', 'L', vk8=cesc)
        call jelira(chs//'.CESC', 'LONMAX', nbcmp)
!
        do kcmp = 1, nbcmp
            nocmp1 = cesc(kcmp)
!
!         -- calcul de kcvrc :
            do kcvrc = 1, nbcvrc
                varc2 = cvrcvarc(kcvrc)
                nocmp2 = cvrccmp(kcvrc)
                if ((exteVariName1 .eq. varc2) .and. (nocmp1 .eq. nocmp2)) goto 4
            end do
            goto 2
!
4           continue
            ASSERT(kcvrc .ge. 1 .and. kcvrc .le. nbcvrc)
!
!           -- boucle sur les mailles :
            nbma = zi(jcesd-1+1)
            ASSERT(nbma .eq. zi(jce1d-1+1))
!
            do ima = 1, nbma
                nbpt = zi(jcesd-1+5+4*(ima-1)+1)
                nbsp = zi(jcesd-1+5+4*(ima-1)+2)
                if (nbsp*nbsp .eq. 0) goto 70
!
                call cesexi('C', jce1d, jce1l, ima, 1, &
                            1, kcvrc, iad1)
!               -- la maille n'est pas concernee par les variables de commande :
                if (iad1 .eq. 0) goto 70
!
!
!               -- la maille porte un element fini qui saurait utiliser
!                  les variables de commande mais elle n'est pas affectee.
!                  on espere que les routines te00ij arreteront en <f> si necessaire.
                if (iad1 .lt. 0) goto 70
!
!               -- On regarde si le champ possede des valeurs sur la maille :
                exival = .false.
                do ipt = 1, nbpt
                    do isp = 1, nbsp
                        call cesexi('C', jcesd, jcesl, ima, ipt, &
                                    isp, kcmp, iad)
                        if (iad .gt. 0) exival = .true.
                    end do
                end do
                if (.not. exival) goto 70
!
!               -- controle du nombre de points :
                ASSERT(nbpt .eq. zi(jce1d-1+5+4*(ima-1)+1))
!
!               -- Controle du nombre de sous-points :
                if (nbsp .ne. zi(jce1d-1+5+4*(ima-1)+2)) then
!                   -- issue23456 : il peut arriver que nbsp=1 mais sans aucune valeur :
                    if (nbsp .eq. 1 .and. .not. exival) goto 70
!

                    nomail = int_to_char8(ima)
                    nute = zi(jmaille-1+ima)
                    call jenuno(jexnum('&CATA.TE.NOMTE', nute), nomte)
                    vali(1) = zi(jce1d-1+5+4*(ima-1)+2)
                    vali(2) = nbsp
                    call utmess('F', 'VARC1_9', sk=nocmp1, ni=2, vali=vali)
                end if
!
!
                do ipt = 1, nbpt
                    do isp = 1, nbsp
                        call cesexi('C', jcesd, jcesl, ima, ipt, &
                                    isp, kcmp, iad)
                        if (iad .gt. 0) then
                            call cesexi('C', jce1d, jce1l, ima, ipt, &
                                        isp, kcvrc, iad1)
                            ASSERT(iad1 .gt. 0)
                            if (cesvi(iad1) .eq. ichs) then
                                valeur = cesv(iad)
                                zl(jce1l-1+iad1) = .true.
                                ce1v(iad1) = valeur
                            end if
                        end if
                    end do
                end do
70              continue
            end do
!
2           continue
        end do
    end do
!
!
!   4. recopie du champ simple dans le champ chvarc
!   -----------------------------------------------------
    call cescel(chvars, modelFED, 'INIT_VARC', nompar, 'NAN', &
                nncp, base, chvarc, 'F', ibid)
!
    dbg = .false.
    if (dbg) call imprsd('CHAMP', chvarc, 6, 'VRCINS/CHVARC')
!
999 continue
    call jedema()
end subroutine
