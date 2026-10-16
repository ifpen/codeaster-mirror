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
! aslint: disable=W1501
!
subroutine dismlg(questi, nomobz, repi, repkz, ierd)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dimge1.h"
#include "asterfort/dismma.h"
#include "asterfort/dismml.h"
#include "asterfort/dismte.h"
#include "asterfort/dismtm.h"
#include "asterfort/exisd.h"
#include "asterfort/ismali.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
!
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi, repkz, nomobz
!
! --------------------------------------------------------------------------------------------------
!
!     --     DISMOI(LIGREL)
!
! --------------------------------------------------------------------------------------------------
!
!    IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE TYPE LIGREL
!    OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: dimge(3)
    aster_logical :: melang
    character(len=4) :: tytm
    character(len=8) :: calcri, mailla, nomacr, cellTypeName, k8bid, typmod2, typmod3
    character(len=16) :: elemTypeName, phenom, modelization, tyvois, formul
    character(len=19) :: nomob
    character(len=32) :: repk
    integer(kind=8) :: jlgrf, iret, nbgrel, igrel, lielSize, elemTypeNume, jsssa, n1
    integer(kind=8) :: ige2, ige1, ige3
    integer(kind=8) :: iexi, iexi2, ico
    integer(kind=8) :: jnomac, nbsm, ism, ibid
    aster_logical :: mail_quad, mail_line, lret
    integer(kind=8) :: ndime
    character(len=8), pointer :: typema(:) => null()
    integer(kind=8), pointer :: nbno(:) => null()
    integer(kind=8), pointer :: liel(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    repk = ' '
    repi = 0
    ierd = 0
!
    nomob = nomobz

    call exisd('LIGREL', nomob, iexi)
    if (iexi .ne. 1) then
        WRITE (6, *) "LIGREL <", nomob, ">"
    end if
    ASSERT(iexi .eq. 1)

    if (questi .eq. 'NOM_MAILLA') then
        call jeveuo(nomob//'.LGRF', 'L', jlgrf)
        repk = zk8(jlgrf-1+1)
    else if (questi .eq. 'JOINTS') then
        call jeveuo(nomob//'.LGRF', 'L', jlgrf)
        repk = zk8(jlgrf-1+4)
    else if (questi .eq. 'PARTITION') then
        call jeveuo(nomob//'.LGRF', 'L', jlgrf)
        repk = zk8(jlgrf-1+2)

    else if (questi .eq. 'BESOIN_VOISIN') then
        repk = 'NON'
        call jeexin(nomob//'.LIEL', iexi)
        if (iexi .gt. 0) then
            call jelira(nomob//'.LIEL', 'NUTIOC', nbgrel)
            do igrel = 1, nbgrel
                call jeveuo(jexnum(nomob//'.LIEL', igrel), 'L', vi=liel)
                call jelira(jexnum(nomob//'.LIEL', igrel), 'LONMAX', lielSize)
                elemTypeNume = liel(lielSize)
                call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
                call teattr('C', 'TYPE_VOISIN', tyvois, iret, typel=elemTypeName)
                if (iret .eq. 0) then
                    repk = 'OUI'
                    exit
                end if
            end do
        end if

    elseif (questi .eq. 'CALC_RIGI') then
        repk = 'NON'
        call jeexin(nomob//'.LIEL', iexi)
        if (iexi .gt. 0) then
            call jelira(nomob//'.LIEL', 'NUTIOC', nbgrel)
            do igrel = 1, nbgrel
                call jeveuo(jexnum(nomob//'.LIEL', igrel), 'L', vi=liel)
                call jelira(jexnum(nomob//'.LIEL', igrel), 'LONMAX', lielSize)
                elemTypeNume = liel(lielSize)
                call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
                call dismte(questi, elemTypeName, repi, calcri, ierd)
                if (calcri .eq. 'OUI') then
                    repk = 'OUI'
                    exit
                end if
            end do
        end if

    else if (questi(1:4) .eq. 'EXI_') then
        repk = 'NON'
        call jeexin(nomob//'.LIEL', iexi)
        if (questi .eq. 'EXI_ELEM') then
            if (iexi .gt. 0) then
                repk = 'OUI'
                goto 99
            end if
        else if (questi .eq. 'EXI_AMOR') then
            if (iexi .gt. 0) then
                call jelira(nomob//'.LIEL', 'NUTIOC', nbgrel)
                do igrel = 1, nbgrel
                    call jeveuo(jexnum(nomob//'.LIEL', igrel), 'L', vi=liel)
                    call jelira(jexnum(nomob//'.LIEL', igrel), 'LONMAX', lielSize)
                    elemTypeNume = liel(lielSize)
                    call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
                    if (lteatt('ABSO', 'OUI', typel=elemTypeName)) then
                        repk = 'OUI'
                    end if
                    if (repk .eq. 'OUI') exit
                end do
                if (repk .eq. 'NON') then
                    call jeexin(nomob//'.SSSA', iexi2)
                    if (iexi2 .gt. 0) then
                        call jelira(nomob//'.SSSA', 'LONMAX', n1)
                        call jeveuo(nomob//'.SSSA', 'L', jsssa)
                        call jeveuo(nomob//'.LGRF', 'L', jlgrf)
                        mailla = zk8(jlgrf-1+1)
                        call jeveuo(mailla//'.NOMACR', 'L', jnomac)
                        nbsm = n1-3
                        do ism = 1, nbsm
                            if (zi(jsssa-1+ism) .eq. 1) then
                                nomacr = zk8(jnomac-1+ism)
                                call dismml(questi, nomacr, ibid, repk, ierd)
                                if (repk .eq. 'OUI') goto 99
                            end if
                        end do
                    end if
                end if
                goto 99
            end if
        else
            if (iexi .gt. 0) then
                call jelira(nomob//'.LIEL', 'NUTIOC', nbgrel)
                do igrel = 1, nbgrel
                    call jeveuo(jexnum(nomob//'.LIEL', igrel), 'L', vi=liel)
                    call jelira(jexnum(nomob//'.LIEL', igrel), 'LONMAX', lielSize)
                    elemTypeNume = liel(lielSize)
                    call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)

                    if (questi .eq. 'EXI_RDM') then
                        if (lteatt('COQUE', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (lteatt('POUTRE', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (lteatt('DISCRET', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_POUTRE') then
                        if (lteatt('POUTRE', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_COQUE') then
                        if (lteatt('COQUE', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_GRILLE') then
                        if (lteatt('GRILLE', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_COQ3D') then
                        call dismte('MODELISATION', elemTypeName, repi, modelization, ierd)
                        if (modelization .eq. 'COQUE_3D') then
                            repk = 'OUI'
                            exit
                        end if

                    else if (questi .eq. 'EXI_PLAQUE') then
                        if (lteatt('PLAQUE', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_TUYAU') then
                        if (lteatt('TUYAU', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_POUX') then
                        if (lteatt('POUX', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_STRX') then
                        if (lteatt('STRX', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_STR2') then
                        if (elemTypeName .eq. 'MECA_POU_D_EM') repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_THM') then
                        call teattr('C', 'TYPMOD2', typmod2, iret, typel=elemTypeName)
                        if (typmod2 .eq. 'THM') then
                            repk = 'OUI'
                            call teattr('C', 'TYPMOD3', typmod3, iret, typel=elemTypeName)
                            if (typmod3 .eq. 'STEADY') then
                                repk = 'OUI_P'
                            end if
                        end if

                    else if (questi .eq. 'EXI_SECH') then
                        if (lteatt('SECH', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_NON_SECH') then
                        if (.not. lteatt('SECH', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') exit

                    else if (questi .eq. 'EXI_HHO') then
                        call teattr('C', 'TYPMOD2', typmod2, iret, typel=elemTypeName)
                        if (typmod2(1:3) .eq. 'HHO') then
                            repk = 'OUI'
                            go to 99
                        else
                            repk = 'NON'
                        end if

                    else if (questi .eq. 'EXI_NO_HHO') then
                        call teattr('C', 'TYPMOD2', typmod2, iret, typel=elemTypeName)
                        if (typmod2(1:3) .eq. 'HHO') then
                            repk = 'NON'
                        else
                            repk = 'OUI'
                            go to 99
                        end if

                    else if (questi .eq. 'EXI_HHO_LINE') then
                        if (lteatt('TYPMOD2', 'HHO', typel=elemTypeName) .or. &
                            lteatt('TYPMOD2', 'HHO_GRAD', typel=elemTypeName)) then
                            call teattr('C', 'FORMULATION', formul, iret, typel=elemTypeName)
                            if (formul .eq. 'HHO_LINE') then
                                repk = 'OUI'
                            else
                                repk = 'NON'
                            end if
                        else
                            repk = 'NON'
                        end if

                    else if (questi .eq. 'EXI_HHO_CSTE') then
                        if (lteatt('TYPMOD2', 'HHO', typel=elemTypeName) .or. &
                            lteatt('TYPMOD2', 'HHO_GRAD', typel=elemTypeName)) then
                            call teattr('C', 'FORMULATION', formul, iret, typel=elemTypeName)
                            if (formul .eq. 'HHO_CSTE') then
                                repk = 'OUI'
                            else
                                repk = 'NON'
                            end if
                        else
                            repk = 'NON'
                        end if

                    else if (questi .eq. 'EXI_HHO_QUAD') then
                        if (lteatt('TYPMOD2', 'HHO', typel=elemTypeName) .or. &
                            lteatt('TYPMOD2', 'HHO_GRAD', typel=elemTypeName)) then
                            call teattr('C', 'FORMULATION', formul, iret, typel=elemTypeName)
                            if (formul .eq. 'HHO_QUAD') then
                                repk = 'OUI'
                            else
                                repk = 'NON'
                            end if
                        else
                            repk = 'NON'
                        end if

                    else if (questi .eq. 'EXI_HHO_CUBI') then
                        if (lteatt('TYPMOD2', 'HHO', typel=elemTypeName) .or. &
                            lteatt('TYPMOD2', 'HHO_GRAD', typel=elemTypeName)) then
                            call teattr('C', 'FORMULATION', formul, iret, typel=elemTypeName)
                            if (formul .eq. 'HHO_CUBI') then
                                repk = 'OUI'
                            else
                                repk = 'NON'
                            end if
                        else
                            repk = 'NON'
                        end if

                    else if (questi .eq. 'EXI_HHO_QUAR') then
                        if (lteatt('TYPMOD2', 'HHO', typel=elemTypeName) .or. &
                            lteatt('TYPMOD2', 'HHO_GRAD', typel=elemTypeName)) then
                            call teattr('C', 'FORMULATION', formul, iret, typel=elemTypeName)
                            if (formul .eq. 'HHO_QUAR') then
                                repk = 'OUI'
                            else
                                repk = 'NON'
                            end if
                        else
                            repk = 'NON'
                        end if

                    else if (questi .eq. 'EXI_AXIS') then
                        if (lteatt('AXIS', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') then
                            exit
                        end if

                    else if (questi .eq. 'EXI_COQSOL') then
                        if ((elemTypeName .eq. 'MESSHELL_SB9' .or. &
                             elemTypeName .eq. 'MESSHELL_SB7')) then
                            repk = 'OUI'
                            exit
                        end if

                    else if (questi .eq. 'EXI_IMPE_ABSO') then
                        if (lteatt('ABSO', 'OUI', typel=elemTypeName) .and. &
                            lteatt('FLUIDE', 'OUI', typel=elemTypeName)) then
                            repk = 'OUI'
                            exit
                        end if

                    else if (questi .eq. 'EXI_CABLE') then
                        if (lteatt('CABLE', 'OUI', typel=elemTypeName)) then
                            repk = 'OUI'
                            exit
                        end if

                    else if (questi .eq. 'EXI_XFEM') then
                        if (lteatt('LXFEM', 'OUI', typel=elemTypeName)) then
                            repk = 'OUI'
                            exit
                        end if

                    else if (questi .eq. 'EXI_INCO') then
                        call dismte('MODELISATION', elemTypeName, repi, modelization, ierd)
                        if (modelization(1:7) .eq. '3D_INCO' .or. &
                            modelization(1:9) .eq. 'AXIS_INCO' .or. &
                            modelization(1:11) .eq. 'D_PLAN_INCO' .or. &
                            modelization(1:12) .eq. '3D_GRAD_INCO' .or. &
                            modelization(1:14) .eq. 'AXIS_GRAD_INCO' .or. &
                            modelization(1:16) .eq. 'D_PLAN_GRAD_INCO') then
                            repk = 'OUI'
                            exit
                        end if

                    else if (questi .eq. 'EXI_C_PLAN') then
                        if (lteatt('C_PLAN', 'OUI', typel=elemTypeName)) repk = 'OUI'
                        if (repk .eq. 'OUI') then
                            exit
                        end if

                    else if (questi .eq. 'EXI_ELTVOL') then
                        call dismte('TYPE_TYPMAIL', elemTypeName, repi, tytm, ierd)
                        if (tytm .eq. 'VOLU') then
                            repk = 'OUI'
                            exit
                        end if

                    else
                        ASSERT(ASTER_FALSE)

                    end if
                end do
            end if
        end if

    else if (questi .eq. 'AXIS') then
        repk = 'NON'
        call jeexin(nomob//'.LIEL', iexi)
        ico = 0
        if (iexi .gt. 0) then
            call jelira(nomob//'.LIEL', 'NUTIOC', nbgrel)
            do igrel = 1, nbgrel
                call jeveuo(jexnum(nomob//'.LIEL', igrel), 'L', vi=liel)
                call jelira(jexnum(nomob//'.LIEL', igrel), 'LONMAX', lielSize)
                elemTypeNume = liel(lielSize)
                call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
                if (lteatt('AXIS', 'OUI', typel=elemTypeName)) then
                    ico = ico+1
                end if
            end do
            if (ico .eq. nbgrel) repk = 'OUI'
        end if

    elseif ((questi .eq. 'NB_SM_MAILLA') .or. &
            (questi .eq. 'NB_SS_ACTI') .or. &
            (questi .eq. 'NB_NL_MAILLA')) then
        call jeexin(nomob//'.SSSA', iexi)
        if (iexi .eq. 0) then
            repi = 0
        else
            call jeveuo(nomob//'.SSSA', 'L', jsssa)
            call jelira(nomob//'.SSSA', 'LONMAX', n1)
            if (questi .eq. 'NB_SM_MAILLA') then
                repi = zi(jsssa-1+n1-2)
            else if (questi .eq. 'NB_SS_ACTI') then
                repi = zi(jsssa-1+n1-1)
            else if (questi .eq. 'NB_NL_MAILLA') then
                repi = zi(jsssa-1+n1)
            end if
        end if

    else if (questi .eq. 'NB_NO_MAILLA') then
        call jeveuo(nomob//'.LGRF', 'L', jlgrf)
        call dismma(questi, zk8(jlgrf), repi, repk, ierd)

    else if (questi .eq. 'NB_MA_MAILLA') then
        call jeveuo(nomob//'.LGRF', 'L', jlgrf)
        call dismma(questi, zk8(jlgrf), repi, repk, ierd)

    else if (questi .eq. 'DIM_GEOM') then
        repi = 0
        ige2 = 0
        call jeexin(nomob//'.LIEL', iexi)
        if (iexi .gt. 0) then
            call jelira(nomob//'.LIEL', 'NUTIOC', nbgrel)
            dimge(1) = 0
            dimge(2) = 0
            dimge(3) = 0
            melang = .false.
            do igrel = 1, nbgrel
                call jeveuo(jexnum(nomob//'.LIEL', igrel), 'L', vi=liel)
                call jelira(jexnum(nomob//'.LIEL', igrel), 'LONMAX', lielSize)
                elemTypeNume = liel(lielSize)
                call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
                call dismte(questi, elemTypeName, ige1, repk, ierd)
                ASSERT((ige1 .ge. 0) .and. (ige1 .le. 3))
                if ((ige2 .eq. 0) .and. (ige1 .ne. 0)) ige2 = ige1
                if ((ige1*ige2 .gt. 0) .and. (ige1 .ne. ige2)) melang = .true.
                if (ige1 .gt. 0) dimge(ige1) = 1
            end do
            if (melang) then
                ige3 = +100*dimge(1)
                ige3 = ige3+10*2*dimge(2)
                ige3 = ige3+1*3*dimge(3)
                ige2 = ige3
            end if
        end if
!        -- SI IL EXISTE DES MACRO-ELEMENTS, IL FAUT EN TENIR COMPTE :
        call jeexin(nomob//'.SSSA', iexi2)
        if (iexi2 .gt. 0) then
            call jelira(nomob//'.SSSA', 'LONMAX', n1)
            call jeveuo(nomob//'.SSSA', 'L', jsssa)
            call jeveuo(nomob//'.LGRF', 'L', jlgrf)
            mailla = zk8(jlgrf-1+1)
            call jeveuo(mailla//'.NOMACR', 'L', jnomac)
            nbsm = n1-3
            do ism = 1, nbsm
                if (zi(jsssa-1+ism) .eq. 1) then
                    nomacr = zk8(jnomac-1+ism)
                    call dismml(questi, nomacr, ige1, repk, ierd)
                    ASSERT(ige1 .ge. 0 .and. ige1 .le. 123)
                    if (ige2 .ne. ige1) then
                        ige2 = dimge1(ige2, ige1)
                    end if
                end if
            end do
        end if
        repi = ige2

    else if (questi .eq. 'NB_GREL') then
        call jeexin(nomob//'.LIEL', iexi)
        if (iexi .gt. 0) then
            call jelira(nomob//'.LIEL', 'NUTIOC', repi)
        else
            repi = 0
        end if

    else if (questi .eq. 'NB_MA_SUP') then
        call jeexin(nomob//'.NEMA', iexi)
        if (iexi .gt. 0) then
            call jelira(nomob//'.NEMA', 'NUTIOC', repi)
        else
            repi = 0
        end if

    else if (questi .eq. 'NB_NO_SUP') then
        call jeveuo(nomob//'.NBNO', 'L', vi=nbno)
        repi = nbno(1)

    else if (questi .eq. 'TYPE_LAGR') then
        call jeveuo(nomob//'.LGRF', 'L', jlgrf)
        repk = zk8(jlgrf-1+3)

    else if (questi .eq. 'PHENOMENE') then
        call jelira(nomob//'.LGRF', 'DOCU', cval=phenom)
        if (phenom(1:4) .eq. 'MECA') then
            repk = 'MECANIQUE'
        else if (phenom(1:4) .eq. 'THER') then
            repk = 'THERMIQUE'
        else if (phenom(1:4) .eq. 'ACOU') then
            repk = 'ACOUSTIQUE'
        else
            WRITE (6, *) "DISMOI: ", nomob, phenom
            ierd = 1
            goto 99
        end if

    else if (questi .eq. 'LINE_QUAD') then
        call jeexin(nomob//'.LIEL', iexi)
        if (iexi .gt. 0) then
            call jelira(nomob//'.LIEL', 'NUTIOC', nbgrel)
            call jeveuo('&CATA.TE.TYPEMA', 'L', vk8=typema)
            mail_quad = .false.
            mail_line = .false.
            do igrel = 1, nbgrel
                call jeveuo(jexnum(nomob//'.LIEL', igrel), 'L', vi=liel)
                call jelira(jexnum(nomob//'.LIEL', igrel), 'LONMAX', lielSize)
                elemTypeNume = liel(lielSize)
                cellTypeName = typema(elemTypeNume)
                call dismtm('DIM_TOPO', cellTypeName, ndime, k8bid, ibid)
                if (ndime .ne. 0) then
                    lret = ismali(cellTypeName)
                    if (lret) mail_line = .true.
                    if (.not. lret) mail_quad = .true.
                end if
            end do
            if (mail_line .and. mail_quad) then
                repk = 'LINE_QUAD'
            else if (mail_line) then
                repk = 'LINE'
            else if (mail_quad) then
                repk = 'QUAD'
            else
                ierd = 1
            end if
        end if

    else
        ierd = 1

    end if
!
99  continue
!
    repkz = repk
    call jedema()
end subroutine
