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

subroutine nmextj(field_type, nb_elem, nb_cmp, list_cmp, type_extr_cmp, type_sele_cmp, &
                  poin_nume, spoi_nume, nb_vale, i_elem, elem_nume, &
                  jcesd, jcesv, jcesl, jcesc, vale_resu)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxliis.h"
#include "asterfort/nmextv.h"
!
!
    character(len=24), intent(in) :: field_type
    integer(kind=8), intent(in) :: nb_elem
    integer(kind=8), intent(in) :: nb_cmp
    character(len=24), intent(in) :: list_cmp
    character(len=8), intent(in) :: type_extr_cmp
    character(len=8), intent(in) :: type_sele_cmp
    integer(kind=8), intent(in) :: poin_nume
    integer(kind=8), intent(in) :: i_elem
    integer(kind=8), intent(in):: elem_nume
    integer(kind=8), intent(in) :: spoi_nume
    integer(kind=8), intent(in) :: jcesd
    integer(kind=8), intent(in) :: jcesv
    integer(kind=8), intent(in) :: jcesl
    integer(kind=8), intent(in) :: jcesc
    integer(kind=8), intent(out) :: nb_vale
    real(kind=8), intent(out) :: vale_resu(*)
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Field extraction datastructure
!
! Extract value(s) at Gauss point
!
! --------------------------------------------------------------------------------------------------
!
! In  i_elem           : index of element
! In  elem_nume        : index of element
! In  poin_nume        : index of point
! In  spoi_nume        : index of subpoint
! In  jcesd            : Jeveux adress to CHAM_ELEM_S.CESD object
! In  jcesv            : Jeveux adress to CHAM_ELEM_S.CESV object
! In  jcesl            : Jeveux adress to CHAM_ELEM_S.CESL object
! In  jcesc            : Jeveux adress to CHAM_ELEM_S.CESC object
! In  field_type       : type of field (name in results datastructure)
! In  nb_cmp           : number of components
! In  list_cmp         : name of object contains list of components
! In  type_extr_cmp    : type of extraction for components
! In  type_sele_cmp    : type of selection for components NOM_CMP or NOM_VARI
! Out vale_resu        : list of result values
! Out nb_vale          : number of result values (one if function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_para_maxi = 20
    character(len=8) :: v_cmp_name(nb_para_maxi)
    real(kind=8) :: v_cmp_vale(nb_para_maxi)
    integer(kind=8) :: nb_cmp_vale, nb_cmp_maxi
    integer(kind=8) :: i_cmp, iret, i_cmp_vale, i_cmp_maxi, i_cmp_r, i_vari
    integer(kind=8) :: iad
    character(len=8) :: cmp_name, vari_name
    character(len=8), pointer :: v_list_cmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    i_cmp_vale = 1
    nb_cmp_maxi = zi(jcesd+4)
    ASSERT(nb_cmp .le. nb_para_maxi)
!
! - Get name of components
!
    call jeveuo(list_cmp, 'L', vk8=v_list_cmp)
    do i_cmp = 1, nb_cmp
        if (type_sele_cmp .eq. 'NOM_CMP') then
            v_cmp_name(i_cmp) = v_list_cmp(i_cmp)
        elseif (type_sele_cmp .eq. 'NOM_VARI') then
            v_cmp_name(i_cmp) = v_list_cmp(nb_elem*(i_cmp-1)+i_elem)
        else
            ASSERT(.false.)
        end if
    end do
!
! - Get value of components
!
    do i_cmp = 1, nb_cmp
        cmp_name = v_cmp_name(i_cmp)
        if (field_type(1:4) .eq. 'VARI') then
            vari_name = cmp_name(2:8)//' '
            call lxliis(vari_name, i_vari, iret)
        else
            i_vari = 0
        end if
        if (field_type(1:4) .eq. 'VARI') then
            i_cmp_r = i_vari
        else
            do i_cmp_maxi = 1, nb_cmp_maxi
                if (cmp_name .eq. zk8(jcesc-1+i_cmp_maxi)) then
                    i_cmp_r = i_cmp_maxi
                end if
            end do
        end if
        call cesexi('C', jcesd, jcesl, elem_nume, poin_nume, &
                    spoi_nume, i_cmp_r, iad)
        if (iad .gt. 0) then
            v_cmp_vale(i_cmp_vale) = zr(jcesv+iad-1)
            i_cmp_vale = i_cmp_vale+1
        end if
    end do
    nb_cmp_vale = i_cmp_vale-1
!
! - Evaluation
!
    call nmextv(nb_cmp_vale, type_extr_cmp, v_cmp_name, v_cmp_vale, nb_vale, &
                vale_resu)
    ASSERT(nb_vale .le. nb_cmp)
!
end subroutine
