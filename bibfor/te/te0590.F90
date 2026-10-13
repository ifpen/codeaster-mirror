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
subroutine te0590(option, nomte)
!
! aslint: disable=C0110
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nifilg.h"
#include "asterfort/nifipd.h"
#include "asterfort/nifism.h"
#include "asterfort/niinit.h"
#include "asterfort/nmtstm.h"
#include "asterfort/tecach.h"
#include "asterfort/tgveri_use.h"
#include "asterfort/tgverm.h"
#include "MeshTypes_type.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_INCO_UPG, AXIS_INCO_UPG, D_PLAN_INCO_UPG
!
! Options: FULL_MECA, FULL_MECA_ELAS
!          RAPH_MECA
!          RIGI_MECA_ELAS, RIGI_MECA_TANG
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    aster_logical, parameter :: lDebug = ASTER_FALSE
    integer(kind=8) :: ndim, nnod, nnog, nnop, npg, nbElrefe
    integer(kind=8) :: icoret, codret, iret
    integer(kind=8) :: iw, ivfd, ivfg, ivfp, idfd, idfg
    integer(kind=8) :: jtab(7), lgpg
    integer(kind=8) :: vu(3, MT_NNOMAX3D), vg(MT_NNOMAX3D), vp(MT_NNOMAX3D), vpi(3, MT_NNOMAX3D)
    integer(kind=8) :: jvGeom, jvMaterc, icontm, ivarim, iuse
    integer(kind=8) :: jvInstmr, jvInstpr, jvDeplmr, jvDeplpr, jvCarcri, ivarix
    integer(kind=8) :: ivectu, icontp, ivarip, imatuu
    integer(kind=8) :: nddl, ia, ja
    character(len=8) :: listElrefe(10), typmod(2)
    aster_logical :: matsym
    character(len=16) :: defoComp, relaComp, typeComp
    aster_logical :: lVect, lMatr, lVari, lSigm, lMatrPred
!     POUR TGVERI
    real(kind=8) :: sdepl(135), svect(135), scont(6*MT_NNOMAX3D)
    real(kind=8) :: epsilo, epsilp, epsilg
    real(kind=8) :: tab_out(3*MT_NNOMAX3D*3*MT_NNOMAX3D)
    real(kind=8), pointer :: varia(:) => null(), smatr(:) => null()
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: na, os, nb, ib, kk
    blas_int :: b_incx, b_incy, b_n
    type(Behaviour_Integ) :: BEHInteg
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    icontp = 1
    ivarip = 1
    imatuu = 1
    ivectu = 1
    ivarix = 1
    codret = 0
    matsym = ASTER_TRUE

! - List of ELREFE
    call elref2(nomte, 10, listElrefe, nbElrefe)
    ASSERT(nbElrefe .ge. 3)

! - Get shape functions
    call elrefe_info(elrefe=listElrefe(3), fami=fami, nno=nnop, jvf=ivfp)
    call elrefe_info(elrefe=listElrefe(2), fami=fami, nno=nnog, jvf=ivfg, jdfde=idfg)
    call elrefe_info(elrefe=listElrefe(1), fami=fami, ndim=ndim, nno=nnod, npg=npg, &
                     jpoids=iw, jvf=ivfd, jdfde=idfd)
    nddl = nnod*ndim+nnog+nnop

! - Modelling
    typmod = ' '
    if (ndim .eq. 2 .and. lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS'
    else if (ndim .eq. 2 .and. lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN'
    else if (ndim .eq. 3) then
        typmod(1) = '3D'
    else
        ASSERT(ASTER_FALSE)
    end if
    typmod(2) = ' '

! - Get index of dof
    call niinit(typmod, ndim, nnod, nnog, nnop, &
                0, vu, vg, vp, vpi)

! - Get input fields
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PINSTMR', 'L', jvInstmr)
    call jevech('PINSTPR', 'L', jvInstpr)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDEPLMR', 'L', jvDeplmr)
    call jevech('PDEPLPR', 'L', jvDeplpr)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nnod, jvGeom, materPara%lcsPara)

! - Get behaviour parameters
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jvCarcri)

! - Get parameters for behaviour
    defoComp = compor(DEFO)
    relaComp = compor(RELA_NAME)
    typeComp = compor(INCRELAS)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHInteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)
    lMatrPred = option .eq. 'RIGI_MECA_TANG'

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jvCarcri), &
                              zr(jvInstmr), zr(jvInstpr), &
                              materPara, BEHInteg)

! - Get output fields
    if (lMatr) then
        call nmtstm(zr(jvCarcri), imatuu, matsym)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', icoret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        call jevech('PVARIMP', 'L', ivarix)
        b_n = to_blas_int(npg*lgpg)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
    end if
!
    call tgveri_use(option, zr(jvCarcri), compor, iuse)
    if (iuse == 1) then
        allocate (varia(2*135*135))
        allocate (smatr(135*135))
    end if
!
100 continue
! - PETITES DEFORMATIONS
    if (defoComp(1:6) .eq. 'PETIT ') then
        call nifipd(BEHInteg, &
                    ndim, nnod, nnog, nnop, npg, &
                    iw, zr(ivfd), zr(ivfg), zr(ivfp), idfd, &
                    vu, vg, vp, zr(jvGeom), typmod, &
                    option, compor, lgpg, zr(jvCarcri), &
                    zr(jvInstmr), zr(jvInstpr), zr(jvDeplmr), zr(jvDeplpr), &
                    zr(icontm), zr(ivarim), zr(icontp), zr(ivarip), &
                    lMatr, lVect, zr(ivectu), zr(imatuu), &
                    codret)

    else if (defoComp .eq. 'GDEF_LOG') then
        call nifilg(BEHInteg, &
                    ndim, nnod, nnog, nnop, npg, &
                    iw, zr(ivfd), zr(ivfg), zr(ivfp), idfd, &
                    vu, vg, vp, zr(jvGeom), typmod, &
                    option, compor, lgpg, zr(jvCarcri), &
                    zr(jvInstmr), zr(jvInstpr), zr(jvDeplmr), zr(jvDeplpr), &
                    zr(icontm), zr(ivarim), zr(icontp), zr(ivarip), lMatr, &
                    lVect, lSigm, lVari, zr(ivectu), zr(imatuu), &
                    matsym, codret)

    else if (defoComp .eq. 'SIMO_MIEHE') then
        call nifism(BEHInteg, &
                    ndim, nnod, nnog, nnop, npg, &
                    iw, zr(ivfd), zr(ivfg), zr(ivfp), idfd, &
                    idfg, vu, vg, vp, zr(jvGeom), &
                    typmod, option, compor, lgpg, &
                    zr(jvCarcri), zr(jvInstmr), zr(jvInstpr), zr(jvDeplmr), zr(jvDeplpr), &
                    zr(icontm), zr(ivarim), zr(icontp), zr(ivarip), &
                    lMatr, lVect, lMatrPred, zr(ivectu), zr(imatuu), &
                    codret)
    else
        call utmess('F', 'ELEMENTS3_16', sk=defoComp)
    end if
!
    if (codret .ne. 0) goto 200
!       Calcul eventuel de la matrice TGTE par PERTURBATION
    call tgverm(option, zr(jvCarcri), compor, nnod, nnog, &
                nnop, zr(jvGeom), ndim, nddl, zr(jvDeplpr), &
                sdepl, vu, vg, vp, zr(ivectu), &
                svect, ndim*2*npg, zr(icontp), scont, npg*lgpg, &
                zr(ivarip), zr(ivarix), zr(imatuu), smatr, matsym, &
                epsilo, epsilp, epsilg, varia, iret)
    if (iret .ne. 0) goto 100
!
200 continue

! - Save return code
    if (lSigm) then
        zi(icoret) = codret
    end if
!
    if (lDebug) then
        if (lMatr) then
            write (6, *) 'MATRICE TANGENTE'
            if (matsym) then
                do ia = 1, nddl
                    write (6, '(108(1X,E11.4))') (zr(imatuu+(ia*(ia-1)/2)+ja-1), ja=1, ia)
                end do
            else
!
! - TERME K:UU      KUU(NDIM,NNO1,NDIM,NNO1)
!
                write (6, *) 'KUU'
                ja = 1
                do na = 1, nnod
                    do ia = 1, ndim
                        os = (vu(ia, na)-1)*nddl
                        do nb = 1, nnod
                            do ib = 1, ndim
                                kk = os+vu(ib, nb)
                                tab_out(ja) = zr(imatuu+kk-1)
                                ja = ja+1
                            end do
                        end do
                    end do
                end do
                do ia = 1, nnod*ndim
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nnod*ndim+ja), ja=1, nnod*ndim)
                end do
!
! - TERME K:GG      KGG(NNO2,NNO2)
!
                write (6, *) 'KGG'
                ja = 1
                do na = 1, nnog
                    os = (vg(na)-1)*nddl
                    do ia = 1, nnog
                        kk = os+vg(ia)
                        tab_out(ja) = zr(imatuu+kk-1)
                        ja = ja+1
                    end do
                end do
                do ia = 1, nnog
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nnog+ja), ja=1, nnog)
                end do
!
! - TERME K:UP      KUP(NDIM,NNO1,NNO3)
!
                write (6, *) 'KUP'
                ja = 1
                do na = 1, nnod
                    do ia = 1, ndim
                        os = (vu(ia, na)-1)*nddl
                        do nb = 1, nnop
                            kk = os+vp(nb)
                            tab_out(ja) = zr(imatuu+kk-1)
                            ja = ja+1
                        end do
                    end do
                end do
                do ia = 1, nnod*ndim
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nnop+ja), ja=1, nnop)
                end do
!
! - TERME K:PU      KPU(NDIM,NNO3,NNO1)
!
                write (6, *) 'KPU'
                ja = 1
                do ia = 1, nnop
                    os = (vp(ia)-1)*nddl
                    do nb = 1, nnod
                        do ib = 1, ndim
                            kk = os+vu(ib, nb)
                            tab_out(ja) = zr(imatuu+kk-1)
                            ja = ja+1
                        end do
                    end do
                end do
                do ia = 1, nnop
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nnod*ndim+ja), ja=1, nnod*ndim)
                end do
!
! - TERME K:UG      KUG(NDIM,NNO1,NNO2)
!
                write (6, *) 'KUG'
                ja = 1
                do na = 1, nnod
                    do ia = 1, ndim
                        os = (vu(ia, na)-1)*nddl
                        do nb = 1, nnog
                            kk = os+vg(nb)
                            tab_out(ja) = zr(imatuu+kk-1)
                            ja = ja+1
                        end do
                    end do
                end do
                do ia = 1, nnod*ndim
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nnog+ja), ja=1, nnog)
                end do
!
! - TERME K:GU      KGU(NDIM,NNO2,NNO1)
!
                write (6, *) 'KGU'
                ja = 1
                do ia = 1, nnog
                    os = (vg(ia)-1)*nddl
                    do nb = 1, nnod
                        do ib = 1, ndim
                            kk = os+vu(ib, nb)
                            tab_out(ja) = zr(imatuu+kk-1)
                            ja = ja+1
                        end do
                    end do
                end do
                do ia = 1, nnog
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nnod*ndim+ja), ja=1, nnod*ndim)
                end do
!
! - TERME K:PG      KPG(NNO3,NNO2)
!
                write (6, *) 'KPG'
                ja = 1
                do ia = 1, nnop
                    os = (vp(ia)-1)*nddl
                    do ib = 1, nnog
                        kk = os+vg(ib)
                        tab_out(ja) = zr(imatuu+kk-1)
                        ja = ja+1
                    end do
                end do
                do ia = 1, nnop
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nnog+ja), ja=1, nnog)
                end do
!
! - TERME K:GP      KPG(NNO2,NNO3)
!
                write (6, *) 'KGP'
                ja = 1
                do ia = 1, nnog
                    os = (vg(ia)-1)*nddl
                    do ib = 1, nnop
                        kk = os+vp(ib)
                        tab_out(ja) = zr(imatuu+kk-1)
                        ja = ja+1
                    end do
                end do
                do ia = 1, nnog
                    write (6, '(108(1X,E11.4))') (tab_out((ia-1)*nnop+ja), ja=1, nnop)
                end do
            end if
        end if
        if (lVect) then
            write (6, *) 'FORCE INTERNE'
            write (6, '(108(1X,E11.4))') (zr(ivectu+ja-1), ja=1, nddl)
        end if
    end if
!
! - Free large arrays
    if (iuse == 1) then
        deallocate (smatr)
        deallocate (varia)
    end if
!
end subroutine
