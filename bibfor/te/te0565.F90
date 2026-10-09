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
subroutine te0565(option, nomte)
!
    use Behaviour_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/enelpg.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ltequa.h"
#include "asterfort/nbsigm.h"
#include "asterfort/reeref.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: nomte, option
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: XFEM
!
! Options: ENEL_ELEM, ENTR_ELEM, ENER_TOTALE, INDIC_ENER, INDIC_SEUIL
!
! --------------------------------------------------------------------------------------------------
!
!  OPTION ENEL_ELEM : CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
!  ================   DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!  EN HPP
!   ENELAS =  SOMME_VOLUME((SIG_T*(1/D)*SIG).DV)
!
!        OU  .SIG       EST LE TENSEUR DES CONTRAINTES DE CAUCHY
!            .D         EST LE TENSEUR DE HOOKE
!
!  EN GRANDES DEFORMATIONS SIMO MIEHE POUR ELAS OU VMIS_ISOT
!   ENERLAS = ENERGIE ELASTIQUE SPECIFIQUE
!           = K(0.5(J^2-1)-lnJ)+0.5mu(tr(J^(-2/3)be)-3)
!           SI PRESENCE DE THERMIQUE, ON AJOUTE UNE CORRECTION
!           SPECIFIQUE PRESENTEE DANS LA DOC R
!  EN GRANDES DEFORMATIONS GDEF_LOG
!   ENERELAS = SOMME_VOLUME((T_T*(1/D)*T).DV)
!        OU  .T       EST LE TENSEUR DES CONTRAINTES DU FORMALISME
!            .D         EST LE TENSEUR DE HOOKE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: famiMater = "XFEM"
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8), parameter :: nbsgm = 6
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8) :: idene1
    integer(kind=8) :: jvDBaseFunc, jvSigm, jvVari, jvGeom, jvMaterc, jvTime
    integer(kind=8) :: jvGaussWeight, jvBaseFunc
    integer(kind=8) :: nbsig, nbvari, ndim, nno
    integer(kind=8) :: npg, iret, i, jtab(7)
    real(kind=8) :: enerElas
    real(kind=8) :: welas, wtotal
    real(kind=8) :: sigmEner(nbsgm)
    real(kind=8) :: time
    real(kind=8) :: f(3, 3), r
    character(len=16) :: relaName, defoComp
    integer(kind=8) :: jpintt, jpmilt, jcnset, jlonch
    integer(kind=8) :: nse, nnop
    character(len=8) :: elrefp, enr
    real(kind=8) :: coorse(3*MT_NNOMAX3D), xg(3), xe(3), ff(MT_NNOMAX3D)
    real(kind=8) :: jac
    integer(kind=8) :: irese, idecpg, idebs, idebv
    integer(kind=8) :: j, in, ino, ise, kpg, n
    aster_logical :: largeStrain, axi
    character(len=8), parameter :: elrese(6) = (/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/)
    character(len=8), parameter :: fami(6) = (/'BID ', 'XINT', 'XINT', 'BID ', 'XINT', 'XINT'/)
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    enerElas = zero
    welas = zero
    wtotal = zero

! - Get parent
    call elref1(elrefp)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nnop)

!   SOUS-ELEMENT DE REFERENCE : RECUP DE NNO, NPG ET IVF
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), nno=nno, npg=npg, &
                     jpoids=jvGaussWeight, jvf=jvBaseFunc, jdfde=jvDBaseFunc)
    axi = lteatt('AXIS', 'OUI')
    nbsig = nbsigm()
    ASSERT(nbsig .le. 6)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(famiMater, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - XFEM
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PLONCHA', 'L', jlonch)
    call teattr('S', 'XFEM', enr, iret)
    if ((iret .eq. 0) .and. ltequa(elrefp, enr)) then
        call jevech('PPMILTO', 'L', jpmilt)
    end if

! - Get stresses
    call jevech('PCONTPR', 'L', jvSigm)

! - Get time
    time = zero
    call tecach('NNO', 'PINSTR', 'L', iret, iad=jvTime)
    if (jvTime .ne. 0) then
        time = zr(jvTime)
    end if

! - Get parameters from behaviour (linear and non-linear cases)
    call behaviourGetParameters(relaName, defoComp)
    if ((defoComp .eq. 'SIMO_MIEHE') .or. (defoComp .eq. 'GDEF_LOG')) then
        largeStrain = ASTER_TRUE
    else
        largeStrain = ASTER_FALSE
    end if

! - Internal state variables (only for non-linear)
    jvVari = 1
    nbvari = 0
    call tecach('ONO', 'PVARIPR', 'L', iret, nval=7, itab=jtab)
    if (iret .eq. 0) then
        jvVari = jtab(1)
        nbvari = max(jtab(6), 1)*jtab(7)
    end if

! - RECUPERATION DU CHAMP POUR STOCKER L'ENERGIE ELSTIQUE
    call jevech('PENERD1', 'E', idene1)
!
! --- RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = zi(jlonch-1+1)

! - BOUCLE D'INTEGRATION SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
!       BOUCLE SUR LES SOMMETS DU SOUS-TRIA (DU SOUS-SEG)
        do in = 1, nno
            ino = zi(jcnset-1+nno*(ise-1)+in)
            do j = 1, ndim
                if (ino .lt. 1000) then
                    coorse(ndim*(in-1)+j) = zr(jvGeom-1+ndim*(ino-1)+j)
                else if (ino .gt. 1000 .and. ino .lt. 2000) then
                    coorse(ndim*(in-1)+j) = zr(jpintt-1+ndim*(ino-1000-1)+j)
                else if (ino .gt. 2000 .and. ino .lt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-2000-1)+j)
                else if (ino .gt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-3000-1)+j)
                end if
            end do
        end do

        do kpg = 1, npg
! --------- Initializations of material parameters on current integration point
            call initParaPoin(kpg, ksp, materPara)

! --------- Coordinates of current Gauss point
            xg = 0.d0
            do i = 1, ndim
                do n = 1, nno
                    xg(i) = xg(i)+zr(jvBaseFunc-1+nno*(kpg-1)+n)*coorse(ndim*(n-1)+i)
                end do
            end do

! --------- Shepe function
            call reeref(elrefp, nnop, zr(jvGeom), xg, ndim, xe, ff)

! --------- Compute jacobian
            if (ndim .eq. 2) then
                call dfdm2d(nno, kpg, jvGaussWeight, jvDBaseFunc, coorse, &
                            jac)
            else if (ndim .eq. 3) then
                call dfdm3d(nno, kpg, jvGaussWeight, jvDBaseFunc, coorse, &
                            jac)
            else
                ASSERT(ASTER_FALSE)
            end if

! --------- CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE):
            if (axi) then
                r = 0.d0
                do n = 1, nnop
                    r = r+ff(n)*zr(jvGeom-1+2*(n-1)+1)
                end do
                ASSERT(r .gt. 0d0)
                jac = jac*r
            end if

! --------- DEBUT DE LA ZONE MEMOIRE DE SIG ET VI CORRESPONDANTE
            idecpg = npg*(ise-1)
            idebs = nbsig*idecpg
            idebv = nbvari*idecpg

! --------- TENSEUR DES CONTRAINTES AU POINT D'INTEGRATION COURANT :
            sigmEner = 0.d0
            do i = 1, nbsig
                sigmEner(i) = zr(jvSigm+idebs+(kpg-1)*nbsig+i-1)
            end do

! --------- CALCUL DU GRADIENT DE LA TRANSFORMATION
            ASSERT(.not. largeStrain)
            f = 0.d0
            do i = 1, 3
                f(i, i) = 1.d0
            end do
            call enelpg(materPara, time, &
                        relaName, defoComp, &
                        f, sigmEner, &
                        nbvari, zr(jvVari+(kpg-1)*nbvari), &
                        enerElas)
            welas = welas+enerElas*jac

        end do
    end do
!
    zr(idene1) = welas
!
end subroutine
