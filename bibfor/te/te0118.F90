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

subroutine te0118(option, nomte)
    implicit none
    character(len=16) :: option, nomte
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dismoi.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/indent.h"
#include "asterfort/iselli.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/provec.h"
#include "asterfort/reeref.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xkamat.h"
#include "asterfort/xnormv.h"
#include "MeshTypes_type.h"
!
!
! Calcul du taux de restitution d'energie elementaire sur les
! elements de bord 3D XFEM
!
! option : 'CALC_G' (charges reelles)
!
! in / option : nom de l'option
! in / nomte  : nom du type element
!
! ======================================================================
!
!   nno_par_max = nombre max de noeuds pour l'element de reference parent
!   -> te dedie aux elements de bords, max atteint en 3D pour les QUAD8
    integer(kind=8) :: nno_par_max
    parameter(nno_par_max=8)
!
!   nno_se_max = nombre max de noeuds pour le sous-element de reference
!   -> te dedie aux elements de bords, max atteint en 3D pour les TRIA6
    integer(kind=8) :: nno_se_max
    parameter(nno_se_max=6)
!
    integer(kind=8) :: ndime, ndim, nnop, nnops, cpt, nno, nnos, npg, nfe, nfh
    integer(kind=8) :: nfiss, ncomp, ncompn, nse, hea_se
    integer(kind=8) :: ipoids, ivf, idfde, ipres, iadzi, iazk24, igeom, idepl
    integer(kind=8) :: jlsn, jlst, jpintt, jcnset, jheavt, jlonch
    integer(kind=8) :: jpmilt, jheavn, iforc
    integer(kind=8) :: ino, j, ise, ifiss, in, inop, kpg, ig
    integer(kind=8) :: nlong_ddl, ddld, ddls, ddlm, deca
    integer(kind=8) :: ithet, igthet
    integer(kind=8) :: ier, iret, irese, jtab(7)
    real(kind=8) :: th1, th2, dth1d1, dth2d2, divt, pres, tcla
    real(kind=8) :: r8pre, sum_teth, sum_forc, r8bit2(2)
    real(kind=8) :: poids, norme, vf
    real(kind=8) :: fk(MT_NNOMAX3D, 3, 3), ka, mu
    integer(kind=8) :: alp, jstno, imate, jbaslo
    real(kind=8) :: coorse(3*nno_se_max), geoloc(2*nno_par_max)
    real(kind=8) :: td1(3), td2(3), nd(3), xg(3), xg_loc(2), he(1), oprim(3)
    real(kind=8) :: depla(3), ff(nno_par_max), coorse_loc(2*nno_par_max)
    real(kind=8) :: dford1(3), dford2(3), forcg(3), dfor(3), forc
    real(kind=8) :: dfdi_loc(nno_par_max, 2)
    character(len=8) :: elrefp, elrese(4), elref, enr, noma
    aster_logical :: axi
    data elrese/'SE2', 'TR3', 'SE3', 'TR6'/
!
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! - prealables
! ----------------------------------------------------------------------
!
! - element de reference parent
    call elref1(elrefp)
    call elrefe_info(fami='RIGI', ndim=ndime, nno=nnop, nnos=nnops)
    ASSERT(nnop .le. nno_par_max)
    ASSERT(ndime .eq. 2)
    axi = lteatt('AXIS', 'OUI')
!
! - dimension de l'espace
    call tecael(iadzi, iazk24, noms=0)
    noma = zk24(iazk24) (1:8)
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
    ASSERT(ndim .eq. 3)
!
! - sous-element de reference
    if (.not. iselli(elrefp)) then
        irese = 2
    else
        irese = 0
    end if
    elref = elrese(ndime+irese)
    call elrefe_info(elrefe=elref, fami='XCON', nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(nno .le. nno_se_max)
!
! - on interdit le multi-heaviside
    call tecach('OOO', 'PHEAVTO', 'L', iret, nval=7, itab=jtab)
    ncomp = jtab(2)
    nfiss = jtab(7)
    ASSERT(nfiss .eq. 1)
!
! - initialisation des dimensions des ddls xfem
! - il ne faut pas appeler xteini car elle ne gere pas les elements de bord
    nfe = 0
    nfh = 0
    call teattr('S', 'XFEM', enr, ier)
    if (enr(1:2) .eq. 'XH') then
        nfh = 1
    end if
    if (enr(1:2) .eq. 'XT' .or. enr(3:3) .eq. 'T') then
        nfe = 1
    end if
    ASSERT((nfe .gt. 0) .or. (nfh .gt. 0))
!
    if (nfe .gt. 0) call jevech('PSTANO', 'L', jstno)
!
! - on interdit tout ddl autre que le deplacement (contact...)
    call tecach('OOO', 'PDEPLAR', 'L', iret, nval=7, itab=jtab)
    ddld = ndim*(1+nfh+nfe)
    nlong_ddl = jtab(2)
    ASSERT(nlong_ddl .eq. nnop*ddld)
    ddls = ddld
    ddlm = ddld
!
! - initialisation terme classique
    tcla = 0.d0
!
! - precision utilisee pour tester la nullite d'un reel
    r8pre = r8prem()
!
! - pas de calcul de G pour les elements ou la valeur de theta est nulle
    call jevech('PTHETAR', 'L', ithet)
    call jevech('PGTHETA', 'E', igthet)
    cpt = 0
    sum_teth = 0.d0
    do inop = 1, nnop
        do j = 1, ndim
            sum_teth = sum_teth+abs(zr(ithet-1+ndim*(inop-1)+j))
        end do
        if (abs(sum_teth) .lt. r8pre) then
            cpt = cpt+1
        end if
    end do
    if (cpt .eq. nnop) then
        goto 999
    end if
!
! ----------------------------------------------------------------------
! - recuperation des entrees / sorties
! ----------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLAR', 'L', idepl)
!
! - si la pression n'est connue sur aucun noeud, on la choisit nulle
    call jevecd('PPRESSR', ipres, 0.d0)
!
! - la programmation relative aux chargements FORCE_FACE n'est pas
! - encore realisee dans ce te, on interdit donc la presence de ce
! - type de charge
    call jevech('PFR2D3D', 'L', iforc)
    cpt = 0
    sum_forc = 0.d0
    do inop = 1, nnop
        do j = 1, ndim
            sum_forc = sum_forc+abs(zr(iforc-1+ndim*(inop-1)+j))
        end do
        if (abs(sum_forc) .lt. r8pre) then
            cpt = cpt+1
        end if
    end do
    if (cpt .ne. nnop) call utmess('F', 'XFEM_98')
!
! - parametres propres a xfem
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PHEAVTO', 'L', jheavt)
    if (nfe .gt. 0) then
        call jevech('PMATERC', 'L', imate)
        call jevech('PBASLOR', 'L', jbaslo)
    end if
!
! - donnees topologiques des fonction Heaviside
    if (enr(1:2) .eq. 'XH') then
        call jevech('PHEA_NO', 'L', jheavn)
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, itab=jtab)
        ncompn = jtab(2)/jtab(3)
    end if
!
! - parametres propres aux elements xfem quadratiques
    if (.not. iselli(elref)) then
        call jevech('PPMILTO', 'L', jpmilt)
    end if
!
! ----------------------------------------------------------------------
! - Definition du repere local de la face supposee plane (element lineaire)
! ----------------------------------------------------------------------
!
! - tangentes et normale a la face
    td1(:) = 0.d0
    td2(:) = 0.d0
    do j = 1, ndim
        td1(j) = zr(igeom+ndim*(2-1)+j-1)-zr(igeom+ndim*(1-1)+j-1)
        td2(j) = zr(igeom+ndim*(3-1)+j-1)-zr(igeom+ndim*(1-1)+j-1)
    end do
!
! - calcul d'une base orthomormee 'Bprime' = (td1, td2, nd)
!   rq : on norme td1 et td2 avant de faire les produits vectoriels
!        pour eviter les pb si on a des "petites" mailles
    call xnormv(ndim, td1, norme)
    ASSERT(norme .gt. r8pre)
    call xnormv(ndim, td2, norme)
    ASSERT(norme .gt. r8pre)
    call provec(td1, td2, nd)
    call xnormv(ndim, nd, norme)
    ASSERT(norme .gt. r8pre)
    call provec(nd, td1, td2)
    call xnormv(ndim, td2, norme)
    ASSERT(norme .gt. r8pre)
!
! - origine 'oprim' du repere local (1er noeud)
    do j = 1, ndim
        oprim(j) = zr(igeom+ndim*(1-1)+j-1)
    end do
!
! - coordonnees des noeuds de l'element parent dans le
! - repere 2D local (oprim, (td1, td2))
    geoloc(:) = 0.d0
    do inop = 1, nnop
        do j = 1, ndim
            geoloc(2*inop-1) = geoloc(2*inop-1) &
                               +(zr(igeom+ndim*(inop-1)+j-1)-oprim(j)) &
                               *td1(j)
            geoloc(2*inop) = geoloc(2*inop) &
                             +(zr(igeom+ndim*(inop-1)+j-1)-oprim(j)) &
                             *td2(j)
        end do
    end do
!
! ----------------------------------------------------------------------
! - Boucle d'integration sur les nse sous-elements
! ----------------------------------------------------------------------
!
! - recuperation de la subdivision des elements en nse sous-elements
    nse = zi(jlonch-1+1)
!
    do ise = 1, nse
!
! ----- coordonnees des noeuds du sous-element dans le repere reel
        coorse(:) = 0.d0
        do in = 1, nno
            ino = zi(jcnset-1+nno*(ise-1)+in)
            do j = 1, ndim
                if (ino .lt. 1000) then
                    coorse(ndim*(in-1)+j) = zr(igeom-1+ndim*(ino-1)+j)
                else if (ino .gt. 1000 .and. ino .lt. 2000) then
                    coorse(ndim*(in-1)+j) = zr(jpintt-1+ndim*(ino-1000- &
                                                              1)+j)
                else if (ino .gt. 2000 .and. ino .lt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-2000- &
                                                              1)+j)
                else if (ino .gt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-3000- &
                                                              1)+j)
                end if
            end do
        end do
!
! ----- coordonnees des noeuds du sous-element dans le repere 2D local
        coorse_loc(:) = 0.d0
        do in = 1, nno
            do j = 1, ndim
                coorse_loc(2*in-1) = coorse_loc(2*in-1) &
                                     +(coorse(ndim*(in-1)+j)-oprim(j)) &
                                     *td1(j)
                coorse_loc(2*in) = coorse_loc(2*in) &
                                   +(coorse(ndim*(in-1)+j)-oprim(j)) &
                                   *td2(j)
            end do
        end do
!
! ----- fonction heaviside constante sur le sous-element et par fissure
        do ifiss = 1, nfiss
            he(ifiss) = zi(jheavt-1+ncomp*(ifiss-1)+ise)
        end do
!
! ----- calcul de l'identifiant du sous-element
        hea_se = xcalc_code(nfiss, he_real=[he])
!
! ----------------------------------------------------------------------
! ----- Boucle sur les points de Gauss du sous-element
! ----------------------------------------------------------------------
!
        do kpg = 1, npg
!
! --------- calcul du poids : poids = poids * jacobien
            call dfdm2d(nno, kpg, ipoids, idfde, coorse_loc, poids)
!
! --------- coordonnees du point de Gauss dans le repere global : xg
            xg(:) = 0.d0
            do in = 1, nno
                vf = zr(ivf-1+nno*(kpg-1)+in)
                do j = 1, ndim
                    xg(j) = xg(j)+vf*coorse(ndim*(in-1)+j)
                end do
            end do
!
! --------- coordonnees du point de Gauss dans le repere
! --------- 2D local (oprim, (td1, td2)) : xg_loc
            xg_loc(:) = 0.d0
            do j = 1, ndim
                xg_loc(1) = xg_loc(1)+(xg(j)-oprim(j))*td1(j)
                xg_loc(2) = xg_loc(2)+(xg(j)-oprim(j))*td2(j)
            end do
!
! --------- derivees des fonctions de formes dans le repere local
            call reeref(elrefp, nnop, geoloc, xg_loc, ndime, &
                        r8bit2, ff, dfdi=dfdi_loc(1:nnop, 1:ndime))
!
! --------- calcul des fonctions d'enrichissement
            if (nfe .gt. 0) then
                call xkamat(zi(imate), ndim, axi, ka, mu, famiz='XCON')
                call xcalfev_wrap(ndim, nnop, zr(jbaslo), zi(jstno), he(1), &
                                  zr(jlsn), zr(jlst), zr(igeom), ka, mu, ff, fk)
            end if
!
! --------- calcul de l'approximation du deplacement
            depla(:) = 0.d0
            do inop = 1, nnop
                call indent(inop, ddls, ddlm, nnops, deca)
                cpt = 0
!               ddls classiques
                do j = 1, ndim
                    cpt = cpt+1
                    depla(j) = depla(j)+ff(inop)*zr(idepl-1+deca+cpt)
                end do
!               ddls heaviside
                do ig = 1, nfh
                    do j = 1, ndim
                        cpt = cpt+1
                        depla(j) = depla(j)+xcalc_heav(zi(jheavn-1+ncompn*(inop-1)+ig), &
                                                       hea_se, &
                                                       zi(jheavn-1+ncompn*(inop-1)+ncompn)) &
                                   *ff(inop)*zr(idepl-1+deca+cpt)
                    end do
                end do
!               ddls enrichis en fond de fissure
                do alp = 1, nfe*ndim
                    cpt = cpt+1
                    do j = 1, ndim
                        depla(j) = depla(j)+fk(inop, alp, j)*zr(idepl-1+deca+cpt)
                    end do
                end do
            end do
!
! --------- calcul de la pression
            pres = 0.d0
            do inop = 1, nnop
                pres = pres+zr(ipres-1+inop)*ff(inop)
            end do
!
! --------- calcul du terme surfacque
            th1 = 0.d0
            th2 = 0.d0
            dth1d1 = 0.d0
            dth2d2 = 0.d0
!
            do j = 1, ndim
                dford1(j) = 0.d0
                dford2(j) = 0.d0
                dfor(j) = 0.d0
                forcg(j) = 0.d0
            end do
!
            do inop = 1, nnop
                do j = 1, ndim
                    th1 = th1+ff(inop)*zr(ithet+3*(inop-1)+j-1)*td1(j)
                    th2 = th2+ff(inop)*zr(ithet+3*(inop-1)+j-1)*td2(j)
                    dth1d1 = dth1d1+zr(ithet-1+3*(inop-1)+j)*td1(j)*dfdi_loc(inop, 1)
                    dth2d2 = dth2d2+zr(ithet-1+3*(inop-1)+j)*td2(j)*dfdi_loc(inop, 2)
                end do
            end do
!
            do j = 1, ndim
                dfor(j) = dfor(j)+dford1(j)*th1+dford2(j)*th2
            end do
!
            divt = dth1d1+dth2d2
!
            do j = 1, ndim
                forc = forcg(j)-pres*nd(j)
                tcla = tcla+poids*(forc*divt+dfor(j))*depla(j)
            end do
!
! ----------------------------------------------------------------------
! ----- Fin boucle sur les points de Gauss du sous-element
! ----------------------------------------------------------------------
!
        end do
!
! ----------------------------------------------------------------------
! - Fin boucle d'integration sur les nse sous-elements
! ----------------------------------------------------------------------
!
    end do
!
999 continue
!
    zr(igthet) = tcla
!
end subroutine
