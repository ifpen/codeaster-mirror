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

subroutine te0036(option, nomte)
    implicit none
!
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN MECANIQUE
!          CORRESPONDANT A UN CHARGEMENT EN PRESSION REPARTIE
!          SUR DES FACES D'ELEMENTS X-FEM
!          (LA PRESSION PEUT ETRE DONNEE SOUS FORME D'UNE FONCTION)
!
!          OPTIONS : 'CHAR_MECA_PRES_R'
!                    'CHAR_MECA_PRES_F'
!                    'CHAR_MECA_FR2D3D'
!                    'CHAR_MECA_FR1D2D'
!                    'CHAR_MECA_FF2D3D'
!                    'CHAR_MECA_FF1D2D'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/dfdm2b.h"
#include "asterfort/dismoi.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/iselli.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/reeref.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/tefrep.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xhmddl.h"
#include "asterfort/xkamat.h"
#include "asterfort/xnormv.h"
#include "asterfort/xteddl.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
    character(len=8) :: nompar(4), noma, elrefp, elrese(4), enr, lag
    character(len=8) :: elref
    integer(kind=8) :: jpintt, jcnset, jheavt, jlonch, jlsn, jlst, k
    integer(kind=8) :: jpmilt, irese, nfiss, jfisno, jtab(7), ncomp, jheavn, jheavs, ncompn
    integer(kind=8) :: jbaslo, imate
    integer(kind=8) :: alp
    integer(kind=8) :: ibid, ier, ndim, nno, nnop, nnops, npg, nnos, kpg
    integer(kind=8) :: ipoids, ivf, idfde, igeom, ipres, itemps, ires, j
    integer(kind=8) :: nfh, nfe, nse, ise
    integer(kind=8) :: in, ino, iadzi, iazk24, jstno
    integer(kind=8) :: iforc, iret, ig, pos, ndime, nddl, ddls, contac
    real(kind=8) :: xg(4), r
    real(kind=8) :: pres, ff(MT_NNOMAX3D), coorse(18), cosa, sina
    real(kind=8) :: nd(3), norme, rb1(3), rb2(3), cisa
    real(kind=8) :: poids, forrep(3), vf, td(3), rbid(1)
    real(kind=8) :: fk(MT_NNOMAX3D, 3, 3), ka, mu
    aster_logical :: lbid, axi, pre1
    integer(kind=8) :: kk, ddlm, nnopm
    data elrese/'SE2', 'TR3', 'SE3', 'TR6'/
!
!-----------------------------------------------------------------------
!     INITIALISATIONS
!-----------------------------------------------------------------------
!
!     ELEMENT DE REFERENCE PARENT
    call elref1(elrefp)
    call elrefe_info(fami='RIGI', ndim=ndime, nno=nnop, nnos=nnops)
    ASSERT(ndime .eq. 1 .or. ndime .eq. 2)
!
    axi = lteatt('AXIS', 'OUI')

    pre1 = .false.
!     SI ON EST DANS LE CAS HM-XFEM, PRE1=.TRUE.
    if (nomte(1:2) .eq. 'HM') then
        pre1 = .true.
    end if
!
!     DIMENSION DE L'ESPACE
    call tecael(iadzi, iazk24, noms=0)
    noma = zk24(iazk24) (1:8)
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
!
!     ATTENTION, NE PAS CONFONDRE NDIM ET NDIME  !!
!     NDIM EST LA DIMENSION DU MAILLAGE
!     NDIME EST DIMENSION DE L'ELEMENT FINI
!     SUR UN ELET DE BORD, ON A :  NDIM = NDIME + 1
!
!     SOUS-ELEMENT DE REFERENCE
    if (.not. iselli(elrefp)) then
        irese = 2
    else
        irese = 0
    end if
    elref = elrese(ndime+irese)
    call elrefe_info(elrefe=elref, fami='RIGI', nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
!     INITIALISATION DES DIMENSIONS DES DDLS X-FEM
!     IL NE FAUT PAS APPELER XTEINI CAR IL NE GERE PAS LES ELEMENTS
!     DE BORD
!      CALL XTEINI(NOMTE,DDLH,NFE,IBID,IBID,IBID,IBID,IBID)
    nfe = 0
    nfh = 0
    nfiss = 1
    call teattr('S', 'XFEM', enr, ier)
    if (enr(1:2) .eq. 'XH') then
! --- NOMBRE DE FISSURES
        call tecach('NOO', 'PHEAVTO', 'L', iret, nval=7, &
                    itab=jtab)
        ncomp = jtab(2)
        nfiss = jtab(7)
        nfh = 1
        if (enr(1:3) .eq. 'XH2') then
            nfh = 2
        else if (enr(1:3) .eq. 'XH3') then
            nfh = 3
        else if (enr(1:3) .eq. 'XH4') then
            nfh = 4
        end if
    end if
!
    if (enr(1:2) .eq. 'XT' .or. enr(3:3) .eq. 'T') then
        nfe = 1
    end if
!
    ASSERT(nfe .gt. 0 .or. nfh .gt. 0)
!
!-----------------------------------------------------------------------
!     RECUPERATION DES ENTREES / SORTIE
!-----------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', igeom)
!
    if (option .eq. 'CHAR_MECA_PRES_R') then
!
!       SI LA PRESSION N'EST CONNUE SUR AUCUN NOEUD, ON LA PREND=0.
        call jevecd('PPRESSR', ipres, 0.d0)
!
    else if (option .eq. 'CHAR_MECA_PRES_F') then
!
        call jevech('PPRESSF', 'L', ipres)
        call jevech('PINSTR', 'L', itemps)
!
    elseif (option .eq. 'CHAR_MECA_FR2D3D' .or.&
     &        option .eq. 'CHAR_MECA_FR1D2D') then
        if (ndim .eq. 3) then
            call tefrep(option, 'PFR2D3D', iforc)
        else if (ndim .eq. 2) then
            call tefrep(option, 'PFR1D2D', iforc)
        end if
!
    elseif (option .eq. 'CHAR_MECA_FF2D3D' .or.&
     &        option .eq. 'CHAR_MECA_FF1D2D') then
!
        if (ndim .eq. 3) then
            call jevech('PFF2D3D', 'L', iforc)
        else if (ndim .eq. 2) then
            call jevech('PFF1D2D', 'L', iforc)
        end if
        call jevech('PINSTR', 'L', itemps)
!
    end if
!
!     PARAMETRES PROPRES A X-FEM
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PHEAVTO', 'L', jheavt)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PSTANO', 'L', jstno)
    if (nfe .gt. 0) then
        call jevech('PBASLOR', 'L', jbaslo)
        call jevech('PMATERC', 'L', imate)
    end if
!     PROPRE AUX ELEMENTS 1D ET 2D (QUADRATIQUES)
    call teattr('S', 'XFEM', enr, ier)
    if (ier .eq. 0 .and. .not. iselli(elref)) call jevech('PPMILTO', 'L', jpmilt)
    if (nfiss .gt. 1) call jevech('PFISNO', 'L', jfisno)
!
!   LECTURE DES DONNES TOPOLOGIQUE DES FONCTIONS HEAVISIDE
    if (enr(1:2) .eq. 'XH' .or. enr(1:2) .eq. 'XT') then
        call jevech('PHEA_NO', 'L', jheavn)
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                    itab=jtab)
        ncompn = jtab(2)/jtab(3)
        call jevech('PHEA_SE', 'L', jheavs)
    end if
!
    call jevech('PVECTUR', 'E', ires)
!
!     RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = zi(jlonch-1+1)
!
!       BOUCLE D'INTEGRATION SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
        coorse(:) = 0.d0
!       BOUCLE SUR LES NOEUDS DU SOUS-TRIA (DU SOUS-SEG)
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
!-----------------------------------------------------------------------
!         BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELT
!-----------------------------------------------------------------------
! -     CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE)
!       if (axi) then
!            r = 0.d0
!            do ino = 1, nnop
!                r = r + ff(ino)*zr(igeom-1+2*(ino-1)+1)
!            end do
!            ASSERT(r.ge.0d0)
!          ATTENTION : LE POIDS N'EST PAS X R
!          CE SERA FAIT PLUS TARD AVEC JAC = JAC X R
!        endif
!
        do kpg = 1, npg
!
!         CALCUL DU POIDS : POIDS = POIDS DE GAUSS * DET(J)
            if (ndime .eq. 2) then
                kk = 2*(kpg-1)*nno
                call dfdm2b(nno, zr(ipoids-1+kpg), zr(idfde+kk), coorse, &
                            poids, nd)
            else if (ndime .eq. 1) then
                kk = (kpg-1)*nno
                call dfdm1d(nno, zr(ipoids-1+kpg), zr(idfde+kk), coorse, rb1, &
                            rb2(1), poids, cosa, sina)
            end if
!
!         COORDONNÉES RÉELLES GLOBALES DU POINT DE GAUSS
            xg(:) = 0.d0
            do j = 1, nno
                vf = zr(ivf-1+nno*(kpg-1)+j)
                do k = 1, ndim
                    xg(k) = xg(k)+vf*coorse(ndim*(j-1)+k)
                end do
            end do
!
            if (ndime .eq. 1) then
                ASSERT(elref(1:2) .eq. 'SE')
                td(:) = 0.d0
                nd(:) = 0.d0
!         CALCUL DE LA NORMALE AU SEGMENT AU POINT DE GAUSS
                nd(1) = cosa
                nd(2) = sina
                call xnormv(2, nd, norme)
                ASSERT(norme .gt. 0.d0)
!         CALCUL DE LA TANGENTE AU SEGMENT AU POINT DE GAUSS
!                do j = 1, nno
!                   vf=zr(idfde-1+nno*(kpg-1)+j)
!                   do k = 1, ndim
!                       td(k)=td(k)+vf*coorse(ndim*(j-1)+k)
!                   end do
!               end do
                td(1) = -sina
                td(2) = cosa
            else if (ndime .eq. 2) then
!         CALCUL DE LA NORMALE A LA FACE AU POINT DE GAUSS
                call xnormv(3, nd, norme)
                ASSERT(norme .gt. 0.d0)
            end if
!
!           CALCUL DES FONCTIONS D'ENRICHISSEMENT
!           -------------------------------------
!
            call reeref(elrefp, nnop, zr(igeom), xg, ndim, rb2, ff)
!
!       FONCTION D'ENRICHISSEMENT AU POINT DE GAUSS ET LEURS DÉRIVÉES
            if (nfe .gt. 0) then
                call xkamat(zi(imate), ndim, axi, ka, mu, famiz='XCON')
                call xcalfev_wrap(ndim, nnop, zr(jbaslo), zi(jstno), real(zi(jheavt-1+ise), 8), &
                                  zr(jlsn), zr(jlst), zr(igeom), ka, mu, ff, fk)
            end if
!
!         CALCUL DES FORCES REPARTIES SUIVANT LES OPTIONS
!         -----------------------------------------------
!
            forrep(:) = 0.d0
            nompar(1) = 'X'
            nompar(2) = 'Y'
            if (ndim .eq. 3) nompar(3) = 'Z'
            if (ndim .eq. 3) nompar(4) = 'INST'
            if (ndim .eq. 2) nompar(3) = 'INST'
!
            if (option .eq. 'CHAR_MECA_PRES_R') then
!
!             CALCUL DE LA PRESSION AUX POINTS DE GAUSS
                pres = 0.d0
                cisa = 0.d0
                do ino = 1, nnop
                    if (ndim .eq. 3) pres = pres+zr(ipres-1+ino)*ff(ino)
                    if (ndim .eq. 2) then
                        pres = pres+zr(ipres-1+2*(ino-1)+1)*ff(ino)
                        cisa = cisa+zr(ipres-1+2*(ino-1)+2)*ff(ino)
                    end if
                end do
!           ATTENTION AU SIGNE : POUR LES PRESSIONS, IL FAUT UN - DVT
!           CAR LE SECOND MEMBRE SERA ECRIT AVEC UN + (VOIR PLUS BAS)
                do j = 1, ndim
                    forrep(j) = -pres*nd(j)
                end do
                if (ndim .eq. 2) then
                    forrep(1) = forrep(1)-cisa*nd(2)
                    forrep(2) = forrep(2)+cisa*nd(1)
                end if
!
            else if (option .eq. 'CHAR_MECA_PRES_F') then
!
!             VALEUR DE LA PRESSION
                xg(ndim+1) = zr(itemps)
!
                call fointe('FM', zk8(ipres), ndim+1, nompar, xg, &
                            pres, ier)
                if (ndim .eq. 2) call fointe('FM', zk8(ipres+1), ndim+1, nompar, xg, &
                                             cisa, ier)
                do j = 1, ndim
                    forrep(j) = -pres*nd(j)
                end do
                if (ndim .eq. 2) then
                    forrep(1) = forrep(1)-cisa*nd(2)
                    forrep(2) = forrep(2)+cisa*nd(1)
                end if
!
            elseif (option .eq. 'CHAR_MECA_FR2D3D' .or.&
     &            option .eq. 'CHAR_MECA_FR1D2D') then
!
                forrep(:) = 0.d0
                do ino = 1, nnop
                    do j = 1, ndim
                        forrep(j) = forrep(j)+zr(iforc-1+ndim*(ino-1)+j) &
                                    *ff(ino)
                    end do
                end do
!
            elseif (option .eq. 'CHAR_MECA_FF2D3D' .or.&
     &            option .eq. 'CHAR_MECA_FF1D2D') then
!
                xg(ndim+1) = zr(itemps)
                do j = 1, ndim
                    call fointe('FM', zk8(iforc-1+j), ndim+1, nompar, xg, &
                                forrep(j), ier)
                end do
!
            end if
            if (axi) then
                r = 0.d0
                do ino = 1, nnop
                    r = r+ff(ino)*zr(igeom-1+2*(ino-1)+1)
                end do
                ASSERT(r .ge. 0d0)
                poids = poids*r
            end if
!
!         CALCUL EFFECTIF DU SECOND MEMBRE
!         --------------------------------
            if (pre1) then
                pos = 0
                do ino = 1, nnop
!
!           TERME CLASSIQUE
                    do j = 1, ndim
                        pos = pos+1
                        zr(ires-1+pos) = zr(ires-1+pos)+forrep(j)*poids*ff(ino)
                    end do
!
!           ON ZAPPE LES TERMES DE PRESSION CLASSIQUE SI ON ES SUR UN NOEUD SOMMET
                    if (ino .le. nnops) pos = pos+1
!           TERME HEAVISIDE
                    do ig = 1, nfh
                        do j = 1, ndim
                            pos = pos+1
                            zr(ires-1+pos) = zr(ires-1+pos)+xcalc_heav( &
                                             zi(jheavn-1+ncompn*(ino-1)+ig) &
                                             , zi(jheavs-1+ise), &
                                             zi(jheavn-1+ncompn*(ino-1)+ncompn)) &
                                             *forrep(j)*poids*ff(ino)
                        end do
!           ON ZAPPE LES TERMES DE PRESSION HEAVISIDE SI ON ES SUR UN NOEUD SOMMET
                        if (ino .le. nnops) then
                            pos = pos+1
                        end if
                    end do
                end do
            else
                pos = 0
                do ino = 1, nnop
!
!           TERME CLASSIQUE
                    do j = 1, ndim
                        pos = pos+1
                        zr(ires-1+pos) = zr(ires-1+pos)+forrep(j)*poids*ff(ino)
                    end do
!
!           TERME HEAVISIDE
                    do ig = 1, nfh
                        do j = 1, ndim
                            pos = pos+1
                            zr(ires-1+pos) = zr(ires-1+pos)+xcalc_heav( &
                                             zi(jheavn-1+ncompn*(ino-1)+ig) &
                                             , zi(jheavs-1+ise), &
                                             zi(jheavn-1+ncompn*(ino-1)+ncompn)) &
                                             *forrep(j)*poids*ff(ino)
                        end do
                    end do
!           TERME SINGULIER
                    do alp = 1, ndim*nfe
                        pos = pos+1
                        do j = 1, ndim
                            zr(ires-1+pos) = zr(ires-1+pos)+fk(ino, alp, j)*forrep(j)*poids
                        end do
                    end do
                end do
            end if
        end do
!
!-----------------------------------------------------------------------
!         FIN DE LA BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELT
!-----------------------------------------------------------------------
!
    end do
!
!     SUPPRESSION DES DDLS SUPERFLUS
    if (pre1) then
        ddls = ndim*(1+nfh)+(1+nfh)
        ddlm = ndim*(1+nfh)
        nnopm = nnop-nnops
        nddl = nnops*ddls+nnopm*ddlm
        contac = 0
        call xhmddl(ndim, nfh, ddls, nddl, nnop, nnops, zi(jstno), .false._1, &
                    option, nomte, rbid, zr(ires), ddlm, nfiss, jfisno, .false._1, contac)
    else
        ddls = ndim*(1+nfh+nfe)
        nddl = nnop*ddls
        call teattr('C', 'XLAG', lag, ibid)
        if (ibid .eq. 0 .and. lag .eq. 'ARETE') then
            nnop = nnos
        end if
        call xteddl(ndim, nfh, nfe, ddls, nddl, &
                    nnop, nnops, zi(jstno), .false._1, lbid, &
                    option, nomte, ddls, &
                    nfiss, jfisno, vect=zr(ires))
    end if
!
!-----------------------------------------------------------------------
!     FIN
!-----------------------------------------------------------------------
!
end subroutine
