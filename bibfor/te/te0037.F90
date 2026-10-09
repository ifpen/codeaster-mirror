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
subroutine te0037(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/iselli.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xhmddl.h"
#include "asterfort/xhmini.h"
#include "asterfort/xjacf2.h"
#include "asterfort/xjacff.h"
#include "asterfort/xkamat.h"
#include "asterfort/xteddl.h"
#include "asterfort/xteini.h"
#include "asterfort/xxmmvd.h"
#include "MeshTypes_type.h"
!
    character(len=16) :: option, nomte
!
!
!.......................................................................
!
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN MECANIQUE
!          CORRESPONDANT A UN CHARGEMENT EN PRESSION REPARTIE
!          SUR LES LEVRES DES FISSURES X-FEM
!          (LA PRESSION PEUT ETRE DONNEE SOUS FORME D'UNE FONCTION)
!
!          OPTIONS : 'CHAR_MECA_PRES_R'
!                    'CHAR_MECA_PRES_F'
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!.......................................................................
!
!
    character(len=8) :: elref, typma, fpg, elc, nompar(4), lag, elrefc, enr, enr2
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfde, jgano, fisno(27, 10)
    integer(kind=8) :: nfh, nfe, singu, ddlc, nnom, ddls, nddl, ier, ddlm
    integer(kind=8) :: igeom, ipres, itemps, ires, iadzi, iazk24, jheavn, ncompn, hea_fa(2)
    integer(kind=8) :: jlst, jptint, jaint, jcface, jlonch, jstno, jbasec, contac
    integer(kind=8) :: i, j, ninter, nface, cface(30, 6), ifa, nfiss, jfisno
    integer(kind=8) :: ibid, ilev, ifiss, ncompc, jtab(7), ncompp, ino
    integer(kind=8) :: imate, jlsn, jbaslo
    integer(kind=8) :: alp
    integer(kind=8) :: nnof, npgf, ipoidf, ivff, idfdef, ipgf, pos, zxain, nptf, ifh
    real(kind=8) :: pres, cisa, forrep(3, 2), ff(MT_NNOMAX3D), jac, nd(3), he(2), mat(1)
    real(kind=8) :: xg(4), dfbid(MT_NNOMAX3D, 3), r27bid(MT_NNOMAX3D), r3bid(3), r
    aster_logical :: pre1, axi
    integer(kind=8) :: compt, nddlm, nddls, nddlp, iret, jheafa, ncomph, ncompb
    real(kind=8) :: thet
    real(kind=8) :: fk(MT_NNOMAX3D, 3, 3), ka, mu
    data he/-1.d0, 1.d0/
!
    call jemarq()
!
!     PAR CONVENTION :
!     LEVRE INFERIEURE (HE=-1) EST LA LEVRE 1, DE NORMALE SORTANTE  ND
!     LEVRE SUPERIEURE (HE=+1) EST LA LEVRE 2, DE NORMALE SORTANTE -ND
!
!-----------------------------------------------------------------------
!     INITIALISATIONS
!-----------------------------------------------------------------------
    zxain = xxmmvd('ZXAIN')
!
    call elref1(elref)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    axi = lteatt('AXIS', 'OUI')
!
    call teattr('C', 'HYDR1', enr2, iret)
    pre1 = (enr2 .eq. '1' .or. enr2 .eq. '2')
!
!-----------------------------------------------------------------------
!     RECUPERATION DES ENTREES / SORTIE
!-----------------------------------------------------------------------
!
    if (option .eq. 'CHAR_MECA_PRES_R') then
!
!       SI LA PRESSION N'EST CONNUE SUR AUCUN NOEUD, ON LA PREND=0.
        call jevecd('PPRESSR', ipres, 0.d0)
        compt = 0
        do i = 1, nno
            thet = abs(zr(ipres-1+(i-1)+1))
            if (thet .lt. r8prem()) compt = compt+1
        end do
        if (compt .eq. nno) goto 999
!
    else if (option .eq. 'CHAR_MECA_PRES_F') then
!
        call jevech('PPRESSF', 'L', ipres)
        call jevech('PINSTR', 'L', itemps)
!
    else
        ASSERT(.false.)
    end if
!
    call jevech('PVECTUR', 'E', ires)
!
!   INITIALISATION DES DIMENSIONS DES DDLS X-FEM
!     SI PRE1=.FALSE. -> MODELISATION MECA XFEM CLASSIQUE
!     SI PRE1=.TRUE.  -> MODELISATION HM XFEM
    if (pre1) then
        call xhmini(nomte, nfh, ddls, ddlm, nddlp, &
                    nfiss, ddlc, contac)
!
        nfe = 0
        nddls = ddls+nddlp+ddlc
        nddlm = ddlm
        nnom = nno-nnos
        nddl = nnos*nddls+nnom*nddlm
    else
        call xteini(nomte, nfh, nfe, singu, ddlc, &
                    nnom, ddls, nddl, ddlm, nfiss, &
                    contac)
    end if
!
    if (nfe .gt. 0) then
        call jevech('PMATERC', 'L', imate)
        call jevech('PBASLOR', 'L', jbaslo)
        call jevech('PLSN', 'L', jlsn)
        call xkamat(zi(imate), ndim, axi, ka, mu, &
                    famiz='XCON')
    end if
!
    call tecael(iadzi, iazk24, noms=0)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3) (1:8)
!
    if (ndim .eq. 3) then
        if (iselli(elref)) then
            elc = 'TR3'
        else
            elc = 'TR3'
        end if
        fpg = 'XCON'
    else if (ndim .eq. 2) then
        if (iselli(elref)) then
            elc = 'SE2'
        else
            elc = 'SE3'
        end if
        fpg = 'MASS'
    end if
!
!     PARAMETRES PROPRES A X-FEM
    call jevech('PLST', 'L', jlst)
    call jevech('PPINTER', 'L', jptint)
    call tecach('OOO', 'PPINTER', 'L', iret, nval=2, &
                itab=jtab)
    ncompp = jtab(2)
    call jevech('PAINTER', 'L', jaint)
    call jevech('PCFACE', 'L', jcface)
    call tecach('OOO', 'PCFACE', 'L', iret, nval=2, &
                itab=jtab)
    ncompc = jtab(2)
    call tecach('OOO', 'PBASECO', 'L', iret, nval=2, &
                itab=jtab)
    ncompb = jtab(2)
    call jevech('PLONGCO', 'L', jlonch)
    call jevech('PSTANO', 'L', jstno)
    call jevech('PBASECO', 'L', jbasec)
    if (nfiss .gt. 1) then
        call jevech('PFISNO', 'L', jfisno)
        call jevech('PHEA_FA', 'L', jheafa)
        call tecach('OOO', 'PHEA_FA', 'L', iret, nval=2, &
                    itab=jtab)
        ncomph = jtab(2)
    end if
!
! --- CONECTIVITÃ~I DES FISSURE ET DES DDL HEAVISIDES
    if (nfiss .eq. 1) then
        do ino = 1, nno
            fisno(ino, 1) = 1
        end do
    else
        do ifh = 1, nfh
            do ino = 1, nno
                fisno(ino, ifh) = zi(jfisno-1+(ino-1)*nfh+ifh)
            end do
        end do
    end if
!
    call jevech('PGEOMER', 'L', igeom)
!
!     RECUPERATION DE LA DEFNITION DES FONCTIONS HEAVISIDE
    hea_fa(1:2) = 0
    call teattr('S', 'XFEM', enr, ier)
    call jevech('PHEA_NO', 'L', jheavn)
    call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                itab=jtab)
    ncompn = jtab(2)/jtab(3)
    if (enr(1:2) .eq. 'XH' .and. nfiss .eq. 1) then
        hea_fa(1) = xcalc_code(1, he_real=[he(1)])
        hea_fa(2) = xcalc_code(1, he_real=[he(2)])
    end if
!
!     RÉCUPÉRATIONS DES DONNÉES SUR LA TOPOLOGIE DES FACETTES
    do ifiss = 1, nfiss
        ninter = zi(jlonch+3*(ifiss-1)-1+1)
        nface = zi(jlonch+3*(ifiss-1)-1+2)
        nptf = zi(jlonch+3*(ifiss-1)-1+3)
        if (nptf .eq. 6 .and. ndim .eq. 3) elc = 'TR6'
        if (ninter .lt. ndim) goto 998
!
        do i = 1, nface
            do j = 1, nptf
                cface(i, j) = zi(jcface-1+ncompc*(ifiss-1)+nptf*(i-1)+j)
            end do
        end do
!
!-----------------------------------------------------------------------
!     BOUCLE SUR LES FACETTES
!-----------------------------------------------------------------------
!
        do ifa = 1, nface
!
            call elrefe_info(elrefe=elc, fami=fpg, nno=nnof, npg=npgf, jpoids=ipoidf, &
                             jvf=ivff, jdfde=idfdef)
!
            ff(:) = 0.d0
!
!       BOUCLE SUR LES POINTS DE GAUSS DES FACETTES
            do ipgf = 1, npgf
!
!         CALCUL DE JAC (PRODUIT DU JACOBIEN ET DU POIDS)
!         ET DES FF DE L'ÉLÉMENT PARENT AU POINT DE GAUSS
!         ET LA NORMALE ND ORIENTÉE DE ESCL -> MAIT
!         ET DE XG : COORDONNEES REELLES DU POINT DE GAUSS
                elrefc = 'NON'
                if (ndim .eq. 3) then
                    call xjacff(elref, elrefc, elc, ndim, fpg, &
                                jptint, ifa, cface, ipgf, nno, &
                                nnos, igeom, jbasec, xg, jac, &
                                ff, r27bid, dfbid, nd, r3bid, &
                                r3bid)
                else if (ndim .eq. 2) then
                    call xjacf2(elref, elrefc, elc, ndim, fpg, &
                                jptint, ifa, cface, nptf, ipgf, &
                                nno, nnos, igeom, jbasec, xg, &
                                jac, ff, r27bid, dfbid, nd, &
                                r3bid)
                end if
!
!         CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE)
                if (axi) then
                    r = 0.d0
                    do ino = 1, nno
                        r = r+ff(ino)*zr(igeom-1+2*(ino-1)+1)
                    end do
                    ASSERT(r .ge. 0d0)
!              ATTENTION : LE POIDS N'EST PAS X R
!              CE SERA FAIT PLUS TARD AVEC JAC = JAC X R
                end if
!
!         CALCUL DES FORCES REPARTIES SUIVANT LES OPTIONS
!         -----------------------------------------------
!
                forrep(:, :) = 0.d0
                nompar(1) = 'X'
                nompar(2) = 'Y'
                if (ndim .eq. 3) nompar(3) = 'Z'
                if (ndim .eq. 3) nompar(4) = 'INST'
                if (ndim .eq. 2) nompar(3) = 'INST'
!
! MODIFIER LE JAC
                if (axi) then
                    jac = jac*r
                end if
!
                if (option .eq. 'CHAR_MECA_PRES_R') then
!
!           CALCUL DE LA PRESSION AUX POINTS DE GAUSS
                    pres = 0.d0
                    cisa = 0.d0
                    do ino = 1, nno
                        if (ndim .eq. 3) pres = pres+zr(ipres-1+ino)*ff(ino)
                        if (ndim .eq. 2) then
                            pres = pres+zr(ipres-1+2*(ino-1)+1)*ff(ino)
                            cisa = cisa+zr(ipres-1+2*(ino-1)+2)*ff(ino)
                        end if
                    end do
!           ATTENTION AU SIGNE : POUR LES PRESSIONS, IL FAUT UN - DVT
!           CAR LE SECOND MEMBRE SERA ECRIT AVEC UN + (VOIR PLUS BAS)
!           ON CALCULE FORREP POUR LES DEUX LEVRES  : 1 = INF ET 2 = SUP
                    do j = 1, ndim
                        forrep(j, 1) = -pres*nd(j)
                        forrep(j, 2) = -pres*(-nd(j))
                    end do
                    if (ndim .eq. 2) then
                        forrep(1, 1) = forrep(1, 1)-cisa*nd(2)
                        forrep(2, 1) = forrep(2, 1)+cisa*nd(1)
                        forrep(1, 2) = forrep(1, 2)-cisa*(-nd(2))
                        forrep(2, 2) = forrep(2, 2)+cisa*(-nd(1))
                    end if
!
                else if (option .eq. 'CHAR_MECA_PRES_F') then
!
!           VALEUR DE LA PRESSION
                    xg(ndim+1) = zr(itemps)
                    call fointe('FM', zk8(ipres), ndim+1, nompar, xg, &
                                pres, ier)
                    if (ndim .eq. 2) call fointe('FM', zk8(ipres+1), ndim+1, nompar, xg, &
                                                 cisa, ier)
                    do j = 1, ndim
                        forrep(j, 1) = -pres*nd(j)
                        forrep(j, 2) = -pres*(-nd(j))
                    end do
                    if (ndim .eq. 2) then
                        forrep(1, 1) = forrep(1, 1)-cisa*nd(2)
                        forrep(2, 1) = forrep(2, 1)+cisa*nd(1)
                        forrep(1, 2) = forrep(1, 2)-cisa*(-nd(2))
                        forrep(2, 2) = forrep(2, 2)+cisa*(-nd(1))
                    end if
                else
                    call utmess('F', 'XFEM_15')
                end if
!
!         CALCUL EFFECTIF DU SECOND MEMBRE SUR LES DEUX LEVRES
                if (pre1) then
                    do ilev = 1, 2
                        pos = 0
                        do ino = 1, nno
!
!               TERME CLASSIQUE
                            do j = 1, ndim
                                pos = pos+1
                                zr(ires-1+pos) = zr(ires-1+pos)+forrep(j, ilev)*jac*ff(ino)
                            end do
!
!               ON ZAPPE LES TERMES DE PRESSION CLASSIQUE SI ON EST SUR UN
!               NOEUD SOMMET
                            if (ino .le. nnos) pos = pos+1
!
!               TERME HEAVISIDE
                            do ifh = 1, nfh
!               EN MULTI-FISSURATION, IL FAUT RECUPERER LES BONNES VALEURS DE HE
                                if (nfiss .gt. 1) then
                                    hea_fa(ilev) = zi(jheafa-1+ncomph*(ifiss-1)+2*(ifa-1)+ilev)
                                end if
                                do j = 1, ndim
                                    pos = pos+1
                                    zr(ires-1+pos) = zr(ires-1+pos)+xcalc_heav(zi(jheavn-1+nco&
                                                     &mpn*(ino-1)+ifh), hea_fa(ilev), zi(jheavn-1+&
                                                     &ncompn*(ino-1)+ncompn))*forrep(j, ilev)*jac&
                                                     &*ff(ino)
                                end do
!               ON ZAPPE LES TERMES DE PRESSION HEAVISIDE SI ON
!               EST SUR UN NOEUD SOMMET
                                if (ino .le. nnos) pos = pos+1
                            end do
!               ON ZAPPE LES TERMES DE CONTACT DU MODELE HM-XFEM
                            if (contac .ge. 2) then
                                if (ino .le. nnos) pos = pos+ddlc
                            end if
                        end do
                    end do
                else
                    do ilev = 1, 2
!
                        pos = 0
                        if (nfe .gt. 0) then
                            if (he(ilev) .gt. 0) then
                                call xcalfev_wrap(ndim, nno, zr(jbaslo), zi(jstno), he(ilev), &
                                                  zr(jlsn), zr(jlst), zr(igeom), ka, mu, &
                                                  ff, fk, face='MAIT')
                            else
                                call xcalfev_wrap(ndim, nno, zr(jbaslo), zi(jstno), he(ilev), &
                                                  zr(jlsn), zr(jlst), zr(igeom), ka, mu, &
                                                  ff, fk, face='ESCL')
                            end if
                        end if
                        do ino = 1, nno
!
!               TERME CLASSIQUE
                            do j = 1, ndim
                                pos = pos+1
                                zr(ires-1+pos) = zr(ires-1+pos)+forrep(j, ilev)*jac*ff(ino)
                            end do
!
!               TERME HEAVISIDE
                            do j = 1, ndim
                                do ifh = 1, nfh
                                    pos = pos+1
                                    zr(ires-1+pos) = zr(ires-1+pos)+xcalc_heav(zi(jheavn-1+nco&
                                                     &mpn*(ino-1)+ifh), hea_fa(ilev), zi(jheavn-1&
                                                     &+ncompn*(ino-1)+ncompn))*forrep(j, ilev)*&
                                                     & jac*ff(ino)
                                end do
                            end do
!
!               TERME SINGULIER
                            do alp = 1, nfe*ndim
                                pos = pos+1
!               PAS DE CISAILLEMENT MODE III
                                if (alp .eq. 3) goto 555
                                do j = 1, ndim
                                    zr(ires-1+pos) = zr(ires-1+pos)+fk(ino, alp, j)*forrep(j, il&
                                                     &ev)*jac
                                end do
555                             continue
                            end do
!
!               ON SAUTE LES POSITIONS DES LAG DE CONTACT FROTTEMENT
!
                            if (contac .eq. 3) then
                                if (ino .le. nnos) pos = pos+ddlc
                            else
                                pos = pos+ddlc
                            end if
!
                        end do
                    end do
                end if
            end do
        end do
998     continue
        jbasec = jbasec+ncompb
        jptint = jptint+ncompp
    end do
!
!     SUPPRESSION DES DDLS SUPERFLUS
    if (pre1) then
        call xhmddl(ndim, nfh, nddls, nddl, nno, nnos, zi(jstno), &
                    ASTER_FALSE, option, nomte, &
                    mat, zr(ires), nddlm, nfiss, jfisno, ASTER_FALSE, contac)
    else
        call teattr('C', 'XLAG', lag, ibid)
        if (ibid .eq. 0 .and. lag .eq. 'ARETE') then
            nno = nnos
        end if
        call xteddl(ndim, nfh, nfe, ddls, nddl, &
                    nno, nnos, zi(jstno), ASTER_FALSE, ASTER_FALSE, &
                    option, nomte, ddlm, nfiss, jfisno, &
                    vect=zr(ires))
    end if
!
!
999 continue
    call jedema()
end subroutine
