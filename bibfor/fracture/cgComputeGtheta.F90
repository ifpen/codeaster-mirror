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
!
subroutine cgComputeGtheta(cgField, cgTheta, cgStudy, cgTable, cgStat)
!
    use calcG_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/cescre.h"
#include "asterfort/cgDiscrField.h"
#include "asterfort/cgHighOrder.h"
#include "asterfort/cgTempNodes.h"
#include "asterfort/chpchd.h"
#include "asterfort/chpver.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcharg.h"
#include "asterfort/getvid.h"
#include "asterfort/glegen.h"
#include "asterfort/gsyste.h"
#include "asterfort/hatSmooth.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/mesomm.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterfort/wkvect.h"
#include "asterfort/exicp.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/jexnum.h"
#include "asterfort/typele.h"
#include "asterfort/nbelem.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "jeveux.h"
!
    type(CalcG_field), intent(in)    :: cgField
    type(CalcG_theta), intent(inout) :: cgTheta
    type(CalcG_Study), intent(inout) :: cgStudy
    type(CalcG_Table), intent(inout) :: cgTable
    type(CalcG_stat), intent(inout) :: cgStat
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Compute G(Theta) in 2D and 3D
!
!----------------------------------------------
    integer(kind=8) :: iret, nsig, ino1, ino2, inga, ibid, i_theta, i, j
    integer(kind=8) :: nchin, iadrt3
    integer(kind=8) :: jcesd, jcesl, jcesd2, jcesl2
    integer(kind=8) :: iad, iel, igr, ima, nbgrel, nel, nncp, nute
    real(kind=8) :: gth(7), som(7)
    real(kind=8) :: s1, s2, s3, sn2, sn1, sn
    character(len=2)  :: codret
    character(len=8)  :: k8b, lpain(50), lpaout(1), model
    character(len=16) :: opti, nomte
    character(len=19) :: chrota, chpesa, cf2d3d, chpres, chvolu, cf1d2d, chepsi, ligrmo, chscpl
    character(len=19) :: chvarc, chvref, chsdeg, modelLigrel, chslag, chscer, chseli, chcer, cheli
    character(len=24) :: chsigi, celmod, sigelno, chtime, chpuls
    character(len=24) :: chgeom, chsig, chgtheta
    character(len=24) :: pavolu, papres, pa2d3d, pepsin, pa1d2d, iscplan
    character(len=24) :: lchin(50), lchout(1)
    aster_logical     :: lfonc, inco
    integer(kind=8), pointer :: v_cesv(:) => null()
    integer(kind=8), pointer :: v_liel(:) => null()
    real(kind=8), pointer :: v_absc(:) => null()
    real(kind=8), pointer :: v_base(:) => null()
    real(kind=8), pointer :: v_basf(:) => null()
    real(kind=8), pointer :: v_cer(:) => null()
    real(kind=8), pointer :: v_eli(:) => null()
    real(kind=8), dimension(cgTheta%nb_theta_field) :: gthi, k1th, k2th, k3th, g1th, g2th, g3th
    real(kind=8), dimension(cgTheta%nnof) :: gs, k1s, k2s, k3s, g1s, g2s, g3s, gis
    real(kind=8) :: finish, start, start0, finish0
    aster_logical :: outExicp

!----------------------------------------------
!
    call cpu_time(start)
!
    call jemarq()
!
!   Initialisation des champs et des paramètres
    chgtheta = '&&cgtheta.CH_G'
    chvarc = '&&cgtheta.VARC'
    chvref = '&&cgtheta.VARC.REF'
    chsigi = '&&cgtheta.CHSIGI'
    celmod = '&&cgtheta.CELMOD'
    sigelno = '&&cgtheta.SIGELNO'
    chvolu = '&&cgtheta.VOLU'
    cf1d2d = '&&cgtheta.1D2D'
    cf2d3d = '&&cgtheta.2D3D'
    chpres = '&&cgtheta.PRES'
    chepsi = '&&cgtheta.EPSI'
    chpesa = '&&cgtheta.PESA'
    chrota = '&&cgtheta.ROTA'
    chtime = '&&cgtheta.CH_INST_R'
    chpuls = '&&cgtheta.PULS'
    chsdeg = '&&cgtheta.CHSDEG'
    chslag = '&&cgtheta.CHSLAG'
    chscer = '&&cgtheta.CHSCER'
    chseli = '&&cgtheta.CHSELI'
    chcer = '&&cgtheta.CHCER'
    cheli = '&&cgtheta.CHELI'
    chscpl = '&&cgtheta.CHSCPL'
    iscplan = '&&cgtheta.ISCPLAN'
    !
    gthi(:) = 0.0
    k1th(:) = 0.0
    k2th(:) = 0.0
    k3th(:) = 0.0
    g1th(:) = 0.0
    g2th(:) = 0.0
    g3th(:) = 0.0
    gs(:) = 0.0
    k1s(:) = 0.0
    k2s(:) = 0.0
    k3s(:) = 0.0
    g1s(:) = 0.0
    g2s(:) = 0.0
    g3s(:) = 0.0
    cgStudy%gth = 0.d0
!
    k8b = '        '
!   Recuperation du champ geometrique
    chgeom = cgStudy%mesh//'.COORDO'
!
!   Recuperation du LIGREL
    model = cgStudy%model(1:8)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
!
!   Recuperation du comportement
    if (cgField%l_incr) then
!
        if (cgField%result_in_type .ne. 'EVOL_NOLI') then
            call utmess('F', 'RUPTURE1_15')
        end if
!
        call rsexch('F', cgField%result_in, 'SIEF_ELGA', cgStudy%nume_ordre, chsig, iret)
    end if
!
!   Cas 2D : creation du champ ISCPLAN
    if (cgField%ndim .eq. 2) then
!       Recuperation de la modelisation
        outExicp = exicp(cgStudy%model, ASTER_TRUE, "", 0)
!
!       Creation d'un champ simple par elements
        call cescre('V', chscpl, 'ELEM', cgStudy%mesh, 'NEUT_I', &
                    0, ' ', [-1], [-1], [-1])

        call jeveuo(chscpl//'.CESD', 'L', jcesd)
        call jeveuo(chscpl//'.CESL', 'E', jcesl)
        call jeveuo(chscpl//'.CESV', 'E', vi=v_cesv)

        call dismoi('NOM_LIGREL', cgStudy%model, 'MODELE', repk=ligrmo)
        call jelira(ligrmo//'.LIEL', 'NMAXOC', nbgrel)
!
        do igr = 1, nbgrel
!
!           Récupération nombre d'éléments dans le groupe
            nel = nbelem(ligrmo, igr)
            call jeveuo(jexnum(ligrmo//'.LIEL', igr), 'L', vi=v_liel)
!
!           Boucle sur les éléments du groupe
            do iel = 1, nel
                ima = v_liel(iel)
                if (ima .lt. 0) cycle
!
                call cesexi('C', jcesd, jcesl, ima, 1, 1, 1, iad)
                iad = abs(iad)
                zl(jcesl-1+iad) = ASTER_TRUE
                if (outExicp) then
                    v_cesv(iad) = 1
                else
                    v_cesv(iad) = 0
                end if

            end do
        end do
!
!       Conversion du champ par élément simple en champ par éléments
        call cescel(chscpl, ligrmo, 'CALC_G', 'ISCPLAN', 'NON', &
                    nncp, 'V', iscplan, 'F', iret)
!
    end if
!   Elements incompressibles
    inco = cgStudy%l_exi_inco
!
!   Recuperation de l'etat initial
    if (cgField%l_incr) then
!
        call getvid('ETAT_INIT', 'SIGM', iocc=1, scal=chsigi, nbret=nsig)
!
!       Verification du type de champ + transfo, si necessaire en champ elno
        if (nsig .ne. 0) then
!
!           chpver renvoit 0 si OK et 1 si PB
            call chpver('C', chsigi(1:19), 'ELNO', 'SIEF_R', ino1)
            call chpver('C', chsigi(1:19), 'NOEU', 'SIEF_R', ino2)
            call chpver('C', chsigi(1:19), 'ELGA', 'SIEF_R', inga)

!           Verification du type de champ
            if (ino1 .eq. 1 .and. ino2 .eq. 1 .and. inga .eq. 1) then
                call utmess('F', 'RUPTURE1_12')
            end if

!           Transformation si champ ELGA
            if (inga .eq. 0) then

!               Traitement du champ pour les elements finis classiques
                call detrsd('CHAMP', celmod)
                call alchml(modelLigrel, 'CALC_G', 'PSIGINR', 'V', celmod, iret, ' ')
                call chpchd(chsigi(1:19), 'ELNO', celmod, 'OUI', 'V', sigelno, model)
                call chpver('F', sigelno(1:19), 'ELNO', 'SIEF_R', ibid)
            end if
        end if
    else
        nsig = 0
    end if
!
!   Recuperation des champs de temperature (T,TREF)
    call vrcins(cgStudy%model, cgStudy%material, k8b, cgStudy%time, chvarc, codret)
    call vrcref(cgStudy%model, cgStudy%material(1:8), k8b, chvref(1:19))
!
!   Traitement des charges
    call gcharg(cgStudy%model, cgStudy%loading, chvolu, cf1d2d, cf2d3d, chpres, &
                chepsi, chpesa, chrota, lfonc, cgStudy%time, cgStudy%nume_ordre)
!
!   Select name of option
    if (lfonc) then
        if (cgField%ndim .eq. 2) then
            pavolu = 'PFFVOLU'
            pa1d2d = 'PFF1D2D'
            papres = 'PPRESSF'
            pepsin = 'PEPSINF'
        else
            pavolu = 'PFFVOLU'
            pa2d3d = 'PFF2D3D'
            papres = 'PPRESSF'
            pepsin = 'PEPSINF'
        end if
!
        if (cgStudy%option .eq. 'K') then
            opti = 'CALC_K_G_F'
        else if (cgStudy%option .eq. 'G' .or. cgStudy%option .eq. 'G_EPSI') then
            opti = 'CALC_G_F'
        else if (cgStudy%option .eq. 'KJ' .or. cgStudy%option .eq. 'KJ_EPSI') then
            opti = 'CALC_G_F'
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        if (cgField%ndim .eq. 2) then
            pavolu = 'PFRVOLU'
            pa1d2d = 'PFR1D2D'
            papres = 'PPRESSR'
            pepsin = 'PEPSINR'
        else
            pavolu = 'PFRVOLU'
            pa2d3d = 'PFR2D3D'
            papres = 'PPRESSR'
            pepsin = 'PEPSINR'
        end if
!
        if (cgStudy%option .eq. 'K') then
            opti = 'CALC_K_G'
        else if (cgStudy%option .eq. 'G' .or. cgStudy%option .eq. 'G_EPSI') then
            opti = 'CALC_G'
        else if (cgStudy%option .eq. 'KJ' .or. cgStudy%option .eq. 'KJ_EPSI') then
            opti = 'CALC_G'
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
!   Création d'un champ constant par élément : LEGENDRE, LAGRANGE

    if (cgtheta%discretization .eq. 'LEGENDRE') then
!
!       Champ scalaire constant par élement :
!            --  degré du polynome
        call cescre('V', chsdeg, 'ELEM', cgStudy%mesh, 'NEUT_I', &
                    0, ' ', [-1], [-1], [-1])
!
        call jeveuo(chsdeg//'.CESD', 'L', jcesd)
        call jeveuo(chsdeg//'.CESL', 'E', jcesl)
        call jeveuo(chsdeg//'.CESV', 'E', vi=v_cesv)
!
    elseif (cgtheta%discretization .eq. 'LINEAIRE') then
!
        call cgTheta%getAbsfon(v_basf)
!
!       Champ scalaire constant par élement :
!            -- abscisse curviligne s(i-1), s(i), s(i+1)
!               pour la fonction de forme courante
        call cescre('V', chslag, 'ELEM', cgStudy%mesh, 'NEUT_R', &
                    0, '', [-1], [-1], [-3])
!
        call jeveuo(chslag//'.CESD', 'L', jcesd)
        call jeveuo(chslag//'.CESL', 'E', jcesl)
        call jeveuo(chslag//'.CESV', 'E', vr=v_absc)
!
    end if
!
    if (cgTheta%form_fiss .eq. 'CERCLE') then
!
!       Champ scalaire constant par élement :
!            --  rayon du cercle
        call cescre('V', chscer, 'ELEM', cgStudy%mesh, 'NEUT_R', &
                    0, ' ', [-1], [-1], [-1])
!
        call jeveuo(chscer//'.CESD', 'L', jcesd2)
        call jeveuo(chscer//'.CESL', 'E', jcesl2)
        call jeveuo(chscer//'.CESV', 'E', vr=v_cer)

    elseif (cgTheta%form_fiss .eq. 'ELLIPSE') then

!       Champ scalaire constant par élement :
!            --  demi grand axe et demi petit axe
        call cescre('V', chseli, 'ELEM', cgStudy%mesh, 'NEUT_R', &
                    0, ' ', [-1], [-1], [-2])
!
        call jeveuo(chseli//'.CESD', 'L', jcesd2)
        call jeveuo(chseli//'.CESL', 'E', jcesl2)
        call jeveuo(chseli//'.CESV', 'E', vr=v_eli)

    end if

    if (cgTheta%form_fiss .eq. 'CERCLE' .or. cgTheta%form_fiss .eq. 'ELLIPSE') then

        call cgHighOrder(cgField, cgTheta, cgStudy, cgStat, chscer, chseli, v_cer, v_eli, &
                         jcesd2, jcesl2, chcer, cheli)
    end if
!
!************************************************
!   Boucle sur le nombre de champ theta
!    --- 1 si 2D
!    --- ndeg+1 ou nnof si 3D
!************************************************
    do i_theta = 1, cgTheta%nb_theta_field
!
!       Declaration des entrees/sorties pour l'appel à calcul (indépendent de i_theta)
        lpaout(1) = 'PGTHETA'
        lchout(1) = chgtheta
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PDEPLAR'
        lchin(2) = cgStudy%depl
        lpain(3) = 'PTHETAR'
        lchin(3) = cgTheta%theta_factors
        lpain(4) = 'PMATERC'
        lchin(4) = cgStudy%mateco
        lpain(5) = 'PVARCPR'
        lchin(5) = chvarc
        lpain(6) = 'PVARCRR'
        lchin(6) = chvref
        lpain(7) = pavolu(1:8)
        lchin(7) = chvolu
        if (cgField%ndim .eq. 2) then
            lpain(8) = pa1d2d(1:8)
            lchin(8) = cf1d2d
        else
            lpain(8) = pa2d3d(1:8)
            lchin(8) = cf2d3d
        end if
        lpain(9) = papres(1:8)
        lchin(9) = chpres
        lpain(10) = 'PPESANR'
        lchin(10) = chpesa
        lpain(11) = 'PROTATR'
        lchin(11) = chrota
        lpain(12) = pepsin(1:8)
        lchin(12) = chepsi
        lpain(13) = 'PCOMPOR'
        lchin(13) = cgField%compor
        lpain(14) = 'PBASLOR'
        lchin(14) = cgTheta%crack//'.BASLOC'
!
        nchin = 14

        if (cgField%ndim .eq. 2) then
            lpain(nchin+1) = 'ISCPLAN'
            lchin(nchin+1) = iscplan
            nchin = nchin+1
        end if
!
        call cgDiscrField(cgField, cgTheta, cgStudy, cgStat, chsdeg, chslag, &
                          v_absc, v_basf, v_cesv, jcesd, jcesl, i_theta, lpain, lchin, nchin)

        if (cgTheta%form_fiss .eq. 'CERCLE') then
            lpain(nchin+1) = 'PCER'
            lchin(nchin+1) = chcer
            nchin = nchin+1
        elseif (cgTheta%form_fiss .eq. 'ELLIPSE') then
            lpain(nchin+1) = 'PELI'
            lchin(nchin+1) = cheli
            nchin = nchin+1
        end if

        if (inco) then
            lpain(nchin+1) = 'PDEPLA'
            lchin(nchin+1) = cgStudy%depl
            nchin = nchin+1
        end if

        if (cgStudy%option .eq. 'K') then
            lpain(nchin+1) = 'PCOURB'
            lchin(nchin+1) = cgTheta%curvature
            lpain(nchin+2) = 'PLSN'
            lchin(nchin+2) = cgTheta%crack//'.LNNO'
            lpain(nchin+3) = 'PLST'
            lchin(nchin+3) = cgTheta%crack//'.LTNO'
            nchin = nchin+3
        end if
!
        if (opti .eq. 'CALC_G_F' .or. opti .eq. 'CALC_K_G_F') then
            call mecact('V', chtime, 'MODELE', modelLigrel, 'INST_R', &
                        ncmp=1, nomcmp='INST', sr=cgStudy%time)
            lpain(nchin+1) = 'PINSTR'
            lchin(nchin+1) = chtime
            nchin = nchin+1
        end if
!
        if (cgStudy%l_modal) then
            call mecact('V', chpuls, 'MODELE', modelLigrel, 'FREQ_R  ', &
                        ncmp=1, nomcmp='FREQ   ', sr=cgStudy%pulse)
            nchin = nchin+1
            lpain(nchin) = 'PPULPRO'
            lchin(nchin) = chpuls
        end if
!
        if (cgField%l_incr) then
            if (opti .eq. 'CALC_G_F' .or. opti .eq. 'CALC_G') then
                lpain(nchin+1) = 'PCONTRR'
                lchin(nchin+1) = chsig
                nchin = nchin+1
            end if
        end if
!
!       Champ de contrainte initiale
        if (nsig .ne. 0) then
            lpain(nchin+1) = 'PSIGINR'
            if (inga .eq. 0) then
!           Champ de contrainte initiale transforme en ELNO
                lchin(nchin+1) = sigelno
            else
!           Champ de contrainte initiale donne par l'utilisateur (NOEUD ou ELNO)
                lchin(nchin+1) = chsigi
            end if
            nchin = nchin+1
        end if
!
        if (cgStudy%option .eq. 'G' .or. cgStudy%option .eq. 'G_EPSI' &
            .or. cgStudy%option .eq. 'KJ' .or. cgStudy%option .eq. 'KJ_EPSI') then
            if (cgStudy%vitesse .ne. ' ') then
                lpain(nchin+1) = 'PVITESS'
                lchin(nchin+1) = cgStudy%vitesse
                lpain(nchin+2) = 'PACCELE'
                lchin(nchin+2) = cgStudy%acce
                nchin = nchin+2
            end if
        end if
!
!       Recuperation des contraintes du resultat pour option G et KJ
        if (cgStudy%option .eq. 'G' .or. cgStudy%option .eq. 'KJ') then
            call rsexch(' ', cgField%result_in, 'SIEF_ELGA', cgStudy%nume_ordre, chsig, iret)
!
            if (iret .ne. 0) then
!               Probleme pour recuperer SIEF_ELGA avec ce numero d'ordre
                call utmess('F', 'RUPTURE0_94', si=cgStudy%nume_ordre)
            end if
!
            lpain(nchin+1) = 'PCONTGR'
            lchin(nchin+1) = chsig
            nchin = nchin+1
        end if
!
!       Sommation des G élémentaires
        call cpu_time(start0)
        call calcul('S', opti, modelLigrel, nchin, lchin, &
                    lpain, 1, lchout, lpaout, 'V', 'OUI')
        call cpu_time(finish0)
        cgStat%cgCmpGtheta_te = cgStat%cgCmpGtheta_te+finish0-start0
        cgStat%nb_cgCmpGtheta_te = cgStat%nb_cgCmpGtheta_te+1
!
!       Somme des G, K1, K2, K3, FIC1, FIC2, FIC3 pour le champ theta actuel
        call cpu_time(start0)
        call mesomm(lchout(1), 7, vr=gth)
        call cpu_time(finish0)
        cgStat%cgCmpGtheta_mes = cgStat%cgCmpGtheta_mes+finish0-start0
        cgStat%nb_cgCmpGtheta_mes = cgStat%nb_cgCmpGtheta_mes+1
!       En 2D, gth contient directement les valeurs à imprimer
!       En 3D, il faut déterminer gth pour chaque noeud du front de fissure
!       (inversion du système linéaire A.G(s)=G(theta) par gsyste plus loin)

        if (cgField%ndim .eq. 3) then
!       En 3D, on stocke pour chacun des theta,  G, K1, K2, K3, FIC1, FIC2, FIC3
!       dans des vecteurs temporaires dédiés
            gthi(i_theta) = gth(1)
            k1th(i_theta) = gth(2)
            k2th(i_theta) = gth(3)
            k3th(i_theta) = gth(4)
            g1th(i_theta) = gth(5)
            g2th(i_theta) = gth(6)
            g3th(i_theta) = gth(7)
        end if
!
    end do
!
!    Cas axis, on normalise par 1/R
    if (cgStudy%l_axis) then
        call cgTheta%getBaseLoc(v_base)
        gth(1:7) = gth(1:7)/v_base(1)
    end if
!
!    Cas 3D, on détermine G(s) et les K(s)
    if (cgField%ndim .eq. 3) then
        if (cgTheta%discretization .eq. 'LINEAIRE') then

!           Correction dans le cas des fonds fermés (LAGRANDE/LINEAIRE) : G(1)=G(N)
            if (cgTheta%l_closed) then
                gthi(cgTheta%nb_theta_field) = gthi(1)
                k1th(cgTheta%nb_theta_field) = k1th(1)
                k2th(cgTheta%nb_theta_field) = k2th(1)
                k3th(cgTheta%nb_theta_field) = k3th(1)
                g1th(cgTheta%nb_theta_field) = g1th(1)
                g2th(cgTheta%nb_theta_field) = g2th(1)
                g3th(cgTheta%nb_theta_field) = g3th(1)
            end if
!
!           CORRECTION VALEURS EXTREMITES
            if (.not. cgTheta%l_closed) then
                if (cgTheta%nb_theta_field .ne. 2) then

                    s1 = v_basf(1)
                    s2 = v_basf(2)
                    s3 = v_basf(3)
                    sn2 = v_basf(cgTheta%nb_theta_field-3)
                    sn1 = v_basf(cgTheta%nb_theta_field-2)
                    sn = v_basf(cgTheta%nb_theta_field-1)

!                    CORRECTION DANS LE CAS LINEAIRE
                    if (.not. cgTheta%milieu) then
                        gthi(1) = gthi(2)*(s2-s1)/(s3-s1)
                        k1th(1) = k1th(2)*(s2-s1)/(s3-s1)
                        k2th(1) = k2th(2)*(s2-s1)/(s3-s1)
                        k3th(1) = k3th(2)*(s2-s1)/(s3-s1)
                        g1th(1) = g1th(2)*(s2-s1)/(s3-s1)
                        g2th(1) = g2th(2)*(s2-s1)/(s3-s1)
                        g3th(1) = g3th(2)*(s2-s1)/(s3-s1)
                        gthi(cgTheta%nb_theta_field) = gthi(cgTheta%nb_theta_field-1) &
                                                       *(sn-sn1)/(sn-sn2)
                        k1th(cgTheta%nb_theta_field) = k1th(cgTheta%nb_theta_field-1) &
                                                       *(sn-sn1)/(sn-sn2)
                        k2th(cgTheta%nb_theta_field) = k2th(cgTheta%nb_theta_field-1) &
                                                       *(sn-sn1)/(sn-sn2)
                        k3th(cgTheta%nb_theta_field) = k3th(cgTheta%nb_theta_field-1) &
                                                       *(sn-sn1)/(sn-sn2)
                        g1th(cgTheta%nb_theta_field) = g1th(cgTheta%nb_theta_field-1) &
                                                       *(sn-sn1)/(sn-sn2)
                        g2th(cgTheta%nb_theta_field) = g2th(cgTheta%nb_theta_field-1) &
                                                       *(sn-sn1)/(sn-sn2)
                        g3th(cgTheta%nb_theta_field) = g3th(cgTheta%nb_theta_field-1) &
                                                       *(sn-sn1)/(sn-sn2)

!                    CORRECTION DANS LE CAS QUADRATIQUE
                    else if (cgTheta%milieu) then
                        gthi(1) = gthi(2)/4.d0
                        k1th(1) = k1th(2)/4.d0
                        k2th(1) = k2th(2)/4.d0
                        k3th(1) = k3th(2)/4.d0
                        g1th(1) = g1th(2)/4.d0
                        g2th(1) = g2th(2)/4.d0
                        g3th(1) = g3th(2)/4.d0
                        gthi(cgTheta%nb_theta_field) = gthi(cgTheta%nb_theta_field-1)/4.d0
                        k1th(cgTheta%nb_theta_field) = k1th(cgTheta%nb_theta_field-1)/4.d0
                        k2th(cgTheta%nb_theta_field) = k2th(cgTheta%nb_theta_field-1)/4.d0
                        k3th(cgTheta%nb_theta_field) = k3th(cgTheta%nb_theta_field-1)/4.d0
                        g1th(cgTheta%nb_theta_field) = g1th(cgTheta%nb_theta_field-1)/4.d0
                        g2th(cgTheta%nb_theta_field) = g2th(cgTheta%nb_theta_field-1)/4.d0
                        g3th(cgTheta%nb_theta_field) = g3th(cgTheta%nb_theta_field-1)/4.d0
                    end if

                end if
            end if
!
            !       On inverse les systèmes linéaires A.G(s)=G(theta)
            call cpu_time(start0)

            !       SYSTEME LINEAIRE:  MATR*GS = GTHI
            call gsyste(cgTheta%matrix, cgTheta%nb_theta_field, cgTheta%nnof, gthi, gs)

            !       SYSTEME LINEAIRE:  MATR*K1S = K1TH
            call gsyste(cgTheta%matrix, cgTheta%nb_theta_field, cgTheta%nnof, k1th, k1s)

            !       SYSTEME LINEAIRE:  MATR*K2S = K2TH
            call gsyste(cgTheta%matrix, cgTheta%nb_theta_field, cgTheta%nnof, k2th, k2s)

            !       SYSTEME LINEAIRE:  MATR*K3S = K3TH
            call gsyste(cgTheta%matrix, cgTheta%nb_theta_field, cgTheta%nnof, k3th, k3s)

            !       SYSTEMES LINEAIRES POUR GIRWIN
            call gsyste(cgTheta%matrix, cgTheta%nb_theta_field, cgTheta%nnof, g1th, g1s)
            call gsyste(cgTheta%matrix, cgTheta%nb_theta_field, cgTheta%nnof, g2th, g2s)
            call gsyste(cgTheta%matrix, cgTheta%nb_theta_field, cgTheta%nnof, g3th, g3s)
!
            if (cgTheta%milieu .and. cgTheta%nb_point_fond .eq. 0 &
                .and. cgTheta%nnof/2-1 .gt. 0) then
                call hatSmooth(cgTheta%nnof, cgTheta%nnof/2+1, v_basf, gs)
                call hatSmooth(cgTheta%nnof, cgTheta%nnof/2+1, v_basf, k1s)
                call hatSmooth(cgTheta%nnof, cgTheta%nnof/2+1, v_basf, k2s)
                call hatSmooth(cgTheta%nnof, cgTheta%nnof/2+1, v_basf, k3s)
                call hatSmooth(cgTheta%nnof, cgTheta%nnof/2+1, v_basf, g1s)
                call hatSmooth(cgTheta%nnof, cgTheta%nnof/2+1, v_basf, g2s)
                call hatSmooth(cgTheta%nnof, cgTheta%nnof/2+1, v_basf, g3s)
            end if
!
            call cpu_time(finish0)
            cgStat%cgCmpGtheta_sys = cgStat%cgCmpGtheta_sys+finish0-start0
            cgStat%nb_cgCmpGtheta_sys = cgStat%nb_cgCmpGtheta_sys+7
!
        else if (cgTheta%discretization .eq. 'LEGENDRE') then
!       On évalue G(s) grâce aux polynômes de Legendre
!       Récupération des valeurs des poluynomes de Legendre pour les  noeuds
!       du fond de fissure
            call wkvect('&&OP0027.LEGENDRE', 'V V R8', &
                        cgTheta%nb_theta_field*cgTheta%nnof, iadrt3)
            call glegen(cgTheta%degree, cgTheta%nnof, cgTheta%lonfis, &
                        cgTheta%absfond, zr(iadrt3))

            do i = 1, cgTheta%nnof
!
                som(:) = 0.d0
!
                do j = 1, cgTheta%nb_theta_field
                    som(1) = som(1)+gthi(j)*zr(iadrt3+(j-1)*cgTheta%nnof+i-1)
                    som(2) = som(2)+k1th(j)*zr(iadrt3+(j-1)*cgTheta%nnof+i-1)
                    som(3) = som(3)+k2th(j)*zr(iadrt3+(j-1)*cgTheta%nnof+i-1)
                    som(4) = som(4)+k3th(j)*zr(iadrt3+(j-1)*cgTheta%nnof+i-1)
                    som(5) = som(5)+g1th(j)*zr(iadrt3+(j-1)*cgTheta%nnof+i-1)
                    som(6) = som(6)+g2th(j)*zr(iadrt3+(j-1)*cgTheta%nnof+i-1)
                    som(7) = som(7)+g3th(j)*zr(iadrt3+(j-1)*cgTheta%nnof+i-1)
                end do
!
                gs(i) = som(1)
                k1s(i) = som(2)
                k2s(i) = som(3)
                k3s(i) = som(4)
                g1s(i) = som(5)
                g2s(i) = som(6)
                g3s(i) = som(7)
                gis(i) = g1s(i)*g1s(i)+g2s(i)*g2s(i)+g3s(i)*g3s(i)
            end do
            call jedetr('&&OP0027.LEGENDRE')
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
!   Ajout des valeurs dans la table de G
!    call cgTempNodes(cgStudy, cgTable)

    do i = 1, cgTheta%nnof
        if (cgField%ndim .eq. 3) then
!       En 2D,  gth est déjà rempli. En 3D, on le remplit pour le noeud courant
            gth(1:7) = [gs(i), k1s(i), k2s(i), k3s(i), g1s(i), g2s(i), g3s(i)]
        end if

!       On recopie gth dans cgStudy%gth avec prise en compte de la symétrie
        if (cgTheta%symech .eq. 'OUI') then
            cgStudy%gth(1:7) = [2.d0*gth(1), 2.d0*gth(2), 0.d0, 0.d0, &
                                2.d0*gth(5), 0.d0, 0.d0]
        else
            cgStudy%gth(1:7) = gth(1:7)
        end if
        call cgTable%addValues(cgField, cgStudy, i)
    end do

!
    call detrsd('CHAMP_GD', chvarc)
    call detrsd('CHAMP_GD', chvref)
    call detrsd('CHAMP_GD', cf1d2d)
    call detrsd('CHAMP_GD', cf2d3d)
    call detrsd('CHAMP_GD', chepsi)
    call detrsd('CHAMP_GD', chpesa)
    call detrsd('CHAMP_GD', chpres)
    call detrsd('CHAMP_GD', chrota)
    call detrsd('CHAMP_GD', chtime)
    call detrsd('CHAMP_GD', chvolu)
    call detrsd('CHAMP_GD', chpuls)
    call detrsd('CHAMP_GD', chgtheta)
    call detrsd('CHAMP_GD', sigelno)
    call detrsd('CHAMP_GD', celmod)
    call detrsd('CHAMP_GD', chsdeg)
    call detrsd('CHAMP_GD', chslag)
    call detrsd('CHAMP_GD', chcer)
    call detrsd('CHAMP_GD', cheli)
    call detrsd('CHAM_ELEM_S', chscpl)
    call detrsd('CHAMP_GD', iscplan)
!
! --- A ne pas supprimer car on supprime un champ externe sinon
!    call detrsd('CHAMP_GD', chsigi)
!
    call jedema()
!
    call cpu_time(finish)
    cgStat%cgCmpGtheta = cgStat%cgCmpGtheta+finish-start
    cgStat%nb_cgCmpGtheta = cgStat%nb_cgCmpGtheta+1
end subroutine
