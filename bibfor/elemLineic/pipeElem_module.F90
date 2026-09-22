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
! aslint: disable=W1306, W0413
!
! ==================================================================================================
!
! Module for pipe elements
!
! ==================================================================================================
!
module pipeElem_module
! ==================================================================================================
    use beamElem_type
    use pipeElem_type
    use beamElem_module
    use loadElemCompute_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: pipeGetDime
    public :: pipeCoorElga, pipeMassIner, pipeForcNoda
    public :: pipeEfgeElgaNlin, pipeEfgeElgaLine, pipeEfgeElno
    public :: pipeLoadPres, pipeLoadLine
    public :: pipeGetProperties, pipeGetSubPoints
    public :: pipeGetRadiusLayer, pipeGetJacobian
    public :: pipeGetElasProp, pipeGetDensity, pipeGetDisp, pipeGetTherProp
    public :: pipeBaseForMatr, pipeBaseForVect
    public :: pipeNMatr, pipeBMatr
    public :: pipeCheckMetric
    private :: pipeGetType, pipeGetBases, pipeGetGeometry
    private :: pipeBMatrStraight, pipeBMatrCurved
    private :: shellBMatrStraight, shellBMatrCurved, pipeGetSection
! ==================================================================================================
    private
#include "asterc/r8pi.h"
#include "asterf_types.h"
#include "asterfort/angvxy.h"
#include "asterfort/assert.h"
#include "asterfort/beamElem_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/normev.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/ppga1d.h"
#include "asterfort/ppgan2.h"
#include "asterfort/prmave.h"
#include "asterfort/promat.h"
#include "asterfort/provec.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/verifg.h"
#include "jeveux.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! pipeCoorElga
!
! Compute coordinates of integration points
!
! In  nbNode           : number of nodes of element
! In  nbFourier        : number of Fourier modes
! In  npg              : number of Gauss points
! In  ipoids           : adress for Gauss weights
! In  ivf              : adress for shape functions
! In  idfde            : adress for derivatives of shape functions
! In  jvGeom           : adress for geometry
! In  jvCoorPg         : adress for coordonnées points de gauss aux sous-points + poids (section)
! In  jvSuppPg         : adress for coordonnées des points de Gauss du support (segment)
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeCoorElga(gauss_support, &
                            nbNode, npg, nbFourier, &
                            ipoids, ivf, idfde, jvGeom, &
                            jvCoorPg, jvSuppPg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: gauss_support
        integer(kind=8), intent(in) :: nbNode, npg, nbFourier
        integer(kind=8), intent(in) :: ipoids, ivf, idfde, jvGeom
        integer(kind=8), intent(in) :: jvCoorPg, jvSuppPg
! ----- Local
        integer(kind=8), parameter :: ndim = 3
        integer(kind=8) :: nbLayer, nbSect
        integer(kind=8) :: nspgLayer, nspgSect, nspg
        real(kind=8) :: weightLayer(2*PIPE_MAX_LAYERS+1)
        real(kind=8) :: weightSect(2*PIPE_MAX_SECTORS+1)
        real(kind=8) :: alpha, layerThickness
        real(kind=8) :: radiusLayer, areaLayer
        real(kind=8) :: yy, zz
        integer(kind=8) :: kpg, kspgLayer, kspgSect
        real(kind=8) :: copg(PIPE_NB_COORELGA, PIPE_MAX_NPG), pgl(3, 3), gm1(3), gm2(3)
        type(pipeElem_Prop) :: pipeElem
        type(sectPipe_Prop) :: sectPipe
        type(beamElem_Prop) :: beamElem
!   ------------------------------------------------------------------------------------------------
!

! ----- Get parameters about layers and sections
        call pipeGetSubPoints(nbLayer, nbSect, &
                              nspgLayer, nspgSect, nspg, &
                              weightLayer, weightSect)

! ----- Get properties of pipe
        call pipeGetProperties(nbNode, nbFourier, pipeElem)
        beamElem = pipeElem%beamElem
        sectPipe = beamElem%sectPipe

! ----- Position et poids des points de gauss de l'élément support (segment)
        call ppga1d(ndim, nbNode, npg, zr(ipoids), zr(ivf), zr(idfde), zr(jvGeom), copg)

        if (gauss_support) then
            do kpg = 1, npg
                zr(jvSuppPg+(kpg-1)*PIPE_NB_COORELGA+0) = copg(1, kpg)
                zr(jvSuppPg+(kpg-1)*PIPE_NB_COORELGA+1) = copg(2, kpg)
                zr(jvSuppPg+(kpg-1)*PIPE_NB_COORELGA+2) = copg(3, kpg)
                zr(jvSuppPg+(kpg-1)*PIPE_NB_COORELGA+3) = copg(4, kpg)
            end do
        end if
!
        gm1 = 0.d0
        alpha = r8pi()/(nbSect)
        layerThickness = sectPipe%thickness/(2.0d0*nbLayer)
        do kpg = 1, npg
!           L'orientation change en fct du point de Gauss dans le cas courbe
            if (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
                pgl = beamElem%pglCell(:, :)
            elseif (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
                pgl = beamElem%pgl(:, :, kpg)
            else
                ASSERT(ASTER_FALSE)
            end if

!           Calcul des coordonnees et stockage. Les sous points sont stockes niveau par niveau.
!           Il y a plusieurs niveaux par couche, en commencant par la section z local = 0 et y >0
            do kspgLayer = 1, 2*nbLayer+1

!               Section concernant le sous-point
                radiusLayer = sectPipe%radiusExt-sectPipe%thickness+(kspgLayer-1)*layerThickness
                areaLayer = radiusLayer*layerThickness*alpha
                do kspgSect = 1, 2*nbSect+1
!                   SUPER IMPORTANT SUPER IMPORTANT SUPER IMPORTANT
!                   La convention des angles de vrilles entre les poutres et tuyaux est différente
!                   Il y radiusPipeMoy un repère indirect pour les tuyaux ==> c'est pas bien
!                       - On décale les angles de 90°.
!                       - Quand tout sera dans jDof'ordre, il faudra calculer correctement yy et zz
!                   A FAIRE DANS : te0478  irmase
!
                    yy = cos(-(kspgSect-1)*alpha-0.5*r8pi())
                    zz = sin(-(kspgSect-1)*alpha-0.5*r8pi())

!                   Position de SP dans la section
                    gm1(2) = radiusLayer*yy
                    gm1(3) = radiusLayer*zz
                    call utpvlg(1, 3, pgl, gm1, gm2)
                    zr(jvCoorPg+((kpg-1)*nspg+(kspgLayer-1)*(2*nbSect+1)+(kspgSect-1))*4+0) = &
                        copg(1, kpg)+gm2(1)
                    zr(jvCoorPg+((kpg-1)*nspg+(kspgLayer-1)*(2*nbSect+1)+(kspgSect-1))*4+1) = &
                        copg(2, kpg)+gm2(2)
                    zr(jvCoorPg+((kpg-1)*nspg+(kspgLayer-1)*(2*nbSect+1)+(kspgSect-1))*4+2) = &
                        copg(3, kpg)+gm2(3)
!                   Pour le poids
                    zr(jvCoorPg+((kpg-1)*nspg+(kspgLayer-1)*(2*nbSect+1)+(kspgSect-1))*4+3) = &
                        copg(4, kpg)*weightSect(kspgSect)*weightLayer(kspgLayer)*areaLayer
                end do
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeMassIner
!
! Compute MASS_INER
!
! In  nbNode           : number of nodes of element
! In  nbFourier        : number of Fourier modes
! In  rho              : density
! In  jvMassIner       : adress for mass/inerty to write
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeMassIner(nbNode, nbFourier, rho, jvMassIner)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nbNode, nbFourier
        real(kind=8), intent(in) :: rho
        integer(kind=8), intent(in) :: jvMassIner
! ----- Local
        real(kind=8) :: xa, xb, xl, iy1, iz1
        real(kind=8) :: pgl(3, 3), pgl3(3, 3)
        real(kind=8) :: cdgLoca(3), cdgGlob(3), cdgGlobT(3)
        real(kind=8) :: mass
        real(kind=8) :: angl(3), coo1(3), coo2(3), coo3(3)
        real(kind=8) :: tau1(3), tau2(3), norm(3), x3(3), y3(3)
        real(kind=8) :: normTau1, normTau2, normNorm
        integer(kind=8) :: iDime, jvGeom
        real(kind=8) :: massInerGlob(6), massInerLoca(6)
        type(pipeElem_Prop) :: pipeElem
        type(sectPipe_Prop) :: sectPipe
        type(beamElem_Prop) :: beamElem
!   ------------------------------------------------------------------------------------------------
!
        call jevech('PGEOMER', 'L', jvGeom)

! ----- Get properties of pipe
        call pipeGetProperties(nbNode, nbFourier, pipeElem)
        beamElem = pipeElem%beamElem
        sectPipe = beamElem%sectPipe

! ----- Principal inerties
        iy1 = r8pi()*(sectPipe%radiusExt**4-sectPipe%radiusInt**4)/4.d0
        iz1 = iy1

! ----- Length of pipe
        xl = beamElem%elemLength
        if (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
            xl = beamElem%thetaElbow*beamElem%radiusElbow
        end if

! ----- Get basis at nodes of pipe
        pgl = 0.d0
        pgl3 = 0.d0
        if (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
            if (nbNode .eq. 4) then
                angl = 0.d0
                do iDime = 1, 3
                    coo1(iDime) = zr(jvGeom-1+iDime)
                    coo2(iDime) = zr(jvGeom-1+3+iDime)
                    coo3(iDime) = (zr(jvGeom-1+6+iDime)+zr(jvGeom-1+9+iDime))*0.5d0
                    tau1(iDime) = coo3(iDime)-coo1(iDime)
                    tau2(iDime) = coo2(iDime)-coo3(iDime)
                    x3(iDime) = coo2(iDime)-coo1(iDime)
                end do
                call normev(tau1, normTau1)
                call normev(tau2, normTau2)
                call provec(tau2, tau1, norm)
                call normev(norm, normNorm)
                call provec(x3, norm, y3)
                call angvxy(x3, y3, angl)
                call matrot(angl, pgl3)
            else
                pgl3 = pipeElem%beamElem%pgl(:, :, 3)
            end if
        elseif (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
            pgl = pipeElem%beamElem%pglCell
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Calcul masse
        mass = rho*sectPipe%area*xl

! ----- Calcul du centre de gravité
        cdgLoca = 0.d0
        xb = 0.d0
        if (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
            xb = 1.d0+ &
                 (sectPipe%radiusMoy**2+sectPipe%thickness**2/4.d0)/ &
                 (2.d0*beamElem%radiusElbow**2)
            cdgLoca(2) = -beamElem%radiusElbow* &
                         (sin(beamElem%thetaElbow/2.d0)/ &
                          (beamElem%thetaElbow/2.d0)*xb-cos(beamElem%thetaElbow/2.d0))
        end if

! ----- From local to global
        cdgGlobT = 0.d0
        if (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
            call utpvlg(1, 3, pgl3, cdgLoca, cdgGlobT)
        elseif (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
            call utpvlg(1, 3, pgl, cdgLoca, cdgGlobT)
        else
            ASSERT(ASTER_FALSE)
        end if
        cdgGlob = 0.d0
        cdgGlob(1) = cdgGlobT(1)+(zr(jvGeom-1+4)+zr(jvGeom-1+1))/2.d0
        cdgGlob(2) = cdgGlobT(2)+(zr(jvGeom-1+5)+zr(jvGeom-1+2))/2.d0
        cdgGlob(3) = cdgGlobT(3)+(zr(jvGeom-1+6)+zr(jvGeom-1+3))/2.d0

! ----- Inertie de l'élément
        if (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
            xa = (sectPipe%area*sectPipe%radiusExt**2+3.0d0*iz1)
            xb = sectPipe%radiusExt*sin(beamElem%thetaElbow/2.d0)/beamElem%thetaElbow/2.d0* &
                 (1.d0+( &
                  sectPipe%radiusMoy*sectPipe%radiusMoy+ &
                  sectPipe%thickness*sectPipe%thickness/4.d0 &
                  )/ &
                  (2.d0*sectPipe%radiusExt**2))
            massInerLoca(1) = rho*xl*(iy1+xa*(0.5d0+ &
                                              sin(beamElem%thetaElbow)/ &
                                              (4.d0*beamElem%thetaElbow)))- &
                              mass*xb*xb
            massInerLoca(2) = 0.d0
            massInerLoca(3) = rho*xl*(iy1+xa*(0.5d0- &
                                              sin(beamElem%thetaElbow)/(4.d0*beamElem%thetaElbow)))
            massInerLoca(4) = 0.d0
            massInerLoca(5) = 0.d0
            massInerLoca(6) = rho*xl*xa-mass*xb*xb
            call utpslg(1, 3, pgl3, massInerLoca, massInerGlob)
        elseif (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
            massInerLoca(1) = rho*(iy1+iz1)*xl
            massInerLoca(2) = 0.d0
            massInerLoca(3) = rho*xl*(iy1+sectPipe%area*xl*xl/12.d0)
            massInerLoca(4) = 0.d0
            massInerLoca(5) = 0.d0
            massInerLoca(6) = rho*xl*(iz1+sectPipe%area*xl*xl/12.d0)
            call utpslg(1, 3, pgl, massInerLoca, massInerGlob)
        else
            ASSERT(ASTER_FALSE)
        end if
        zr(jvMassIner) = mass
        zr(jvMassIner+1) = cdgGlob(1)
        zr(jvMassIner+2) = cdgGlob(2)
        zr(jvMassIner+3) = cdgGlob(3)
        zr(jvMassIner+4) = massInerGlob(1)
        zr(jvMassIner+5) = massInerGlob(3)
        zr(jvMassIner+6) = massInerGlob(6)
        zr(jvMassIner+7) = massInerGlob(2)
        zr(jvMassIner+8) = massInerGlob(4)
        zr(jvMassIner+9) = massInerGlob(5)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetSection
!
! Get section of pipe
!
! Out sectPipe         : properties of pipe's section
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetSection(sectPipe)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(sectPipe_Prop), intent(out) :: sectPipe
! ----- Local
        integer(kind=8), parameter :: nbCara = 2
        character(len=8), parameter :: caraName(nbCara) = (/'R1 ', 'EP1'/)
        real(kind=8) :: caraVale(nbCara)
!   ------------------------------------------------------------------------------------------------
!
        call poutre_modloc('CAGEP1', caraName, nbCara, lvaleur=caraVale)
        sectPipe%radiusExt = caraVale(1)
        sectPipe%thickness = caraVale(2)
        sectPipe%radiusMoy = sectPipe%radiusExt-sectPipe%thickness/2.d0
        sectPipe%radiusInt = sectPipe%radiusExt-sectPipe%thickness
        sectPipe%area = r8pi()*(sectPipe%radiusExt**2-sectPipe%radiusInt**2)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetType
!
! Get type of pipe
!
! In  jvPipe           : adress for parameters of pipe
! In  nbNode           : number of nodes of element
! Out pipeType         : type of pipe (straight or elbow)
! Out lModiMetric      : using MODI_METRIQUE
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetType(jvPipe, nbNode, pipeType, lModiMetric)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvPipe
        integer(kind=8), intent(in) :: nbNode
        integer(kind=8), intent(out) :: pipeType
        aster_logical, intent(out) :: lModiMetric
! ----- Local
        integer(kind=8) :: pipeTypeInfo, iCmp
!   ------------------------------------------------------------------------------------------------
!
        pipeType = PIPE_TYPE_UNDEF
        if (nbNode .eq. 3) then
            iCmp = 9
        else if (nbNode .eq. 4) then
            iCmp = 12
        else
            ASSERT(ASTER_FALSE)
        end if
        pipeTypeInfo = nint(zr(jvPipe-1+iCmp+1))
        lModiMetric = pipeTypeInfo .lt. 10
        if (pipeTypeInfo .ge. 10) then
            pipeTypeInfo = pipeTypeInfo-10
        end if
        if (pipeTypeInfo .eq. 0) then
            pipeType = PIPE_TYPE_STRAIGHT
        elseif (pipeTypeInfo .eq. 1) then
            pipeType = PIPE_TYPE_ELBOW
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetBases
!
! Get basis at nodes of pipe
!
! In  jvPipe           : adress for parameters of pipe
! In  nbNode           : number of nodes of element
! In  pipeType         : type of pipe
! Out pgl              : matrix for basis definition (straight pipe)
! Out pgl1             : matrix for basis definition (elbow pipe, at node 1)
! Out pgl2             : matrix for basis definition (elbow pipe, at node 2)
! Out pgl3             : matrix for basis definition (elbow pipe, at node 3)
! Out pgl4             : matrix for basis definition (elbow pipe, at node 4, only for SEG4 support)
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetBases(jvPipe, nbNode, pipeType, &
                            pgl, pgl1, pgl2, pgl3, pgl4)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvPipe, nbNode, pipeType
        real(kind=8), intent(out) :: pgl(3, 3)
        real(kind=8), intent(out) :: pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), pgl4(3, 3)
! ----- Local
        integer(kind=8) :: iCmp
        real(kind=8) :: angl1(3), angl2(3), angl3(3), angl4(3)
!   ------------------------------------------------------------------------------------------------
!
        pgl = 0.d0
        pgl1 = 0.d0
        pgl2 = 0.d0
        pgl3 = 0.d0
        pgl4 = 0.d0
        if (pipeType .eq. PIPE_TYPE_STRAIGHT) then
            do iCmp = 1, 3
                angl1(iCmp) = zr(jvPipe-1+iCmp)
            end do
            call matrot(angl1, pgl)
        elseif (pipeType .eq. PIPE_TYPE_ELBOW) then
            do iCmp = 1, 3
                angl1(iCmp) = zr(jvPipe-1+iCmp)
                angl2(iCmp) = zr(jvPipe-1+iCmp+3)
                angl3(iCmp) = zr(jvPipe-1+iCmp+6)
            end do
            call matrot(angl1, pgl1)
            call matrot(angl2, pgl2)
            call matrot(angl3, pgl3)
            if (nbNode .eq. 4) then
                do iCmp = 1, 3
                    angl4(iCmp) = zr(jvPipe-1+iCmp+9)
                end do
                call matrot(angl4, pgl4)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetGeometry
!
! Get geometry of pipe
!
! In  jvPipe           : adress for parameters of pipe
! In  nbNode           : number of nodes of element
! Out thetaElbow       : angle of elbow
! Out radiusElbow      : radius of elbow
! Out omegaElbow       : ANGLE ENTRE N ET LA GENERATRICE
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetGeometry(jvPipe, nbNode, &
                               thetaElbow, radiusElbow, omegaElbow, pipeLength)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvPipe
        integer(kind=8), intent(in) :: nbNode
        real(kind=8), intent(out) :: thetaElbow, radiusElbow, omegaElbow, pipeLength
! ----- Local
        integer(kind=8) :: iCmp
!   ------------------------------------------------------------------------------------------------
!
        if (nbNode .eq. 3) then
            iCmp = 9
        else if (nbNode .eq. 4) then
            iCmp = 12
        else
            ASSERT(ASTER_FALSE)
        end if
        pipeLength = zr(jvPipe-1+iCmp+2)
        radiusElbow = zr(jvPipe-1+iCmp+3)
        thetaElbow = zr(jvPipe-1+iCmp+4)
        omegaElbow = zr(jvPipe-1+iCmp+5)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetElasProp
!
! Get elasticity parameters for pipes
!
! Use temperature at point => give kpg AND kspg
! Use given temperature => give temp
!
! In  jvMaterCode      : adress for material parameters
! In  temp             : current temperature
! In  kpg              : current Gauss point
! In  kspg             : current Gauss sub-point
! Out c                : elasticity matrix for pipe
! Out cisail           : shear constant
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetElasProp(jvMaterCode, &
                               kpg_, kspg_, temp_, &
                               c_, cisail_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvMaterCode
        integer(kind=8), optional, intent(in) :: kpg_, kspg_
        real(kind=8), optional, intent(in) :: temp_
        real(kind=8), optional, intent(out) :: c_(PIPE_TENS_SIZE, PIPE_TENS_SIZE)
        real(kind=8), optional, intent(out) :: cisail_
! ----- Local
        character(len=4), parameter :: fami = "RIGI"
        integer(kind=8), parameter :: nbProp = 2, nbPara = 1
        character(len=8), parameter :: paraName = 'TEMP'
        character(len=16), parameter :: propName(nbProp) = (/'E ', 'NU'/)
        real(kind=8) :: propVale(nbProp)
        integer(kind=8) :: propError(nbProp)
        character(len=16) :: elasKeyword
        aster_logical :: lTempAtPg
        real(kind=8) :: youngModulus, poissonRatio
        real(kind=8) :: beta, g, cisail, c(PIPE_TENS_SIZE, PIPE_TENS_SIZE)
!   ------------------------------------------------------------------------------------------------
!
        c = 0.d0
        cisail = 0.d0

! ----- Get status of temperature
        lTempAtPg = ASTER_FALSE
        if (present(kpg_)) then
            lTempAtPg = ASTER_TRUE
            ASSERT(present(kspg_))
        else
            ASSERT(present(temp_))
        end if

! ----- Type of elasticity
        call rccoma(zi(jvMaterCode), 'ELAS', 1, elasKeyword)
        if (elasKeyword .ne. 'ELAS') then
            call utmess('F', 'PIPE1_46', sk=elasKeyword)
        end if

! ----- Get parameters of elasticity
        if (lTempAtPg) then
            call rcvalb(fami, kpg_, kspg_, '+', &
                        zi(jvMaterCode), ' ', elasKeyword, &
                        0, ' ', [0.d0], &
                        nbProp, propName, propVale, &
                        propError, 1)
        else
            call rcvala(zi(jvMaterCode), ' ', elasKeyword, &
                        nbPara, paraName, [temp_], &
                        nbProp, propName, propVale, &
                        propError, 1)
        end if
        youngModulus = propVale(1)
        poissonRatio = propVale(2)
        beta = 1.d0/(1.d0-poissonRatio**2)
        cisail = youngModulus/(1.d0+poissonRatio)
        g = 1.d0/(2.d0*(1.d0+poissonRatio))
        c(1, 1) = beta
        c(1, 2) = poissonRatio*beta
        c(1, 3) = 0.d0
        c(1, 4) = 0.d0

        c(2, 1) = poissonRatio*beta
        c(2, 2) = beta
        c(2, 3) = 0.d0
        c(2, 4) = 0.d0

        c(3, 1) = 0.d0
        c(3, 2) = 0.d0
        c(3, 3) = g
        c(3, 4) = 0.d0

        c(4, 1) = 0.d0
        c(4, 2) = 0.d0
        c(4, 3) = 0.d0
        c(4, 4) = g

        c = youngModulus*c

        if (present(cisail_)) then
            cisail_ = cisail
        end if
        if (present(c_)) then
            c_ = c
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetTherProp
!
! Get thermal parameters for pipes
!
! In  jvMaterCode      : adress for material parameters
! In  temp             : current temperature
! Out c                : elasticity matrix for pipe
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetTherProp(jvMaterCode, temp, c)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvMaterCode
        real(kind=8), intent(in) :: temp
        real(kind=8), intent(out) :: c(2, 2)
! ----- Local
        character(len=4), parameter :: fami = "RIGI"
        integer(kind=8), parameter :: nbProp = 2, nbPara = 1
        character(len=8), parameter :: paraName = 'TEMP'
        character(len=16), parameter :: propName(nbProp) = (/'E ', 'NU'/)
        real(kind=8) :: propVale(nbProp)
        integer(kind=8) :: propError(nbProp)
        character(len=16) :: elasKeyword
        real(kind=8) :: youngModulus, poissonRatio
        real(kind=8) :: beta
!   ------------------------------------------------------------------------------------------------
!
        c = 0.d0

! ----- Type of elasticity
        call rccoma(zi(jvMaterCode), 'ELAS', 1, elasKeyword)
        if (elasKeyword .ne. 'ELAS') then
            call utmess('F', 'PIPE1_46', sk=elasKeyword)
        end if

! ----- Get parameters of elasticity
        call rcvalb(fami, 1, 1, '+', &
                    zi(jvMaterCode), ' ', elasKeyword, &
                    nbPara, paraName, [temp], &
                    nbProp, propName, propVale, &
                    propError, 1)
        youngModulus = propVale(1)
        poissonRatio = propVale(2)
        beta = 1.d0/(1.d0-poissonRatio**2)
!
        c(1, 1) = youngModulus*beta
        c(1, 2) = youngModulus*poissonRatio*beta
        c(2, 1) = youngModulus*poissonRatio*beta
        c(2, 2) = youngModulus*beta
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetSubPoints
!
! Get sub-points for pipe elements
!
! Out nbLayer          : number of layers in layers of pipe
! Out nbSect           : number of sectors in sectors of pipe
! Out nspgLayer        : number of sub-points in layers of pipe
! Out nspgSect         : number of sub-points in sectors of pipe
! Out nspg             : total number of sub-points
! Out weightLayer      : weight of integration scheme for layers
! Out weightSect       : weight of integration scheme for sectors
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetSubPoints(nbLayer_, nbSect_, &
                                nspgLayer_, nspgSect_, nspg_, &
                                weightLayer_, weightSect_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), optional, intent(out) :: nbLayer_, nbSect_
        integer(kind=8), optional, intent(out) :: nspgLayer_, nspgSect_, nspg_
        real(kind=8), optional, intent(out) :: weightLayer_(2*PIPE_MAX_LAYERS+1)
        real(kind=8), optional, intent(out) :: weightSect_(2*PIPE_MAX_SECTORS+1)
! ----- Local
        integer(kind=8) :: jvSupPoint
        integer(kind=8) :: iLayer, iSect
        integer(kind=8) :: nbLayer, nbSect
        integer(kind=8) :: nspgLayer, nspgSect, nspg
        real(kind=8) :: weightLayer(2*PIPE_MAX_LAYERS+1), weightSect(2*PIPE_MAX_SECTORS+1)
!   ------------------------------------------------------------------------------------------------
!
        call jevech('PNBSP_I', 'L', jvSupPoint)

! ----- Get parameters about layers and sections
        nbLayer = zi(jvSupPoint-1+1)
        nbSect = zi(jvSupPoint-1+2)
        ASSERT(nbLayer .le. PIPE_MAX_LAYERS)
        ASSERT(nbSect .le. PIPE_MAX_SECTORS)
        nspgLayer = 2*nbLayer+1
        nspgSect = 2*nbSect+1
        ASSERT(nspgLayer .le. 2*PIPE_MAX_LAYERS+1)
        ASSERT(nspgSect .le. 2*PIPE_MAX_SECTORS+1)
        nspg = nspgLayer*nspgSect

! ----- For DEBUG
        if (PIPE_DEBUG) then
            WRITE (6, *) "Pipe - Integration"
            WRITE (6, *) " nbLayer :", nbLayer
            WRITE (6, *) " nbSect :", nbSect
            WRITE (6, *) " nspgLayer :", nspgLayer
            WRITE (6, *) " nspgSect :", nspgSect
            WRITE (6, *) " nspg :", nspg
        end if

! ----- Get integration scheme for sub-points (layers/sectors)
        weightLayer(1) = 1.d0/3.d0
        do iLayer = 1, nbLayer-1
            weightLayer(2*iLayer) = 4.d0/3.d0
            weightLayer(2*iLayer+1) = 2.d0/3.d0
        end do
        weightLayer(2*nbLayer) = 4.d0/3.d0
        weightLayer(2*nbLayer+1) = 1.d0/3.d0
        weightSect(1) = 1.d0/3.d0
        do iSect = 1, nbSect-1
            weightSect(2*iSect) = 4.d0/3.d0
            weightSect(2*iSect+1) = 2.d0/3.d0
        end do
        weightSect(2*nbSect) = 4.d0/3.d0
        weightSect(2*nbSect+1) = 1.d0/3.d0

        if (present(weightLayer_)) then
            weightLayer_ = weightLayer
        end if
        if (present(weightSect_)) then
            weightSect_ = weightSect
        end if
        if (present(nbLayer_)) then
            nbLayer_ = nbLayer
        end if
        if (present(nbSect_)) then
            nbSect_ = nbSect
        end if
        if (present(nspgLayer_)) then
            nspgLayer_ = nspgLayer
        end if
        if (present(nspgSect_)) then
            nspgSect_ = nspgSect
        end if
        if (present(nspg_)) then
            nspg_ = nspg
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! shellBMatrStraight
!
! Compute B matrix (u => epsi) for straight pipe - Shell part
!
! In  pipeElem         : properties of pipe element element
! In  nbDof            : number of DOF in element
! In  kpg              : current Gauss point
! In  phi              : coordinate of sub-point in section
! In  zeta             : Coordinate of sub-point in thickness
! In  radiusLayer      : radius of current layer
! In  ff               : shape functions on beam
! In  df1              : first derivative of shape functions on beam
! In  df2              : second derivative of shape functions on beam
! IO  b                : B matrix (u => epsi)
!
! --------------------------------------------------------------------------------------------------
    subroutine shellBMatrStraight(pipeElem, nbDof, kpg, &
                                  phi, zeta, radiusLayer, &
                                  ff, df1, df2, &
                                  b)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbDof, kpg
        real(kind=8), intent(in) :: phi, zeta, radiusLayer
        real(kind=8), intent(in) :: ff(*), df1(*), df2(*)
        real(kind=8), intent(inout) :: b(PIPE_TENS_SIZE, nbDof)
! ----- Local
        integer(kind=8) :: iNode, iBloc, iFourier, iColon, nbNode, nbFourier
        real(kind=8) :: hk, dhk, d2hk, cosmfi, sinmfi, cosfi, sinfi
        real(kind=8) :: radiusPipeMoy, elemLength
!   ------------------------------------------------------------------------------------------------
!
        cosfi = cos(phi)
        sinfi = sin(phi)
        nbFourier = pipeElem%nbFourier
        nbNode = pipeElem%nbNode
        elemLength = pipeElem%beamElem%elemLength
        radiusPipeMoy = pipeElem%beamElem%sectPipe%radiusMoy

        do iNode = 1, nbNode
            hk = ff(nbNode*(kpg-1)+iNode)
            dhk = df1(nbNode*(kpg-1)+iNode)*(2.d0/elemLength)
            d2hk = df2(nbNode*(kpg-1)+iNode)*(2.d0/elemLength)*(2.d0/elemLength)

            iBloc = (9+6*(nbFourier-1))*(iNode-1)

            do iFourier = 2, nbFourier
                iColon = iBloc+6+3*(iFourier-2)
                cosmfi = cos(iFourier*phi)
                sinmfi = sin(iFourier*phi)

                b(1, iColon+1) = dhk*cosmfi
                b(1, iColon+2) = 0.d0
                b(1, iColon+3) = -zeta*d2hk*cosmfi

                b(2, iColon+1) = 0.d0
                b(2, iColon+2) = (iFourier/radiusLayer)*hk*cosmfi*(1.d0+zeta/radiusPipeMoy)
                b(2, iColon+3) = (1.d0/radiusLayer)*hk*cosmfi* &
                                 (1.d0+zeta*iFourier*iFourier/radiusPipeMoy)

                b(3, iColon+1) = -(iFourier/radiusLayer)*hk*sinmfi
                b(3, iColon+2) = dhk*sinmfi*(1+zeta/radiusPipeMoy)
                b(3, iColon+3) = zeta*iFourier*dhk*sinmfi*(1.d0/radiusPipeMoy+1.d0/radiusLayer)

                b(4, iColon+1) = 0.d0
                b(4, iColon+2) = 0.d0
                b(4, iColon+3) = 0.d0

            end do
            do iFourier = 2, nbFourier
                iColon = iBloc+6+3*(nbFourier-1)+3*(iFourier-2)
                cosmfi = cos(iFourier*phi)
                sinmfi = sin(iFourier*phi)

                b(1, iColon+1) = dhk*sinmfi
                b(1, iColon+2) = 0.d0
                b(1, iColon+3) = -zeta*d2hk*sinmfi

                b(2, iColon+1) = 0.d0
                b(2, iColon+2) = -(iFourier/radiusLayer)*hk*sinmfi*(1.d0+zeta/radiusPipeMoy)
                b(2, iColon+3) = (1.d0/radiusLayer)*hk*sinmfi* &
                                 (1.d0+zeta*iFourier*iFourier/radiusPipeMoy)

                b(3, iColon+1) = (iFourier/radiusLayer)*hk*cosmfi
                b(3, iColon+2) = dhk*cosmfi*(1.d0+zeta/radiusPipeMoy)
                b(3, iColon+3) = -zeta*iFourier*dhk*cosmfi*(1.d0/radiusPipeMoy+1.d0/radiusLayer)

                b(4, iColon+1) = 0.d0
                b(4, iColon+2) = 0.d0
                b(4, iColon+3) = 0.d0

            end do

            iColon = iBloc+6*(nbFourier-1)+6
            b(1, iColon+1) = -zeta*d2hk
            b(2, iColon+1) = hk/radiusLayer
            b(3, iColon+1) = 0.d0
            b(4, iColon+1) = 0.d0

            b(1, iColon+2) = -zeta*d2hk*cosfi
            b(2, iColon+2) = (2.d0/radiusLayer)*hk*cosfi*(1.d0+zeta/radiusPipeMoy)
            b(3, iColon+2) = dhk*sinfi*(1.d0+2.d0*zeta/radiusPipeMoy+zeta/radiusLayer)
            b(4, iColon+2) = 0.d0

            b(1, iColon+3) = -zeta*d2hk*sinfi
            b(2, iColon+3) = (2.d0/radiusLayer)*hk*sinfi*(1.d0+zeta/radiusPipeMoy)
            b(3, iColon+3) = -dhk*cosfi*(1.d0+2.d0*zeta/radiusPipeMoy+zeta/radiusLayer)
            b(4, iColon+3) = 0.d0

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! shellBMatrCurved
!
! Compute B matrix (u => epsi) for curved pipe - Shell part
!
! In  pipeElem         : properties of pipe element element
! In  nbDof            : number of DOF in element
! In  kpg              : current Gauss point
! In  npg              : number of Gauss points
! In  xpg              : coordinates of Gauss points
! In  phi              : coordinate of sub-point in section
! In  zeta             : Coordinate of sub-point in thickness
! In  radiusLayer      : radius of current layer
! In  ff               : shape functions on beam
! In  df1              : first derivative of shape functions on beam
! In  df2              : second derivative of shape functions on beam
! IO  b                : B matrix (u => epsi)
!
! --------------------------------------------------------------------------------------------------
    subroutine shellBMatrCurved(pipeElem, nbDof, &
                                kpg, npg, xpg, &
                                phi, zeta, radiusLayer, &
                                ff, df1, df2, &
                                b)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbDof, kpg, npg
        real(kind=8), intent(in) :: xpg(npg)
        real(kind=8), intent(in) :: phi, zeta, radiusLayer
        real(kind=8), intent(in) :: ff(*), df1(*), df2(*)
        real(kind=8), intent(inout) :: b(PIPE_TENS_SIZE, nbDof)
! ----- Local
        integer(kind=8) :: iNode, iBloc, iFourier, iColon, nbNode, nbFourier
        real(kind=8) :: hk, dhk, d2hk, cosfi, sinfi
        real(kind=8) :: cosmte, sinmte, coste, sinte
        real(kind=8) :: denr, dena, dent, tk(PIPE_MAX_NODE), ck, sk
        real(kind=8) :: radiusElbow, omegaElbow, thetaElbow
        real(kind=8) :: elemLength, radiusPipeMoy
!   ------------------------------------------------------------------------------------------------
!
        nbNode = pipeElem%nbNode
        nbFourier = pipeElem%nbFourier
        elemLength = pipeElem%beamElem%elemLength
        radiusElbow = pipeElem%beamElem%radiusElbow
        thetaElbow = pipeElem%beamElem%thetaElbow
        omegaElbow = pipeElem%beamElem%omegaElbow
        tk = pipeElem%beamElem%tk
        radiusPipeMoy = pipeElem%beamElem%sectPipe%radiusMoy

        coste = cos(phi-omegaElbow)
        sinte = sin(phi-omegaElbow)
        cosfi = cos(phi)
        sinfi = sin(phi)
        denr = radiusElbow+radiusLayer*sinfi
        dena = radiusElbow+radiusPipeMoy*sinfi
        dent = (1.d0/denr+radiusPipeMoy/(radiusLayer*dena))

        do iNode = 1, nbNode

            hk = ff(nbNode*(kpg-1)+iNode)
            dhk = df1(nbNode*(kpg-1)+iNode)*(2.d0/thetaElbow)
            d2hk = df2(nbNode*(kpg-1)+iNode)*(2.d0/thetaElbow)*(2.d0/thetaElbow)
            ck = cos((1.d0+xpg(kpg))*thetaElbow/2.d0-tk(iNode))
            sk = sin((1.d0+xpg(kpg))*thetaElbow/2.d0-tk(iNode))

            iBloc = (9+6*(nbFourier-1))*(iNode-1)

            do iFourier = 2, nbFourier
                iColon = iBloc+6+3*(iFourier-2)
                cosmte = cos(iFourier*(phi-omegaElbow))
                sinmte = sin(iFourier*(phi-omegaElbow))
                b(1, iColon+1) = dhk*cosmte/dena
                b(1, iColon+2) = hk*sinmte*cosfi*(1.d0+zeta/radiusPipeMoy)/denr
                b(1, iColon+3) = -zeta*d2hk*cosmte/(denr*dena)+ &
                                 hk*cosmte*sinfi/denr+ &
                                 zeta*iFourier*hk*sinmte*cosfi/(radiusPipeMoy*denr)

                b(2, iColon+1) = 0.d0
                b(2, iColon+2) = (iFourier/radiusLayer)*hk*cosmte*(1.d0+zeta/radiusPipeMoy)
                b(2, iColon+3) = (1.d0/radiusLayer)*hk*cosmte* &
                                 (1.d0+zeta*iFourier*iFourier/radiusPipeMoy)

                b(3, iColon+1) = -(iFourier/radiusLayer)*hk*sinmte- &
                                 hk*cosmte*cosfi/denr- &
                                 zeta*hk*cosfi*sinfi*cosmte/dena*dent+ &
                                 zeta*hk*(cosmte*cosfi-iFourier*sinmte*sinfi)/(radiusLayer*dena)
                b(3, iColon+2) = dhk*sinmte*(1+zeta/radiusPipeMoy)/denr
                b(3, iColon+3) = zeta*iFourier*dhk*sinmte*(1.d0/(radiusPipeMoy*denr)+ &
                                                           1.d0/(radiusLayer*dena))+ &
                                 zeta*dhk*cosfi*cosmte/dena*dent

                b(4, iColon+1) = 0.d0
                b(4, iColon+2) = 0.d0
                b(4, iColon+3) = 0.d0

            end do

            do iFourier = 2, nbFourier
                iColon = iBloc+6+3*(nbFourier-1)+3*(iFourier-2)
                cosmte = cos(iFourier*(phi-omegaElbow))
                sinmte = sin(iFourier*(phi-omegaElbow))

                b(1, iColon+1) = dhk*sinmte/dena
                b(1, iColon+2) = hk*cosmte*cosfi*(1.d0+zeta/radiusPipeMoy)/denr
                b(1, iColon+3) = -zeta*d2hk*sinmte/(denr*dena)+ &
                                 hk*sinmte*sinfi/denr- &
                                 zeta*iFourier*hk*cosmte*cosfi/(radiusPipeMoy*denr)

                b(2, iColon+1) = 0.d0
                b(2, iColon+2) = -(iFourier/radiusLayer)*hk*sinmte*(1.d0+zeta/radiusPipeMoy)
                b(2, iColon+3) = (1.d0/radiusLayer)*hk*sinmte* &
                                 (1.d0+zeta*iFourier*iFourier/radiusPipeMoy)

                b(3, iColon+1) = (iFourier/radiusLayer)*hk*cosmte-hk*sinmte*cosfi/denr- &
                                 zeta*hk*cosfi*sinfi*sinmte/dena*dent+ &
                                 zeta*hk*(sinmte*cosfi+iFourier*cosmte*sinfi)/(radiusLayer*dena)
                b(3, iColon+2) = dhk*cosmte*(1.d0+zeta/radiusPipeMoy)/denr
                b(3, iColon+3) = -zeta*iFourier*dhk*cosmte*(1.d0/(radiusPipeMoy*denr)+ &
                                                            1.d0/(radiusLayer*dena))+ &
                                 zeta*dhk*cosfi*sinmte/dena*dent

                b(4, iColon+1) = 0.d0
                b(4, iColon+2) = 0.d0
                b(4, iColon+3) = 0.d0

            end do

            iColon = iBloc+6*(nbFourier-1)+6
            b(1, iColon+1) = hk*sinfi/denr-zeta*d2hk/(dena*denr)
            b(2, iColon+1) = hk/radiusLayer
            b(3, iColon+1) = zeta*dhk*cosfi/dena*dent
            b(4, iColon+1) = 0.d0

            b(1, iColon+2) = -zeta*d2hk*coste/(denr*dena)+ &
                             2.d0*hk*cosfi* &
                             sinte/denr*zeta/radiusPipeMoy+hk*(coste*sinfi+cosfi*sinte)/denr
            b(2, iColon+2) = (2.d0/radiusLayer)*hk*coste*(1.d0+zeta/radiusPipeMoy)
            b(3, iColon+2) = dhk*sinte* &
                             ((1.d0+2.d0*zeta/radiusPipeMoy)/denr+zeta/(radiusLayer*dena))+ &
                             zeta*dhk*coste*cosfi/dena*dent
            b(4, iColon+2) = 0.d0
            !
            b(1, iColon+3) = -zeta*d2hk*sinte/(dena*denr)+ &
                             hk*(sinfi*sinte-cosfi*coste)/denr- &
                             2.d0*zeta*hk*cosfi*coste/(radiusPipeMoy*denr)
            b(2, iColon+3) = (2.d0/radiusLayer)*hk*sinte*(1.d0+zeta/radiusPipeMoy)
            b(3, iColon+3) = -dhk*coste* &
                             ((1.d0+2.d0*zeta/radiusPipeMoy)/denr+zeta/(radiusLayer*dena))+ &
                             zeta*dhk*sinte*cosfi/dena*dent
            b(4, iColon+3) = 0.d0
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeNMatr
!
! Get N matrix of shape functions and its transpose (for shell part)
!
! In  pipeElem         : properties of pipe element element
! In  nbDof            : number of DOF in element
! In  kpg              : current Gauss point
! In  phi              : coordinate of sub-point in section
! In  ff               : shape functions on beam
! IO  nvec             : vector of discrete disp/rota
! IO  tnvec              transposed vector of discrete disp/rota
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeNMatr(pipeElem, nbDof, kpg, &
                         phi, &
                         ff, nvec, tnvec)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbDof, kpg
        real(kind=8), intent(in) :: phi
        real(kind=8), intent(in) :: ff(*)
        real(kind=8), intent(inout) :: nvec(PIPE_NBDOF_BEAM, nbDof), tnvec(nbDof, PIPE_NBDOF_BEAM)
! ----- Local
        integer(kind=8) :: iNode, iBloc, iFourier, iColon, nbNode, nbFourier
        real(kind=8) :: hk, cosmfi, sinmfi, cosfi, sinfi
!   ------------------------------------------------------------------------------------------------
!
        nbNode = pipeElem%nbNode
        nbFourier = pipeElem%nbFourier
        cosfi = cos(phi)
        sinfi = sin(phi)
        do iNode = 1, nbNode
            hk = ff(nbNode*(kpg-1)+iNode)
            iBloc = (9+6*(nbFourier-1))*(iNode-1)

            do iFourier = 2, nbFourier
                iColon = iBloc+6+6*(iFourier-2)
                cosmfi = cos(iFourier*phi)
                sinmfi = sin(iFourier*phi)
                nvec(4, iColon+1) = hk*cosmfi
                nvec(4, iColon+4) = hk*sinmfi
                nvec(5, iColon+2) = hk*sinmfi
                nvec(5, iColon+5) = hk*cosmfi
                nvec(6, iColon+3) = hk*cosmfi
                nvec(6, iColon+6) = hk*sinmfi

                tnvec(iColon+1, 4) = hk*cosmfi
                tnvec(iColon+4, 4) = hk*sinmfi
                tnvec(iColon+2, 5) = hk*sinmfi
                tnvec(iColon+5, 5) = hk*cosmfi
                tnvec(iColon+3, 6) = hk*cosmfi
                tnvec(iColon+6, 6) = hk*sinmfi
            end do
            iColon = iBloc+6*(nbFourier-1)+6
            nvec(5, iColon+2) = hk*sinfi
            nvec(5, iColon+3) = -hk*cosfi
            nvec(6, iColon+1) = hk
            nvec(6, iColon+2) = hk*cosfi
            nvec(6, iColon+3) = hk*sinfi
!
            tnvec(iColon+2, 5) = hk*sinfi
            tnvec(iColon+3, 5) = -hk*cosfi
            tnvec(iColon+1, 6) = hk
            tnvec(iColon+2, 6) = hk*cosfi
            tnvec(iColon+3, 6) = hk*sinfi
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeBMatr
!
! Compute B matrix (u => epsi)
!
! In  pipeElem         : properties of pipe element
! In  nbDof            : number of DOF in element
! In  npg              : number of Gauss points
! In  xpg              : coordinates of Gauss points
! In  phi              : coordinate of sub-point in section
! In  zeta             : Coordinate of sub-point in thickness
! In  radiusLayer      : radius of current layer
! In  kpg              : index of Gauss point
! In  ff               : shape functions on beam
! In  df1              : first derivative of shape functions on beam
! In  df2              : second derivative of shape functions on beam
! Out b                : B matrix (u => epsi)
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeBMatr(pipeElem, nbDof, &
                         npg, xpg, &
                         phi, zeta, radiusLayer, &
                         kpg, ff, df1, df2, &
                         b_, bT_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbDof
        real(kind=8), intent(in) :: phi, zeta, radiusLayer
        integer(kind=8), intent(in) :: npg
        real(kind=8), intent(in) :: xpg(npg)
        integer(kind=8), intent(in) :: kpg
        real(kind=8), intent(in) :: ff(*), df1(*), df2(*)
        real(kind=8), optional, intent(out) :: b_(PIPE_TENS_SIZE, nbDof)
        real(kind=8), optional, intent(out) :: bT_(nbDof, PIPE_TENS_SIZE)
! ----- Locals
        integer(kind=8) :: iDof, iTens
        real(kind=8) :: b(PIPE_TENS_SIZE, nbDof)
!   ------------------------------------------------------------------------------------------------
!
        if (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
! --------- For straight pipe
            call pipeBMatrStraight(pipeElem, nbDof, &
                                   phi, zeta, radiusLayer, &
                                   kpg, ff, df1, df2, &
                                   b)

        else if (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
! --------- For elbow pipe
            call pipeBMatrCurved(pipeElem, nbDof, &
                                 npg, xpg, &
                                 phi, zeta, radiusLayer, &
                                 kpg, ff, df1, df2, &
                                 b)
        else
            ASSERT(ASTER_FALSE)
        end if

        if (present(b_)) then
            b_ = b
        end if

! ----- Transpose of B matrix
        if (present(bT_)) then
            bT_ = 0.d0
            do iTens = 1, PIPE_TENS_SIZE
                do iDof = 1, nbDof
                    bT_(iDof, iTens) = b(iTens, iDof)
                end do
            end do
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeBMatrStraight
!
! Compute B matrix (u => epsi) for straight pipe
!
! In  pipeElem         : properties of pipe element
! In  nbDof            : number of DOF in element
! In  phi              : coordinate of sub-point in section
! In  zeta             : Coordinate of sub-point in thickness
! In  radiusLayer      : radius of current layer
! In  kpg              : index of Gauss point
! In  ff               : shape functions on beam
! In  df1              : first derivative of shape functions on beam
! In  df2              : second derivative of shape functions on beam
! IO  b                : B matrix (u => epsi)
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeBMatrStraight(pipeElem, nbDof, &
                                 phi, zeta, radiusLayer, &
                                 kpg, ff, df1, df2, &
                                 b)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbDof
        real(kind=8), intent(in) :: phi, zeta, radiusLayer
        integer(kind=8), intent(in) :: kpg
        real(kind=8), intent(in) :: ff(*), df1(*), df2(*)
        real(kind=8), intent(inout) :: b(PIPE_TENS_SIZE, nbDof)
! ----- Local
        integer(kind=8) :: shiftDOF
!   ------------------------------------------------------------------------------------------------
!
        shiftDOF = 9+6*(pipeElem%nbFourier-1)

! ----- Compute B matrix (u => epsi) for straight pipe - Beam part
        call beamBMatrStraight(pipeElem%beamElem, nbDof, kpg, &
                               phi, radiusLayer, shiftDOF, &
                               ff, df1, b)

! ----- Compute B matrix (u => epsi) for straight pipe - Shell part
        call shellBMatrStraight(pipeElem, nbDof, kpg, &
                                phi, zeta, radiusLayer, &
                                ff, df1, df2, &
                                b)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeBMatrCurved
!
! Compute B matrix (u => epsi) for curved pipe (elbow)
!
! In  pipeElem         : properties of pipe element
! In  nbDof            : number of DOF in element
! In  npg              : number of Gauss points
! In  xpg              : coordinates of Gauss points
! In  phi              : coordinate of sub-point in section
! In  zeta             : Coordinate of sub-point in thickness
! In  radiusLayer      : radius of current layer
! In  kpg              : index of Gauss point
! In  ff               : shape functions on beam
! In  df1              : first derivative of shape functions on beam
! In  df2              : second derivative of shape functions on beam
! IO  b                : B matrix (u => epsi)
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeBMatrCurved(pipeElem, nbDof, &
                               npg, xpg, &
                               phi, zeta, radiusLayer, &
                               kpg, ff, df1, df2, &
                               b)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbDof
        real(kind=8), intent(in) :: phi, zeta, radiusLayer
        integer(kind=8), intent(in) :: npg
        real(kind=8), intent(in) :: xpg(npg)
        integer(kind=8), intent(in) :: kpg
        real(kind=8), intent(in) :: ff(*), df1(*), df2(*)
        real(kind=8), intent(inout) :: b(PIPE_TENS_SIZE, nbDof)
! ----- Local
        integer(kind=8) :: nbFourier, shiftDOF
!   ------------------------------------------------------------------------------------------------
!
        nbFourier = pipeElem%nbFourier
        shiftDOF = 9+6*(nbFourier-1)

! ----- Compute B matrix (u => epsi) for elbow - Beam part
        call beamBMatrCurved(pipeElem%beamElem, nbDof, &
                             npg, kpg, xpg, &
                             phi, radiusLayer, shiftDOF, &
                             ff, df1, b)

! ----- Compute B matrix (u => epsi) for elbow - Shell part
        call shellBMatrCurved(pipeElem, nbDof, &
                              kpg, npg, xpg, &
                              phi, zeta, radiusLayer, &
                              ff, df1, df2, &
                              b)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeLoadPres
!
! Compute load for internal pressure in pipe element
!
! In  typeScal         : type of scalar for load (real, complex or function)
! In  pipeElem         : properties of pipe element
! In  nbNode           : number of nodes of element
! In  nbFourier        : number of Fourier modes
! In  nbSect           : number of sectors in section of pipe
! In  npg              : number of Gauss points
! In  nspgSect         : number of sub-points in sectors of pipe
! In  jvLoad           : adress to data from load
! In  jvTime           : adress to time parameters
! In  jvGeom           : adress to initial coordinates of nodes
! In  iPoids           : adress to weight of Causs
! In  lAbsCurv         : flag for curvilinear coordinates
! In  absCurv          : value of curvilinear coordinates at nodes
! In  ff               : shape functions
! In  xpg              : coordinates of Gauss points
! In  weightSect       : weight of integration scheme for sectors
! In  nbDof            : number of DOF in element
! In  jvVect           : adress to load vector
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeLoadPres(typeScal, pipeElem, &
                            nbNode, nbFourier, nbSect, &
                            npg, nspgSect, &
                            jvLoad, jvTime, jvGeom, ipoids, &
                            lAbsCurv, absCurv, ff, &
                            xpg, weightSect, &
                            nbDof, jvVect)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: typeScal
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbNode, nbFourier, nbSect
        integer(kind=8), intent(in) :: npg, nspgSect
        integer(kind=8), intent(in) :: jvLoad, jvTime, jvGeom, ipoids
        aster_logical, intent(in) :: lAbsCurv
        real(kind=8), intent(in) :: absCurv(nbNode)
        real(kind=8), intent(in) :: ff(*)
        real(kind=8), intent(in) :: xpg(npg), weightSect(nspgSect)
        integer(kind=8), intent(in) :: nbDof
        integer(kind=8), intent(in) :: jvVect
! ----- Local
        integer(kind=8) :: iNode, kpg, kspgSect
        real(kind=8) :: hk, preskpg(npg), ck, sk, cosfi, sinfi
        real(kind=8) :: poids, phi, presNode(nbNode)
        integer(kind=8) :: indxDX, indxDY, indxDZ, indxWO, indxWI1, indxWO1
        real(kind=8) :: elemLength
        type(sectPipe_Prop) :: sectPipe
        type(beamElem_Prop) :: beamElem
!   ------------------------------------------------------------------------------------------------
!
        beamElem = pipeElem%beamElem
        sectPipe = beamElem%sectPipe

! ----- Evaluation of pressure at nodes
        call evalPresAtNodes(typeScal, &
                             nbNode, &
                             jvLoad, jvTime, jvGeom, &
                             lAbsCurv, absCurv, &
                             presNode)

! ----- Evaluation of pressure at Gauss points (interpolation)
        do kpg = 1, npg
            presKpg(kpg) = 0.d0
            do iNode = 1, nbNode
                hk = ff(nbNode*(kpg-1)+iNode)
                presKpg(kpg) = hk*presNode(iNode)+presKpg(kpg)
            end do
        end do

! ----- Evaluation of pressure on pipe
        do iNode = 1, nbNode
            indxDX = jvVect-1+(6+6*(nbFourier-1)+3)*(iNode-1)+1
            indxDY = jvVect-1+(6+6*(nbFourier-1)+3)*(iNode-1)+2
            indxDZ = jvVect-1+(6+6*(nbFourier-1)+3)*(iNode-1)+3
            indxWO = jvVect-1+(6+6*(nbFourier-1)+3)*(iNode-1)+1+6+6*(nbFourier-1)
            indxWI1 = jvVect-1+(6+6*(nbFourier-1)+3)*(iNode-1)+2+6+6*(nbFourier-1)
            indxWO1 = jvVect-1+(6+6*(nbFourier-1)+3)*(iNode-1)+3+6+6*(nbFourier-1)
            do kpg = 1, npg
                hk = ff(nbNode*(kpg-1)+iNode)
                ck = 1.d0
                sk = 0.d0
                if (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
                    ck = cos((1.d0+xpg(kpg))*beamElem%thetaElbow/2.d0-beamElem%tk(iNode))
                    sk = sin((1.d0+xpg(kpg))*beamElem%thetaElbow/2.d0-beamElem%tk(iNode))
                end if

                do kspgSect = 1, nspgSect
                    if (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
                        poids = zr(ipoids-1+kpg)* &
                                weightSect(kspgSect)*(beamElem%elemLength/2.d0)* &
                                2.d0*r8pi()/(2.d0*nbSect)*beamElem%sectPipe%radiusInt
                        zr(indxWO) = zr(indxWO)+hk*poids*presKpg(kpg)

                    elseif (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
                        phi = (kspgSect-1)*2.d0*r8pi()/(2.d0*nbSect)
                        cosfi = cos(phi)
                        sinfi = sin(phi)
                        elemLength = beamElem%thetaElbow*( &
                                     beamElem%radiusElbow+beamElem%sectPipe%radiusInt*sinfi)
                        poids = zr(ipoids-1+kpg)* &
                                weightSect(kspgSect)*(elemLength/2.d0)* &
                                2.d0*r8pi()/(2.d0*nbSect)*beamElem%sectPipe%radiusInt
                        zr(indxDX) = zr(indxDX)+hk*poids*presKpg(kpg)*sinfi*sk
                        zr(indxDY) = zr(indxDY)-hk*poids*presKpg(kpg)*sinfi*ck
                        zr(indxDZ) = zr(indxDZ)-hk*poids*presKpg(kpg)*cosfi
                        zr(indxWO) = zr(indxWO)+hk*poids*presKpg(kpg)
                        zr(indxWI1) = zr(indxWI1)+hk*poids*presKpg(kpg)*cos(phi-beamElem%omegaElbow)
                        zr(indxWO1) = zr(indxWO1)+hk*poids*presKpg(kpg)*sin(phi-beamElem%omegaElbow)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                end do
            end do
        end do

! ----- Change base
        call pipeBaseForVect('LG', pipeElem, nbDof, zr(jvVect))
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetDensity
!
! Get density for pipes
!
! In  jvMaterCode      : adress for material parameters
! Out rho              : density
! In  temp             : given temperature
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetDensity(jvMaterCode, rho, temp_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvMaterCode
        real(kind=8), optional, intent(in) :: temp_
        real(kind=8), intent(out) :: rho
! ----- Local
        integer(kind=8), parameter :: nbProp = 1, nbPara = 1
        character(len=8), parameter :: paraName = 'TEMP'
        character(len=16), parameter :: propName(nbProp) = (/'RHO'/)
        real(kind=8) :: propVale(nbProp)
        integer(kind=8) :: propError(nbProp)
        character(len=16) :: elasKeyword
        integer(kind=8), parameter:: kpg = 1, spt = 1
        character(len=8), parameter :: fami = 'FPG1', poum = '+'
!   ------------------------------------------------------------------------------------------------
!
        rho = 0.d0
        if (present(temp_)) then
            call rcvala(zi(jvMaterCode), ' ', 'ELAS', &
                        nbPara, paraName, [temp_], &
                        nbProp, propName, propVale, &
                        propError, 1)
            rho = propVale(1)
        else
            call rccoma(zi(jvMaterCode), 'ELAS', 1, elasKeyword)
            if (elasKeyword .eq. 'ELAS' .or. &
                elasKeyword .eq. 'ELAS_ISTR' .or. &
                elasKeyword .eq. 'ELAS_ORTH') then
                call rcvalb(fami, kpg, spt, poum, zi(jvMaterCode), &
                            ' ', elasKeyword, &
                            0, ' ', [0.d0], &
                            nbProp, propName, propVale, propError, &
                            1)
                rho = propVale(1)
            else
                call utmess('F', 'PIPE1_50')
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeLoadLine
!
! Compute load for lineic force
!
! In  typeScal         : type of scalar for load (real, complex or function)
! In  lGravity         : flag if lineic force is gravity
! In  pipeElem         : properties of pipe element
! In  nbNode           : number of nodes of element
! In  nbFourier        : number of Fourier modes
! In  npg              : number of Gauss points
! In  nspgSect         : number of sub-points in sectors of pipe
! In  nspgLayer        : number of sub-points in layers of pipe
! In  nbSect           : number of sectors in section of pipe
! In  nbLayer          : number of layers in section of pipe
! In  jvLoad           : adress to data from load
! In  jvTime           : adress to time parameters
! In  jvGeom           : adress to initial coordinates of nodes
! In  iPoids           : adress to weight of Causs
! In  rho              : density
! In  ff               : shape functions
! In  xpg              : coordinates of Gauss points
! In  weightSect       : weight of integration scheme for sectors
! In  weightLayer      : weight of integration scheme for layers
! In  nbDof            : number of DOF in element
! In  jvVect           : adress to load vector
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeLoadLine(typeScal, lGravity, pipeElem, &
                            nbNode, nbFourier, &
                            nbSect, nbLayer, &
                            npg, nspgSect, nspgLayer, &
                            jvLoad, jvTime, jvGeom, ipoids, &
                            rho, &
                            ff, weightSect, weightLayer, &
                            nbDof, jvVect)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: typeScal
        aster_logical, intent(in) :: lGravity
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbNode, nbFourier
        integer(kind=8), intent(in) :: nbSect, nbLayer
        integer(kind=8), intent(in) :: npg, nspgSect, nspgLayer
        integer(kind=8), intent(in) :: jvLoad, jvTime, jvGeom, ipoids
        real(kind=8), intent(in) :: rho
        real(kind=8), intent(in) :: ff(*), weightSect(nspgSect), weightLayer(nspgLayer)
        integer(kind=8), intent(in) :: nbDof, jvVect
! ----- Local
        aster_logical :: lWind, lGlob, lCplxRealPart
        integer(kind=8) :: iNode, iComp, nbComp, kpg, indxDof, kspgSect, kspgLayer, iDof
        real(kind=8) :: hk, phi, radiusLayer
        real(kind=8) :: jacobi, poids
        real(kind=8) :: lineLoadGlob(6), lineLoadLoca(6)
        real(kind=8) :: loadVect(nbDof)
        type(sectPipe_Prop) :: sectPipe
        type(beamElem_Prop) :: beamElem
!   ------------------------------------------------------------------------------------------------
!
        beamElem = pipeElem%beamElem
        sectPipe = beamElem%sectPipe

! ----- Checks
        if (.not. lGravity) then
            if (typeScal .eq. 'R') then
                lGlob = zr(jvLoad+6) .eq. 0.d0
                if (.not. lGlob) then
                    call utmess('F', 'PIPE1_45')
                end if
            end if
            if (typeScal .eq. 'F') then
                lWind = zk8(jvLoad+6) .eq. 'VENT'
                lGlob = zk8(jvLoad+6) .eq. 'GLOBAL'
                if (lWind) then
                    call utmess('F', 'PIPE1_44')
                end if
                if (.not. lGlob) then
                    call utmess('F', 'PIPE1_45')
                end if
            end if
        end if

! ----- Evaluation of  components of lineic force in global base
        if (typeScal .eq. 'R') then
            lCplxRealPart = ASTER_TRUE
            call evalLineLoad(typeScal, lCplxRealPart, lGravity, &
                              jvLoad, jvTime, jvGeom, &
                              rho, lineLoadGlob)
            if (.not. lGravity) then
                lineLoadGlob = lineLoadGlob/sectPipe%area
            end if
        end if

! ----- Complex algebra case
        nbComp = 1
        if (typeScal .eq. 'C') then
            nbComp = 2
        end if

! ----- Compute load
        loadVect = 0.d0
        do iComp = 1, nbComp
            lCplxRealPart = iComp .eq. 1

            do iNode = 1, nbNode
! ------------- Evaluation of  components of lineic force in global base
                if (typeScal .ne. 'R') then
                    call evalLineLoad(typeScal, lCplxRealPart, lGravity, &
                                      jvLoad, jvTime, jvGeom, &
                                      rho, lineLoadGlob, iNode)
! ----------------- By sector area
                    if (.not. lGravity) then
                        lineLoadGlob = lineLoadGlob/sectPipe%area
                    end if
                end if

! ------------- Change base
                lineLoadLoca = 0.d0
                if (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
                    call utpvgl(1, 6, beamElem%pglCell, lineLoadGlob, lineLoadLoca)
                else
                    if (iNode .eq. 1) then
                        call utpvgl(1, 6, beamElem%pgl(:, :, 1), lineLoadGlob, lineLoadLoca)
                    elseif (iNode .eq. 2) then
                        call utpvgl(1, 6, beamElem%pgl(:, :, 2), lineLoadGlob, lineLoadLoca)
                    elseif (iNode .eq. 3) then
                        call utpvgl(1, 6, beamElem%pgl(:, :, 3), lineLoadGlob, lineLoadLoca)
                    elseif (iNode .eq. 4) then
                        call utpvgl(1, 6, beamElem%pgl(:, :, 4), lineLoadGlob, lineLoadLoca)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                end if

! ------------- Compute load
                indxDof = (9+6*(nbFourier-1))*(iNode-1)
                do kpg = 1, npg
                    hk = ff(nbNode*(kpg-1)+iNode)
                    do kspgLayer = 1, nspgLayer
! --------------------- Radius of layer
                        call pipeGetRadiusLayer(pipeElem, kspgLayer, nbLayer, radiusLayer)

                        do kspgSect = 1, nspgSect
! ------------------------- Coordinate of sub-point in section
                            phi = (kspgSect-1)*2.d0*r8pi()/(2.d0*nbSect)

! ------------------------- Compute jacobian
                            call pipeGetJacobian(pipeElem, nbSect, nbLayer, phi, radiusLayer, &
                                                 jacobi)

! ------------------------- Compute integral
                            poids = zr(ipoids-1+kpg)* &
                                    weightLayer(kspgLayer)*weightSect(kspgSect)* &
                                    jacobi

! ------------------------- Compute vector
                            loadVect(indxDof+1) = loadVect(indxDof+1)+poids*hk*lineLoadLoca(1)
                            loadVect(indxDof+2) = loadVect(indxDof+2)+poids*hk*lineLoadLoca(2)
                            loadVect(indxDof+3) = loadVect(indxDof+3)+poids*hk*lineLoadLoca(3)

                        end do
                    end do
                end do
            end do

! --------- Change to global base
            call pipeBaseForVect('LG', pipeElem, nbDof, loadVect)

! --------- Save values in vector
            if (typeScal .eq. 'R' .or. typeScal .eq. 'F') then
                ASSERT(iComp .eq. 1)
                ASSERT(nbComp .eq. 1)
                do iDof = 1, nbDof
                    zr(jvVect-1+iDof) = loadVect(iDof)
                end do
            else
                if (lCplxRealPart) then
                    do iDof = 1, nbDof
                        zc(jvVect-1+iDof) = loadVect(iDof)
                    end do
                else
                    do iDof = 1, nbDof
                        zc(jvVect-1+iDof) = dcmplx(dble(zc(jvVect-1+iDof)), dble(loadVect(iDof)))
                    end do
                end if
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetProperties
!
! Get properties of pipe
!
! In  nbNode           : number of nodes of element
! In  nbFourier        : number of Fourier modes
! Out pipeElem         : properties of pipe element
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetProperties(nbNode, nbFourier, pipeElem)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nbNode, nbFourier
        type(pipeElem_Prop), intent(out) :: pipeElem
! ----- Local
        integer(kind=8) :: pipeType, jvPipe
        aster_logical :: lModiMetric
        real(kind=8) :: thetaElbow, radiusElbow, omegaElbow, pipeLength
        real(kind=8) :: pgl(3, 3)
        real(kind=8) :: pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), pgl4(3, 3)
        type(sectPipe_Prop) :: sectPipe
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(nbNode .le. PIPE_MAX_NODE)
        call jevech('PCAORIE', 'L', jvPipe)

! ----- General properties
        pipeElem%nbNode = nbNode
        pipeElem%nbFourier = nbFourier

! ----- Get geometric properties of section
        call pipeGetSection(sectPipe)
        pipeElem%beamElem%sectType = BEAM_SECT_PIPE
        pipeElem%beamElem%sectPipe = sectPipe

! ----- Get type of pipe
        call pipeGetType(jvPipe, nbNode, pipeType, lModiMetric)
        pipeElem%pipeType = pipeType
        pipeElem%lModiMetric = lModiMetric

! ----- Get geometry of pipe
        call pipeGetGeometry(jvPipe, nbNode, &
                             thetaElbow, radiusElbow, omegaElbow, pipeLength)

! ----- Set properties of beam
        pipeElem%beamElem%nbNode = nbNode
        pipeElem%beamElem%elemLength = pipeLength
        if (pipeType .eq. PIPE_TYPE_STRAIGHT) then
            pipeElem%beamElem%beamType = BEAM_TYPE_STRAIGHT
        elseif (pipeType .eq. PIPE_TYPE_ELBOW) then
            pipeElem%beamElem%beamType = BEAM_TYPE_ELBOW
            pipeElem%beamElem%thetaElbow = thetaElbow
            pipeElem%beamElem%radiusElbow = radiusElbow
            pipeElem%beamElem%omegaElbow = omegaElbow
            if (nbNode .eq. 3) then
                pipeElem%beamElem%tk(1) = 0.d0
                pipeElem%beamElem%tk(2) = thetaElbow
                pipeElem%beamElem%tk(3) = thetaElbow/2.d0
            else if (nbNode .eq. 4) then
                pipeElem%beamElem%tk(1) = 0.d0
                pipeElem%beamElem%tk(2) = thetaElbow
                pipeElem%beamElem%tk(3) = thetaElbow/3.d0
                pipeElem%beamElem%tk(4) = 2.d0*thetaElbow/3.d0
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Get basis at nodes of pipe
        call pipeGetBases(jvPipe, nbNode, pipeType, &
                          pgl, pgl1, pgl2, pgl3, pgl4)
        pipeElem%beamElem%pglCell = pgl
        pipeElem%beamElem%pgl(:, :, 1) = pgl1
        pipeElem%beamElem%pgl(:, :, 2) = pgl2
        pipeElem%beamElem%pgl(:, :, 3) = pgl3
        pipeElem%beamElem%pgl(:, :, 4) = pgl4
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetRadiusLayer
!
! Get radius of current layer (with or without MODI_METRIQUE)
!
! In  pipeElem         : properties of pipe element
! In  kspgLayer        : current sub-point in layer of section of pipe
! In  nbLayer          : number of layers in section of pipe
! Out radiusLayer      : radius of current layer (with or without MODI_METRIQUE)
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetRadiusLayer(pipeElem, kspgLayer, nbLayer, radiusLayer)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: kspgLayer, nbLayer
        real(kind=8), intent(out) :: radiusLayer
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(pipeElem%beamElem%sectType .eq. BEAM_SECT_PIPE)
        radiusLayer = pipeElem%beamElem%sectPipe%radiusMoy
        if (pipeElem%lModiMetric) then
            radiusLayer = radiusLayer+ &
                          (kspgLayer-1)*pipeElem%beamElem%sectPipe%thickness/(2.d0*nbLayer)- &
                          pipeElem%beamElem%sectPipe%thickness/2.d0
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeBaseForMatr
!
! Change base for matrix from local element base to global base
! Only for beam dof.
!
! In  pipeElem         : properties of pipe element
! In  nbDof            : number of DOF in element
! IO  matr             : matrix to project
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeBaseForMatr(pipeElem, nbDof, matr)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbDof
        real(kind=8), intent(inout) :: matr(nbDof, nbDof)
! ----- Locals
        integer(kind=8) :: iDof, jDof, lDof, iNode, nbNode
        real(kind=8) :: p(nbDof, nbDof), ktemp(nbDof, nbDof)
!   ------------------------------------------------------------------------------------------------
!
        nbNode = pipeElem%nbNode
        p = 0.d0
        do iDof = 1, nbDof
            do jDof = 1, nbDof
                if (iDof .eq. jDof) then
                    p(iDof, jDof) = 1.d0
                end if
            end do
        end do

        if (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
            do iNode = 1, nbNode
                do iDof = 1, 3
                    do jDof = 1, 3
                        p((iNode-1)*nbDof/nbNode+iDof, (iNode-1)*nbDof/nbNode+jDof) = &
                            pipeElem%beamElem%pglCell(iDof, jDof)
                        p((iNode-1)*nbDof/nbNode+3+iDof, (iNode-1)*nbDof/nbNode+3+jDof) = &
                            pipeElem%beamElem%pglCell(iDof, jDof)
                    end do
                end do
            end do
        else if (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
            do iDof = 1, 3
                do jDof = 1, 3
                    p(iDof, jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 1)
                    p(3+iDof, 3+jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 1)
                end do
            end do
            do iDof = 1, 3
                do jDof = 1, 3
                    p(nbDof/nbNode+iDof, nbDof/nbNode+jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 2)
                    p(nbDof/nbNode+3+iDof, nbDof/nbNode+3+jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 2)
                end do
            end do
            do iDof = 1, 3
                do jDof = 1, 3
                    p(2*nbDof/nbNode+iDof, 2*nbDof/nbNode+jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 3)
                    p(2*nbDof/nbNode+3+iDof, 2*nbDof/nbNode+3+jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 3)
                end do
            end do
            if (nbNode .eq. 4) then
                do iDof = 1, 3
                    do jDof = 1, 3
                        p(3*nbDof/nbNode+iDof, 3*nbDof/nbNode+jDof) = &
                            pipeElem%beamElem%pgl(iDof, jDof, 4)
                        p(3*nbDof/nbNode+3+iDof, 3*nbDof/nbNode+3+jDof) = &
                            pipeElem%beamElem%pgl(iDof, jDof, 4)
                    end do
                end do
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- [KTEMP] = [P]^T * [K]
        ktemp = 0.d0
        do iDof = 1, nbDof
            do jDof = 1, nbDof
                do lDof = 1, nbDof
                    ktemp(iDof, jDof) = ktemp(iDof, jDof)+p(lDof, iDof)*matr(lDof, jDof)
                end do
            end do
        end do

! ----- K = [P]^T * [K] * [P] = [KTEMP] * [P]
        matr = 0.d0
        do iDof = 1, nbDof
            do jDof = 1, nbDof
                do lDof = 1, nbDof
                    matr(iDof, jDof) = matr(iDof, jDof)+ktemp(iDof, lDof)*p(lDof, jDof)
                end do
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeBaseForVect
!
! Change base for vector from local element base to global base
! Only for beam dof.
!
! In  code             : LG => local to global
!                        GL => global to local
! In  pipeElem         : properties of pipe element
! In  nbDof            : number of DOF in element
! IO  vect             : vector
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeBaseForVect(code, pipeElem, nbDof, vect)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=2), intent(in) :: code
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbDof
        real(kind=8), intent(inout) :: vect(nbDof)
! ----- Locals
        integer(kind=8) :: iDof, jDof, iNode, nbNode
        real(kind=8) :: p(nbDof, nbDof), vtemp(nbDof)
!   ------------------------------------------------------------------------------------------------
!
        nbNode = pipeElem%nbNode
        p = 0.d0
        do iDof = 1, nbDof
            do jDof = 1, nbDof
                if (iDof .eq. jDof) then
                    p(iDof, jDof) = 1.d0
                end if
            end do
        end do

        if (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
            do iNode = 1, nbNode
                do iDof = 1, 3
                    do jDof = 1, 3
                        p((iNode-1)*nbDof/nbNode+iDof, (iNode-1)*nbDof/nbNode+jDof) = &
                            pipeElem%beamElem%pglCell(iDof, jDof)
                        p((iNode-1)*nbDof/nbNode+3+iDof, (iNode-1)*nbDof/nbNode+3+jDof) = &
                            pipeElem%beamElem%pglCell(iDof, jDof)
                    end do
                end do
            end do
        else if (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
            do iDof = 1, 3
                do jDof = 1, 3
                    p(iDof, jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 1)
                    p(3+iDof, 3+jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 1)
                end do
            end do
            do iDof = 1, 3
                do jDof = 1, 3
                    p(nbDof/nbNode+iDof, nbDof/nbNode+jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 2)
                    p(nbDof/nbNode+3+iDof, nbDof/nbNode+3+jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 2)
                end do
            end do
            do iDof = 1, 3
                do jDof = 1, 3
                    p(2*nbDof/nbNode+iDof, 2*nbDof/nbNode+jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 3)
                    p(2*nbDof/nbNode+3+iDof, 2*nbDof/nbNode+3+jDof) = &
                        pipeElem%beamElem%pgl(iDof, jDof, 3)
                end do
            end do
            if (nbNode .eq. 4) then
                do iDof = 1, 3
                    do jDof = 1, 3
                        p(3*nbDof/nbNode+iDof, 3*nbDof/nbNode+jDof) = &
                            pipeElem%beamElem%pgl(iDof, jDof, 4)
                        p(3*nbDof/nbNode+3+iDof, 3*nbDof/nbNode+3+jDof) = &
                            pipeElem%beamElem%pgl(iDof, jDof, 4)
                    end do
                end do
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

        vtemp = 0.d0
        if (code .eq. 'LG') then
! -------- {VTEMP} = [P]^T * {V}
            do iDof = 1, nbDof
                do jDof = 1, nbDof
                    vtemp(iDof) = vtemp(iDof)+p(jDof, iDof)*vect(jDof)
                end do
            end do
        else if (code .eq. 'GL') then
! -------- {VTEMP} = [P] * {V}
            do iDof = 1, nbDof
                do jDof = 1, nbDof
                    vtemp(iDof) = vtemp(iDof)+p(iDof, jDof)*vect(jDof)
                end do
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
        vect(1:nbDof) = vtemp(1:nbDof)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetJacobian
!
! Compute jacobian for current pipe element
!
! In  pipeElem         : properties of pipe element
! In  nbLayer          : number of layers in section of pipe
! In  nbSect           : number of sectors in section of pipe
! In  phi              : coordinate of sub-point in section
! In  radiusLayer      : radius of current layer
! Out jacobi           : jacobian
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetJacobian(pipeElem, nbSect, nbLayer, phi, radiusLayer, &
                               jacobi)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbSect, nbLayer
        real(kind=8), intent(in) :: phi, radiusLayer
        real(kind=8), intent(out) :: jacobi
!   ------------------------------------------------------------------------------------------------
!
        jacobi = 0.d0
        if (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
            jacobi = (pipeElem%beamElem%elemLength/2.d0)* &
                     pipeElem%beamElem%sectPipe%thickness* &
                     2.d0*r8pi()/(4.d0*nbLayer*nbSect)*radiusLayer

        elseif (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
            jacobi = (pipeElem%beamElem%thetaElbow* &
                      (pipeElem%beamElem%radiusElbow+radiusLayer*sin(phi))/2.d0)* &
                     pipeElem%beamElem%sectPipe%thickness* &
                     2.d0*r8pi()/(4.d0*nbLayer*nbSect)*radiusLayer
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetDisp
!
! Get displacements (in local base)
!
! In  pipeElem         : properties of pipe element
! In  nbDof            : number of DOF in element
! In  dispGlob         : displacements in global base
! Out dispLoca         : displacements in local base
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetDisp(pipeElem, nbDof, dispGlob, &
                           dispLoca)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: nbDof
        real(kind=8), intent(in) :: dispGlob(nbDof)
        real(kind=8), intent(out) :: dispLoca(nbDof)
! ----- Local
        integer(kind=8) :: nbNode
        real(kind=8) :: dispTemp(nbDof)
!   ------------------------------------------------------------------------------------------------
!
        dispLoca = 0.d0
        nbNode = pipeElem%nbNode

! ----- Move displacements in local base
        dispTemp = dispGlob
        call pipeBaseForVect('GL', pipeElem, nbDof, dispTemp)
        dispLoca = dispTemp
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeForcNoda
!
! Compute vector of nodal forces
!
! In  pipeElem         : properties of pipe element
! In  npg              : number of Gauss points
! In  xpg              : coordinates of Gauss points
! In  weightPg         : weight of Gauss points
! In  nbSect           : number of sectors in section of pipe
! In  nbLayer          : number of layers in section of pipe
! In  nspgSect         : number of sub-points in sectors of pipe
! In  nspgLayer        : number of sub-points in layers of pipe
! In  weightSect       : weight of integration scheme for sectors
! In  weightLayer      : weight of integration scheme for layers
! In  ff               : shape functions on beam
! In  df1              : first derivative of shape functions on beam
! In  df2              : second derivative of shape functions on beam
! In  nbDof            : number of DOF in element
! Out forcNoda         : vector of nodal forces
! In  jvSigm           : address to variable stress
! In  sigmRefe         : constant stress
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeForcNoda(pipeElem, &
                            npg, xpg, weightPg, &
                            nbSect, nbLayer, &
                            nspgSect, nspgLayer, &
                            weightSect, weightLayer, &
                            ff, df1, df2, &
                            nbDof, forcNoda, &
                            jvSigm_, sigmRefe_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: npg
        real(kind=8), intent(in) :: xpg(npg), weightPg(npg)
        integer(kind=8), intent(in) :: nbLayer, nbSect
        integer(kind=8), intent(in) :: nspgSect, nspgLayer
        real(kind=8), intent(in) :: weightSect(nspgSect), weightLayer(nspgLayer)
        real(kind=8), intent(in) :: ff(*), df1(*), df2(*)
        integer(kind=8), intent(in) :: nbDof
        real(kind=8), intent(out) :: forcNoda(nbDof)
        integer(kind=8), optional, intent(in) :: jvSigm_
        real(kind=8), optional, intent(in) :: sigmRefe_
! ----- Local
        type(sectPipe_Prop) :: sectPipe
        type(beamElem_Prop) :: beamElem
        integer(kind=8) :: iret, iDof, iRefe
        integer(kind=8) :: kpg, kspg, kspgLayer, kspgSect
        real(kind=8) :: phi, zeta
        real(kind=8) :: radiusLayer, sigmTemp(PIPE_TENS_SIZE), sigm(PIPE_TENS_SIZE)
        real(kind=8) :: bT(nbDof, PIPE_TENS_SIZE), prodMatrVect(nbDof)
        real(kind=8) :: jacobi, poids
!   ------------------------------------------------------------------------------------------------
!
        forcNoda = 0.d0
        beamElem = pipeElem%beamElem
        sectPipe = beamElem%sectPipe

! ----- Loop on Gauss points (on segment)
        kspg = 0
        do kpg = 1, npg
            bT = 0.d0
! --------- Loop on sub-points in thickness
            do kspgLayer = 1, nspgLayer
! ------------- Radius of pipe
                call pipeGetRadiusLayer(pipeElem, kspgLayer, nbLayer, radiusLayer)

! ------------- Coordinate of sub-point in thickness
                zeta = (kspgLayer-1)*pipeElem%beamElem%sectPipe%thickness/(2.d0*nbLayer)- &
                       pipeElem%beamElem%sectPipe%thickness/2.d0

! ------------- Loop on sub-points in section
                do kspgSect = 1, nspgSect
                    kspg = kspg+1

! ----------------- Coordinate of sub-point in section
                    phi = (kspgSect-1)*2.d0*r8pi()/(2.d0*nbSect)

! ----------------- Compute B^T matrix (u => epsi)
                    call pipeBMatr(pipeElem, nbDof, &
                                   npg, xpg, &
                                   phi, zeta, radiusLayer, &
                                   kpg, ff, df1, df2, &
                                   bT_=bT)

! ----------------- Compute jacobian
                    call pipeGetJacobian(pipeElem, nbSect, nbLayer, phi, radiusLayer, &
                                         jacobi)

! ----------------- Compute integral
                    poids = weightPg(kpg)* &
                            weightLayer(kspgLayer)*weightSect(kspgSect)* &
                            jacobi

! ----------------- Compute nodal force
                    if (present(jvSigm_)) then
                        sigm(1) = zr(jvSigm_-1+6*(kspg-1)+1)
                        sigm(2) = zr(jvSigm_-1+6*(kspg-1)+2)
                        sigm(3) = zr(jvSigm_-1+6*(kspg-1)+4)
                        sigm(4) = zr(jvSigm_-1+6*(kspg-1)+5)
                        call prmave(0, bT, nbDof, nbDof, PIPE_TENS_SIZE, &
                                    sigm, PIPE_TENS_SIZE, prodMatrVect, nbDof, &
                                    iret)
                        ASSERT(iret .eq. 0)
                        do iDof = 1, nbDof
                            forcNoda(iDof) = forcNoda(iDof)+prodMatrVect(iDof)*poids
                        end do
                    elseif (present(sigmRefe_)) then
                        sigmTemp = 0.d0
                        do iRefe = 1, 4
                            sigmTemp(iRefe) = sigmRefe_
                            call prmave(0, bT, nbDof, nbDof, PIPE_TENS_SIZE, &
                                        sigmTemp, PIPE_TENS_SIZE, prodMatrVect, nbDof, &
                                        iret)
                            ASSERT(iret .eq. 0)
                            sigmTemp(iRefe) = 0.d0
                            do iDof = 1, nbDof
                                forcNoda(iDof) = forcNoda(iDof)+abs(prodMatrVect(iDof)*poids)
                            end do
                        end do
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                end do
            end do
        end do

! ----- Change base
        call pipeBaseForVect('LG', pipeElem, nbDof, forcNoda)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeEfgeElno
!
! Compute EFGE_ELNO from EFGE_ELGA
!
! In  pipeElem         : properties of pipe element
! In  jgano            : adress to Gauss to point matrix
! In  nbNode           : number of nodes of element
! In  npg              : number of Gauss points
! In  ff               : shape functions on beam
! In  efgeElga         : generalized forces at Gauss points (ELGA)
! Out efgeElno         : generalized forces at nodes in cell (ELNO)
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeEfgeElno(pipeElem, &
                            jgano, nbNode, npg, xpg, ff, &
                            efgeElga, efgeElno)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: jgano, nbNode, npg
        real(kind=8), intent(in) :: xpg(npg), ff(*)
        real(kind=8), intent(in) :: efgeElga(PIPE_TENS_SIZE, PIPE_NBDOF_BEAM)
        real(kind=8), intent(out) :: efgeElno(PIPE_NBDOF_BEAM*nbNode)
! ----- Local
        type(sectPipe_Prop) :: sectPipe
        type(beamElem_Prop) :: beamElem
        integer(kind=8) :: kpg, iNode, iDofBeam
        integer(kind=8) :: ih, ip, i1, i2
        real(kind=8) :: hk(PIPE_MAX_NODE, PIPE_MAX_NPG)
        real(kind=8) :: co(PIPE_MAX_NPG, PIPE_MAX_NODE), si(PIPE_MAX_NPG, PIPE_MAX_NODE)
        real(kind=8) :: cp(2, 2), cv(2, 2)
        real(kind=8) :: fno(PIPE_NBDOF_BEAM)
        real(kind=8) :: alphaf, betaf, alpham, betam
        real(kind=8) :: xa, xb, xc, xd
        real(kind=8) :: efge(PIPE_NBDOF_BEAM), valePg(PIPE_MAX_NPG)

!   ------------------------------------------------------------------------------------------------
!
        efgeElno = 0.d0
        beamElem = pipeElem%beamElem
        sectPipe = beamElem%sectPipe

! ----- Matrix of shape functions for each node
        hk = 0.d0
        do iNode = 1, nbNode
            do kpg = 1, npg
                hk(iNode, kpg) = ff(nbNode*(kpg-1)+iNode)
            end do
        end do

        if ((nbNode .eq. 3) .and. (npg .eq. 3)) then
! !      POUR NE PAS SUPPRIMER LA SAVANTE PROGRAMMATION DE PATRICK
            do kpg = 1, npg
                do iNode = 1, nbNode
                    if (pipeElem%pipeType .eq. PIPE_TYPE_STRAIGHT) then
                        co(kpg, iNode) = 1.d0
                        si(kpg, iNode) = 0.d0
                    elseif (pipeElem%pipeType .eq. PIPE_TYPE_ELBOW) then
                        co(kpg, iNode) = &
                            cos((1.d0+xpg(kpg))*beamElem%thetaElbow/2.d0-beamElem%tk(iNode))
                        si(kpg, iNode) = &
                            sin((1.d0+xpg(kpg))*beamElem%thetaElbow/2.d0-beamElem%tk(iNode))
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                end do
            end do
            do iNode = 1, nbNode
                if (iNode .eq. 1) then
                    ih = 2
                    ip = 1
                    i1 = 1
                    i2 = 3
                else if (iNode .eq. 2) then
                    ih = 1
                    ip = 2
                    i1 = 3
                    i2 = 1
                else
                    do iDofBeam = 1, PIPE_NBDOF_BEAM
                        fno(iDofBeam) = efgeElga(2, iDofBeam)
                    end do
                    goto 380
                end if
                cp(1, 1) = co(1, ih)*co(1, 3)+si(1, ih)*si(1, 3)
                cp(1, 2) = -co(1, ih)*si(1, 3)+si(1, ih)*co(1, 3)
                cp(2, 1) = -cp(1, 2)
                cp(2, 2) = cp(1, 1)
                cv(1, 1) = co(3, ih)*co(3, 3)+si(3, ih)*si(3, 3)
                cv(1, 2) = -co(3, ih)*si(3, 3)+si(3, ih)*co(3, 3)
                cv(2, 1) = -cp(1, 2)
                cv(2, 2) = cp(1, 1)
                alphaf = hk(ih, 3)*(co(1, ih)*efgeElga(1, 1)+ &
                                    si(1, ih)*efgeElga(1, 2))- &
                         hk(ih, 3)*hk(3, 1)*(cp(1, 1)*efgeElga(2, 1)+ &
                                             cp(1, 2)*efgeElga(2, 2))- &
                         hk(ih, 1)*(co(3, ih)*efgeElga(3, 1)+ &
                                    si(3, ih)*efgeElga(3, 2))+ &
                         hk(ih, 1)*hk(3, 3)*(cv(1, 1)*efgeElga(2, 1)+ &
                                             cv(1, 2)*efgeElga(2, 2))
                betaf = hk(ih, 3)*(-si(1, ih)*efgeElga(1, 1)+ &
                                   co(1, ih)*efgeElga(1, 2))- &
                        hk(ih, 3)*hk(3, 1)*(cp(2, 1)*efgeElga(2, 1)+ &
                                            cp(2, 2)*efgeElga(2, 2))- &
                        hk(ih, 1)*(-si(3, ih)*efgeElga(3, 1)+ &
                                   co(3, ih)*efgeElga(3, 2))+ &
                        hk(ih, 1)*hk(3, 3)*(cv(2, 1)*efgeElga(2, 1)+ &
                                            cv(2, 2)*efgeElga(2, 2))
                alpham = hk(ih, 3)*(co(1, ih)*efgeElga(1, 4)+ &
                                    si(1, ih)*efgeElga(1, 5))- &
                         hk(ih, 3)*hk(3, 1)*(cp(1, 1)*efgeElga(2, 4)+ &
                                             cp(1, 2)*efgeElga(2, 5))- &
                         hk(ih, 1)*(co(3, ih)*efgeElga(3, 4)+ &
                                    si(3, ih)*efgeElga(3, 5))+ &
                         hk(ih, 1)*hk(3, 3)*(cv(1, 1)*efgeElga(2, 4)+ &
                                             cv(1, 2)*efgeElga(2, 5))
                betam = hk(ih, 3)*(-si(1, ih)*efgeElga(1, 4)+ &
                                   co(1, ih)*efgeElga(1, 5))- &
                        hk(ih, 3)*hk(3, 1)*(cp(2, 1)*efgeElga(2, 4)+ &
                                            cp(2, 2)*efgeElga(2, 5))- &
                        hk(ih, 1)*(-si(3, ih)*efgeElga(3, 4)+ &
                                   co(3, ih)*efgeElga(3, 5))+ &
                        hk(ih, 1)*hk(3, 3)*(cv(2, 1)*efgeElga(2, 4)+ &
                                            cv(2, 2)*efgeElga(2, 5))
                cp(1, 1) = co(1, ih)*co(1, ip)+si(1, ih)*si(1, ip)
                cp(1, 2) = -co(1, ih)*si(1, ip)+si(1, ih)*co(1, ip)
                cp(2, 1) = -cp(1, 2)
                cp(2, 2) = cp(1, 1)
                cv(1, 1) = co(3, ih)*co(3, ip)+si(3, ih)*si(3, ip)
                cv(1, 2) = -co(3, ih)*si(3, ip)+si(3, ih)*co(3, ip)
                cv(2, 1) = -cp(1, 2)
                cv(2, 2) = cp(1, 1)
                xa = hk(ip, 1)*hk(ih, 3)*cp(1, 1)-hk(ip, 3)*hk(ih, 1)*cv(1, 1)
                xb = hk(ip, 1)*hk(ih, 3)*cp(1, 2)-hk(ip, 3)*hk(ih, 1)*cv(1, 2)
                xc = hk(ip, 1)*hk(ih, 3)*cp(2, 1)-hk(ip, 3)*hk(ih, 1)*cv(2, 1)
                xd = hk(ip, 1)*hk(ih, 3)*cp(2, 2)-hk(ip, 3)*hk(ih, 1)*cv(2, 2)
                efge = 0.d0
                efge(1) = (xd*alphaf-xb*betaf)/(xa*xd-xb*xc)
                efge(2) = (-xc*alphaf+xa*betaf)/(xa*xd-xb*xc)
                efge(3) = (hk(ih, i2)*efgeElga(i1, 3)- &
                           hk(ih, i1)*efgeElga(i2, 3)- &
                           efgeElga(2, 3)*(hk(3, i1)*hk(ih, i2)- &
                                           hk(3, i2)*hk(ih, i1)))/(hk(1, 1)*hk(2, 3)- &
                                                                   hk(1, 3)*hk(2, 1))
                efge(4) = (xd*alpham-xb*betam)/(xa*xd-xb*xc)
                efge(5) = (-xc*alpham+xa*betam)/(xa*xd-xb*xc)
                efge(6) = (hk(ih, i2)*efgeElga(i1, 6)- &
                           hk(ih, i1)*efgeElga(i2, 6)- &
                           efgeElga(2, 6)*(hk(3, i1)*hk(ih, i2)- &
                                           hk(3, i2)*hk(ih, i1)))/(hk(1, 1)*hk(2, 3)- &
                                                                   hk(1, 3)*hk(2, 1))
380             continue
                do iDofBeam = 1, PIPE_NBDOF_BEAM
                    efgeElno(PIPE_NBDOF_BEAM*(iNode-1)+iDofBeam) = efge(iDofBeam)
                end do
            end do
        else
            do iDofBeam = 1, PIPE_NBDOF_BEAM
                valePg = 0.d0
                valePg(1:npg) = efgeElga(1:npg, iDofBeam)
                efge = 0.d0
                call ppgan2(jgano, 1, 1, valePg, efge)
                do iNode = 1, nbNode
                    efgeElno(PIPE_NBDOF_BEAM*(iNode-1)+iDofBeam) = efge(iNode)
                end do
            end do
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeCheckMetric
!
! Check properties of pipe for modiMetric
!
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeCheckMetric(nomteZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: nomteZ
! ----- Local
        integer(kind=8) :: jvPipe, jvCodret, jvCheck, jvIndicR
        integer(kind=8) :: nbFourier, nbNode
        integer(kind=8) :: pipeType
        real(kind=8) :: criteria, tole
        type(sectPipe_Prop) :: sectPipe
        aster_logical :: lModiMetric
!   ------------------------------------------------------------------------------------------------
!
        call elrefe_info(fami='RIGI', nno=nbNode)
        call jevech('PCAORIE', 'L', jvPipe)
        call jevech('PCHCKPR', 'L', jvCheck)
        call jevech('PCODRET', "E", jvCodret)
        call jevech('PINDICR', "E", jvIndicR)
        nbFourier = 3
        if (nomteZ .eq. 'MET6SEG3') then
            nbFourier = 6
        end if
        call pipeGetType(jvPipe, nbNode, pipeType, lModiMetric)
        call pipeGetSection(sectPipe)
        tole = zr(jvCheck-1+1)
        criteria = sectPipe%thickness/sectPipe%radiusMoy
        zr(jvIndicR-1+1) = criteria
        if (lModiMetric) then
            if (criteria .gt. tole) then
                zi(jvCodret-1+1) = 2
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeEfgeElgaNlin
!
! Compute generalized forces at Gauss points (EFGE_ELGA) - Non-linear case (stresses are given)
!
! In  pipeElem         : properties of pipe element
! In  jvSigm           : adress for stresses
! In  nbSect           : number of sectors in section of pipe
! In  nbLayer          : number of layers in section of pipe
! In  weightSect       : weight of integration scheme for sectors
! In  weightLayer      : weight of integration scheme for layers
! Out efgeElga         : generalized forces at Gauss points (ELGA)
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeEfgeElgaNlin(pipeElem, &
                                jvSigm, npg, &
                                nbSect, nbLayer, &
                                nspgSect, nspgLayer, &
                                weightSect, weightLayer, &
                                efgeElga)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8), intent(in) :: jvSigm
        integer(kind=8), intent(in) :: npg
        integer(kind=8), intent(in) :: nbSect, nbLayer
        integer(kind=8), intent(in) :: nspgSect, nspgLayer
        real(kind=8), intent(in) :: weightSect(nspgSect), weightLayer(nspgLayer)
        real(kind=8), intent(out) :: efgeElga(PIPE_MAX_NPG, PIPE_NBDOF_BEAM)
! ----- Local
        integer(kind=8) :: kpg, kspg, kspgLayer, kspgSect
        type(sectPipe_Prop) :: sectPipe
        type(beamElem_Prop) :: beamElem
        real(kind=8) :: radiusLayer, phi, poids
        real(kind=8) :: efge(PIPE_NBDOF_BEAM), sigm(PIPE_TENS_SIZE)
!   ------------------------------------------------------------------------------------------------
!
        efgeElga = 0.d0
        beamElem = pipeElem%beamElem
        sectPipe = beamElem%sectPipe

! ----- Loop on Gauss points (on segment)
        kspg = 0
        do kpg = 1, npg
! --------- Loop on sub-points in thickness
            efge = 0.d0
            do kspgLayer = 1, nspgLayer
! ------------- Radius of pipe
                call pipeGetRadiusLayer(pipeElem, kspgLayer, nbLayer, radiusLayer)

! ------------- Loop on sub-points in section
                do kspgSect = 1, nspgSect
                    kspg = kspg+1

! ----------------- Coordinate of sub-point in section
                    phi = (kspgSect-1)*2.d0*r8pi()/(2.d0*nbSect)

! ----------------- Get stresses
                    sigm(1) = zr(jvSigm-1+6*(kspg-1)+1)
                    sigm(2) = zr(jvSigm-1+6*(kspg-1)+2)
                    sigm(3) = zr(jvSigm-1+6*(kspg-1)+4)
                    sigm(4) = zr(jvSigm-1+6*(kspg-1)+5)

! ----------------- Compute quadrature integral
                    poids = weightLayer(kspgLayer)*weightSect(kspgSect)* &
                            sectPipe%thickness*2.d0*r8pi()/(4.d0*nbLayer*nbSect)*radiusLayer

! ----------------- Compute value
                    efge(1) = efge(1)+poids*sigm(1)
                    efge(2) = efge(2)-poids*(sin(phi)*sigm(4)+cos(phi)*sigm(3))
                    efge(3) = efge(3)+poids*(sin(phi)*sigm(3)-cos(phi)*sigm(4))
                    efge(4) = efge(4)-poids*sigm(3)*radiusLayer
                    efge(5) = efge(5)-poids*sigm(1)*radiusLayer*cos(phi)
                    efge(6) = efge(6)+poids*sigm(1)*radiusLayer*sin(phi)

                end do
            end do
            efgeElga(kpg, :) = efge(:)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeEfgeElgaLine
!
! Compute generalized forces at Gauss points (EFGE_ELGA) - Linear case
!
! In  pipeElem         : properties of pipe element
! In  jvMaterCode      : adress for material parameters
! In  jvDisp           : adress for displacements
! In  nbDof            : number of DOF in element
! In  nspg             : number of Gauss sub-points
! In  npg              : number of Gauss points
! In  xpg              : coordinates of Gauss points
! In  nbSect           : number of sectors in section of pipe
! In  nbLayer          : number of layers in section of pipe
! In  nspgSect         : number of sub-points in sectors of pipe
! In  nspgLayer        : number of sub-points in layers of pipe
! In  weightSect       : weight of integration scheme for sectors
! In  weightLayer      : weight of integration scheme for layers
! In  ff               : shape functions on beam
! In  df1              : first derivative of shape functions on beam
! In  df2              : second derivative of shape functions on beam
! Out efgeElga         : generalized forces at Gauss points (ELGA)
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeEfgeElgaLine(pipeElem, &
                                jvMaterCode, jvDisp, &
                                nbDof, nspg, npg, xpg, &
                                nbSect, nbLayer, &
                                nspgSect, nspgLayer, &
                                weightSect, weightLayer, &
                                ff, df1, df2, &
                                efgeElga)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(pipeElem_Prop), intent(in) :: pipeElem
        integer(kind=8) :: jvMaterCode, jvDisp
        integer(kind=8), intent(in) :: nbDof
        integer(kind=8), intent(in) :: npg, nspg
        real(kind=8), intent(in) :: xpg(npg)
        real(kind=8), intent(in) :: ff(*), df1(*), df2(*)
        integer(kind=8), intent(in) :: nbSect, nbLayer
        integer(kind=8), intent(in) :: nspgSect, nspgLayer
        real(kind=8), intent(in) :: weightSect(nspgSect), weightLayer(nspgLayer)
        real(kind=8), intent(out) :: efgeElga(PIPE_TENS_SIZE, PIPE_NBDOF_BEAM)
! ----- Local
        type(sectPipe_Prop) :: sectPipe
        type(beamElem_Prop) :: beamElem
        integer(kind=8) :: iret
        integer(kind=8) :: iDofBeam, kpg, kspg, kspgLayer, kspgSect
        real(kind=8) :: meanTemp, elasMatr(PIPE_TENS_SIZE, PIPE_TENS_SIZE)
        real(kind=8) :: dispLoca(nbDof)
        real(kind=8) :: sigmTher(2), epsiTher, sigm(PIPE_TENS_SIZE)
        real(kind=8) :: radiusLayer, phi, zeta, poids
        real(kind=8) :: b(PIPE_TENS_SIZE, nbDof)
        real(kind=8) :: prodMatrMatr(nbDof, PIPE_TENS_SIZE)
        real(kind=8) :: efge(PIPE_NBDOF_BEAM)
!   ------------------------------------------------------------------------------------------------
!
        efgeElga = 0.d0
        beamElem = pipeElem%beamElem
        sectPipe = beamElem%sectPipe

! ----- Compute mean temperature (on all point and "sous-point" gauss)
        iret = 0
        call moytem('RIGI', npg, nspg, '+', meanTemp, iret)
        if (iret .ne. 0) then
            meanTemp = 0.d0
        end if

! ----- Get elastic properties
        call jevech('PMATERC', 'L', jvMaterCode)
        call pipeGetElasProp(jvMaterCode, temp_=meanTemp, &
                             c_=elasMatr)

! ----- Get displacements  (in local base)
        call jevech('PDEPLAR', 'L', jvDisp)
        call pipeGetDisp(pipeElem, nbDof, zr(jvDisp), dispLoca)

        do kpg = 1, npg
! --------- Compute thermal strain and stress
            call verifg('RIGI', kpg, nspg, '+', zi(jvMaterCode), epsiTher)
            sigmTher(1) = (elasMatr(1, 1)+elasMatr(1, 2))*epsiTher
            sigmTher(2) = (elasMatr(2, 1)+elasMatr(2, 2))*epsiTher

! --------- Loop on sub-points in thickness
            efge = 0.d0
            do kspgLayer = 1, nspgLayer
! ------------- Radius of pipe
                call pipeGetRadiusLayer(pipeElem, kspgLayer, nbLayer, radiusLayer)

! ------------- Coordinate of sub-point in thickness
                zeta = (kspgLayer-1)*pipeElem%beamElem%sectPipe%thickness/(2.d0*nbLayer)- &
                       pipeElem%beamElem%sectPipe%thickness/2.d0

! ------------- Loop on sub-points in section
                do kspgSect = 1, nspgSect
                    kspg = kspg+1

! ----------------- Coordinate of sub-point in section
                    phi = (kspgSect-1)*2.d0*r8pi()/(2.d0*nbSect)

! ----------------- Compute B matrix (u => epsi)
                    call pipeBMatr(pipeElem, nbDof, &
                                   npg, xpg, &
                                   phi, zeta, radiusLayer, &
                                   kpg, ff, df1, df2, &
                                   b)

! ----------------- Product [C] [B]
                    call promat(elasMatr, PIPE_TENS_SIZE, PIPE_TENS_SIZE, PIPE_TENS_SIZE, &
                                b, PIPE_TENS_SIZE, PIPE_TENS_SIZE, nbDof, &
                                prodMatrMatr)

! ----------------- Product {sigm} = [C] [B] {U}
                    call prmave(0, prodMatrMatr, PIPE_TENS_SIZE, PIPE_TENS_SIZE, nbDof, &
                                dispLoca, nbDof, sigm, PIPE_TENS_SIZE, iret)
                    ASSERT(iret .eq. 0)

! ----------------- Compute quadrature integral
                    poids = weightLayer(kspgLayer)*weightSect(kspgSect)* &
                            sectPipe%thickness*2.d0*r8pi()/(4.d0*nbLayer*nbSect)*radiusLayer

! ----------------- Compute value
                    efge(1) = efge(1)+poids*(sigm(1)-sigmTher(1))
                    efge(2) = efge(2)-poids*(sin(phi)*sigm(4)+cos(phi)*sigm(3))
                    efge(3) = efge(3)+poids*(sin(phi)*sigm(3)-cos(phi)*sigm(4))
                    efge(4) = efge(4)-poids*sigm(3)*radiusLayer
                    efge(5) = efge(5)-poids*(sigm(1)-sigmTher(1))*radiusLayer*cos(phi)
                    efge(6) = efge(6)+poids*(sigm(1)-sigmTher(1))*radiusLayer*sin(phi)

                end do
            end do
            do iDofBeam = 1, PIPE_NBDOF_BEAM
                efgeElga(kpg, iDofBeam) = efge(iDofBeam)
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! pipeGetDime
!
! Get main dimensions of pipe element
!
! In  nomte            : type of finite element
! In  fami             : name of integration scheme
! Out nbNode           : number of nodes
! Out nbFourier        : number of Fourier modes
! Out nbDof            : number of DOF in element
!
! --------------------------------------------------------------------------------------------------
    subroutine pipeGetDime(nomteZ, famiZ, &
                           nbNode_, nbFourier_, nbDof_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: nomteZ, famiZ
        integer(kind=8), optional, intent(out) :: nbNode_, nbFourier_, nbDof_
! ----- Local
        integer(kind=8) :: nbNode, nbFourier, nbDof
!   ------------------------------------------------------------------------------------------------
!
        call elrefe_info(fami=famiZ, nno=nbNode)
        nbFourier = 3
        if (nomteZ .eq. 'MET6SEG3') then
            nbFourier = 6
        end if
        nbDof = nbNode*(6+3+6*(nbFourier-1))
        ASSERT(nbDof .le. PIPE_MAX_DOF)

        if (nomteZ .eq. 'MET3SEG3') then
            ASSERT(nbDof .eq. 63)
        else if (nomteZ .eq. 'MET6SEG3') then
            ASSERT(nbDof .eq. 117)
        else if (nomteZ .eq. 'MET3SEG4') then
            ASSERT(nbDof .eq. 84)
        else
            ASSERT(ASTER_FALSE)
        end if

        if (present(nbNode_)) nbNode_ = nbNode
        if (present(nbDof_)) nbDof_ = nbDof
        if (present(nbFourier_)) nbFourier_ = nbFourier
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module pipeElem_module
