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
subroutine assmam(jvBase, matrAsseZ, &
                  nbMatrElem, listMatrElem, coefMatrElem, &
                  numeDofZ, motcle, matrScalType)
!
!----------------------------------------------------------------------
!  attention : cette routine ne doit pas etre appellee directement :
!              il faut appeler son "chapeau" : asmatr
!----------------------------------------------------------------------
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/cheksd.h"
#include "asterfort/asmpi_barrier.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/assma1.h"
#include "asterfort/assma2.h"
#include "asterfort/assma3.h"
#include "asterfort/crelil.h"
#include "asterfort/detrsd.h"
#include "asterfort/digdel.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jaexin.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedbg2.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jerazo.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/nbno.h"
#include "asterfort/typmat.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpr.h"
#include "asterfort/uttcpu.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerobj.h"
#include "asterfort/nbddlMaxMa.h"
#include "asterfort/getDistributionParameters.h"
#include "asterfort/isParallelMatrix.h"
!
    character(len=1), intent(in) :: jvBase
    integer(kind=8), intent(in) :: nbMatrElem
    character(len=*), intent(in) :: matrAsseZ, listMatrElem(nbMatrElem), numeDofZ
    integer(kind=8), intent(in) :: matrScalType
    real(kind=8), intent(in) :: coefMatrElem(*)
    character(len=4), intent(in) :: motcle
!-----------------------------------------------------------------------
! Assemblage Morse avec preconditionnement des matr_elem de mailles
! "Lagrange".
!-----------------------------------------------------------------------
! int k* base   : base sur laquelle on veut creer la matr_asse
! out k* matrAsseZ  : l'objet matrAsseZ de type matr_asse est cree et rempli
! in  k* matrAsseZ  : nom de l'objet de type matr_asse a creer
! in  i  nbMatrElem  : nombre de matr_elem  de la liste listMatrElem
! in  k* listMatrElem : liste des matr_elem
! in  i  coefMatrElem : liste des coefficients multiplicateurs des matr_elem
! in  k* numeDofZ     : nom du numero_ddl
! in  k4 motcle : 'ZERO' ou 'cumu'
!                 'ZERO':si un objet de nom matrAsseZ et de type
!                        matr_asse existe on l'ecrase
!                 'CUMU':si un objet de nom matrAsseZ et de type
!                        matr_asse existe on l'enrichi
! in  i   matrScalType  : type (r/c) de la matr_asse
!                          1 --> reelles
!                          2 --> complexes
!-----------------------------------------------------------------------
    aster_logical, parameter :: dbg = ASTER_FALSE
    character(len=16) :: optio, optio2
    character(len=1) :: typsca
    character(len=2) :: tt
    character(len=8) ::  nogdco, nogdsi, mesh, mesh2, model, mo2
    character(len=8) :: symel, kempic
    character(len=14) :: numeDof, nu14
    character(len=19) :: matrAsse, mat19, resu, matel, ligre1
    character(len=24) :: cMeshName
    character(len=1) :: matsym
    character(len=3) :: matd, answer
    real(kind=8) :: c1, temps(7)
    aster_logical :: acreer, cumul, ldistme, lmatd, l_parallel_matrix
    aster_logical :: lmasym, lmesym, ldgrel, lparallel_mesh, lligrel_cp
    integer(kind=8) :: admodl, i, nbi1mx
    integer(kind=8) :: jdesc
    integer(kind=8) :: jadli, jadne, jnueq, jnulo1
    integer(kind=8) :: jposd1, jtmp2, lgtmp2
    integer(kind=8) :: ibid, iconx1, iconx2, idbgav
    integer(kind=8) :: jprn1, jprn2, jresl, maxDDLMa
    integer(kind=8) :: iel, ier, ifm, igr
    integer(kind=8) :: ilima, ilimo, ilinu
    integer(kind=8) :: imat, iresu
    integer(kind=8) :: iret, itbloc
    integer(kind=8) :: jrefa, jsmdi, jsmhc, jvalm(2)
    integer(kind=8) :: lcmodl, mode, n1, nbelm
    integer(kind=8) :: nblc, nbnomx, nbnoss, nbresu
    integer(kind=8) :: ncmp, nbvel, nec, nel, nequ, nbproc
    integer(kind=8) :: niv, nlili, nmxcmp, nnoe
    integer(kind=8) :: nugd, rang, ieq, idia, ellagr, iexi, ilinu2
    mpi_int :: msize
    integer(kind=8), pointer :: smde(:) => null()
    character(len=24), pointer :: noli(:) => null()
    character(len=24), pointer :: relr(:) => null()
    character(len=24), pointer :: tco(:) => null()
    character(len=24), pointer :: v_cmno(:) => null()
    integer(kind=8), pointer :: assma3_tab1(:) => null(), assma3_tab2(:) => null()
    integer(kind=8), pointer :: numsd(:) => null()
    integer(kind=8), pointer :: v_crco(:) => null()
    integer(kind=8), pointer :: v_refp(:) => null()

!-----------------------------------------------------------------------
!     FONCTIONS FORMULES :
!-----------------------------------------------------------------------

#define zzngel(ili) zi(jadli+3*(ili-1))
#define zznelg(ili,igrel) zi(zi(jadli+3*(ili-1)+2)+igrel)- \
    zi(zi(jadli+3*(ili-1)+2)+igrel-1)-1
#define zzliel(ili,igrel,iel) zi(zi(jadli+3*(ili-1)+1)-1+ \
    zi(zi(jadli+3*(ili-1)+2)+igrel-1)+iel-1)

!----------------------------------------------------------------------
! Gestion du parallisme :
! -----------------------
! La routine assemble des matr_elem pour en faire une matr_asse.
! Si les matr_elem sont "distribues" (c'est a dire que chaque processeur
! ne connait qu'une partie des matrices elementaires), on va produire
! une matr_asse "distribuee" (de contenu different sur chaque processeur).
!
! Une matr_asse "distribuee" peut etre de petite taille (MATR_DISTRIBUEE='OUI')
! ou non (MATR_DISTRIBUEE='NON').
!
! 2 booleens pilotent le parallelisme de cette routine :
!  ldistme : il existe au moins un matr_elem distribue
!            => la matr_asse produite sera "distribuee"
!  lmatd : l'utilisateur a demande MATR_DISTRIBUEE='OUI'
!
! On verifie que :
!   lmatd=.T.    => ldistme=.T.
!
!
! Precisions :
! ------------
!  Le booleen lmatd sert essentiellement a determiner le stockage de la matr_asse
!  qu'il faut utiliser : un stockage global ou un stockage local.
!
!  La decision de produire une matr_asse distribuee est prise des que l'on
!  trouve un (ou plusieurs) resuelem distribue(s) dans les matr_elem a assembler.
!
!  Quand une matr_asse est distribuee, la matrice "totale" peut etre obtenue
!  en faisant (au moins par la pensee) une "simple" somme des matr_asse possedees
!  par les differents processeurs. Il est donc fondamental qu'une matrice elementaire
!  (ou la matrice d'un macro-element) ne soit assemblee que sur un seul processeur.
!
!  Si un resuelem est distribue, on peut recuperer la partition attachee a ce
!  resuelem. C'est cette parttion qui servira pour l'assemblage : chaque processeur
!  n'assemble que "ses" elements dans la partition.
!  Si plusieurs resuelem sont distribues, on verifie que leurs partitions sont
!  identiques. Sinon : erreur <F>.
!  Le nombre de processeurs lors de l'assemblage doit etre identique a celui de la
!  partition.
!
!  Les matrices liees aux macro-elements sont attachees aux matr_elem.
!  Elles sont toujours calculees (et identiques) sur TOUS les processeurs
!  (operateur MACR_ELEM_STAT).
!  Si la matr_asse est distribuee, c'est le processeur 0 (et lui seul) qui va assembler
!  les matrices des macro-elements.
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call jedbg2(idbgav, 0)
    call infniv(ifm, niv)

    matrAsse = matrAsseZ
    cMeshName = ' '

! - Get numbering
    numeDof = numeDofZ
    if (dbg) call cheksd(numeDof, 'SD_NUME_DDL', iret)

! - Get model
    call dismoi('NOM_MODELE', numeDof, 'NUME_DDL', repk=model)

! - Get mesh
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NOM_MAILLA', numeDof, 'NUME_DDL', repk=mesh2)
    ASSERT(mesh .eq. mesh2)
    call dismoi('PARALLEL_MESH', mesh, 'MAILLAGE', repk=answer)
    lparallel_mesh = (answer .eq. 'OUI')
    if (.not. lparallel_mesh) then
        call asmpi_barrier()
    end if

    call uttcpu('CPU.CALC.1', 'DEBUT', ' ')
    call uttcpu('CPU.ASSE.1', 'DEBUT', ' ')
    call uttcpu('CPU.ASSE.2', 'DEBUT', ' ')

    call dismoi('NB_NO_SS_MAX', mesh, 'MAILLAGE', repi=nbnoss)
    call dismoi('NOM_GD', numeDof, 'NUME_DDL', repk=nogdco)
    call dismoi('NOM_GD_SI', nogdco, 'GRANDEUR', repk=nogdsi)
    call dismoi('NB_CMP_MAX', nogdsi, 'GRANDEUR', repi=nmxcmp)
    ncmp = nmxcmp
    call dismoi('NUM_GD_SI', nogdsi, 'GRANDEUR', repi=nugd)
    nec = nbec(nugd)
    call jeveuo(jexatr('&CATA.TE.MODELOC', 'LONCUM'), 'L', lcmodl)
    call jeveuo(jexnum('&CATA.TE.MODELOC', 1), 'L', admodl)
    call jeexin(mesh//'.CONNEX', iret)
    if (iret .gt. 0) then
        call jeveuo(mesh//'.CONNEX', 'L', iconx1)
        call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', iconx2)
    else
        iconx1 = 0
        iconx2 = 0
    end if
!   ellagr : 0 : il n'existe pas d'element de lagrange
!            1 : il existe des elements de lagrange
    ellagr = 0

!   -- calcul de lmatd et jnueq :
!   -----------------------------
    call dismoi('MATR_DISTRIBUEE', numeDof, 'NUME_DDL', repk=matd)
    lmatd = (matd .eq. 'OUI')
    if (lmatd) then
        call jeveuo(numeDof//'.NUML.NUEQ', 'L', jnueq)
    else
        call jeveuo(numeDof//'.NUME.NUEQ', 'L', jnueq)
    end if

!   -- calcul de :
!   --------------
!     lmasym: .true   : matrice assemblee symetrique
!             .false. : matrice assemblee non-symetrique
!     acreer: .true.  : il faut creer la matr_asse
!             .false. : la matr_asse existe deja
!     cumul : .true.  : on accumule dans la matr_asse
!             .false. : on remet la matr_asse a zero
!                       (elle doit exister)
!     tt  : tt(1) : type (r/c) de ce que l'on assemble
!           tt(2) : type (r/c) de la sd_matr_asse
!------------------------------------------------------------
    matsym = typmat(nbMatrElem, listMatrElem)
    ASSERT(matsym .eq. 'S' .or. matsym .eq. 'N')
    lmasym = (matsym .eq. 'S')
    if (motcle(1:4) .eq. 'ZERO') then
        cumul = .false.
        acreer = .true.
    else if (motcle(1:4) .eq. 'CUMU') then
        cumul = .true.
        acreer = .false.
        call jelira(matrAsse//'.VALM', 'NMAXOC', nblc)
        ASSERT(nblc .eq. 1 .or. nblc .eq. 2)
        if (nblc .eq. 2) lmasym = .false.
    else
        ASSERT(.false.)
    end if
    ASSERT((matrScalType .eq. 1) .or. (matrScalType .eq. 2))
    if (matrScalType .eq. 1) then
        tt = '?R'
    else if (matrScalType .eq. 2) then
        tt = '?C'
    end if

    call jelira(numeDof//'.NUME.REFN', 'LONMAX', n1)
    ASSERT(n1 .eq. 5)

! - Get parameters for distribution of elementary matrices
    call getDistributionParameters(nbMatrElem, listMatrElem, &
                                   ldistme, ldgrel, &
                                   rang, nbproc, &
                                   numsd)

    if (lmatd) then
        ASSERT(ldistme)
    end if

!   -- allocation des objets .NUMLOX et .POSDDX:
!   ----------------------------------------------
!   50 est suppose etre le + gd nombre de noeuds d'une maille
!      standard (jusqu'a present : 27 (hexa27))
    nbnomx = max(nbnoss, 50)
    call wkvect('&&ASSMAM.NUMLO1', 'V V I', 2*nbnomx, jnulo1)
    call wkvect('&&ASSMAM.POSDD1', 'V V I', nbnomx*nmxcmp, jposd1)
!
!   -- allocation d'un objet de travail utilise dans asretm :
!      ce vecteur est agrandi si necessaire dans asretm
    lgtmp2 = 400
    call wkvect('&&ASSMAM.TMP2', 'V V I', lgtmp2, jtmp2)

    if (acreer) then
        call detrsd('MATR_ASSE', matrAsse)
    else
        mat19 = matrAsse
        nu14 = numeDof
        call jeveuo(mat19//'.REFA', 'L', jrefa)
        ASSERT(zk24(jrefa-1+2) (1:14) .eq. nu14)
        call jedetr(mat19//'.REFA')
    end if
!  -- calcul d un repertoire,temporaire, matrAsse.lili a partir
!     de la liste de matrices elementaires

    call crelil('F', nbMatrElem, listMatrElem, matrAsse//'.LILI', 'V', &
                '&MAILLA', matrAsse, ibid, mesh, ibid, &
                ibid, ilimo, nlili, nbelm, nume_=numeDof)
    call jeveuo(matrAsse//'.ADLI', 'E', jadli)
    call jeveuo(matrAsse//'.ADNE', 'E', jadne)

    if (niv .ge. 2) then
        call uttcpu('CPU.ASSMAM', 'INIT ', ' ')
        call uttcpu('CPU.ASSMAM', 'DEBUT', ' ')
    end if
!
! --- On alloue ici les tableaux de travail pour assma3
!     l'augmentation du temps d'assemblage est négligeable
!
    maxDDLMa = nbddlMaxMa(numeDof, matrAsse, nbMatrElem, listMatrElem)
! Le x2 est par sécurité - certains processeur peuvent ne pas avoir de ddls à assembler
    nbi1mx = max(1, 2*maxDDLMa*maxDDLMa)
    call wkvect('&&ASSMAM.TAB1', 'V V I', nbi1mx, vi=assma3_tab1)
    call wkvect('&&ASSMAM.TAB2', 'V V I', nbi1mx, vi=assma3_tab2)
!
!   -- calcul de mat19, nu14, jsmhc, jsmdi, ... :
!   ----------------------------------------------
    mat19 = matrAsse
    nu14 = numeDof

    call jeveuo(nu14//'.SMOS.SMHC', 'L', jsmhc)
    call jeveuo(nu14//'.SMOS.SMDI', 'L', jsmdi)
    call jeexin(nu14//'.NUME.REFP', iret)
    if (iret .ne. 0) then
        call jeveuo(nu14//'.NUME.REFP', 'L', vi=v_refp)
    end if
    if (lmatd) then
        call jeveuo(nu14//'.NUML.PRNO', 'L', jprn1)
        call jeveuo(jexatr(nu14//'.NUML.PRNO', 'LONCUM'), 'L', jprn2)
    else
        call jeveuo(nu14//'.NUME.PRNO', 'L', jprn1)
        call jeveuo(jexatr(nu14//'.NUME.PRNO', 'LONCUM'), 'L', jprn2)
    end if

!   -- creation et remplissage de .REFA
!   -------------------------------------
    call wkvect(mat19//'.REFA', jvBase//' V K24', 20, jrefa)
    zk24(jrefa-1+1) = mesh
    zk24(jrefa-1+2) = nu14
    zk24(jrefa-1+8) = 'ASSE'
    if (lmasym) then
        zk24(jrefa-1+9) = 'MS'
    else
        zk24(jrefa-1+9) = 'MR'
    end if
    zk24(jrefa-1+10) = 'NOEU'

    call jeveuo(nu14//'.SMOS.SMDE', 'L', vi=smde)
    nequ = smde(1)
    itbloc = smde(2)
    ASSERT(smde(3) .eq. 1)
    if (lmasym) then
        nblc = 1
    else
        nblc = 2
    end if

!   -- allocation (ou non) de .VALM :
!   ---------------------------------
    if (acreer) then
        call jecrec(mat19//'.VALM', jvBase//' V '//tt(2:2), 'NU', 'DISPERSE', 'CONSTANT', &
                    nblc)
        call jeecra(mat19//'.VALM', 'LONMAX', itbloc)
        do i = 1, nblc
            call jecroc(jexnum(mat19//'.VALM', i))
        end do
    else
        if (.not. cumul) then
            do i = 1, nblc
                call jerazo(jexnum(mat19//'.VALM', i), itbloc, 1)
            end do
        end if
    end if

!   -- mise en memoire des 1 (ou 2) blocs de .VALM :
!   ------------------------------------------------
    call jeveuo(jexnum(mat19//'.VALM', 1), 'E', jvalm(1))
    call jelira(jexnum(mat19//'.VALM', 1), 'TYPE', cval=typsca)
    ASSERT(tt(2:2) .eq. typsca)
    if (.not. lmasym) then
        call jeveuo(jexnum(mat19//'.VALM', 2), 'E', jvalm(2))
    else
        jvalm(2) = 0
    end if

!   3. boucle sur les matr_elem
!   =============================
    do imat = 1, nbMatrElem
        c1 = coefMatrElem(imat)
        matel = listMatrElem(imat)
        call dismoi('NOM_MODELE', matel, 'MATR_ELEM', repk=mo2)
        call dismoi('SUR_OPTION', matel, 'MATR_ELEM', repk=optio)

        if (imat .eq. 1) then
            optio2 = optio
        else
            if (optio2 .ne. optio) optio2 = '&&MELANGE'
        end if

        if (mo2 .ne. model) then
            call utmess('F', 'ASSEMBLA_5')
        end if

!       3.1 traitement des macro-elements :
!       ----------------------------------
        call assma2(ldistme, lmasym, tt, nu14, ncmp, matel, &
                    c1, jvalm, jtmp2, lgtmp2)

!       3.2 traitement des elements finis classiques
!       -------------------------------------------
        call jeexin(matel//'.RELR', iret)
        if (iret .eq. 0) goto 80

        call jelira(matel//'.RELR', 'LONUTI', nbresu)
        if (nbresu .gt. 0) call jeveuo(matel//'.RELR', 'L', vk24=relr)

!       -- boucle sur les resu_elem
!       ============================
        do iresu = 1, nbresu
            resu = relr(iresu) (1:19)
            call jeexin(resu//'.DESC', ier)
            if (ier .eq. 0) goto 70

            call dismoi('MPI_COMPLET', resu, 'RESUELEM', repk=kempic)
            if (kempic .eq. 'NON') then
                ASSERT(ldistme)
            end if

!           -- parfois, certains resuelem sont == 0.
            if (zerobj(resu//'.RESL')) goto 70

!           -- nom du ligrel
            call jeveuo(resu//'.NOLI', 'L', vk24=noli)
            ligre1 = noli(1) (1:19)

            call jenonu(jexnom(matrAsse//'.LILI', ligre1), ilima)
            call jenonu(jexnom(nu14//'.NUME.LILI', ligre1), ilinu)
            lligrel_cp = .false.
            call jeexin(noli(1) (1:19)//'._TCO', ier)
            if (ier .ne. 0) then
                call jeveuo(noli(1) (1:19)//'._TCO', "L", vk24=tco)
                lligrel_cp = (tco(1) .eq. 'LIGREL_CP')
            end if
            if (lligrel_cp) then
                call asmpi_info(size=msize)
                if (msize > 1) then
                    call jeveuo(jexnum(nu14//'.NUME.CRCO', ilinu), 'L', vi=v_crco)
                end if
                call jeveuo(ligre1//".CMNO", 'L', vk24=v_cmno)
                cMeshName = v_cmno(1)
                ilinu2 = v_refp(ilinu)
            else
                ilinu2 = ilinu
            end if

            call dismoi('TYPE_SCA', resu, 'RESUELEM', repk=typsca)
            ASSERT(typsca .eq. 'R' .or. typsca .eq. 'C')
            tt(1:1) = typsca
            call dismoi('TYPE_MATRICE', resu, 'RESUELEM', repk=symel)
            ASSERT(symel(1:1) .eq. 'S' .or. symel(1:1) .eq. 'N')
            lmesym = (symel(1:1) .eq. 'S')
            if (lmasym) then
                ASSERT(lmesym)
            end if
!
            call jeveuo(resu//'.DESC', 'L', jdesc)
!
!           -- boucle sur les grels du ligrel
!           ==================================
            do igr = 1, zzngel(ilima)
                if (ldgrel .and. mod(igr, nbproc) .ne. rang) goto 60

!               -- il se peut que le grel igr soit vide :
                call jaexin(jexnum(resu//'.RESL', igr), iexi)
                if (iexi .eq. 0) goto 60

                mode = zi(jdesc+igr+1)
                if (mode .gt. 0) then
                    nnoe = nbno(mode)
                    ASSERT(nnoe .le. nbnomx)
                    nbvel = digdel(mode)
!                           -- nombre d'elements du grel igr du ligrel ligre1/ilima
                    nel = zznelg(ilima, igr)
                    call jeveuo(jexnum(resu//'.RESL', igr), 'L', jresl)

!                   boucle sur les elements du grel
!                   -------------------------------
                    do iel = 1, nel
                        call assma3(lmasym, lmesym, tt, igr, iel, &
                                    c1, rang, jnueq, numsd, jresl, &
                                    nbvel, nnoe, ldistme, ldgrel, &
                                    ilima, jadli, jadne, jprn1, jprn2, &
                                    jnulo1, jposd1, admodl, &
                                    lcmodl, mode, nec, nmxcmp, ncmp, &
                                    jsmhc, jsmdi, iconx1, iconx2, jtmp2, &
                                    lgtmp2, jvalm, ilinu2, ellagr, &
                                    nbi1mx, assma3_tab1, assma3_tab2, &
                                    v_crco, lligrel_cp)
                    end do
                    call jelibe(jexnum(resu//'.RESL', igr))
                end if
60              continue
            end do
70          continue
        end do
80      continue
    end do

!   -- mise a jour de REFA(4)
    call jeveuo(mat19//'.REFA', 'E', jrefa)
    if (acreer) then
        zk24(jrefa-1+4) = optio2
    else
        if (zk24(jrefa-1+4) .ne. optio2) zk24(jrefa-1+4) = '&&MELANGE'
    end if
    ! NSELLENET ce lien est-il legitime ??? !!!!!!!!!!!!!
    zk24(jrefa-1+5) = cMeshName
    ! NSELLENET !!!!!!!!!!!!!

    if (niv .ge. 2) then
        call uttcpu('CPU.ASSMAM', 'FIN', ' ')
        call uttcpr('CPU.ASSMAM', 7, temps)
        write (ifm, '(A44,D11.4,D11.4,D11.4)') &
            'TEMPS CPU/SYS/ELAPSED ASSEMBLAGE M        : ', &
            temps(5), temps(6), temps(7)
    end if

    if (.not. lmasym) then
!       -- par prudence, on affecte aux termes diagonaux du bloc inferieur
!          les valeurs des termes diagonaux du bloc superieur
        do ieq = 1, nequ
            idia = zi(jsmdi+ieq-1)
            zr(jvalm(2)+idia-1) = zr(jvalm(1)+idia-1)
        end do
    end if

    l_parallel_matrix = isParallelMatrix(mat19)
!   -- il faut communiquer ellagr entre les procs :
    if (ldistme .or. l_parallel_matrix) then
        call asmpi_comm_vect('MPI_MAX', 'I', sci=ellagr)
    end if

!   -- mise a l'echelle des coef. de lagrange si necessaire :
    if (ellagr .gt. 0) call assma1(mat19, ldistme, l_parallel_matrix)

    if (.not. ldistme) then
        zk24(jrefa-1+11) = 'MPI_COMPLET'
    else
        if (lmatd) then
            zk24(jrefa-1+11) = 'MATR_DISTR'
        else
            zk24(jrefa-1+11) = 'MPI_INCOMPLET'
        end if
    end if
!
    call jedetr('&&ASSMAM.TAB1')
    call jedetr('&&ASSMAM.TAB2')
    call jedetr(matrAsse//'.ADNE')
    call jedetr(matrAsse//'.ADLI')
    call jedetr('&&ASSMAM.NUMLO1')
    call jedetr('&&ASSMAM.POSDD1')
    call jedetr('&&ASSMAM.TMP2')
    call jedbg2(ibid, idbgav)
    if (dbg) call cheksd(matrAsse, 'SD_MATR_ASSE', iret)

    if (.not. lparallel_mesh) then
        call asmpi_barrier()
    end if
    call uttcpu('CPU.CALC.1', 'FIN', ' ')
    call uttcpu('CPU.ASSE.1', 'FIN', ' ')
    call uttcpu('CPU.ASSE.2', 'FIN', ' ')
    call jedema()
end subroutine
