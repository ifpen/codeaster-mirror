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

subroutine promor(nuz, base, printz)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/infniv.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/nbec.h"
#include "asterfort/utmess.h"
#include "asterfort/uttrii_i4.h"
#include "asterfort/wkvect.h"
#include "asterc/vector_vector_add_values.h"
#include "asterc/vector_vector_delete.h"
#include "asterc/vector_vector_get_data.h"
#include "asterc/vector_vector_get_data_sizes.h"
#include "asterc/vector_vector_new.h"
#include "asterc/vector_vector_resize.h"
#include "asterc/vector_vector_sort_unique.h"
!
    character(len=*) :: nuz
    character(len=1) :: base
    aster_logical, optional :: printz
!     ------------------------------------------------------------------
!     CALCUL DE LA STRUCTURE COMPACTE D'UNE MATRICE
!     ------------------------------------------------------------------
! IN  K*14 NU     : NOM DE LA SD_UME_DDL A COMPLETER.
! OUT K*14 NU     : L'OBJET NU EST COMPLETE DES OBJETS .SMOS.XXXX
! IN  K1  BASE    : BASE DE CREATION DU STOCKAGE
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    character(len=8) :: ma, mo, exiele, partit
!----------------------------------------------------------------------
    character(len=14) :: nu
    aster_logical :: ldist, ldgrel, lmadis, printt, lligrel_cp
    character(len=19) :: nomlig
    character(len=24) :: crco
    integer(kind=8) :: iconx2, ili, iel
    integer(kind=8) :: idprn1, iexi
    integer(kind=8) :: idprn2, ifm, niv, iret, nnoe, jnueq
    integer(kind=8) :: vali(3), neqx, iilib, igr, numa, k1, n1, iad1, nddl1, n12
    integer(kind=8) :: iddl, iamail, ncoef, jsmde, igd, nbss
    integer(kind=8) :: iadequ, nlili, nequ, mxddlt
    integer(kind=8) :: ima, nddlt, jalm, nel, nec, nbsma
    integer(kind=8) ::  rang, imd, ili2, vec_ptr, vecSize, position
!
    real(kind=8) :: valr(2), rcoef, requ
    integer(kind=8) :: nbproc
    integer(kind=8), pointer :: maille(:) => null()
    integer(kind=8), pointer :: adne(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: adli(:) => null()
    character(len=24), pointer :: prtk(:) => null()
    character(len=24), pointer :: tco(:) => null()
    integer(kind=8), pointer :: sssa(:) => null()
    integer(kind=8), pointer :: v_refp(:) => null()
    integer(kind=8), pointer :: v_crco(:) => null()
    integer(kind=8), pointer :: v_smdi(:) => null()
    integer(kind=4), pointer :: v_smhc(:) => null()
    mpi_int :: mrank, msize
!
!-----------------------------------------------------------------------
!     FONCTIONS LOCALES D'ACCES AUX DIFFERENTS CHAMPS DES
!     S.D. MANIPULEES DANS LE SOUS PROGRAMME
!-----------------------------------------------------------------------
!---- FONCTION D ACCES AU CHAMP CONNEX DE LA S.D. MAILLA DE TYPE
!     MAILLAGE
!     ZZCONX(IMAIL,J) = NUMERO DANS LA NUMEROTATION DU MAILLAGE
!         DU NOEUD J DE LA MAILLE IMAIL
#define zzconx(imail,j) connex(zi(iconx2+imail-1)+j-1)
!
!---- NBRE DE NOEUDS DE LA MAILLE IMAIL DU MAILLAGE
!
#define zznbne(imail) zi(iconx2+imail)-zi(iconx2+imail-1)
!
!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS LIEL DES S.D. LIGREL
!     REPERTORIEES DANS LE REPERTOIRE TEMPORAIRE .MATAS.LILI
!     ZZLIEL(ILI,IGREL,J) =
!      SI LA JIEME MAILLE DU LIEL IGREL DU LIGREL ILI EST:
!          -UNE MAILLE DU MAILLAGE : SON NUMERO DANS LE MAILLAGE
!          -UNE MAILLE TARDIVE : -POINTEUR DANS LE CHAMP .NEMA
!
#define zzliel(ili,igrel,j) zi(adli(1+3*(ili-1)+1)-1+ \
    zi(adli(1+3*(ili-1)+2)+igrel-1)+j-1)
!
!---- NBRE DE GROUPES D'ELEMENTS (DE LIEL) DU LIGREL ILI
!
#define zzngel(ili) adli(1+3*(ili-1))
!
!---- NBRE DE NOEUDS DE LA MAILLE TARDIVE IEL ( .NEMA(IEL))
!     DU LIGREL ILI REPERTOIRE .LILI
!     (DIM DU VECTEUR D'ENTIERS .LILI(ILI).NEMA(IEL) )
!
#define zznsup(ili,iel) zi(adne(1+3*(ili-1)+2)+iel)- \
    zi(adne(1+3*(ili-1)+2)+iel-1)-1
!
!---- NBRE D ELEMENTS DU LIEL IGREL DU LIGREL ILI DU REPERTOIRE TEMP.
!     .MATAS.LILI(DIM DU VECTEUR D'ENTIERS .LILI(ILI).LIEL(IGREL) )
!
#define zznelg(ili,igrel) zi(adli(1+3*(ili-1)+2)+igrel)- \
    zi(adli(1+3*(ili-1)+2)+igrel-1)-1
!
!---- NBRE D ELEMENTS SUPPLEMENTAIRE (.NEMA) DU LIGREL ILI DU
!     REPERTOIRE TEMPORAIRE .MATAS.LILI
!
!
!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS NEMA DES S.D. LIGREL
!     REPERTORIEES DANS LE REPERTOIRE TEMPO. .MATAS.LILI
!
#define zznema(ili,iel,j) zi(adne(1+3*(ili-1)+1)-1+ \
    zi(adne(1+3*(ili-1)+2)+iel-1)+j-1)
!
!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS PRNO DES S.D. LIGREL
!     REPERTORIEES DANS NU.LILI DE LA S.D. NUME_DDL ET A LEURS ADRESSES
!     ZZPRNO(ILI,NUNOEL,1) = NUMERO DE L'EQUATION ASSOCIEES AU 1ER DDL
!                            DU NOEUD NUNOEL DANS LA NUMEROTATION LOCALE
!                            AU LIGREL ILI DE .LILI
!     ZZPRNO(ILI,NUNOEL,2) = NOMBRE DE DDL PORTES PAR LE NOEUD NUNOEL
!     ZZPRNO(ILI,NUNOEL,2+1) = 1ER CODE
!     ZZPRNO(ILI,NUNOEL,2+NEC) = NEC IEME CODE
!
#define zzprno(ili,nunoel,l) zi(idprn1-1+zi(idprn2+ili-1)+ \
    (nunoel-1)*(nec+2)+l-1)
!----------------------------------------------------------------------
!
    call infniv(ifm, niv)
    call jemarq()
    nu = nuz
    printt = ASTER_TRUE
    if (present(printz)) then
        printt = printz
    end if
!
!
    call dismoi('NOM_MODELE', nu, 'NUME_DDL', repk=mo)
    call dismoi('NUM_GD_SI', nu, 'NUME_DDL', repi=igd)
    call dismoi('NOM_MAILLA', nu, 'NUME_DDL', repk=ma)
!
!---- QUEL TYPE DE PARTITION ?
!     LDIST=.TRUE.  : LES CALCULS ELEMENTAIRES SONT DISTRIBUES
!     LDGREL=.TRUE. : DISTRIBUTION DE TYPE 'GROUP_ELEM'
    ldist = .false.
    ldgrel = .false.
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)

    if (mo .ne. ' ') then
        call dismoi('NOM_LIGREL', mo, 'MODELE', repk=nomlig)
        call dismoi('PARTITION', nomlig, 'LIGREL', repk=partit)
!
        call exisd('PARTITION', partit, iexi)
        if (iexi .ne. 0) then
            ldist = .true.
            call jeveuo(partit//'.PRTK', 'L', vk24=prtk)
            ldgrel = prtk(1) .eq. 'SOUS_DOMAINE' .or. prtk(1) .eq. 'GROUP_ELEM'
            if (.not. ldgrel) then
                call jeveuo(partit//'.NUPR', 'L', vi=maille)
            end if
        end if
    end if
!
    call jeexin(ma//'.CONNEX', iret)
    if (iret .gt. 0) then
        call jeveuo(ma//'.CONNEX', 'L', vi=connex)
        call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', iconx2)
    end if
!
    if (mo .eq. ' ') then
        nbss = 0
    else
        call dismoi('NB_SS_ACTI', mo, 'MODELE', repi=nbss)
        if (nbss .gt. 0) then
            call dismoi('NB_SM_MAILLA', mo, 'MODELE', repi=nbsma)
            call jeveuo(mo//'.MODELE    .SSSA', 'L', vi=sssa)
        end if
    end if
!
!
    call jeveuo(nu//'     .ADNE', 'E', vi=adne)
    call jeveuo(nu//'     .ADLI', 'E', vi=adli)
!
    crco = nu//'.NUME.CRCO'
!
!     -- CAS MATR_DISTRIBUE='OUI' => LMADIS=.TRUE.
    call jeexin(nu//'.NUML.DELG', imd)
    lmadis = (imd .ne. 0)
    if (.not. lmadis) then
        call jeveuo(nu//'.NUME.NEQU', 'L', iadequ)
        call jelira(nu//'.NUME.PRNO', 'NMAXOC', nlili)
        call jeveuo(nu//'.NUME.PRNO', 'L', idprn1)
        call jeveuo(jexatr(nu//'.NUME.PRNO', 'LONCUM'), 'L', idprn2)
        call jeveuo(nu//'.NUME.NUEQ', 'L', jnueq)
    else
        call jeveuo(nu//'.NUML.NEQU', 'L', iadequ)
        call jelira(nu//'.NUME.PRNO', 'NMAXOC', nlili)
        call jeveuo(nu//'.NUML.PRNO', 'L', idprn1)
        call jeveuo(jexatr(nu//'.NUML.PRNO', 'LONCUM'), 'L', idprn2)
        call jeveuo(nu//'.NUML.NUEQ', 'L', jnueq)
    end if
    call jeexin(nu//'.NUME.REFP', iret)
    if (iret .ne. 0) then
        call jeveuo(nu//'.NUME.REFP', 'L', vi=v_refp)
    end if
!
    nec = nbec(igd)
    nequ = zi(iadequ)
!
    call vector_vector_new(vec_ptr)
    call vector_vector_resize(vec_ptr, nequ)
!
!     -- ALLOCATION DU VECTEUR &&PROMOR.ANCIEN.LM
!     CE VECTEUR SERA AGRANDI SI NECESSAIRE
    mxddlt = 100
    call wkvect('&&PROMOR.ANCIEN.LM', 'V V S', mxddlt, jalm)
!
!     -- ALLOCATION DE .SMOS.SMDI
    call jeexin(nu//'.SMOS', iret)
    if (iret == 0) then
        call detrsd('STOC_MORSE', nu//'.SMOS')
    end if
!
!---- INITIALISATION DES TABLEAUX POUR LE STOCKAGE MORSE
!     ATTENTION:   PENDANT LA CONSTRUCTION DE LA STRUCTURE CHAINEE
!                  (SMDI,SMHC,ISUIV) DE LA MATRICE ON A
!     ZI(JSMDI+.): POINTEUR DEBUT DE CHAINE
!     ZI(JSMHC-1+II) : MAILLON QUI CONTIENT L'INDICE COLONNE
!                       DE LA CHAINE II
!     ISUIV(II)     : MAILLON SUIVANT DE LA MEME CHAINE.
!
!     NEQX   : COMPTEUR DU NOMBRE D'EQUATION (CONTROLE)
    neqx = 0
!
!     IILIB  : 1-ERE PLACE LIBRE
    iilib = 1
!
!     -- BOUCLE SUR LES LIGREL DE NU//'.NUME.LILI' :
!     -----------------------------------------------
    do ili = 2, nlili
        call jenuno(jexnum(nu//'.NUME.LILI', ili), nomlig)
        call jeexin(nomlig//'._TCO', iret)
        lligrel_cp = .false.
        if (iret .ne. 0) then
            call jeveuo(nomlig//'._TCO', "L", vk24=tco)
            lligrel_cp = (tco(1) .eq. 'LIGREL_CP')
        end if
        if (lligrel_cp) then
            ili2 = v_refp(ili)
            call jeexin(jexnum(crco, ili), iret)
            if (iret .ne. 0) then
                call jeveuo(jexnum(crco, ili), 'L', vi=v_crco)
            end if
        else
            ili2 = ili
        end if
        call dismoi('EXI_ELEM', nomlig, 'LIGREL', repk=exiele)
!
        if (nomlig(1:8) .eq. mo) then
            call dismoi('NB_SS_ACTI', mo, 'MODELE', repi=nbss)
        else
            nbss = 0
        end if
        if (exiele .eq. 'NON') goto 90
!
!       1. TRAITEMENT DES ELEMENTS FINIS CLASSIQUES:
!       --------------------------------------------
        do igr = 1, zzngel(ili)
            if (lmadis) then
                if (ldgrel .and. mod(igr, nbproc) .ne. rang) goto 80
            end if
            nel = zznelg(ili, igr)
            do iel = 1, nel
                nddlt = 0
                numa = zzliel(ili, igr, iel)
!
                if (numa .gt. 0) then
!                   -- MAILLES DU MAILLAGE :
                    if (lmadis .and. ldist .and. .not. ldgrel) then
                        if (maille(numa) .ne. rang) goto 70
                    end if
!
                    nnoe = zznbne(numa)
                    do k1 = 1, nnoe
                        n1 = zzconx(numa, k1)
                        iad1 = zzprno(1, n1, 1)
                        nddl1 = zzprno(1, n1, 2)
                        if (mxddlt .lt. (nddlt+nddl1)) then
                            mxddlt = 2*(nddlt+nddl1)
                            call juveca('&&PROMOR.ANCIEN.LM', mxddlt)
                            call jeveuo('&&PROMOR.ANCIEN.LM', 'E', jalm)
                        end if
                        do iddl = 1, nddl1
                            zi4(jalm+nddlt+iddl-1) = int(zi(jnueq-1+iad1+iddl-1), kind=4)
                        end do
                        nddlt = nddlt+nddl1
                    end do
!
                else
!                   -- MAILLES TARDIVES :
                    if (lmadis .and. ldist .and. .not. ldgrel) then
                        if (rang .ne. 0) goto 70
                    end if
!
                    numa = -numa
                    nnoe = zznsup(ili, numa)
                    do k1 = 1, nnoe
                        n1 = zznema(ili, numa, k1)
                        if (n1 .lt. 0) then
                            n1 = -n1
                            if (lligrel_cp) then
                                n12 = v_crco(n1)
                            else
                                n12 = n1
                            end if
                            iad1 = zzprno(ili2, n12, 1)
                            nddl1 = zzprno(ili2, n12, 2)
                        else
                            iad1 = zzprno(1, n1, 1)
                            nddl1 = zzprno(1, n1, 2)
                        end if
                        if (mxddlt .lt. (nddlt+nddl1)) then
                            mxddlt = 2*(nddlt+nddl1)
                            call juveca('&&PROMOR.ANCIEN.LM', mxddlt)
                            call jeveuo('&&PROMOR.ANCIEN.LM', 'E', jalm)
                        end if
                        do iddl = 1, nddl1
                            zi4(jalm+nddlt+iddl-1) = int(zi(jnueq-1+iad1+iddl-1), kind=4)
                        end do
                        nddlt = nddlt+nddl1
                    end do
                end if
!
!       -- TRI EN ORDRE CROISSANT POUR L'INSERTION DES COLONNES
                ASSERT(nddlt .le. mxddlt)
                call uttrii_i4(zi4(jalm), nddlt)
!
!       -- INSERTION DES COLONNES DE L'ELEMENT DANS
!           LA STRUCTURE CHAINEE
                do iddl = 0, nddlt-1
                    position = zi4(jalm+iddl)-1
                    call vector_vector_add_values(vec_ptr, position, iddl+1, zi4(jalm))
                end do
70              continue
            end do
80          continue
        end do
!
!       3. TRAITEMENT DES SOUS-STRUCTURES STATIQUES :
!       ---------------------------------------------
90      continue
        if (nbss .gt. 0) then
            do ima = 1, nbsma
                if (sssa(ima) .eq. 0) goto 130
                nddlt = 0
                call jeveuo(jexnum(ma//'.SUPMAIL', ima), 'L', iamail)
                call jelira(jexnum(ma//'.SUPMAIL', ima), 'LONMAX', nnoe)
                do k1 = 1, nnoe
                    n1 = zi(iamail-1+k1)
                    ASSERT(n1 .ne. 0)
                    iad1 = zzprno(1, n1, 1)
                    nddl1 = zzprno(1, n1, 2)
                    if (mxddlt .lt. (nddlt+nddl1)) then
                        mxddlt = 2*(nddlt+nddl1)
                        call juveca('&&PROMOR.ANCIEN.LM', mxddlt)
                        call jeveuo('&&PROMOR.ANCIEN.LM', 'E', jalm)
                    end if
                    do iddl = 1, nddl1
                        zi4(jalm+nddlt+iddl-1) = int(zi(jnueq-1+iad1+iddl-1), kind=4)
                    end do
                    nddlt = nddlt+nddl1
                end do
!
                ASSERT(nddlt .le. mxddlt)
                call uttrii_i4(zi4(jalm), nddlt)
                do iddl = 0, nddlt-1
                    position = zi4(jalm+iddl)-1
                    call vector_vector_add_values(vec_ptr, position, iddl+1, zi4(jalm))
                end do
130             continue
            end do
        end if
    end do
!
!     DESIMBRIQUATION DE CHAINES POUR OBTENIR LA STRUCTURE COMPACTE
!     (SMDI,SMHC) DE LA MATRICE
    call vector_vector_sort_unique(vec_ptr)
    call vector_vector_get_data_sizes(vec_ptr, vecSize, ncoef)
    call wkvect(nu//'.SMOS.SMDI', base//' V I', nequ, vi=v_smdi)
    call wkvect(nu//'.SMOS.SMHC', base//' V S', ncoef, vi4=v_smhc)
    call vector_vector_get_data(vec_ptr, v_smdi, v_smhc)
    call vector_vector_delete(vec_ptr)
    v_smdi(1) = 1
!   Le contenu du SMDI ne correspond pas à la position du premier terme
!   d'une ligne donnée dans le SMHC mais à la position du terme diagonal
!   de la ligne dans le SMHC, c'est pourquoi on doit réaliser la boucle
!   suivante
    do iddl = 2, nequ-1
        v_smdi(iddl) = v_smdi(iddl+1)-1
    end do
    v_smdi(nequ) = ncoef
!
!     -- CREATION ET REMPLISSAGE DE .SMDE
    call wkvect(nu//'.SMOS.SMDE', base//' V I', 6, jsmde)
    zi(jsmde-1+1) = nequ
    zi(jsmde-1+2) = ncoef
    zi(jsmde-1+3) = 1
    call jedetr('&&PROMOR.ANCIEN.LM   ')
!
    if (niv .ge. 1 .and. printt .and. .not. isParallelMesh(ma)) then
        vali(1) = nequ
        vali(2) = ncoef
        vali(3) = 2*ncoef-nequ
        rcoef = ncoef
        requ = nequ
        valr(1) = (100.d0*(2.d0*rcoef-requ))/(requ*requ)
        call utmess('I', 'FACTOR_2', ni=3, vali=vali, sr=valr(1))
    end if
!
    call jedema()
!
end subroutine
