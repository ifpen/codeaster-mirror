/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2026 - EDF - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#ifndef ASTER_FORT_CALCUL_H_
#define ASTER_FORT_CALCUL_H_

#include "aster.h"

/* ******************************************************
 *
 * Interfaces of fortran subroutines called from C/C++.
 *
 * ******************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#define CALL_ADDMATRASSE( a, b, c, d, e ) CALLOOPPO( ADDMATRASSE, addmatrasse, a, b, c, d, e )
void DEFSSPPS( ADDMATRASSE, addmatrasse, const char *, STRING_SIZE, const char *, STRING_SIZE,
               const ASTERDOUBLE *, const ASTERDOUBLE *, const char *, STRING_SIZE );

#define CALLO_ASASVE( a, b, c, d, e ) CALLOOOOO( ASASVE, asasve, a, b, c, d, e )
void DEFSSSSS( ASASVE, asasve, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
               STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_ASCOVA( a, b, c, d, e, f, g, h )                                                     \
    CALLOOOOPOOO( ASCOVA, ascova, a, b, c, d, e, f, g, h )
void DEFSSSSPSSS( ASCOVA, ascova, const char *, STRING_SIZE, const char *, STRING_SIZE,
                  const char *, STRING_SIZE, const char *, STRING_SIZE, const ASTERDOUBLE *,
                  const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_ASCAVC_WRAP( a, b, c, d, e, f )                                                      \
    CALLOOOPOO( ASCAVC_WRAP, ascavc_wrap, a, b, c, d, e, f )
void DEFSSSPSS( ASCAVC_WRAP, ascavc_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                const char *, STRING_SIZE, const ASTERDOUBLE *, const char *, STRING_SIZE,
                const char *, STRING_SIZE );

#define CALL_ASMATR( a, b, c, d, e, f, g, h, i )                                                   \
    CALLPSSSSSSPS( ASMATR, asmatr, a, b, c, d, e, f, g, h, i )
void DEFPSSSSSSPS( ASMATR, asmatr, ASTERINTEGER *, const char *, STRING_SIZE, const char *,
                   STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                   STRING_SIZE, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
                   STRING_SIZE );

#define CALL_ASSVEC( a, b, c, d, e, f, g )                                                         \
    CALLSSPSPSP( ASSVEC_WRAP, assvec_wrap, a, b, c, d, e, f, g )
void DEFSSPSPSP( ASSVEC_WRAP, assvec_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 ASTERINTEGER *, const char *, STRING_SIZE, ASTERDOUBLE *, const char *,
                 STRING_SIZE, ASTERINTEGER * );

#define CALL_ASSVECWITHMASK( a, b, c, d, e, f, g, h, i )                                           \
    CALLSSPSPSPSP( ASSVECWITHMASK, assvecwithmask, a, b, c, d, e, f, g, h, i )
void DEFSSPSPSPSP( ASSVECWITHMASK, assvecwithmask, const char *, STRING_SIZE, const char *,
                   STRING_SIZE, ASTERINTEGER *, const char *, STRING_SIZE, ASTERDOUBLE *,
                   const char *, STRING_SIZE, ASTERINTEGER *, const char *, STRING_SIZE,
                   const ASTERLOGICAL * );

#define CALL_ASSMIV( a, b, c, d, e, f, g ) CALLSSPSPSP( ASSMIV, assmiv, a, b, c, d, e, f, g )
void DEFSSPSPSP( ASSMIV, assmiv, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 ASTERINTEGER *, const char *, STRING_SIZE, ASTERDOUBLE *, const char *,
                 STRING_SIZE, ASTERINTEGER * );

#define CALLO_AP_ASSEMBLY_VECTOR( a ) CALLO( AP_ASSEMBLY_VECTOR, ap_assembly_vector, a )
void DEFS( AP_ASSEMBLY_VECTOR, ap_assembly_vector, const char *, STRING_SIZE );

#define CALLO_CONLAG( a, b ) CALLOP( CONLAG, conlag, a, b )
void DEFSP( CONLAG, conlag, const char *, STRING_SIZE, ASTERDOUBLE * );

#define CALLO_CORICHWRITE( a, b ) CALLOP( CORICHWRITE, corichwrite, a, b )
void DEFSP( CORICHWRITE, corichwrite, const char *, STRING_SIZE, const ASTERINTEGER * );

#define CALLO_CRESOL_WRAP( a, b, c ) CALLOOO( CRESOL_WRAP, cresol_wrap, a, b, c )
void DEFSSS( CRESOL_WRAP, cresol_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
             const char *, STRING_SIZE );

#define CALLO_CHPCHD( a, b, c, e, f, h, g ) CALLOOOOOOO( CHPCHD, chpchd, a, b, c, e, f, h, g )
void DEFSSSSSSS( CHPCHD, chpchd, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                 STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                 STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_MATRIX_FACTOR( a, b, c, d, e, f, g )                                                 \
    CALLOOPOOPP( MATRIX_FACTOR, matrix_factor, a, b, c, d, e, f, g )
void DEFSSPSSPP( MATRIX_FACTOR, matrix_factor, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 ASTERINTEGER *, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 ASTERINTEGER *, ASTERINTEGER * );

#define CALL_MATR_ASSE_SYME( a ) CALLO( MATR_ASSE_SYME, matr_asse_syme, a )
void DEFS( MATR_ASSE_SYME, matr_asse_syme, const char *, STRING_SIZE );

#define MATR_ASSE_COMPUTE_KINEMATIC_RHS( a, b, c )                                                 \
    CALLOOO( MATR_ASSE_COMPUTE_KINEMATIC_RHS, matr_asse_compute_kinematic_rhs, a, b, c )
void DEFSSS( MATR_ASSE_COMPUTE_KINEMATIC_RHS, matr_asse_compute_kinematic_rhs, const char *,
             STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_NMDOCH_WRAP( a, b ) CALLOO( NMDOCH_WRAP, nmdoch_wrap, a, b )
void DEFSS( NMDOCH_WRAP, nmdoch_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_NTDOCH_WRAP( a, b ) CALLOO( NTDOCH_WRAP, ntdoch_wrap, a, b )
void DEFSS( NTDOCH_WRAP, ntdoch_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_ACDOCH_WRAP( a, b ) CALLOO( ACDOCH_WRAP, acdoch_wrap, a, b )
void DEFSS( ACDOCH_WRAP, acdoch_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_NUMCIMA( a, b, c, d ) CALLOOOO( NUMCIMA, numcima, a, b, c, d )
void DEFSSSS( NUMCIMA, numcima, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
              STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_NUMERO_WRAP( a, b, c, d, e, f )                                                      \
    CALLOOOOOP( NUMERO_WRAP, numero_wrap, a, b, c, d, e, f )
void DEFSSSSSP( NUMERO_WRAP, numero_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                ASTERLOGICAL * );

#define CALLO_NUMER3_WRAP( a, b, c, d, e, f, g )                                                   \
    CALLOOOOOOP( NUMER3_WRAP, numer3_wrap, a, b, c, d, e, f, g )
void DEFSSSSSSP( NUMER3_WRAP, numer3_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 const char *, STRING_SIZE, ASTERLOGICAL * );

#define CALLO_NUME_DDL_MATR( a, b, c, d ) CALLOOPP( NUME_DDL_MATR, nume_ddl_matr, a, b, c, d )
void DEFSSPP( NUME_DDL_MATR, nume_ddl_matr, const char *, STRING_SIZE, const char *, STRING_SIZE,
              ASTERINTEGER *, ASTERLOGICAL * );

#define CALLO_NUME_DDL_CHAMELEM( a, b, c, d )                                                      \
    CALLOOOP( NUME_DDL_CHAMELEM, nume_ddl_chamelem, a, b, c, d )
void DEFSSSP( NUME_DDL_CHAMELEM, nume_ddl_chamelem, const char *, STRING_SIZE, const char *,
              STRING_SIZE, const char *, STRING_SIZE, ASTERLOGICAL * );

#define CALL_MVMULT( a, b, c ) CALLOOO( MVMULT, mvmult, a, b, c )
void DEFSSS( MVMULT, mvmult, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
             STRING_SIZE );

#define CALL_REDETR( a ) CALLO( REDETR, redetr, a )
void DEFS( REDETR, redetr, const char *, STRING_SIZE );

#define CALLO_RESOUD( a, b, c, d, e, f, g, h, i, j, k, l, m, n )                                   \
    CALLOOOOPOOOPPOPPP( RESOUD, resoud, a, b, c, d, e, f, g, h, i, j, k, l, m, n )
void DEFSSSSPSSSPPSPPP( RESOUD, resoud, const char *, STRING_SIZE, const char *, STRING_SIZE,
                        const char *, STRING_SIZE, const char *, STRING_SIZE, ASTERINTEGER *,
                        const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                        STRING_SIZE, ASTERDOUBLE *, ASTERDOUBLE *, const char *, STRING_SIZE,
                        ASTERLOGICAL *, ASTERINTEGER *, ASTERINTEGER * );

#define CALLO_SDMPIC( a, b ) CALLOO( SDMPIC, sdmpic, a, b )
void DEFSS( SDMPIC, sdmpic, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VECHME_WRAP( a, b, c, d, e, f, g, h, i, l, m, n )                                    \
    CALLOOOOPPPOOOOO( VECHME_WRAP, vechme_wrap, a, b, c, d, e, f, g, h, i, l, m, n )
void DEFSSSSPPPSSSSS( VECHME_WRAP, vechme_wrap, const char *, STRING_SIZE, const char *,
                      STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                      const ASTERDOUBLE *, const ASTERDOUBLE *, const ASTERDOUBLE *, const char *,
                      STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                      const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_RSINCHPRE( a, b, c, d ) CALLOOOP( RSINCHPRE, rsinchpre, a, b, c, d )
void DEFSSSP( RSINCHPRE, rsinchpre, const char *, STRING_SIZE, const char *, STRING_SIZE,
              const char *, STRING_SIZE, const ASTERINTEGER * );

#define CALLO_RSINCH( a, b, c, d, e, f, g, h, i, l, m, n )                                         \
    CALLOOOPOOOPOPOP( RSINCH, rsinch, a, b, c, d, e, f, g, h, i, l, m, n )
void DEFSSSPSSSPSPSP( RSINCH, rsinch, const char *, STRING_SIZE, const char *, STRING_SIZE,
                      const char *, STRING_SIZE, const ASTERDOUBLE *, const char *, STRING_SIZE,
                      const char *, STRING_SIZE, const char *, STRING_SIZE, const ASTERINTEGER *,
                      const char *, STRING_SIZE, const ASTERDOUBLE *, const char *, STRING_SIZE,
                      const ASTERINTEGER * );

#define CALLO_VEDIME( a, b, c, d, e, f ) CALLOOOPOO( VEDIME, vedime, a, b, c, d, e, f )
void DEFSSSPSS( VEDIME, vedime, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                STRING_SIZE, const ASTERDOUBLE *, const char *, STRING_SIZE, const char *,
                STRING_SIZE );

#define CALLO_VEBUME( a, b, c, d, e, f ) CALLOOOOPO( VEBUME, vebume, a, b, c, d, e, f )
void DEFSSSSPS( VEBUME, vebume, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                STRING_SIZE, const char *, STRING_SIZE, const ASTERDOUBLE *, const char *,
                STRING_SIZE );

#define CALLO_VELAME( a, b, c, d, e ) CALLOOOOO( VELAME, velame, a, b, c, d, e )
void DEFSSSSS( VELAME, velame, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
               STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VRCINS_WRAP( a, b, c, d, e, f, g )                                                   \
    CALLOOOPOOO( VRCINS_WRAP, vrcins_wrap, a, b, c, d, e, f, g )
void DEFSSSPSSS( VRCINS_WRAP, vrcins_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
                 const char *, STRING_SIZE, const ASTERDOUBLE *, const char *, STRING_SIZE,
                 const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VRCREF( a, b, c, d, e ) CALLOOOOO( VRCREF, vrcref, a, b, c, d, e )
void DEFSSSSS( VRCREF, vrcref, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
               STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VTCREB_WRAP( a, b, c, d ) CALLOOOO( VTCREB_WRAP, vtcreb_wrap, a, b, c, d )
void DEFSSSS( VTCREB_WRAP, vtcreb_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
              const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_VTCOPY( a, b, c ) CALLOOP( VTCOPY, vtcopy, a, b, c )
void DEFSSP( VTCOPY, vtcopy, const char *, STRING_SIZE, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_NMDOCC( a, b, c, d, e, f ) CALLOOPOOP( NMDOCC, nmdocc, a, b, c, d, e, f )
void DEFSSPSSP( NMDOCC, nmdocc, const char *, STRING_SIZE, const char *, STRING_SIZE,
                ASTERLOGICAL *, const char *, STRING_SIZE, const char *, STRING_SIZE,
                ASTERLOGICAL * );

#define CALLO_NXDOCC( a, b, c ) CALLOOO( NXDOCC, nxdocc, a, b, c )
void DEFSSS( NXDOCC, nxdocc, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
             STRING_SIZE );

#define CALLO_NMDOCR( a, b, c ) CALLOOO( NMDOCR, nmdocr, a, b, c )
void DEFSSS( NMDOCR, nmdocr, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
             STRING_SIZE );

#define CALLO_NMDOCM( a, b, c ) CALLOOO( NMDOCM, nmdocm, a, b, c )
void DEFSSS( NMDOCM, nmdocm, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
             STRING_SIZE );

#define CALLO_AFVARC( a, b, c ) CALLOOO( AFVARC, afvarc, a, b, c )
void DEFSSS( AFVARC, afvarc, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
             STRING_SIZE );

#define CALLO_CMTREF( a, b ) CALLOO( CMTREF, cmtref, a, b )
void DEFSS( CMTREF, cmtref, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_CESVAR( a, b, c, d ) CALLOOOO( CESVAR, cesvar, a, b, c, d )
void DEFSSSS( CESVAR, cesvar, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
              STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_ALCHML( a, b, c, d, e, f, g ) CALLOOOOOPO( ALCHML, alchml, a, b, c, d, e, f, g )
void DEFSSSSSPS( ALCHML, alchml, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                 STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE, ASTERINTEGER *,
                 const char *, STRING_SIZE );

#define CALLO_CALCUL( a, b, c, d, e, f, g, h, i, j, k )                                            \
    CALLSSSPSSPSSSS( CALCUL_CWRAP, calcul_cwrap, a, b, c, d, e, f, g, h, i, j, k )
void DEFSSSPSSPSSSS( CALCUL_CWRAP, calcul_cwrap, const char *, STRING_SIZE, const char *,
                     STRING_SIZE, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
                     STRING_SIZE, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
                     STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                     const char *, STRING_SIZE );

#define CALLO_CHECKSUPERELEMENT( a, b ) CALLOO( CHECKSUPERELEMENT, checksuperelement, a, b )
void DEFSS( CHECKSUPERELEMENT, checksuperelement, const char *, STRING_SIZE, const char *,
            STRING_SIZE );

#define CALLO_GETERRORCODE( a, b ) CALLOP( GETERRORCODE, geterrorcode, a, b )
void DEFSP( GETERRORCODE, geterrorcode, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_SS2MME( a, b, c ) CALLOOO( SS2MME, ss2mme, a, b, c )
void DEFSSS( SS2MME, ss2mme, const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
             STRING_SIZE );

#define CALLO_ME2MME_EVOL( a, b, c, d, e, f, g, h, i, j, k, l, m, n )                              \
    CALLOOOOPOPOOPPPOO( ME2MME_EVOL, me2mme_evol, a, b, c, d, e, f, g, h, i, j, k, l, m, n )
void DEFSSSSPSPSSPPPSS( ME2MME_EVOL, me2mme_evol, const char *, STRING_SIZE, const char *,
                        STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE,
                        ASTERINTEGER *, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
                        STRING_SIZE, const char *, STRING_SIZE, ASTERDOUBLE *, ASTERDOUBLE *,
                        ASTERDOUBLE *, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALLO_HASBEHAVIOURFEATURE( a, b, c, d, e )                                                 \
    CALLOOOOP( HASBEHAVIOURFEATURE, hasbehaviourfeature, a, b, c, d, e )
void DEFSSSSP( HASBEHAVIOURFEATURE, hasbehaviourfeature, const char *, STRING_SIZE, const char *,
               STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE, ASTERLOGICAL * );

#define CALL_GRANDEURID( a, b, c ) CALLOOP( GRANDEURID, grandeurid, a, b, c )
void DEFSSP( GRANDEURID, grandeurid, const char *, STRING_SIZE, const char *, STRING_SIZE,
             ASTERINTEGER * );

#define CALL_GET_LOCAL_MODE( a, b, c, d ) CALLOOPP( GET_LOCAL_MODE, get_local_mode, a, b, c, d )
void DEFSSPP( GET_LOCAL_MODE, get_local_mode, const char *, STRING_SIZE, const char *, STRING_SIZE,
              const ASTERINTEGER *, ASTERINTEGER * );

#define CALL_DIGDEL_WRAP( a, b ) CALLPP( DIGDEL_WRAP, digdel_wrap, a, b )
void DEFPP( DIGDEL_WRAP, digdel_wrap, const ASTERINTEGER *, ASTERINTEGER * );

#define CALL_CHCKVARI( a, b, c, d ) CALLOOOO( CHCKVARI, chckvari, a, b, c, d )
void DEFSSSS( CHCKAVRI, chckvari, const char *, STRING_SIZE, const char *, STRING_SIZE,
              const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALL_COMPAREFIELDSHAPE( a, b, c, d, e )                                                    \
    CALLOOPOP( COMPAREFIELDSHAPE, comparefieldshape, a, b, c, d, e )
void DEFSSPSP( COMPAREFIELDSHAPE, comparefieldshape, const char *, STRING_SIZE, const char *,
               STRING_SIZE, ASTERLOGICAL *, const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_VARINONU_WRAP( a, b, c, d, e )                                                       \
    CALLOOOOO( VARINONU_WRAP, varinonu_wrap, a, b, c, d, e )
void DEFSSSSS( VARINONU_WRAP, varinonu_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
               const char *, STRING_SIZE, const char *, STRING_SIZE, const char *, STRING_SIZE );

#ifdef __cplusplus
}
#endif

/* FIN ASTER_FORT_CALCUL_H */
#endif
