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

#ifndef ASTER_FORT_UTILS_H_
#define ASTER_FORT_UTILS_H_

#include "aster.h"

/* ******************************************************
 *
 * Interfaces of fortran subroutines called from C/C++.
 *
 * ******************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#define CALL_AFFICH( a, b ) CALLOO( AFFICH, affich, a, b )
extern void DEFSS( AFFICH, affich, const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALL_CODENT( a, b, c ) CALLPSS( CODENT, codent, a, b, c )
#define CALLO_CODENT( a, b, c ) CALLPOO( CODENT, codent, a, b, c )
extern void DEFPSS( CODENT, codent, ASTERINTEGER *, const char *, STRING_SIZE, const char *,
                    STRING_SIZE );

#define CALL_CODLET_WRAP( a, b, c, d ) CALLPSSS( CODLET_WRAP, codlet_wrap, a, b, c, d )
#define CALLO_CODLET_WRAP( a, b, c, d ) CALLPOOO( CODLET_WRAP, codlet_wrap, a, b, c, d )
extern void DEFPSSS( CODLET_WRAP, codlet_wrap, ASTERINTEGER *, const char *, STRING_SIZE,
                     const char *, STRING_SIZE, const char *, STRING_SIZE );

#define CALL_DISMOI( a, b, c, d, e, f, g ) CALLSSSPSSP( DISMOI, dismoi, a, b, c, d, e, f, g )
#define CALLO_DISMOI( a, b, c, d, e, f, g ) CALLOOOPOOP( DISMOI, dismoi, a, b, c, d, e, f, g )
extern void DEFSSSPSSP( DISMOI, dismoi, const char *, STRING_SIZE, const char *, STRING_SIZE,
                        const char *, STRING_SIZE, ASTERINTEGER *, const char *, STRING_SIZE,
                        const char *, STRING_SIZE, ASTERINTEGER * );

#define CALLO_CARCHA( a, b, c, d, e ) CALLOOOOO( CARCHA, carcha, a, b, c, d, e )
extern void DEFSSSSS( CARCHA, carcha, const char *, STRING_SIZE, const char *, STRING_SIZE,
                      const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                      STRING_SIZE );

#define CALL_FCLOSE( a ) CALLP( FCLOSE, fclose, a )
extern void DEFP( FCLOSE, fclose, ASTERINTEGER * );

#define CALL_INFMAJ_EXT( a ) CALLP( INFMAJ_EXT, infmaj_ext, a )
extern void DEFP( INFMAJ_EXT, infmaj_ext, const ASTERINTEGER *const );

#define CALL_ISDECO( a, b, c ) CALLPPP( ISDECO, isdeco, a, b, c )
void DEFPPP( ISDECO, isdeco, ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER * );

#define CALL_ISNNEM() CALL0( ISNNEM, isnnem )
extern ASTERINTEGER DEF0( ISNNEM, isnnem );

#define CALL_MATFPE( a ) CALLP( MATFPE, matfpe, a )
extern void DEFP( MATFPE, matfpe, ASTERINTEGER * );

#define CALL_NBEC( a ) CALLP( NBEC, nbec, a )
ASTERINTEGER DEFP( NBEC, nbec, const ASTERINTEGER *const );

#define CALL_R8PI() CALL0( R8PI, r8pi )
extern ASTERDOUBLE DEF0( R8PI, r8pi );

#define CALL_R8VIDE() CALL0( R8VIDE, r8vide )
extern ASTERDOUBLE DEF0( R8VIDE, r8vide );

#define CALL_ULOPEN( a, b, c, d, e ) CALLPSSSS( ULOPEN, ulopen, a, b, c, d, e )
extern void DEFPSSSS( ULOPEN, ulopen, ASTERINTEGER *, char *, STRING_SIZE, char *, STRING_SIZE,
                      char *, STRING_SIZE, char *, STRING_SIZE );

#define CALL_UTGTME( a, b, c, d ) CALLPSPP( UTGTME, utgtme, a, b, c, d )
extern void DEFPSPP( UTGTME, utgtme, ASTERINTEGER *, char *, STRING_SIZE, ASTERDOUBLE *,
                     ASTERINTEGER * );

#define CALL_UTNCMP( a, b, c ) CALLSPS( UTNCMP, utncmp, a, b, c )
extern void DEFSPS( UTNCMP, utncmp, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
                    STRING_SIZE );

#define CALL_UTNCMP2( a, b, c, d ) CALLSPSS( UTNCMP2, utncmp2, a, b, c, d )
extern void DEFSPSS( UTNCMP2, utncmp2, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
                     STRING_SIZE, const char *, STRING_SIZE );

#define CALL_UTNCMP3( a, b, c, d ) CALLSPSS( UTNCMP3, utncmp3, a, b, c, d )
extern void DEFSPSS( UTNCMP3, utncmp3, const char *, STRING_SIZE, ASTERINTEGER *, const char *,
                     STRING_SIZE, const char *, STRING_SIZE );

#define CALL_UTMESS( cod, idmess ) CALLSS( UTMESS_CWRAP, utmess_cwrap, cod, idmess )
extern void DEFSS( UTMESS_CWRAP, utmess_cwrap, char *, STRING_SIZE, char *, STRING_SIZE );

/* particulier car on fixe les longueurs des chaines valk */
#define VALK_SIZE 128
extern void DEFSSPSPPPPPS( UTMESS_CORE, utmess_core, char *, STRING_SIZE, char *, STRING_SIZE,
                           ASTERINTEGER *, char *, STRING_SIZE, ASTERINTEGER *, ASTERINTEGER *,
                           ASTERINTEGER *, ASTERDOUBLE *, ASTERINTEGER *, char *, STRING_SIZE );
#ifdef ASTER_STRLEN_AT_END
#define CALL_UTMESS_CORE( cod, idmess, nk, valk, ni, vali, nr, valr, nexcep, fname )               \
    F_FUNC( UTMESS_CORE, utmess_core )                                                             \
    ( cod, idmess, nk, valk, ni, vali, nr, valr, nexcep, fname, strlen( cod ), strlen( idmess ),   \
      VALK_SIZE, strlen( fname ) )
#else
#define CALL_UTMESS_CORE( cod, idmess, nk, valk, ni, vali, nr, valr, nexcep, fname )               \
    F_FUNC(UTMESS_CORE, utmess_core)(cod, strlen(cod), idmess, strlen(idmess), nk, \
                                     valk, VALK_SIZE, ni, vali, nr, valr,
                                     nexcep, fname, strlen(fname))
#endif

#define CALL_UTPTME( a, b, c ) CALLSPP( UTPTME, utptme, a, b, c )
extern void DEFSPP( UTPTME, utptme, char *, STRING_SIZE, ASTERDOUBLE *, ASTERINTEGER * );

#define CALL_UTTCSM( a ) CALLP( UTTCSM, uttcsm, a )
extern void DEFP( UTTCSM, uttcsm, ASTERDOUBLE * );

#define CALL_VERIGD_WRAP( a, b, c, d ) CALLSSPP( VERIGD_WRAP, verigd_wrap, a, b, c, d )
#define CALLO_VERIGD_WRAP( a, b, c, d ) CALLOOPP( VERIGD_WRAP, verigd_wrap, a, b, c, d )
void DEFSSPP( VERIGD_WRAP, verigd_wrap, const char *, STRING_SIZE, const char *, STRING_SIZE,
              ASTERINTEGER *, ASTERINTEGER * );

#define CALL_UTEREF2( a, b, c, d, e, f, g, h, i, l, m, n, o, p, q )                                \
    CALLPPPSSSPPPPPPPPP( UTEREF2, uteref2, a, b, c, d, e, f, g, h, i, l, m, n, o, p, q )
#define CALLO_UTEREF2( a, b, c, d, e, f, g, h, i, l, m, n, o, p, q )                               \
    CALLPPPOOOPPPPPPPPP( UTEREF2, uteref2, a, b, c, d, e, f, g, h, i, l, m, n, o, p, q )
extern void DEFPPPSSSPPPPPPPPP( UTEREF2, uteref2, ASTERINTEGER *, ASTERINTEGER *, ASTERLOGICAL *,
                                const char *, STRING_SIZE, const char *, STRING_SIZE, const char *,
                                STRING_SIZE, ASTERINTEGER *, ASTERINTEGER *, ASTERINTEGER *,
                                ASTERINTEGER *, ASTERINTEGER *, ASTERDOUBLE *, ASTERDOUBLE *,
                                ASTERDOUBLE *, ASTERINTEGER * );

/* intrinsics */
#define CALL_ABORTF() CALL0( ABORTF, abortf )
extern void DEF0( ABORTF, abortf );

#ifdef __cplusplus
}
#endif

/* FIN ASTER_FORT_UTILS_H */
#endif
