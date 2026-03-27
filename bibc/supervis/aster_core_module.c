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

/*
 * Module for the core functions of Code_Aster :
 *  - give access to the main object (for fortran calls) :
 *      - JDC : current JDC object,
 *      - CoreOptions : command line options + basic informations
 *      - MessageLog : utility to print message,
 *      - aster_core : pure python module that define functions more easily
 *        than here.
 *  - give informations about the execution.
 */
#include "aster.h"

#include "aster_core_module.h"

#include "aster_exceptions.h"
#include "aster_fort_mpi.h"
#include "aster_fort_utils.h"
#include "aster_module.h"
#include "aster_mpi.h"
#include "aster_utils.h"
#include "shared_vars.h"

#include <signal.h>

/*! aster_core C module */
static PyObject *aster_core = (PyObject *)0;

static char register_jdc_doc[] = "Enregistre des objets globaux.";

static PyObject *register_jdc( PyObject *self, PyObject *args ) {
    /*
     * Register the Python instances for usage from fortran/libaster
     */
    PyObject *val;
    PyObject *coreopts, *msglog;
    if ( !PyArg_ParseTuple( args, "OO:register_jdc", &coreopts, &msglog ) )
        return NULL;
    register_sh_params( coreopts );
    register_sh_msglog( msglog );

    // Add some wrappers for convenience
    val = PyObject_GetAttrString( coreopts, "get_option" );
    if ( PyObject_SetAttrString( aster_core, "get_option", val ) < 0 )
        MYABORT( "erreur lors de l'initialisation de 'aster_core.get_option'." );

    Py_DECREF( val );
    Py_INCREF( Py_None );
    return Py_None;
}

/*
 * Functions based on JDC object.
 */
ASTERINTEGER DEFS( JDCGET, jdcget, char *attr, STRING_SIZE l_attr ) {
    /*
     * Permet de récupérer la valeur entière d'un attribut du jdc.
     */
    PyObject *val;
    ASTERINTEGER value;

    val = PyObject_CallMethod( get_sh_params(), "get_option", "s#", attr, l_attr );
    if ( val == NULL ) {
        MYABORT( "erreur dans JDCGET" );
    }
    if ( !PyLong_Check( val ) )
        MYABORT( "Seuls les attributs de type entier peuvent etre recuperes !" );

    value = (ASTERINTEGER)PyLong_AsLong( val );

    Py_XDECREF( val );
    return value;
}

void DEFSP( JDCSET, jdcset, char *attr, STRING_SIZE l_attr, ASTERINTEGER *value ) {
    /*
     * Permet de positionner la valeur entière d'un attribut du jdc à `value`.
     */
    PyObject *res;

    res = PyObject_CallMethod( get_sh_params(), "set_option", "s#l", attr, l_attr, (long)*value );
    if ( res == NULL )
        MYABORT( "erreur dans JDCSET" );
    Py_XDECREF( res );
}

PyObject *GetJdcAttr( _IN char *attribut ) {
    /*
     * Retourne un attribut du 'jdc' en tant que PyObject.
     *
     * Return value: New reference.
     */
    PyObject *objattr;
    objattr = PyObject_GetAttrString( get_sh_params(), attribut );
    /* traiter l'erreur "objattr == NULL" dans l'appelant */
    return objattr;
}

static double _cache_tpmax = -1.;

void reset_tpmax() {
    /*
     * Reset le cache de tpmax.
     * La valeur est mise en cache pour éviter le passage au Python à chaque
     * appel de uttrst/uttcp0/uttcpu.
     */
    _cache_tpmax = -1;
}

double get_tpmax() {
    /*
     * Retourne le temps maximum autorisé pour l'exécution
     */
    int iret = 0;
    double tpmax;
    if ( _cache_tpmax < 0 ) {
        tpmax = asterc_getopt_double( "tpmax", &iret );
        if ( iret == 0 ) {
            _cache_tpmax = tpmax;
        }
    }
    return _cache_tpmax;
}

void DEFP( RDTMAX, rdtmax, _IN ASTERDOUBLE *tsub ) {
    /*
     * Réduit le temps maximum de l'exécution : tpmax = tpmax - tsub
     */
    PyObject *res;

    res = PyObject_CallMethod( get_sh_params(), "sub_tpmax", "d", (double)( *tsub ) );
    if ( res == NULL )
        MYABORT( "erreur dans RDTMAX" );
    // reset du cache
    reset_tpmax();
    Py_DECREF( res );
    return;
}

/*
 * Functions based on CoreOpts object.
 */
PyObject *asterc_getopt( _IN char *option ) {
    /*
     * Interface Fortran/Python pour récupérer une option de la ligne de commande.
     * Retourne :
     *  iret = 0 : tout est ok
     *  iret > 0   erreur
     *      iret = 1 : longueur de valk insuffisante, valeur tronquée
     *      iret = 4 : option inexistante, type incorrect.
     */
    PyObject *res;

    res = PyObject_CallMethod( get_sh_params(), "get_option", "s", option );
    if ( !res )
        MYABORT( "erreur lors de l'appel a la methode CoreOptions.get_option" );

    return res;
}

long asterc_getopt_long( _IN char *option, _OUT int *iret ) {
    /*
     * Interface C/Python pour récupérer une option de la ligne de commande.
     * Retourne la valeur et un code retour :
     *  iret = 0 : tout est ok
     *  iret = 4 : option inexistante, type incorrect.
     */
    PyObject *res;
    long value = 0;
    *iret = 4;
    res = asterc_getopt( option );
    if ( PyLong_Check( res ) ) {
        value = PyLong_AsLong( res );
        *iret = 0;
    } else if ( PyLong_Check( res ) ) {
        value = PyLong_AsLong( res );
        *iret = 0;
    }
    Py_DECREF( res );
    return value;
}

double asterc_getopt_double( _IN char *option, _OUT int *iret ) {
    /*
     * Interface C/Python pour récupérer une option de la ligne de commande.
     * Retourne la valeur et un code retour :
     *  iret = 0 : tout est ok
     *  iret = 4 : option inexistante, type incorrect.
     */
    PyObject *res;
    double value = 0.;
    *iret = 4;
    res = asterc_getopt( option );
    if ( PyFloat_Check( res ) ) {
        value = PyFloat_AsDouble( res );
        *iret = 0;
    } else if ( PyLong_Check( res ) ) {
        value = (double)PyLong_AsLong( res );
        *iret = 0;
    }
    Py_DECREF( res );
    return value;
}

char *asterc_getopt_string( _IN char *option, _OUT int *iret ) {
    /*
     * Interface C/Python pour récupérer une option de la ligne de commande.
     * Retourne la valeur et un code retour :
     *  iret = 0 : tout est ok
     *  iret = 4 : option inexistante, type incorrect.
     */
    PyObject *res;
    char *value = NULL;
    const char *stmp = NULL;
    STRING_SIZE lv;
    *iret = 4;
    res = asterc_getopt( option );
    if ( PyUnicode_Check( res ) ) {
        stmp = PyUnicode_AsUTF8( res );
        lv = strlen( stmp );
        value = MakeFStrFromCStr( stmp, strlen( stmp ) );
        *iret = 0;
    }
    Py_DECREF( res );
    return value;
}

void DEFSPP( GTOPTI, gtopti, _IN char *option, STRING_SIZE lopt, _OUT ASTERINTEGER *vali,
             _OUT ASTERINTEGER *iret ) {
    /*
     * Interface C/Python pour récupérer une option de la ligne de commande.
     * Retourne la valeur et un code retour :
     *  iret = 0 : tout est ok
     *  iret = 4 : option inexistante, type incorrect.
     */
    long value;
    int ret = 0;
    char *sopt;
    sopt = MakeCStrFromFStr( option, lopt );
    value = asterc_getopt_long( sopt, &ret );

    *vali = (ASTERINTEGER)value;
    *iret = (ASTERINTEGER)ret;
    FreeStr( sopt );
}

void DEFSPP( GTOPTR, gtoptr, _IN char *option, STRING_SIZE lopt, _OUT ASTERDOUBLE *valr,
             _OUT ASTERINTEGER *iret ) {
    /*
     * Interface C/Python pour récupérer une option de la ligne de commande.
     * Retourne la valeur et un code retour :
     *  iret = 0 : tout est ok
     *  iret = 4 : option inexistante, type incorrect.
     */
    double value;
    int ret = 0;
    char *sopt;
    sopt = MakeCStrFromFStr( option, lopt );
    value = asterc_getopt_double( sopt, &ret );

    *valr = (ASTERDOUBLE)value;
    *iret = (ASTERINTEGER)ret;
    FreeStr( sopt );
}

void DEFSSP( GTOPTK, gtoptk, _IN char *option, STRING_SIZE lopt, _OUT char *valk, STRING_SIZE lvalk,
             _OUT ASTERINTEGER *iret ) {
    /*
     * Interface C/Python pour récupérer une option de la ligne de commande.
     * Retourne la valeur et un code retour :
     *  iret = 0 : tout est ok
     *  iret = 1 : longueur de valk insuffisante, valeur tronquée
     *  iret = 4 : option inexistante, type incorrect.
     */
    char *value;
    int ret = 0;
    char *sopt;
    sopt = MakeCStrFromFStr( option, lopt );
    value = asterc_getopt_string( sopt, &ret );

    if ( ret == 0 ) {
        if ( strlen( value ) > lvalk )
            ret = 1;
        CopyCStrToFStr( valk, value, lvalk );
    }
    *iret = (ASTERINTEGER)ret;
    FreeStr( sopt );
    FreeStr( value );
}

static char get_mem_stat_doc[] = "Interface d'appel a la routine fortran UTGTME.\n";

static PyObject *asterc_get_mem_stat( PyObject *self, PyObject *args ) {
    /*
     *  Interface d'appel à la routine fortran UTGTME
     */
    PyObject *t_valres, *t_res;
    int inbpar;
    ASTERINTEGER nbpar, codret;
    char *nompar;
    ASTERDOUBLE *valres;
    /* doit impérativement correspondre aux longueurs des chaines de caractères
     * fortran */
    STRING_SIZE long_nompar = 8;
    void *malloc( size_t size );

    /* Conversion en tableaux de chaines */
    inbpar = (int)PyTuple_Size( args );
    nbpar = (ASTERINTEGER)inbpar;
    nompar = MakeTabFStr( inbpar, long_nompar );
    convertxt( inbpar, args, nompar, long_nompar );

    /* allocation des variables de sortie */
    valres = (ASTERDOUBLE *)malloc( nbpar * sizeof( ASTERDOUBLE ) );

    CALL_UTGTME( &nbpar, nompar, valres, &codret );

    t_valres = MakeTupleFloat( (long)inbpar, valres );

    /* retour de la fonction */
    t_res = PyTuple_New( 2 );
    PyTuple_SetItem( t_res, 0, t_valres );
    PyTuple_SetItem( t_res, 1, PyLong_FromLong( (long)codret ) );

    FreeStr( nompar );
    free( valres );

    return t_res;
}

/*
 * Functions based on MessageLog object.
 */
void DEFSPSPSPPPPS( UTPRIN, utprin, _IN char *typmess, _IN STRING_SIZE ltype,
                    _IN ASTERINTEGER *exc_typ, _IN char *idmess, _IN STRING_SIZE lidmess,
                    _IN ASTERINTEGER *nbk, _IN char *valk, _IN STRING_SIZE lvk,
                    _IN ASTERINTEGER *nbi, _IN ASTERINTEGER *vali, _IN ASTERINTEGER *nbr,
                    _IN ASTERDOUBLE *valr, _IN char *fname, _IN STRING_SIZE lfn ) {
    /*
     * Fortran/Python interface to print the messages.
     *
     * WARNING: In the case that the error indicator has already been set, we must
     * restore it after PyObject_CallMethod.
     */
    char test = ' ';
    char *kvar;
    int i, iret, iexc = 0;
    PyObject *args, *kwargs, *pyfname, *pFunc, *res;
    PyObject *tup_valk, *tup_vali, *tup_valr;
    PyObject *etype, *eval, *etb;

    if ( PyErr_Occurred() ) {
        iexc = 1;
        PyErr_Fetch( &etype, &eval, &etb );
    }

    tup_valk = PyTuple_New( (Py_ssize_t)*nbk );
    for ( i = 0; i < *nbk; i++ ) {
        kvar = valk + i * lvk;
        char *copyChaine = MakeCStrFromFStr( kvar, lvk );
        PyTuple_SetItem( tup_valk, i, PyUnicode_FromString( copyChaine ) );
        FreeStr( copyChaine );
    }

    tup_vali = PyTuple_New( (Py_ssize_t)*nbi );
    for ( i = 0; i < *nbi; i++ ) {
        PyTuple_SetItem( tup_vali, i, PyLong_FromLong( (long)vali[i] ) );
    }

    tup_valr = PyTuple_New( (Py_ssize_t)*nbr );
    for ( i = 0; i < *nbr; i++ ) {
        PyTuple_SetItem( tup_valr, i, PyFloat_FromDouble( (double)valr[i] ) );
    }

    args = Py_BuildValue( "s#s#OOOi", typmess, ltype, idmess, lidmess, tup_valk, tup_vali, tup_valr,
                          (int)*exc_typ );
    kwargs = PyDict_New();
    pyfname = PyUnicode_FromStringAndSize( fname, lfn );
    iret = PyDict_SetItemString( kwargs, "files", pyfname );
    if ( iret != 0 ) {
        MYABORT( "error the given filename in utprin" );
    }

    res = PyObject_Call( get_sh_msglog(), args, kwargs );
    if ( !res ) {
        MYABORT( "erreur lors de l'appel a MessageLog" );
    }
    if ( iexc == 1 ) {
        PyErr_Restore( etype, eval, etb );
    }

    Py_DECREF( pyfname );
    Py_DECREF( args );
    Py_XDECREF( kwargs );
    Py_DECREF( tup_valk );
    Py_DECREF( tup_vali );
    Py_DECREF( tup_valr );
    Py_DECREF( res );
}

void DEFPP( CHKMSG, chkmsg, _IN ASTERINTEGER *info_alarm, _OUT ASTERINTEGER *iret ) {
    /*
     * Interface Fortran/Python pour la vérification que tout s'est bien
     * passé, destinée à etre appelée dans FIN ou au cours d'une commande.
     * Argument IN :
     *    info_alarm = 1  on vérifie si les alarmes ignorées ont été émises ou
     * non. = 0  on ne fait pas cette vérif Retourne : iret = 0 : tout est ok iret
     * > 0   erreur
     */
    PyObject *res;

    res = PyObject_CallMethod( get_sh_msglog(), "check_counter", "i", (int)*info_alarm );
    if ( !res )
        MYABORT( "erreur lors de l'appel a la methode MessageLog.check_counter" );
    *iret = (ASTERINTEGER)PyLong_AsLong( res );

    Py_DECREF( res );
}

void DEFSS( UTALRM, utalrm, _IN char *is_bool, _IN STRING_SIZE lbool, _IN char *idmess,
            _IN STRING_SIZE lidm ) {
    /* Interface Fortran/Python pour masquer et rétablir l'affichage d'une alarme.
     *
     * call utalrm('OFF', 'CALCULEL5_7') == MasquerAlarme('CALCULEL5_7')
     * call utalrm('ON', 'CALCULEL5_7') == RetablirAlarme('CALCULEL5_7')
     */
    char *onoff, *s_id;
    PyObject *res;

    onoff = MakeCStrFromFStr( is_bool, lbool );
    s_id = MakeCStrFromFStr( idmess, lidm );

    if ( !strcmp( onoff, "OFF" ) ) {
        res = PyObject_CallMethod( get_sh_msglog(), "disable_alarm", "sO", s_id, Py_True );
    } else {
        res = PyObject_CallMethod( get_sh_msglog(), "reset_alarm", "sO", s_id, Py_True );
    }

    Py_DECREF( res );
    FreeStr( onoff );
    FreeStr( s_id );
}

void DEFP( GTALRM, gtalrm, _OUT ASTERINTEGER *nb ) {
    /* Interface Fortran/Python pour obtenir si des alarmes ont été émises.
     */
    PyObject *res;
    res = PyObject_CallMethod( get_sh_msglog(), "get_info_alarm_nb", "" );
    if ( !res )
        MYABORT( "erreur lors de l'appel a la methode 'get_info_alarm_nb'" );
    *nb = (ASTERINTEGER)PyLong_AsLong( res );
    Py_DECREF( res );
}

/*
 * Functions defined in aster_core
 */
void DEFP( PRHEAD, prhead, _IN ASTERINTEGER *part ) {
    /*
     * Interface Fortran/Python pour l'affichage des informations systèmes
     * en début d'exécution
     * Voir help(aster_core.print_header)
     */
    PyObject *res;
    PyObject *func = GetJdcAttr( "print_header" );
    res = PyObject_CallFunction( func, "i", (int)( *part ) );
    if ( !res )
        MYABORT( "erreur lors de l'appel a la fonction E_Global.print_header" );
    Py_DECREF( func );
    Py_DECREF( res );
}

void DEFSSP( CHEKSD, cheksd, _IN char *nomsd, _IN STRING_SIZE lnom, _IN char *typsd,
             _IN STRING_SIZE ltyp, _OUT ASTERINTEGER *iret ) {
    /*
        Interface Fortran/C pour vérifier que la structure de données `nomsd`
        est conforme au type `typsd`.

        Exemple d'appel :
            call cheksd('MA', 'sd_maillage', iret)
    */
    PyObject *res;
    PyObject *func = GetJdcAttr( "checksd" );
    res = PyObject_CallFunction( func, "s#s#", nomsd, lnom, typsd, ltyp );
    if ( !res )
        MYABORT( "erreur lors de l'appel a la methode CHECKSD" );
    *iret = (ASTERINTEGER)PyLong_AsLong( res );
    Py_DECREF( func );
    Py_DECREF( res );
}

void DEFSSPPPPPPPPPPPP( TESTRESU_PRINT, testresu_print, _IN char *refer, _IN STRING_SIZE lref,
                        _IN char *legend, _IN STRING_SIZE lleg, _IN ASTERINTEGER *llab,
                        _IN ASTERINTEGER *skip, _IN ASTERINTEGER *rela, _IN ASTERDOUBLE *tole,
                        _IN ASTERINTEGER *typ, _IN ASTERDOUBLE *refr, _IN ASTERDOUBLE *valr,
                        _IN ASTERINTEGER *refi, _IN ASTERINTEGER *vali, _IN ASTERDOUBLE *refc,
                        _IN ASTERDOUBLE *valc, _IN ASTERDOUBLE *compare ) {
    /*
        Interface Fortran/C pour imprimer le résultat d'un TEST_RESU

        def testresu_print(type_ref, legend, label, skip, relative,
                           tole, ref, val, compare=1.):
    */
    PyObject *res, *func, *args, *kwargs, *ref, *val, *comp = NULL;
    int ityp;

    func = PyObject_GetAttrString( get_sh_params(), "testresu_print" );
    /* positional arguments */
    args = Py_BuildValue( "s#s#llld", refer, lref, legend, lleg, (long)( *llab ), (long)( *skip ),
                          (long)( *rela ), (double)( *tole ) );
    /* keyword arguments */
    ityp = (int)( *typ );
    switch ( ityp ) {
    case 1:
        ref = PyFloat_FromDouble( (double)( *refr ) );
        val = PyFloat_FromDouble( (double)( *valr ) );
        break;
    case 2:
        ref = PyLong_FromLong( (long)( *refi ) );
        val = PyLong_FromLong( (long)( *vali ) );
        break;
    case 3:
        ref = PyComplex_FromDoubles( (double)( *refc ), (double)( *( refc + 1 ) ) );
        val = PyComplex_FromDoubles( (double)( *valc ), (double)( *( valc + 1 ) ) );
        break;
    }
    kwargs = PyDict_New();
    PyDict_SetItemString( kwargs, "ref", ref );
    PyDict_SetItemString( kwargs, "val", val );
    if ( (float)( *compare ) != 1. ) {
        comp = PyFloat_FromDouble( (double)( *compare ) );
        PyDict_SetItemString( kwargs, "compare", comp );
        Py_DECREF( comp );
    }

    res = PyObject_Call( func, args, kwargs );
    if ( !res ) {
        MYABORT( "erreur lors de l'appel a testresu_print" );
    }

    Py_DECREF( res );
    Py_XDECREF( func );
    Py_XDECREF( args );
    Py_XDECREF( kwargs );
    Py_XDECREF( ref );
    Py_XDECREF( val );
}

/*
 * Functions to communicate the execution status in parallel
 */

static PyObject *aster_mpi_warn( PyObject *self, PyObject *args ) {
    ASTERINTEGER iexc = 1;
    /* call ASMPI_WARN */
    try {
        CALL_ASMPI_WARN( &iexc );
    }
    exceptAll { raiseException(); }
    endTry();
    Py_INCREF( Py_None );
    return Py_None;
}

/*
 * Methods of the aster_core module.
 */
static PyMethodDef methods[] = { { "register", register_jdc, METH_VARARGS, register_jdc_doc },
                                 { "get_mem_stat", asterc_get_mem_stat, METH_VARARGS,
                                   get_mem_stat_doc },
                                 { "MPI_Warn", aster_mpi_warn, METH_VARARGS },
                                 // { "get_option",  ... } : method added in register_jdc
                                 { NULL, NULL, 0, NULL } };

static struct PyModuleDef aster_core_def = {
    PyModuleDef_HEAD_INIT, "aster_core", NULL, -1, methods, NULL, NULL, NULL, NULL
};

#ifndef ASTER_WITHOUT_PYMOD
PyObject *PyInit_aster_core( void ) {
    aster_core = PyModule_Create( &aster_core_def );
    return aster_core;
}
#endif
