/**
 * @file FortranInterface.cxx
 * @brief Python bindings for Fortran interface.
 * @author Mathieu Courtois
 * @section LICENCE
 *   Copyright (C) 1991 - 2026  EDF www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "PythonBindings/FortranInterface.h"

#include "aster_mpi.h"
#include "aster_pybind.h"

void exportFortranToPython( py::module_ &mod ) {

    // These functions are for internal use.
    mod.def( "jeveux_init", &jeveux_init, R"(
Initialize the memory manager (Jeveux).

Arguments:
    mpi_comm (int): Identifier of MPI communicator (from ``py2f()``).
        )",
             py::arg( "mpi_comm" ) );

    mod.def( "jeveux_finalize", &jeveux_finalize, R"(
Finalize the memory manager (Jeveux).
        )" );

    mod.def( "call_oper", &call_oper, R"(
Call a Fortran operator ('op' subroutine).

Arguments:
    syntax (CommandSyntax): Object containing the user syntax.
    jxveri (int): If non null `JXVERI` is called after calling the operator.
        )",
             py::arg( "syntax" ), py::arg( "jxveri" ) );

    mod.def( "call_oper_init", &call_oper_init, R"(
Execute initializations before and after operator but without executing any
operator.
        )" );

    mod.def( "cmd_ctxt_enter", &call_cmd_ctxt_enter, R"(
Call Fortran 'cmd_ctxt_enter' subroutine.
        )" );

    mod.def( "cmd_ctxt_exit", &call_cmd_ctxt_exit, R"(
Call Fortran 'cmd_ctxt_exit' subroutine.
        )" );

    mod.def( "write", &call_print, R"(
Print a string using the fortran subroutine.

Arguments:
    text (str): Text to be printed.
        )",
             py::arg( "text" ) );

    mod.def( "affich", &call_affich, R"(
Print a string using the fortran subroutine on an internal file.

Arguments:
    code (str): Code name of the internal file : 'MESSAGE' or 'CODE'.
    text (str): Text to be printed.
        )",
             py::arg( "code" ), py::arg( "text" ) );

    mod.def( "jeveux_status", &get_sh_jeveux_status, R"(
Return the status of the Jeveux memory manager.

Returns:
    int: 0 after initialization and after shutdown, 1 during the execution.
        )" );

    mod.def( "jeveux_delete", &jeveux_delete, R"(
Force manual deletion of Jeveux objects.

Warning: Use only for debugging usages, it is dangerous for objects integrity
and cpu consuming.

Arguments:
    prefix (str): Root name of the Jeveux datastructure.
        )",
             py::arg( "prefix" ) );

    mod.def( "deleteTemporaryObjects", deleteTemporaryObjects, R"(
Delete temporary Jeveux objects
        )" );

    mod.def( "deleteCachedObjects", deleteCachedObjects, R"(
Delete temporary and cached Jeveux objects (temporary, matrix, base, ...)
        )" );

    mod.def( "onFatalError", &onFatalError,
             R"(
Get/set the behavior in case of error.

Arguments:
    value (str, optional): Set the new behavior in case of error (one of "ABORT",
        "EXCEPTION", "EXCEPTION+VALID" or "INIT"). If `value` is not provided,
        the current behavior is returned.

Returns:
    str: Current value
    )",
             py::arg( "value" ) = "" );

    mod.def( "matfpe", &call_matfpe,
             R"(
Enable or disable floating point exceptions.

Arguments:
    value (int): -1 to disable the FPE interception, 1 to enable FPE detection.
    )",
             py::arg( "value" ) );

    mod.def( "fe_invalid", &aster_fe_invalid,
             R"(
Enable or disable FE_INVALID exception.

Arguments:
    value (int): -1 to disable the interception, 1 to enable detection.
    )",
             py::arg( "value" ) );

    mod.def( "reset_tpmax", &reset_tpmax, R"(
Reset the cache value of tpmax (for internal use only).
        )" );

    mod.def( "asmpi_set", &asmpi_set, R"(
Set the current MPI communicator.

Arguments:
    comm (int): id of the communicator.
        )" );
    mod.def( "asmpi_get", &asmpi_get, R"(
Get the current MPI communicator.

Returns:
    comm (int): id of the communicator.
        )" );
    mod.def( "asmpi_free", &asmpi_free, R"(
Free the MPI communicator in argument.

Arguments:
    comm (int): id of the communicator.
        )" );
    mod.def( "asmpi_split", &asmpi_split, R"(
Split the MPI communicator in argument.

Arguments:
    comm (int): id of the parent communicator to split.
    color (int): color to which the calling process will belong.
    name (str): name of the new communicator.

Returns:
    comm (int) : id of the communicator.
        )" );
    mod.def( "asmpi_info", &asmpi_info, R"(
Return the rank and size of the MPI communicator.

Arguments:
    comm (int): id of the communicator.

Returns:
    rank (int) : rank of the communicator.
    size (int) : size of the communicator.
        )" );
};
