# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

import os
from pathlib import PureWindowsPath
from subprocess import PIPE, Popen

from waflib import Configure, Errors, Logs, Utils


def options(self):
    self.load("python")  # --nopyc/--nopyo are enabled below


def configure(self):
    self.check_python()
    self.check_numpy()
    self.check_asrun()
    self.check_mpi4py()


###############################################################################
@Configure.conf
def check_python(self):
    self.load("python")
    self.check_python_version((3, 5, 0))
    self.check_python_headers()
    if self.env.CC_IS_INTEL:
        self.env["LIB_PYEMBED"] = list(set(self.env["LIB_PYEMBED"]))
        # Best is to clear PYEMBED and PYEXT {c/cxx}flags
        for lang in ("CFLAGS", "CXXFLAGS"):
            for feat in ("PYEMBED", "PYEXT"):
                self.env[lang + "_" + feat] = []
    for lang in ("CFLAGS", "CXXFLAGS", "LINKFLAGS"):
        for feat in ("PYEMBED", "PYEXT"):
            self.remove_flags(lang + "_" + feat, "^-O[0-9]")
            self.remove_flags(lang + "_" + feat, "^-flto")
            self.remove_flags(lang + "_" + feat, "^-ffat\-lto")
    cfgext = self.env["CONFIG_PARAMETERS"].get("cfgext", "yaml")
    if cfgext in ("", "yaml"):
        try:
            self.check_python_module("yaml")
        except Errors.ConfigurationError:
            cfgext = "json"
    self.env["CFG_EXT"] = cfgext


@Configure.conf
def check_numpy(self):
    if not self.env["PYTHON"]:
        self.fatal("load python tool first")
    self.check_python_module("numpy")
    self.check_numpy_headers()


@Configure.conf
def check_numpy_headers(self):
    if not self.env["PYTHON"]:
        self.fatal("load python tool first")
    self.start_msg("Checking for numpy includes")
    # retrieve includes dir from numpy module
    cmd = self.env.PYTHON + ["-c", "\nimport numpy\nprint(numpy.get_include())"]
    numpy_includes = self.cmd_and_log(cmd, shell=False).strip()
    self.end_msg(numpy_includes)
    self.start_msg("Checking for numpy arrayobject.h")
    # Bad path formating on msys2
    if self.is_defined("ASTER_PLATFORM_MINGW") and not self.is_defined("ASTER_PLATFORM_MSYS2"):
        incs = PureWindowsPath(numpy_includes)
        parts = list(incs.parts)
        if incs.anchor:
            parts[0] = incs.root
        for i, sub in enumerate(parts):
            if sub == "lib":
                parts[i] = "Lib"
        numpy_includes = PureWindowsPath(*parts).as_posix()
    # check the given includes dirs
    self.check(
        feature="c",
        header_name="Python.h numpy/arrayobject.h",
        includes=[numpy_includes],
        use="PYEMBED",
        uselib_store="NUMPY",
        errmsg="Could not find the numpy development headers",
    )
    self.end_msg(numpy_includes)


@Configure.conf
def check_asrun(self):
    if not self.env["PYTHON"]:
        self.fatal("load python tool first")
    try:
        self.check_python_module("asrun")
    except Errors.WafError:
        # optional
        pass


@Configure.conf
def check_mpi4py(self):
    if not self.env.BUILD_MPI:
        return
    if not self.env["PYTHON"]:
        self.fatal("load python tool first")
    try:
        self.check_python_module("mpi4py")
    except Errors.ConfigurationError:
        if self.options.with_py_mpi4py:
            raise


@Configure.conf
def check_optimization_python(self):
    self.setenv("debug")
    self.env["PYC"] = self.env["PYO"] = 0
    self.setenv("release")
    self.env["PYC"] = self.env["PYO"] = 0


def _get_default_pythonpath(python):
    """Default sys.path should be added into PYTHONPATH"""
    env = os.environ.copy()
    env["PYTHONPATH"] = ""
    proc = Popen([python, "-c", "import sys; print(sys.path)"], stdout=PIPE, env=env)
    system_path = eval(proc.communicate()[0])
    return system_path


check_cfg_old = getattr(Configure.ConfigurationContext, "check_cfg")


@Configure.conf
def check_cfg(self, *k, **kw):
    ret = check_cfg_old(self, *k, **kw)
    if not self.env.CC_IS_INTEL and "clang" not in self.env.CC_NAME.lower():
        return ret

    to_be_removed = ["-ffat-lto-objects", "-flto", "-flto-partition=none", "-fuse-linker-plugin"]
    if kw["uselib_store"] in ("PYEMBED", "PYEXT"):
        use = kw["uselib_store"]
        keys = [i for i in self.env.keys() if use in i]
        for key in keys:
            inval = Utils.to_list(self.env[key])
            Logs.debug(f"{key}: {inval}")
            if type(inval) in (int, str):
                continue
            self.env[key] = [arg for arg in inval if arg not in to_be_removed]
            if self.env[key] != inval:
                Logs.debug("CHANGED: %s : %s", key, self.env[key])
    return ret
