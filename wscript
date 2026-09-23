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

"""
Build script for code_aster


Note:
- All defines conditionning the compilation must be set using `conf.define()`.
  They will be exported into `asterc_config.h`/`asterf_config.h`.

- If some of them are also required during the build step, another variable
  must be passed using `env` (ex. BUILD_MED)
"""

import os
import os.path as osp
import sys
import sysconfig

from waflib import Build, Configure, Logs, Utils
from waflib.Tools.c_config import DEFKEYS
from waflib.Tools.fc import fc

from waftools.wafutils import remove_previous

top = "."
out = "build"
install_suffix = os.environ.get("WAF_DEFAULT_VARIANT") or os.environ.get("WAF_SUFFIX", "mpi")
default_prefix = "../install/%s" % install_suffix

if sys.version_info < (3, 11):
    Logs.error("Python 3.11 or newer is required.")
    sys.exit(1)


def options(self):
    orig_get_usage = self.parser.get_usage

    def _usage():
        return orig_get_usage() + os.linesep.join((
            "",
            "Environment variables:",
            "  CC             : C compiler",
            "  FC             : Fortran compiler",
            "  CXX            : C++ compiler",
            "  DEFINES        : extra preprocessor defines",
            "  LINKFLAGS      : extra linker options",
            "  CFLAGS         : extra C compilation options",
            "  CXXFLAGS       : extra C++ compilation options",
            "  FCFLAGS        : extra Fortran compilation options",
            "  {C,CXX,FC}FLAGS_ASTER_DEBUG may be used to add options only in debug mode",
            '  LIBPATH_x, LIB_x, INCLUDES_x, PYPATH_x : paths for component "x" for libs, '
            "includes, python modules",
            "  CONFIG_PARAMETERS_name=value: extra configuration parameters (for config.yaml/json)",
            "  WAFBUILD_ENV   : environment file to be included in runtime environment file",
            "  PREFIX         : default installation prefix to be used, "
            "if no --prefix option is given.",
            "  ASTER_BLAS_INT_SIZE  : kind of integers to use in the fortran blas/lapack "
            "calls (4 or 8, default is 4)",
            "  ASTER_MUMPS_INT_SIZE : kind of integers to use in the fortran mumps calls "
            " (4 or 8, default is 4)",
            "  CATALO_CMD     : command line used to build the elements catalog. "
            "It is just inserted before the executable "
            "(may define additional environment variables or a wrapper that takes "
            "all arguments, see catalo/wscript)",
            "",
        ))

    self.parser.get_usage = _usage

    self.load("use_config")
    self.load("gnu_dirs")

    # see waflib/Tools/gnu_dirs.py for the group name
    group = self.get_option_group("Installation prefix")
    descr = group.description or ""
    # replace path in description
    new_descr = descr.replace("/usr/local", default_prefix)
    new_descr += (
        ". Using 'waf_variant', '%s' will be automatically replaced by 'variant'." % install_suffix
    )
    group.description = new_descr
    # change default value for '--prefix'
    option = self.parser.get_option("--prefix")
    if option:
        self.parser.remove_option("--prefix")
    group.add_option(
        "--prefix",
        dest="prefix",
        default=None,
        help="installation prefix [default: %r]" % default_prefix,
    )

    group = self.get_option_group("Build and installation options")
    group.add_option(
        "--fast",
        dest="custom_fc_sig",
        action="store_true",
        default=False,
        help="use fast algorithm based on modification time to "
        "check for dependencies of fortran sources",
    )
    group.add_option(
        "--safe",
        dest="custom_fc_sig",
        action="store_false",
        help="use safe algorithm based on content to check for "
        "implicit dependencies of fortran sources",
    )
    group.add_option(
        "--coverage",
        dest="enable_coverage",
        action="store_true",
        help="enable options for coverage (actually only for debug build)",
    )
    group.add_option(
        "--enable-asan",
        dest="enable_asan",
        action="store_true",
        default=False,
        help="enable address sanitizer in debug mode",
    )
    group.add_option(
        "--disable-asan",
        dest="enable_asan",
        action="store_false",
        help="disable address sanitizer in debug mode",
    )
    group.add_option(
        "--mingw-cross-compilation",
        dest="dest_mingw",
        action="store_true",
        default=os.environ.get("ENABLE_MINGW", "0") != "0",
        help="enable cross compilation for mingw",
    )
    # MSVC Options
    group.add_option(
        "--site-packages",
        dest="spdir",
        default="auto",
        help="Python site-packages directory (absolute path or relative to PREFIX), "
        "special values: 'auto' default Python site-packages, 'no' direcly in lib.",
    )
    group.add_option(
        "--bindir",
        dest="bindir",
        default=None,
        help="Binary directory for DLL files on Windows conda/rattler builds [default: PREFIX/bin]",
    )
    group.add_option(
        "--enable-aster-subdir",
        dest="use_aster_subdir",
        action="store_true",
        help="Enable the '/aster' subdirectory in installation paths "
        "(may be useful for some packaging)",
    )
    group.add_option(
        "--disable-aster-subdir",
        dest="use_aster_subdir",
        action="store_false",
        default=False,
        help="Disable the '/aster' subdirectory in installation paths "
        "(this is the default, useful for conda/rattler builds where files _"
        "should go directly to lib)",
    )
    group = self.add_option_group("code_aster options")

    self.load("parallel", tooldir="waftools")
    self.load("python_cfg", tooldir="waftools")
    self.load("mathematics", tooldir="waftools")
    self.load("med_cfg", tooldir="waftools")
    self.load("metis", tooldir="waftools")
    self.load("parmetis", tooldir="waftools")
    self.load("mumps", tooldir="waftools")
    self.load("scotch", tooldir="waftools")
    self.load("petsc", tooldir="waftools")
    self.load("runtest", tooldir="waftools")
    self.recurse("libs")
    self.recurse("bibfor")
    self.recurse("code_aster")
    self.recurse("run_aster")
    self.recurse("bibcxx")
    self.recurse("bibc")
    self.recurse("mfront")
    self.recurse("i18n")
    self.recurse("data")
    self.recurse("doc")
    self.recurse("astest")

    group.add_option(
        "-E",
        "--embed-all",
        dest="embed_all",
        action="store_true",
        default=False,
        help="activate all embed-* options (except embed-python)",
    )
    group.add_option(
        "--enable-all",
        dest="enable_all",
        action="store_true",
        default=None,
        help="activate all 'enable-*' options, means that all prerequisites are required"
        "(same as ENABLE_ALL=1 environment variable)",
    )
    group.add_option(
        "--no-enable-all",
        dest="enable_all",
        action="store_false",
        help="try to build with some missing prerequisites (same as "
        "ENABLE_ALL=0 environment variable)",
    )


@Configure.conf
def all_components(self):
    components = [
        "HDF5",
        "MED",
        "MFRONT",
        "MGIS",
        "METIS",
        "PARMETIS",
        "SCOTCH",
        "MUMPS",
        "PETSC",
        "PETSC4PY",
        "MATH",
        "NUMPY",
        "PYBIND11",
        "OPENMP",
        "MPI",
        "Z",
        "CXX",
        "ASRUN",
        "PYTHON",
    ]
    return components


def configure(self):
    opts = self.options
    self.setenv("default")
    self.load("official_platforms", tooldir="waftools")

    # add environment variables into `self.env`
    self.add_os_flags("FC")
    self.add_os_flags("CC")
    self.add_os_flags("CXX")
    self.add_os_flags("CFLAGS")
    self.add_os_flags("CXXFLAGS")
    self.add_os_flags("FCFLAGS")
    self.add_os_flags("LINKFLAGS")
    self.add_os_flags("DEFINES")
    self.add_os_flags("WAFBUILD_ENV")
    self.add_os_flags("CFLAGS_ASTER_DEBUG")
    self.add_os_flags("CXXFLAGS_ASTER_DEBUG")
    self.add_os_flags("FCFLAGS_ASTER_DEBUG")

    all_prods = self.all_components()
    if opts.enable_all is None:
        enabled_comp = [getattr(opts, "enable_" + comp.lower(), None) for comp in all_prods]
        # if a product is explicitly disabled, can not enable all!
        enable_all = os.environ.get("ENABLE_ALL") != "0" and False not in enabled_comp
        opts.enable_all = enable_all
    for comp in all_prods:
        self.add_os_flags("LIBPATH_" + comp)
        self.add_os_flags("LIB_" + comp)
        self.add_os_flags("INCLUDES_" + comp)
        self.add_os_flags("PYPATH_" + comp)
        if not opts.enable_all:
            continue
        if not opts.parallel and comp in ["PARMETIS", "PETSC", "MPI"]:
            continue
        opt = "enable_" + comp.lower()
        if hasattr(opts, opt):
            setattr(opts, opt, True)

    self.env["CONFIG_PARAMETERS"] = {}
    for key in os.environ.keys():
        if not key.startswith("CONFIG_PARAMETERS_"):
            continue
        name = key.split("CONFIG_PARAMETERS_")[1]
        self.env["CONFIG_PARAMETERS"][name] = os.environ[key]

    # compute default prefix
    if self.env.PREFIX in ("", "/"):
        self.env.PREFIX = osp.abspath(default_prefix)
        if "PREFIX_ROOT" in os.environ:
            self.env.PREFIX = osp.join(os.environ["PREFIX_ROOT"], install_suffix)
    self.msg("Setting prefix to", self.env.PREFIX)

    # wine wrapper must be set up before checking compilers
    self.env.DEST_MINGW = opts.dest_mingw
    self.load("ext_aster", tooldir="waftools")
    self.load("use_config")
    self.load("gnu_dirs")
    self.env["CODEASTERPATH"] = self.path.find_dir("code_aster").abspath()
    self.env["CATALOPATH"] = self.path.find_dir("catalo").abspath()

    self.load("parallel", tooldir="waftools")
    self.check_platform()
    self.set_installdirs()
    self.load("python_cfg", tooldir="waftools")
    # keep compatibility for as_run
    if self.get_define("ASTER_HAVE_MPI"):
        self.env.ASRUN_MPI_VERSION = 1

    self.recurse("libs")
    # bib* configure functions may add options required by prerequisites
    self.recurse("bibfor")
    self.recurse("bibcxx")
    self.recurse("bibc")
    self.check_asan()

    self.load("mathematics", tooldir="waftools")
    self.load("med_cfg", tooldir="waftools")
    self.load("metis", tooldir="waftools")
    self.load("parmetis", tooldir="waftools")
    self.load("mumps", tooldir="waftools")
    self.load("scotch", tooldir="waftools")
    self.load("petsc", tooldir="waftools")
    self.load("runtest", tooldir="waftools")

    self.recurse("code_aster")
    self.recurse("run_aster")
    self.recurse("extern")
    self.recurse("mfront")
    self.recurse("i18n")
    self.recurse("data")
    self.recurse("doc")
    self.recurse("astest")

    # variants
    self.check_optimization_options()
    self.write_config_headers()


def build(self):
    env = self.env.derive()

    fc._use_custom_sig = self.options.custom_fc_sig
    # shared the list of dependencies between bibc/bibfor
    # the order may be important
    if not self.variant:
        self.fatal(
            'Call "waf build_debug" or "waf build_release", and read '
            "the comments in the wscript file!"
        )
    if self.variant == "release" and env["CFLAGS_ASAN"]:
        self.fatal(
            "The project was configured with '--enable-asan' option (AddressSanitizer). "
            "Only the 'debug' variant is relevant."
        )
    if self.cmd.startswith("install"):
        # because we can't know which files are obsolete `rm *.py{,c,o}`
        # Only clean code_aster and run_aster subdirectories
        code_aster_dir = self.root.find_node(osp.join(env.SITEPACKAGESDIR, "code_aster"))
        run_aster_dir = self.root.find_node(osp.join(env.SITEPACKAGESDIR, "run_aster"))
        if code_aster_dir:
            remove_previous(code_aster_dir, ["**/*.py", "**/*.pyc", "**/*.pyo"])
        if run_aster_dir:
            remove_previous(run_aster_dir, ["**/*.py", "**/*.pyc", "**/*.pyo"])
        remove_previous(
            self.root.find_node(env.ASTERDATADIR), ["datg/**/*", "materiau/**/*", "tests_data/**/*"]
        )

    self.load("ext_aster", tooldir="waftools")
    self.recurse("libs")
    self.recurse("bibfor")
    self.recurse("code_aster")
    self.recurse("run_aster")
    self.recurse("bibcxx")
    self.recurse("bibc")
    self.recurse("mfront")
    self.recurse("i18n")
    self.recurse("catalo")
    self.recurse("data")
    self.recurse("astest")
    self.install_as(
        osp.join(env.SITEPACKAGESDIR, "code_aster", "Utilities", "aster_config.py"),
        ["aster_config.py"],
    )


def build_elements(self):
    self.recurse("catalo")


def init(self):
    from waflib.Build import BuildContext, CleanContext, InstallContext, UninstallContext

    _all = (
        BuildContext,
        CleanContext,
        InstallContext,
        UninstallContext,
        TestContext,
        I18NContext,
        DocContext,
    )
    for x in ["debug", "release"]:
        for y in _all:
            name = y.__name__.replace("Context", "").lower()

            class tmp(y):
                cmd = name + "_" + x
                variant = x

    # default to release
    for y in _all:

        class tmp(y):
            variant = os.environ.get("WAF_DEFAULT_VARIANT") or "release"


def all(self):
    from waflib import Options

    lst = ["install_release", "install_debug"]
    Options.commands = lst + Options.commands


class BuildElementContext(Build.BuildContext):
    """execute the build for elements catalog only using an installed Aster (also performed at install, for internal use only)"""

    cmd = "_buildelem"
    fun = "build_elements"


def runtest(self):
    self.load("runtest", tooldir="waftools")


class TestContext(Build.BuildContext):
    """facility to execute a testcase"""

    cmd = "test"
    fun = "runtest"


def update_i18n(self):
    self.recurse("i18n")


class I18NContext(Build.BuildContext):
    """build the i18n files"""

    cmd = "i18n"
    fun = "update_i18n"


def build_doc(self):
    self.recurse("doc")


class DocContext(Build.BuildContext):
    """build the documentation files"""

    cmd = "doc"
    fun = "build_doc"


@Configure.conf
def reset_msg(self):
    """Reset message level"""
    if self.in_msg:
        self.end_msg("cancelled", "YELLOW")
    self.in_msg = 0


@Configure.conf
def set_installdirs(self):
    # set the installation subdirectories
    # Check if --libdir was explicitly provided and override LIBDIR if needed
    # The gnu_dirs module defaults to lib64 on 64-bit systems, but conda/rattler
    # builds typically expect lib to be used
    if self.options.libdir:
        # User explicitly set libdir, use it
        self.env.LIBDIR = self.options.libdir

    # Check if --bindir was explicitly provided for Windows DLL installation
    if self.options.bindir:
        self.env.BINDIR = self.options.bindir
    else:
        # Default to PREFIX/bin if not specified
        self.env.BINDIR = osp.join(self.env.PREFIX, "bin")

    # For conda/rattler builds, we don't want the /aster subdirectory for libraries
    # but we DO want it for data files (tests, profile.sh, etc.)
    if self.options.use_aster_subdir:
        norm_lib = lambda path: osp.normpath(osp.join(path, "aster"))
    else:
        norm_lib = lambda path: osp.normpath(path)
    norm_data = lambda path: osp.normpath(osp.join(path, "aster"))

    self.env["ASTERLIBDIR"] = norm_lib(self.env.LIBDIR)
    self.env["ASTERBINDIR"] = self.env.BINDIR
    self.env["ASTERINCLUDEDIR"] = norm_lib(self.env.INCLUDEDIR)
    self.env["ASTERDATADIR"] = norm_data(self.env.DATADIR)
    if not self.env.LOCALEDIR:
        self.env.LOCALEDIR = osp.join(self.env.PREFIX, "share", "locale")
    self.env["ASTERLOCALEDIR"] = norm_data(self.env.LOCALEDIR)
    # set relative paths for profile.sh
    for var in ("LIBDIR", "DATADIR", "LOCALEDIR"):
        self.env["RELATIVE_" + var] = osp.relpath(self.env["ASTER" + var], self.env["PREFIX"])
    # use --site-packages by default
    if self.options.spdir == "auto":
        scheme = "posix_user" if self.env.ASTER_PLATFORM_POSIX else "nt_user"
        self.env.SITEPACKAGESDIR = sysconfig.get_path(
            "platlib", scheme=scheme, vars={"userbase": self.env.PREFIX}
        )
    elif self.options.spdir == "no":
        self.env.SITEPACKAGESDIR = self.env.ASTERLIBDIR
    else:
        self.env.SITEPACKAGESDIR = osp.join(self.env.PREFIX, self.options.spdir)
    self.env["RELATIVE_SITEPACKAGESDIR"] = osp.relpath(self.env.SITEPACKAGESDIR, self.env["PREFIX"])
    # for aster_config.py
    self.define("ASTER_BINDIR", self.env["ASTERBINDIR"])
    self.define("ASTER_LIBDIR", self.env["ASTERLIBDIR"])
    self.define("ASTER_DATADIR", self.env["ASTERDATADIR"])
    self.define("ASTER_SITEPACKAGESDIR", self.env.SITEPACKAGESDIR)
    self.msg("Setting bindir to", self.env["ASTERBINDIR"])
    self.msg("Setting libdir to", self.env["ASTERLIBDIR"])
    self.msg("Setting resource dir to", self.env["ASTERDATADIR"])
    self.msg("Setting python dir to", self.env["SITEPACKAGESDIR"])


@Configure.conf
def check_platform(self):
    self.start_msg("Getting platform")
    # convert waf (sys.platform) to code_aster terminology
    os_name = self.env.DEST_OS
    if os_name == "cygwin":
        os_name = "linux"
    elif os_name == "win32":
        os_name = "mingw32"
    if "64" in self.env.DEST_CPU:
        if os_name.endswith("32"):
            os_name = os_name[:-2]
        os_name += "64"
        self.define("ASTER_HAVE_64_BITS", 1)
    plt = "ASTER_PLATFORM_" + os_name.upper()
    if not os_name.startswith("mingw"):
        self.define("ASTER_PLATFORM_POSIX", 1)
        self.env.ASTER_PLATFORM_POSIX = True
        self.undefine("ASTER_PLATFORM_MINGW")
    else:
        self.define("ASTER_PLATFORM_MINGW", 1)
        self.env.ASTER_PLATFORM_MINGW = True
        self.undefine("ASTER_PLATFORM_POSIX")
    self.env.ASTER_PLATFORM = plt
    if os.getenv("MSYSTEM"):
        ## Define and additional variable for MSYS2
        self.define("ASTER_PLATFORM_MSYS2", 1)
        self.env.ASTER_PLATFORM_MSYS2 = True

    self.define(plt, 1)
    self.end_msg(plt)


@Configure.conf
def check_optimization_options(self):
    # adapt the environment of the build variants
    self.setenv("debug", env=self.all_envs["default"])
    self.setenv("release", env=self.all_envs["default"])
    # these functions must switch between each environment
    self.check_optimization_cflags()
    self.check_optimization_cxxflags()
    self.check_optimization_fcflags()
    self.check_optimization_python()
    self.check_variant_vars()


@Configure.conf
def check_variant_vars(self):
    self.setenv("debug")
    self.env["ASTER_BEHAVIOUR_LIB"] = "AsterMFrOfficialDebug"
    self.define("ASTER_BEHAVIOUR_LIB", self.env["ASTER_BEHAVIOUR_LIB"])

    self.setenv("release")
    self.env["ASTER_BEHAVIOUR_LIB"] = "AsterMFrOfficial"
    self.define("ASTER_BEHAVIOUR_LIB", self.env["ASTER_BEHAVIOUR_LIB"])


# same idea than waflib.Tools.c_config.write_config_header
# but defines are not removed from `env`
# XXX see write_config_header(remove=True/False) + format Fortran ?
class ConfigHelper:
    def __init__(self, language):
        self._lang = language

    def cmt(self, text):
        return {"C": "/* {0} */", "Fortran": "! {0}", "Python": "# {0}"}[self._lang].format(text)

    @property
    def filename(self):
        return {"C": "asterc_config.h", "Fortran": "asterf_config.h", "Python": "aster_config.py"}[
            self._lang
        ]

    @property
    def header(self):
        return [self.cmt("WARNING! Automatically generated by `waf configure`!")]

    @property
    def guard_begin(self):
        if self._lang in ("C", "Fortran"):
            guard = Utils.quote_define_name(self.filename)
            return ["#ifndef {0}_".format(guard), self.define(guard).strip() + "_", ""]
        elif self._lang == "Python":
            return ["config = dict("]
        return [""]

    @property
    def guard_end(self):
        if self._lang in ("C", "Fortran"):
            return ["#endif", ""]
        elif self._lang == "Python":
            return [")"]
        return [""]

    def support(self, var):
        """Tell if the language supports the variable name."""
        if var.startswith("HAVE_"):
            # do not keep variables automatically set by waf/check_cfg
            return False
        # ASTER_C_* defines will be used if language='C', not 'Fortran'.
        if self._lang != "C" and var.startswith("ASTER_C_"):
            return False
        return True

    def define(self, var, value=""):
        if self._lang == "Python":
            if not value.isdigit():
                value = "r" + value
            fmt = "   {0}=" + "{1}," if value else "True"
        else:
            fmt = "#define {0} {1}"
        return fmt.format(var, value)

    def undefine(self, var):
        fmt = self.cmt("#undef {0}")
        return fmt.format(var)


@Configure.conf
def write_config_headers(self):
    # Write both xxxx_config.h files for C and Fortran,
    # then remove entries from DEFINES
    for variant in ("debug", "release"):
        self.setenv(variant)
        self.write_config_h("Fortran", variant)
        self.write_config_h("C", variant)
        self.write_config_h("Python", variant)
        for key in self.env[DEFKEYS]:
            self.undefine(key)
        self.env[DEFKEYS] = []


@Configure.conf
def write_config_h(self, language, variant, configfile=None, env=None):
    # Write a configuration header containing defines
    self.start_msg("Write config file")
    cfg = ConfigHelper(language)
    configfile = configfile or cfg.filename
    env = env or self.env
    lst = cfg.header
    lst.append("")
    lst.extend(cfg.guard_begin)
    lst.extend(self.get_config_h(cfg))
    lst.extend(cfg.guard_end)

    node = self.bldnode or self.path.get_bld()
    node = node.make_node(osp.join(variant, configfile))
    node.parent.mkdir()
    node.write("\n".join(lst))
    incpath = node.parent.abspath()
    env.append_unique("INCLUDES", incpath)
    # config files are not removed on "waf clean"
    env.append_unique(Build.CFG_FILES, [node.abspath()])
    self.end_msg(node.bldpath())


@Configure.conf
def get_config_h(self, cfg):
    # Create the contents of a ``config.h`` file from the defines
    # set in conf.env.define_key / conf.env.include_key. No include guards are added.
    lst = []
    for x in self.env[DEFKEYS]:
        if not cfg.support(x):
            continue
        if self.is_defined(x):
            lst.append(cfg.define(x, self.get_define(x)))
        else:
            lst.append(cfg.undefine(x))
    lst.append("")
    return lst
