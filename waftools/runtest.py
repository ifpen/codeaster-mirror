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

import json
import os
import os.path as osp
import platform
import tempfile
from configparser import ConfigParser
from fnmatch import fnmatchcase
from glob import glob
from pathlib import Path
from subprocess import Popen

from waflib import Errors, Logs, TaskGen

try:
    import yaml
except ImportError:
    yaml = None


def options(self):
    """To get the names of the testcases"""
    group = self.get_option_group("code_aster options")
    group.add_option(
        "-n",
        "--name",
        dest="testname",
        action="append",
        default=None,
        help="name of testcases to run (as_run must be in PATH)",
    )
    group.add_option(
        "--outputdir",
        action="store",
        default=None,
        metavar="DIR",
        help="directory to store the output files. A default "
        "value can be stored in `~/.config/aster/config.yaml`",
    )
    group.add_option(
        "--exectool",
        dest="exectool",
        action="store",
        default=None,
        help="run a testcase by passing additional arguments "
        '(possible values are "env" + those '
        "defined in the configuration)",
    )
    group.add_option(
        "--time_limit",
        dest="time_limit",
        action="store",
        default=None,
        help="override the time limit of the testcase",
    )


def configure(self):
    """Store developer preferences"""
    prefs = {}
    key = "outputdir"

    fcfg = Path.home() / ".config" / "aster" / "config.yaml"
    self.start_msg("Reading user prefs from %s" % fcfg)
    if not fcfg.is_file() or not yaml:
        self.end_msg("not found")
        fcfg = fcfg.with_name("config.json")
        self.start_msg("Reading user prefs from %s" % fcfg)
    value = self.options.outputdir or ""
    if value:
        pass
    elif fcfg.is_file():
        with fcfg.open("rb") as fobj:
            if fcfg.suffix == ".yaml":
                content = yaml.load(fobj, Loader=yaml.Loader)
            else:
                content = json.load(fobj)
        params = _filter(content, "server", "name", platform.node())
        params.update(_filter(content, "version", "path", self.env["PREFIX"]))
        value = params.get(key, "")
    else:
        self.end_msg("not found")
        fcfg = Path.home() / ".gitconfig"
        self.start_msg("Reading user prefs from %s" % fcfg)
        if fcfg.is_file():
            # strict=False will ignore duplicated sections (may at least occur with fetch)
            cfg = ConfigParser(strict=False)
            cfg.read(str(fcfg))
            value = cfg.get("aster", key, fallback="")
        else:
            self.end_msg("not found")
            self.start_msg("No resource file found, empty content")

    dkey = "PREFS_{}".format(key.upper())
    self.env[dkey] = value
    if value:
        prefs[key] = value
    self.end_msg(repr(prefs))


@TaskGen.feature("test")
def runtest(self):
    """Run a testcase by calling as_run"""
    opts = self.options
    run_aster = osp.join(self.env["BINDIR"], "run_aster")
    if not osp.isfile(run_aster):
        Logs.error("'run_aster' not found, please check your $PATH")
        return
    args = []
    if opts.exectool == "env":
        args.append("--env")
        workdir = tempfile.mkdtemp(prefix="runtest_")
        args.extend(["--workdir", workdir])
    elif opts.exectool is not None:
        args.append("--exectool=%s" % opts.exectool)
    if opts.time_limit:
        args.append("--time_limit={0}".format(opts.time_limit))
    dtmp = opts.outputdir or self.env["PREFS_OUTPUTDIR"] or tempfile.mkdtemp(prefix="runtest_")
    try:
        os.makedirs(dtmp)
    except (OSError, IOError):
        pass
    Logs.info("destination of output files: %s" % dtmp)
    status = 0
    if not opts.testname:
        raise Errors.WafError("no testcase name provided, use the -n option")
    for test in opts.testname:
        export = test + ".export"
        exp = glob("astest/" + export) + glob("../validation/astest/" + export)
        if not exp:
            raise FileNotFoundError(test + ".export")
        cmd = [run_aster, "--test"]
        cmd.extend(args)
        cmd.append(osp.abspath(exp[0]))
        Logs.info("running %s in '%s'" % (test, self.variant))
        ext = "." + osp.basename(self.env["PREFIX"]) + "." + self.variant
        out = osp.join(dtmp, osp.basename(test) + ext) + ".output"
        err = osp.join(dtmp, osp.basename(test) + ext) + ".error"
        Logs.info("`- command: %s" % (" ".join(cmd)))
        Logs.info("`- output in %s" % out)
        # do not run from source directory to import installed files
        current = os.getcwd()
        os.chdir(dtmp)
        with open(out, "w") as fobj, open(err, "w") as ferr:
            proc = Popen(cmd, stdout=fobj, stderr=ferr, bufsize=1)
        retcode = proc.wait()
        os.chdir(current)
        with open(out, "rb") as fobj:
            btext = fobj.read()
        text = btext.decode("utf8", "replace")
        if retcode == 2 and "NOOK_TEST_RESU" in text:
            retcode = "nook"
        if retcode == 0:
            func = Logs.info
        else:
            func = Logs.error
            status += 1
        func("`- exit %s" % retcode)
    if status != 0:
        raise Errors.WafError("testcase failed")


# function extracted from run_aster/config.py
def _filter(content, section, filter_key, filter_value):
    """Filter content by keeping sections that match the filter.

    Arguments:
        content (dict): file content with optional "server" and
            "version" list.

    Returns:
        dict: Version parameters for the current server and version.
    """
    params = {}
    candidates = content.get(section, [])
    if not isinstance(candidates, list):
        candidates = [candidates]
    for cfg in candidates:
        if not isinstance(cfg, dict):
            print("warning: dict expected for %r, not: %s", section, cfg)
            continue
        if not fnmatchcase(filter_value, cfg.get(filter_key, "")):
            continue
        config = cfg.get("config", {})
        if not isinstance(config, dict):
            print("warning: dict expected for 'config', not: %s", config)
            continue
        params.update(config)
    return params
