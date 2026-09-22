#!/usr/bin/env python3
# coding: utf-8

"""
Generate automodule blocks.
"""

import argparse
import os
import os.path as osp
from collections import OrderedDict

EPILOG = """
EXAMPLES:
    Generate a file like ``supervis.rst``::

        python3 generate_rst.py --manual code_aster/__init__.py code_aster/Supervis/*.py

    Generate files for *DataStructure* and derivated subclasses::

        python3 generate_rst.py --objects
"""

SECTIONS = ("DataStructure",)

automodule_block = """.. automodule:: {0}
   :show-inheritance:
   :members:
   :special-members: __init__
""".format

autoclass_block = """.. autoclass:: code_aster.Objects.{0}
   :show-inheritance:
   :members:
""".format

auto_documentation = """.. AUTOMATICALLY CREATED BY generate_rst.py - DO NOT EDIT MANUALLY!

{link}

{intro}

{content}
""".format

title_ds = """
################################################################################
Subclasses of :py:class:`~code_aster.Objects.{0}`
################################################################################
""".format

title_ds_alone = """
********************************************************************************
:py:class:`~code_aster.Objects.{0}` object
********************************************************************************
""".format

subtitle = """
********************************************************************************
:py:class:`~code_aster.Objects.{0}` object
********************************************************************************
""".format

table_decl_type = """
.. list-table::
   :widths: 40 25
   :header-rows: 1

   * - C++/Python object
     - type in the commands catalog [#f1]_"""

table_line_type = """   * - :py:class:`~code_aster.Objects.{0}`
     - {1}""".format

table_footer_type = """
.. rubric:: Footnotes

.. [#f1] If empty, it means that ``.getType()`` can not be called without a full execution,
         see the related header file.
"""

table_decl = """
.. list-table::
   :widths: 40
   :header-rows: 1

   * - Class name"""

table_line = """   * - :py:class:`~code_aster.Objects.{0}`""".format


def automodule(filename):
    """Return a block to document a module."""
    name = osp.splitext(filename)[0].replace("/", ".")
    return automodule_block(name)


def _build_text():
    """Generate sphinx blocks for all libaster objects."""
    import code_aster.Objects as OBJ

    pyb_instance = OBJ.DataStructure.mro()[1]
    # pyb_enum = OBJ.Physics.mro()[1]
    pyb_enum = object()

    # sections: directly derivated from pybind11 instance
    sections = [OBJ.DataStructure]
    addsect = []
    for name, obj in list(OBJ.__dict__.items()):
        if not isinstance(obj, type):
            continue
        # if pyb_instance in obj.mro():
        # some objects are missed? or included via dependencies?
        # check AcousticDirichletBC for example
        if obj.mro()[1] is pyb_instance:
            # print("1:", name, obj.mro())
            addsect.append((name, obj))
    for _, obj in sorted(addsect):
        sections.append(obj)
    # sections.append(pyb_enum)
    sections.append(Exception)
    # print(len(sections), "sections")

    # dict of names of subclasses
    dictobj = OrderedDict([(i, []) for i in sections])
    for name, obj in list(OBJ.__dict__.items()):
        # if obj is not OBJ.Material:
        #     continue
        if not isinstance(obj, type) or issubclass(
            obj, (OBJ.UnavailableObject, OBJ.InternalStateBuilder, OBJ.NamedTuple)
        ):
            continue
        found = False
        for subtyp in sections:
            if issubclass(obj, subtyp) or obj is subtyp:
                dictobj[subtyp].append((name, obj))
                found = True
                # print("Found:", name ,">>>", subtyp)
                break
        if not found and not issubclass(obj, OBJ.PyDataStructure):
            raise KeyError("pybind11 class not found: {0}".format(obj.mro()))

    dictdesc = OrderedDict()
    dicttabl = OrderedDict()
    for subtyp, couples in list(dictobj.items()):
        typename = subtyp.__name__
        # put subclass first
        try:
            couples.remove((typename, subtyp))
        except ValueError:
            if subtyp not in (pyb_enum, Exception):
                print(subtyp, typename)
                print([i[0] for i in couples])
                raise
        couples.sort()
        if subtyp not in (pyb_enum, Exception):
            couples.insert(0, (typename, subtyp))
        objs = [i[0] for i in couples]

        lines = []
        ltab = []
        if len(objs) > 1 and typename in SECTIONS:
            lines.append(title_ds(typename))
        else:
            lines.append(title_ds_alone(typename))
        nberr = 0
        for name, astyp in couples:
            if typename in SECTIONS:
                lines.append(subtitle(name))
                try:
                    inst = astyp()
                    typ_extr = "``" + inst.getType().lower() + "``"
                except TypeError:
                    typ_extr = ""
                    nberr += 1
                ltab.append(table_line_type(name, typ_extr))
            else:
                ltab.append(table_line(name))
            lines.append(autoclass_block(name))
        print(f"{nberr} errors / {len(couples)} types")

        dictdesc[typename] = os.linesep.join(lines)
        if ltab:
            if typename in SECTIONS:
                ltab.insert(0, table_decl_type)
                ltab.append(table_footer_type)
            dicttabl[typename] = os.linesep.join(ltab)
    # raise ValueError("DEBUG: dictdesc keys:", list(dictdesc.keys()))
    return dictdesc, dicttabl


def all_objects(destdir):
    """Generate sphinx blocks for all libaster objects."""
    descr, tables = _build_text()

    # generate a page for each of the section
    for sect in SECTIONS:
        cnt = descr.pop(sect)
        with open(osp.join(destdir, f"objects_{sect}.rst"), "w") as fobj:
            params = dict(link=f".. _devguide-objects_{sect}:", content=cnt, intro="")
            fobj.write(auto_documentation(**params))
        cnt = tables.pop(sect)
        with open(osp.join(destdir, f"table_{sect}.rst"), "w") as fobj:
            params = dict(link="", content=cnt, intro="")
            fobj.write(auto_documentation(**params))

    # generate a page for all other classes
    with open(osp.join(destdir, "objects_Others.rst"), "w") as fobj:
        params = dict(
            link=".. _devguide-objects_Others:",
            content=os.linesep.join(list(descr.values())),
            intro="""
####################################
Index of all other available objects
####################################

Documentation of all other types.
""",
        )
        fobj.write(auto_documentation(**params))

    # raise ValueError(list(tables.values())[1])
    ltabs = list(tables.values())
    ltabs.insert(0, table_decl)
    with open(osp.join(destdir, "table_Others.rst"), "w") as fobj:
        params = dict(link="", content=os.linesep.join(ltabs), intro="")
        fobj.write(auto_documentation(**params))


def main():
    default_dest = os.getcwd()
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        action="store_const",
        dest="action",
        const="objects",
        help="for C++ only objects (needs to import libaster)",
    )
    parser.add_argument(
        "--manual", action="store_const", dest="action", const="manual", help="for only few objects"
    )
    parser.add_argument(
        "-d",
        "--destdir",
        action="store",
        metavar="DIR",
        default=default_dest,
        help="directory where `rst` files will be written",
    )
    parser.add_argument("file", metavar="FILE", nargs="*", help="file to analyse")
    args = parser.parse_args()

    if args.action == "objects":
        all_objects(args.destdir)
    elif args.action == "manual":
        for name in args.file:
            print(automodule(name))


if __name__ == "__main__":
    import code_aster

    code_aster.rc.restart = False
    from code_aster import CA

    all_objects(os.environ["AUTODOC_DESTDIR"])
    CA.close()
