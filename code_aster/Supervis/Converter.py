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
This module defines a converter object used to automatically adapt the
syntax of commands with trivial conversions.
"""

from typing import Any

from ..Objects import Function
from ..Utilities import logger


def type_fixer(value: Any, idx: int, valid_type: list[type]) -> tuple:
    """Try to fix a TypeError exception.

    Arguments:
        value (misc): value passed to the simple keyword, may be a list/tuple.
        idx (int): index of the value that raises the type error.
        valid_types (list[type]): list of expected types.

    Returns:
        tuple: the first element (type *bool*) tells if the value has changed,
        the second element is the value.
    """
    logger.debug("FIXER: value=%s, idx=%d, types=%s", value, idx, valid_type)
    returned = value
    changed = False
    as_list = isinstance(returned, (list, tuple))
    if as_list:
        value = returned[idx]
    typenames = []
    for typ in valid_type:
        try:
            typenames.append(typ.getType())
        except AttributeError:
            pass
    # convert float as constant function
    if isinstance(value, (float, int)) and "FONCTION_SDASTER" in typenames:
        changed = True
        value = constant_function(value)
        if as_list:
            returned[idx] = value
        else:
            returned = value
    return changed, returned


def constant_function(value: float) -> Function:
    """Create a constant function."""
    func = Function()
    func.setValues([1.0], [value])
    func.setResultName("TOUTRESU")
    func.setParameterName("TOUTPARA")
    func.setInterpolation("LIN LIN")
    func.setExtrapolation("CC")
    func.setAsConstant()
    return func
