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
:py:class:`MechanicalLoad` --- Assignment of mechanical load
*******************************************************************
"""

from libaster import MechanicalLoadComplex, MechanicalLoadFunction, MechanicalLoadReal

from ..Utilities import injector
from ..Objects.Serialization import InternalStateBuilder


class MechanicalLoadStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *MechanicalLoad*."""

    def save(self, obj):
        """Return the internal state of a *MechanicalLoad* to be pickled.

        Arguments:
            obj (*MechanicalLoad*): The *MechanicalLoad* object to be pickled.

        Returns:
            *InternalStateBuilder*: The internal state itself.
        """
        super().save(obj)
        self._st["fed"] = obj.getFiniteElementDescriptor()
        return self

    def restore(self, obj):
        """Restore the *MechanicalLoad* content from the previously saved internal
        state.

        Arguments:
            obj (*MechanicalLoad*): The *MechanicalLoad* object to be restored.
        """
        super().restore(obj)
        if self._st["fed"]:
            obj.setFiniteElementDescriptor(self._st["fed"])


@injector(MechanicalLoadReal)
class ExtendedMechanicalLoadReal:
    cata_sdj = "SD.sd_char_meca.sd_char_meca"
    internalStateBuilder = MechanicalLoadStateBuilder

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a MechanicalLoadReal
        object during unpickling.
        """
        return (self.getName(), self.getModel())


@injector(MechanicalLoadFunction)
class ExtendedMechanicalLoadFunction:
    cata_sdj = "SD.sd_char_meca.sd_char_meca"
    internalStateBuilder = MechanicalLoadStateBuilder

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a MechanicalLoadFunction
        object during unpickling.
        """
        return (self.getName(), self.getModel())


@injector(MechanicalLoadComplex)
class ExtendedMechanicalLoadComplex:
    cata_sdj = "SD.sd_char_meca.sd_char_meca"
    internalStateBuilder = MechanicalLoadStateBuilder

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a MechanicalLoadComplex
        object during unpickling.
        """
        return (self.getName(), self.getModel())
