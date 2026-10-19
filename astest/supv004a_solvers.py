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

import unittest
from unittest.mock import MagicMock

from code_aster.Commands import DEFI_LIST_ENTI, DEFI_LIST_INST, DEFI_LIST_REEL
from code_aster import CA
from code_aster.Solvers import (
    ConvergenceManager,
    PhysicalState,
    ProblemType,
    StorageManager,
    TimeStepper,
)


class TestTimeStepper(unittest.TestCase):
    """Check for TimeStepper."""

    def setUp(self):
        self.list0 = DEFI_LIST_REEL(VALE=0.0)
        self.listr = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=10.0, PAS=1.0))
        self.listi = DEFI_LIST_INST(DEFI_LIST=_F(VALE=(1.0, 2.0)))

    def test00_init(self):
        stp = TimeStepper([0.0, 1.0, 2.0, 3.0])
        self.assertSequenceEqual(stp._times, [1.0, 2.0, 3.0])
        self.assertEqual(stp.size(), 3)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0], initial=1.0)
        self.assertSequenceEqual(stp._times, [2.0, 3.0])
        self.assertEqual(stp.size(), 2)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0])
        self.assertSequenceEqual(stp._times, [1.0, 2.0, 3.0])
        self.assertEqual(stp.size(), 3)
        self.assertEqual(stp.remaining(), stp.size())

        stp.setInitial(1.0)
        self.assertSequenceEqual(stp._times, [2.0, 3.0])
        self.assertEqual(stp.size(), 2)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0], initial=0.0, final=2.5)
        self.assertSequenceEqual(stp._times, [1.0, 2.0, 2.5])
        self.assertEqual(stp.size(), 3)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0])
        stp.setInitial(0.0)
        stp.setFinal(2.5)
        self.assertSequenceEqual(stp._times, [1.0, 2.0, 2.5])
        self.assertEqual(stp.size(), 3)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0])
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        self.assertAlmostEqual(stp.getInitial(), stp.getFinal())
        self.assertSequenceEqual(stp._times, [])
        self.assertEqual(stp.size(), 0)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0], initial=3.0)
        self.assertAlmostEqual(stp.getInitial(), 3.0)
        self.assertAlmostEqual(stp.getInitial(), stp.getFinal())
        self.assertSequenceEqual(stp._times, [])
        self.assertEqual(stp.size(), 0)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([0.0, 1.0, 2.0, 3.0], initial=2.0, final=2.5)
        self.assertSequenceEqual(stp._times, [2.5])
        self.assertEqual(stp.size(), 1)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([1.0, 2.0, 3.0], final=0.5)
        self.assertSequenceEqual(stp._times, [0.5])
        self.assertEqual(stp.size(), 1)
        self.assertEqual(stp.remaining(), stp.size())

        stp = TimeStepper([], final=1.0)
        self.assertSequenceEqual(stp._times, [1.0])
        self.assertEqual(stp.size(), 1)
        self.assertEqual(stp.remaining(), stp.size())

        with self.assertRaisesRegex(ValueError, "ordered"):
            TimeStepper([0.0, 1.0, 3.0, 2.0])

    def test01_initial(self):
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitial(0.0)
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        self.assertEqual(stp.size(), 3)
        self.assertEqual(stp.remaining(), stp.size())
        self.assertTrue(stp.isInitialStep())
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)
        stp.completed()
        self.assertFalse(stp.isInitialStep())
        self.assertEqual(stp.remaining(), 2)
        self.assertFalse(stp.isFinished())
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        stp.completed()
        self.assertEqual(stp.remaining(), 1)
        self.assertFalse(stp.isFinished())
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)
        stp.completed()
        self.assertEqual(stp.remaining(), 0)
        self.assertTrue(stp.isFinished())
        with self.assertRaisesRegex(IndexError, "no more timesteps"):
            stp.completed()

        eps = 1.0e-3
        stp = TimeStepper([0.25, 1.0, 2.0], eps)
        stp.setInitial(0.25 + eps * 0.999)
        self.assertLess(stp.getInitial() - 0.25, eps)
        self.assertSequenceEqual(stp._times, [1.0, 2.0])
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        self.assertEqual(stp.size(), 2)
        self.assertTrue(stp.isInitialStep())
        stp.completed()
        self.assertEqual(stp.remaining(), 1)
        self.assertFalse(stp.isFinished())
        self.assertFalse(stp.isInitialStep())
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)
        self.assertSequenceEqual(stp._times, [1.0, 2.0])

        stp.setInitial(0.3)
        self.assertSequenceEqual(stp._times, [1.0, 2.0])
        self.assertAlmostEqual(stp.getInitial(), 0.3)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        self.assertEqual(stp.size(), 2)
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)

        stp.setInitial(0.1)
        self.assertSequenceEqual(stp._times, [1.0, 2.0])
        self.assertAlmostEqual(stp.getInitial(), 0.1)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        self.assertEqual(stp.size(), 2)
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)

        stp.setInitial(0.1 + eps * 0.999)
        self.assertSequenceEqual(stp._times, [1.0, 2.0])
        stp.completed()
        self.assertEqual(stp.size(), 2)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitial(1.0)
        self.assertAlmostEqual(stp.getInitial(), 1.0)
        self.assertEqual(stp.size(), 1)
        self.assertSequenceEqual(stp._times, [2.0])
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)

        stp.setInitial(2.5)
        self.assertAlmostEqual(stp.getInitial(), 2.5)
        self.assertEqual(stp.size(), 1)
        self.assertFalse(stp.isFinished())
        stp.completed()
        self.assertTrue(stp.isFinished())

    def test02_initial(self):
        stp = TimeStepper([2.0, 4.0, 6.0, 8.0, 10.0])
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        self.assertEqual(stp.size(), 5)
        stp.setInitial(2.0)
        self.assertAlmostEqual(stp.getInitial(), 2.0)
        self.assertEqual(stp.size(), 4)

        stp = TimeStepper([2.0, 4.0, 6.0, 8.0, 10.0])
        stp.setInitial(3.0)
        self.assertAlmostEqual(stp.getInitial(), 3.0)
        self.assertEqual(stp.size(), 4)

        stp = TimeStepper([2.0, 4.0, 6.0, 8.0, 10.0])
        stp.setInitial(1.0)
        self.assertAlmostEqual(stp.getInitial(), 1.0)
        self.assertEqual(stp.size(), 5)

        stp = TimeStepper([2.0, 4.0, 6.0, 8.0, 10.0])
        stp.setInitial(2.0)
        stp.setInitial(1.0)
        self.assertEqual(stp.size(), 4)

    def test03_final(self):
        eps = 1.0e-3
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0], eps)
        stp.setInitial(0.0)
        stp.setFinal(2.0 - eps * 0.999)
        self.assertLess(stp.getFinal() - 2.0, eps)
        self.assertEqual(stp.size(), 3)
        stp.setFinal(2.0 + eps * 0.999)
        self.assertLess(stp.getFinal() - 2.0, eps)
        self.assertEqual(stp.size(), 3)

        stp.setFinal(1.9)
        self.assertAlmostEqual(stp.getFinal(), 1.9)
        self.assertEqual(stp.size(), 3)
        self.assertSequenceEqual(stp._times, [0.25, 1.0, 1.9])
        for _ in range(3):
            stp.completed()
        with self.assertRaisesRegex(IndexError, "no more timesteps"):
            stp.completed()
        with self.assertRaisesRegex(IndexError, "no more timesteps"):
            step = stp.getCurrent()

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitial(0.0)
        stp.setFinal(2.5)
        self.assertEqual(stp.size(), 4)
        self.assertSequenceEqual(stp._times, [0.25, 1.0, 2.0, 2.5])
        for _ in range(3):
            stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.5)
        stp.completed()
        self.assertTrue(stp.isFinished())

        stp = TimeStepper([0.0, 0.25, 1.0, 2.0])
        stp.setInitial(0.0)
        self.assertEqual(stp.size(), 3)
        stp.setFinal(0.8)
        self.assertEqual(stp.size(), 2)
        stp.setFinal()
        self.assertEqual(stp.size(), 2)

    def test04_basic(self):
        stp = TimeStepper([0.0, 0.25, 1.0, 2.0, 3.0])
        self.assertEqual(stp.size(), 4)
        self.assertFalse(stp.isFinished())

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.25)
        delta_t = stp.getIncrement()
        self.assertEqual(delta_t, 0.25)
        stp.completed()

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        delta_t = stp.getIncrement()
        self.assertAlmostEqual(delta_t, 0.75)
        stp.completed()

        step = stp.getCurrent()
        self.assertAlmostEqual(step, 2.0)
        delta_t = stp.getIncrement()
        self.assertAlmostEqual(delta_t, 1.0)
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 3.0)
        delta_t = stp.getIncrement()
        self.assertAlmostEqual(delta_t, 1.0)
        stp.completed()

        other = stp.copy()
        self.assertEqual(other.size(), 4)
        self.assertFalse(other.isFinished())
        step = other.getCurrent()
        self.assertAlmostEqual(step, 0.25)
        delta_t = other.getIncrement()
        self.assertEqual(delta_t, 0.25)

        self.assertTrue(stp.isFinished())
        with self.assertRaisesRegex(IndexError, "no more timesteps"):
            stp.completed()

    def test05_cmp(self):
        stp = TimeStepper([1.0], epsilon=0.01)
        self.assertEqual(stp.cmp(1.0, 1.0001), 0)
        self.assertEqual(stp.cmp(1.0, 1.1), -1)
        self.assertEqual(stp.cmp(1.2, 1.1), 1)

    def test06_epsilon(self):
        values = [1.0, 2.0, 2.001, 2.002, 3.0]
        stp = TimeStepper(values)
        self.assertEqual(stp.size(), 5)
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            TimeStepper(values, epsilon=0.01)

        list2 = DEFI_LIST_REEL(VALE=values)
        stp = TimeStepper.from_keywords(LIST_INST=list2, INST_INIT=None, PRECISION=1.0e-6)
        self.assertEqual(stp.size(), 5)
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            TimeStepper.from_keywords(LIST_INST=list2, INST_INIT=None, PRECISION=1.0e-2)

    def test07_meca_statique(self):
        stp = TimeStepper([0.0], initial=None)
        self.assertEqual(stp.size(), 1)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 0.0)
        self.assertIsNone(stp.getPrevious())
        self.assertIsNone(stp.getIncrement())
        stp.completed()
        self.assertTrue(stp.isFinished())

    def test08_ther_lineaire(self):
        stp = TimeStepper.from_keywords(LIST_INST=self.list0, INST_INIT=None, PRECISION=1.0e-6)
        self.assertEqual(stp.size(), 1)
        self.assertAlmostEqual(stp.getInitial(), None)
        self.assertAlmostEqual(stp.getFinal(), 0.0)

        stp = TimeStepper.from_keywords(LIST_INST=self.listr, INST_FIN=5.0, PRECISION=1.0e-6)
        self.assertEqual(stp.size(), 5)
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        self.assertAlmostEqual(stp.getFinal(), 5.0)
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        stp.completed()
        self.assertEqual(stp.remaining(), 4)

        stp = TimeStepper([-0.5, 5.5, 11.5], initial=-0.5)
        self.assertEqual(stp.size(), 2)
        self.assertAlmostEqual(stp.getInitial(), -0.5)
        self.assertAlmostEqual(stp.getFinal(), 11.5)

        stp = TimeStepper.from_keywords(
            LIST_INST=self.listr, NUME_INST_INIT=0, NUME_INST_FIN=1, PRECISION=1.0e-6
        )
        self.assertEqual(stp.size(), 1)
        self.assertAlmostEqual(stp.getInitial(), 0.0)
        self.assertAlmostEqual(stp.getFinal(), 1.0)

    def test09_from_timeslist(self):
        stp = TimeStepper.from_keywords(LIST_INST=self.listi)
        self.assertEqual(stp.size(), 1)
        self.assertAlmostEqual(stp.getInitial(), 1.0)
        self.assertAlmostEqual(stp.getFinal(), 2.0)

        list2 = self.listi.copy()
        list2.setValues([-1.0, 0.0, 1.0])
        stp = TimeStepper.from_keywords(LIST_INST=list2)
        self.assertEqual(stp.size(), 2)
        self.assertAlmostEqual(stp.getInitial(), -1.0)
        self.assertAlmostEqual(stp.getFinal(), 1.0)

    def test10_skip(self):
        stp = TimeStepper([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 2.0)
        stp._skip_before(4.5)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 5.0)

    def test20_event(self):
        stp = TimeStepper([1.0, 1.1, 2.0])
        stp.register_event(TimeStepper.Interrupt(TimeStepper.Error()))
        with self.assertRaisesRegex(CA.ConvergenceError, "MESSAGEID"):
            stp.failed(CA.ConvergenceError("MESSAGEID"))

        stp = TimeStepper([1.0, 1.1, 2.0])
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 1.1)

    def test21_split(self):
        stp = TimeStepper([0.0, 1.0, 1.1, 2.0, 3.0])
        stp._maxLevel = 3
        self.assertEqual(stp.size(), 4)
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        stp.register_event(TimeStepper.Split(TimeStepper.Error(), nbSubSteps=2, minStep=0.05))
        # print("\n+ split #1")
        stp.failed(CA.ConvergenceError("MESSAGEID"))
        # [0.5, 1.0, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 5)
        self.assertAlmostEqual(stp.getCurrent(), 0.5)
        # print("+ split #2")
        stp.failed(CA.ConvergenceError("MESSAGEID"))
        # [0.25, 0.5, 1.0, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 6)
        self.assertAlmostEqual(stp.getCurrent(), 0.25)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 0.5)
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        # print("+ split #1")
        stp.failed(CA.ConvergenceError("MESSAGEID"))
        # [0.25, 0.5, 0.75, 1.0, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 7)
        self.assertAlmostEqual(stp.getCurrent(), 0.75)
        # print("+ split #2")
        stp.failed(CA.ConvergenceError("MESSAGEID"))
        # [0.25, 0.5, 0.625, 0.75, 1.0, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 8)
        self.assertAlmostEqual(stp.getCurrent(), 0.625)
        # print("+ split #3")
        stp.failed(CA.ConvergenceError("MESSAGEID"))
        # [0.25, 0.5, 0.5625, 0.625, 0.75, 1.0, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 9)
        self.assertAlmostEqual(stp.getCurrent(), 0.5625)
        # print("+ split #4")
        with self.assertRaisesRegex(CA.SolverError, "max.*subdivision"):
            stp.failed(CA.ConvergenceError("MESSAGEID"))
        stp.completed()
        stp.completed()
        stp.completed()
        stp.completed()
        self.assertAlmostEqual(stp.getCurrent(), 1.1)
        # print("+ split #1")
        stp.failed(CA.ConvergenceError("MESSAGEID"))
        # [0.25, 0.5, 0.5625, 0.625, 0.75, 1.0, 1.05, 1.1, 2.0, 3.0]
        self.assertEqual(stp.size(), 10)
        self.assertAlmostEqual(stp.getCurrent(), 1.05)
        # print("+ split #2")
        with self.assertRaisesRegex(CA.SolverError, "trop petit"):
            stp.failed(CA.ConvergenceError("MESSAGEID"))

    def test22_finalize(self):
        stp = TimeStepper([0.0, 1.0, 1.1, 2.0, 3.0])
        nb_max_steps = 2
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        self.assertFalse(stp.isFinished())
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.1)
        self.assertFalse(stp.isFinished())
        done = stp.size() - stp.remaining()
        self.assertEqual(done, 1)
        self.assertGreaterEqual(done + 1, nb_max_steps)
        stp.setFinal(step, current=step)
        self.assertEqual(stp.size(), 2)
        self.assertSequenceEqual(stp._times, [1.0, 1.1])
        stp.completed()
        self.assertEqual(stp.remaining(), 0)
        self.assertTrue(stp.isFinished())

    def test23_event_finalize(self):
        stp = TimeStepper([0.0, 1.0, 1.1, 2.0, 3.0])
        nb_max_steps = 2
        stp.register_event(TimeStepper.Finalize(TimeStepper.MaximumNbOfSteps(nb_max_steps)))
        phys = MagicMock(name="phys_state")
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.0)
        self.assertTrue(stp.check_event(phys))
        self.assertFalse(stp.isFinished())
        stp.completed()
        step = stp.getCurrent()
        self.assertAlmostEqual(step, 1.1)
        self.assertTrue(stp.check_event(phys))
        self.assertFalse(stp.isFinished())
        stp.completed()
        self.assertEqual(stp.remaining(), 0)
        self.assertTrue(stp.isFinished())

    def test30_autosplit(self):
        residuals = [
            0.54505899475475184,
            0.28564664846422366,
            0.54818981044729176,
            0.28745504803237071,
            0.54742632449177908,
            0.28701074524171433,
            0.54761426225949084,
        ]
        xdet, xa0, xa1 = TimeStepper.AutoSplit._extrapol(residuals)
        self.assertAlmostEqual(xdet, 10.0202568758)
        self.assertAlmostEqual(xa0, 68.5865097516)
        self.assertAlmostEqual(xa1, 2.8476138960)
        crit = dict(RESI_GLOB_RELA=1.0e-10, ITER_GLOB_MAXI=10)
        nbSteps, ratio = TimeStepper.AutoSplit._splittingRatio(residuals, crit)
        self.assertEqual(nbSteps, 4)
        self.assertAlmostEqual(ratio, 0.14285714285)

    def test40_auto(self):
        times = [-1.0, 0.0, 1.0, 10.0]
        stp = TimeStepper(times, initial=None)
        stp.setInitial(-1.0)
        on_error = TimeStepper.Error()
        stp._maxLevel = 3
        split = TimeStepper.Split(on_error, nbSubSteps=4, minStep=1.0e-2)
        stp.register_event(split)

        stp.register_event(TimeStepper.Finalize(TimeStepper.MaximumNbOfSteps(100000)))
        always = TimeStepper.Always()
        mult = TimeStepper.AdaptConst(always, factor=2.0)
        stp.register_event(mult)

        self.assertAlmostEqual(stp.getInitial(), -1.0)

        self.assertAlmostEqual(stp.getCurrent(), 0.0)
        self.assertEqual(stp.size(), 3)
        self.assertTrue(all(stp._forced))
        self.assertAlmostEqual(stp.getIncrement(), 1.0)
        self.assertAlmostEqual(stp.getNextIncrement(), 1.0)
        do_step(stp, True, "t=0")

        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        self.assertEqual(stp.size(), 3)
        self.assertTrue(all(stp._forced))
        self.assertAlmostEqual(stp.getIncrement(), 1.0)
        self.assertAlmostEqual(stp.getNextIncrement(), 9.0)
        do_step(stp, False, "error at t=1.0")

        self.assertAlmostEqual(stp.getCurrent(), 0.25)
        self.assertEqual(stp.size(), 6)
        self.assertEqual(stp._forced.count(False), 3)
        self.assertAlmostEqual(stp.getIncrement(), 0.25)
        self.assertAlmostEqual(stp.getNextIncrement(), 0.75)
        do_step(stp, True, "ok at t=0.25")

        self.assertAlmostEqual(stp.getCurrent(), 0.5)
        self.assertEqual(stp.size(), 6)
        self.assertEqual(stp._forced.count(False), 2)
        self.assertAlmostEqual(stp.getIncrement(), 0.25)
        self.assertAlmostEqual(stp.getNextIncrement(), 0.5)
        do_step(stp, True, "ok at t=0.5 accel limited to dt=0.25 to maximize next step")

        self.assertAlmostEqual(stp.getCurrent(), 1.0)
        # 0.75 has been removed
        self.assertEqual(stp.size(), 5)
        self.assertTrue(all(stp._forced))
        self.assertAlmostEqual(stp.getIncrement(), 0.5)
        self.assertAlmostEqual(stp.getNextIncrement(), 9.0)
        do_step(stp, True, "ok at t=1.0 accel with dt=0.5")

        self.assertAlmostEqual(stp.getCurrent(), 2.0)
        # 2.0 has been inserted
        self.assertEqual(stp.size(), 6)
        self.assertEqual(stp._forced.count(False), 1)
        self.assertAlmostEqual(stp.getIncrement(), 1.0)
        self.assertAlmostEqual(stp.getNextIncrement(), 8.0)
        do_step(stp, True, "ok at t=2.0 accel with dt=1.0")

        self.assertAlmostEqual(stp.getCurrent(), 4.0)
        # 4.0 has been inserted
        self.assertEqual(stp.size(), 7)
        self.assertEqual(stp._forced.count(False), 1)
        self.assertAlmostEqual(stp.getIncrement(), 2.0)
        self.assertAlmostEqual(stp.getNextIncrement(), 6.0)
        do_step(stp, True, "ok at t=4.0 accel with dt=2.0")

        self.assertAlmostEqual(stp.getCurrent(), 6.0)
        # 6.0 has been inserted
        self.assertEqual(stp.size(), 8)
        self.assertEqual(stp._forced.count(False), 1)
        self.assertAlmostEqual(stp.getIncrement(), 2.0)
        self.assertAlmostEqual(stp.getNextIncrement(), 4.0)
        do_step(stp, True, "ok at t=6.0 accel with dt=2.0")

        self.assertAlmostEqual(stp.getCurrent(), 10.0)
        self.assertEqual(stp.size(), 8)
        self.assertTrue(all(stp._forced))
        self.assertAlmostEqual(stp.getIncrement(), 4.0)
        self.assertIsNone(stp.getNextIncrement())
        do_step(stp, True, "ok at t=10.0 accel with dt=4.0")

        self.assertTrue(stp.isFinished())
        with self.assertRaisesRegex(IndexError, "no more timesteps"):
            do_step(stp, True, "ended")

    def test41_auto(self):

        stp = TimeStepper([0.0, 20.0])
        # ECHEC, /4
        on_error = TimeStepper.Error()
        stp._maxLevel = 3
        split = TimeStepper.Split(on_error, nbSubSteps=4, minStep=1.0e-2)
        stp.register_event(split)
        # ADAPTATION, FIXE +75%
        stp.register_event(TimeStepper.Finalize(TimeStepper.MaximumNbOfSteps(100000)))
        always = TimeStepper.Always()
        mult = TimeStepper.AdaptConst(always, factor=2.0)
        stp.register_event(mult)

        # "PAS INIT"
        self.assertTrue(stp.isInitialStep())
        self.assertAlmostEqual(stp.getInitial(), 0.0)

        dt0 = 1.0
        self.assertLessEqual(stp.getInitial() + dt0, stp.getCurrent())
        stp._insert(0, stp.getInitial() + dt0)
        do_step(stp, False, "t=1.0 fails")
        do_step(stp, True, "ok at t=0.25")
        do_step(stp, True, "ok at t=0.75")
        do_step(stp, True, "ok at t=1.75")
        do_step(stp, True, "ok at t=3.75")
        do_step(stp, True, "ok at t=7.75")
        maxdt = (20.0 - 7.75) / 3.0
        do_step(stp, True, f"ok at t={7.75 + maxdt}")
        do_step(stp, True, "ok at t=20.0")
        self.assertTrue(stp.isFinished())
        self.assertSequenceEqual(stp._times, [0.25, 0.75, 1.75, 3.75, 7.75, 7.75 + maxdt, 20.0])


def do_step(stp, ok: bool, label: str = "", verbose: bool = False):
    if verbose:
        print(f"+ step at {stp.getCurrent()}: {label}")
        print(f"  increment: {stp.getIncrement()}, next incr: {stp.getNextIncrement()}")
        print(f"  list: {stp._times}")
        print(f"  forced: {stp._forced}")
    if ok:
        state = MagicMock(name="phys_state")
        chk = stp.check_event(state)
        assert chk, "should not occur here!"
        stp.completed()
    else:
        stp.failed(CA.ConvergenceError("MESSAGEID"))
        return False
    return True


class TestPhysicalState(unittest.TestCase):
    """Check for PhysicalState"""

    class FakeField:
        def __init__(self, value):
            self.value = value

        def __add__(self, other):
            add = self.copy()
            add.value += other.value
            return add

        def copy(self):
            return TestPhysicalState.FakeField(self.value)

        def getValues(self):
            return [self.value]

        def setValues(self, value):
            self.value = value

    def test01_stash(self):
        phys = PhysicalState(size=2, pb_type=ProblemType.MecaStat)
        # []
        self.assertEqual(len(phys._stack), 0)
        self.assertIsNone(phys.time_prev)
        self.assertIsNone(phys.U_t)
        self.assertIsNone(phys._stash)

        phys.stash()
        # set values at t=0
        phys.time_prev = 0.0
        phys.time_step = 0.0
        # can not use the setters with FakeField object
        phys.current.U_t = TestPhysicalState.FakeField(99.0)
        phys.current.deltaU = TestPhysicalState.FakeField(1.0)
        self.assertEqual(phys.time_prev, 0.0)
        self.assertEqual(phys.U_t.value, 99.0)
        self.assertEqual(phys.deltaU.value, 1.0)
        self.assertEqual(phys.U.value, 100.0)
        self.assertIsNotNone(phys._stash)
        self.assertIsNone(phys._stash._time_prev)
        self.assertIsNone(phys._stash.U_t)
        # valid this state
        phys.commit()
        self.assertEqual(phys.U_t.value, 100.0)
        self.assertIsNotNone(phys.deltaU)
        self.assertEqual(phys.deltaU.value, 0.0)
        # [phys_t0]
        self.assertEqual(len(phys._stack), 1)

        phys.stash()
        phys.time_prev = 1.0
        phys.current.U_t = TestPhysicalState.FakeField(123.0)
        self.assertEqual(phys.time_prev, 1.0)
        self.assertEqual(phys.time_curr, 1.0)
        self.assertEqual(phys.U_t.value, 123.0)
        self.assertIsNotNone(phys._stash)
        self.assertEqual(phys._stash._time_prev, 0.0)
        self.assertEqual(phys._stash.U_t.value, 100.0)

        phys.revert()
        self.assertEqual(phys.time_prev, 0.0)
        self.assertEqual(phys.U_t.value, 100.0)
        self.assertIsNone(phys._stash)

        phys.time_prev = 1.0
        phys.current.deltaU = TestPhysicalState.FakeField(99.0)
        self.assertEqual(phys.time_prev, 1.0)
        # valid this state
        phys.commit()
        self.assertEqual(phys.U_t.value, 199.0)
        # [phys_t0, phys_t1]
        self.assertEqual(len(phys._stack), 2)

    def test02_stack(self):
        phys = PhysicalState(size=3, pb_type=ProblemType.MecaStat)
        # []
        self.assertEqual(len(phys._stack), 0)
        self.assertIsNone(phys.time_prev)
        self.assertIsNone(phys.U_t)
        self.assertIsNone(phys._stash)

        # set values at t=0
        phys.time_prev = 0.0
        phys.time_step = 0.0
        # can not use the setters with FakeField object
        phys.current.U_t = TestPhysicalState.FakeField(99.0)
        phys.current.deltaU = TestPhysicalState.FakeField(1.0)
        phys.commit()
        self.assertEqual(phys.U_t.value, 100.0)
        # [phys_t0]
        self.assertEqual(len(phys._stack), 1)

        phys.time_prev = 1.0
        phys.current.deltaU = TestPhysicalState.FakeField(33.0)
        self.assertEqual(phys.time_prev, 1.0)
        phys.commit()
        self.assertEqual(phys.U_t.value, 133.0)
        # [phys_t0, phys_t1]
        self.assertEqual(len(phys._stack), 2)

        phys.time_prev = 2.0
        phys.current.deltaU = TestPhysicalState.FakeField(67.0)
        self.assertEqual(phys.time_prev, 2.0)
        phys.commit()
        self.assertEqual(phys.U_t.value, 200.0)
        # [phys_t0, phys_t1, phys_t2]
        self.assertEqual(len(phys._stack), 3)

        phys.time_prev = 3.0
        phys.current.deltaU = TestPhysicalState.FakeField(22.0)
        self.assertEqual(phys.time_prev, 3.0)
        phys.commit()
        self.assertEqual(phys.U_t.value, 222.0)
        # [phys_t1, phys_t2, phys_t3]
        self.assertEqual(len(phys._stack), 3)

        # phys.debugPrint("\n")
        # phys.current.debugPrint("current: ")


class TestConvergenceManager(unittest.TestCase):
    """Check for ConvergenceManager"""

    def test01_param(self):
        resi = ConvergenceManager.Parameter.factory("RESI_GLOB_RELA", 1.0e-6)
        iter = ConvergenceManager.Parameter.factory("ITER_GLOB_MAXI", 10)
        with self.assertRaisesRegex(TypeError, "unknown parameter"):
            ConvergenceManager.Parameter.factory("PARA", 0.0)

        self.assertFalse(resi.isSet())
        self.assertFalse(resi.isConverged())
        self.assertFalse(resi.isFinished())
        resi.value = 1.123e-4
        self.assertTrue(resi.isSet())
        self.assertFalse(resi.isConverged())
        self.assertFalse(resi.isFinished())
        resi.reset()
        self.assertFalse(resi.isSet())
        resi.value = 1.123e-7
        self.assertTrue(resi.isSet())
        self.assertTrue(resi.isConverged())
        self.assertTrue(resi.isFinished())

        self.assertFalse(iter.isSet())
        self.assertTrue(iter.isConverged())
        self.assertFalse(iter.isFinished())
        iter.value = 0
        self.assertTrue(iter.isSet())
        self.assertTrue(iter.isConverged())
        self.assertFalse(iter.isFinished())
        iter.minValue = 1
        self.assertFalse(iter.isConverged())
        self.assertFalse(iter.isFinished())
        iter.value = 1
        self.assertTrue(iter.isConverged())
        self.assertFalse(iter.isFinished())
        iter.value = 9
        self.assertTrue(iter.isConverged())
        self.assertFalse(iter.isFinished())
        iter.value = 10
        self.assertTrue(iter.isConverged())
        self.assertTrue(iter.isFinished())

    def test02_conv(self):
        conv = ConvergenceManager()
        resi_rela = conv.setdefault("RESI_GLOB_RELA", 1.0e-6)
        self.assertAlmostEqual(resi_rela.reference, 1.0e-6)

        self.assertIs(conv.get("PARAMETER"), ConvergenceManager.undef)
        resi_rela = conv.setdefault("RESI_GLOB_RELA", 1.0e99)
        self.assertAlmostEqual(resi_rela.reference, 1.0e-6)

        # done by evalNormResidual
        resi_rela.value = 1.123e-4
        conv.setdefault("RESI_GLOB_MAXI", 1.0e99).value = 9.87e3
        self.assertFalse(conv.isConverged())
        self.assertFalse(conv.isFinished())

        iter = conv.setdefault("ITER_GLOB_MAXI", 20)
        iter.value = 15
        self.assertFalse(conv.isConverged())
        self.assertFalse(conv.isFinished())
        iter.value = 20
        self.assertFalse(conv.isConverged())
        self.assertTrue(conv.isFinished())

        # test copy
        self.assertEqual(len(conv._param), 3)
        params = conv.getParameters()
        self.assertEqual(len(conv._param), 3)
        refk = ["ITER_GLOB_MAXI", "RESI_GLOB_MAXI", "RESI_GLOB_RELA"]
        keys = sorted(list(params.keys()))
        self.assertSequenceEqual(refk, keys)
        for key in refk:
            ref = conv._param.get(key)
            copy = params[key]
            self.assertIsNot(ref, copy)
            self.assertAlmostEqual(ref.reference, copy.reference)
            self.assertAlmostEqual(ref.value, copy.value)

    def test03_resi_geom(self):
        conv = ConvergenceManager()
        resi_rela = conv.setdefault("RESI_GLOB_RELA", 1.0e-6)
        resi_geom = conv.setdefault("RESI_GEOM", 1.0e-6)

        conv.initialize()  # RESI_GEOM will be ignored if not defined
        resi_rela.value = 1.0e-8
        self.assertFalse(conv.isConverged())

        conv.initialize("RESI_GEOM")  # RESI_GEOM will be initialized to -1.0
        resi_rela.value = 1.0e-8
        self.assertFalse(conv.isConverged())

        conv.initialize("RESI_GEOM")
        resi_rela.value = 1.0e-8
        resi_geom.value = 1.0e-8
        self.assertTrue(conv.isConverged())

        conv.initialize("RESI_GEOM")
        resi_rela.value = 1.0e-8
        resi_geom.value = ConvergenceManager.undef
        self.assertFalse(conv.isConverged())


class TestStorageManager(unittest.TestCase):
    """Check for StorageManager"""

    def setUp(self):
        # logger.setLevel(INFO)
        self.result = MagicMock(name="result")
        self.result.getNumberOfIndexes.return_value = 0
        self.prob = MagicMock(name="phys_pb")
        self.state = MagicMock(name="phys_state")
        self.field = MagicMock(name="field")

    def test01_mnl_basics(self):
        # calculation at 0., 1. stored at 0, 1
        store = StorageManager(self.result)
        self.assertIs(store.getResult(), self.result)
        self.assertEqual(store._last_idx, -999)
        self.assertEqual(store._stor_idx, -1)
        self.assertFalse(self.result.clear.called)

        self.assertTrue(store._to_be_stored(0, 0.0))
        self.assertTrue(store._to_be_stored(1, 1.0))
        self.assertFalse(store._has_successor(0))
        self.assertFalse(store._has_successor(1))

        # initial state
        done = store.storeState(0, 0.0, self.prob, self.state)
        self.assertTrue(done)
        self.assertEqual(store._last_idx, 0)
        self.assertEqual(store._stor_idx, 0)
        self.assertFalse(store._has_successor(0))

        # first step
        done = store.storeState(1, 1.0, self.prob, self.state)
        self.assertTrue(done)
        self.assertEqual(store._last_idx, 1)
        self.assertEqual(store._stor_idx, 1)

    def test01_mnl_restart(self):
        self.result.getNumberOfIndexes.return_value = 2
        last_index = 1

        # restart at 1. already (re)stored at 1
        # calculation at 2., 3., 4. stored at -, 2, -
        # then force storage of 4. at 3
        store = StorageManager(self.result, PAS_ARCH=2, reused=True)
        self.assertIs(store.getResult(), self.result)
        self.assertEqual(store._last_idx, -999)
        self.assertEqual(store._stor_idx, -1)

        # first storage at index 'last_index + 1'
        store.setFirstStorageIndex(storage_index=last_index + 1)
        self.assertTrue(self.result.clear.called)
        # same indexes than at the end of the previous calculation
        self.assertEqual(store._last_idx, -999)
        self.assertEqual(store._stor_idx, 1)

        self.assertTrue(store._to_be_stored(0, 1.0))
        self.assertFalse(store._to_be_stored(1, 2.0))
        self.assertTrue(store._to_be_stored(2, 3.0))
        self.assertFalse(store._to_be_stored(3, 4.0))
        self.assertFalse(store._has_successor(0))
        self.assertFalse(store._has_successor(1))
        self.assertFalse(store._has_successor(2))
        self.assertFalse(store._has_successor(3))

        # initial state 0 at 1.0, already stored
        done = store.storeState(0, 1.0, self.prob, self.state)
        self.assertFalse(done)
        self.assertEqual(store._last_idx, 0)
        self.assertEqual(store._stor_idx, 1)

        # next step 1 at 2.0, rules not checked
        done = store.storeState(1, 2.0, self.prob, self.state)
        self.assertFalse(done)
        self.assertEqual(store._last_idx, 0)
        self.assertEqual(store._stor_idx, 1)

        # next step 2 at 3.0
        done = store.storeState(2, 3.0, self.prob, self.state)
        self.assertTrue(done)
        self.assertEqual(store._last_idx, 2)
        self.assertEqual(store._stor_idx, 2)

        # next step 3 at 4.0, rules not checked
        done = store.storeState(3, 4.0, self.prob, self.state)
        self.assertFalse(done)
        self.assertEqual(store._last_idx, 2)
        self.assertEqual(store._stor_idx, 2)

        # last step 3 at 4.0
        done = store.storeState(3, 4.0, self.prob, self.state, ignore_policy=True)
        self.assertTrue(done)
        self.assertEqual(store._last_idx, 3)
        self.assertEqual(store._stor_idx, 3)


class TestResult(unittest.TestCase):
    """Check for static methods of Result"""

    def setUp(self):
        self.listr = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=10.0, PAS=1.0))

    def test01_indexes(self):
        params = {"NUME_ORDRE": (1, 2, 3, 4, 5), "INST": (1.0, 2.0, 3.0, 4.0, 5.0)}
        # TOUT_ORDRE
        res = CA.Result.getIndexesFromKeywords(params, {"TOUT_ORDRE": "OUI"})
        self.assertSequenceEqual(res, (1, 2, 3, 4, 5))
        # NUME_ORDRE
        res = CA.Result.getIndexesFromKeywords(params, {"NUME_ORDRE": [1, 3, 5]})
        self.assertSequenceEqual(res, (1, 3, 5))
        res = CA.Result.getIndexesFromKeywords(params, {"NUME_ORDRE": [0, 1, 3, 5, 8, 12]})
        self.assertSequenceEqual(res, (1, 3, 5))
        # LIST_ORDRE
        lord = DEFI_LIST_ENTI(VALE=(2, 4))
        res = CA.Result.getIndexesFromKeywords(params, {"LIST_ORDRE": lord})
        self.assertSequenceEqual(res, (2, 4))
        # INST
        res = CA.Result.getIndexesFromKeywords(params, {"INST": [1.0, 3.0, 5.0]})
        self.assertSequenceEqual(res, (1, 3, 5))
        # LIST_INST
        res = CA.Result.getIndexesFromKeywords(params, {"LIST_INST": self.listr})
        self.assertSequenceEqual(res, (1, 2, 3, 4, 5))


if __name__ == "__main__":
    CA.init()
    unittest.main()
