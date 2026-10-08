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


from code_aster.Commands import *
from code_aster import CA

import os
import json

from math import fsum, pi, sqrt, log

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True

except ImportError:
    HAS_MATPLOTLIB = False


# Convergence du couplage CG(P2)-HHO(k=1) sur une interface coincidente.
# Solution manufacturee ; erreur discrete en energie et ordre log-log.

CA.init("--test")
test = CA.TestCase()

with open("fort.10") as datafile:
    parameters = json.load(datafile)

E = 2000.0
GAMMA_FACTOR = 0.1  # COEF_PENA = GAMMA_FACTOR * E
NU = 0.30
MU = E / (2.0 * (1.0 + NU))
LAM = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))

GROUP_CG = ("Haut", "LIAISON_Haut", "Pres", "Sym_Haut", "Droit_Haut")
GROUP_HHO = ("Bas", "LIAISON_Bas", "Encast", "Sym_Bas", "Droit_Bas")
GROUP_CG_MODEL = GROUP_CG
GROUP_HHO_MODEL = GROUP_HHO
GROUP_VOL = ("Haut", "Bas")
GROUP_BOUNDARY_CG = ("Pres", "Sym_Haut", "Droit_Haut")
GROUP_BOUNDARY_HHO = ("Encast", "Sym_Bas", "Droit_Bas")
GROUP_VOL_CG = "Haut"
GROUP_VOL_HHO = "Bas"
GROUP_ORIE = (
    "LIAISON_Haut",
    "LIAISON_Bas",
    "Encast",
    "Pres",
    "Sym_Haut",
    "Sym_Bas",
    "Droit_Haut",
    "Droit_Bas",
)

BASE_UNIT = 20
NLEVELS = parameters["number_of_refinements"]


def diameter_cell(coor, nodes):
    nb_nodes = len(nodes)
    pts = np.array([coor.getNode(node).getValues() for node in nodes])

    return max(
        np.linalg.norm(pts[i] - pts[j]) for i in range(nb_nodes) for j in range(i + 1, nb_nodes)
    )


def diameter(mesh):
    mesh_lin = mesh.convertToLinear()
    nbCells = mesh_lin.getNumberOfCells()
    coor = mesh_lin.getCoordinates()
    connec = mesh_lin.getConnectivity()

    sum = 0.0
    for c_id in range(nbCells):
        sum += diameter_cell(coor, connec[c_id])

    diam_avg = sum / nbCells
    return diam_avg


def local_order(h, e):
    """Ordres entre deux maillages successifs, du grossier au fin."""
    return [float("nan")] + [
        np.log(e[i] / e[i - 1]) / np.log(h[i] / h[i - 1]) for i in range(1, len(h))
    ]


def _unique_nodes_from_groups(mesh, groups):
    connec = mesh.getConnectivity()
    nodes = set()
    for group in groups:
        for cell_id in mesh.getCells(group):
            for node in connec[cell_id]:
                nodes.add(node)
    return list(nodes)


def _coord_bounds(mesh, groups):
    coor = mesh.getCoordinates()
    nodes = _unique_nodes_from_groups(mesh, groups)
    xs, ys = [], []
    for node in nodes:
        point = coor.getNode(node)
        xs.append(point.x())
        ys.append(point.y())
    return min(xs), max(xs), min(ys), max(ys)


def build_manufactured_formulas(mesh):
    xs0, xs1, _, _ = _coord_bounds(mesh, ("Sym_Haut", "Sym_Bas"))
    _, _, yb0, yb1 = _coord_bounds(mesh, ("Encast",))
    _, _, yt0, yt1 = _coord_bounds(mesh, ("Pres",))
    xmin, xmax, _, _ = _coord_bounds(mesh, GROUP_VOL)

    x_sym = 0.5 * (xs0 + xs1)
    y_bot = 0.5 * (yb0 + yb1)
    y_top = 0.5 * (yt0 + yt1)
    lx = max(abs(xmax - x_sym), abs(xmin - x_sym))
    alpha = pi / lx
    beta = 2.0 * alpha

    s = f"(X-({x_sym:.17e}))"
    B = f"((Y-({y_bot:.17e}))*(Y-({y_top:.17e})))"
    Bp = f"(2*Y-({(y_bot + y_top):.17e}))"
    Bpp = "2.0"
    A = f"({s}*sin(({alpha:.17e})*{s}))"
    Ap = f"(sin(({alpha:.17e})*{s})+({alpha:.17e})*{s}*cos(({alpha:.17e})*{s}))"
    App = f"(2*({alpha:.17e})*cos(({alpha:.17e})*{s})-({alpha:.17e})*({alpha:.17e})*{s}*sin(({alpha:.17e})*{s}))"
    C = f"({s}*sin(({beta:.17e})*{s}))"
    Cp = f"(sin(({beta:.17e})*{s})+({beta:.17e})*{s}*cos(({beta:.17e})*{s}))"
    Cpp = f"(2*({beta:.17e})*cos(({beta:.17e})*{s})-({beta:.17e})*({beta:.17e})*{s}*sin(({beta:.17e})*{s}))"

    ux = f"({A})*({B})"
    uy = f"({C})*({B})"
    fx = (
        f"-(({2 * MU + LAM:.17e})*({App})*({B})"
        f" + ({LAM + MU:.17e})*({Cp})*({Bp})"
        f" + ({MU:.17e})*({A})*({Bpp}))"
    )
    fy = (
        f"-(({LAM + MU:.17e})*({Ap})*({Bp})"
        f" + ({MU:.17e})*({Cpp})*({B})"
        f" + ({2 * MU + LAM:.17e})*({C})*({Bpp}))"
    )

    return (
        FORMULE(VALE=ux, NOM_PARA=("X", "Y")),
        FORMULE(VALE=uy, NOM_PARA=("X", "Y")),
        FORMULE(VALE=fx, NOM_PARA=("X", "Y")),
        FORMULE(VALE=fy, NOM_PARA=("X", "Y")),
    )


def solve_on_mesh(mesh, order):
    """Resout le probleme et renvoie h, l'erreur et les deux energies discretes."""
    UX_EX, UY_EX, FX_EX, FY_EX = build_manufactured_formulas(mesh)

    mesh = MODI_MAILLAGE(reuse=mesh, MAILLAGE=mesh, ORIE_PEAU=_F(GROUP_MA_PEAU=GROUP_ORIE))

    model = AFFE_MODELE(
        MAILLAGE=mesh,
        AFFE=(
            _F(GROUP_MA=GROUP_CG_MODEL, PHENOMENE="MECANIQUE", MODELISATION="D_PLAN"),
            _F(
                GROUP_MA=GROUP_HHO_MODEL,
                PHENOMENE="MECANIQUE",
                MODELISATION="D_PLAN_HHO",
                FORMULATION=order,
            ),
        ),
    )
    material = DEFI_MATERIAU(ELAS=_F(E=E, NU=NU, RHO=1.0))
    fieldmat = AFFE_MATERIAU(MAILLAGE=mesh, MODELE=model, AFFE=_F(TOUT="OUI", MATER=material))

    diri_cg = AFFE_CHAR_CINE(MODELE=model, MECA_IMPO=_F(GROUP_MA=GROUP_BOUNDARY_CG, DX=0.0, DY=0.0))
    diri_hho = AFFE_CHAR_CINE(
        MODELE=model, MECA_IMPO=_F(GROUP_MA=GROUP_BOUNDARY_HHO, DX=0.0, DY=0.0)
    )
    load_cg = AFFE_CHAR_MECA_F(
        MODELE=model, FORCE_INTERNE=_F(GROUP_MA=GROUP_VOL_CG, FX=FX_EX, FY=FY_EX)
    )
    load_hho = AFFE_CHAR_MECA_F(
        MODELE=model, FORCE_INTERNE=_F(GROUP_MA=GROUP_VOL_HHO, FX=FX_EX, FY=FY_EX)
    )
    raccord = AFFE_CHAR_MECA(
        MODELE=model,
        LIAISON_MASSIF=_F(
            GROUP_MA_ESCL="LIAISON_Haut",
            GROUP_MA_MAIT="LIAISON_Bas",
            METHODE="NITSCHE",
            COEF_PENA=GAMMA_FACTOR * E,
            APPARIEMENT="MORTAR",
            TYPE_APPA="FORCEBRUTE",
            APPA_TOLE=1.0e-8,
            AIRE_TOLE=1.0e-8,
        ),
    )

    phys_pb = CA.PhysicalProblem(model, fieldmat)
    phys_pb.addLoad(load_cg)
    phys_pb.addLoad(load_hho)
    phys_pb.addLoad(raccord)
    phys_pb.addDirichletBC(diri_cg)
    phys_pb.addDirichletBC(diri_hho)
    phys_pb.computeDOFNumbering()

    disc = CA.DiscreteComputation(phys_pb)
    me_rigi = disc.getLinearStiffnessMatrix(with_dual=True, assembly=False)
    me_cpl = disc.getMechanicalLinearCouplingMatrix(None)
    A = CA.AssemblyMatrixDisplacementReal(phys_pb)
    A.assemble([me_rigi, me_cpl], phys_pb.getListOfLoads())
    # Rigidite CG/HHO, sans les termes de couplage Nitsche.
    A_rigi = CA.AssemblyMatrixDisplacementReal(phys_pb)
    A_rigi.assemble([me_rigi], phys_pb.getListOfLoads())

    rhs = disc.getVolumetricForces()
    diri = disc.getDirichletBC()
    solver = CA.MumpsSolver()
    solver.factorize(A)
    u_h = solver.solve(rhs, diri)

    # Reduction HHO et interpolation nodale CG de la solution exacte.
    hho = CA.HHO(phys_pb)
    u_hho = hho.projectOnHHOSpace([UX_EX, UY_EX])
    h_vals, (h_nodes, h_cmps) = u_hho.getValuesWithDescription([], [])
    hho_map = {(n, c): v for n, c, v in zip(h_nodes, h_cmps, h_vals)}

    (g_nodes, g_cmps), g_rows = u_h.getDescription().getDOFsWithDescription([], [], True)
    cg_node_set = set(_unique_nodes_from_groups(mesh, GROUP_CG_MODEL))
    coor = mesh.getCoordinates()

    vals = [0.0] * u_h.size()
    for node, cmp, row in zip(g_nodes, g_cmps, g_rows):
        key = (node, cmp)
        if key in hho_map:
            vals[row] = hho_map[key]
        elif node in cg_node_set and cmp in ("DX", "DY"):
            p = coor.getNode(node)
            vals[row] = UX_EX(p.x(), p.y()) if cmp == "DX" else UY_EX(p.x(), p.y())
    u_I = u_h.copy()
    u_I.setValues(vals)

    e = u_I - u_h
    ev = e.getValues()
    qH1 = fsum(a * b for a, b in zip((A_rigi * e).getValues(), ev))

    # Penalisation : A_Gamma(2*gamma) - A_Gamma(gamma).
    raccord2 = AFFE_CHAR_MECA(
        MODELE=model,
        LIAISON_MASSIF=_F(
            GROUP_MA_ESCL="LIAISON_Haut",
            GROUP_MA_MAIT="LIAISON_Bas",
            METHODE="NITSCHE",
            COEF_PENA=2.0 * GAMMA_FACTOR * E,
            APPARIEMENT="MORTAR",
            TYPE_APPA="FORCEBRUTE",
            APPA_TOLE=1.0e-8,
            AIRE_TOLE=1.0e-8,
        ),
    )
    pp2 = CA.PhysicalProblem(model, fieldmat)
    pp2.addLoad(load_cg)
    pp2.addLoad(load_hho)
    pp2.addLoad(raccord2)
    pp2.addDirichletBC(diri_cg)
    pp2.addDirichletBC(diri_hho)
    pp2.computeDOFNumbering()
    me_cpl2 = CA.DiscreteComputation(pp2).getMechanicalLinearCouplingMatrix(None)
    A_cg = CA.AssemblyMatrixDisplacementReal(phys_pb)
    A_cg.assemble([me_cpl], phys_pb.getListOfLoads())
    A_cg2 = CA.AssemblyMatrixDisplacementReal(phys_pb)
    A_cg2.assemble([me_cpl2], phys_pb.getListOfLoads())
    pen_e = [a2 - a1 for a1, a2 in zip((A_cg * e).getValues(), (A_cg2 * e).getValues())]
    qB = fsum(p * b for p, b in zip(pen_e, ev))

    # ||I_h u - u_h||_Vh^2 : rigidite CG/HHO et penalisation d'interface.
    qVh = qH1 + qB
    err = sqrt(qVh) if qVh > 0.0 else float("nan")

    # Taille caracteristique pour le domaine rectangulaire 2D.
    h = diameter(mesh)
    return h, err, qH1, qB


error = {}
conv_order = {}

for order in ("LINEAIRE", "QUADRATIQUE"):
    MESH = LIRE_MAILLAGE(UNITE=BASE_UNIT)
    if order == "LINEAIRE":
        MESH = MESH.convertToQuadratic()
    else:
        MESH = MESH.convertToQuadratic()

    MESH = CREA_MAILLAGE(MAILLAGE=MESH, MODI_HHO=_F(GROUP_MA=GROUP_HHO))

    error[order] = {"h": [], "H1": []}
    hs, errs, qH1s, qBs = [], [], [], []
    for i_reff in range(NLEVELS):
        print("\nCONVERGENCE_CG_HHO niveau %d" % i_reff)
        h, err, qH1, qB = solve_on_mesh(MESH, order)
        hs.append(h)
        errs.append(err)
        qH1s.append(qH1)
        qBs.append(qB)
        error[order]["h"].append(h)
        error[order]["H1"].append(err)

        if i_reff < NLEVELS - 1:
            MESH = MESH.refine(1)

    hs = np.array(hs)
    errs = np.array(errs)
    qH1s = np.array(qH1s)
    qBs = np.array(qBs)
    qVhs = qH1s + qBs
    oVh = local_order(hs, errs)

    # Conditions de definition des logarithmes et des ordres locaux.
    if not (np.all(np.isfinite(hs)) and np.all(hs > 0.0)):
        raise ValueError("Tailles de maillage non finies ou non positives.")
    if not np.all(np.diff(hs) < 0.0):
        raise ValueError("Les maillages doivent etre ordonnes du grossier au fin.")
    if not (np.all(np.isfinite(errs)) and np.all(errs > 0.0)):
        raise ValueError("Erreurs non finies ou non positives dans la norme Vh.")

    # Validite de la norme et convergence observee.
    test.assertTrue(bool(np.all(np.isfinite(qH1s))))
    test.assertTrue(bool(np.all(np.isfinite(qBs))))
    test.assertTrue(bool(np.all(qVhs > 0.0)))
    test.assertTrue(bool(np.all(qBs >= -1.0e-8 * np.maximum(1.0, np.abs(qVhs)))))
    test.assertTrue(bool(np.all(np.isfinite(oVh[1:]))))
    test.assertTrue(bool(errs[-1] < errs[0]))

    conv_order[order] = {}
    for norm in ("H1",):
        xlog = [log(h) for h in error[order]["h"]]
        ylog = [log(y) for y in error[order][norm]]
        A = np.vstack([xlog, np.ones(len(xlog))]).T
        m, c = np.linalg.lstsq(A, ylog, rcond=None)[0]
        conv_order[order][norm] = [m, c]

    print("\n%-13s %-13s %-8s" % ("h", "ERR_Vh", "pVh"))
    for i in range(len(hs)):
        cvorder = "-" if i == 0 else "%.3f" % oVh[i]
        print("%-13.5e %-13.6e %-8s" % (hs[i], errs[i], cvorder))
    print("\nORDRE_GLOBAL_Vh %.6f" % conv_order[order]["H1"][0])

# compute convergence order

for order in ("LINEAIRE", "QUADRATIQUE"):

    print("Convergence order %s -> H1-norm: %f" % (order, conv_order[order]["H1"][0]))

    # test convergence order
    test.assertAlmostEqual(conv_order[order]["H1"][0], parameters["H1"][order], delta=1e-4)

# plot figure with matplotlib
if HAS_MATPLOTLIB and os.getenv("DISPLAY"):

    # disable floating point exceptions from matplotlib
    with CA.disable_fpe():

        ylim = {"H1": [1e-7, 1.0]}
        the_conv = {
            "LINEAIRE": {"L2": {"o": 3, "c": 2e-2}, "H1": {"o": 2, "c": 0.2}},
            "QUADRATIQUE": {"L2": {"o": 4, "c": 1e-3}, "H1": {"o": 3, "c": 0.02}},
        }
        # plot the data
        for norm in ("H1",):
            plt.plot(1, 1, 1)

            plt.xscale("log")
            plt.yscale("log")

            # set the limits
            plt.xlim([0.02, 2])
            plt.ylim(ylim[norm])
            plt.xlabel("mesh-size")
            plt.ylabel("%s-error" % norm)

            plt.plot(
                error["LINEAIRE"]["h"],
                error["LINEAIRE"][norm],
                marker="o",
                color="tab:blue",
                label="k=1, computed",
            )
            m, c = conv_order["LINEAIRE"][norm]
            plt.plot(
                error["LINEAIRE"]["h"],
                [
                    the_conv["LINEAIRE"][norm]["c"] * h ** the_conv["LINEAIRE"][norm]["o"]
                    for h in error["LINEAIRE"]["h"]
                ],
                "k--",
                label="k=1, theorical",
                color="tab:blue",
            )

            plt.plot(
                error["QUADRATIQUE"]["h"],
                error["QUADRATIQUE"][norm],
                marker="o",
                color="tab:orange",
                label="k=2, computed",
            )
            m, c = conv_order["QUADRATIQUE"][norm]
            plt.plot(
                error["QUADRATIQUE"]["h"],
                [
                    the_conv["QUADRATIQUE"][norm]["c"] * h ** the_conv["QUADRATIQUE"][norm]["o"]
                    for h in error["LINEAIRE"]["h"]
                ],
                "k--",
                label="k=2, theorical",
                color="tab:orange",
            )

            plt.legend()
            plt.title("%s convergence error for HHO" % norm)
            plt.show()

            # save plot
            # savedir = "/tmp/" # or os.getcwd()
            # plt.savefig(os.path.join(savedir, "%s_error.png" % norm))
            plt.clf()

test.printSummary()
FIN()
