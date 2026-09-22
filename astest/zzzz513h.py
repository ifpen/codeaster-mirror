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
import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True

except ImportError:
    HAS_MATPLOTLIB = False

from math import sqrt, log

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

###################################################################################
#
#   Analytical solution
#   Linear elasticity - axisymmetric
#
#   Weak form: (grad u, grad v) = (f,v)
#   HHO unknowns : huT = (uT, udT)
#
####################################################################################


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

    diam = -1.0
    for c_id in range(nbCells):
        diam = max(diam, diameter_cell(coor, connec[c_id]))

    return diam


# number of refinement
nb_reff = 6

E = 200000.0
Nu = 0.3

lamb = E * Nu / (1 + Nu) / (1 - 2 * Nu)
mu = E / 2 / (1 + Nu)

# define analytical solution
u_R = FORMULE(VALE="sin(pi*X)*sin(pi*Y)", NOM_PARA=("X", "Y"))
u_Z = FORMULE(VALE="cos(pi*X)*cos(pi*Y)", NOM_PARA=("X", "Y"))

# define load function
f_R = FORMULE(
    VALE="((2*pi*pi*mu*X*X*sin(pi*X) - pi*X*(lamb + 2*mu)*cos(pi*X) + (lamb + 2*mu)*sin(pi*X))*sin(pi*Y))/(X*X)",
    NOM_PARA=("X", "Y"),
    lamb=lamb,
    mu=mu,
)
f_Z = FORMULE(
    VALE="(pi*(-lamb*sin(pi*X) + 2*pi*mu*X*cos(pi*X))*cos(pi*Y))/X",
    NOM_PARA=("X", "Y"),
    lamb=lamb,
    mu=mu,
)

# error save
error = {}
conv_order = {}

# initial_mesh - domain [0,1]^2 - 4 quads
mesh0_quad = CA.Mesh.buildSquare(refine=1)
# split in triangular mesh
mesh0_tri = CREA_MAILLAGE(MAILLAGE=mesh0_quad, MODI_MAILLE=_F(TOUT="OUI", OPTION="QUAD_TRIA3"))
# convert for hho-cells
mesh0_hho = CREA_MAILLAGE(MAILLAGE=mesh0_tri, MODI_HHO=_F(TOUT="OUI"))


for order in ("LINEAIRE", "QUADRATIQUE"):
    error[order] = {"h": [], "L2": [], "H1": []}
    mesh = mesh0_hho
    for i_reff in range(nb_reff):
        ## DEFINE PROBLEM
        # create mesh - refine previous mesh
        mesh = mesh.refine(1)

        # define material
        coeff = DEFI_MATERIAU(ELAS=_F(E=E, NU=Nu, RHO=1.0), HHO=_F(COEF_STAB=2 * mu))

        # apply material on mesh
        mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

        # define finite element model
        model = AFFE_MODELE(
            MAILLAGE=mesh,
            AFFE=_F(TOUT="OUI", MODELISATION="AXIS_HHO", FORMULATION=order, PHENOMENE="MECANIQUE"),
        )

        # define Dirichlet BC
        bc = AFFE_CHAR_CINE_F(
            MODELE=model, MECA_IMPO=_F(GROUP_MA=("RIGHT", "LEFT", "TOP", "BOTTOM"), DX=u_R, DY=u_Z)
        )

        # define external load
        load = AFFE_CHAR_MECA_F(
            MODELE=model, FORCE_INTERNE=_F(GROUP_MA=("SURFACE"), FX=f_R, FY=f_Z)
        )

        ## COMPUTE DISCRETE SOLUTION
        # define physical problem
        phys_pb = CA.PhysicalProblem(model, mater)
        phys_pb.addLoad(load)
        phys_pb.addDirichletBC(bc)

        # compute DOF numbering
        phys_pb.computeDOFNumbering()

        # create discrete computation
        disc_comp = CA.DiscreteComputation(phys_pb)

        # compute rigidity matrix: (lambda * GkT(huT), GkT(hvT))_T + lambda * stab(huT, hvT)
        rigidity = disc_comp.getLinearStiffnessMatrix(assembly=True)

        # compute external load: (f, vT)_T
        rhs = disc_comp.getVolumetricForces()

        # compute Dirichlet BC to apply
        diriBCs = disc_comp.getDirichletBC()

        # define linear solver - MUMPS
        mySolver = CA.MumpsSolver()

        # factorize and solve
        mySolver.factorize(rigidity)
        u_hho = mySolver.solve(rhs, diriBCs)

        ## COMPUTE ERROR
        # create hho handler
        hho = CA.HHO(phys_pb)

        # Project analytical solution on HHO space
        u_proj = hho.projectOnHHOSpace([u_R, u_Z])

        # compute difference
        u_diff = u_hho - u_proj

        # to compute norm
        norm_L2 = disc_comp.getMassMatrix(assembly=True)
        norm_H1 = rigidity

        # compute L2 and H1-errors
        error[order]["h"].append(diameter(mesh))
        error[order]["L2"].append(sqrt((norm_L2 * u_diff).dot(u_diff)))
        error[order]["H1"].append(sqrt((norm_H1 * u_diff).dot(u_diff)))

# compute convergence order
for order in ("LINEAIRE", "QUADRATIQUE"):
    conv_order[order] = {}
    for norm in ("L2", "H1"):
        xlog = [log(h) for h in error[order]["h"]]
        ylog = [log(y) for y in error[order][norm]]
        A = np.vstack([xlog, np.ones(len(xlog))]).T
        m, c = np.linalg.lstsq(A, ylog, rcond=None)[0]
        conv_order[order][norm] = [m, c]

    print(
        "Convergence order %s -> L2-norm: %f and H1-norm: %f"
        % (order, conv_order[order]["L2"][0], conv_order[order]["H1"][0])
    )

    # test convergence order
    test.assertAlmostEqual(
        conv_order[order]["L2"][0],
        {"LINEAIRE": 2.7141075455980594, "QUADRATIQUE": 3.918285655874547}[order],
        delta=1e-4,
    )
    test.assertAlmostEqual(
        conv_order[order]["H1"][0],
        {"LINEAIRE": 1.9407444513567884, "QUADRATIQUE": 2.9844370844900094}[order],
        delta=1e-4,
    )

# plot figure with matplotlib
if HAS_MATPLOTLIB and os.getenv("DISPLAY"):

    # disable floating point exceptions from matplotlib
    with CA.disable_fpe():

        ylim = {"L2": [1e-8, 5e-2], "H1": [1e-5, 8e1]}
        the_conv = {
            "LINEAIRE": {"L2": 3, "H1": 2},
            "QUADRATIQUE": {"L2": 4, "H1": 3},
            "CUBIQUE": {"L2": 5, "H1": 4},
        }

        def compute_slope(h, y):
            xlog = np.log(h)
            ylog = np.log(y)
            A = np.vstack([xlog, np.ones(len(xlog))]).T
            m, _ = np.linalg.lstsq(A, ylog, rcond=None)[0]
            return m

        def reference_curve(h, y0, p):
            h = np.asarray(h, dtype=float)
            h0 = h[0]
            return y0 * (h / h0) ** p

        # plot the data
        for norm in ("L2", "H1"):

            fig, ax = plt.subplots()

            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim([7e-3, 0.5])
            ax.set_ylim(ylim[norm])

            ax.set_xlabel(r"$h$")
            ax.set_ylabel(rf"$\|u - u_h\|_{{{norm}}}$")

            # ----- LINEAIRE -----
            h = np.asarray(error["LINEAIRE"]["h"], dtype=float)
            y = np.asarray(error["LINEAIRE"][norm], dtype=float)

            slope = compute_slope(h, y)
            p_th = the_conv["LINEAIRE"][norm]

            ax.plot(h, y, "o-", color="tab:blue", label=rf"$k=1$ (slope $\approx {slope:.2f}$)")

            ax.plot(
                h,
                reference_curve(h, y[0], p_th),
                "--",
                color="tab:blue",
                label=rf"$\mathcal{{O}}(h^{{{p_th}}})$",
            )

            # ----- QUADRATIQUE -----
            h = np.asarray(error["QUADRATIQUE"]["h"], dtype=float)
            y = np.asarray(error["QUADRATIQUE"][norm], dtype=float)

            slope = compute_slope(h, y)
            p_th = the_conv["QUADRATIQUE"][norm]

            ax.plot(h, y, "s-", color="tab:orange", label=rf"$k=2$ (slope $\approx {slope:.2f}$)")

            ax.plot(
                h,
                reference_curve(h, y[0], p_th),
                "--",
                color="tab:orange",
                label=rf"$\mathcal{{O}}(h^{{{p_th}}})$",
            )

            ax.legend(loc="lower right")
            fig.tight_layout()

            # --- EXPORT PGF ---
            try:
                fig.savefig(f"convergence_{norm}.pgf")
                fig.savefig(f"convergence_{norm}.pdf")
            except AttributeError as exc:
                print(f"PGF/PDF export disabled: {exc}")

            plt.close(fig)

# close
CA.close()
