globals().clear()

import os
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


# --------------------------------------------------
# imports input file: replace with the desired input file located in input directory
# --------------------------------------------------
from input.buckling import gmshfile, materials, solid, dbc, dbcn, tbc, LambdaT, \
                           hist_mid_node, hist_disp_node, toll, maxiter, H, L



gmshfile = "input/" + gmshfile

import readgmsh
from B3elemental import B3_Fe_surf
from T6elemental import T6_g2n
from T6elementalNL import T6_KeNL, T6_SgNL
import elasticity_postgmsh


ndof = 2
Sdim = 3
ne = 6

analysis, T6, B3, nodes = readgmsh.readgmsh(
    gmshfile, materials, solid, dbc, dbcn, tbc
)


# --------------------------------------------------
# Sets global numbering of unknowns
# --------------------------------------------------
neq = 0
for e in range(analysis.NT6):
    for k in range(ne):
        n = T6.nodes[e, k]
        for d in range(ndof):
            if nodes.dof[n, d] == -1:
                nodes.dof[n, d] = neq
                neq += 1

print("Number of elements:", analysis.NT6)
print("Number of nodes:", analysis.NN)
print("Number of unknowns:", neq)


# --------------------------------------------------
# Assemble reference external surface loads once
# --------------------------------------------------
print("Assembling surface loads...")

Fref = np.zeros(neq, dtype=float)

for e in range(analysis.NB3):
    conn = B3.nodes[e, :]
    Xe = nodes.coor[conn, :]
    Ge = nodes.dof[conn, :].reshape(-1)
    loadxy = B3.loadxy[e, :]
    press = B3.press[e]

    Fe = B3_Fe_surf(Xe, loadxy, press)

    Le0 = Ge >= 0
    Fref[Ge[Le0]] += Fe[Le0]


# --------------------------------------------------
# Initialize nonlinear analysis
# --------------------------------------------------
UD = nodes.U.copy()       # prescribed displacement pattern read from input
nodes.U[:, :] = 0.0       # current total displacement starts from zero

#LambdaT = np.asarray(LambdaT, dtype=float).reshape(-1)
numstep = len(LambdaT) - 1
hist = np.zeros((numstep, 3), dtype=float)       # lambda, iterations, residual
hist_out = np.zeros((numstep, 2), dtype=float)   # mid-span displacement, load estimate

print("Nonlinear elastic analysis")
print("Total steps:", numstep)


# --------------------------------------------------
# Loading sequence and Newton iterations
# --------------------------------------------------
for step in range(numstep):
    lambda_old = LambdaT[step]
    lambda_new = LambdaT[step + 1]
    dlambda = lambda_new - lambda_old

    # Same convention as geomNL.m:
    # nodes(n).U = nodes(n).U + Dlambda * nodes(n).UD
    nodes.U += dlambda * UD

    # Same convention as geomNL.m: Fext = F, not lambda_new * F.
    Fext = Fref.copy()

    iter_count = 0
    resid = 1.0
    residref = 1.0

    while resid > toll * residref:
        iter_count += 1

        # --------------------------------------------------
        # Assemble tangent matrix and internal force vector
        # --------------------------------------------------
        Fint = np.zeros(neq, dtype=float)
        rows_all = []
        cols_all = []
        data_all = []

        for e in range(analysis.NT6):
            conn = T6.nodes[e, :]
            Xe = nodes.coor[conn, :]
            Ge = nodes.dof[conn, :].reshape(-1)
            Ue = nodes.U[conn, :].reshape(-1)
            mate = materials[T6.material[e]]

            Ke, Finte = T6_KeNL(Xe, mate, Ue)

            Le0 = Ge >= 0
            Ie = Ge[Le0]
            nfree = Ie.size

            if nfree:
                Fint[Ie] += Finte[Le0]
                Ke00 = Ke[Le0][:, Le0]
                rows_all.append(np.repeat(Ie, nfree))
                cols_all.append(np.tile(Ie, nfree))
                data_all.append(Ke00.ravel())

        rows = np.concatenate(rows_all)
        cols = np.concatenate(cols_all)
        data = np.concatenate(data_all)
        K = coo_matrix((data, (rows, cols)), shape=(neq, neq)).tocsc()
        K.sum_duplicates()

        # --------------------------------------------------
        # Newton correction
        # --------------------------------------------------
        R = -Fext - Fint
        resid = float(np.linalg.norm(R))
        if iter_count == 1:
            residref = max(resid, 1.0)

        if resid <= toll * residref:
            break

        Du = spsolve(K, -R)

        free = nodes.dof >= 0
        nodes.U[free] += Du[nodes.dof[free]]

        print(
            f"Step: {step + 1:4d}  "
            f"Iter: {iter_count:3d}  "
            f"Lambda: {lambda_new:.6g}  "
            f"Residuum: {resid:.6e}"
        )

        if iter_count >= maxiter:
            raise RuntimeError(
                f"Newton did not converge at step {step + 1}; "
                f"residual = {resid:.6e}, tolerance = {toll * residref:.6e}."
            )

    hist[step, :] = [lambda_new, iter_count, resid]

    # --------------------------------------------------
    # Force computation for buckling history
    # --------------------------------------------------
    work = 0.0

    for e in range(analysis.NT6):
        conn = T6.nodes[e, :]
        Xe = nodes.coor[conn, :]
        Ue = nodes.U[conn, :].reshape(-1)
        mate = materials[T6.material[e]]

        Ke, Finte = T6_KeNL(Xe, mate, Ue)
        work += float(Ue @ Finte)

    Ud = nodes.U[hist_disp_node - 1, 1]
    if abs(Ud) > 0.0:
        hist_out[step, :] = [nodes.U[hist_mid_node - 1, 0], work / Ud]
    else:
        hist_out[step, :] = [nodes.U[hist_mid_node - 1, 0], np.nan]


# --------------------------------------------------
# Stress evaluation first at Gauss nodes and then extrapolates to nodes
# --------------------------------------------------
Sn = np.zeros((analysis.NN, Sdim), dtype=float)
counter = np.zeros(analysis.NN, dtype=int)

for e in range(analysis.NT6):
    conn = T6.nodes[e, :]
    Xe = nodes.coor[conn, :]
    Ue = nodes.U[conn, :].reshape(-1)
    mate = materials[T6.material[e]]

    GLg, Sg = T6_SgNL(Xe, mate, Ue)
    Sne = T6_g2n(Sg)
    Sn[conn, :] += Sne
    counter[conn] += 1

visited = counter > 0
Sn[visited, :] = Sn[visited, :] / counter[visited].reshape(-1, 1)


# --------------------------------------------------
# Post processing in Gmsh
# --------------------------------------------------
elasticity_postgmsh.postgmsh(gmshfile, analysis, nodes, Sn)


# --------------------------------------------------
# Post-processing: compares with buckling force from beam theory
# --------------------------------------------------
os.makedirs("output", exist_ok=True)
np.savetxt(
    "output/buckling_history.txt",
    hist_out,
    header="mid_span_displacement  compressive_load",
)

import matplotlib.pyplot as plt

plt.figure()
plt.title("Buckling of doubly clamped beam")
plt.xlabel("Mid-span deflection")
plt.ylabel("Compressive load")
plt.plot(np.abs(hist_out[:, 0]), np.abs(hist_out[:, 1]), "ok-", label="Mesh1")

Jbeam = H**3 / 12.0
force = 4.0 * np.pi**2 * Jbeam / L**2
plt.plot(
    np.abs(hist_out[:, 0]),
    force * np.ones(numstep),
    "k-.",
    label="Theory",
)
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig("output/buckling_comparison.png", dpi=200)
plt.show()

print("Buckling history written to output/buckling_history.txt")
print("Buckling plot written to output/buckling_comparison.png")
print("Done")
