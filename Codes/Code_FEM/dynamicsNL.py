import time
import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import eigsh
from sksparse.cholmod import cholesky

# --------------------------------------------------
# --------------------------------------------------
from input.dynamics_tyreincl import gmshfile, materials, solid, dbc, dbcn, tbc, bf, tf, Dt, output_interval, v0, Lambda, mu
# --------------------------------------------------
# --------------------------------------------------

gmshfile="input/" + gmshfile
import readgmsh
analysis, T6, B3, nodes = readgmsh.readgmsh(gmshfile, materials, solid, dbc, dbcn, tbc)

from T6elemental import T6_Ke, T6_Me, T6_Fe
from T6elementalNL import T6_FeNL, T6_FeNL_numba
from B3elemental import B3_Fe_contact_friction
import dynamics_postgmsh

# --------------------------------------------------
# problem constants
# --------------------------------------------------
ndof = 2
DG = 2
ne = 6
nen_dof = ne * ndof                 # 12 for a T6 plane element
max_nnz_per_elem = nen_dof * nen_dof

# --------------------------------------------------
# global numbering of unknowns
# --------------------------------------------------
neq = 0
for e in range(analysis.NT6):
    conn = T6.nodes[e, :]
    for n in conn:
        for d in range(ndof):
            if nodes.dof[n, d] == -1:
                nodes.dof[n, d] = neq
                neq += 1

analysis.neq = neq
print(f"Number of unknowns: {neq}")

# --------------------------------------------------
# assemble K, M, F with preallocated triplet storage
# --------------------------------------------------
print(f"Number of elements: {analysis.NT6}")
print("Assembling matrices and body-force vector...")

t0 = time.perf_counter()
nelem = analysis.NT6
rows_all = []
cols_all = []
dataK_all = []
dataM_all = []

F = np.zeros(neq, dtype=np.float64)
ptr = 0
percold = -1

for e in range(nelem):
    perc = 20 * ((5 * (e + 1)) // nelem) if nelem > 0 else 100
    if perc != percold:
        percold = perc
        print(f"Assembling: {percold}%")

    conn = T6.nodes[e, :]
    Xe = nodes.coor[conn, :DG]
    Ge = nodes.dof[conn, :ndof].reshape(-1)
    mate = materials[T6.material[e]]

    Ke = T6_Ke(Xe, mate)
    Me = T6_Me(Xe, mate)
    Fe = T6_Fe(Xe, bf)

    Le0 = Ge >= 0
    Ie = Ge[Le0]  # 0-based unknown indices
    nfree = Ie.size
    if nfree:
        # Build (row,col,val) triplets for this element block
        Ke00 = Ke[Le0][:,Le0]  # Stiffness block for free dofs
        Me00 = Me[Le0][:,Le0]  # Stiffness block for free dofs
        rows = np.repeat(Ie, nfree)
        cols = np.tile(Ie, nfree)
        dataK = Ke00.ravel()
        dataM = Me00.ravel()
        rows_all.append(rows)
        cols_all.append(cols)
        dataK_all.append(dataK)
        dataM_all.append(dataM)

        F[Ie] += Fe[Le0]

rows_all = np.concatenate(rows_all) 
cols_all = np.concatenate(cols_all) 
dataK_all = np.concatenate(dataK_all) 
dataM_all = np.concatenate(dataM_all) 

K = coo_matrix((dataK_all, (rows_all, cols_all)), shape=(neq, neq)).tocsc()
M = coo_matrix((dataM_all, (rows_all, cols_all)), shape=(neq, neq)).tocsc()
K.sum_duplicates()
M.sum_duplicates()
print(f"Assembly done in {time.perf_counter() - t0:.3f} s")

# --------------------------------------------------
# solution phase
# --------------------------------------------------
print("Computing maximum stable time step estimate...")
lam_max = eigsh(K, k=1, M=M, which="LM", return_eigenvectors=False, tol=1e-6)[0]
Dt_stable = 2.0 / np.sqrt(lam_max)
print(f"Dt_stable = {Dt_stable}")

# --------------------------------------------------
# Time marching
# --------------------------------------------------
print("Time marching solution")
nstep = int(np.floor(tf / Dt))
Dstep = int(np.floor(output_interval / Dt))
Dstep = max(Dstep, 1)

# initial conditions
U = np.zeros(neq, dtype=np.float64)
Ud = np.zeros(neq, dtype=np.float64)
for n in range(analysis.NN):
    dof = int(nodes.dof[n, 1])
    if dof >= 0:
        Ud[dof] = v0

Um1 = U - Dt * Ud

# --------------------------------------------------
# initial output (placeholder)
# --------------------------------------------------
post_flag = 0
dynamics_postgmsh.postgmsh(gmshfile, analysis, nodes, post_flag, Dt=Dt, step=0)

print("Factoring mass matrix...")
factorM = cholesky(M)

percold = -1
for step in range(1, nstep + 1):
    perc = 10 * ((10 * step) // nstep) if nstep > 0 else 100
    if perc != percold:
        percold = perc
        probe_val = U[nodes.dof[3-1, 2-1]]  # node 3 in gmsh!
        print(f"Time history: {percold}%  Displ: {probe_val}")
        
 # internal forces update to account for rotation (SV-K model)
    Fg = F.copy()
    for e in range(nelem):   # loop over elements
        mate = materials[T6.material[e]]
        conn = T6.nodes[e, :]
        Ge = nodes.dof[conn, :ndof].reshape(-1)
        Xe = nodes.coor[conn, :DG]
        Ue = U[Ge]
#        Finte = T6_FeNL(Xe, mate, Ue)
        Finte = T6_FeNL_numba(Xe, mate['E'], mate['nu'], Ue)  # super fast
        Fg[Ge] += Finte

    # adds contact assuming no dbc constraints
    ne=3
    for e in range(analysis.NB3):
        conn = B3.nodes[e, :]
        Xe = nodes.coor[conn, :DG]
        Ge = nodes.dof[conn, :ndof].reshape(-1)
        Fe = B3_Fe_contact_friction(Xe, U[Ge].reshape(ne, ndof), Lambda, mu)        
        Fg[Ge] += Fe
    
    # performs the transition U_n -> U_{n+1}
    Up1 = (Dt ** 2) * factorM(Fg) + 2.0 * U - Um1
    Um1 = U
    U = Up1

    if step % Dstep == 0:
        nodes.U[:] = 0.0
        dof = nodes.dof
        mask = dof >= 0
        nodes.U[mask] = U[dof[mask]]

        post_flag = 1
        dynamics_postgmsh.postgmsh(gmshfile, analysis, nodes, post_flag, Dt=Dt, step=step)

# --------------------------------------------------
# cleanup / final post (placeholder)
# --------------------------------------------------
del K, M
print("Launching GMSH")
print("...........................")
post_flag = 2
dynamics_postgmsh.postgmsh(gmshfile, analysis, nodes, post_flag, Dt=Dt, step=step)
