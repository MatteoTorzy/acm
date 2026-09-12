
globals().clear()

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

# --------------------------------------------------
# imports input file: replace "modal_beamCC" with the desired input file located in input directory
# --------------------------------------------------

# --------------------------------------------------
# --------------------------------------------------
from input.modal_beamCC import gmshfile, materials, solid, dbc, dbcn, tbc
# --------------------------------------------------
# --------------------------------------------------

gmshfile="input/" + gmshfile
import readgmsh
analysis,T6,B3,nodes = readgmsh.readgmsh(gmshfile, materials, solid, dbc, dbcn, tbc)

ndof = 2   # number of dof per node
Sdim = 3   # number of output stress comp

# --------------------------------------------------
# Sets global numbering of unknowns
# --------------------------------------------------
ne = 6     # number of nodes per T6
neq = 0
for e in range(analysis.NT6):  # only nodes associated to elements are active
    for k in range(ne):
        n = T6.nodes[e, k]
        for d in range(ndof):
            if nodes.dof[n, d] == -1:
                nodes.dof[n, d] = neq  # dofs are 0-based numbered
                neq += 1
print("Number of unknowns:", neq)


# --------------------------------------------------
# Assemble stiffness matrix
# --------------------------------------------------
print("Assembling stiffness matrix...")

F = np.zeros(neq)
from T6elemental import T6_Ke, T6_Me

rows_all = []
cols_all = []
dataK_all = []
dataM_all = []

for e in range(analysis.NT6):
    conn = T6.nodes[e, :]
    Xe = nodes.coor[conn, :]
    Ge = nodes.dof[conn, :].reshape(-1)
    UDe = nodes.U[conn, :].reshape(-1)
    mate = materials[T6.material[e]]
    
    Ke = T6_Ke(Xe, mate)
    Me = T6_Me(Xe, mate)

    Le0 = Ge >= 0    
    Ie = Ge[Le0]  
    Ke00 = Ke[Le0][:,Le0]  # Stiffness block for free dofs
    Me00 = Me[Le0][:,Le0]  # Stiffness block for free dofs
    nfree = Ie.size
    if nfree:
        # Build (row,col,val) triplets for this element block
        rows = np.repeat(Ie, nfree)
        cols = np.tile(Ie, nfree)
        dataK = Ke00.ravel()
        dataM = Me00.ravel()
        rows_all.append(rows)
        cols_all.append(cols)
        dataK_all.append(dataK)
        dataM_all.append(dataM)

rows_all = np.concatenate(rows_all) 
cols_all = np.concatenate(cols_all) 
dataK_all = np.concatenate(dataK_all) 
dataM_all = np.concatenate(dataM_all) 
K = coo_matrix((dataK_all, (rows_all, cols_all)), shape=(neq, neq)).tocsc()
M = coo_matrix((dataM_all, (rows_all, cols_all)), shape=(neq, neq)).tocsc()

# --------------------------------------------------
# Compute eigenvalues
# --------------------------------------------------
print("Computing eigenvalues...")

nmodes = 10
eigvals, V = eigsh(K, k=nmodes, M=M, sigma=0.0, which="LM")

omega = np.sqrt(eigvals)
print(omega)

import modal_postgmsh
modal_postgmsh.postgmsh(gmshfile, analysis, nodes, V, ndof=ndof)
print("Done")


