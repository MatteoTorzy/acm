
globals().clear()

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


# --------------------------------------------------
# imports input file: replace "elasticity_prova" with the desired input file located in input directory
# --------------------------------------------------
from input.example import gmshfile, materials, solid, dbc, dbcn, tbc
#from input.elasticity_plate_hole import gmshfile, materials, solid, dbc, dbcn, tbc
#from input.elasticity_plate_ellipse import gmshfile, materials, solid, dbc, dbcn, tbc
#from input.elasticity_plate_compress import gmshfile, materials, solid, dbc, dbcn, tbc

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
from T6elemental import T6_Ke

rows_all = []
cols_all = []
data_all = []

for e in range(analysis.NT6):
    conn = T6.nodes[e, :]
    Xe = nodes.coor[conn, :]
    Ge = nodes.dof[conn, :].reshape(-1)
    UDe = nodes.U[conn, :].reshape(-1)
    mate = materials[T6.material[e]]

    Ke = T6_Ke(Xe, mate)  # elemental stiffness matrix

    Le0 = Ge >= 0     # mask
    Ie = Ge[Le0]      # 0-based unknown indices
    Ke00 = Ke[Le0,:][:,Le0]  # Stiffness block for free dofs
    nfree = Ie.size
    if nfree:
        rows = np.repeat(Ie, nfree)
        cols = np.tile(Ie, nfree)
#        data = Ke00.ravel()
        data = Ke00.reshape(-1)
        rows_all.append(rows)
        cols_all.append(cols)
        data_all.append(data)

    LeD = Ge < 0    # rhs contribution
    if LeD.size > 0 and nfree:
       F[Ie] -= Ke[Le0,:][:,LeD] @ UDe[LeD]  # Ke[Le0][:,LeD] is submatrix

rows_all = np.concatenate(rows_all) 
cols_all = np.concatenate(cols_all) 
data_all = np.concatenate(data_all) 
K = coo_matrix((data_all, (rows_all, cols_all)), shape=(neq, neq)).tocsc()
K.sum_duplicates()  # important if many entries repeat (they do in FEM)


# --------------------------------------------------
# Sparsity visualization
# --------------------------------------------------
if False :
    from scipy.sparse.csgraph import reverse_cuthill_mckee
    import matplotlib.pyplot as plt
   
    perm = reverse_cuthill_mckee(K)
    Kp = K[perm, :][:, perm]
    plt.figure()
    plt.spy(K, markersize=1)
    plt.title("Original")
    plt.figure()
    plt.spy(Kp, markersize=1)
    plt.title("Reordered")
    plt.show()


# --------------------------------------------------
# Assemble surface loads (B3)
# --------------------------------------------------
print("Assembling surface loads...")

from B3elemental import B3_Fe_surf

for e in range(analysis.NB3):

    conn = B3.nodes[e, :]
    Xe = nodes.coor[conn, :]
    Ge = nodes.dof[conn, :].reshape(-1)
    loadxy = B3.loadxy[e, :]
    press = B3.press[e]

    Fe = B3_Fe_surf(Xe, loadxy, press)  # elemental nodal force array

    Le0 = Ge >= 0   # creates a boolean mask
    F[Ge[Le0]] += Fe[Le0]    
    

# --------------------------------------------------
# Solve system
# --------------------------------------------------
print("Solving system...")
U = spsolve(K, F)


# --------------------------------------------------
# Scatter solution back to nodes
# --------------------------------------------------
mask = nodes.dof >= 0
nodes.U[mask] = U[nodes.dof[mask]]  # dofs are 0-based numbered


# --------------------------------------------------
# Stress evaluation first at Gauss nodes and then exptrapolates to nodes for gmsh rendering
# --------------------------------------------------
from T6elemental import T6_Sg, T6_g2n

Sn = np.zeros((analysis.NN,Sdim),float)   
counter = np.zeros(analysis.NN,int)

ne = 6
for e in range(analysis.NT6):
    conn = T6.nodes[e, :]                    # 0-based connectivity (6,)
    Xe = nodes.coor[conn, :]                 # (6, DG)
    Ue = nodes.U[conn, :].reshape(-1)        # (ne*ndof,)
    mate = materials[T6.material[e]]         # material properties
    Sg = T6_Sg(Xe, mate, Ue)                 # (3, Sdim) stresses at Gauss points
    Sne = T6_g2n(Sg)                         # (6, Sdim) extrapolated to nodes
    Sn[conn, :] += Sne                       # accumulate contributions
    counter[conn] += 1

# naive average (assumes every node was visited at least once)
Sn = Sn / counter.reshape(-1, 1)


# --------------------------------------------------
# Post processing
# --------------------------------------------------
import elasticity_postgmsh
elasticity_postgmsh.postgmsh(gmshfile,analysis,nodes,Sn)


print("Done")


