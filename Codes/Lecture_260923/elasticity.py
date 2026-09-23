# preamble - imports
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


# --------------------------------------------------
# imports input file: replace "elasticity_prova" with the desired input file located in input directory
# --------------------------------------------------
#from input.elasticity_plate_hole import gmshfile, materials, solid, dbc, dbcn, tbc
#from input.elasticity_plate_ellipse import gmshfile, materials, solid, dbc, dbcn, tbc
#from input.elasticity_plate_compress import gmshfile, materials, solid, dbc, dbcn, tbc
from input.example import gmshfile, materials, solid, dbc, dbcn, tbc

#gmshfile --> supplied mesh file
#materials, solid, dbc, dbcn, tbc --> material and boundary-condition definitions
#you can open input.py 

gmshfile="input/" + gmshfile #create path for location of gmshfile

import readgmsh

#Mesh parsing and material/BC assignment: input dictionaries become mesh, material and BC arrays
analysis,T6,B3,nodes = readgmsh.readgmsh(gmshfile, materials, solid, dbc, dbcn, tbc) #pass all input variables to readgmsh

#data structures required for the finite element analysis:
#analysis: N_nodes, N_T6, N_B3 (loaded)
#T6: connectivity, material (0-based) for each T6
#B3: connectivity, xyBC, press_BC for active B3
#nodes: U (dof values), nodal coor, dof (numbering)

# print(analysis)
# print(T6.nodes[:2])
# print(T6.material[:2])
# print(np.unique(nodes.dof,return_counts=True))

ndof = 2   # number of dof per node
Sdim = 3   # number of output stress comp

# --------------------------------------------------
# Sets global numbering of unknowns: each active free displacement receives an equation number
# --------------------------------------------------
ne = 6     # number of nodes per T6
neq = 0
    
#progressive dof numbering (nodes.dof)
for e in range(analysis.NT6): #loop over T6 elements
    for k in range(ne): #loop over nodes per element
        n = T6.nodes[e, k] #get node number from connectivity
        for d in range(ndof): #loop over dofs per per node
            # initially nodes.dof = -1 everywhere
            if nodes.dof[n, d] == -1: #skip marked constrained dofs (-2) and already numbered dofs
                nodes.dof[n, d] = neq  # dofs are 0-based numbered
                neq += 1 #progressive numbering
#all -1 entries should have been removed at the end of numbering

print("Number of unknowns:", neq)

#12,796 total displacement components
#10 quadratic segments on either vertical side contain 21 consrained node each
#their x-constrained components gives 42 conditions
#fixing the y component of a mesh node1 gives one condition 
#the system has 12,796 - 43 = 12,753 unknowns

# --------------------------------------------------
# Assemble stiffness matrix
# --------------------------------------------------
print("Assembling stiffness matrix...")

F = np.zeros(neq) #initialize rhs
from T6elemental import T6_Ke

rows_all = [] #empty list of global row indices for stiffness contributions
cols_all = [] #empty list of paired global row indices for stiffness contributions
data_all = [] #empty list of local stiffness contributions

#loop over T6 elements
for e in range(analysis.NT6):
    conn = T6.nodes[e, :] # 0-based connectivity (6,)
    Xe = nodes.coor[conn, :] # (6, 2) #element nodal coordinates
    Ge = nodes.dof[conn, :].reshape(-1) #flatten array # (12,) of element global dofs
    UDe = nodes.U[conn, :].reshape(-1) #flatten array # (12,) of element dofs values
    mate = materials[T6.material[e]] #element material id

    Ke = T6_Ke(Xe, mate)  # element stiffness matrix
    
    #Ge is used to separate free and prescribed degrees of freedom
    #contributions associated with unknowns are added to the global stiffness matrix
    #the effect of prescribed displacements is transferred to rhs
    
    #only free dofs enter in the stiffness assembly
    Le0 = Ge >= 0     # Boolean mask for negative values (constrianed dofs)
    Ie = Ge[Le0]      # 0-based unknown (reduced) indices
    
    # stiffness block for free dofs from a two-step slicing
    Ke00 = Ke[Le0,:][:,Le0]  #extract square sub-matrix selecting only Le0 rows and colums
    nfree = Ie.size
    if nfree: #expand each local block into coordinate triplets, e.g., for Ie=[2,7], the triplets have index pairs (2,2),(2,7),(7,2),(7,7)
        rows = np.repeat(Ie, nfree) #repeat each dof in Ie a given number of times
        cols = np.tile(Ie, nfree) #repeat Ie a given number of times 
        data = Ke00.reshape(-1) #Flattening of Ke00 (row-major flattening)
        #consistent with the orders of rows and cols
        
        rows_all.append(rows)
        cols_all.append(cols)
        data_all.append(data) #different elements sharing a node contribute shared global entries!


    LeD = Ge < 0    # rhs contribution (displacement induced) - selects excluded positions
    
    #Kff*Uf + Kfd*Ud = ff --> Kff*Uf = ff - Kfd*Ud (account imposed displacements into other dofs)
    
    if LeD.size > 0 and nfree: #if only constrained dofs, this is not activated
       F[Ie] -= Ke[Le0,:][:,LeD] @ UDe[LeD] # - Kfd*Ud 
       
    #######################################
    
    # #optional for body forces:
    # from T6elemental import T6_Fe
    # bf = np.array([0.0, -1.0]) # example force density
    # Fe_body = T6_Fe(Xe, bf)
    # F[Ie] += Fe_body[Le0]
    # #body forces should then be included also in the external-force vector for reaction recovery
    
    #######################################
    
# from lists to arrays:
rows_all = np.concatenate(rows_all) 
cols_all = np.concatenate(cols_all) 
data_all = np.concatenate(data_all) 


#creates global K from coordinate trplets
K = coo_matrix((data_all, (rows_all, cols_all)), shape=(neq, neq)).tocsc() 
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html 

# this format creates the global K without allocating a huge matrix full of zeros
# the .tocsc() method convert the ‘ijv’ or ‘triplet’ format into CSC (compressed sparse column), 
# sparse format reorganize the data in memory to allow for much faster computations. 

# the COO format is fast to build matrices by appending new entries 
# it allows no reordering and duplicated coordinates. 
# however, operations would require every time to scan the entire (unordered) list
# which is absouletely inefficient

# In CSC elements are ordered column by column in memory to speed up linear algebra
# it is a compression tool that stores starting and ending columns positions

K.sum_duplicates()  # important if many entries repeat (they do in FEM) - standard FEM assembly
# elements sharing the same node give rise to different triplets for the same dof 
# after .tocsc() duplicates are positioned side by side
# .sum_duplicates() adds duplcates in a single cell

# --------------------------------------------------
# Sparsity visualization: display original and reordered sparsity patterns
# the local element matrix is small and dense
# the assembled matrix is large and sparse
# this because local shape functions overlaps
# global shape functions have compact support
# nodes not sharing an element have no direct stiffness coupling
# therefore global entry can be nonzero only if its two DOFs share an element
# --------------------------------------------------
if False :

    # print(f'Potential matrix entries: {neq**2}')
    # print('Contributions before merging:', len(data_all))
    # print('Stored entries in CSC:', K.size)
    # print(f'Stored-entry fraction: {K.size / neq**2 * 100:.2f}%')
    
    # dense_memory = neq * neq * K.dtype.itemsize
    # sparse_memory = K.data.nbytes + K.indices.nbytes + K.indptr.nbytes
    # saving = (1 - (sparse_memory / dense_memory)) * 100
    # print(f'Memory dense matrix: {dense_memory / 1024**2:.2f} MB')
    # print(f'Memory sparse matrix: {dense_memory / 1024**2:.2f} MB')
    # print(f'RAM saving:       {saving:.2f}%')
    
    from scipy.sparse.csgraph import reverse_cuthill_mckee
    import matplotlib.pyplot as plt
   
    #reverse Cuthill-McKee renumbering changes the ordering to optimize the skyline
    #it does not change the number of stored entries
    
    #Kp is only plotted: the solve still uses K
    #sparse solvers may apply their own ordering
    
    #solving with Kp requires permuting RHS as F[perm] and restore U
    
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
# Assemble surface loads (B3): boundary force contributions are added to the RHS
# --------------------------------------------------
print("Assembling surface loads...")

from B3elemental import B3_Fe_surf

#loop over B3 elements counted as "active" 
for e in range(analysis.NB3):

    conn = B3.nodes[e, :] # 0-based connectivity (3,)
    Xe = nodes.coor[conn, :] # (3, 2) #element nodal coordinates
    Ge = nodes.dof[conn, :].reshape(-1) #flatten array # (6,) of element global dofs
    loadxy = B3.loadxy[e, :] #xyBC traction components
    press = B3.press[e] #press_BC

    Fe = B3_Fe_surf(Xe, loadxy, press)  # elemental nodal force array

    #Only the contributions associated with free dofs are added to the rhs
    Le0 = Ge >= 0   # Boolean mask for negative values (constrianed dofs)
    F[Ge[Le0]] += Fe[Le0] # prescribed components are absent from the rhs
    
# --------------------------------------------------
# Solve system
# --------------------------------------------------
print("Solving system...")
U = spsolve(K, F)
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve.html
# sparse solve uses factorization and ordering without forming an inverse

# --------------------------------------------------
# Scatter solution back to nodes
# --------------------------------------------------

# U contains only free displacements.
# nodes.U contains the complete nodal field
# create a mask to update free dofs without changing prescribed ones
mask = nodes.dof >= 0
nodes.U[mask] = U[nodes.dof[mask]]

# --------------------------------------------------
# Stress evaluation first at Gauss nodes and then exptrapolates to nodes for gmsh rendering
# --------------------------------------------------
from T6elemental import T6_Sg, T6_g2n

Sn = np.zeros((analysis.NN,Sdim),float) #initialize stress vector at all nodes (NN,3)
counter = np.zeros(analysis.NN,int) #counter for elems contributing to each nodal stress

#loop over T6 elements
for e in range(analysis.NT6):
    conn = T6.nodes[e, :]                    # 0-based connectivity (6,)
    Xe = nodes.coor[conn, :]                 # (6, 2) #element nodal coordinates
    Ue = nodes.U[conn, :].reshape(-1)        # (12,) - flattened element dof values
    mate = materials[T6.material[e]]         # material properties
    Sg = T6_Sg(Xe, mate, Ue)                 # (3, Sdim) stresses at Gauss points
    Sne = T6_g2n(Sg)                         # (6, Sdim) extrapolated to nodes
    Sn[conn, :] += Sne                       # accumulate contributions from T6
    # Sn[conn, :] = Sn[conn, :]+  Sne 
    counter[conn] += 1                       # increase counter for current T6 nodes

# arithmetic average
Sn = Sn / counter.reshape(-1, 1) #last dimension must be 1
#smooths nodal field without ensuring local equilibrium

# --------------------------------------------------
# Post processing
# --------------------------------------------------
import elasticity_postgmsh
elasticity_postgmsh.postgmsh(gmshfile,analysis,nodes,Sn)

# What happens if scaling every traction and prescribed displacement by alpha?
# What happens if scaling E?

print("Done")

# The procedure generates three files:
# • out.msh: contains nodal results (displacements and stresses)
# • opti.geo: defines visualization options
# • post.msh: a script that merges the mesh and results


