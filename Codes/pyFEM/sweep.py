import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import gmsh                     
from pathlib import Path
import readgmsh
import elasticity_postgmsh
import re
import matplotlib.pyplot as plt

from input.elasticity_plate_ellipse import materials, solid, dbc, dbcn, tbc
  
R1_values = np.array([0.5, 0.25, 0.1, 0.05, 0.025, 0.01])
mesh_scale = 1. # repeat with reduced mesh to check sensitivity and assess convergence
sigma_tip = []

geo_file = Path("input/elasticity_plate_ellipse.geo")
original_geo_text = geo_file.read_text()

for R1 in R1_values:
    gmshfile = f"input/ellipse_{R1:g}.msh"
    
    geo_text = original_geo_text
    
    parameters = {
    "R1": R1,
    "lc1": 0.5 * mesh_scale,
    "lc2": 0.1 * mesh_scale,
    "lc3": 0.001 * mesh_scale,
    }

    for name, value in parameters.items():
        geo_text, replacements = re.subn(
            rf"^[ \t]*{re.escape(name)}\s*=\s*[^;]+;",
            f"{name} = {value:.16g};",
            geo_text,
            count=1,
            flags=re.MULTILINE,
        )

    geo_file.write_text(geo_text)
    
    # Regenerate the mesh using the modified geometry
    gmsh.initialize()

    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(str(geo_file))
    R2 = gmsh.parser.getNumber("R2")[0]

    gmsh.model.mesh.generate(2)
    gmsh.write(gmshfile)
    gmsh.finalize()
        

    # Arrays are created on every iteration.
    analysis,T6,B3,nodes = readgmsh.readgmsh(
        gmshfile, materials, solid, dbc, dbcn, tbc)
        
    ndof = 2   # number of dof per node
    Sdim = 3   # number of output stress comp
    
    
    # --------------------------------------------------
    # Sets global numbering of unknowns
    # --------------------------------------------------
    ne = 6     # number of nodes per T6
    neq = 0
    for e in range(analysis.NT6):  
        for k in range(ne):
            n = T6.nodes[e, k]
            for d in range(ndof):
                if nodes.dof[n, d] == -1:
                    nodes.dof[n, d] = neq  
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
        Ke00 = Ke[Le0,:][:,Le0] 
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
           F[Ie] -= Ke[Le0,:][:,LeD] @ UDe[LeD]  
    
    rows_all = np.concatenate(rows_all) 
    cols_all = np.concatenate(cols_all) 
    data_all = np.concatenate(data_all) 
    K = coo_matrix((data_all, (rows_all, cols_all)), shape=(neq, neq)).tocsc()
    K.sum_duplicates()  # important if many entries repeat (they do in FEM)
    
    
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
    nodes.U[mask] = U[nodes.dof[mask]]  
    
    
    # --------------------------------------------------
    # Stress evaluation first at Gauss nodes and then exptrapolates to nodes for gmsh rendering
    # --------------------------------------------------
    from T6elemental import T6_Sg, T6_g2n
    
    Sn = np.zeros((analysis.NN,Sdim),float)   
    counter = np.zeros(analysis.NN,int)
    
    ne = 6
    for e in range(analysis.NT6):
        conn = T6.nodes[e, :]                    
        Xe = nodes.coor[conn, :]                 
        Ue = nodes.U[conn, :].reshape(-1)        
        mate = materials[T6.material[e]]         
        Sg = T6_Sg(Xe, mate, Ue)                 
        Sne = T6_g2n(Sg)                         
        Sn[conn, :] += Sne                       
        counter[conn] += 1
    
    # naive average (assumes every node was visited at least once)
    Sn = Sn / counter.reshape(-1, 1)
    
    
    # --------------------------------------------------
    # Post processing
    # --------------------------------------------------
    #elasticity_postgmsh.postgmsh(gmshfile,analysis,nodes,Sn)
    
    # record sigma_xx at the two vertical tips (0, +/-R2).
    tips = np.isclose(nodes.coor[:, 0], 0.0, atol=1e-10) & \
           np.isclose(np.abs(nodes.coor[:, 1]), R2, atol=1e-10)
    sigma_tip.append(np.max(Sn[tips, 0]))
    print(f"Tip sigma_xx = {sigma_tip[-1]:.6g}")

# compare all cases after the loop
np.savetxt("output/stress_vs_R1.csv", np.column_stack((R1_values, sigma_tip)),
           delimiter=",", header="R1_semi_minor,sigma_xx_tip", comments="")


plt.figure()
plt.plot(R1_values, sigma_tip, "*-", c='goldenrod', mfc='black', ms=20, linewidth=2)
plt.gca().invert_xaxis()
plt.xlabel("Semi-minor axis R1 (decreasing to the right)",fontsize=14)
plt.ylabel("Recovered tip stress $\sigma_{xx}$",fontsize=14)
#plt.grid(True)
plt.title("Tip stress against hole opening",fontsize=14)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig("output/stress_vs_R1.png", dpi=180)
plt.show()
