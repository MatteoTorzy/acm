

# COPY AT THE END OF ELASTICITY

# --------------------------------------------------
# Reaction force from imposed displacement, energy method
# --------------------------------------------------
# delta is the imposed displacement whose conjugate resultant is required.
# Use the same value imposed in the input file, for example -0.1.
delta = -0.1
work = 0.0
ne = 6
for e in range(analysis.NT6):
    conn = T6.nodes[e, :]
    Xe = nodes.coor[conn, :]
    Ue = nodes.U[conn, :].reshape(-1)
    mate = materials[T6.material[e]]

    Ke = T6_Ke(Xe, mate)
    work += Ue @ Ke @ Ue

reaction_force_energy = work / delta
print("Reaction force, energy method:", reaction_force_energy)

# Exact Cartesian solution for the compression benchmark without friction.
# Change L and H to match the specimen dimensions used in the mesh.
mate =  materials[T6.material[0]]  # ges material of first element (assume all equal)
E, nu = mate['E'], mate['nu']
L = 10.0
H = 5.0
epsyy = delta / L
sigmayy = E / (1.0 - nu**2) * epsyy
reaction_force_exact = sigmayy * H
print("Reaction force, exact no-friction solution:", reaction_force_exact)
