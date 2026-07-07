import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

# --------------------------------------------------
# imports input file: replace with the desired input file located in input directory
# --------------------------------------------------
# Example usage:
# from input.elasticity_plate_hole_ellipse import gmshfile, material, solid, dbc, dbcn, tbc
# gmshfile = "input/" + gmshfile
#
# import elasticity_readgmsh
# analysis, T6, B3, nodes = elasticity_readgmsh.readgmsh(gmshfile, material, solid, dbc, dbcn, tbc)



# gmshfile = "input/" + gmshfile
# import elasticity_readgmsh
# analysis, T6, B3, nodes = elasticity_readgmsh.readgmsh(gmshfile, material, solid, dbc, dbcn, tbc)



def solve_modal(analysis, T6, nodes, material, neigs=10):
    ndof = 2
    ne = 6

    # --------------------------------------------------
    # Sets global numbering of unknowns
    # --------------------------------------------------
    neq = 0
    for e in range(analysis.NT6):
        for k in range(ne):
            n = T6.nodes[e, k]
            for d in range(ndof):
                if nodes.dof[n, d] == 0:
                    neq += 1
                    nodes.dof[n, d] = neq
    analysis.neq = neq
    print("Number of unknowns:", neq)

    # --------------------------------------------------
    # Assemble stiffness and mass matrices
    # --------------------------------------------------
    print("Assembling stiffness and mass matrices...")

    rows_k = []
    cols_k = []
    data_k = []

    rows_m = []
    cols_m = []
    data_m = []

    for e in range(analysis.NT6):
        conn = T6.nodes[e, :]
        Xe = nodes.coor[conn, :]
        Ge = nodes.dof[conn, :].reshape(-1)
        mat_id = int(T6.material[e])
        mate = material[mat_id - 1]

        Ke = T6_Ke(Xe, mate)
        Me = T6_Me(Xe, mate)

        Le0 = np.where(Ge > 0)[0]
        Ie = Ge[Le0].astype(np.int64) - 1
        nfree = Ie.size

        if nfree:
            Ke00 = Ke[np.ix_(Le0, Le0)]
            Me00 = Me[np.ix_(Le0, Le0)]

            rows = np.repeat(Ie, nfree)
            cols = np.tile(Ie, nfree)

            rows_k.append(rows)
            cols_k.append(cols)
            data_k.append(Ke00.ravel())

            rows_m.append(rows)
            cols_m.append(cols)
            data_m.append(Me00.ravel())

    if rows_k:
        rows_k = np.concatenate(rows_k)
        cols_k = np.concatenate(cols_k)
        data_k = np.concatenate(data_k)
        K = coo_matrix((data_k, (rows_k, cols_k)), shape=(neq, neq)).tocsr()

        rows_m = np.concatenate(rows_m)
        cols_m = np.concatenate(cols_m)
        data_m = np.concatenate(data_m)
        M = coo_matrix((data_m, (rows_m, cols_m)), shape=(neq, neq)).tocsr()
    else:
        K = coo_matrix((neq, neq)).tocsr()
        M = coo_matrix((neq, neq)).tocsr()

    # --------------------------------------------------
    # Solve generalized eigenproblem
    # --------------------------------------------------
    print("Eig solver")

    if neq < 2:
        raise ValueError("At least 2 free DOFs are required for eigsh.")

    nmodes = min(neigs, neq - 1)
    eigvals, V = eigsh(K, k=nmodes, M=M, which="SM")

    eigvals = np.real(eigvals)
    eigvals[eigvals < 0.0] = 0.0
    idx = np.argsort(eigvals)
    eigvals = eigvals[idx]
    V = V[:, idx]

    omega = np.sqrt(eigvals)
    freq = omega / (2.0 * np.pi)

    # --------------------------------------------------
    # Scatter modal vectors back to nodes
    # --------------------------------------------------
    modes = np.zeros((analysis.NN, ndof, nmodes), dtype=float)
    for n in range(analysis.NN):
        for d in range(ndof):
            dof = nodes.dof[n, d]
            if dof > 0:
                modes[n, d, :] = V[dof - 1, :]
            else:
                modes[n, d, :] = nodes.U[n, d]

    return K, M, eigvals, omega, freq, V, modes


if __name__ == "__main__":
    import elasticity_readgmsh
    from T6elemental import T6_Ke

    # --------------------------------------------------
    # imports input file: uncomment one case below
    # --------------------------------------------------
    # from input.elasticity_prova import gmshfile, material, solid, dbc, dbcn, tbc
    # from input.elasticity_plate_compress import gmshfile, material, solid, dbc, dbcn, tbc
    # from input.elasticity_plate_hole import gmshfile, material, solid, dbc, dbcn, tbc
    from input.elasticity_plate_hole_ellipse import gmshfile, material, solid, dbc, dbcn, tbc

    gmshfile = "input/" + gmshfile
    analysis, T6, B3, nodes = elasticity_readgmsh.readgmsh(gmshfile, material, solid, dbc, dbcn, tbc)

    K, M, eigvals, omega, freq, V, modes = solve_modal(
        analysis, T6, nodes, material, neigs=10
    )

    print("Eigenvalues:")
    print(eigvals)
    print("Circular frequencies:")
    print(omega)
    print("Frequencies [Hz]:")
    print(freq)
