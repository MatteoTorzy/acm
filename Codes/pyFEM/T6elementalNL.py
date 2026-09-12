import numpy as np


def T6_KeNL(X, mate, Ue):
    """
    T6 plane-strain geometrically nonlinear element.

    Saint-Venant--Kirchhoff material, Total Lagrangian formulation.
    
        Finte = - integral(Bnl.T @ SK) dV
        KTe   =   integral(Bnl.T @ A @ Bnl + Kg) dV

    Parameters
    ----------
    X : (6, 2) array_like
        Reference nodal coordinates.
    mate : dict-like
        Material dictionary with keys 'E' and 'nu'.
    Ue : (12,) array_like
        Element displacement vector [u1,v1,u2,v2,...,u6,v6].

    Returns
    -------
    KTe : (12, 12) ndarray
        Consistent tangent stiffness matrix.
    Finte : (12,) ndarray
        Internal force vector, with the same sign as the MATLAB code.
    """

    E = mate["E"]
    nu = mate["nu"]

    A = (E / ((1.0 + nu) * (1.0 - 2.0 * nu))) * np.array([
        [1.0 - nu, nu, 0.0],
        [nu, 1.0 - nu, 0.0],
        [0.0, 0.0, (1.0 - 2.0 * nu) / 2.0],
    ], dtype=float)

    a_gauss = (1.0 / 6.0) * np.array([
        [4.0, 1.0, 1.0],
        [1.0, 4.0, 1.0],
        [1.0, 1.0, 4.0],
    ], dtype=float)
    w_gauss = np.array([1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0], dtype=float)

    Ke = np.zeros((12, 12), dtype=float)
    Finte = np.zeros(12, dtype=float)

    for g in range(3):
        a = a_gauss[g, :]

        # Derivatives of T6 shape functions with respect to (a1, a2).
        # Shape: (6 nodes, 2 parent coordinates).
        D = np.array([
            [4.0 * a[0] - 1.0, 0.0],
            [0.0, 4.0 * a[1] - 1.0],
            [-4.0 * a[2] + 1.0, -4.0 * a[2] + 1.0],
            [4.0 * a[1], 4.0 * a[0]],
            [-4.0 * a[1], 4.0 * (a[2] - a[1])],
            [4.0 * (a[2] - a[0]), -4.0 * a[0]],
        ], dtype=float)

        Jmat = X.T @ D
        detJ = Jmat[0, 0] * Jmat[1, 1] - Jmat[0, 1] * Jmat[1, 0]
        invJ = (1.0 / detJ) * np.array([
            [Jmat[1, 1], -Jmat[0, 1]],
            [-Jmat[1, 0], Jmat[0, 0]],
        ], dtype=float)

        G = D @ invJ  # dN/dX, shape (6, 2)

        # Gradient operator giving gradU = [u_X, u_Y, v_X, v_Y].
        GM = np.zeros((4, 12), dtype=float)
        for i in range(6):
            GM[0, 2 * i] = G[i, 0]
            GM[1, 2 * i] = G[i, 1]
            GM[2, 2 * i + 1] = G[i, 0]
            GM[3, 2 * i + 1] = G[i, 1]

        gradU = GM @ Ue
        ux, uy, vx, vy = gradU

        # Green-Lagrange strain in engineering vector form [E11, E22, 2E12].
        GL = np.array([
            ux + 0.5 * ux**2 + 0.5 * vx**2,
            vy + 0.5 * uy**2 + 0.5 * vy**2,
            uy + vx + ux * uy + vx * vy,
        ], dtype=float)

        # Second Piola-Kirchhoff stress vector [S11, S22, S12].
        S = A @ GL

        # Linear and displacement-dependent strain-displacement matrices.
        B = np.array([
            GM[0, :],
            GM[3, :],
            GM[1, :] + GM[2, :],
        ], dtype=float)

        Bu = np.array([
            ux * GM[0, :] + vx * GM[2, :],
            uy * GM[1, :] + vy * GM[3, :],
            ux * GM[1, :] + vx * GM[3, :] + uy * GM[0, :] + vy * GM[2, :],
        ], dtype=float)

        B += Bu

        # Geometric stiffness contribution.
        Bg = np.array([
            S[0] * GM[0, :] + S[2] * GM[1, :],
            S[2] * GM[0, :] + S[1] * GM[1, :],
            S[0] * GM[2, :] + S[2] * GM[3, :],
            S[2] * GM[2, :] + S[1] * GM[3, :],
        ], dtype=float)

        weight = detJ * w_gauss[g]
        Finte -= B.T @ S * weight
        Ke += (B.T @ A @ B + GM.T @ Bg) * weight

    return Ke, Finte


def T6_SgNL(X, mate, Ue):
    """
    Gauss-point Green-Lagrange strain and second Piola-Kirchhoff stress.

    Returns
    -------
    Eg : (3, 3) ndarray
        Rows are [E11, E22, 2E12].
    Sg : (3, 3) ndarray
        Rows are [S11, S22, S12].
    """
    X = np.asarray(X, dtype=float)
    Ue = np.asarray(Ue, dtype=float).reshape(12)

    E = mate["E"]
    nu = mate["nu"]
    A = (E / ((1.0 + nu) * (1.0 - 2.0 * nu))) * np.array([
        [1.0 - nu, nu, 0.0],
        [nu, 1.0 - nu, 0.0],
        [0.0, 0.0, (1.0 - 2.0 * nu) / 2.0],
    ], dtype=float)

    a_gauss = (1.0 / 6.0) * np.array([
        [4.0, 1.0, 1.0],
        [1.0, 4.0, 1.0],
        [1.0, 1.0, 4.0],
    ], dtype=float)

    GLg = np.zeros((3, 3), dtype=float)
    Sg = np.zeros((3, 3), dtype=float)

    for g in range(3):
        a = a_gauss[g, :]
        D = np.array([
            [4.0 * a[0] - 1.0, 0.0],
            [0.0, 4.0 * a[1] - 1.0],
            [-4.0 * a[2] + 1.0, -4.0 * a[2] + 1.0],
            [4.0 * a[1], 4.0 * a[0]],
            [-4.0 * a[1], 4.0 * (a[2] - a[1])],
            [4.0 * (a[2] - a[0]), -4.0 * a[0]],
        ], dtype=float)
        Jmat = X.T @ D
        detJ = Jmat[0, 0] * Jmat[1, 1] - Jmat[0, 1] * Jmat[1, 0]
        invJ = (1.0 / detJ) * np.array([
            [Jmat[1, 1], -Jmat[0, 1]],
            [-Jmat[1, 0], Jmat[0, 0]],
        ], dtype=float)
        G = D @ invJ

        GM = np.zeros((4, 12), dtype=float)
        for i in range(6):
            GM[0, 2 * i] = G[i, 0]
            GM[1, 2 * i] = G[i, 1]
            GM[2, 2 * i + 1] = G[i, 0]
            GM[3, 2 * i + 1] = G[i, 1]

        ux, uy, vx, vy = GM @ Ue
        GLg[g, :] = np.array([
            ux + 0.5 * ux**2 + 0.5 * vx**2,
            vy + 0.5 * uy**2 + 0.5 * vy**2,
            uy + vx + ux * uy + vx * vy,
        ], dtype=float)
        Sg[g, :] = A @ GLg[g, :]

    return GLg, Sg



def T6_FeNL(X, mate, Ue):
    """
    Nonlinear internal force vector for a T6 triangle
    using Green-Lagrange strain and 2nd Piola-Kirchhoff stress.

    Parameters
    ----------
    X : (6, 2) array_like
        Element nodal coordinates [[x1,y1], ..., [x6,y6]].
    mate : dict
        Material properties with keys:
        - 'E'  : Young's modulus
        - 'nu' : Poisson ratio
    Ue : (12,) or (12,1) array_like
        Element displacement vector:
        [u1, v1, u2, v2, ..., u6, v6]^T

    Returns
    -------
    Fe : (12,) ndarray
        Element internal force vector.
    """

    E, nu = mate['E'], mate['nu']
    Finte = np.zeros(12, dtype=float)

    # six-point Gauss rule
    a1 = 0.445948490915965
    a2 = 0.091576213509771
    a_gauss = np.array([
        [a1, a1, 1.0 - 2.0 * a1],
        [a1, 1.0 - 2.0 * a1, a1],
        [1.0 - 2.0 * a1, a1, a1],
        [a2, a2, 1.0 - 2.0 * a2],
        [a2, 1.0 - 2.0 * a2, a2],
        [1.0 - 2.0 * a2, a2, a2]
    ], dtype=float)

    w1 = 0.111690794839005
    w2 = 0.054975871827661
    w_gauss = np.array([w1, w1, w1, w2, w2, w2], dtype=float)

    # plane strain constitutive matrix
    A = (E / ((1.0 + nu) * (1.0 - 2.0 * nu))) * np.array([
        [1.0 - nu, nu, 0.0],
        [nu, 1.0 - nu, 0.0],
        [0.0, 0.0, (1.0 - 2.0 * nu) / 2.0]
    ], dtype=float)

    for i in range(len(w_gauss)):
        a = a_gauss[i, :]

        D = np.array([
            [4.0 * a[0] - 1.0, 0.0],
            [0.0, 4.0 * a[1] - 1.0],
            [-4.0 * a[2] + 1.0, -4.0 * a[2] + 1.0],
            [4.0 * a[1], 4.0 * a[0]],
            [-4.0 * a[1], 4.0 * (a[2] - a[1])],
            [4.0 * (a[2] - a[0]), -4.0 * a[0]]
        ], dtype=float)

        F = X.T @ D
        J = F[0, 0] * F[1, 1] - F[0, 1] * F[1, 0]

        invF = (1.0 / J) * np.array([
            [ F[1, 1], -F[0, 1]],
            [-F[1, 0],  F[0, 0]]
        ], dtype=float)

        G = D @ invF

        GM = np.zeros((4, 12), dtype=float)
        for i in range(6):
            GM[0, 2 * i] = G[i, 0]
            GM[1, 2 * i] = G[i, 1]
            GM[2, 2 * i + 1] = G[i, 0]
            GM[3, 2 * i + 1] = G[i, 1]

        gradU = GM @ Ue

        GL = np.array([
            gradU[0] + 0.5 * gradU[0]**2 + 0.5 * gradU[2]**2,
            gradU[3] + 0.5 * gradU[1]**2 + 0.5 * gradU[3]**2,
            gradU[1] + gradU[2] + gradU[0] * gradU[1] + gradU[2] * gradU[3]
        ], dtype=float)

        S = A @ GL

        B = np.array([
            GM[0, :],
            GM[3, :],
            GM[1, :] + GM[2, :]
        ], dtype=float)

        Bu = np.array([
            gradU[0] * GM[0, :] + gradU[2] * GM[2, :],
            gradU[1] * GM[1, :] + gradU[3] * GM[3, :],
            gradU[0] * GM[1, :] + gradU[2] * GM[3, :] +
            gradU[1] * GM[0, :] + gradU[3] * GM[2, :]
        ], dtype=float)

        B += Bu

        Finte -= B.T @ S * J * w_gauss[i]

    return Finte

