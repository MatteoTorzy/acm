import numpy as np


def T6_Ke(X, mate):
    """
    Parameters
    ----------
    X : (6,2) array_like
        Nodal coordinates [ [x1,y1], ... , [x6,y6] ].
    mate : array_like length 2
        [E, nu]

    Returns
    -------
    Ke : (12,12) ndarray
        Element stiffness matrix (6 nodes * 2 dof/node).
    """
    E, nu = mate['E'], mate['nu']

    # matrix of elastic constants for plane strain
    A = (E / ((1 + nu) * (1 - 2 * nu))) * np.array([
        [1 - nu, nu,      0.0],
        [nu,     1 - nu,  0.0],
        [0.0,    0.0, (1 - 2 * nu) / 2]
    ], dtype=float)

    # Gauss abscissae and weights
    a_gauss = (1/6) * np.array([
        [4, 1, 1],
        [1, 4, 1],
        [1, 1, 4]
    ], dtype=float)
    w_gauss = np.array([1/6, 1/6, 1/6], dtype=float)

    Ke = np.zeros((12, 12), dtype=float)

    for g in range(3):
        a = a_gauss[g, :]  # a = [a1, a2, a3]

        D = np.array([
            [4*a[0] - 1, 0.0],
            [0.0, 4*a[1] - 1],
            [-4*a[2] + 1, -4*a[2] + 1],
            [4*a[1], 4*a[0]],
            [-4*a[1],4*(a[2] - a[1])],
            [4*(a[2] - a[0]), -4*a[0]],
        ], dtype=float)
        F = X.T @ D
        J = F[0, 0]*F[1, 1] - F[0, 1]*F[1, 0]
        invF = (1.0 / J) * np.array([
            [ F[1, 1], -F[0, 1]],
            [-F[1, 0],  F[0, 0]]
        ], dtype=float)
        G = D @ invF

        B = np.zeros((3, 12), dtype=float)
        for i in range(6):
            B[0, 2*i]     = G[i, 0]
            B[1, 2*i + 1] = G[i, 1]
            B[2, 2*i]     = G[i, 1]
            B[2, 2*i + 1] = G[i, 0]

        Ke += (B.T @ A @ B) * (J * w_gauss[g])

    return Ke



def T6_Me(X, mate):
    """
    Elemental consistent mass matrix for a T6 triangle.

    Parameters
    ----------
    X : (6, 2) array_like
        Element nodal coordinates.
    mate : sequence
        Material properties. The third entry is the density rho.

    Returns
    -------
    Me : (12, 12) ndarray
        Element mass matrix.
    """
    rho = mate['rho']

    # Six-point Gauss rule, matching the MATLAB routine.
    a1 = 0.445948490915965
    a2 = 0.091576213509771
    a_gauss = np.array([
        [a1, a1, 1.0 - 2.0 * a1],
        [a1, 1.0 - 2.0 * a1, a1],
        [1.0 - 2.0 * a1, a1, a1],
        [a2, a2, 1.0 - 2.0 * a2],
        [a2, 1.0 - 2.0 * a2, a2],
        [1.0 - 2.0 * a2, a2, a2],
    ], dtype=float)
    w1 = 0.111690794839005
    w2 = 0.054975871827661
    w_gauss = np.array([w1, w1, w1, w2, w2, w2], dtype=float)

    Me = np.zeros((12, 12), dtype=float)

    for g in range(len(w_gauss)):
        a = a_gauss[g, :]

        D = np.array([
            [4.0 * a[0] - 1.0, 0.0],
            [0.0, 4.0 * a[1] - 1.0],
            [-4.0 * a[2] + 1.0, -4.0 * a[2] + 1.0],
            [4.0 * a[1], 4.0 * a[0]],
            [-4.0 * a[1], 4.0 * (a[2] - a[1])],
            [4.0 * (a[2] - a[0]), -4.0 * a[0]],
        ], dtype=float)

        F = X.T @ D
        J = F[0, 0] * F[1, 1] - F[0, 1] * F[1, 0]

        NL = np.array([
            a[0] * (2.0 * a[0] - 1.0),
            a[1] * (2.0 * a[1] - 1.0),
            a[2] * (2.0 * a[2] - 1.0),
            4.0 * a[0] * a[1],
            4.0 * a[1] * a[2],
            4.0 * a[0] * a[2],
        ], dtype=float)

        N = np.zeros((2, 12), dtype=float)
        for i in range(6):
            N[0, 2 * i] = NL[i]
            N[1, 2 * i + 1] = NL[i]

        Me += rho * (N.T @ N) * J * w_gauss[g]

    return Me



def T6_Fe(X, bf):
    """
    Nodal forces due to body forces for a T6 triangular element.

    Parameters
    ----------
    X : array_like, shape (6, 2)
        Nodal coordinates of the 6-node triangle.
    bf : array_like, shape (2,) or (2, 1)
        Body force vector [bx, by].

    Returns
    -------
    Fe : ndarray, shape (12,)
        Element nodal force vector.
    """
    a_gauss = (1.0 / 6.0) * np.array([
        [4.0, 1.0, 1.0],
        [1.0, 4.0, 1.0],
        [1.0, 1.0, 4.0]
    ])
    w_gauss = np.array([1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0])

    X = np.asarray(X, dtype=float)
    bf = np.asarray(bf, dtype=float).reshape(2)
    Fe = np.zeros(12)

    for g in range(3):
        a = a_gauss[g, :]

        # Derivative of shape functions w.r.t. a1, a2
        D = np.array([
            [4 * a[0] - 1, 0],
            [0, 4 * a[1] - 1],
            [-4 * a[2] + 1, -4 * a[2] + 1],
            [4 * a[1], 4 * a[0]],
            [-4 * a[1], 4 * (a[2] - a[1])],
            [4 * (a[2] - a[0]), -4 * a[0]]
        ])

        F = X.T @ D
        J = F[0, 0] * F[1, 1] - F[0, 1] * F[1, 0]

        # Shape functions
        NL = np.array([
            a[0] * (2 * a[0] - 1),
            a[1] * (2 * a[1] - 1),
            a[2] * (2 * a[2] - 1),
            4 * a[0] * a[1],
            4 * a[1] * a[2],
            4 * a[0] * a[2]
        ])

        N = np.array([
            [NL[0], 0,     NL[1], 0,     NL[2], 0,     NL[3], 0,     NL[4], 0,     NL[5], 0],
            [0,     NL[0], 0,     NL[1], 0,     NL[2], 0,     NL[3], 0,     NL[4], 0,     NL[5]]
        ])

        Fe += N.T @ bf * J * w_gauss[g]

    return Fe


def T6_Sg(X, mate, Ue):
    """
    Stresses at Gauss points for a T6 triangle 

    Parameters
    ----------
    X : (6,2) array_like
        Nodal coordinates [[x1,y1],...,[x6,y6]]
    mate : array_like length 2
        [E, nu]
    Ue : (12,) or (12,1) array_like
        Element displacement vector [u1,v1,u2,v2,...,u6,v6]^T

    Returns
    -------
    Sg : (3,3) ndarray
        Stresses at 3 Gauss points; each row is [sigma_x, sigma_y, tau_xy].
    """
    E, nu = mate['E'], mate['nu']
    Ue = Ue.reshape(-1)

    # constitutive matrix plane strain
    A = (E / ((1 + nu) * (1 - 2 * nu))) * np.array([
        [1 - nu, nu,      0.0],
        [nu,     1 - nu,  0.0],
        [0.0,    0.0, (1 - 2 * nu) / 2]
    ], dtype=float)

    # Gauss abscissae (exactly as MATLAB)
    a_gauss = (1/6) * np.array([
        [4, 1, 1],
        [1, 4, 1],
        [1, 1, 4]
    ], dtype=float)

    Sg = np.zeros((3, 3), dtype=float)

    for g in range(3):
        a = a_gauss[g, :]  # [a1, a2, a3]

        # D is (6,2) after transpose
        D = np.array([
            [4*a[0] - 1,                 0.0],
            [0.0,                 4*a[1] - 1],
            [-4*a[2] + 1,         -4*a[2] + 1],
            [4*a[1],              4*a[0]],
            [-4*a[1],         4*(a[2] - a[1])],
            [4*(a[2] - a[0]),     -4*a[0]],
        ], dtype=float)

        # F = X' * D  -> (2,6) @ (6,2) = (2,2)
        F = X.T @ D

        # J = det(F)
        J = F[0, 0]*F[1, 1] - F[0, 1]*F[1, 0]

        # invF = inv(F) (explicit, to match MATLAB)
        invF = (1.0 / J) * np.array([
            [ F[1, 1], -F[0, 1]],
            [-F[1, 0],  F[0, 0]]
        ], dtype=float)

        # G = D * invF -> (6,2)
        G = D @ invF

        # Build B (3 x 12)
        B = np.zeros((3, 12), dtype=float)
        for i in range(6):
            B[0, 2*i]     = G[i, 0]
            B[1, 2*i + 1] = G[i, 1]
            B[2, 2*i]     = G[i, 1]
            B[2, 2*i + 1] = G[i, 0]

        # Stress at Gauss point: sigma = A * B * Ue
        Sg[g, :] = (A @ (B @ Ue))

    return Sg



def T6_g2n(qG):
    """
    Extrapolation from Gauss points to nodes for T6 element.

    Parameters
    ----------
    qG : ndarray of shape (3, m)
        Values at the 3 Gauss points.

    Returns
    -------
    qN : ndarray of shape (6, m)
        Extrapolated values at the 6 nodes.
    """

    # Value of 'a' at nodes (3 x 6)
    aN = np.array([
        [1.0, 0.0, 0.0, 0.5, 0.0, 0.5],
        [0.0, 1.0, 0.0, 0.5, 0.5, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.5, 0.5]
    ])

    # Mg coefficients (3 x 3)
    N = (1.0 / 3.0) * np.array([
        [5.0, -1.0, -1.0],
        [-1.0, 5.0, -1.0],
        [-1.0, -1.0, 5.0]
    ])

    # Mg evaluated at nodes (3 x 6)
    fN = N @ aN

    # Linear extrapolation to nodes (6 x m)
    qN = fN.T @ qG

    return qN


