
from numba import njit
@njit(cache=True)
def T6_FeNL_numba(X, E, nu, Ue):
    """
    Numba-ready nonlinear internal force vector for a T6 triangle
    using Green-Lagrange strain and 2nd Piola-Kirchhoff stress.

    Parameters
    ----------
    X : (6, 2) float64 ndarray
        Element nodal coordinates.
    E : float
        Young's modulus.
    nu : float
        Poisson ratio.
    Ue : (12,) float64 ndarray
        Element displacement vector [u1,v1,u2,v2,...,u6,v6].

    Returns
    -------
    Fe : (12,) float64 ndarray
        Element internal force vector.
    """
    Fe = np.zeros(12, dtype=np.float64)


    A1 = 0.445948490915965 
    A2 = 0.091576213509771

    A_GAUSS = np.array([
        [A1, A1, 1.0 - 2.0 * A1],
        [A1, 1.0 - 2.0 * A1, A1],
        [1.0 - 2.0 * A1, A1, A1],
        [A2, A2, 1.0 - 2.0 * A2],
        [A2, 1.0 - 2.0 * A2, A2],
        [1.0 - 2.0 * A2, A2, A2],
        ], dtype=np.float64)

    W_GAUSS = np.array([
        0.111690794839005,
        0.111690794839005,
        0.111690794839005,
        0.054975871827661,
        0.054975871827661,
        0.054975871827661,
        ], dtype=np.float64)
    
    # Plane strain constitutive matrix entries
    c = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
    A00 = c * (1.0 - nu)
    A01 = c * nu
    A11 = c * (1.0 - nu)
    A22 = c * (1.0 - 2.0 * nu) / 2.0

    for ig in range(6):
        a0 = A_GAUSS[ig, 0]
        a1 = A_GAUSS[ig, 1]
        a2 = A_GAUSS[ig, 2]
        w = W_GAUSS[ig]

        # D = dN/da, shape (6,2)
        D00 = 4.0 * a0 - 1.0
        D01 = 0.0
        D10 = 0.0
        D11 = 4.0 * a1 - 1.0
        D20 = -4.0 * a2 + 1.0
        D21 = -4.0 * a2 + 1.0
        D30 = 4.0 * a1
        D31 = 4.0 * a0
        D40 = -4.0 * a1
        D41 = 4.0 * (a2 - a1)
        D50 = 4.0 * (a2 - a0)
        D51 = -4.0 * a0

        # F = X^T D, shape (2,2)
        F00 = (
            X[0, 0] * D00 + X[1, 0] * D10 + X[2, 0] * D20 +
            X[3, 0] * D30 + X[4, 0] * D40 + X[5, 0] * D50
        )
        F01 = (
            X[0, 0] * D01 + X[1, 0] * D11 + X[2, 0] * D21 +
            X[3, 0] * D31 + X[4, 0] * D41 + X[5, 0] * D51
        )
        F10 = (
            X[0, 1] * D00 + X[1, 1] * D10 + X[2, 1] * D20 +
            X[3, 1] * D30 + X[4, 1] * D40 + X[5, 1] * D50
        )
        F11 = (
            X[0, 1] * D01 + X[1, 1] * D11 + X[2, 1] * D21 +
            X[3, 1] * D31 + X[4, 1] * D41 + X[5, 1] * D51
        )

        J = F00 * F11 - F01 * F10
        invJ = 1.0 / J

        # invF
        invF00 =  F11 * invJ
        invF01 = -F01 * invJ
        invF10 = -F10 * invJ
        invF11 =  F00 * invJ

        # G = D * invF, shape (6,2)
        Gx0 = D00 * invF00 + D01 * invF10
        Gy0 = D00 * invF01 + D01 * invF11
        Gx1 = D10 * invF00 + D11 * invF10
        Gy1 = D10 * invF01 + D11 * invF11
        Gx2 = D20 * invF00 + D21 * invF10
        Gy2 = D20 * invF01 + D21 * invF11
        Gx3 = D30 * invF00 + D31 * invF10
        Gy3 = D30 * invF01 + D31 * invF11
        Gx4 = D40 * invF00 + D41 * invF10
        Gy4 = D40 * invF01 + D41 * invF11
        Gx5 = D50 * invF00 + D51 * invF10
        Gy5 = D50 * invF01 + D51 * invF11

        # gradU = GM * Ue
        # gradU = [u,x ; u,y ; v,x ; v,y]
        g0 = (
            Gx0 * Ue[0] + Gx1 * Ue[2] + Gx2 * Ue[4] +
            Gx3 * Ue[6] + Gx4 * Ue[8] + Gx5 * Ue[10]
        )
        g1 = (
            Gy0 * Ue[0] + Gy1 * Ue[2] + Gy2 * Ue[4] +
            Gy3 * Ue[6] + Gy4 * Ue[8] + Gy5 * Ue[10]
        )
        g2 = (
            Gx0 * Ue[1] + Gx1 * Ue[3] + Gx2 * Ue[5] +
            Gx3 * Ue[7] + Gx4 * Ue[9] + Gx5 * Ue[11]
        )
        g3 = (
            Gy0 * Ue[1] + Gy1 * Ue[3] + Gy2 * Ue[5] +
            Gy3 * Ue[7] + Gy4 * Ue[9] + Gy5 * Ue[11]
        )

        # Green-Lagrange strain
        GL0 = g0 + 0.5 * g0 * g0 + 0.5 * g2 * g2
        GL1 = g3 + 0.5 * g1 * g1 + 0.5 * g3 * g3
        GL2 = g1 + g2 + g0 * g1 + g2 * g3

        # S = A * GL
        S20 = A00 * GL0 + A01 * GL1
        S21 = A01 * GL0 + A11 * GL1
        S22 = A22 * GL2


        scale = J * w

        for a in range(6):
            iu = 2 * a
            iv = iu + 1

            if a == 0:
                gx = Gx0
                gy = Gy0
            elif a == 1:
                gx = Gx1
                gy = Gy1
            elif a == 2:
                gx = Gx2
                gy = Gy2
            elif a == 3:
                gx = Gx3
                gy = Gy3
            elif a == 4:
                gx = Gx4
                gy = Gy4
            else:
                gx = Gx5
                gy = Gy5

            # Bt row contributions for displacement dof u_a
            bt0_u = gx + g0 * gx
            bt1_u = g1 * gy
            bt2_u = gy + g0 * gy + g1 * gx

            # Bt row contributions for displacement dof v_a
            bt0_v = g2 * gx
            bt1_v = gy + g3 * gy
            bt2_v = gx + g2 * gy + g3 * gx

            Fe[iu] -= (bt0_u * S20 + bt1_u * S21 + bt2_u * S22) * scale
            Fe[iv] -= (bt0_v * S20 + bt1_v * S21 + bt2_v * S22) * scale

    return Fe
