import numpy as np


def B3_Fe_surf(X, loadxy, press):

    a_gauss = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    w_gauss = np.array([1.0, 1.0])
    Fe = np.zeros(6)

    for g in range(2):
        a = a_gauss[g]
        D = np.array([a - 0.5, a + 0.5, -2*a])  # deriv. of shape functions 
        F = D @ X                               
        J = np.sqrt(F[0]**2 + F[1]**2)          # jacobian

        TD = loadxy + press / J * np.array([F[1], -F[0]])  # total force
        
        NL = np.array([          
            0.5 * a * (a - 1),
            0.5 * a * (1 + a),
            1 - a**2
        ])  # list of B3 shape functions
        N = np.array([           
            [NL[0], 0,      NL[1], 0,      NL[2], 0],
            [0,     NL[0],  0,     NL[1],  0,     NL[2]]
        ])  # (2,6) matrix of shape functions

        Fe = Fe + N.T @ TD * J * w_gauss[g]

    return Fe


import numpy as np




def B3_Fe_contact(X, U, lambda_contact):
    """
    Nodal forces due to contact with rigid surface in y = 0.
    Parameters
    ----------
    X : (3, 2) ndarray
        Initial nodal coordinates of the B3 boundary element.
    U : (3, 2) ndarray
        Current nodal displacements of the same element.
    lambda_contact : float
        Contact penalty coefficient.

    Returns
    -------
    Fe : (6,) ndarray
        Element nodal force vector due to contact.
    """

    a_gauss = np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
    w_gauss = np.array([1.0, 1.0])
    Fe = np.zeros(6, dtype=float)
    Y = X + U  # actual coordinates

    for g in range(2):
        a = a_gauss[g]

        NL = np.array([    
            0.5 * a * (a - 1.0),
            0.5 * a * (1.0 + a),
            1.0 - a ** 2,
        ])    # vector of shape functions
        y2 = NL @ Y[:, 1]

        if y2 < 0.0:
        # derivative of shape functions
            D = np.array([a - 0.5, a + 0.5, -2.0 * a])
            F = D @ X
            J = np.sqrt(F[0] ** 2 + F[1] ** 2)
            N = np.array([
                [NL[0], 0.0,   NL[1], 0.0,   NL[2], 0.0],
                [0.0,   NL[0], 0.0,   NL[1], 0.0,   NL[2]],
                ])    # matrix of shape functions
     
            TD = lambda_contact * np.array([0.0, -y2])
            Fe = Fe + N.T @ TD * J * w_gauss[g]

    return Fe



def B3_Fe_contact_friction(X, U, lambda_contact, mu_s):
    """
    Nodal forces due to contact with a rigid line y = -x, including friction,

    Parameters
    ----------
    X : (3, 2) ndarray
        Initial nodal coordinates of the B3 boundary element.
    U : (3, 2) ndarray
        Current nodal displacements of the same element.
    lambda_contact : float
        Contact penalty coefficient.
    mu_s : float
        Friction coefficient.

    Returns
    -------
    Fe : (6,) ndarray
        Element nodal force vector due to contact and friction.
    """

    m = -1.0/2.0
#    m = 0.0
#    mu_s=0.0
    
    alpha = np.arctan(m)

    a_gauss = np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
    w_gauss = np.array([1.0, 1.0])

    Fe = np.zeros(6, dtype=float)
    Y = X + U  # actual coordinates

    for g in range(2):
        a = a_gauss[g]
        NL = np.array([
            0.5 * a * (a - 1.0),
            0.5 * a * (1.0 + a),
            1.0 - a ** 2,
        ])
        Y_point = NL @ Y
        d_crit = (Y_point[1] - m * Y_point[0]) / np.sqrt(1.0 + m ** 2)

        if d_crit < 0.0:
            D = np.array([a - 0.5, a + 0.5, -2.0 * a])
            F = D @ X
            J = np.sqrt(F[0] ** 2 + F[1] ** 2)
            N = np.array([
                [NL[0], 0.0, NL[1], 0.0, NL[2], 0.0],
                [0.0, NL[0], 0.0, NL[1], 0.0, NL[2]],
                ])

            d_crit = abs(d_crit)
            TD = d_crit * lambda_contact * np.array([
                -np.sin(alpha),
                np.cos(alpha),
            ])

            friction = -mu_s * d_crit * lambda_contact * np.array([
                np.cos(alpha),
                np.sin(alpha),
            ])

            TD = TD + friction
            Fe = Fe + N.T @ TD * J * w_gauss[g]

    return Fe
