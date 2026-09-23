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
    
    #computes the stiffness of a six-node triangular element
    
    E, nu = mate['E'], mate['nu'] #get material properties

    # matrix of elastic constants for plane strain - strains into stresses
    A = (E / ((1 + nu) * (1 - 2 * nu))) * np.array([
        [1 - nu, nu,      0.0],
        [nu,     1 - nu,  0.0],
        [0.0,    0.0, (1 - 2 * nu) / 2]
    ], dtype=float)

    # Gauss abscissae and weights
    a_gauss = (1/6) * np.array([ #permutations of 2/3, 1/6, 1/6
        [4, 1, 1],
        [1, 4, 1],
        [1, 1, 4]
    ], dtype=float) # area coordinates (a1,a2,1−a1−a2) for three points
    
    w_gauss = np.array([1/6, 1/6, 1/6], dtype=float) # sum to 1/2 (parent triangle area)
    
    # for an affine T6 with constant material, B is linear in position
    # B^TAB is quadratic; the three-point rule is exact
    # for curved geometry, it is approximate

    # initialize local stiffness matrix - 12 displ components on T6
    Ke = np.zeros((12, 12), dtype=float)
    
    #1. evaluate interpolation quantities at Gauss points
    #2. account for the physical geometry
    #3. accumulate the integral

    for g in range(3):
        a = a_gauss[g, :]  # a = [a1, a2, a3] at current Gauss point
    
        #D[i,:] = [dNi/da1, dNi/da2] is a 6 by 2 matrix with derivatives in reference triangle
        #N1 = a1(2a1-1)
        #N2 = a2(2a2-1)
        #N3 = (1-a1-a2)(1-2a1-2a2)
        #N4 = 4a1a2
        #N5 = 4a2(1-a1-a2)
        #N6 = 4a1(1-a1-a2)
        
        #the shape functions tell us how each node influences a point in the element
        #their derivatives tell us how nodal displacements produce deformation at that point
        
        D = np.array([
            [4*a[0] - 1, 0.0],
            [0.0, 4*a[1] - 1],
            [-4*a[2] + 1, -4*a[2] + 1],
            [4*a[1], 4*a[0]],
            [-4*a[1],4*(a[2] - a[1])],
            [4*(a[2] - a[0]), -4*a[0]],
        ], dtype=float)
        #note this does not depend on physical coordinates X
        #it depends only on the shape functions and the current Gauss point
                
        # F describes how reference coordinates map to physical coordinates 
        # (dx=F*da -> F_ij = dx_i/da_j)
        # local stretching, rotation, and distortion of the reference geometry into the physical triangle
        
        F = X.T @ D #geometry Jacobian matrix, with rows [dx/da1, dx/da2] and [dy/da1, dy/da2]
          
        #determinant of the Jacobian matrix converts the area dΩ=∣detF∣da1​da2​
        J = F[0, 0]*F[1, 1] - F[0, 1]*F[1, 0]
        
        invF = (1.0 / J) * np.array([ #F^-1
            [ F[1, 1], -F[0, 1]],
            [-F[1, 0],  F[0, 0]]
        ], dtype=float)
        
        #derivatives of the shape functions wrt the physical coordinates
        G = D @ invF #G[i,:] = [dNi/dx, dNi/dy]
            
        B = np.zeros((3, 12), dtype=float) #Allocate B 
        
        # B matrix converts 12 nodal displacements into 3 strains - (3,12)
        # Bi = [[Ni,x, 0],     #ϵxx​=∑​Ni,x * ​ui
        #       [ 0, Ni,y],    #ϵyy​=∑​Ni,y * ​vi
        #       [Ni,y, Ni,x]]  #γxy​=∑​(Ni,y * ​ui ​+ Ni,x​ * vi​)
        
        # B = [B1 B2 B3 B4 B5 B6]  --> ϵ = Bue
        
        for i in range(6):
            B[0, 2*i]     = G[i, 0] #ϵxx - Ni,x
            B[1, 2*i + 1] = G[i, 1] #ϵyy - Ni,y
            B[2, 2*i]     = G[i, 1] #γxy - Ni,y
            B[2, 2*i + 1] = G[i, 0] #γxy - Ni,x
            
        # add contributions from Gauss points
        Ke += (B.T @ A @ B) * (J * w_gauss[g]) # δueT​[∫​​BT * A * B * dΩ]ue - depends on strains
        
    return Ke



# def T6_Me(X, mate):
#     """
#     Elemental consistent mass matrix for a T6 triangle.

#     Parameters
#     ----------
#     X : (6, 2) array_like
#         Element nodal coordinates.
#     mate : sequence
#         Material properties. The third entry is the density rho.

#     Returns
#     -------
#     Me : (12, 12) ndarray
#         Element mass matrix.
#     """
#     rho = mate['rho']

#     # Six-point Gauss rule, matching the MATLAB routine.
#     a1 = 0.445948490915965
#     a2 = 0.091576213509771
#     a_gauss = np.array([
#         [a1, a1, 1.0 - 2.0 * a1],
#         [a1, 1.0 - 2.0 * a1, a1],
#         [1.0 - 2.0 * a1, a1, a1],
#         [a2, a2, 1.0 - 2.0 * a2],
#         [a2, 1.0 - 2.0 * a2, a2],
#         [1.0 - 2.0 * a2, a2, a2],
#     ], dtype=float)
#     w1 = 0.111690794839005
#     w2 = 0.054975871827661
#     w_gauss = np.array([w1, w1, w1, w2, w2, w2], dtype=float)

#     Me = np.zeros((12, 12), dtype=float)

#     for g in range(len(w_gauss)):
#         a = a_gauss[g, :]

#         D = np.array([
#             [4.0 * a[0] - 1.0, 0.0],
#             [0.0, 4.0 * a[1] - 1.0],
#             [-4.0 * a[2] + 1.0, -4.0 * a[2] + 1.0],
#             [4.0 * a[1], 4.0 * a[0]],
#             [-4.0 * a[1], 4.0 * (a[2] - a[1])],
#             [4.0 * (a[2] - a[0]), -4.0 * a[0]],
#         ], dtype=float)

#         F = X.T @ D
#         J = F[0, 0] * F[1, 1] - F[0, 1] * F[1, 0]

#         NL = np.array([
#             a[0] * (2.0 * a[0] - 1.0),
#             a[1] * (2.0 * a[1] - 1.0),
#             a[2] * (2.0 * a[2] - 1.0),
#             4.0 * a[0] * a[1],
#             4.0 * a[1] * a[2],
#             4.0 * a[0] * a[2],
#         ], dtype=float)

#         N = np.zeros((2, 12), dtype=float)
#         for i in range(6):
#             N[0, 2 * i] = NL[i]
#             N[1, 2 * i + 1] = NL[i]

#         Me += rho * (N.T @ N) * J * w_gauss[g]

#     return Me



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
    bf = np.asarray(bf, dtype=float).reshape(2) #different input required!
    # bf is a force density
    
    # initialize equivalent nodal loads - 12 displ components on T6
    Fe = np.zeros(12)

    for g in range(3):
        
        # repeat the geometry map
        
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
        
        ###### up to here identical to T6_Ke

        # Shape functions - use the quadratic Ni values to form N
        NL = np.array([
            a[0] * (2 * a[0] - 1),
            a[1] * (2 * a[1] - 1),
            a[2] * (2 * a[2] - 1),
            4 * a[0] * a[1],
            4 * a[1] * a[2],
            4 * a[0] * a[2]
        ]) # list of T6 shape functions evaluated at Gauss point

        N = np.array([
            [NL[0], 0,     NL[1], 0,     NL[2], 0,     NL[3], 0,     NL[4], 0,     NL[5], 0],
            [0,     NL[0], 0,     NL[1], 0,     NL[2], 0,     NL[3], 0,     NL[4], 0,     NL[5]]
        ])

        #δuT​[∫​​NT * b dΩ] - depends on displacements
        Fe += N.T @ bf * J * w_gauss[g] #add contributions from Gauss points
        
        #For a straight T6 under constant b, each corner-node force is zero and each
        #midside-node force is A*b/3.

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
    
    #gathers the full 12-component element displacement
    #repeats A, D, F, J, G and B as above
    #computes Sg[g,:] = A @ (B @ Ue) at each of three Gauss points
    #Sg has shape (3,3): three Gauss points by three in-plane stress components

    E, nu = mate['E'], mate['nu'] #get material properties
    Ue = Ue.reshape(-1) #already okay
    
    
    ######################### computation of A, D, F, J, G and B identical to T6_Ke

    # constitutive matrix plane strain
    A = (E / ((1 + nu) * (1 - 2 * nu))) * np.array([
        [1 - nu, nu,      0.0],
        [nu,     1 - nu,  0.0],
        [0.0,    0.0, (1 - 2 * nu) / 2]
    ], dtype=float)

    # Gauss abscissae
    a_gauss = (1/6) * np.array([
        [4, 1, 1],
        [1, 4, 1],
        [1, 1, 4]
    ], dtype=float)

    Sg = np.zeros((3, 3), dtype=float) #initialize stress field at 3 Gauss points

    for g in range(3):
        a = a_gauss[g, :]  # a = [a1, a2, a3] at current Gauss point

        #D[i,:] = [dNi/da1, dNi/da2] is a 6 by 2 matrix with derivatives in reference triangle
        D = np.array([
            [4*a[0] - 1,                 0.0],
            [0.0,                 4*a[1] - 1],
            [-4*a[2] + 1,         -4*a[2] + 1],
            [4*a[1],              4*a[0]],
            [-4*a[1],         4*(a[2] - a[1])],
            [4*(a[2] - a[0]),     -4*a[0]],
        ], dtype=float)

        # F = X' * D  -> (2,6) @ (6,2) = (2,2)
        F = X.T @ D #geometry Jacobian matrix, with rows [dx/da1, dx/da2] and [dy/da1, dy/da2]

        # J = det(F)
        J = F[0, 0]*F[1, 1] - F[0, 1]*F[1, 0] #dΩ=∣detF∣da1​da2​

        invF = (1.0 / J) * np.array([
            [ F[1, 1], -F[0, 1]],
            [-F[1, 0],  F[0, 0]]
        ], dtype=float) #F^-1

        # G = D * invF -> (6,2)
        G = D @ invF #G[i,:] = [dNi/dx, dNi/dy]

        # B matrix converts 12 nodal displacements into 3 strains - (3,12)
        # Bi = [[Ni,x, 0],     #ϵxx​=∑​Ni,x * ​ui
        #       [ 0, Ni,y],    #ϵyy​=∑​Ni,y * ​vi
        #       [Ni,y, Ni,x]]  #γxy​=∑​(Ni,y * ​ui ​+ Ni,x​ * vi​)
        
        # B = [B1 B2 B3 B4 B5 B6]
        
        B = np.zeros((3, 12), dtype=float)
        for i in range(6):
            B[0, 2*i]     = G[i, 0] #ϵxx - Ni,x
            B[1, 2*i + 1] = G[i, 1] #ϵyy - Ni,y
            B[2, 2*i]     = G[i, 1] #γxy - Ni,y
            B[2, 2*i + 1] = G[i, 0] #γxy - Ni,x

        ######################### compute stress at Gauss point: sigma = A * B * Ue
        
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

    # qG[i,:] = [sigma_xi, sigma_yi, tau_xyi] # at Gauss points
    
    # assumes a linear variation in barycentric coordinates (T6)
    # and evaluates stress field at the six nodes
    
    # linear extrapolation is applied independently to each component

    # Value of 'a' at nodes (3 x 6) - barycentric coordinates of element nodes
    aN = np.array([
        [1.0, 0.0, 0.0, 0.5, 0.0, 0.5],
        [0.0, 1.0, 0.0, 0.5, 0.5, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.5, 0.5]
    ])

    # for the Gauss points we impose qG = H * c
    
    # H collects the barycentric coordinates at the Gauss points (3,3)
    # c are the interpolation (unknown) coefficients (3,1)
    
    ########## c= H^-1 * qG ########## 
    
    # call H^-1 as N:
        
    # Mg coefficients (3 x 3)
    N = (1.0 / 3.0) * np.array([
        [5.0, -1.0, -1.0],
        [-1.0, 5.0, -1.0],
        [-1.0, -1.0, 5.0]
    ]) #this is H^-1
    
    # Extrapolation at nodes:
    # qN = aN^T * c 
    #    = aN^T * H^-1 * qG
    #    = aN^T * N * qG
    
    #since N is symmetric:
    #    (N * aN)​^T = aN^T * ​N^T = aN^T * ​N

    # Mg evaluated at nodes (3 x 6)
    fN = N @ aN

    # Linear extrapolation to nodes (6 x m)
    qN = fN.T @ qG

    return qN


