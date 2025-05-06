# Author: Jaafar Alaswad (jalaswad@bu.edu)

import numpy as np

def newton_raphson_solver(E, width, height, length, lambda_, mu,
                           coords, connect, prescribed_dofs,
                           numel, ne, ne_L, ngp_c, ngp_l,
                           n_load_steps, max_iterations, tolerance,
                           ANS_membrane, ANS_shear, ANS_curvature):
    """
    Newton-Raphson solver for incremental loading of a nonlinear 3D cantilever beam.

    Parameters
    ----------
    E, width, height, length : Material and geometric properties
    lambda_, mu : Lamé parameters
    coords : np.ndarray
        Initial nodal coordinates.
    connect : np.ndarray
        Element connectivity.
    prescribed_dofs : array-like
        Degrees of freedom with prescribed displacements.
    numel : int
        Number of elements.
    ne : int
        Number of nodes per element.
    ne_L : int
        Number of nodes per element in axial direction.
    ngp_c, ngp_l : int
        Number of Gauss points in cross-sectional and axial directions.
    n_load_steps : int
        Number of load steps.
    max_iterations : int
        Maximum Newton-Raphson iterations per load step.
    tolerance : float
        Convergence tolerance.

    Returns
    -------
    u : np.ndarray
        Final displacement vector.
    """
    ndof = coords.shape[0] * 3   # total number of nodes: 3 DOFs per node
    ndof_e = ne * 3              # 3 DOFs per element node

    # Initialize displacement vector
    u = np.zeros((ndof, 1))

    # Initialize convergence history storage
    ResNorm = np.zeros((n_load_steps, max_iterations))      # Euclidean norm
    EnerNorm = np.zeros((n_load_steps, max_iterations))     # Energy norm

    displacements_all = []  # Store solution per load step to visualize later

    # Incremental force method loop
    for inc in range(1, n_load_steps + 1):

        # Newton-Raphson loop per increment
        for itn in range(1, max_iterations + 1):
    
            # Initialize internal force, external force, stiffness
            F_int = np.zeros((ndof, 1))
            F_ext = np.zeros((ndof, 1))

            K = np.zeros((ndof, ndof))

            # External force computation (follwer load updated every iteration)
            F_ext, ke_l = compute_external_force(E, length, height, width,
                                                 n_load_steps, inc,
                                                 coords, connect, ne, ne_L, numel,
                                                 ngp_c, ndof, u)

            # Loop over all elements
            for ie in range(numel):
                nx0_e = np.zeros((ndof_e, 3))   # Nodal position of the current element in the reference configuration
                ue = np.zeros((ndof_e, 3))      # Initialize nodal displacements

                # Loop over nodes of the element
                for i in range(ne):
                    num = connect[ie, i]            # Global node number for local node i
                    nx0_e[i, :] = coords[num, :]    # Extract reference coordinates from global array
                    ue[i, :] = u[3*num:3*num+3, 0]  # Extract displacement components for global node
                ke, fe = element_3D_nonlinear(ne, ne_L, nx0_e, ue, lambda_, mu, ngp_c, ngp_l, ANS_membrane, ANS_shear, ANS_curvature)

                # Assemble ke_l for loaded last element only
                if ie == numel - 1:
                    ke += ke_l
                
                # Assembly of global stiffness matrix and internal force vector
                for i in range(ne):
                    ii = connect[ie, i]
                    for j in range(ne):
                        jj = connect[ie, j]
                        K[3*ii:3*ii+3, 3*jj:3*jj+3] += ke[3*i:3*i+3, 3*j:3*j+3]
                    F_int[3*ii:3*ii+3] += fe[3*i:3*i+3]


            # Residual
            R = F_int - F_ext
            
            # Solve linear system
            du, R_r = solve_linear_system(K, R, prescribed_dofs)

            # Update displacements
            u += du

            # Residual and Energy norms
            ResNorm[inc-1, itn-1] = np.linalg.norm(R_r)
            EnerNorm[inc-1, itn-1] = abs(R.T @ du)[0, 0]

            # Print current load step and Newton-Raphson iteration,
            # along with residual and energy norms for monitoring convergence
            print(f"Load Step {inc}, Iteration {itn}: Residual Norm = {ResNorm[inc-1, itn-1]:10.6e}, Energy Norm = {EnerNorm[inc-1, itn-1]:10.6e}")

            # Convergence check
            if EnerNorm[inc-1, itn-1] < tolerance:
                print("=====")
                break

        # Store displacements for visualization later
        displacements_all.append(u.copy())

    return u, displacements_all



def compute_external_force(E, length, height, width, n_load_steps, inc,
                                 coords, connect, ne, ne_L, numel, ngp_c, ndof, u):
    """
    Apply constant conservative surface traction (shear) on last element's surface.

    Parameters
    ----------
    (Same as original bending version)

    Returns
    -------
    F_ext : np.ndarray
        External force vector (ndof x 1).
    ke_l : np.ndarray
        Zero matrix (no load stiffness for conservative shear).
    """

    ndof_e = 3 * ne
    F_ext = np.zeros((ndof, 1))
    ke_l = np.zeros((ndof_e, ndof_e))  # Just zero here - linearization of ext. virtual work = 0

    list_e_surf = [numel - 1]  # Only last element is loaded

    for ie in list_e_surf:
        fe_ext = np.zeros((ndof_e, 1))

        # Apply force at top 4 nodes (z-direction)
        force = (1.0 / length**2) * (inc / n_load_steps)

        # Load last 4 nodes
        for k in range(4):
            m = ne - 1 - k
            fe_ext[3 * m + 2, 0] = force

        # Assemble into global vector
        for m in range(ne):
            mm = connect[ie, m]
            F_ext[3 * mm:3 * mm + 3, 0] += fe_ext[3 * m:3 * m + 3, 0]

    return F_ext, ke_l


def element_3D_nonlinear(ne, ne_L, nx0_e, ue, lambda_, mu, ngp_c, ngp_l, ANS_membrane, ANS_shear, ANS_curvature):
    """
    Compute element tangent stiffness matrix and internal force vector for Saint Venant-Kirchhoff (SVK) material,
    using ANS method to alleviate membrane, transverse shear, and curvature-thickness locking.

    Parameters
    ----------
    ne : int
        Number of nodes per element.
    ne_L : int
        Number of nodes along the beam length.
    nx0_e : np.ndarray
        Initial element nodal coordinates (ne, 3).
    ue : np.ndarray
        Element nodal displacements (ne, 3).
    lambda_ : float
        First Lamé parameter.
    mu : float
        Second Lamé parameter.
    ngp_c : int
        Number of Gauss points in cross-sectional directions.
    ngp_l : int
        Number of Gauss points along beam length.

    Returns
    -------
    ke : np.ndarray
        Element stiffness matrix (3*ne x 3*ne).
    fe : np.ndarray
        Element internal force vector (3*ne x 1).
    """

    ndof_e = 3 * ne     # Degrees of freedom per element

    # Initialize element material and geometric stiffnesses and internal force vector
    ke_m = np.zeros((ndof_e, ndof_e))
    ke_g = np.zeros((ndof_e, ndof_e))
    fe = np.zeros((ndof_e, 1))

    # Gauss points
    eta, wgp_c = gauss(ngp_c)
    zta = eta
    xi, wgp_l = gauss(ngp_l)

    # Current configuration
    qcur = nx0_e + ue

    # Loop over Gauss points
    for k1 in range(ngp_l):
        for k2 in range(ngp_c):
            for k3 in range(ngp_c):

                # Evaluate shape functions derivative at Gauss points
                _, dN = shape_functions_3D(ne_L, xi[k1], eta[k2], zta[k3])

                # Covariant basis in reference coordinates
                H = np.zeros((3, 3))
                for m in range(ne):
                    H[:, 0] += dN[:, m] * nx0_e[m, 0]
                    H[:, 1] += dN[:, m] * nx0_e[m, 1]
                    H[:, 2] += dN[:, m] * nx0_e[m, 2]

                # Compute determinant of Jacobian
                detJ = np.linalg.det(H)
                
                # Construct covariant basis vectors in reference configuration
                G_1 = H[:, 0]
                G_2 = H[:, 1]
                G_3 = H[:, 2]

                # Construct covariant metric tensor in reference configuration
                Gco = np.array([
                    [np.dot(G_1, G_1), np.dot(G_1, G_2), np.dot(G_1, G_3)],
                    [np.dot(G_2, G_1), np.dot(G_2, G_2), np.dot(G_2, G_3)],
                    [np.dot(G_3, G_1), np.dot(G_3, G_2), np.dot(G_3, G_3)]
                ])

                # Invert to get contravariant metric tensor in reference configuration
                Gcontra = np.linalg.inv(Gco)
                
                # Construct contravariant basis vectors in reference configuration - No need when F not calculated!
                # G1 = Gcontra[0, 0]*G_1 + Gcontra[0, 1]*G_2 + Gcontra[0, 2]*G_3
                # G2 = Gcontra[1, 0]*G_1 + Gcontra[1, 1]*G_2 + Gcontra[1, 2]*G_3
                # G3 = Gcontra[2, 0]*G_1 + Gcontra[2, 1]*G_2 + Gcontra[2, 2]*G_3

                # Covariant basis in current configuration
                H = np.zeros((3, 3))
                for m in range(ne):
                    H[0, :] += dN[:, m] * qcur[m, 0]
                    H[1, :] += dN[:, m] * qcur[m, 1]
                    H[2, :] += dN[:, m] * qcur[m, 2]

                # Construct covariant basis vectors in current configuration
                g_1 = H[:, 0]
                g_2 = H[:, 1]
                g_3 = H[:, 2]

                # Construct covariant metric tensor in current configuration
                gco = np.array([
                    [np.dot(g_1, g_1), np.dot(g_1, g_2), np.dot(g_1, g_3)],
                    [np.dot(g_2, g_1), np.dot(g_2, g_2), np.dot(g_2, g_3)],
                    [np.dot(g_3, g_1), np.dot(g_3, g_2), np.dot(g_3, g_3)]
                ])

                # Deformation gradient - no need!
                # F = np.outer(g_1, G1) + np.outer(g_2, G2) + np.outer(g_3, G3)

                # Green-Lagrange strain tensor
                E = 0.5 * (gco - Gco)

                # 2nd Piola-Kirchhoff stress (SVK- material)
                S = np.zeros((3, 3))
                for I in range(3):
                    for J in range(3):
                        for K in range(3):
                            for L in range(3):
                                S[I, J] += compute_C_SVK_component(Gcontra, lambda_, mu, I, J, K, L) * E[K, L]

                # Cast stress in Voigt notation
                Svgt = np.array([
                    S[0, 0],
                    S[1, 1],
                    S[2, 2],
                    S[0, 1],
                    S[0, 2],
                    S[1, 2]
                ]).reshape(6, 1)

                # Initialize standard strain-displacement matrix Be
                Be = np.zeros((6, ndof_e))

                # Calculate standard strain-displacement matrix Be
                for m in range(ne):
                    B = np.zeros((6, 3))
                    B[0, :] = g_1 * dN[0, m]
                    B[1, :] = g_2 * dN[1, m]
                    B[2, :] = g_3 * dN[2, m]
                    B[3, :] = g_2 * dN[0, m] + g_1 * dN[1, m]
                    B[4, :] = g_3 * dN[0, m] + g_1 * dN[2, m]
                    B[5, :] = g_3 * dN[1, m] + g_2 * dN[2, m]
                    Be[:, 3*m:3*m+3] = B

                # ANS - refer to accompanying document for details
                # Modify standard Be using ANS for membrane and transverse shear locking
                if ANS_membrane or ANS_shear:
                    
                    ntying1 = ne_L - 1      # Number of tying points based on reduced integration
                    xibar1 = np.array(gauss(ntying1)).flatten()     # Tying points in xi
                    Nmat1 = lagrange_shape_1D(ntying1, xi[k1])      # Evaluate 1D Lagrange shape functions for tying points

                    M1 = np.zeros((ntying1, ntying1))   # Initialize Matrix M (Equation 8)
                    for i in range(ntying1):            # Calculate Matrix M (Equation 8)
                        M1[i, :] = lagrange_shape_1D(ntying1, xibar1[i])

                    if ANS_membrane:    # Initialize modified 1st row for membrane
                        Bebar_1 = np.zeros((ntying1, 3 * ne))

                    if ANS_shear:       # Initialize modified 4th and 5th rows for transverse shear
                        Bebar_4 = np.zeros((ntying1, 3 * ne))
                        Bebar_5 = np.zeros((ntying1, 3 * ne))

                    for k in range(ntying1):
                        _, dNbar1 = shape_functions_3D(ne_L, xibar1[k], eta[k2], zta[k3])   # Evaluate shape functions at tying points

                        # Covariant basis in current configuration at tying points
                        Pbar1 = np.zeros((3, 3))
                        for m in range(ne):
                            Pbar1[0, :] += dNbar1[0:3, m] * qcur[m, 0]
                            Pbar1[1, :] += dNbar1[0:3, m] * qcur[m, 1]
                            Pbar1[2, :] += dNbar1[0:3, m] * qcur[m, 2]

                        # Construct covariant basis vectors at tying points
                        gbar_1 = Pbar1[:, 0]
                        gbar_2 = Pbar1[:, 1]
                        gbar_3 = Pbar1[:, 2]

                        # Construct the ANS modified B matrix for membrane locking
                        if ANS_membrane:
                            for m in range(ne):
                                Bbar1 = np.zeros(3)
                                Bbar1[0] = gbar_1[0] * dNbar1[0, m]
                                Bbar1[1] = gbar_1[1] * dNbar1[0, m]
                                Bbar1[2] = gbar_1[2] * dNbar1[0, m]
                                Bebar_1[k, 3 * m : 3 * m + 3] = Bbar1
                        
                        # Construct the ANS modified B matrix for transverse shear locking
                        if ANS_shear:
                            for m in range(ne):
                                Bbar4 = np.zeros(3)
                                Bbar4[0] = gbar_2[0] * dNbar1[0, m] + gbar_1[0] * dNbar1[1, m]
                                Bbar4[1] = gbar_2[1] * dNbar1[0, m] + gbar_1[1] * dNbar1[1, m]
                                Bbar4[2] = gbar_2[2] * dNbar1[0, m] + gbar_1[2] * dNbar1[1, m]
                                Bebar_4[k, 3 * m : 3 * m + 3] = Bbar4

                                Bbar5 = np.zeros(3)
                                Bbar5[0] = gbar_3[0] * dNbar1[0, m] + gbar_1[0] * dNbar1[2, m]
                                Bbar5[1] = gbar_3[1] * dNbar1[0, m] + gbar_1[1] * dNbar1[2, m]
                                Bbar5[2] = gbar_3[2] * dNbar1[0, m] + gbar_1[2] * dNbar1[2, m]
                                Bebar_5[k, 3 * m : 3 * m + 3] = Bbar5

                        # Replace 1st row for membrane - Beware Python indexing!
                        if ANS_membrane:
                            Row1 = Nmat1 @ np.linalg.solve(M1, Bebar_1)
                            Be[0, :] = Row1

                        # Replace 4th and 5th rows for transverse shear - Beware Python indexing!
                        if ANS_shear:
                            Row4 = Nmat1 @ np.linalg.solve(M1, Bebar_4)
                            Be[3, :] = Row4
                            Row5 = Nmat1 @ np.linalg.solve(M1, Bebar_5)
                            Be[4, :] = Row5
                
                # Modify standard Be using ANS for curvature-thickness locking
                if ANS_curvature:

                    ntying2 = ne_L      # Number of tying points based on equidistant nodes
                    xibar2 = np.linspace(-1, 1, ntying2)  # Use equally spaced tying points

                    Nmat2 = lagrange_shape_1D(ntying2, xi[k1])  # Tying points in xi
                    M2 = np.eye(ntying2)    # Matrix M collapses to identity - euaidistant nodes

                    # Initialize modified 2nd and 3rd rows of B matrix
                    Bebar_2 = np.zeros((ntying2, 3 * ne))
                    Bebar_3 = np.zeros((ntying2, 3 * ne))

                    for k in range(ntying2):
                        _, dNbar2 = shape_functions_3D(ne_L, xibar2[k], eta[k2], zta[k3])   # Evaluate shape functions derivatives at tying points

                        # Covariant basis in current configuration at tying points
                        Pbar2 = np.zeros((3, 3))
                        for m in range(ne):
                            Pbar2[0, :] += dNbar2[0:3, m] * qcur[m, 0]
                            Pbar2[1, :] += dNbar2[0:3, m] * qcur[m, 1]
                            Pbar2[2, :] += dNbar2[0:3, m] * qcur[m, 2]

                        # Construct covariant vectors at tying points in current configuration
                        gbar_1 = Pbar2[:, 0]
                        gbar_2 = Pbar2[:, 1]
                        gbar_3 = Pbar2[:, 2]

                        # Calculate 2nd row of B matrix
                        for m in range(ne):
                            Bbar2 = np.zeros(3)
                            Bbar2[0] = gbar_2[0] * dNbar2[1, m]
                            Bbar2[1] = gbar_2[1] * dNbar2[1, m]
                            Bbar2[2] = gbar_2[2] * dNbar2[1, m]
                            Bebar_2[k, 3 * m : 3 * m + 3] = Bbar2

                        # Calculate 3rd row of B matrix
                        for m in range(ne):
                            Bbar3 = np.zeros(3)
                            Bbar3[0] = gbar_3[0] * dNbar2[2, m]
                            Bbar3[1] = gbar_3[1] * dNbar2[2, m]
                            Bbar3[2] = gbar_3[2] * dNbar2[2, m]
                            Bebar_3[k, 3 * m : 3 * m + 3] = Bbar3

                    # Calculate and replace 2nd row of B matrix - Beware Python indexing!
                    Row2 = Nmat2 @ np.linalg.solve(M2, Bebar_2)
                    Be[1, :] = Row2

                    # Calculate and replace 3rd row of B matrix - Beware Python indexing!
                    Row3 = Nmat2 @ np.linalg.solve(M2, Bebar_3)
                    Be[2, :] = Row3


                # Internal force vector
                fe += Be.T @ Svgt * detJ * wgp_l[k1] * wgp_c[k2] * wgp_c[k3]

                # Material stiffness matrix
                matC = np.zeros((6, 6))
                for I in range(3):
                    for J in range(I, 3):
                        for K in range(3):
                            for L in range(K, 3):
                                row = voigt_index(I, J)
                                col = voigt_index(K, L)
                                matC[row, col] = compute_C_SVK_component(Gcontra, lambda_, mu, I, J, K, L)
                                if row != col:
                                    matC[col, row] = matC[row, col]

                # Material stiffness contribution
                ke_m += Be.T @ matC @ Be * detJ * wgp_l[k1] * wgp_c[k2] * wgp_c[k3]
                
                # Geometric stiffness contribution
                
                # Initialize sub-matrices
                ke_g11 = np.zeros((ndof_e, ndof_e))
                ke_g12 = np.zeros((ndof_e, ndof_e))
                ke_g13 = np.zeros((ndof_e, ndof_e))
                ke_g21 = np.zeros((ndof_e, ndof_e))
                ke_g22 = np.zeros((ndof_e, ndof_e))
                ke_g23 = np.zeros((ndof_e, ndof_e))
                ke_g31 = np.zeros((ndof_e, ndof_e))
                ke_g32 = np.zeros((ndof_e, ndof_e))
                ke_g33 = np.zeros((ndof_e, ndof_e))

                # Initialize tying arrays for membrane (Equation 19)
                if ANS_membrane:
                    tying_dNxi_dNxi = np.zeros((ntying1, 1))

                # Initialize tying arrays for transverse shear (Equation 19)
                if ANS_shear:
                    tying_dNxi_dNeta = np.zeros((ntying1, 1))
                    tying_dNeta_dNxi = np.zeros((ntying1, 1))
                    tying_dNxi_dNzta = np.zeros((ntying1, 1))
                    tying_dNzta_dNxi = np.zeros((ntying1, 1))
                
                # Initialize tying arrays for curvature-thickness (Equation 19)
                if ANS_curvature:
                    tying_dNeta_dNeta = np.zeros((ntying2, 1))
                    tying_dNzta_dNzta = np.zeros((ntying2, 1))

                # Loop over element nodes
                for K in range(ne):
                    for L in range(ne):
                        
                        # Lump constants and define identity matrix
                        factor = detJ * wgp_l[k1] * wgp_c[k2] * wgp_c[k3]
                        I3 = np.eye(3)

                        # In-plane shear contributions lways calculated the usual way (never modified by ANS)
                        ke_g23[3*K:3*K+3, 3*L:3*L+3] += dN[1, K] * S[1, 2] * I3 * dN[2, L] * factor
                        ke_g32[3*K:3*K+3, 3*L:3*L+3] += dN[2, K] * S[2, 1] * I3 * dN[1, L] * factor

                        # --- Membrane contribution to stress (Equation 14) ---
                        if ANS_membrane:    # calculated stresses in a modified way if alleviated
                            for k in range(ntying1):
                                _, dNbar1 = shape_functions_3D(ne_L, xibar1[k], eta[k2], zta[k3])
                                tying_dNxi_dNxi[k, 0] = dNbar1[0, K] * dNbar1[0, L]
                            proj_dNxi_dNxi = Nmat1 @ np.linalg.inv(M1) @ tying_dNxi_dNxi
                            Sxi_xi = S[0, 0] * proj_dNxi_dNxi
                        else:   # Calculate stresses the usual way if not alleviated
                            Sxi_xi = S[0, 0] * dN[0, K] * dN[0, L]
                        ke_g11[3*K:3*K+3, 3*L:3*L+3] += Sxi_xi * I3 * factor    # Calculate membrane geometric contribution (Equations 18 and 19)

                        # --- Transverse shear contribution to stress (Equation 14) ---
                        if ANS_shear:   # calculated stresses in a modified way if alleviated
                            for k in range(ntying1):
                                _, dNbar1 = shape_functions_3D(ne_L, xibar1[k], eta[k2], zta[k3])
                                tying_dNxi_dNeta[k, 0] = dNbar1[0, K] * dNbar1[1, L]
                                tying_dNeta_dNxi[k, 0] = dNbar1[1, K] * dNbar1[0, L]
                                tying_dNxi_dNzta[k, 0] = dNbar1[0, K] * dNbar1[2, L]
                                tying_dNzta_dNxi[k, 0] = dNbar1[2, K] * dNbar1[0, L]

                            proj_dNxi_dNeta = Nmat1 @ np.linalg.inv(M1) @ tying_dNxi_dNeta
                            proj_dNeta_dNxi = Nmat1 @ np.linalg.inv(M1) @ tying_dNeta_dNxi
                            proj_dNxi_dNzta = Nmat1 @ np.linalg.inv(M1) @ tying_dNxi_dNzta
                            proj_dNzta_dNxi = Nmat1 @ np.linalg.inv(M1) @ tying_dNzta_dNxi

                            Sxi_eta = S[0, 1] * proj_dNxi_dNeta
                            Seta_xi = S[1, 0] * proj_dNeta_dNxi
                            Sxi_zta = S[0, 2] * proj_dNxi_dNzta
                            Szta_xi = S[2, 0] * proj_dNzta_dNxi
                        else:   # Calculate stresses the usual way if not alleviated
                            Sxi_eta = S[0, 1] * dN[0, K] * dN[1, L]
                            Seta_xi = S[1, 0] * dN[1, K] * dN[0, L]
                            Sxi_zta = S[0, 2] * dN[0, K] * dN[2, L]
                            Szta_xi = S[2, 0] * dN[2, K] * dN[0, L]

                        # Calculate transvere shear geometric contributions (Equations 18 and 19)
                        ke_g12[3*K:3*K+3, 3*L:3*L+3] += Sxi_eta * I3 * factor
                        ke_g21[3*K:3*K+3, 3*L:3*L+3] += Seta_xi * I3 * factor
                        ke_g13[3*K:3*K+3, 3*L:3*L+3] += Sxi_zta * I3 * factor
                        ke_g31[3*K:3*K+3, 3*L:3*L+3] += Szta_xi * I3 * factor

                        # --- Curvature-thickness contribution to stress (Equation 14) ---
                        if ANS_curvature:
                            for k in range(ntying2):
                                _, dNbar2 = shape_functions_3D(ne_L, xibar2[k], eta[k2], zta[k3])
                                tying_dNeta_dNeta[k, 0] = dNbar2[1, K] * dNbar2[1, L]
                                tying_dNzta_dNzta[k, 0] = dNbar2[2, K] * dNbar2[2, L]

                            proj_dNeta_dNeta = Nmat2 @ np.linalg.inv(M2) @ tying_dNeta_dNeta
                            proj_dNzta_dNzta = Nmat2 @ np.linalg.inv(M2) @ tying_dNzta_dNzta

                            Seta_eta = S[1, 1] * proj_dNeta_dNeta
                            Szta_zta = S[2, 2] * proj_dNzta_dNzta
                        else:   # Calculate stresses the usual way if not alleviated
                            Seta_eta = S[1, 1] * dN[1, K] * dN[1, L]
                            Szta_zta = S[2, 2] * dN[2, K] * dN[2, L]

                        # Calculate curvature-thickness geometric contributions (Equations 18 and 19)
                        ke_g22[3*K:3*K+3, 3*L:3*L+3] += Seta_eta * I3 * factor
                        ke_g33[3*K:3*K+3, 3*L:3*L+3] += Szta_zta * I3 * factor


                # Sum all parts to form final geometric stiffness matrix (Equation 16)
                ke_g += (
                    ke_g11 + ke_g12 + ke_g13 +
                    ke_g21 + ke_g22 + ke_g23 +
                    ke_g31 + ke_g32 + ke_g33
)

    # Add material and geometric parts of element stiffness matrix
    ke = ke_m + ke_g
    return ke, fe



def gauss(ngp):
    """
    Returns Gauss-Legendre quadrature points and weights for [-1, 1].

    Parameters
    ----------
    ngp : int
        Number of Gauss points (1 to 10 supported).

    Returns
    -------
    xi : np.ndarray
        Coordinates of Gauss points.
    w : np.ndarray
        Associated weights.
    """
    if not (1 <= ngp <= 10):
        raise ValueError("Number of Gauss points must be between 1 and 10.")

    data = {
        1: ([0.0], [2.0]),
        2: ([-1/np.sqrt(3), 1/np.sqrt(3)], [1.0, 1.0]),
        3: ([-np.sqrt(3/5), 0.0, np.sqrt(3/5)], [5/9, 8/9, 5/9]),
        4: (
            [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116],
            [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451]
        ),
        5: (
            [-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459],
            [0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851]
        ),
        6: (
            [-0.9324695142, -0.6612093865, -0.2386191861,
              0.2386191861,  0.6612093865,  0.9324695142],
            [0.1713244924, 0.3607615730, 0.4679139346,
             0.4679139346, 0.3607615730, 0.1713244924]
        ),
        7: (
            [-0.9491079123, -0.7415311856, -0.4058451514,
              0.0,
              0.4058451514, 0.7415311856, 0.9491079123],
            [0.1294849662, 0.2797053915, 0.3818300505,
             0.4179591837,
             0.3818300505, 0.2797053915, 0.1294849662]
        ),
        8: (
            [-0.9602898565, -0.7966664774, -0.5255324099, -0.1834346425,
              0.1834346425,  0.5255324099,  0.7966664774,  0.9602898565],
            [0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834,
             0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363]
        ),
        9: (
            [-0.9681602395, -0.8360311073, -0.6133714327, -0.3242534234, 0.0,
              0.3242534234, 0.6133714327, 0.8360311073, 0.9681602395],
            [0.0812743884, 0.1806481607, 0.2606106964, 0.3123470770, 0.3302393550,
             0.3123470770, 0.2606106964, 0.1806481607, 0.0812743884]
        ),
        10: (
            [-0.9739065285, -0.8650633667, -0.6794095683, -0.4333953941, -0.1488743390,
              0.1488743390, 0.4333953941, 0.6794095683, 0.8650633667, 0.9739065285],
            [0.0666713443, 0.1494513492, 0.2190863625, 0.2692667193, 0.2955242247,
             0.2955242247, 0.2692667193, 0.2190863625, 0.1494513492, 0.0666713443]
        )
    }

    xi, w = data[ngp]
    return np.array(xi), np.array(w)


def solve_linear_system(K, R, prescribed_dofs):
    """
    Solves the linear system K d = -R with Dirichlet boundary conditions.

    Parameters
    ----------
    K : np.ndarray
        Global stiffness matrix (ndof x ndof).
    R : np.ndarray
        Residual force vector (ndof x 1).
    prescribed_dofs : array-like
        List or array of constrained DOFs (Dirichlet boundary conditions).

    Returns
    -------
    d : np.ndarray
        Solution displacement vector (ndof x 1).
    R_r : np.ndarray
        Reduced residual vector after removing prescribed DOFs.
    """
    ndof = K.shape[0]  # Total number of DOFs

    # All DOFs initially
    all_dofs = np.arange(ndof)

    # Free DOFs (not prescribed)
    free_dofs = np.setdiff1d(all_dofs, prescribed_dofs)

    # Initialize displacement vector
    d = np.zeros((ndof, 1))

    # Reduced system solver
    R_r = R[free_dofs]
    s = -np.linalg.solve(K[np.ix_(free_dofs, free_dofs)], R_r)      # Comment: you may want to use sparse solver for larger problems!

    # Set displacements
    d[free_dofs, 0] = s[:, 0]

    return d, R_r


def compute_C_SVK_component(Gcontra, lambda_, mu, I, J, K, L):
    """
    Computes one component of the Saint Venant–Kirchhoff material elasticity tensor.

    Parameters
    ----------
    Gcontra : np.ndarray
        Contravariant metric tensor at reference configuration (3x3 array).
    lambda_ : float
        First Lamé parameter.
    mu : float
        Second Lamé parameter (shear modulus).
    I, J, K, L : int
        Tensor indices (0-based indexing).

    Returns
    -------
    comp : float
        Component value of the reduced elasticity tensor.
    """
    comp = (lambda_ * Gcontra[I, J] * Gcontra[K, L] +
            mu * (Gcontra[I, K] * Gcontra[J, L] + Gcontra[I, L] * Gcontra[J, K]))
    return comp


def shape_functions_3D(ne_L, xi, eta, zta):
    """
    Compute shape functions and their derivatives for a 3D brick element
    with 2x2 nodes per cross-section and ne_L nodes in axial direction.

    Parameters
    ----------
    ne_L : int
        Number of nodes along the xi (length) direction.
    xi, eta, zta : float
        Natural coordinates in the reference cube [-1, 1]^3.

    Returns
    -------
    N : (4*ne_L,) array
        Shape functions evaluated at (xi, eta, zta).
    dN : (3, 4*ne_L) array
        Derivatives of shape functions with respect to xi, eta, zta.
    """
    # Nodes in each direction
    xi_nodes = np.linspace(-1, 1, ne_L)
    eta_nodes = [-1, 1]
    zta_nodes = [-1, 1]

    # 1D shape functions and derivatives in xi
    N_xi = np.array([
        np.prod([(xi - xi_nodes[m]) / (xi_nodes[k] - xi_nodes[m])
                 for m in range(ne_L) if m != k])
        for k in range(ne_L)
    ])

    dN_xi = np.zeros(ne_L)
    for i in range(ne_L):
        for j in range(ne_L):
            if j != i:
                term = 1.0 / (xi_nodes[i] - xi_nodes[j])
                for k in range(ne_L):
                    if k != i and k != j:
                        term *= (xi - xi_nodes[k]) / (xi_nodes[i] - xi_nodes[k])
                dN_xi[i] += term

    # Bilinear shape functions in (eta, zta)
    N_eta = [0.5 * (1 - eta), 0.5 * (1 + eta)]
    dN_eta = [-0.5, 0.5]

    N_zta = [0.5 * (1 - zta), 0.5 * (1 + zta)]
    dN_zta = [-0.5, 0.5]

    ne = 4 * ne_L
    N = np.zeros(ne)
    dN = np.zeros((3, ne))  # rows: d/dxi, d/deta, d/dzta

    index = 0
    for k, Nxi in enumerate(N_xi):
        for j, Nzt in enumerate(N_zta):
            for i, Net in enumerate(N_eta):
                N[index] = Net * Nzt * Nxi
                dN[0, index] = N_eta[i] * N_zta[j] * dN_xi[k]
                dN[1, index] = dN_eta[i] * N_zta[j] * Nxi
                dN[2, index] = N_eta[i] * dN_zta[j] * Nxi
                index += 1

    return N, dN


def voigt_index(i, j):
    """
    Return the Voigt notation index corresponding to the tensor component (i, j).

    Voigt notation maps 3×3 symmetric tensor indices to a 6-element vector:
        (0,0) → 0   normal xx
        (1,1) → 1   normal yy
        (2,2) → 2   normal zz
        (0,1),(1,0) → 3   shear xy
        (0,2),(2,0) → 4   shear xz
        (1,2),(2,1) → 5   shear yz
    """
    if i == j:
        return i
    elif (i, j) in [(0, 1), (1, 0)]:
        return 3
    elif (i, j) in [(0, 2), (2, 0)]:
        return 4
    elif (i, j) in [(1, 2), (2, 1)]:
        return 5


def lagrange_shape_1D(ne_L, xi):
    """
    Compute Lagrangian shape functions at a point xi for a 1D element
    with ne_L nodes along the xi direction.

    Parameters
    ----------
    ne_L : int
        Number of nodes in xi direction (order + 1).
    xi : float
        Evaluation point in the reference interval [-1, 1].

    Returns
    -------
    N_xi : np.ndarray
        Lagrangian shape function values at xi.
    """
    xi_nodes = np.linspace(-1, 1, ne_L)
    N_xi = np.zeros(ne_L)

    for k in range(ne_L):
        product = 1.0
        for m in range(ne_L):
            if m != k:
                product *= (xi_nodes[m] - xi) / (xi_nodes[m] - xi_nodes[k])
        N_xi[k] = product

    return N_xi