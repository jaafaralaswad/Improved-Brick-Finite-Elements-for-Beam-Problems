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
    ndof = coords.shape[0] * 3  # 3 DOFs per node
    ndof_e = ne * 3              # 3 DOFs per element node

    # Initialize displacement vector
    u = np.zeros((ndof, 1))

    # Initialize history storage
    ResNorm = np.zeros((n_load_steps, max_iterations))
    EnerNorm = np.zeros((n_load_steps, max_iterations))

    displacements_all = []  # Store solution per load step

    # Incremental force method loop
    for inc in range(1, n_load_steps + 1):

        for itn in range(1, max_iterations + 1):
    
            # Initialize internal force, external force, stiffness
            F_int = np.zeros((ndof, 1))
            F_ext = np.zeros((ndof, 1))

            K = np.zeros((ndof, ndof))

            # External force computation (incremental)
            F_ext, ke_l = compute_external_force(E, length, height, width,
                                                 n_load_steps, inc,
                                                 coords, connect, ne, ne_L, numel,
                                                 ngp_c, ndof, u)

            # Element loop
            for ie in range(numel):
                nx0_e = np.zeros((ndof_e, 3))
                ue = np.zeros((ndof_e, 3))

                for i in range(ne):
                    num = connect[ie, i]
                    nx0_e[i, :] = coords[num, :]
                    ue[i, :] = u[3*num:3*num+3, 0]
                ke, fe = element_3D_nonlinear(ne, ne_L, nx0_e, ue, lambda_, mu, ngp_c, ngp_l, ANS_membrane, ANS_shear, ANS_curvature)

                # Assemble ke_l for last element
                if ie == numel - 1:
                    ke += ke_l
                
                # Assembly
                for i in range(ne):
                    ii = connect[ie, i]
                    for j in range(ne):
                        jj = connect[ie, j]
                        K[3*ii:3*ii+3, 3*jj:3*jj+3] += ke[3*i:3*i+3, 3*j:3*j+3]
                    F_int[3*ii:3*ii+3] += fe[3*i:3*i+3]


            # Residual
            R = F_int - F_ext
            
            # Solve
            du, R_r = solve_linear_system(K, R, prescribed_dofs)

            # Update
            u += du

            # Residual and Energy norms
            ResNorm[inc-1, itn-1] = np.linalg.norm(R_r)
            EnerNorm[inc-1, itn-1] = abs(R.T @ du)[0, 0]

            print(f"Load Step {inc}, Iteration {itn}: Residual Norm = {ResNorm[inc-1, itn-1]:10.6e}, Energy Norm = {EnerNorm[inc-1, itn-1]:10.6e}")

            # Convergence check
            if EnerNorm[inc-1, itn-1] < tolerance:
                print("=====")
                break

        displacements_all.append(u.copy())

    return u, displacements_all



def compute_external_force(E, length, height, width, n_load_steps, inc,
                            coords, connect, ne, ne_L, numel, ngp_c, ndof, u):
    """
    Compute the external force vector and load stiffness matrix for surface traction.

    Parameters
    ----------
    E : float
        Young's modulus.
    length : float
        Initial length of the beam.
    height : float
        Initial height of the cross-section.
    width : float
        Initial width of the cross-section.
    n_load_steps : int
        Total number of load steps.
    inc : int
        Current load step.
    coords : np.ndarray
        Initial nodal coordinates (n_nodes, 3).
    connect : np.ndarray
        Element connectivity array (n_elements, n_nodes_per_element).
    ne : int
        Number of nodes per element.
    ne_L : int
        Number of nodes along the length.
    numel : int
        Number of elements.
    ngp_c : int
        Number of Gauss points per cross-sectional direction.
    ndof : int
        Total number of degrees of freedom.
    u : np.ndarray
        Current displacement vector.

    Returns
    -------
    F_ext : np.ndarray
        External force vector (ndof x 1).
    ke_l : np.ndarray
        Load stiffness matrix (element level).
    """
    ndof_e = 3 * ne
    list_e_surf = [numel - 1]  # Loading only the last element

    # Gauss points and weights
    XI, wgp = gauss(ngp_c)
    xi_fixed = 1.0  # fixed natural coordinate

    eta = XI
    zta = XI

    # Initialize force and stiffness
    F_ext = np.zeros((ndof, 1))
    ke_l = np.zeros((ndof_e, ndof_e))

    # Current nodal coordinates
    qcur = coords + u.reshape(-1, 3)

    for ie in list_e_surf:

        fe_ext = np.zeros((ndof_e, 1))

        for k2 in range(ngp_c):
            for k3 in range(ngp_c):

                # Shape functions and derivatives
                N, dN = shape_functions_3D(ne_L, xi_fixed, eta[k2], zta[k3])

                # Build Ne matrix
                Ne = np.zeros((3, ndof_e))
                Z = 0.0
                for m in range(ne):
                    mm = connect[ie, m]
                    Ne[0, 3*m] = N[m]
                    Ne[1, 3*m+1] = N[m]
                    Ne[2, 3*m+2] = N[m]
                    Z += N[m] * coords[mm, 2]

                # Surface traction calculation
                Ie = height**3 * width / 12
                M = (np.pi * E * Ie / length) * (inc / n_load_steps)  # half circle
                marm = Z - height/2
                sf = (-12 * M / height**3) * marm

                # Convective vectors
                comp1 = np.zeros((3, 1))
                comp2 = np.zeros((3, 1))
                for m in range(ne):
                    mm = connect[ie, m]
                    comp1 += (dN[1, m] * qcur[mm, :]).reshape(3,1)
                    comp2 += (dN[2, m] * qcur[mm, :]).reshape(3,1)

                # Normal vector (cross product)
                normal = np.zeros((3, 1))
                normal[0] = comp1[1]*comp2[2] - comp1[2]*comp2[1]
                normal[1] = comp1[2]*comp2[0] - comp1[0]*comp2[2]
                normal[2] = comp1[0]*comp2[1] - comp1[1]*comp2[0]

                # N_bar matrices
                N_bar_eta = np.zeros((3, 3))
                N_bar_eta[0, 1] = comp1[2]
                N_bar_eta[1, 0] = -comp1[2]
                N_bar_eta[0, 2] = -comp1[1]
                N_bar_eta[2, 0] = comp1[1]
                N_bar_eta[1, 2] = comp1[0]
                N_bar_eta[2, 1] = -comp1[0]

                N_bar_zta = np.zeros((3, 3))
                N_bar_zta[0, 1] = comp2[2]
                N_bar_zta[1, 0] = -comp2[2]
                N_bar_zta[0, 2] = -comp2[1]
                N_bar_zta[2, 0] = comp2[1]
                N_bar_zta[1, 2] = comp2[0]
                N_bar_zta[2, 1] = -comp2[0]

                # Load stiffness matrix
                for K in range(ne):
                    for L in range(ne):
                        ke_l[3*K:3*K+3, 3*L:3*L+3] -= (
                            sf * N[K] * (dN[1, L]*N_bar_zta - dN[2, L]*N_bar_eta) * wgp[k2] * wgp[k3]
                        )

                # External force contribution
                fe_ext += (Ne.T @ (sf * normal)) * wgp[k2] * wgp[k3]

        # Assembly into global force vector
        for m in range(ne):
            mm = connect[ie, m]
            F_ext[3*mm:3*mm+3, 0] += fe_ext[3*m:3*m+3, 0]

    return F_ext, ke_l










def element_3D_nonlinear(ne, ne_L, nx0_e, ue, lambda_, mu, ngp_c, ngp_l, ANS_membrane, ANS_shear, ANS_curvature):
    """
    Compute element tangent stiffness matrix and internal force vector for Saint Venant-Kirchhoff (SVK) material.

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

    ndof_e = 3 * ne

    # Initialize
    ke_m = np.zeros((ndof_e, ndof_e))
    ke_g = np.zeros((ndof_e, ndof_e))
    fe = np.zeros((ndof_e, 1))

    # Gauss points
    eta, wgp_c = gauss(ngp_c)
    zta = eta
    xi, wgp_l = gauss(ngp_l)

    # Current configuration
    qcur = nx0_e + ue

    for k1 in range(ngp_l):
        for k2 in range(ngp_c):
            for k3 in range(ngp_c):

                N, dN = shape_functions_3D(ne_L, xi[k1], eta[k2], zta[k3])

                # Compute covariant basis in reference
                H = np.zeros((3, 3))
                for m in range(ne):
                    H[:, 0] += dN[:, m] * nx0_e[m, 0]
                    H[:, 1] += dN[:, m] * nx0_e[m, 1]
                    H[:, 2] += dN[:, m] * nx0_e[m, 2]

                detJ = np.linalg.det(H)

                G_1 = H[:, 0]
                G_2 = H[:, 1]
                G_3 = H[:, 2]

                Gco = np.array([
                    [np.dot(G_1, G_1), np.dot(G_1, G_2), np.dot(G_1, G_3)],
                    [np.dot(G_2, G_1), np.dot(G_2, G_2), np.dot(G_2, G_3)],
                    [np.dot(G_3, G_1), np.dot(G_3, G_2), np.dot(G_3, G_3)]
                ])

                Gcontra = np.linalg.inv(Gco)
                
                G1 = Gcontra[0, 0]*G_1 + Gcontra[0, 1]*G_2 + Gcontra[0, 2]*G_3
                G2 = Gcontra[1, 0]*G_1 + Gcontra[1, 1]*G_2 + Gcontra[1, 2]*G_3
                G3 = Gcontra[2, 0]*G_1 + Gcontra[2, 1]*G_2 + Gcontra[2, 2]*G_3

                # Covariant basis in current configuration
                H = np.zeros((3, 3))
                for m in range(ne):
                    H[0, :] += dN[:, m] * qcur[m, 0]
                    H[1, :] += dN[:, m] * qcur[m, 1]
                    H[2, :] += dN[:, m] * qcur[m, 2]

                g_1 = H[:, 0]
                g_2 = H[:, 1]
                g_3 = H[:, 2]

                gco = np.array([
                    [np.dot(g_1, g_1), np.dot(g_1, g_2), np.dot(g_1, g_3)],
                    [np.dot(g_2, g_1), np.dot(g_2, g_2), np.dot(g_2, g_3)],
                    [np.dot(g_3, g_1), np.dot(g_3, g_2), np.dot(g_3, g_3)]
                ])

                F = np.outer(g_1, G1) + np.outer(g_2, G2) + np.outer(g_3, G3)

                # Green-Lagrange strain
                E = 0.5 * (gco - Gco)

                # 2nd Piola-Kirchhoff stress (SVK)
                S = np.zeros((3, 3))
                for I in range(3):
                    for J in range(3):
                        for K in range(3):
                            for L in range(3):
                                S[I, J] += compute_C_SVK_component(Gcontra, lambda_, mu, I, J, K, L) * E[K, L]

                # Voigt notation
                Svgt = np.array([
                    S[0, 0],
                    S[1, 1],
                    S[2, 2],
                    S[0, 1],
                    S[0, 2],
                    S[1, 2]
                ]).reshape(6, 1)

                # Strain-displacement matrix Be
                Be = np.zeros((6, ndof_e))

                for m in range(ne):
                    B = np.zeros((6, 3))
                    B[0, :] = g_1 * dN[0, m]
                    B[1, :] = g_2 * dN[1, m]
                    B[2, :] = g_3 * dN[2, m]
                    B[3, :] = g_2 * dN[0, m] + g_1 * dN[1, m]
                    B[4, :] = g_3 * dN[0, m] + g_1 * dN[2, m]
                    B[5, :] = g_3 * dN[1, m] + g_2 * dN[2, m]
                    Be[:, 3*m:3*m+3] = B


                # Modify B using ANS for membrane and transverse shear locking
                if ANS_membrane or ANS_shear:

                    ntying1 = ne_L - 1
                    xibar1 = gauss(ntying1)
                    Nmat1 = lagrange_shape_1D(ntying1, xi[k1])

                    M1 = np.zeros((ntying1, ntying1))
                    for i in range(ntying1):
                        M1[i, :] = lagrange_shape_1D(ntying1, xibar1[i])

                    if ANS_membrane:
                        Bebar_1 = np.zeros((ntying1, 3 * ne))

                    if ANS_shear:
                        Bebar_4 = np.zeros((ntying1, 3 * ne))
                        Bebar_5 = np.zeros((ntying1, 3 * ne))

                    for k in range(ntying1):
                        Nbar1, dNbar1 = shape_functions_3D(ne_L, xibar1[k], eta[k2], zta[k3])

                        # Covariant tangents
                        Pbar1 = np.zeros((3, 3))
                        for m in range(ne):
                            Pbar1[0, :] += dNbar1[0:3, m] * qcur[m, 0]
                            Pbar1[1, :] += dNbar1[0:3, m] * qcur[m, 1]
                            Pbar1[2, :] += dNbar1[0:3, m] * qcur[m, 2]
                    
                        gbar_1 = Pbar1[:, 0]
                        gbar_2 = Pbar1[:, 1]
                        gbar_3 = Pbar1[:, 2]

                        if ANS_membrane:
                            for m in range(ne):
                                Bbar1 = np.zeros(3)
                                Bbar1[0] = gbar_1[0] * dNbar1[0, m]
                                Bbar1[1] = gbar_1[1] * dNbar1[0, m]
                                Bbar1[2] = gbar_1[2] * dNbar1[0, m]
                                Bebar_1[k, 3 * m : 3 * m + 3] = Bbar1
                                
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

                        # Final row updates
                        if ANS_membrane:
                            Row1 = Nmat1 @ np.linalg.solve(M1, Bebar_1)
                            Be[0, :] = Row1

                        if ANS_shear:
                            Row4 = Nmat1 @ np.linalg.solve(M1, Bebar_4)
                            Be[3, :] = Row4
                            Row5 = Nmat1 @ np.linalg.solve(M1, Bebar_5)
                            Be[4, :] = Row5
                

                # Modify B using ANS for curvature-thickness locking
                if ANS_curvature:

                    ntying2 = ne_L
                    xibar2 = np.linspace(-1, 1, ntying2)  # Use equally spaced tying points

                    Nmat2 = lagrange_shape_1D(ntying2, xi[k1])
                    M2 = np.eye(ntying2)

                    Bebar_2 = np.zeros((ntying2, 3 * ne))
                    Bebar_3 = np.zeros((ntying2, 3 * ne))

                    for k in range(ntying2):
                        Nbar2, dNbar2 = shape_functions_3D(ne_L, xibar2[k], eta[k2], zta[k3])

                        Pbar2 = np.zeros((3, 3))
                        for m in range(ne):
                            Pbar2[0, :] += dNbar2[0:3, m] * qcur[m, 0]
                            Pbar2[1, :] += dNbar2[0:3, m] * qcur[m, 1]
                            Pbar2[2, :] += dNbar2[0:3, m] * qcur[m, 2]

                        gbar_1 = Pbar2[:, 0]
                        gbar_2 = Pbar2[:, 1]
                        gbar_3 = Pbar2[:, 2]

                        # Calculate 2nd row
                        for m in range(ne):
                            Bbar2 = np.zeros(3)
                            Bbar2[0] = gbar_2[0] * dNbar2[1, m]
                            Bbar2[1] = gbar_2[1] * dNbar2[1, m]
                            Bbar2[2] = gbar_2[2] * dNbar2[1, m]
                            Bebar_2[k, 3 * m : 3 * m + 3] = Bbar2

                        # Calculate 3rd row
                        for m in range(ne):
                            Bbar3 = np.zeros(3)
                            Bbar3[0] = gbar_3[0] * dNbar2[2, m]
                            Bbar3[1] = gbar_3[1] * dNbar2[2, m]
                            Bbar3[2] = gbar_3[2] * dNbar2[2, m]
                            Bebar_3[k, 3 * m : 3 * m + 3] = Bbar3

                    # Replace 2nd and 3rd rows of B matrix
                    Row2 = Nmat2 @ np.linalg.solve(M2, Bebar_2)
                    Be[1, :] = Row2

                    Row3 = Nmat2 @ np.linalg.solve(M2, Bebar_3)
                    Be[2, :] = Row3




                # Internal force
                fe += Be.T @ Svgt * detJ * wgp_l[k1] * wgp_c[k2] * wgp_c[k3]

                # Material stiffness matrix (C)
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

                # Preallocate tying arrays as needed
                if ANS_membrane:
                    tying_dNxi_dNxi = np.zeros((ntying1, 1))

                if ANS_shear:
                    tying_dNxi_dNeta = np.zeros((ntying1, 1))
                    tying_dNeta_dNxi = np.zeros((ntying1, 1))
                    tying_dNxi_dNzta = np.zeros((ntying1, 1))
                    tying_dNzta_dNxi = np.zeros((ntying1, 1))
                
                if ANS_curvature:
                    tying_dNeta_dNeta = np.zeros((ntying2, 1))
                    tying_dNzta_dNzta = np.zeros((ntying2, 1))

                for K in range(ne):
                    for L in range(ne):

                        factor = detJ * wgp_l[k1] * wgp_c[k2] * wgp_c[k3]
                        I3 = np.eye(3)

                        # Always calculated (never modified)
                        ke_g23[3*K:3*K+3, 3*L:3*L+3] += dN[1, K] * S[1, 2] * I3 * dN[2, L] * factor
                        ke_g32[3*K:3*K+3, 3*L:3*L+3] += dN[2, K] * S[2, 1] * I3 * dN[1, L] * factor

                        # --- Membrane (xi-xi) ---
                        if ANS_membrane:
                            for k in range(ntying1):
                                _, dNbar1 = shape_functions_3D(ne_L, xibar1[k], eta[k2], zta[k3])
                                tying_dNxi_dNxi[k, 0] = dNbar1[0, K] * dNbar1[0, L]
                            proj_dNxi_dNxi = Nmat1 @ np.linalg.inv(M1) @ tying_dNxi_dNxi
                            Sxi_xi = S[0, 0] * proj_dNxi_dNxi
                        else:
                            Sxi_xi = S[0, 0] * dN[0, K] * dN[0, L]
                        ke_g11[3*K:3*K+3, 3*L:3*L+3] += Sxi_xi * I3 * factor

                        # --- Shear terms (modified if ANS_shear) ---
                        if ANS_shear:
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
                        else:
                            Sxi_eta = S[0, 1] * dN[0, K] * dN[1, L]
                            Seta_xi = S[1, 0] * dN[1, K] * dN[0, L]
                            Sxi_zta = S[0, 2] * dN[0, K] * dN[2, L]
                            Szta_xi = S[2, 0] * dN[2, K] * dN[0, L]

                        ke_g12[3*K:3*K+3, 3*L:3*L+3] += Sxi_eta * I3 * factor
                        ke_g21[3*K:3*K+3, 3*L:3*L+3] += Seta_xi * I3 * factor
                        ke_g13[3*K:3*K+3, 3*L:3*L+3] += Sxi_zta * I3 * factor
                        ke_g31[3*K:3*K+3, 3*L:3*L+3] += Szta_xi * I3 * factor

                        # --- Curvature terms (modified if ANS_curvature) ---
                        if ANS_curvature:
                            for k in range(ntying2):
                                _, dNbar2 = shape_functions_3D(ne_L, xibar2[k], eta[k2], zta[k3])
                                tying_dNeta_dNeta[k, 0] = dNbar2[1, K] * dNbar2[1, L]
                                tying_dNzta_dNzta[k, 0] = dNbar2[2, K] * dNbar2[2, L]

                            proj_dNeta_dNeta = Nmat2 @ np.linalg.inv(M2) @ tying_dNeta_dNeta
                            proj_dNzta_dNzta = Nmat2 @ np.linalg.inv(M2) @ tying_dNzta_dNzta

                            Seta_eta = S[1, 1] * proj_dNeta_dNeta
                            Szta_zta = S[2, 2] * proj_dNzta_dNzta
                        else:
                            Seta_eta = S[1, 1] * dN[1, K] * dN[1, L]
                            Szta_zta = S[2, 2] * dN[2, K] * dN[2, L]

                        ke_g22[3*K:3*K+3, 3*L:3*L+3] += Seta_eta * I3 * factor
                        ke_g33[3*K:3*K+3, 3*L:3*L+3] += Szta_zta * I3 * factor












                # Sum all parts to form final geometric stiffness matrix
                ke_g += (
                    ke_g11 + ke_g12 + ke_g13 +
                    ke_g21 + ke_g22 + ke_g23 +
                    ke_g31 + ke_g32 + ke_g33
)






    ke = ke_m + ke_g
    return ke, fe























def gauss(ngp):
    """
    Returns Gauss points and weights for 1D numerical integration.

    Parameters
    ----------
    ngp : int
        Number of Gauss points (1 to 10 supported).

    Returns
    -------
    XI : np.ndarray
        Coordinates of Gauss points.
    ALPHAG : np.ndarray
        Weights associated with Gauss points.
    """
    XI = np.zeros(ngp)
    ALPHAG = np.zeros(ngp)

    if ngp == 1:
        XI[0] = 0.0
        ALPHAG[0] = 2.0

    elif ngp == 2:
        XI[:] = [-1/np.sqrt(3), 1/np.sqrt(3)]
        ALPHAG[:] = [1.0, 1.0]

    elif ngp == 3:
        XI[:] = [-np.sqrt(3/5), 0.0, np.sqrt(3/5)]
        ALPHAG[:] = [5/9, 8/9, 5/9]

    elif ngp == 4:
        XI[:] = [-0.861136311594053, -0.339981043584856,
                  0.339981043584856,  0.861136311594053]
        ALPHAG[:] = [0.347854845137454, 0.652145154862546,
                     0.652145154862546, 0.347854845137454]

    elif ngp == 5:
        XI[:] = [-0.906179845938664, -0.538469310105683, 0.0,
                  0.538469310105683,  0.906179845938664]
        ALPHAG[:] = [0.236926885056189, 0.478628670499366,
                     0.568888888888889, 0.478628670499366,
                     0.236926885056189]

    elif ngp == 6:
        XI[:] = [-0.932469514203152, -0.661209386466265, -0.238619186083197,
                  0.238619186083197,  0.661209386466265,  0.932469514203152]
        ALPHAG[:] = [0.171324492379170, 0.360761573048139, 0.467913934572691,
                     0.467913934572691, 0.360761573048139, 0.171324492379170]

    elif ngp == 7:
        XI[:] = [-0.949107912342759, -0.741531185599394, -0.405845151377397,
                  0.0,
                  0.405845151377397,  0.741531185599394,  0.949107912342759]
        ALPHAG[:] = [0.129484966168870, 0.279705391489277, 0.381830050505119,
                     0.417959183673469,
                     0.381830050505119, 0.279705391489277, 0.129484966168870]

    elif ngp == 8:
        XI[:] = [-0.960289856497536, -0.796666477413627, -0.525532409916329,
                 -0.183434642495650,
                  0.183434642495650,  0.525532409916329,  0.796666477413627,
                  0.960289856497536]
        ALPHAG[:] = [0.101228536290376, 0.222381034453374, 0.313706645877887,
                     0.362683783378362,
                     0.362683783378362, 0.313706645877887, 0.222381034453374,
                     0.101228536290376]

    elif ngp == 9:
        XI[:] = [-0.968160239507626, -0.836031107326636, -0.613371432700590,
                 -0.324253423403809, 0.0,
                  0.324253423403809,  0.613371432700590,  0.836031107326636,
                  0.968160239507626]
        ALPHAG[:] = [0.081274388361574, 0.180648160694857, 0.260610696402935,
                     0.312347077040003, 0.330239355001260,
                     0.312347077040003, 0.260610696402935, 0.180648160694857,
                     0.081274388361574]

    elif ngp == 10:
        XI[:] = [-0.973906528517172, -0.865063366688985, -0.679409568299024,
                 -0.433395394129247, -0.148874338981631,
                  0.148874338981631,  0.433395394129247,  0.679409568299024,
                  0.865063366688985,  0.973906528517172]
        ALPHAG[:] = [0.066671344308688, 0.149451349150581, 0.219086362515982,
                     0.269266719309996, 0.295524224714753,
                     0.295524224714753, 0.269266719309996, 0.219086362515982,
                     0.149451349150581, 0.066671344308688]

    else:
        raise ValueError("Number of Gauss points must be between 1 and 10.")

    return XI, ALPHAG




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

    # Reduced system solve
    R_r = R[free_dofs]
    s = -np.linalg.solve(K[np.ix_(free_dofs, free_dofs)], R_r)

    # Set displacements
    d[free_dofs, 0] = s[:, 0]
    d[prescribed_dofs, 0] = 0.0

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
    Computes shape functions and their derivatives for a nonlinear 3D brick element.

    Parameters
    ----------
    ne_L : int
        Number of nodes along the length (xi) direction.
    xi, eta, zta : float
        Parametric coordinates.

    Returns
    -------
    N : np.ndarray
        Shape function values at (xi, eta, zta), shape (4*ne_L, ).
    dN : np.ndarray
        Derivatives of shape functions w.r.t (xi, eta, zta), shape (3, 4*ne_L).
        dN[0, :] = dN/dxi, dN[1, :] = dN/deta, dN[2, :] = dN/dzta
    """
    # Parametric coordinates of nodes in cross-sectional directions (eta and zta)
    eta_nodes = [-1, 1]
    zta_nodes = [-1, 1]

    # Linear Lagrangian shape functions in eta and zta
    N_eta = [0.5 * (1 - eta), 0.5 * (1 + eta)]
    N_zta = [0.5 * (1 - zta), 0.5 * (1 + zta)]

    # Parametric coordinates of nodes along xi direction
    xi_nodes = np.linspace(-1, 1, ne_L)

    # 1D Lagrangian shape functions in xi
    N_xi = np.zeros(len(xi_nodes))
    for k in range(len(xi_nodes)):
        product = 1.0
        for m in range(len(xi_nodes)):
            if m != k:
                product *= (xi_nodes[m] - xi) / (xi_nodes[m] - xi_nodes[k])
        N_xi[k] = product

    # Combine shape functions in all directions
    N = np.zeros(4 * len(xi_nodes))
    num = 0
    for k in range(len(xi_nodes)):
        for j in range(2):
            for i in range(2):
                N[num] = N_eta[i] * N_zta[j] * N_xi[k]
                num += 1

    # Initialize derivatives dN (3 rows: dN/dxi, dN/deta, dN/dzta)
    dN = np.zeros((3, 4 * len(xi_nodes)))

    # Derivatives w.r.t eta (row 1)
    num = 0
    for k in range(len(xi_nodes)):
        for j in range(2):
            for i in range(2):
                dN[1, num] = ((-0.5 if eta_nodes[i] == -1 else 0.5) * N_zta[j] * N_xi[k])
                num += 1

    # Derivatives w.r.t zta (row 2)
    num = 0
    for k in range(len(xi_nodes)):
        for j in range(2):
            for i in range(2):
                dN[2, num] = ((-0.5 if zta_nodes[j] == -1 else 0.5) * N_eta[i] * N_xi[k])
                num += 1

    # Derivatives w.r.t xi (row 0)
    dXi = np.zeros(len(xi_nodes))
    for i in range(len(xi_nodes)):
        dl = 0.0
        for j in range(len(xi_nodes)):
            if j != i:
                dl_temp = 1.0 / (xi_nodes[i] - xi_nodes[j])
                for k in range(len(xi_nodes)):
                    if k != i and k != j:
                        dl_temp *= (xi - xi_nodes[k]) / (xi_nodes[i] - xi_nodes[k])
                dl += dl_temp
        dXi[i] = dl

    num = 0
    num1 = 0
    for k in range(len(xi_nodes)):
        for j in range(2):
            for i in range(2):
                dN[0, num] = dXi[num1] * N_eta[i] * N_zta[j]
                num += 1
        num1 += 1

    return N, dN


def voigt_index(i, j):
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
