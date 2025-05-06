# Author: Jaafar Alaswad (jalaswad@bu.edu)

import numpy as np

def validate_inputs(width, height, length, E, nu, numel, ne_L, ngp_c, ngp_l, n_load_steps, max_iterations, tolerance):
    """
    Validate the input parameters for the finite element simulation.
    
    Raises:
        ValueError: If any input is outside the allowed range.
    """
    # Geometry checks
    if width <= 0 or height <= 0 or length <= 0:
        raise ValueError("All geometry dimensions (width, height, length) must be positive.")

    # Material properties
    if E <= 0:
        raise ValueError("Young's modulus 'E' must be positive.")
    if not (0 <= nu < 0.5):
        raise ValueError("Poisson's ratio 'nu' must be in the range [0, 0.5).")

    # Discretization parameters
    if numel < 1:
        raise ValueError("'numel' (number of elements) must be at least 1.")
    if ne_L < 2:
        raise ValueError("'ne_L' (nodes per element in axial direction) must be at least 2.")
    if not (1 <= ngp_c <= 10):
        raise ValueError("'ngp_c' (Gauss points per cross-sectional direction) must be in [1, 10].")
    if not (1 <= ngp_l <= 10):
        raise ValueError("'ngp_l' (Gauss points in axial direction) must be in [1, 10].")

    # Solver parameters
    if n_load_steps < 1:
        raise ValueError("'n_load_steps' must be at least 1.")
    if max_iterations < 1:
        raise ValueError("'max_iterations' must be at least 1.")
    if tolerance <= 0:
        raise ValueError("'tolerance' must be positive.")


def create_mesh(width, height, length, numel, ne_L):
    """
    Parameters:
    width, height: Cross-sectional dimensions
    length: Beam length
    numel: Total number of elements
    ne_L: Number of nodes per element in the length direction

    Returns:
    coords: Nodal coordinates
    nnode_W, nnode_H, nnode_L: Number of nodes in each direction
    ndof: Total number of degrees of freedom
    connect: Element connectivity matrix
    """

    ne = 4 * ne_L           # Total number of nodes per element (solid-beam formulation)

    # Number of nodes per direction
    nnode_W = 2
    nnode_H = 2
    nnode_L = numel + 1 + (ne_L - 2) * numel

    # Total number of nodes and degrees of freedom
    nnode = nnode_W * nnode_H * nnode_L
    ndof = 3 * nnode

    # Element sizes
    elen_W = width
    elen_H = height
    elen_L = length / numel
    nelen_L = elen_L / (ne_L - 1) # Distance between nodes in the length direction

    # Initialize nodal coordinate vector
    coords = np.zeros(ndof)

    # Global node counter
    num = 0

    # Loop over nodes in length direction
    for i in range(nnode_L):
        x_L = nelen_L * i
        for j in range(nnode_H):
            x_H = elen_H * j
            for k in range(nnode_W):
                x_W = elen_W * k

                coords[3*num] = x_L       # x coordinate (along beam)
                coords[3*num + 1] = x_W   # y coordinate
                coords[3*num + 2] = x_H   # z coordinate
                num += 1
    coords = coords.reshape(-1, 3)

    # Connectivity matrix
    connect = np.zeros((numel, ne), dtype=int)
    for i in range(numel):
        for j in range(ne):
            connect[i, j] = j + i * nnode_W * nnode_H + i * (ne_L - 2) * nnode_W * nnode_H

    return coords, nnode_W, nnode_H, nnode_L, ndof, connect, ne







def impose_supports(prescribed_nodes):
    """
    Generate list of constrained degrees of freedom (DOFs) based on prescribed nodes.

    Parameters:
        prescribed_nodes (list or array): Nodes where displacement is prescribed (0-based indexing).

    Returns:
        prescribed_dofs (np.ndarray): Array of constrained DOF indices.
    """
    # Number of prescribed nodes
    num_prescribed_nodes = len(prescribed_nodes)

    # Total number of constrained degrees of freedom (3 dofs per node)
    num_constrained_dofs = 3 * num_prescribed_nodes

    # Degrees of freedom with prescribed displacements
    prescribed_dofs = np.zeros(num_constrained_dofs, dtype=int)

    # List of constrained dofs
    for i, node in enumerate(prescribed_nodes):
        prescribed_dofs[3*i    ] = 3*node     # x DOF
        prescribed_dofs[3*i + 1] = 3*node + 1 # y DOF
        prescribed_dofs[3*i + 2] = 3*node + 2 # z DOF

    return prescribed_dofs


def compute_lame_parameters(E, nu):
    """
    Calculate the Lamé parameters (lambda and mu) given Young's modulus (E) and Poisson's ratio (nu).

    Parameters:
        E (float): Young's modulus
        nu (float): Poisson's ratio
    
    Returns:
        tuple: A tuple containing lambda_ (first Lamé parameter) and mu (second Lamé parameter, shear modulus)
    """
    lambda_ = E * nu / ((1 + nu) * (1 - 2 * nu))  # First Lamé parameter
    mu = E / (2 * (1 + nu))                       # Second Lamé parameter (Shear modulus)
    
    return lambda_, mu