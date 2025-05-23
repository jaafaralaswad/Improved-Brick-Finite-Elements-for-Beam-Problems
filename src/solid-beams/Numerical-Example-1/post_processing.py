# Author: Jaafar Alaswad (jalaswad@bu.edu)

import numpy as np
import matplotlib.pyplot as plt

def compute_tip_displacement(u, coords, connect, shape_func_3D, ne_L):
    """
    Compute interpolated displacement at the center of the tip of the beam.

    Parameters
    ----------
    u : np.ndarray
        Global displacement vector of shape (n_dofs,).
    coords : np.ndarray
        Nodal coordinates of shape (n_nodes, 3).
    connect : np.ndarray
        Element connectivity array of shape (n_elements, nodes_per_element).
    shape_func_3D : callable
        Shape function evaluator: shape_func_3D(ne_L, xi, eta, zeta) → (N, dN)
    ne_L : int
        Number of nodes per element in axial direction.

    Returns
    -------
    u_center : np.ndarray of shape (3,)
        Interpolated displacement at tip center.
    """

    # Reshape global displacement vector into (n_nodes, 3) form
    u_nodes = u.flatten().reshape(-1, 3)

    # Get the global node indices of the last element (assumed to be the tip)
    tip_node_indices = connect[-1]

    # Extract displacement vectors for the tip element nodes
    u_tip_nodes = u_nodes[tip_node_indices]

    # Evaluate shape functions at the center of the tip face:
    # xi = 1.0 (beam tip), eta = 0, zeta = 0 → center of the cross-section
    N, _ = shape_func_3D(ne_L, 1.0, 0.0, 0.0)

    # Interpolate displacement at the center using the shape functions
    u_center = N @ u_tip_nodes

    return u_center[0]


def plot_tip_displacement_x(tip_displacements, length, save_path="tip_displacement_x.png"):
    """
    Plot normalized X displacement at the beam tip over normalized load.

    Parameters
    ----------
    tip_displacements : list of np.ndarray
        List of tip displacements (each a 3-element vector).
    length : float
        Beam length (used for normalization).
    save_path : str
        File path to save the figure.
    """


    disp_array = np.array(tip_displacements)
    x_disp = np.abs(disp_array)   # Absolute value of x displacement

    # Add zero at the beginning for undeformed "step"
    x_disp = np.insert(x_disp, 0, 0.0)

    # Normalize by beam length
    x_disp_normalized = x_disp / length
    
    # Compute normalized load parameter k over the same number of steps
    n_steps = len(x_disp_normalized)
    k_normalized = np.linspace(0, 0.5, n_steps)  # k = ML / (2πEI), ranges from 0 to 0.5

    # analytical solution from Euler's elastica for bending into half a circle
    analytical = np.array([
        0, 0.016368357, 0.064510716, 0.141606309, 0.243173271,
        0.363380228, 0.495448848, 0.632116989, 0.766127679,
        0.890707595, 1
    ])
    analytical_k = np.linspace(0, 0.5, len(analytical))

    # Plot settings
    plt.figure(figsize=(6, 4))
    plt.plot(k_normalized, x_disp_normalized, marker='o', linestyle='-', color='blue', label='Numerical Solution')
    plt.plot(analytical_k, analytical, linestyle='-', color='black', label='Analytical Solution')
    plt.xlabel(r"$k = \dfrac{ML}{2\pi EI}$", fontsize=12)
    plt.ylabel(r"$\dfrac{u}{L}$", fontsize=12)
    plt.grid(True)
    plt.legend()
    plt.xlim(0, 0.5)
    plt.ylim(0, 1.2)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.show()