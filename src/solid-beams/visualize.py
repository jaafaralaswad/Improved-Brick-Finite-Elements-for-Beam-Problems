import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def plot_mesh_3D(fname: str, coords: np.ndarray, connect: np.ndarray, 
                 prescribed_nodes: np.ndarray = None, show_node_numbers: bool = False, 
                 elev: float = 20, azim: float = -60):
    """
    Plots a 3D mesh based on nodal coordinates.

    Parameters
    ----------
    fname : str
        Filename to save the plot.
    coords : np.ndarray of shape (n_nodes, 3)
        Array of nodal coordinates.
    connect : np.ndarray
        Element connectivity array.
    prescribed_nodes : np.ndarray, optional
        List of node indices with prescribed displacements (highlighted in red).
    show_node_numbers : bool, optional
        Whether to display node indices next to the points (default is False).
    elev : float, optional
        Elevation angle for the 3D plot view.
    azim : float, optional
        Azimuth angle for the 3D plot view.
    """
    coords = coords.reshape(-1, 3)

    fig = plt.figure(figsize=(8, 6), dpi=150)
    ax = fig.add_subplot(111, projection='3d')

    xs, ys, zs = coords[:, 0], coords[:, 1], coords[:, 2]

    ax.set_xticks(np.unique(np.round(xs, decimals=6)))
    ax.set_yticks(np.unique(np.round(ys, decimals=6)))
    ax.set_zticks(np.unique(np.round(zs, decimals=6)))

    # Scatter plot: different colors for clamped and free nodes
    if prescribed_nodes is None:
        prescribed_nodes = []
    prescribed_nodes = np.array(prescribed_nodes, dtype=int)

    all_nodes = np.arange(coords.shape[0])
    free_nodes = np.setdiff1d(all_nodes, prescribed_nodes)

    ax.scatter(xs[free_nodes], ys[free_nodes], zs[free_nodes], 
               c='blue', s=10, depthshade=True, label='Free Nodes')
    ax.scatter(xs[prescribed_nodes], ys[prescribed_nodes], zs[prescribed_nodes], 
               c='red', s=20, depthshade=True, label='Clamped Nodes')

    # Draw edges                
    for element in connect:
        element = element.astype(int)
        nslices = len(element) // 4

        for s in range(nslices):
            i = s * 4
            n0, n1, n2, n3 = element[i:i+4]

            # Connect slice
            ax.plot(*zip(coords[n0], coords[n1]), color='gray', linewidth=1)
            ax.plot(*zip(coords[n0], coords[n2]), color='gray', linewidth=1)
            ax.plot(*zip(coords[n1], coords[n3]), color='gray', linewidth=1)
            ax.plot(*zip(coords[n2], coords[n3]), color='gray', linewidth=1)

            # Connect to next slice
            if s < nslices - 1:
                j = (s + 1) * 4
                m0, m1, m2, m3 = element[j:j+4]
                ax.plot(*zip(coords[n0], coords[m0]), color='gray', linewidth=1)
                ax.plot(*zip(coords[n1], coords[m1]), color='gray', linewidth=1)
                ax.plot(*zip(coords[n2], coords[m2]), color='gray', linewidth=1)
                ax.plot(*zip(coords[n3], coords[m3]), color='gray', linewidth=1)

    # Show node numbers if requested
    if show_node_numbers:
        for i, (x, y, z) in enumerate(coords):
            ax.text(x, y, z, str(i), fontsize=8, color='black',
                    verticalalignment='bottom', horizontalalignment='right')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(elev=elev, azim=azim)

    # Set aspect ratio
    xrange = np.ptp(xs)
    yrange = np.ptp(ys)
    zrange = np.ptp(zs)
    ax.set_box_aspect([xrange, yrange, zrange])

    ax.grid(True)
    ax.set_title('Reference Configuration', fontsize=14, fontweight='bold')

    # Add legend
    ax.legend()

    plt.tight_layout()
    plt.savefig(fname, dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    plt.close(fig)





def plot_deformed_mesh_3D(fname: str, coords: np.ndarray, connect: np.ndarray, 
                          u: np.ndarray, scale: float = 1.0, 
                          prescribed_nodes: np.ndarray = None, show_node_numbers: bool = False, 
                          elev: float = 20, azim: float = -60):
    """
    Plots the deformed 3D mesh.

    Parameters
    ----------
    fname : str
        Filename to save the plot.
    coords : np.ndarray of shape (n_nodes, 3)
        Original nodal coordinates.
    connect : np.ndarray
        Element connectivity array.
    u : np.ndarray of shape (n_dofs,)
        Global displacement vector.
    scale : float, optional
        Scale factor for visualizing displacements (default is 1.0).
    prescribed_nodes : np.ndarray, optional
        List of node indices with prescribed displacements (highlighted in red).
    show_node_numbers : bool, optional
        Whether to display node indices (default is False).
    elev : float, optional
        Elevation angle for 3D view.
    azim : float, optional
        Azimuth angle for 3D view.
    """
    coords = coords.reshape(-1, 3)
    u = u.reshape(-1, 3)  # Reshape u to (n_nodes, 3)

    # Compute deformed coordinates
    coords_def = coords + scale * u

    fig = plt.figure(figsize=(8, 6), dpi=150)
    ax = fig.add_subplot(111, projection='3d')

    xs, ys, zs = coords_def[:, 0], coords_def[:, 1], coords_def[:, 2]

    ax.set_xticks(np.unique(np.round(xs, decimals=6)))
    ax.set_yticks(np.unique(np.round(ys, decimals=6)))
    ax.set_zticks(np.unique(np.round(zs, decimals=6)))

    if prescribed_nodes is None:
        prescribed_nodes = []
    prescribed_nodes = np.array(prescribed_nodes, dtype=int)

    all_nodes = np.arange(coords.shape[0])
    free_nodes = np.setdiff1d(all_nodes, prescribed_nodes)

    ax.scatter(xs[free_nodes], ys[free_nodes], zs[free_nodes], 
               c='blue', s=10, depthshade=True, label='Free Nodes')
    ax.scatter(xs[prescribed_nodes], ys[prescribed_nodes], zs[prescribed_nodes], 
               c='red', s=20, depthshade=True, label='Clamped Nodes')

    # Draw edges
    for element in connect:
        element = element.astype(int)
        nslices = len(element) // 4

        for s in range(nslices):
            i = s * 4
            n0, n1, n2, n3 = element[i:i+4]

            ax.plot(*zip(coords_def[n0], coords_def[n1]), color='gray', linewidth=1)
            ax.plot(*zip(coords_def[n0], coords_def[n2]), color='gray', linewidth=1)
            ax.plot(*zip(coords_def[n1], coords_def[n3]), color='gray', linewidth=1)
            ax.plot(*zip(coords_def[n2], coords_def[n3]), color='gray', linewidth=1)

            if s < nslices - 1:
                j = (s + 1) * 4
                m0, m1, m2, m3 = element[j:j+4]
                ax.plot(*zip(coords_def[n0], coords_def[m0]), color='gray', linewidth=1)
                ax.plot(*zip(coords_def[n1], coords_def[m1]), color='gray', linewidth=1)
                ax.plot(*zip(coords_def[n2], coords_def[m2]), color='gray', linewidth=1)
                ax.plot(*zip(coords_def[n3], coords_def[m3]), color='gray', linewidth=1)

    if show_node_numbers:
        for i, (x, y, z) in enumerate(coords_def):
            ax.text(x, y, z, str(i), fontsize=8, color='black',
                    verticalalignment='bottom', horizontalalignment='right')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.view_init(elev=elev, azim=azim)

    # Set aspect ratio
    xrange = np.ptp(xs)
    yrange = np.ptp(ys)
    zrange = np.ptp(zs)
    ax.set_box_aspect([xrange, yrange, zrange])

    ax.grid(True)
    ax.set_title('Deformed Configuration', fontsize=14, fontweight='bold')
    ax.legend()

    plt.tight_layout()
    plt.savefig(fname, dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.show()
    plt.close(fig)