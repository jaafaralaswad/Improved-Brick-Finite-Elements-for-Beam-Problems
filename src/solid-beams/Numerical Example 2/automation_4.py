# This automation code was created using ChatGPT

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from pre_processing import validate_inputs, create_mesh, impose_supports, compute_lame_parameters
from processing import newton_raphson_solver, shape_functions_3D
from post_processing import compute_tip_displacement

# ==============================
#       SETUP FIGURE DIRECTORY
# ==============================

figures_dir = Path(__file__).resolve().parent / "generated-figures"
figures_dir.mkdir(parents=True, exist_ok=True)

# ==============================
#       COMMON PARAMETERS
# ==============================

width = 1.0
height = 1.0
E = 12
nu = 0.0
numel = 3
ne_L = 3
ngp_c = 2
ngp_l = ne_L
n_load_steps = 10
max_iterations = 20
tolerance = 1e-15
visualize_reference = False
visualize_final = False

# ==============================
#       ANS CONFIGURATIONS
# ==============================

cases = [
    {"label": "a = 100 - No locking alleviated", "length": 300, "ANS_membrane": False, "ANS_shear": False, "ANS_curvature": False},
    {"label": "a = 100 - Only membrane locking alleviated", "length": 300, "ANS_membrane": True, "ANS_shear": False, "ANS_curvature": False},
    {"label": "a = 100 - Only transverse shear locking alleviated", "length": 300, "ANS_membrane": False, "ANS_shear": True, "ANS_curvature": False},
    {"label": "a = 100 - Only curvature-thickness locking alleviated", "length": 300, "ANS_membrane": False, "ANS_shear": False, "ANS_curvature": True},
    {"label": "a = 100 - All locking modes alleviated", "length": 300, "ANS_membrane": True, "ANS_shear": True, "ANS_curvature": True},
]

# ==============================
#       RUN SIMULATIONS
# ==============================

tip_displacements_all = []

for case in cases:
    length = case["length"]
    validate_inputs(width, height, length, E, nu, numel, ne_L, ngp_c, ngp_l,
                    n_load_steps, max_iterations, tolerance)

    coords, nnode_W, nnode_H, nnode_L, ndof, connect, ne = create_mesh(width, height, length, numel, ne_L)

    prescribed_nodes = [0, 1, 2, 3]
    prescribed_dofs = impose_supports(prescribed_nodes)
    lambda_, mu = compute_lame_parameters(E, nu)

    u, displacements_all = newton_raphson_solver(
        E, width, height, length, lambda_, mu,
        coords, connect, prescribed_dofs,
        numel, ne, ne_L, ngp_c, ngp_l,
        n_load_steps, max_iterations, tolerance,
        case["ANS_membrane"], case["ANS_shear"], case["ANS_curvature"]
    )

    tip_displacements = [
        compute_tip_displacement(u_i, coords, connect, shape_functions_3D, ne_L)
        for u_i in displacements_all
    ]

    tip_displacements_all.append((case["label"], tip_displacements, length))

# ==============================
#       PLOT RESULTS
# ==============================

plt.figure(figsize=(7, 5))

markers = {
    "a = 100 - No locking alleviated": "o",
    "a = 100 - Only membrane locking alleviated": "s",
    "a = 100 - Only transverse shear locking alleviated": "^",
    "a = 100 - Only curvature-thickness locking alleviated": "X",
    "a = 100 - All locking modes alleviated": "d",
}

linestyles = {
    "a = 100 - No locking alleviated": "-",
    "a = 100 - Only membrane locking alleviated": "--",
    "a = 100 - Only transverse shear locking alleviated": "-.",
    "a = 100 - Only curvature-thickness locking alleviated": (0, (3, 1, 1, 1)),
    "a = 100 - All locking modes alleviated": ":",
}

for label, tips, length in tip_displacements_all:
    disp_array = np.array(tips)
    z_disp = disp_array[:, 2]  # Use signed value, not abs
    z_disp = np.insert(z_disp, 0, 0.0)
    z_disp_normalized = z_disp / length
    n_steps = len(z_disp_normalized)
    k_normalized = np.linspace(0, 4.0, n_steps)

    plt.plot(
        k_normalized, z_disp_normalized,
        marker=markers[label],
        linestyle=linestyles[label],
        linewidth=2,
        label=label
    )

plt.xlabel(r"$k = \dfrac{P L^2}{EI}$", fontsize=12)
plt.ylabel(r"Normalized vertical tip displacement $w / L$", fontsize=12)
plt.title("Influence of Indiviual Locking Modes", fontsize=14)
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig(figures_dir / "locking_alleviation_comparison.png", dpi=300)
plt.close()