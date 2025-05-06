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

figures_dir = Path(__file__).resolve().parents[2] / "generated-figures"
figures_dir.mkdir(parents=True, exist_ok=True)

# ==============================
#       COMMON PARAMETERS
# ==============================

width = 1.0
height = 0.5
length = 20.0
E = 1.2e7
nu = 0.0
numel = 20
ne_L = 2
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
    {"label": "No locking alleviated", "ANS_membrane": False, "ANS_shear": False, "ANS_curvature": False},
    {"label": "Only membrane locking alleviated", "ANS_membrane": True, "ANS_shear": False, "ANS_curvature": False},
    {"label": "Only transverse shear locking alleviated", "ANS_membrane": False, "ANS_shear": True, "ANS_curvature": False},
    {"label": "Only curvature-thickness locking alleviated", "ANS_membrane": False, "ANS_shear": False, "ANS_curvature": True},
    {"label": "All locking modes alleviated", "ANS_membrane": True, "ANS_shear": True, "ANS_curvature": True},
]

# ==============================
#       RUN SIMULATIONS
# ==============================

tip_displacements_all = []

for case in cases:
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

    tip_displacements_all.append((case["label"], tip_displacements))

# ==============================
#       PLOT RESULTS
# ==============================

plt.figure(figsize=(7, 5))

# Match styles by label string
markers = {
    "No locking alleviated": "d",
    "Only membrane locking alleviated": "o",
    "Only transverse shear locking alleviated": "s",
    "Only curvature-thickness locking alleviated": "^",
    "All locking modes alleviated": "X"
}

linestyles = {
    "No locking alleviated": ":",
    "Only membrane locking alleviated": "-",
    "Only transverse shear locking alleviated": "--",
    "Only curvature-thickness locking alleviated": "-.",
    "All locking modes alleviated": (0, (3, 1, 1, 1))
}

for label, tips in tip_displacements_all:
    disp_array = np.array(tips)
    x_disp = np.abs(disp_array[:, 0])
    x_disp = np.insert(x_disp, 0, 0.0)
    x_disp_normalized = x_disp / length
    n_steps = len(x_disp_normalized)
    k_normalized = np.linspace(0, 0.5, n_steps)

    plt.plot(
        k_normalized, x_disp_normalized,
        marker=markers[label],
        linestyle=linestyles[label],
        linewidth=2,
        label=label
    )

# Analytical solution from Euler's elastica
analytical = np.array([
    0, 0.016368357, 0.064510716, 0.141606309, 0.243173271,
    0.363380228, 0.495448848, 0.632116989, 0.766127679,
    0.890707595, 1.0
])
plt.plot(k_normalized, analytical, linestyle='-', color='black',
         linewidth=2, label='Analytical solution')

plt.xlabel(r"$k = \dfrac{ML}{2\pi EI}$", fontsize=12)
plt.ylabel(r"Normalized tip displacement $u_x / L$", fontsize=12)
plt.title("Influence of Individual Locking Modes", fontsize=14)
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.savefig(figures_dir / "compare_ans_with_analytical.png", dpi=300)
plt.close()
