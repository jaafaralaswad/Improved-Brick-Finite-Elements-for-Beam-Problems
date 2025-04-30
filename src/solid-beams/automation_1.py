import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from pre_processing import validate_inputs, create_mesh, impose_supports, compute_lame_parameters
from processing import newton_raphson_solver
from post_processing import compute_tip_displacement
from processing import shape_functions_3D

# ==============================
#       SETUP FIGURE DIRECTORY
# ==============================

figures_dir = Path(__file__).resolve().parents[2] / "generated-figures"
figures_dir.mkdir(exist_ok=True)

# ==============================
#       FIXED PARAMETERS
# ==============================

width = 1
height = 1
E = 1.2e7
nu = 0.0
numel = 10
ne_L = 2
ngp_c = 2
ngp_l = ne_L
n_load_steps = 10
max_iterations = 20
tolerance = 1e-15
prescribed_nodes = [0, 1, 2, 3]

# ==============================
#       SIMULATION CONFIGS
# ==============================

configs = [
    {"label": "No ANS",         "length": 20, "ANS_membrane": False, "ANS_shear": False, "ANS_curvature": False},
    {"label": "ANS Membrane",   "length": 20, "ANS_membrane": True,  "ANS_shear": False, "ANS_curvature": False},
    {"label": "ANS Shear",      "length": 20, "ANS_membrane": False, "ANS_shear": True,  "ANS_curvature": False},
    {"label": "ANS Curvature",  "length": 20, "ANS_membrane": False, "ANS_shear": False, "ANS_curvature": True},
    {"label": "All ANS",        "length": 20, "ANS_membrane": True,  "ANS_shear": True,  "ANS_curvature": True},
]

# ==============================
#       RUN ALL SIMULATIONS
# ==============================

displacement_curves = []

for config in configs:
    length = config["length"]

    # Validate input parameters
    validate_inputs(width, height, length, E, nu, numel, ne_L, ngp_c, ngp_l,
                    n_load_steps, max_iterations, tolerance)

    # Mesh and DOFs
    coords, nnode_W, nnode_H, nnode_L, ndof, connect, ne = create_mesh(width, height, length, numel, ne_L)
    prescribed_dofs = impose_supports(prescribed_nodes)

    # Material
    lambda_, mu = compute_lame_parameters(E, nu)

    # Solve
    u, displacements_all = newton_raphson_solver(
        E, width, height, length, lambda_, mu,
        coords, connect, prescribed_dofs,
        numel, ne, ne_L, ngp_c, ngp_l,
        n_load_steps, max_iterations,
        tolerance,
        config["ANS_membrane"],
        config["ANS_shear"],
        config["ANS_curvature"]
    )

    # Compute and store X-displacement at tip
    tip_displacements = [
        compute_tip_displacement(u_i, coords, connect, shape_functions_3D, ne_L)
        for u_i in displacements_all
    ]
    disp_x = np.array([d[0] for d in tip_displacements])
    displacement_curves.append((config["label"], disp_x))

# ==============================
#       ANALYTICAL SOLUTION
# ==============================

k_analytical = np.array([
    0.00, 0.05, 0.10, 0.15, 0.20,
    0.25, 0.30, 0.35, 0.40, 0.45,
    0.50, 0.55, 0.60, 0.65, 0.70,
    0.75, 0.80, 0.85, 0.90, 0.95,
    1.00
])/2

uL_analytical = np.array([
    0.000000000,
    0.004107265,
    0.016368357,
    0.036602238,
    0.064510716,
    0.099683684,
    0.141606309,
    0.189668042,
    0.243173271,
    0.301353415,
    0.363380228,
    0.428380067,
    0.495448848,
    0.563667407,
    0.632116989,
    0.699894561,
    0.766127679,
    0.829988630,
    0.890707595,
    0.947584593,
    1.000000000
])

# ==============================
#       PLOT COMPARISON
# ==============================

plt.figure(figsize=(6.5, 4.5))

# Adjusted styles for 5 simulations
offsets = [-0.012, -0.006, 0.0, 0.006, 0.012]
linestyles = ['-', '--', '-.', ':', (0, (1, 1))]
markers = ['o', 's', '^', 'D', 'v']
colors = ['C0', 'C1', 'C2', 'C3', 'C4']

for i, (label, disp_x) in enumerate(displacement_curves):
    u_norm = np.abs(disp_x) / length
    k = np.linspace(0, 0.5, len(u_norm)) + offsets[i]
    plt.plot(k, u_norm, marker=markers[i], linestyle=linestyles[i],
             color=colors[i], label=label)

# Analytical solution
plt.plot(k_analytical, uL_analytical, color='black', linewidth=2.5, label='Analytical', zorder=10)

plt.xlabel(r"$k = \dfrac{ML}{2\pi EI}$", fontsize=12)
plt.ylabel(r"$\dfrac{u}{L}$", fontsize=12)
plt.title("Normalized Tip Displacement Comparison", fontsize=12)
plt.xlim(0, 0.5)
plt.ylim(0, 1.2)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(figures_dir / "tip_displacement_comparison_normalized.png", dpi=300)
plt.show()
