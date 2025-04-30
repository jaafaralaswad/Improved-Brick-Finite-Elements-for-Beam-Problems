import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

from pre_processing import validate_inputs, create_mesh, impose_supports, compute_lame_parameters
from processing import newton_raphson_solver, compute_external_force, element_3D_nonlinear, gauss, solve_linear_system, compute_C_SVK_component, shape_functions_3D, voigt_index
from post_processing import compute_tip_displacement, plot_tip_displacement_x
import visualize as viz

# ==============================
#       SETUP FIGURE DIRECTORY
# ==============================

# Set directory to save figures and plots
figures_dir = Path(__file__).resolve().parents[2] / "generated-figures"
figures_dir.mkdir(exist_ok=True)

# ==============================
#       DEFINE PROBLEM SETUP
# ==============================

# Geometry of the rectangular cantilever beam
width = 1.0     # Beam width (cross-sectional)
height = 0.5    # Beam height (cross-sectional)
length = 20.0   # Beam length (axial direction)

# Material properties (Saint Venant–Kirchhoff model)
E = 1.2e7        # Young's modulus
nu = 0.0        # Poisson's ratio (set to 0 to avoid Poisson and volumetric locking)

# Finite element discretization
numel = 10              # Number of elements along the beam length
ne_L = 2               # Number of nodes per element in axial direction (minimum = 2)

# Integration points
ngp_c = 2       # Number of Gauss points in each cross-sectional direction
ngp_l = ne_L    # Number of Gauss points in axial direction (≤ 10)

# Locking alleviation techniques (Assumed Natural Strain method)
ANS_membrane = False    # Alleviate membrane locking
ANS_shear = False       # Alleviate transverse shear locking
ANS_curvature = False   # Alleviate curvature-thickness locking

# Incremental loading and Newton-Raphson settings
n_load_steps = 10       # Number of load steps
max_iterations = 20     # Maximum Newton-Raphson iterations per load step
tolerance = 1e-15       # Convergence tolerance

# Visualization options
visualize_reference = True      # Visualize reference configuration
visualize_final = True          # Visualize final configuration
plot_displacement = True        # Plot horizontal tip displacement 

# ==============================
#           PRE-PROCESSING
# ==============================

# Pre-processing: validate inputs, generate mesh, apply boundary conditions, compute material properties, and visualize reference configuration

# ------------------------------
# Validate input parameters
# ------------------------------
validate_inputs(width, height, length, E, nu, numel, ne_L, ngp_c, ngp_l, n_load_steps, max_iterations, tolerance)

# ------------------------------
# Generate mesh (nodes, connectivity, DOFs)
# ------------------------------
coords, nnode_W, nnode_H, nnode_L, ndof, connect, ne = create_mesh(width, height, length, numel, ne_L)

# ------------------------------
# Apply displacement boundary conditions (supports)
# ------------------------------
prescribed_nodes = [0, 1, 2, 3]  # Nodes to be clamped
prescribed_dofs = impose_supports(prescribed_nodes) # Corresponding constrained DOFs

# ------------------------------
# Compute Lamé parameters (material properties)
# ------------------------------
lambda_, mu = compute_lame_parameters(E, nu) # Calculate Lame parameters

# ------------------------------
# Visualize reference (undeformed) configuration
# ------------------------------
if visualize_reference:
    img_fname = figures_dir / "reference_configuration.png"
    viz.plot_mesh_3D(str(img_fname), coords, connect, prescribed_nodes)


# ==============================
#             SOLVING THE SYSTEM
# ==============================

u, displacements_all = newton_raphson_solver(
    E, width, height, length, lambda_, mu,
    coords, connect, prescribed_dofs,
    numel, ne, ne_L, ngp_c, ngp_l,
    n_load_steps, max_iterations,
    tolerance, ANS_membrane,
    ANS_shear, ANS_curvature
)

# ==============================
#                POST-PROCESSING
# ==============================

# Interpolate the displacements at the center of the tip of the beam
if plot_displacement:
    # Tip displacements for every load step
    tip_displacements = [
        compute_tip_displacement(u_i, coords, connect, shape_functions_3D, ne_L)
        for u_i in displacements_all
    ]
    print(tip_displacements)
    # Plot X displacement curve and save it
    plot_tip_displacement_x(tip_displacements, length, save_path=figures_dir / "tip_displacement_x.png")


# ================================
# VISUALIZE DEFORMED CONFIGURATION
# ================================
if visualize_final:
    img_fname_def = figures_dir / "final_configuration.png"
    viz.plot_deformed_mesh_3D(str(img_fname_def), coords, connect, u, scale=1.0, prescribed_nodes=prescribed_nodes)