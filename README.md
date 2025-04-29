![Python Version](https://img.shields.io/badge/python-3.12-blue)
![OS](https://img.shields.io/badge/os-ubuntu%20%7C%20macos%20%7C%20windows-blue)
![License](https://img.shields.io/badge/license-MIT-green)


# Improved Brick Finite Elements for Beam Problems

## Table of Contents

- [Introduction](#introduction)
- [Conda Environment, Installation, and Testing](#conda-environment-installation-and-testing)
- [The Boundary Value Problem](#the-boundary-value-problem)
- [Nonlinear Finite Element Procedure](#nonlinear-finite-element-procedure)
- [Locking Modes](#locking-modes)
- [Assumed Natural Strain Method](#assumed-natural-strain-method)
- [Problem Setup and Usage Instructions](#problem-setup-and-usage-instructions)
- [Results](#results)
- [More Information](#more-information)

## Motivation

It is common to model slender structures using classical **beam finite elements**, which are efficient and widely used in engineering applications.

<p align="center">
  <img src="README_figures/beam_elements.jpg" alt="Beam elements" width="400"/>
</p>

However, **1D beam theories** come with important limitations: they struggle to capture **advanced material behaviors**, cannot represent **cross-sectional deformations**, and pose challenges for **multiphysics coupling** and **contact interactions**. These limitations stem from the **inherent assumptions of 1D beam theory** and are **not artifacts of the finite element method**.


To overcome these limitations, this project adopts a **solid-beam modeling approach** using **3D brick finite elements**. These formulations offer greater flexibility and accuracy in capturing the full 3D deformation of beam-like structures.

<p align="center">
  <img src="README_figures/brick_elements.jpg" alt="Solid-beam elements" width="400"/>
</p>


However, **low-order brick elements** often suffer from **locking**, leading to inaccurate results. This project develops **brick finite element beams**, combining **Lagrange shape functions** with the **Assumed Natural Strain (ANS)** method to alleviate locking. We apply the code to a **benchmark problem** and investigate the **effects of different locking modes**. We also compare **ANS** to other approaches like **h-refinement**, **p-refinement**, and **reduced integration**.

Below are two videos comparing the analytical and finite element solutions for a cantilever beam subjected to a tip moment.

The **first video** shows the **analytical solution**, where the beam forms a full circle under a tip moment of $M = 2π·EI / L$.

The **second video** presents the **finite element simulation** using a coarse 3D mesh with four fully integrated trilinear elements. Due to **locking**, the beam fails to deform into a full circle and appears overly stiff.

Click on each thumbnail below to watch the video on **YouTube**:

[![Cantilever Beam Forming a Full Circle – Analytical Solution](https://img.youtube.com/vi/4MB-QnyLL30/0.jpg)](https://www.youtube.com/watch?v=4MB-QnyLL30)

[![Beam Fails to Curve Fully – FEM with Locking (4 Elements, Full Integration)](https://img.youtube.com/vi/REFapxlKZWc/0.jpg)](https://www.youtube.com/watch?v=REFapxlKZWc)


## Conda Environment, Install, and Testing


## The Boundary Value Problem

In order to construct a complete BVP, we need to bring together kinematics, balance relations and material law, in addition to the boundary conditions. We formulate the problem using **convective coordinates** for better structural mechanical interpretation and later implementation flexibility.

### 1. Kinematics

We distinguish between the **reference configuration** $\mathcal{B}_0$ with material coordinates $\mathbf{X}$, and the **current configuration** $\mathcal{B}$ with spatial coordinates $\mathbf{x}$. The displacement is given by $\mathbf{u} = \mathbf{x} - \mathbf{X}$.

- The **deformation gradient** is computed via the chain rule:  
  $\mathbf{F} = \dfrac{\partial \mathbf{x}}{\partial \mathbf{X}} = \dfrac{\partial \mathbf{x}}{\partial \xi^i} \otimes \dfrac{\partial \xi^i}{\partial \mathbf{X}} = \mathbf{g}_i \otimes \mathbf{G}^i$  
  where:  
  - $\xi^i \in \{\xi, \eta, \zeta\}$ are **convective (natural) coordinates**  
  - $\mathbf{g}_i$ are **covariant base vectors** in the current configuration  
  - $\mathbf{G}^i$ are **contravariant base vectors** in the reference configuration

- The **Green-Lagrange strain tensor** is given by  
  $\mathbf{E} = \dfrac{1}{2}(\mathbf{C} - \mathbf{I})$, with $\mathbf{C} = \mathbf{F}^\top \mathbf{F}$.

- In the **contravariant basis** $\{\mathbf{G}^i\}$, the strain tensor is written as  
  $\mathbf{E} = E^{ij} \, \mathbf{G}_i \otimes \mathbf{G}_j$, with $E^{ij} = \dfrac{1}{2}(g^{ij} - G^{ij})$.

> Metric tensors $g^{ij} = \mathbf{g}_i \cdot \mathbf{g}_j$ and $G^{ij} = \mathbf{G}_i \cdot \mathbf{G}_j$ allow transformation between covariant and contravariant components.

- In **Voigt notation**, we exploit the symmetry of the strain tensor to write it as  
$\hat{\mathbf{E}} = \left[ \overline{E}{\xi \xi},\ \overline{E}{\eta \eta},\ \overline{E}{\zeta \zeta},\ 2\overline{E}{\xi \eta},\ 2\overline{E}{\eta \zeta},\ 2\overline{E}{\xi \zeta} \right]^\mathrm{T}$.



### 2. Balance Relations

- The **strong form** of the equilibrium equations is:  
  $\text{Div}(\mathbf{P}) + \rho_0 (\mathbf{b} - \ddot{\mathbf{x}}) = 0 \quad \text{in } \mathcal{B}_0$  
  where:  
  - $\mathbf{P}$ is the **first Piola-Kirchhoff stress**  
  - $\mathbf{b}$ is the **body force** per unit reference volume  
  - $\rho_0$ is the **reference density**

> For static problems, inertial effects are neglected: $\ddot{\mathbf{x}} = 0$.

### 3. Material Law: Saint-Venant Kirchhoff Material

- The **Helmholtz free energy** is defined in the **reference configuration** as  
  $\psi_0(C) = \dfrac{1}{8} \lambda (I_C - 3)^2 + \dfrac{1}{4} \mu (I_C^2 - 2I_C - 2II_{C} + 3)$  
  where $I_C = \text{tr}(\mathbf{C})$, $II_{C} = \text{tr}(\mathbf{C}^2)$.

- The **second Piola-Kirchhoff stress** can be calculated as  
  $\mathbf{S} = \lambda \text{tr}(\mathbf{E}) \mathbf{I} + 2\mu \mathbf{E}$.

- The 4th order **elasticity tensor** in the **contravariant basis** is given by  
  $\mathbb{C}^{ijkl} = \lambda G^{ij} G^{kl} + \mu (G^{ik} G^{jl} + G^{il} G^{jk})$.

> This model is linear in strain and suitable for small-to-moderate strains.

### 4. Boundary Conditions

The problem is completed with:

- **Displacement boundary conditions** on $\Gamma_u$:  
  $\mathbf{u} = \bar{\mathbf{u}} \quad \text{on } \Gamma_u$

- **Traction boundary conditions** on $\Gamma_t$:  
  $\mathbf{P} \mathbf{N} = \bar{\mathbf{t}} \quad \text{on } \Gamma_t$

## Nonlinear Finite Element Procedure

We are working in **convective curvilinear coordinates**, so while we follow the general framework from *Nonlinear Finite Element Methods* by **Wriggers (2008)**, some steps are **modified** to suit this coordinate system. Below, we summarize only the **key results** and **necessary adjustments**.

The weak form for the problem neglecting inertial effects and body forces is given by

$\int_{\Omega_0} \mathbf{S} : \delta \mathbf{E} \mathrm{d}V - \int_{\partial \Omega_0} \mathbf{t}_0 \cdot \delta \mathbf{x} \mathrm{d}A = 0$

  where:  
  - $\mathbf{S}$ is **2nd Piola-Kirchhoff stress tensor**.
  - $\mathbf{t}_0$ is **surface traction** in the reference configuration.

One can write:

$g(\mathbf{x}, \delta \mathbf{x}) = g^{\mathrm{int}}(\mathbf{x}, \delta \mathbf{x}) - g^{\mathrm{ext}}(\mathbf{x}, \delta \mathbf{x}) = 0$

where:
  - $g^{\mathrm{int}}(\mathbf{x}, \delta \mathbf{x}) := \int_{\Omega_0} \mathbf{S} : \delta \mathbf{E} \mathrm{d}V$ is the **internal virtual work**.
  - $g^{\mathrm{ext}}(\mathbf{x}, \delta \mathbf{x}) := \int_{\partial \Omega_0} \mathbf{t}_0 \cdot \delta \mathbf{x} \mathrm{d}A$ is the **external virtual work**.

The **lienarization** of the **internal virtual work** is given by

$\Delta g^{\mathrm{int}} = \int_{\Omega_0} \Delta \mathbf{E} : \mathbb{C}(\mathbf{E}) : \delta \mathbf{E} \mathrm{d}V + \int_{\Omega_0} \mathbf{S} : \Delta \delta \mathbf{E} \mathrm{d}V$

where:
- $\int_{\Omega_0} \Delta \mathbf{E} : \mathbb{C}(\mathbf{E}) : \delta \mathbf{E} \mathrm{d}V$ is the material part.
- $\int_{\Omega_0} \mathbf{S} : \Delta \delta \mathbf{E} \mathrm{d}V$ is the geometric part.

The **linearization** of the **external virtual work** is zero when the applied load is **conservative**. However, if the load is **non-conservative** (e.g., a **follower load**), the linearization is generally **non-zero** and must be taken account for. This leads to the following consequences:

- The **weak form** cannot be derived from a strain energy **functional**, since a **non-conservative load** has no potential.
- The **load must be updated** during every **Newton-Raphson iteration** to stay consistent.
- An additional (third) **contribution to the tangent matrix** must be added — but **only for elements directly loaded**. Because of this, the tangent matrix is **no longer symmetric**.

For details about the **linearization of the external virtual work**, see *Wriggers (2008), Section 4.2.5*. We follow the procedure **exactly as outlined** there.

In this project, we model **bending moment** acting at the tip of the beam as a **follower load**, where $\hat{p} = -\left(\dfrac{12M}{I}\right)\zeta, \quad \zeta \in \left[-\dfrac{h}{2}, \dfrac{h}{2}\right]$. We **apply** it as a **first Piola-Kirchhoff stress tensor**, meaning it is defined over the **undeformed area**. The load has a **constant magnitude**; only **direction changes** with deformation.


<p align="center">
  <img src="README_figures/bending-traction.jpg" alt="Bending traction" width="400"/>
</p>


We discretize the domain into $n_{\text{el}}$ finite elements.

Within the **isoparametric concept**, both the geometry and displacements are approximated using the same **shape functions** $N_I(\xi, \eta, \zeta)$:

$$
\mathbf{u}(\xi, \eta, \zeta) \approx \sum_{I=1}^{n_e} N_I(\xi, \eta, \zeta) \mathbf{u}_I
$$

The mapping between parametric $(\xi, \eta, \zeta)$ and physical coordinates $(X_1, X_2, X_3)$ is given by the **Jacobian**:

$$ \begin{bmatrix}
\frac{\partial (\cdot)}{\partial \xi} \\
\frac{\partial (\cdot)}{\partial \eta} \\
\frac{\partial (\cdot)}{\partial \zeta}
\end{bmatrix} =
\mathbf{J}
\begin{bmatrix}
\frac{\partial (\cdot)}{\partial X_1} \\
\frac{\partial (\cdot)}{\partial X_2} \\
\frac{\partial (\cdot)}{\partial X_3}
\end{bmatrix}
$$

with:

$$
\mathbf{J} =
\begin{bmatrix}
\frac{\partial X_1}{\partial \xi} & \frac{\partial X_2}{\partial \xi} & \frac{\partial X_3}{\partial \xi} \\
\frac{\partial X_1}{\partial \eta} & \frac{\partial X_2}{\partial \eta} & \frac{\partial X_3}{\partial \eta} \\
\frac{\partial X_1}{\partial \zeta} & \frac{\partial X_2}{\partial \zeta} & \frac{\partial X_3}{\partial \zeta}
\end{bmatrix}
$$

Inverted:

$$
\begin{bmatrix}
\frac{\partial N_I}{\partial X_1} \\
\frac{\partial N_I}{\partial X_2} \\
\frac{\partial N_I}{\partial X_3}
\end{bmatrix}
= \mathbf{J}^{-1}
\begin{bmatrix}
\frac{\partial N_I}{\partial \xi} \\
\frac{\partial N_I}{\partial \eta} \\
\frac{\partial N_I}{\partial \zeta}
\end{bmatrix}
$$

The variation of the Green-Lagrange strain tensor $\delta \mathbf{E}$ in **Voigt notation**:

$$
\delta \widehat{\mathbf{E}} =
\begin{bmatrix}
\delta E_{\xi\xi} \\
\delta E_{\eta\eta} \\
\delta E_{\zeta\zeta} \\
2\delta E_{\xi\eta} \\
2\delta E_{\eta\zeta} \\
2\delta E_{\xi\zeta}
\end{bmatrix}
\approx
\sum_{I=1}^{n_e}
\mathbf{B}_I \, \delta \mathbf{u}_I
$$

where each $\mathbf{B}_I$ matrix is:

$$
\mathbf{B}_I = \begin{bmatrix}
\left[ N_{I,\xi} \, \mathbf{g}_\xi^T \right] \\
\left[ N_{I,\eta} \, \mathbf{g}_\eta^T \right] \\
\left[ N_{I,\zeta} \, \mathbf{g}_\zeta^T \right] \\
\left[ N_{I,\xi} \, \mathbf{g}_\eta^T + N_{I,\eta} \, \mathbf{g}_\xi^T \right] \\
\left[ N_{I,\eta} \, \mathbf{g}_\zeta^T + N_{I,\zeta} \, \mathbf{g}_\eta^T \right] \\
\left[ N_{I,\xi} \, \mathbf{g}_\zeta^T + N_{I,\zeta} \, \mathbf{g}_\xi^T \right]
\end{bmatrix}
$$

- $\mathbf{g}_\xi$, $\mathbf{g}_\eta$, $\mathbf{g}_\zeta$: covariant base vectors in current configuration.



## Locking modes




## Assumed Natural Strain Method

Here, we adopt the **Assumed Natural Strain (ANS)** method as outlined in:  

*Caseiro, J.F., Valente, R.F., Reali, A., Kiendl, J., Auricchio, F., & Alves de Sousa, R.*  
["On the Assumed Natural Strain method to alleviate locking in solid-shell NURBS-based finite elements."](https://doi.org/10.1007/s00466-014-0978-4 ) *Computational Mechanics*, **53**, 1341–1353 (2014).

However, we **adapt the formulation for beam problems**, **replace NURBS with Lagrange polynomials**, and **extend it to geometrically nonlinear analyses**.

## Problem Setup and Usage Instructions

The problem we solve here is a rectangular **cantiliever beam** subjected to **bending moment applied at its tip**.  

<p align="center">
  <img src="README_figures/cantilever.jpg" alt="Cantilever beam" width="400"/>
</p>


All the user needs to do is define the parameters in the `DEFINE PROBLEM SETUP` block from `main.py`.


<p align="center">
  <img src="README_figures/user-parameters.jpg" alt="User parameters" width="800"/>
</p>


The parameters are:

- **width**: Width of the beam (cross-sectional).
- **height**: Height of the beam (cross-sectional).
- **length**: Length of the beam (axial direction).
- **E**: Young’s modulus.
- **nu**: Poisson’s ratio. Set to `0.0` to eliminate Poisson and volumetric locking.

- **numel**: Number of finite elements along the beam length. Increase this value to perform **h-refinement**. Each cross-section contains only a single element.
- **ne_L**: Number of nodes per element along the beam axis.  
  - `ne_L = 2` → linear shape functions in the axial direction  
  - `ne_L = 3` → quadratic shape functions, and so on  
  Increase this value to perform **p-refinement**. Shape functions in the cross-sectional directions are always linear.

- **ngp_c**: Number of Gauss points in each cross-sectional direction. Use `ngp_c = 2` for full integration. Setting `ngp_c = 1` often produces inaccurate results.
- **ngp_l**: Number of Gauss points in the length direction.  
- `ngp_l = ne_L` → **full integration**  
- `ngp_l = ne_L - 1` → **reduced integration** in the **axial direction**  
**Note:** No **stabilization techniques** are included in the current implementation.

- **ANS_membrane**: Enables/disables the Assumed Natural Strain (ANS) method to alleviate **membrane locking**.
- **ANS_shear**: Enables/disables ANS to alleviate **transverse shear locking**.
- **ANS_curvature**: Enables/disables ANS to alleviate **curvature-thickness locking**.

- **n_load_steps**: Number of load steps to incrementally apply the external load.
- **max_iterations**: Maximum number of Newton–Raphson iterations allowed per load step.
- **tolerance**: Convergence tolerance on the **energy norm**.

Other modifications, such as changing the load and boundary conditions, must be made manually within their respective functions.

## Results


## More Information

The current code serves as a platform for further development. Additional validation can be performed by solving more benchmark problems and conducting more thorough comparisons between different locking alleviation methods. **Poisson locking** and **volumetric locking** are suppressed by using a **linear material model** and setting the **Poisson's ratio** to zero.   The formulation can be augmented with methods such as the **Enhanced Assumed Strain (EAS)** technique to properly address these effects when present.

The code is already structured to support these **extensions**. We had intended to develop these **implementations** within **FEniCSx** to establish a more **comprehensive** and more **usable** **framework** for beam problems, especially that, to the best of our knowledge, these important techniques are not yet available there, and in general, they lack beginner-friendly documentation, despite their significance for addressing locking phenomena in finite element simulations. However, **time limitations** have prevented us from pursuing this effort at present.
