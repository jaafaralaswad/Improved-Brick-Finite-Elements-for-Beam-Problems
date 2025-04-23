![Python Version](https://img.shields.io/badge/python-3.12-blue)
![OS](https://img.shields.io/badge/os-ubuntu%20%7C%20macos%20%7C%20windows-blue)
![License](https://img.shields.io/badge/license-MIT-green)


# Improved Brick Finite Elements for Beam Problems

## Table of Contents

- [Introduction](#introduction)
- [Conda Environment, Installation, and Testing](#conda-environment-installation-and-testing)
- [The Boundary Value Problem](#BVP)
- [Locking Modes](#locking-modes)
- [Assumed Natural Strain Method](#ANS)
- [Results](#results)

## Introduction


## Conda Environment, Install, and Testing


## The Boundary Value Problem
This code solves a 3D finite deformation elasticity problem using **Saint-Venant Kirchhoff hyperelasticity** in the **Lagrangian** description, employing **convective coordinates**.

This code solves a 3D nonlinear finite deformation elasticity problem in the Lagrangian framework using **Saint-Venant Kirchhoff hyperelasticity**. All quantities are formulated in the **reference configuration** using **convective coordinates** for better mechanical interpretation and implementation flexibility.

### 1. Kinematics

- We distinguish between the **reference configuration** $ \mathcal{B}_0 $ with material coordinates $ \mathbf{X} $, and the **current configuration** $ \mathcal{B} $ with spatial coordinates $ \mathbf{x} $. The displacement is:
  $$
  \mathbf{u} = \mathbf{x} - \mathbf{X}
  $$

- The **deformation gradient** is computed via the chain rule:
  \[
  \mathbf{F} = \frac{\partial \mathbf{x}}{\partial \mathbf{X}} = \frac{\partial \mathbf{x}}{\partial \xi^i} \otimes \frac{\partial \xi^i}{\partial \mathbf{X}} = \mathbf{g}_i \otimes \mathbf{G}^i
  \]
  where:
  - \( \xi^i \in \{\xi, \eta, \zeta\} \) are **convective (natural) coordinates**
  - \( \mathbf{g}_i \) are **covariant base vectors** in the current configuration
  - \( \mathbf{G}^i \) are **contravariant base vectors** in the reference configuration

- The **Green-Lagrange strain tensor** is:
  \[
  \mathbf{E} = \frac{1}{2}(\mathbf{C} - \mathbf{I}), \quad \text{with } \mathbf{C} = \mathbf{F}^\top \mathbf{F}
  \]

- In **contravariant basis** \( \{ \mathbf{G}^i \} \), the strain tensor is:
  \[
  \mathbf{E} = E^{ij} \, \mathbf{G}_i \otimes \mathbf{G}_j, \quad \text{with } E^{ij} = \frac{1}{2}(g^{ij} - G^{ij})
  \]

> Metric tensors \( g^{ij} = \mathbf{g}_i \cdot \mathbf{g}_j \), and \( G^{ij} = \mathbf{G}_i \cdot \mathbf{G}_j \) allow transformation between covariant and contravariant bases.

### 2. Balance Relations

- In the reference configuration, the **strong form** of the equilibrium equations is:
  \[
  \text{Div}(\mathbf{P}) + \rho_0 (\mathbf{b} - \ddot{\mathbf{x}}) = 0 \quad \text{in } \mathcal{B}_0
  \]
  where:
  - \( \mathbf{P} \) is the **first Piola-Kirchhoff stress**
  - \( \mathbf{b} \): body force per unit reference volume
  - \( \rho_0 \): reference density

> For static problems, inertial effects are neglected: \( \ddot{\mathbf{x}} = 0 \).

### 3. Constitutive Model: Saint-Venant Kirchhoff Material

- The Helmholtz free energy density is:
  \[
  \psi(C) = \frac{1}{8} \lambda (I_C - 3)^2 + \frac{1}{4} \mu (I_C^2 - 2I_C - 2I_{II} + 3)
  \]
  where \( I_C = \text{tr}(\mathbf{C}) \), \( I_{II} = \text{tr}(\mathbf{C}^2) \)

- The **second Piola-Kirchhoff stress** is:
  \[
  \mathbf{S} = \lambda \, \text{tr}(\mathbf{E}) \, \mathbf{I} + 2\mu \mathbf{E}
  \]

- In **contravariant basis**, the components are:
  \[
  S^{ij} = \lambda \, \text{tr}(\mathbf{E}) G^{ij} + 2\mu \, E^{ij}
  \]

- The **elasticity tensor** is:
  \[
  \mathbb{C}^{ijkl} = \lambda G^{ij} G^{kl} + \mu (G^{ik} G^{jl} + G^{il} G^{jk})
  \]

> This model is linear in strain, and suitable for small-to-moderate strains. All tensors are projected in the contravariant basis \( \{ \mathbf{G}^i \otimes \mathbf{G}^j \} \).

### 4. Boundary Conditions

The problem is completed with:

- **Displacement boundary conditions** on \( \Gamma_u \):
  \[
  \mathbf{u} = \bar{\mathbf{u}} \quad \text{on } \Gamma_u
  \]

- **Traction boundary conditions** on \( \Gamma_t \):
  \[
  \mathbf{P} \mathbf{N} = \bar{\mathbf{t}} \quad \text{on } \Gamma_t
  \]

  ## Locking modes

To be written.


## Assumed Natural Strain Method

To be written.


## More Information

To be written
