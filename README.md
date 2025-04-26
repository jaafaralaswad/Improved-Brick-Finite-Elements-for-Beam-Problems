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
- [Results](#results)

## Introduction

This project focuses on modeling beam structures using 3D brick finite elements, an approach that avoids many limitations of classical 1D beam theories. Traditional 1D beam theories have difficulty handling **advanced material behaviors**, cannot represent **cross-sectional deformations**, and pose challenges for **multiphysics coupling** and **contact interactions**. In contrast, so-called **solid-beam** formulations based on 3D brick elements are capable of surpassing these limitations.

However, low-order brick elements often suffer from **locking**, leading to inaccurate results. This project develops brick finite elements beams, combining Lagrange shape functions with the **Assumed Natural Strain (ANS)** method to alleviate locking.



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

### 3. Constitutive Model: Saint-Venant Kirchhoff Material

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

## Locking modes

To be written.


## Assumed Natural Strain Method

Here, we adopt the **Assumed Natural Strain (ANS)** method as outlined in:  

*Caseiro, J.F., Valente, R.F., Reali, A., Kiendl, J., Auricchio, F., & Alves de Sousa, R.*  
["On the Assumed Natural Strain method to alleviate locking in solid-shell NURBS-based finite elements."](https://doi.org/10.1007/s00466-014-0978-4 ) *Computational Mechanics*, **53**, 1341–1353 (2014).

However, we **adapt the formulation for beam problems**, **replace NURBS with Lagrange polynomials**, and **extend it to geometrically nonlinear analyses**.


## More Information

To be written
