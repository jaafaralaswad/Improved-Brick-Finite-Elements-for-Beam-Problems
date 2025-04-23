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

### 1. Kinematics

We distinguish between the reference configuration (`$\mathcal{B}_0$`) with material coordinates (`$\mathbf{X}$`), and the current configuration (`$\mathcal{B}$`) with spatial coordinates (`$\mathbf{x}$`). The displacement is:

```
u = x - X
```

The deformation gradient is computed via the chain rule:

```
F = ∂x/∂X = ∂x/∂ξᵢ ⊗ ∂ξᵢ/∂X = gᵢ ⊗ Gⁱ
```

where:

- `$ξ^i ∈ \{ξ, η, ζ\}$` are **convective (natural) coordinates**
- `$g_i$` are **covariant base vectors** in the current configuration
- `$G^i$` are **contravariant base vectors** in the reference configuration

The Green-Lagrange strain tensor is:

```
E = 1/2 (C - I),     where C = Fᵀ F
```

In contravariant basis `{G^i}`, the strain tensor is:

```
E = Eⁱʲ Gᵢ ⊗ Gⱼ,     with   Eⁱʲ = 1/2 (gⁱʲ - Gⁱʲ)
```

Metric tensors are:

```
gⁱʲ = gᵢ · gⱼ,    Gⁱʲ = Gᵢ · Gⱼ
```

These allow for transformation between covariant and contravariant components.

---

### 2. Balance Relations

In the reference configuration, the strong form of the equilibrium equations is:

```
Div(P) + ρ₀ (b - ẍ) = 0    in  𝔅₀
```

where:

- `$P$` is the **first Piola-Kirchhoff stress**
- `$b$`: body force per unit reference volume
- `$ρ₀$`: reference mass density

For static problems, inertial effects are neglected: `$ẍ = 0$`.

---

### 3. Constitutive Model: Saint-Venant Kirchhoff Material

The Helmholtz free energy density is:

```
ψ(C) = (1/8) λ (I_C - 3)² + (1/4) μ (I_C² - 2 I_C - 2 I_II + 3)
```

where:

- `$I_C = tr(C)$`
- `$I_{II} = tr(C^2)$`

The **second Piola-Kirchhoff stress** is:

```
S = λ tr(E) I + 2μ E
```

In contravariant basis, the components are:

```
Sⁱʲ = λ tr(E) Gⁱʲ + 2μ Eⁱʲ
```

The elasticity tensor is:

```
Cⁱʲᵏˡ = λ Gⁱʲ Gᵏˡ + μ (Gⁱᵏ Gʲˡ + Gⁱˡ Gʲᵏ)
```

This model is linear in strain and suitable for small-to-moderate deformations. All stress and elasticity tensors are projected in the contravariant basis `{Gⁱ ⊗ Gʲ}`.

---

### 4. Boundary Conditions

The boundary value problem is completed with:

- **Displacement conditions** on `$\Gamma_u$`:
  
  ```
  u = ū     on  Γ_u
  ```

- **Traction conditions** on `$\Gamma_t$`:

  ```
  P N = t̄     on  Γ_t
  ```


  ## Locking modes

To be written.


## Assumed Natural Strain Method

To be written.


## More Information

To be written
