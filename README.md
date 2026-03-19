
# Consensus-based Distributed Optimization for Multi-agent Systems over Multiplex Networks

<p align="center">
  <a href="https://ieeexplore.ieee.org/abstract/document/10772699">
    <img src="https://img.shields.io/badge/IEEE-10772699-blue?style=flat-square&logo=ieee" alt="IEEE Paper"/>
  </a>
  <a href="https://arxiv.org/abs/2304.01875">
    <img src="https://img.shields.io/badge/arXiv-2304.01875-b31b1b?style=flat-square&logo=arxiv" alt="arXiv"/>
  </a>
  <img src="https://img.shields.io/badge/Julia-1.8%2B-9558B2?style=flat-square&logo=julia" alt="Julia"/>
  <img src="https://img.shields.io/badge/License-MIT-green?style=flat-square" alt="License"/>
  <img src="https://img.shields.io/badge/Status-Published-brightgreen?style=flat-square" alt="Status"/>
</p>


---

## 📄 Paper

**Title:** Consensus-based Distributed Optimization for Multi-agent Systems over Multiplex Networks

**Authors:** Christian D. Rodríguez-Camargo, Andres F. Urquijo-Rodríguez, Eduardo Mojica-Nava

**Journal:** IEEE Transactions on Network Science and Engineering (2024)

**DOI:** [10.1109/TNSE.2024.10772699](https://doi.org/10.1109/TNSE.2024.10772699)

**Preprint:** [arXiv:2304.01875](https://arxiv.org/abs/2304.01875)

> **Abstract:** Multilayer networks provide a more comprehensive framework for exploring real-world and engineering systems than traditional single-layer networks. This paper proposes two algorithms for distributed optimization problems in multiplex networks using the supra-Laplacian matrix and its diffusion dynamics: a distributed saddle-point algorithm and a distributed gradient descent algorithm. By relating consensus and diffusion dynamics, we obtain the multiplex supra-Laplacian matrix and extend classical distributed optimization to the multiplex setting. Convergence is analyzed with theoretical guarantees. Numerical examples validate both algorithms and explore the impact of interlayer diffusion on consensus time. A coordinated dispatch application for interdependent energy–gas infrastructure networks is also presented.

---

## 🗂️ Repository Structure

```
.
├── README.md
├── LICENSE
├── Project.toml                   # Julia package dependencies
├── Manifest.toml
├── consensus_bilayer.jl           # Bilayer multiplex: saddle-point & gradient descent (Figs. 3–5)
└── consensus_three_layer.jl       # Three-layer energy–gas dispatch (Fig. 7)
```

---

## 🔬 Overview

This repository contains the Julia code used to produce the numerical results in the paper. The implementation solves two distributed continuous-time optimization algorithms for **multiplex networks** — a class of multilayer networks where the same set of nodes is connected across multiple types of interaction layers.

### Algorithms implemented

**1. Distributed saddle-point dynamics** (`consensus_bilayer.jl` → `integrate_system`)

Solves the augmented supra-Lagrangian saddle-point problem. The second-order ODE system implemented is:

$$\dot{y} = v, \qquad \dot{v} = -(\mathcal{L} + \mathbb{I}) v - \mathcal{L}^2 y$$

where $\mathcal{L}$ is the supra-Laplacian matrix of the multiplex network. Integration uses a fixed-step **RK4** solver with a discrete termination callback triggered when $\sum_{i < j}|y_i - y_j| \leq 10^{-2}$.

**2. Distributed gradient descent** (`consensus_bilayer.jl` → `integrate_system_gd`)

A soft-penalty-based variant with a time-decaying gain $\varsigma(t) = \frac{1}{t + \theta}$ (default $\theta = 10000$):

$$\dot{y} = -(\mathcal{L} + \varsigma(t) \mathbb{I}) y - \varsigma(t) \nabla\tilde{f}(y)$$

Convergence is terminated when $\sum_{i < j}|y_i - y_j| \leq 10^{-3}$.

**3. Coordinated dispatch for energy–gas multiplex network** (`consensus_three_layer.jl`)

Three-layer multiplex network (fully connected / mixed / ring topologies) modeling a multi-energy system. The RHS includes the power and gas balance constraints:

$$\dot{y} = v, \qquad \dot{v} = -(\mathcal{L} + \mathbb{I}) v - \mathcal{L}^2 q - \alpha(\mathbf{q}) \tilde{\mathbf{1}}$$

where $\alpha(\mathbf{q})$ encodes the active power and inelastic gas demand residuals via the fuel conversion factor $\phi = 0.7$. Integration uses **AutoTsit5(Rosenbrock23())**, an adaptive stiff/non-stiff solver.

---

## ⚙️ Installation

### Prerequisites

- [Julia](https://julialang.org/downloads/) ≥ 1.8
- A working Python environment with `matplotlib` installed (required by `PyPlot` / `PyCall`)
- (Optional) A LaTeX distribution for rendered axis labels (`text.usetex = true` in `consensus_three_layer.jl`)

### Setup

Clone the repository and instantiate the Julia environment:

```bash
git clone https://github.com/alfagalileo/Consensus-based-Distributed-Optimization.git
cd Consensus-based-Distributed-Optimization
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

### Dependencies

All dependencies are captured in `Project.toml`. The main packages are:

| Package | Role |
|---|---|
| `OrdinaryDiffEq` | ODE solvers (RK4, AutoTsit5, Rosenbrock23) |
| `DiffEqCallbacks` | Discrete termination callbacks |
| `Graphs` | Graph construction and Laplacian matrix assembly |
| `BlockArrays` | Supra-Laplacian block matrix construction |
| `LinearAlgebra`, `SparseArrays` | Dense and sparse matrix operations |
| `PyPlot`, `PyCall` | Plotting via Matplotlib |
| `IterTools` | Combinatorial subset iteration |
| `JSON`, `DelimitedFiles` | Data export |


