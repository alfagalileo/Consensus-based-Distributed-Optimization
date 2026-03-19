
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


---

## 🚀 Usage

### Bilayer network — saddle-point algorithm

Runs the saddle-point ODE on a 2-layer multiplex network (ring + complete graph topologies) and plots the consensus dynamics. Reproduces **Figures 3 and 4**.

```julia
include("consensus_bilayer.jl")

# Saddle-point: 5 nodes/layer, D^[1]=D^[2]=1.0, D^[1,2]=0.1, t_max=130
integrate_system(5, 1.0, 1.0, 0.1; to=130.0)

# Gradient descent on the same network
integrate_system_gd(5, 1.0, 1.0, 0.1)
```

Function signatures:

```julia
integrate_system(n, d1, d2, dx; to=5000.0)
# n    – number of nodes per layer
# d1   – intralayer diffusion constant D^[1] (layer 1: ring graph)
# d2   – intralayer diffusion constant D^[2] (layer 2: complete graph)
# dx   – interlayer diffusion constant D^[1,2]
# to   – maximum integration time

integrate_system_gd(n, d1, d2, dx)
# same arguments; uses time-varying gain ζ(t) = 1/(t + 10000)
```

### Consensus time sweep — critical phenomena

Sweeps $D_x \in [0.5, 1.5]$ over 1000 values for network sizes $N \in \{11, 15, 21, 31\}$ using multithreading. Reproduces **Figure 5**.

```julia
include("consensus_bilayer.jl")
main()   # results optionally exported to data.json (uncomment lines 170–172)
```

Run Julia with multiple threads to accelerate the sweep:

```bash
julia --project=. --threads=auto consensus_bilayer.jl
```

### Three-layer energy–gas dispatch

Runs the coordinated dispatch simulation on a 3-layer multiplex network of 7 generators per layer with parameters from the paper: $D^{[T]}=D^{[G]}=0.2$, $D^{[K]}=0.8$, $D_x=0.6$, $\phi=0.7$, $P_D = G_D = 100$. Reproduces **Figure 7**.

```julia
include("consensus_three_layer.jl")

# Single run (paper default parameters)
integrate_system(7, [0.2, 0.8, 0.2], 0.6; to=30.0)

# Sweep D_x for N=11 and export results to data_threelayers.json
main()
```

Function signature:

```julia
integrate_system(n, d, dx; to=5000.0)
# n   – number of nodes per layer
# d   – vector of intralayer diffusion constants [D^[T], D^[K], D^[G]]
# dx  – interlayer diffusion constant D^[α,β] (uniform across layer pairs)
# to  – maximum integration time
```

---

## 📊 Reproducing Paper Figures

| Figure | Description | File | Entry point |
|--------|-------------|------|-------------|
| Fig. 3 | Saddle-point consensus, 2-layer / 3-node network | `consensus_bilayer.jl` | `integrate_system(5, 1.0, 1.0, 0.1; to=130.0)` |
| Fig. 4 | Gradient descent with varying $\theta$ | `consensus_bilayer.jl` | `integrate_system_gd(...)` with different `ζ(t)` |
| Fig. 5 | Consensus time $t_c$ vs. $D_x$ — critical phenomena | `consensus_bilayer.jl` | `main()` |
| Fig. 7 | Energy–gas dispatch consensus (3 layers, 7 nodes) | `consensus_three_layer.jl` | `integrate_system(7, [0.2,0.8,0.2], 0.6; to=30.0)` |

---

## 📐 Key Parameters

| Symbol | Variable in code | Description |
|--------|-----------------|-------------|
| $D^{[\alpha]}$ | `d1`, `d2`, `d[i]` | Intralayer diffusion constants per layer |
| $D^{[\alpha,\beta]}$ | `dx` | Interlayer diffusion constant |
| $N$ | `n` | Number of nodes per layer |
| $\varsigma(t)$ | `ζ(t)` | Time-decaying gain for gradient descent: $1/(t+\theta)$ |
| $\phi$ | `ϕ` | Fuel conversion factor (energy–gas dispatch), default `0.7` |
| $P_D$, $G_D$ | `Pd`, `Gd` | Electrical power and inelastic gas demand, default `100.0` |
| $\epsilon$ | convergence threshold | $\sum_{i<j}\|y_i - y_j\| \leq 10^{-2}$ (saddle-point), $10^{-3}$ (GD, dispatch) |

---

## 📝 Citation

If you use this code in your research, please cite:

```bibtex
@article{rodriguez2024consensus,
  title   = {Consensus-based Distributed Optimization for Multi-agent Systems
             over Multiplex Networks},
  author  = {Rodr{\'i}guez-Camargo, Christian D. and
             Urquijo-Rodr{\'i}guez, Andres F. and
             Mojica-Nava, Eduardo},
  journal = {IEEE Transactions on Network Science and Engineering},
  year    = {2024},
  doi     = {10.1109/TNSE.2024.10772699},
  note    = {arXiv preprint arXiv:2304.01875}
}
```

---

## 🤝 Acknowledgments

This work was partially supported by:

- **Minciencias** Grant CT 542-2020 — *Programa de Investigación en Tecnologías Emergentes para Microrredes Eléctricas Inteligentes con Alta Penetración de Energías Renovables*
- **Uniminuto** VIII Convocatoria para el Desarrollo y Fortalecimiento de los Grupos de Investigación (code C119-173)
- **Uniminuto** Convocatoria de investigación para prototipado de tecnologías que promueven el cuidado o la restauración del medioambiente (code CPT123-200-5220)
- **EPSRC** (Grants No. EP/R513143/1 and No. EP/T517793/1)

---

## 👥 Authors

| Name | Affiliation | Contact |
|------|-------------|---------|
| **Christian D. Rodríguez-Camargo** | AMOPP Group, University College London & PAAS-UN, Universidad Nacional de Colombia | christian.rodriguez-camargo.21@ucl.ac.uk |
| **Andres F. Urquijo-Rodríguez** | Pontificia Universidad Javeria & Universidad Nacional de Colombia | afurquijor@unal.edu.co |
| **Eduardo Mojica-Nava** | Dept. of Electrical and Electronics Engineering, Universidad Nacional de Colombia | eamojican@unal.edu.co |

---

## 📬 Contact

For questions about the code or the paper, please open an [issue](../../issues) or contact the corresponding authors listed above.

---


