
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
