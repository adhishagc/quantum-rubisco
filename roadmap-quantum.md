# Quantum Computing Roadmap: Rubisco Catalysis & Kinetic Optimization

[![Status: Conceptual Roadmap](https://img.shields.io/badge/Status-Research%20Roadmap-blue.svg)](#)
[![Focus: FTQC%20%7C%20NISQ%20%7C%20QML](https://img.shields.io/badge/Focus-FTQC%20%7C%20NISQ%20%7C%20QML-purple.svg)](#)
[![Target: Rubisco%20Transition%20State](https://img.shields.io/badge/Target-Transition%20State%20Discrimination-green.svg)](#)

This document establishes a technical roadmap for applying **Fault-Tolerant Quantum Computing (FTQC) Resource Estimation**, **NISQ hybrid quantum chemistry**, **Quantum Machine Learning (QML)**, and **Quantum Optimization** to the Rubisco kinetic datasets and catalytic mechanisms documented in this repository.

---

## Table of Contents
1. [Biochemical Motivation & Quantum Advantage](#1-biochemical-motivation--quantum-advantage)
2. [Research Tracks](#2-research-tracks)
   - [Track 1: FTQC Resource Estimation for Transition-State Discrimination](#track-1-ftqc-resource-estimation-for-transition-state-discrimination)
   - [Track 2: NISQ Hybrid Chemistry (QM/MM + VQE on the Triplet Spin-Crossing Complex)](#track-2-nisq-hybrid-chemistry-qmmm--vqe-on-the-triplet-spin-crossing-complex)
   - [Track 3: Quantum Machine Learning (QML) on Empirical Kinetic Datasets](#track-3-quantum-machine-learning-qml-on-empirical-kinetic-datasets)
   - [Track 4: Quantum Combinatorial Optimization (QUBO/QAOA) for Directed Evolution](#track-4-quantum-combinatorial-optimization-quboqaoa-for-directed-evolution)
3. [Phased Implementation Timeline](#3-phased-implementation-timeline)
4. [Software Ecosystem & Toolchains](#4-software-ecosystem--toolchains)
5. [Key References](#5-key-references)

---

## 1. Biochemical Motivation & Quantum Advantage

### The Core Problem: Discrimination at the Transition State
As established by the empirical findings in this repository (*Flamholz et al., Biochemistry 2019*), Rubisco's primary evolutionary constraint is **Mechanistic Proposal #2**: the simultaneous coupling of carboxylation and oxygenation catalytic efficiencies ($k_{\text{cat},C}/K_C$ vs $k_{\text{cat},O}/K_O$, $R = 0.94$).

Specificity ($S_{C/O}$) is dictated by the difference in activation free energy barriers between carboxylation and oxygenation:

$$\Delta \Delta G^\ddagger = \Delta G^\ddagger_{\text{oxygenation}} - \Delta G^\ddagger_{\text{carboxylation}} = RT \ln(S_{C/O})$$

Across the empirical dataset:
- **C3 Plants** ($S_{C/O} \approx 90\text{--}105$): $\Delta \Delta G^\ddagger \approx 2.65\text{--}2.75\text{ kcal/mol}$
- **Red Algae** ($S_{C/O} \approx 140\text{--}166$): $\Delta \Delta G^\ddagger \approx 2.95\text{--}3.05\text{ kcal/mol}$
- **Cyanobacteria** ($S_{C/O} \approx 40\text{--}50$): $\Delta \Delta G^\ddagger \approx 2.20\text{--}2.35\text{ kcal/mol}$

The entire biological difference between crop plants and elite red algal Rubiscos rests on a sub-kcal/mol free energy window:
$$\delta (\Delta \Delta G^\ddagger) \approx 0.3\text{--}0.7\text{ kcal/mol}$$

### Why Classical Electronic Structure Fails
1. **Multi-Reference Open-Shell Character:** Ground-state $\text{O}_2$ is a triplet ($^3\Sigma_g^-$), whereas the RuBP 2,3-enediolate is a closed-shell singlet. The oxygenation pathway involves electron transfer forming a transient radical pair ($\text{RuBP}^{\bullet+} / \text{O}_2^{\bullet-}$) and spin-orbit intersystem crossing.
2. **DFT Errors:** Density Functional Theory (DFT) exhibits severe self-interaction and spin-contamination errors ($\pm 3\text{--}6\text{ kcal/mol}$), dwarfing the subtle physical barrier that distinguishes Rubisco variants.
3. **Classical Multi-Reference Limits:** Classical CASPT2 and DMRG become computationally intractable when the active space exceeds $\sim 30\text{--}40$ spatial orbitals, preventing simultaneous simulation of the substrate, the enediolate core, the catalytic $\text{Mg}^{2+}$ ion, and the primary coordination shell.

---

## 2. Research Tracks

```
                        ┌────────────────────────────────────────────────────────┐
                        │      Rubisco Quantum Computing Research Roadmap       │
                        └──────────────────────────┬─────────────────────────────┘
                                                   │
         ┌─────────────────────┬───────────────────┴─────────────────┬─────────────────────┐
         ▼                     ▼                                     ▼                     ▼
┌─────────────────┐   ┌─────────────────┐                   ┌─────────────────┐   ┌─────────────────┐
│     Track 1     │   │     Track 2     │                   │     Track 3     │   │     Track 4     │
│  FTQC Resource  │   │  NISQ Hybrid    │                   │   QML Kinetic   │   │  QUBO Enzyme    │
│   Estimation    │   │ QM/MM Simulation│                   │   Regression    │   │     Design      │
└─────────────────┘   └─────────────────┘                   └─────────────────┘   └─────────────────┘
```

---

### Track 1: FTQC Resource Estimation for Transition-State Discrimination
*Primary Objective:* Quantify the exact quantum computing hardware requirements (logical qubits, T-states, runtime, surface code cycles) to calculate Rubisco's carboxylation vs. oxygenation transition states within chemical accuracy ($< 1\text{ kcal/mol}$).

#### 1. Active Space Formulations
Using high-resolution crystallographic coordinates (e.g., PDB: `8RUC`, `1WDD`):
- **Minimal Active Space $(22e, 22o)$ (44 spin-orbitals):**  
  RuBP C1–C3 fragments, incoming $\text{CO}_2$ or $\text{O}_2$, and coordinating $\text{Mg}^{2+}$ d-orbitals.
- **Moderate Active Space $(40e, 40o)$ (80 spin-orbitals):**  
  Enediolate core + $\text{Mg}^{2+}$ + carbamylated Lys201 + gas substrate + key catalytic residues (His294, Asp203, Glu204).
- **Extended Active Space $(64e, 64o)$ (128 spin-orbitals):**  
  Full first-coordination sphere, bridging water molecules, and second-shell electrostatic donor residues.

#### 2. Compilation & Block-Encoding Paradigms
Benchmark state-of-the-art Hamiltonian representations:
- **Double Factorization (DF):** Rank reduction of the two-electron integral tensor.
- **Tensor Hypercontraction (THC):** Low-rank factorization reducing T-complexity to $\mathcal{O}(N)$.
- **Qubitization & Quantum Phase Estimation (QPE):** Compute spectral precision $\epsilon = 0.5\text{ mHa}$ ($0.31\text{ kcal/mol}$) to reliably resolve the $S_{C/O}$ gap between lineages.

#### 3. Architectural Benchmarking
- Evaluate on:
  - Rotated Surface Codes at physical error rates $p = 10^{-3}$ and $p = 10^{-4}$.
  - Magic State Distillation (15-to-1 or 20-to-4 T-factories).
  - Surface code cycle time assumption: $1\ \mu\text{s}$ (superconducting) vs $100\ \mu\text{s}$ (trapped ion).

---

### Track 2: NISQ Hybrid Chemistry (QM/MM + VQE on the Triplet Spin-Crossing Complex)
*Primary Objective:* Simulate simplified active-site models on near-term hardware and state-vector simulators to model the $\text{RuBP}-\text{O}_2$ initial encounter complex.

#### 1. Classical Embedding Framework
- **Layer 1 (MM):** Entire hexadecameric Rubisco enzyme scaffold + solvent box modeled with CHARMM/AMBER force fields.
- **Layer 2 (DFT Embedding):** Surrounding active-site pocket (100–150 atoms) providing electrostatic embedding.
- **Layer 3 (Quantum Core):** 8-to-16 qubit active space representing the reacting enediolate $\text{C2}=\text{C3}$ bond and ground-state $\text{O}_2$.

#### 2. Algorithmic Implementation
- **ADAPT-VQE:** Dynamically grows compact problem-tailored ansatzes with shallow circuit depth, mitigating noise.
- **Quantum Subspace Expansion (QSE):** Computes singlet and triplet excited states to identify the potential energy crossing point along the $\text{C2--O}$ bond-formation coordinate.
- **Error Mitigation:** Apply zero-noise extrapolation (ZNE), readout error mitigation, and symmetry verification (particle number and total spin conservation $S^2, S_z$).

---

### Track 3: Quantum Machine Learning (QML) on Empirical Kinetic Datasets
*Primary Objective:* Leverage the curated empirical dataset ([DatasetS2_RubiscoKinetics_Merged.csv](file:///Users/adhishagammanpila/Documents/Engineering/rubisco/data/DatasetS2_RubiscoKinetics_Merged.csv) and [DatasetS4_RubiscoKineticsFull_Merged.csv](file:///Users/adhishagammanpila/Documents/Engineering/rubisco/data/DatasetS4_RubiscoKineticsFull_Merged.csv)) to test whether Quantum Kernels offer expressivity advantages on small, highly non-linear biological datasets ($N \approx 300$).

#### 1. Feature Representation
- **Sequence Embeddings:** Extract representations from pre-trained protein language models (ESM-2, ProtTrans) for the Rubisco large subunit (RbcL).
- **Dimensionality Reduction:** Compress representations via PCA / UMAP to $d \in [8, 16]$ dimensions.
- **Targets:** Continuous kinetic parameters: $\ln(k_{\text{cat},C})$, $\ln(S_{C/O})$, and $\ln(K_C)$.

#### 2. Model Architectures
- **Quantum Kernel Ridge Regression (QKRR):**
  - Construct parameterized quantum feature maps (ZZ-feature maps, Hamiltonian evolution circuits).
  - Evaluate quantum kernel matrix:
    $$K_{ij} = |\langle 0^{\otimes n} | U^\dagger(x_j) U(x_i) | 0^{\otimes n} \rangle|^2$$
  - Compare prediction accuracy and generalization metrics against classical RBF, linear, and Gaussian Process regressors.
- **Projected Quantum Kernels:** Mitigate quantum kernel exponential concentration (barren plateaus in kernel space) by projecting into reduced physical subspace observables.

---

### Track 4: Quantum Combinatorial Optimization (QUBO/QAOA) for Directed Evolution
*Primary Objective:* Formulate multi-residue Rubisco enzyme engineering as a constrained combinatorial optimization problem.

#### 1. Multi-Objective Objective Function
$$\min_{x} \left[ - w_1 \cdot \Delta k_{\text{cat},C}(x) + w_2 \cdot \max(0, S_{\text{target}} - S_{C/O}(x)) + w_3 \cdot \Delta \Delta G_{\text{fold}}(x) \right]$$

Subject to:
- Maximum mutation budget: $\sum_i x_i \le K$ (e.g., $K \le 5$ simultaneous mutations).
- Fold stability constraint: $\Delta \Delta G_{\text{fold}}(x) \le 0\text{ kcal/mol}$.

#### 2. Execution Platforms
- **Quantum Annealing (D-Wave Advantage):** Embed dense epistasis graphs into Pegasus/Zephyr topologies using minor-embedding algorithms.
- **QAOA (Gate-Based):** Run depth $p = 1\text{--}3$ QAOA circuits on simulators and gate-based quantum processors.

---

## 3. Phased Implementation Timeline

| Phase | Milestone | Primary Deliverable | Tools / Platforms |
|:---:|:---|:---|:---|
| **Phase 1** *(Month 1–3)* | **Kinetic Inversion & Active Space Setup** | Map dataset $S_{C/O}$ values to $\Delta \Delta G^\ddagger$; generate PySCF molecular orbitals for active-site transition states. | PySCF, OpenFermion, PDB |
| **Phase 2** *(Month 4–6)* | **FTQC Resource Estimation Benchmark** | First publication on full-scale logical qubit, T-count, and runtime resource estimation for Rubisco $\text{CO}_2/\text{O}_2$ discrimination. | Azure Quantum QRE, Qualtran |
| **Phase 3** *(Month 7–9)* | **QML Regression Benchmark** | Benchmark Quantum Kernel Ridge Regression against classical GPR on Dataset S2/S4 kinetic prediction. | PennyLane, Qiskit Machine Learning |
| **Phase 4** *(Month 10–12)* | **NISQ VQE Simulation & QUBO Optimization** | ADAPT-VQE simulation of the $\text{RuBP}-\text{O}_2$ triplet model; QUBO multi-point mutant design on D-Wave. | Qiskit Nature, D-Wave Ocean SDK |

---

## 4. Software Ecosystem & Toolchains

- **Quantum Chemistry & Electronic Structure:**
  - [PySCF](https://github.com/pyscf/pyscf): Molecular orbital generation, CASSCF active space selection, integral export.
  - [OpenFermion](https://github.com/quantumlib/OpenFermion): Electronic structure Hamiltonian mapping (Jordan-Wigner, Bravyi-Kitaev).
- **Fault-Tolerant Resource Estimation:**
  - [Qualtran](https://github.com/quantumlib/qualtran) (Google Quantum AI): Detailed fault-tolerant algorithm compilation and gate synthesis.
  - [Azure Quantum Resource Estimator](https://learn.microsoft.com/en-us/azure/quantum/user-guide-resource-estimator) (Microsoft): Space-time volume, physical qubit, and T-factory estimation.
  - [pyLIQTR](https://github.com/isi-usc-edu/pyLIQTR): Qubitization and quantum signal processing circuits.
- **NISQ Simulation & QML:**
  - [PennyLane](https://github.com/PennyLaneAI/pennylane): Differentiable quantum programming and quantum kernel workflows.
  - [Qiskit Nature](https://github.com/qiskit-community/qiskit-nature): VQE, ADAPT-VQE, and active-space reduction.
- **Optimization:**
  - [D-Wave Ocean SDK](https://github.com/dwavesystems/dwave-ocean-sdk): Constrained Quadratic Models (CQM) and QUBO solvers.

---

## 5. Key References

1. **Rubisco Kinetic Dataset & Trade-Off Analysis:**
   - Flamholz, A. I. et al. *Revisiting Trade-offs between Rubisco Kinetic Parameters.* **Biochemistry 2019**, 58 (33), 3365–3376. [DOI: 10.1021/acs.biochem.9b00237](https://doi.org/10.1021/acs.biochem.9b00237).
2. **Quantum Chemistry Resource Estimation Benchmarks:**
   - Reiher, M. et al. *Elucidating reaction mechanisms on quantum computers.* **PNAS 2017**, 114 (29), 7555–7560. [DOI: 10.1073/pnas.1619152114](https://doi.org/10.1073/pnas.1619152114).
   - Lee, S. et al. *Even More Efficient Quantum Computations of Chemistry Through Tensor Hypercontraction.* **PRX Quantum 2021**, 2, 030305. [DOI: 10.1103/PRXQuantum.2.030305](https://doi.org/10.1103/PRXQuantum.2.030305).
   - von Burg, V. et al. *Quantum computing enhanced computational catalysis.* **Phys. Rev. Research 2021**, 3, 033055. [DOI: 10.1103/PhysRevResearch.3.033055](https://doi.org/10.1103/PhysRevResearch.3.033055).
3. **Transition-State Theory of Rubisco:**
   - Tcherkez, G. G. B. et al. *Despite slow catalysis and confused substrate specificity, all ribulose bisphosphate carboxylases may be nearly perfectly optimized.* **PNAS 2006**, 103 (19), 7246–7251. [DOI: 10.1073/pnas.0600605103](https://doi.org/10.1073/pnas.0600605103).
   - Savir, Y. et al. *The catalytic beauty of RuBisCO: Carboxylation, oxygenation, and evolutionary tradeoffs.* **PLOS Comput. Biol. 2010**, 6 (9), e1000923. [DOI: 10.1371/journal.pcbi.1000923](https://doi.org/10.1371/journal.pcbi.1000923).
