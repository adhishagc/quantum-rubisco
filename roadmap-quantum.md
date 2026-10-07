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

## 1. Biochemical Motivation & The Quantum Bottleneck

### Where the Empirical Data Demonstrates >1D Space
A longstanding dogma in plant biochemistry (*Savir et al. 2010*, *Tcherkez et al. 2006*) held that Rubisco kinetics is strictly 1-dimensional—that all enzymes are trapped along a rigid 1D Pareto frontier where optimizing one parameter deterministically degrades others.

The expanded dataset compiled in this repository refutes this 1D constraint with quantitative statistical evidence:

1. **Principal Component Analysis Dimensionality Drop ([Figure S8 PCA](notebooks/Figure%20S8%20PCA.ipynb)):**
   - **Historical Dataset (Savir et al. 2010, $N \approx 18$ Form I Rubiscos):**
     - **PC1 explained 91.0% of total variance** (PC2 explained only 6.2%). Because $>91\%$ of variance was captured by a single axis, previous authors concluded Rubisco was rigidly 1-dimensional.
   - **Expanded Form I Dataset ($N > 200$ Rubiscos):**
     - **PC1 drops to 71.3%** of variance.
     - **PC2 increases to 14.5%** of variance.
     - **PC3 increases to 9.8%** of variance.
     - **PC4 increases to 4.4%** of variance.
   - **Full Dataset (Form I, II, III combined):**
     - **PC1 explains only 61.6%**, while **PC2 explains 27.3%**.
   - *Conclusion:* Between **28.7% and 38.4% of total kinetic variance** lies orthogonal to the primary axis. Rubisco is not trapped on a 1D curve; it explores an open multi-dimensional volume.

2. **Decoupling of Carboxylation Rate and Specificity ([Figure 5](notebooks/Figure5%20Previous%20kcatC%20Correlations.ipynb)):**
   - If Rubisco were 1D, knowing $k_{\text{cat},C}$ would deterministically predict specificity $S_{C/O}$ ($R^2 \approx 1.0$).
   - In the expanded Form I dataset, log-scale Pearson $R = -0.56 \implies R^2 \approx 0.31$.
   - **$k_{\text{cat},C}$ explains only 31% of the variance in $S_{C/O}$**—leaving **69% of the variation unaccounted for by any 1D trade-off line**.
   - Across all Rubisco forms globally, $R = 0.03 \implies R^2 \approx 0.001$ (effectively zero correlation).

3. **Covariance Residuals & Multi-Dimensional Spread ([Figure S9](notebooks/FigureS9%20Residual%20Analysis.ipynb)):**
   - Orthogonal regression residuals show significant non-zero second eigenvalues in log-space covariance matrices ($[0.086, 0.020]$ for $k_{\text{cat},C}\text{--}K_C$ and $[0.034, 0.006]$ for $k_{\text{cat},C}\text{--}S_{C/O}$), proving genuine physical dispersion beyond experimental measurement noise.

---

### What Exact Bottleneck Requires Quantum Computing?

> **Clarification:** Fitting the CSV data, computing PCA, or running ODR regressions does **not** require quantum computing; classical computers execute these in milliseconds.  
> The **quantum bottleneck** lies in calculating the **underlying transition-state potential energy surfaces** that generated these kinetic values in nature.

#### 1. The Multi-Reference Open-Shell Triplet $\text{O}_2$ Problem
The core finding of this repository is **Mechanistic Proposal #2**: Rubisco is physically constrained by the challenge of discriminating $\text{CO}_2$ from $\text{O}_2$ at the addition transition state.
- **Carboxylation ($\text{CO}_2$ Addition):** A closed-shell singlet intermediate (RuBP 2,3-enediolate) attacks closed-shell electrophilic $\text{CO}_2$. Single-reference methods (like DFT) capture this reasonably well.
- **Oxygenation ($\text{O}_2$ Addition):** Ground-state molecular oxygen is an open-shell **triplet** ($^3\Sigma_g^-$ with 2 unpaired electrons, $S=1$). To react with singlet RuBP ($S=0$), the system must undergo **radical single-electron transfer** to form a transient superoxide radical pair ($\text{RuBP}^{\bullet+} / \text{O}_2^{\bullet-}$) followed by spin-orbit intersystem crossing.

Classical single-reference Density Functional Theory (DFT) fails dramatically for open-shell diradicals:
- It suffers from severe **spin contamination** ($\langle S^2 \rangle \ne 0$) and self-interaction errors.
- Typical DFT errors in activation barriers ($\Delta G^\ddagger$) for oxygenation are **$\pm 3\text{ to }6\text{ kcal/mol}$**.

#### 2. The Sub-kcal/mol Scale of the Biological Signal
From transition-state theory:
$$\Delta \Delta G^\ddagger = \Delta G^\ddagger_{\text{oxygenation}} - \Delta G^\ddagger_{\text{carboxylation}} = RT \ln(S_{C/O})$$

Using empirical values from this repository:
- **Red Algae** ($S_{C/O} \approx 160$): $\Delta \Delta G^\ddagger = 3.00\text{ kcal/mol}$
- **C3 Crop Plants** ($S_{C/O} \approx 95$): $\Delta \Delta G^\ddagger = 2.70\text{ kcal/mol}$
- **Cyanobacteria** ($S_{C/O} \approx 45$): $\Delta \Delta G^\ddagger = 2.25\text{ kcal/mol}$

The entire biological difference between a standard crop Rubisco and an elite red algal Rubisco is **only $\approx 0.30\text{ kcal/mol}$** ($\sim 0.5\text{ mHa}$)!

| Computational Method | Typical Error in Barrier | Can it resolve the biological difference ($0.30\text{ kcal/mol}$)? |
|:---|:---:|:---:|
| **Classical DFT** | $\pm 3\text{--}6\text{ kcal/mol}$ | **No** (Error is 10–20× larger than the target signal) |
| **Classical CASPT2 / DMRG** | $\pm 1\text{--}2\text{ kcal/mol}$ | **No** (Limited to active spaces $\le 30\text{--}35$ spatial orbitals) |
| **Fault-Tolerant Quantum (QPE)** | **$< 0.5\text{ kcal/mol}$** | **Yes** (Achieves true chemical accuracy) |

#### 3. The Exponential Wall of Classical Active Spaces
To resolve open-shell transition states without DFT, quantum chemists use complete active space methods (CASSCF/CASPT2). However, the classical configuration space scales combinatorially:
$$\dim(\mathcal{H}) = \binom{2 N_{\text{orbitals}}}{N_{\text{electrons}}}$$

Modeling Rubisco's transition state accurately requires an active space spanning:
- The enediolate $\text{C1--C3}$ skeleton
- The substrate ($\text{O}_2$ or $\text{CO}_2$)
- The catalytic $\text{Mg}^{2+}$ ion
- Coordinating residues (carbamylated Lys201, His294, Asp203, Glu204)

This yields an active space of at least **$(40e, 40o)$ to $(64e, 64o)$**:
- A $(40e, 40o)$ active space contains over **$10^{22}$ Slater determinants**, exceeding the memory and compute capacity of any classical supercomputer.
- A **Fault-Tolerant Quantum Computer** using Quantum Phase Estimation (QPE) scales **polynomially** ($\mathcal{O}(N^3\text{ to } N^4)$ with Tensor Hypercontraction), making this active space tractable and enabling the first first-principles calculation of Rubisco's selectivity barrier.

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
