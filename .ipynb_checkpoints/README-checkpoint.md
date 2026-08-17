

[![CI](https://github.com/Abiyzelalem27/DeutschJozsaNoise/actions/workflows/python_CI.yml/badge.svg)](https://github.com/Abiyzelalem27/DeutschJozsaNoise/actions/workflows/python_CI.yml)

[![codecov](https://codecov.io/github/Abiyzelalem27/DeutschJozsaNoise/graph/badge.svg)](https://codecov.io/github/Abiyzelalem27/DeutschJozsaNoise)

# Quantum Algorithms 

A Classical simulator for quantum circuits and quantum algorithms, with a focus on the Deutsch–Jozsa algorithm and quantum noise analysis.

## Features

* Brute-force state-vector simulation of the Deutsch–Jozsa algorithm
* Constant-0, constant-1, and balanced quantum-oracle implementations
* Coherent X-, Y-, and Z-axis rotation-error analysis
* Depolarizing-noise analysis at four circuit locations (E1–E4)
* Sensitivity, scalability, and shot-based performance evaluation

## Noise Model 

<p align="center">
  <img src="images/dja_noise_model.png" width="850">
</p>

<p align="center">
<b>Figure:</b> Modified Deutsch–Jozsa circuit with localized noise channels
applied at different stages of the computation.
</p>

### Constant-0 Function — 9 Qubits

<p align="center">
  <img src="images/rotation_success_constant_0_n9_X_q4.png"
       alt="Constant-0 function with X-axis rotation error" width="31%">
  <img src="images/rotation_success_constant_0_n9_Y_q4.png"
       alt="Constant-0 function with Y-axis rotation error" width="31%">
  <img src="images/rotation_success_constant_0_n9_Z_q4.png"
       alt="Constant-0 function with Z-axis rotation error" width="31%">
</p>

<p align="center">
  <b>Figure 1:</b> Success probability for the constant-0 function under
  X-, Y-, and Z-axis rotation errors. The rotation error is applied to
  qubit 4 in the nine-qubit case.
</p>

### Balanced Function — 9 Qubits

<p align="center">
  <img src="images/rotation_success_balanced_n9_X_q4.png"
       alt="Balanced function with X-axis rotation error" width="31%">
  <img src="images/rotation_success_balanced_n9_Y_q4.png"
       alt="Balanced function with Y-axis rotation error" width="31%">
  <img src="images/rotation_success_balanced_n9_Z_q4.png"
       alt="Balanced function with Z-axis rotation error" width="31%">
</p>

<p align="center">
  <b>Figure 2:</b> Success probability for the balanced function under
  X-, Y-, and Z-axis rotation errors. The rotation error is applied to
  qubit 4 in the nine-qubit case.
</p>
---
### Constant-1 Function — 9 Qubits

<p align="center">
  <img src="images/rotation_success_constant_1_n9_X_q4.png"
       alt="Constant-1 with X-axis rotation error" width="31%">
  <img src="images/rotation_success_constant_1_n9_Y_q4.png"
       alt="Constant-1 with Y-axis rotation error" width="31%">
  <img src="images/rotation_success_constant_1_n9_Z_q4.png"
       alt="Constant-1 with Z-axis rotation error" width="31%">
</p>

<p align="center">
  <b>Figure 2:</b> Success probability for the constant-1 function under
  X-, Y-, and Z-axis rotation errors. The error is applied to qubit 4
  in the nine-qubit case.
</p>

---

### Rotation-Error Scalability

<p align="center">
  <img src="images/scalability_balanced_E3_after_oracle_X_90deg.png"
       alt="Balanced-function scalability with error after oracle" width="47%">
  <img src="images/scalability_balanced_E4_after_final_H_X_90deg.png"
       alt="Balanced-function scalability with error after final Hadamard" width="47%">
</p>

<p align="center">
  <b>Figure 4:</b> Scalability of the balanced-function implementation
  under a 90-degree X-rotation error introduced after the oracle (E3)
  and after the final Hadamard layer (E4).
</p>

---

## Depolarizing-Noise Results

The following figure shows the effect of depolarizing noise introduced at different stages of the Deutsch–Jozsa algorithm.

### Constant-0 Function — 9 Qubits

<p align="center">
  <img src="images/depolarizing_constant0_9qubits.png"
       alt="Effect of depolarizing noise on the Deutsch–Jozsa algorithm for the constant-0 function"
       width="850">
</p>

<p align="center">
  <b>Figure 5:</b> Success probability for the constant-0 function under
  depolarizing noise introduced before the first Hadamard layer (E1),
  after the first Hadamard layer (E2), after the oracle (E3), and after
  the final Hadamard layer (E4). The ideal circuit maintains unit success
  probability.
</p>

## Repository Structure

```text
DeutschJozsaNoise/
│
├── quantum_algorithms/
│   ├── __init__.py
│   ├── operators.py
│   ├── black_box.py
│   ├── deutsch_jozsa.py
│   ├── depolarizing.py
│   └── plotting.py
│
├── Rotation_gate/   # notebook for  rotation-error analysis, visualizations, and results
├── Depolarizing/  # notebook for depolarizing-noise analysis, visualizations, and results
├── Internship/  # Internship report and PowerPoint presentation
├── images/
│
├── README.md
├── pyproject.toml
├── LICENSE
└── .gitignore
```

## Example

```python
import quantum_algorithms as qa

state = qa.deutsch_jozsa.deutsch_jozsa(
    3,
    qa.deutsch_jozsa.f_balanced_parity
)

print(state)
```


## 📚 References

* **Quantum Information and Quantum Simulation (QIQS) Group**  
  Friedrich Schiller University Jena, Germany  
  Research group and academic environment associated with the course.  
  Website: [qiqs-jena.de](https://qiqs-jena.de/)

* **Nielsen & Chuang**  
  *Quantum Computation and Quantum Information*

 

## 👤 Author

* **Abiy Zelalem Tegegne**

GitHub: [Abiyzelalem27](https://github.com/Abiyzelalem27) 


