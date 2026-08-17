

[![CI](https://github.com/Abiyzelalem27/DeutschJozsaNoise/actions/workflows/python_CI.yml/badge.svg)](https://github.com/Abiyzelalem27/DeutschJozsaNoise/actions/workflows/python_CI.yml)

[![codecov](https://codecov.io/github/Abiyzelalem27/DeutschJozsaNoise/graph/badge.svg)](https://codecov.io/github/Abiyzelalem27/DeutschJozsaNoise)

# Quantum Algorithms 

A brute-force classical simulator for quantum circuits and quantum algorithms, with a focus on the Deutsch–Jozsa algorithm and quantum noise analysis.


## Features

* Quantum circuit simulation on classical hardware
* Deutsch–Jozsa algorithm implementation
* Oracle-based black-box functions
* Quantum measurement simulation 
* Single-qubit rotation gate errors
* Noise analysis at different stages of the algorithm
* Statistical measurement experiments


## Rotation-Gate Error Results

The following figures show the effect of single-qubit rotation errors on
the success probability of the Deutsch–Jozsa algorithm.

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

## Noise Model 

<p align="center">
  <img src="images/dja_noise_model.png" width="850">
</p>

<p align="center">
<b>Figure:</b> Modified Deutsch–Jozsa circuit with localized noise channels
applied at different stages of the computation.
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
├── Rotation_gate/
│   ├── Rotational_Error_Analysis.ipynb
│   ├── Rotational_Error_Analysis_and_Sensitivity.ipynb
│   ├── Rotational_Error_Performance_and_Scalability.ipynb
│   └── Rotational_Error_Scalability_and_Average_Success.ipynb
│
├── Depolarizing/
│   ├── Depolarizing_Noise_Analysis.ipynb
│   ├── Depolarizing_Noise_Analysis_and_Sensitivity.ipynb
│   ├── Depolarizing_Noise_Performance_and_Scalability.ipynb
│   └── dja_depolarizing_results_shots_1024.csv
│
├── Internship/
│   └── Final_Internship_Report.pdf
│
├── images/
│   ├── deutsch_jozsa_circuit.png
│   ├── dja_noise_model.png
│   └── oracle_operator.png
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


