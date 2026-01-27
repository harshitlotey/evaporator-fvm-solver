# Discretized Finite-Volume Evaporator Solver
 
> **[Read the full re port here (PDF)](./report.pdf)**

A 1D numerical simulation of a fin-and-tube evaporator using a discretized control-volume approach. This solver models refrigerant phase change, pressure drop, and heat transfer row-by-row, integrating established empirical correlations for high-fidelity performance prediction.

### Overview

This project implements a **Finite Volume Method (FVM)** approach to simulate the steady-state performance of a multi-circuit evaporator. Unlike lumped-parameter models, this solver discretizes the heat exchanger tubes into finite segments, solving mass, momentum, and energy conservation equations for each control volume.

Fluid & Operational Scope: Built on the CoolProp wrapper, the solver supports a comprehensive library of thermodynamic fluids (refrigerants). It is specifically engineered for subcritical heat exchange applications, accurately modeling the non-linear thermodynamic evolution across the saturation dome—capturing the complex transition from the two-phase mixture region into the superheated vapor regime.

### Key Features

* **Discretized Domain:** Divides the coil into user-defined segments to capture local property variations (quality, enthalpy, temperature).
* **Physics-Based Correlations:**
    * **Boiling Heat Transfer:** Implements **Kandlikar’s Correlation (1990)** to dynamically switch between nucleate boiling and convective boiling dominance based on local Boiling number ($Bo$) and Convection number ($Co$).
    * **Two-Phase Pressure Drop:** Uses the **Friedel Correlation (1979)** for frictional pressure drop, accounting for liquid/vapor interaction.
    * **Air-Side Heat Transfer:** Uses **Wang et al. (2000)** correlations for plain fin-and-tube geometry to determine Colburn $j$ - factor and friction factor $f$.


* **Robust Numerical Solving:**
    * Uses **Brent’s Method** (`scipy.optimize.root_scalar`) to iteratively solve the energy balance at the tube wall for every segment ($Q_{air}=Q_{ref}$).
    * Implements a **row-by-row air temperature propagation** scheme, calculating the degradation of air temperature as it moves through the coil depth.

### Governing Equations & Methodology

The solver iterates through the tube length $L$, updating the refrigerant state at each step $i \to i+1$:

**1. Energy Balance (Enthalpy Update):**

$$H_{i+1} = H_i + \frac{dQ_{segment}}{\dot{m}_{ref}}$$

**2. Wall Temperature Convergence:**
For each segment, the solver minimizes the residual function to find the wall temperature $T_{wall}$:

$$R(T_{wall}) = h_{air}A_{o}(T_{air} - T_{wall}) - h_{ref}A_{i}(T_{wall} - T_{ref})$$

Where $h_{ref}$ is derived from Kandlikar's correlation:

$$h_{ref} = \max(h_{nucleate}, h_{convective})$$

**3. Pressure Drop:**

$$P_{i+1} = P_i - (\Delta P_{friction} + \Delta P_{momentum})$$

### Installation

1. Clone the repository:
```bash
git clone https://github.com/harshitlotey/evaporator-fvm-solver.git

```

2. Install dependencies:
```bash
pip install -r requirements.txt

```
*Dependencies:* `CoolProp`, `numpy`, `scipy`, `matplotlib`

###  Usage

1. **Configure the Inputs:** Modify `evapProps.json` to define your geometry and operating conditions:
```json
{
    "Refrigerant": "R32",
    "No. of Rows": 2,
    "Tube OD": 7.0,
    "Inlet Air Temp": 27.0,
    ...
}

```

2. **Run the Solver:**
```bash
python3 evap.py

```
3. **Output:** The script generates console output for capacity of the HX and the output temperature of the refrigerant along with plots:
    * **P-h Diagram:** Visualizing the evaporation process.
    * **Quality vs. Length:** Tracking the phase change progress.
    * **Temperature Profile:** Refrigerant temperature along the coil length.

### Author

**Harshit Dhiman**
*R&D Engineer | Heat Transfer Specialist*
[LinkedIn](https://www.linkedin.com/in/harshitlotey) | [Email](mailto:harshitlotey@gmail.com)
