# Advanced Numerical Methods: CPU Cooling Simulation

A MATLAB-based numerical simulation package for modeling 2D transient heat transfer within an active-cooled aluminum CPU heatsink.

The project solves the transient heat conduction equation to evaluate the thermal distribution and maximum operating temperatures of a heatsink under load. It evaluates custom finite difference solvers against an industry-standard finite element reference, focusing on a complex, state-dependent "Smart Fan" boundary condition.

## Project Overview
The project was developed and evaluated in three sequential phases:

### Phase 1: Explicit FTCS Solver (`solve_Explicit.m`)
* Implements the Forward-Time Central-Space (FTCS) finite difference method.
* Highly optimized utilizing MATLAB matrix vectorization rather than nested spatial loops.
* Automatically manages strict time-step constraints (CFL conditions) to guarantee numerical stability.
* Accurately captures the thermal evolution from initial 1D transient heating to a fully 2D steady-state distribution.

### Phase 2: Implicit ADI Solver (`solve_ADI.m`)
* Employs the Alternating Direction Implicit (ADI) method, splitting each time step into two half-steps.
* Solves the implicit matrices efficiently by pre-allocating and factorizing the tridiagonal coefficient matrices prior to the time-stepping loop.
* Utilizes the Thomas algorithm for rapid matrix inversion at each time step.

### Phase 3: MATLAB PDE Toolbox Validation (`main_Reference.m`)
* Serves as the mathematical finite-element benchmark for the custom FVM solvers.
* Generates an unstructured triangular mesh (`generateMesh`) and solves the exact same geometry, material properties, and non-linear boundary conditions using `createpde`.

## Solver & Modeling Notes
* **2D Planar Assumption:** The heatsink is modeled with a uniform cross-section ($20\text{mm} \times 20\text{mm}$), allowing the problem to be solved dynamically in the 2D $xy$-plane without volume equations.
* **CPU Heat Generation:** Modeled as a constant Neumann heat flux boundary condition ($q'' = 200,000 \text{ W/m}^2$) applied to the bottom edge.
* **Smart Cooling Fan (Top Boundary):** The top edge implements a dynamic, non-linear Robin boundary condition. A proportional controller ($K=100$) monitors the maximum surface temperature. If the temperature exceeds the target of $90^\circ\text{C}$, the controller scales the convective heat transfer coefficient ($h_{\text{top}}$) from an idle $1000 \text{ W/m}^2\text{K}$ up to a maximum physical limit of $4000 \text{ W/m}^2\text{K}$. At steady state ($t \approx 50\text{s}$), $h_{\text{top}}$ stabilizes around $3540 \text{ W/m}^2\text{K}$.

## Mesh Independency & Error Analysis
To ensure the accuracy of the custom finite difference solvers, a rigorous mesh independency and relative error analysis was conducted against the MATLAB PDE Toolbox reference.

The validation confirms critical numerical behaviors:
1. **Mesh Sensitivity:** A coarse mesh ($2 \times 4$) exhibited unacceptably high initial errors. Conversely, an overly fine mesh ($80 \times 40$) caused the **ADI solver to become numerically unstable**. This instability stems from the Thomas algorithm's sensitivity to sharp temperature gradients combined with the non-linear boundary logic.
2. **Optimal Configuration:** An intermediate mesh ($10 \times 20$) mitigated the Fourier number-induced initial error spikes and provided an optimal balance, keeping the relative error across solvers strictly under $2\%$.
3. **Conclusion:** For this specific non-linear boundary problem, the **Explicit FTCS** method proved to be the most efficient and robust choice, successfully managing stability without the boundary gradient issues seen in the implicit approach.

## Repository Structure
```text
Advanced-Numerical-Methods-CPU-Cooling/
├── .gitignore
├── README.md
├── LICENSE
├── main_Unified.m         # Executes and plots all solvers side-by-side
├── main_Explicit.m        # Standalone execution for Explicit FTCS
├── main_ADI.m             # Standalone execution for ADI
├── main_Reference.m       # Standalone execution for PDE Toolbox
├── solve_Explicit.m       # Vectorized FTCS computational core
└── solve_ADI.m            # Pre-factorized ADI computational core
```

## Requirements and Installation

### Prerequisites
- MATLAB R2020a or higher
- Partial Differential Equation (PDE) Toolbox (required for the Reference solver)

### Setup

1. Clone the repository:

```bash
git clone https://github.com/Alireza-Kianifard/Advanced-Numerical-Methods-CPU-Cooling.git
cd Advanced-Numerical-Methods-CPU-Cooling
```

2. Open MATLAB and navigate to the cloned repository directory.

3. Run the unified simulation (or any of the other "main" codes) from the MATLAB command window to see the side-by-side comparison:
```bash
matlab
run main_Unified.m
```