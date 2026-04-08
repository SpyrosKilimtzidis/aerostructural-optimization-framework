# Aerostructural Optimization Framework

## Overview

This repository presents a high-fidelity **aerostructural optimization framework** for the preliminary design of composite UAV wings. The framework couples **Computational Fluid Dynamics (CFD)** and **Finite Element Method (FEM)** analyses within a **Surrogate-Based Optimization (SBO)** environment.

The methodology enables the simultaneous optimization of aerodynamic and structural design variables, achieving high-performance configurations while significantly reducing computational cost through surrogate modeling techniques.

---

## Key Features

* Coupled **CFD–FEM aerostructural analysis**
* **Parametric wing geometry generation**
* Automated **mesh generation and simulation workflow**
* **Kriging surrogate modeling**
* **Expected Improvement (EI)** infill strategy
* Multi-variable optimization including:

  * Aerodynamic parameters (aspect ratio, taper ratio, sweep, twist)
  * Structural parameters (skin thickness, spars, ribs, layout)

---

## Methodology

### Aerodynamic Analysis

* High-fidelity CFD simulations using RANS equations
* Spalart–Allmaras turbulence model
* Low-Reynolds-number regime suitable for UAV applications

### Structural Analysis

* FEM model developed in MSC Nastran
* Shell (CQUAD4) and beam (CBEAM) elements
* Composite material modeling (Hexcel IM7/8552)
* Evaluation of:

  * Strength (Failure Index)
  * Stiffness (buckling constraints)

### Coupling Strategy

* Pressure field mapping from CFD to FEM
* Automated geometry and mesh generation via parametric scripts
* Fully integrated aerostructural workflow

### Optimization Framework

* Surrogate-Based Optimization (SBO)
* Kriging models for objective and constraints
* Expected Improvement (EI) for adaptive sampling
* Design space exploration using Latin Hypercube Sampling (LHS)

---

## Optimization Problem

**Objective:**

* Maximize UAV range

**Design Variables:**

* Aerodynamic: aspect ratio, taper ratio, sweep angle, twist
* Structural: thicknesses, spar locations, rib spacing, stringer spacing

**Constraints:**

* Lift coefficient requirement
* Structural strength (Failure Index)
* Buckling eigenvalue (stability)
* Safety factors for structural components

---

## Workflow

1. Define design variables and bounds
2. Generate parametric wing geometry
3. Perform CFD analysis
4. Map aerodynamic loads to FEM model
5. Run structural and buckling analyses
6. Train surrogate models (Kriging)
7. Select new design points via EI criterion
8. Iterate until convergence

---

## Example Usage

```matlab
% Example workflow

% Define design variables
params = defineDesignVariables();

% Generate geometry
generateWingGeometry(params);

% Run CFD analysis
runCFD();

% Map loads and run FEM
runFEM();

% Perform optimization step
optimizeSBO();
```

---

## Repository Structure

```
src/
  geometry/        % Parametric geometry generation
  fem/             % FEM model creation and analysis
  aero/            % CFD setup and execution
  optimization/    % Surrogate modeling and SBO

examples/
  run_case.m       % Example optimization workflow

docs/
  figures/         % Diagrams and workflow figures
```

---

## Related Publication

This repository supports the methodology presented in:

**Kilimtzidis et al. (2026)**
*Coupled Aerostructural Optimization of a Composite Low Reynolds Wing Using Surrogate Modeling Techniques*
Journal: *Drones*

---

## Applications

* UAV wing design
* Aerostructural optimization
* Multidisciplinary Design Optimization (MDO)
* Preliminary aircraft design

---

## Requirements

* MATLAB
* MSC Patran / Nastran
* ANSYS Fluent (or equivalent CFD solver)

---

## Author

Spyridon Kilimtzidis
PhD in Aerostructural Optimization

---

## License

MIT License
