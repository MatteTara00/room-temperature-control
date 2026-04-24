# Room Temperature Control using LMIs

## Overview

This project focuses on the control of temperature in a four-room building using Linear Matrix Inequalities (LMIs). The system is modeled as a multi-input multi-output (MIMO) dynamic system, where each room is treated as a subsystem with its own temperature dynamics.

Different control architectures are investigated, including centralized, decentralized, and distributed approaches, and their performance is evaluated in both continuous and discrete time domains.

---

## Contributions

* Mathematical modeling of a multi-room thermal system
* Linearization around equilibrium operating conditions
* Design of state-feedback controllers using LMIs
* Implementation of centralized, decentralized, and distributed control architectures
* Stability and performance analysis in continuous and discrete time
* H₂ optimization and sector-based LMI design

---

## System Description

The system consists of:

* **4 states** → temperature of each room
* **4 inputs** → heating/cooling actuators
* **4 outputs** → temperature measurements

The building is modeled as a set of interconnected subsystems, where thermal coupling between rooms is taken into account.

As observed from system analysis, the coupling between rooms is relatively low, enabling the use of decentralized and distributed control strategies.

---

## Control Architectures

The project explores different control configurations:

### Centralized Control

* Single controller with full system knowledge
* Best performance, highest complexity

### Decentralized Control

* One controller per room
* No communication between subsystems

### Distributed Control

* Controllers communicate with neighboring subsystems
* Two configurations based on system interconnections

According to the analysis, all architectures are feasible and stable, with no fixed modes detected .

---

## LMI-Based Control Design

The control design is based on solving different LMI problems:

### Continuous-Time LMIs

* Stability conditions
* Performance constraints (eigenvalue placement)
* Sector constraints
* H₂ optimization

### Discrete-Time LMIs

* Stability
* Performance
* H₂ minimization

All LMI problems were found to be feasible across the considered architectures .

---

## Implementation

The entire project is implemented in **MATLAB**, including:

* System modeling
* LMI formulation and solution
* Controller synthesis
* Simulation and performance evaluation

---

## Workflow

1. System modeling and linearization
2. Subsystem decomposition
3. Stability analysis
4. LMI-based controller design
5. Simulation and comparison of architectures
6. Performance evaluation

---

## Requirements

* MATLAB
* Control System Toolbox
* Optimization / LMI Toolbox

---

## Usage

1. Open MATLAB
2. Load the project folder
3. Run the main scripts for:

   * system definition
   * controller synthesis
   * simulation

---

## Notes

* The system exhibits low inter-room coupling
* All control architectures are feasible
* Distributed control provides a good trade-off between performance and scalability

---
