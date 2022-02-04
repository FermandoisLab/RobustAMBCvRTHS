# Adaptive model-based compensation for linear case

## Description

In this example the adaptive model-based compensator is designed, calibrated and applied to the RTHS benchmark. A calibrated gain matrix Gamma is provided to run several simulations with linear systems includding uncertainty. However, the adaptive gain matrix Gamma can be recalculated and analized using the calibration procedure.

## Files 

### AMB_L_1_vRTHS.m  (Matlab script)

In this script the problem is defined and the virtual RTHS can be run multiple times for considering uncertainty.

Definitions:
  - Simulation parameters.
  - Transfer system dynamics.
  - Compensator: initial parameters, filters, adaptive gain matrix (pre-calibrated matrix is provided).
  - Earthquake.
  - Reference structure.
Simulations:
  For each simulation experimental and numerical substructure are defined and simulated.
  - Displacement plots.
  - Adaptive parameters plots.

### AMB_L_2_calibration.m  (Matlab script)

In this script the calibration simulation to evaluate different adaptive gains is defined.

Definitions:
  - Calibration plants.
  - Calibration earthquake.
  - Numerical substructures.
Simulation:
  - Calibration simulation example with predefined gains.

Function definition:
  -Matlab function with adaptive gains coefficients as inputs and mean J2 error for N calibration simulations as outputs (fun()).

### AMB_L_3_sensitivy_inputs.m  (Matlab script)

A sensitivity analysis for a fixed gain matrix is carried out.

Function definition:
  - Matlab function with adaptive gains coefficients as inputs and J2 errors and simulation parameters as outputs (plants parameters, eartquakes scales, numerical structures) (fun2()).

Plots:
  - J2 vs calibration parameters.
  - Histograms.

### AMB_L_4_sensitivity_gains.m  (Matlab script)

A sensitivity analysis for R2 indicator for different adaptive gains matrix is carried out.

Definitions:
- Adaptive gains range to be evaluated.

Evaluation:
- Combinations of Gamma are evaluted using fun().
- R2 plots.

### AMB_L_5_optimization.m  (Matlab script)

Particle swarm optimization is implemented to find the optimal gains and the neighborhood is analyzed evaluating a grid of Gamma in different projections of the search space.

- Particle swarm optimization for a bounded search space.
- R2 evaluation around the optimal obtained value. (This may take time).
- R2 plots in color maps.

### AMB_L_6_R2function.m  (Matlab function)

Matlab function with R2 as outptut and simulation parameters as inputs. The functions runs N simulations with random parameters and a fixed gain matrix diag(10.^x).

### AMB_L_7_J2andInputs.m  (Matlab function)

Matlab function with J2 errors and simulated calibration parameters as outptuts and simulation parameters as inputs. The functions runs N simulations with random parameters and a fixed gain matrix diag(10.^x).

### AMB_L_8_vRTHS_simulation.sim (Simulink model)

Simulink model for vRTHS problem.

### AMB_L_9_calibration_simulation.sim (Simulink model)

Simulink model for the calibration simulations.

### Others

ElCentroAccelNoScaling.mat (Matlab file). Acceleration and time data of El Centro earthquake.

KobeAccelNoScaling.mat (Matlab file). Acceleration and time data of Kobe earthquake.

MauleAccelNoScaling.mat (Matlab file). Acceleration and time data of El Maule earthquake.

## Instructions

vRTHS with the predefined gain matrix:
1. Run *AMB_L_1_vRTHS.m*. (The adaptive gain matrix can be modified by the user)

To recalculate gain matrix:
1. Run *AMB_L_2_calibration.m* to define the calibration simulation.
2. (Optional) Run *AMB_L_3_sensitivy_inputs.m* to analyze the calibration simulation parameters.
3. (Optional) Run *AMB_L_4_sensitivity_gains.m* to analyze random gain matrices.
4. Run *AMB_L_5_optimization.m* to obtain a new adaptive gain matrix and the color maps for R2.
5. Update the variable adaptivegain in *AMB_L_1_vRTHS.m* and run to prove the new gain matrix.



