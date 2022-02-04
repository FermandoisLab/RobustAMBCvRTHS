# Adaptive model-based compensation for non-linear case

## Description

In this example the adaptive model-based compensator is applied to the RTHS benchmark with a modified experimental substructure. Bouc-wen with degradation models are implemented to consider complex non-linear behavior in the control plant. In this example, the goal is to prove the adaptive gain matrix calibrated for linear cases but now applied to non-linear case.

## Files 

### AMB_NL_1_vRTHS.m  (Matlab script)

In this script the problem is defined and the virtual RTHS can be run.

Definitions:
  - Simulation parameters.
  - Transfer system dynamics.
  - Compensator: initial parameters, filters, adaptive gain matrix (pre-calibrated matrix is provided).
  - Earthquake.
  - Reference structure.
  - Experimental substructure (Bouc-Wen parameters).
  - Numerical substructure.
Simulation:
  - Displacement plots.
  - Force plots.
  - Adaptive parameters plots.

### AMB_NL_2_simulation.sim (Simulink model)

Simulink model for vRTHS problem.

### Others

ElCentroAccelNoScaling.mat (Matlab file). Acceleration and time data of El Centro earthquake.

KobeAccelNoScaling.mat (Matlab file). Acceleration and time data of Kobe earthquake.

### Freq_Resp_Tong.m (Matlab function)

Matlab function wiht frequency evaluation index FEI (Guo et al 2014)

## Instructions

1. Run *AMB_NL_1_vRTHS.m*.




