# AdaptiveCompensation_CGGalmez

## Description

This repository contains five applications of adaptive dynamic compensators for virtual real-time hybrid simulation problems. The adaptive compensators corresponds to two variations of adaptive feedforward controllers. The Adaptive model-based compensator (AMB) consists in the inverse of the control plant in continuous form, but implemented with finite difference to achieve a causal controller. The AMB utilizes gradient adaptive law to identify the control plant model parameters and therefore requires the calibration of an adaptive gain matrix Gamma. On the other hand, discrete adaptive model-based compensator (dAMB) is also initially formulated from a model of the control plant, but the compensator is defined as a FIR filter and recursive least square error algorithm is utilized to update the FIR filter directly.

<img src="figures/ControlArchitecture.jpg" alt="Compensation" width="800"/>

## Requirements

- Matlab R2021a or superior

## Folders description

A brief description of each application is presented and description of each file and instructions to execute are provided in the readme file inside the respective folders.

### 1. AMB_Linear

In this example the adaptive model-based compensator is designed, calibrated and applied to the RTHS benchmark. A calibrated gain matrix Gamma is provided to run several simulations with linear systems includding uncertainty. However, the adaptive gain matrix Gamma can be recalculated and analized using the calibration procedure.

### 2. AMB_NonLinear 

In this example the adaptive mode-based compensator is applied to a modification of the RTHS benchmark, where the linear experimental substructure is replaced by non-linear models. Modified Bouc-Wen models are utilized to consider hysteretic models with degradation. This example does not include gains calibration since the goal is to prove the gains obtained for the linear case.

### 3. dAMB_NonLinear

In this example the discrete adaptive model-based compensator is applied to the RTHS benchamark with non-linear experimental substructure. In this case, degradation is not considered because for dynamic compensation the hysteretic case without degradation is more difficult to compensate. For the dAMB compensator, a .pdf file is provided with the formulation and equations.

### 4. dAMB+Feedback 

In this example the discrete adaptive model-based compensator is combined with a feedback controller to improve the tracking in pressence of non-linearities in the control plant. The vRTHS problem consists in the RTHS benchmark with non-linear experimental substructure. To design the feedback controller, parametric and non parametric uncertainty is considered in the control plant. Loop shaping and nyquist plots are utilized to compare alternatives of linear feedback controllers. The feedback controllers improves the tracking for the non-linear case, but a bad design of feedback controller can deteriorate tracking or even can cause instability.

### 5. maRTHS with dAMB

In this example a multi-axial RTHS problem is defined. Consists in 4DOF frame where one column corresponds to the experimental substructure with displacement and rotation as DOF to be imposed. The servo-hydraulic actuators are modeled to considere actuator-specimen interaction and actuator-actuator interaction. The compensation consists in a decentralized compensation, where displacement and rotation are transformed into local actuator target displacements using kinematic transformations. Each actuator is compensated using discrete-adaptive model based compensation. This example provide evidence of complex issues for maRTHS such as stability problems and actuators interaction.

### 6. dAMB (compensator only)

This folder contains the discrete adaptive model based compensator code and simulink block.
