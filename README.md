# LCR-Circuit-Modernization

## Overview

This repository documents the evolution of a computational physics study on Driven LCR-Series Circuits. Originally developed as a dissertation in April 2025 using procedural Fortran 90, this project is currently being modernized into a high-performance Python/NumPy framework.

The primary goal is to solve the governing second-order differential equation for an AC-driven LCR circuit using the 4th-Order Runge-Kutta (RK4) numerical method and to benchmark the performance gains of modern vectorization.

## Technical Background

The dynamic behavior of the circuit is modeled by Kirchhoff’s Voltage Law (KVL):

$$L\frac{d^2q}{dt^2} + R\frac{dq}{dt} + \frac{1}{C}q = V_0\sin(\omega t)$$

In this project, the 2nd-order ODE is decomposed into a system of two coupled 1st-order equations to allow for high-precision integration via the RK4 algorithm.

 $\frac{dq}{dt} = i$

 $\frac{di}{dt} = \frac{1}{L} \left( V_0\sin(\omega t) - Ri - \frac{q}{C} \right)$

This method captures transient responses and steady-state resonance with a validated accuracy of $10^{-6}$.

## Methodology

This project utilizes the 4th-Order Runge-Kutta (RK4) algorithm, which provides a local truncation error of $O(h^5)$ and a global error of $O(h^4)$. This allows for high-precision modeling of complex damping states:

- Underdamped: Characterized by oscillating transients before reaching steady-state resonance.

- Overdamped: A non-oscillatory return to equilibrium.

- Critically Damped: The fastest return to steady state without oscillation.

## Repository Structure

v1_legacy_fortran/
Contains the original research artifacts from the 2025 dissertation:

- Source Code: Procedural Fortran 90 implementations of the RK4 solver.

- Visualization: GnuPlot scripts used for transient and steady-state analysis.

- Documentation: The original unpublished dissertation PDF.

v2_modern_python/ (In Progress)
The current "March 2026" milestone focused on:

- Vectorized Solvers: Re-implementing the RK4 engine using NumPy for multi-parameter sweeps.  

- Benchmarking: Side-by-side execution time comparisons between compiled Fortran and interpreted (vectorized) Python.

- Advanced Analysis: 3D resonance surfaces and Phase Space ($q$ vs. $i$) mapping.

## Modernization Goals

Performance Auditing: Measuring the runtime overhead of Python vs. Fortran for large-scale simulations ($N > 10^6$ steps).

Numerical Validation: Cross-verifying custom RK4 results against ` scipy.integrate.solve_ivp. `

Phase Space Analysis: Visualizing energy dissipation and attractor stability in complex damping scenarios.

## How to Use

1. Legacy: View original Fortran logic and GnuPlot visualizations in the `v1_legacy_fortran/` directory.

2. Modern: Run the Python engine (requires NumPy and Matplotlib) found in `v2_modern_python/`.
