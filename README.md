# 1D Viscous Burgers' Equation Simulation
> **Note:** completed in March 2026.

Simulating 1D viscous Burgers' equation in python using finite differences, then plotting it against the exact solution to see how well the math holds up.

# Visualizations

<p align="center">
  <img src="burgers_nx200-3D.png" alt="3D plot of burgers equation" width="59%">
  <br>
  <em>Figure 1: 3D evolution of u(x,t) showing non linear wave steepening into a sharp viscous shock ν = 0.07.</em>
</p>

<br>

<p align="center">
  <img src="burgers_nx200.png" alt="2D numerical vs analytical comparison" width="59%">
  <br>
  <em>Figure 2: Validation of the finite difference solver against the analytical solution via Cole-Hopf. </em>
</p>




## Overview




Burgers' equation combines two main physics behaviors: non linear convection (wave steepening) and diffusion (smoothing out gradients). 

* **Numerical Model:** I used a 1st order upwind scheme for the advection term $u \frac{\partial u}{\partial x}$ and a 2nd order central difference for the viscosity term $\nu \frac{\partial^2 u}{\partial x^2}$ with periodic boundaries.
* **Cole-Hopf?** Non linear PDEs are unusually hard to solve analytically. So using cole-hopf transformation converts the burgers' into a linear heat equation, which happens to have a known exact solution!!. Then Using `sympy` to handle that transformation gives a clean ground truth solution to check the finite difference code against.

## Setup & Running

```bash
pip install numpy sympy matplotlib
python burgers_simulation.py
