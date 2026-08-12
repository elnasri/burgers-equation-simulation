# 1D Viscous Burgers' Equation Simulation
> **Note:** Project completed in March 2026.

Simulating the 1D viscous Burgers' equation in Python using finite differences, then plotting it against the exact solution to see how well the math holds up.

Burgers' equation combines two main physics behaviors: non linear convection (wave steepening) and diffusion (smoothing out gradients). 

* **Numerical Model:** I used a standard 1st order upwind scheme for the advection term $u \frac{\partial u}{\partial x}$ and a 2nd order central difference for the viscosity term $\nu \frac{\partial^2 u}{\partial x^2}$ with periodic boundaries.
* **Cole-Hopf?** Non linear PDEs are unusually hard to solve analytically. The Cole-Hopf transformation converts the burgers' into a linear heat equation, which happens to have a known exact solution!!.Then Using `sympy` to handle that transformation gives a clean ground truth solution to check the finite difference code against.

## Setup & Running

```bash
pip install numpy sympy matplotlib
python burgers_simulation.py
