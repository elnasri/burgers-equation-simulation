# 1D Viscous Burgers' Equation Simulation

Simulating 1D viscous Burgers' equation in python using finite differences, then plotting it against the exact solution to see how well the math holds up.

# Visualizations

<p align="center">
  <img src="burgers_nu0.5_nx400.png" alt="3D plot of smooth burgers equation solution" width="55%">
  <br>
  <em>Figure 1: A Smooth 3D surface plot of the Burgers' equation solution field u(x,t) for ν = 0.5 and Nx = 400.</em>
</p>

<p align="center">
  <img src="burgers_nx200-3D.png" alt="3D plot of burgers equation" width="55%">
  <br>
  <em>Figure 2: 3D evolution of u(x,t) showing non linear wave steepening into a sharp viscous shock ν = 0.07.</em>
</p>

<br>

<p align="center">
  <img src="burgers_nx200.png" alt="2D numerical vs analytical comparison" width="59%">
  <br>
  <em>Figure 3: Validation of the finite difference solver against the analytical solution via Cole-Hopf. </em>
</p>




## Overview




Burgers' equation combines two main physical behaviors: nonlinear convection (the wave steepening) and diffusion (smoothing of gradient). Burgers' equation is classified as a parabolic PDE, and is solved here using finite differences

* **Numerical Model:** Discretised using a 1st-order backward upwind scheme for advection ($u \frac{\partial u}{\partial x}$) and a 2nd-order central difference for diffusion ($\nu \frac{\partial^2 u}{\partial x^2}$) with periodic boundaries. Backward differencing is used because $u(x,t) > 0$ across the domain.

* **Analytical Validation via the Cole Hopf Transformation:** Non linear PDEs are unusually difficult to solve analytically. The Cole-Hopf transformation reduces the nonlinear Burgers' equation to the linear heat equation, which happens to have a known closed-form solution. Then using `sympy` to handle that transformation gives a clean ground truth solution to check the finite difference code against.


## Grid Convergence Study

I Ran a 3 grid GCI check (N = 400, 200, 100) on the higher viscosity solver (`burgers_higher_viscosity_simulation.py`, ν = 0.5) to check if the solver is actually converging the way it should.
Here, the solution variable $\phi$ represents the global velocity $L_2$ norm: $\phi = \sqrt{\frac{1}{N}\sum u_i^2}$


| Mesh Level | Grid Size ($N$) | Solution ($\phi$) |
| :--- | :--- | :--- |
| **Fine ($N_1$)** | 400 | 4.094330 |
| **Medium ($N_2$)** | 200 | 4.081411 |
| **Coarse ($N_3$)** | 100 | 4.057673 |

* **Order of Convergence ($p$):** 0.8778
* **Extrapolated Limit ($\phi_{\text{ext}}$):** 4.109753
* **Fine Grid GCI ($GCI_{\text{fine}}$):** 0.4709%

$p \approx 0.88$, close to 1st order. This matches theoretical expectations since advection uses a 1st-order upwind scheme while diffusion uses a 2nd-order central difference scheme, causing the lower-order advection error to dominate. The GCI on the finest grid is well under 1%.
 
 ## Flow Regimes

Varying viscosity ($\nu$) alters the balance between advection and diffusion, represented by the Reynolds number:

* **High Re ($Re \approx 90, \nu = 0.07$):** Low diffusion allows nonlinear wave steepening into a sharp viscous shock (figures 2 & 3)
* **Low Re ($Re \approx 12.5, \nu = 0.5$):** High diffusion smooths out gradients over time, providing a stable baseline for the GCI study (figure 1)


## Setup & Running

```bash
pip install numpy sympy matplotlib
python burgers_simulation.py                         # Low viscosity shock case (v = 0.07)
python burgers_higher_viscosity_simulation.py         # High viscosity case (v = 0.5) with GCI study
