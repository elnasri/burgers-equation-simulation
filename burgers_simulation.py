import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


# Exact solution using the Cole-Hopf transformation
x, t, nu = sp.symbols("x t nu")

phi = (
    sp.exp(-(x - 4 * t) ** 2 / (4 * nu * (t + 1)))
    + sp.exp(-(x - 4 * t - 2 * sp.pi) ** 2 / (4 * nu * (t + 1)))
)

phi_x = phi.diff(x)
u_expr = -2 * nu * (phi_x / phi) + 4
u_func = sp.lambdify((t, x, nu), u_expr, modules="numpy")


# Parameters
nx = 200
nt = 100
nu_vis = 0.07

dx = 2 * np.pi / (nx - 1)
dt = dx * nu_vis

x_grid = np.linspace(0, 2 * np.pi, nx)
u = np.array([u_func(0, x0, nu_vis) for x0 in x_grid])


# Time marching
for n in range(nt):
    u_old = u.copy()

    u[1:-1] = (
        u_old[1:-1]
        - u_old[1:-1] * dt / dx * (u_old[1:-1] - u_old[:-2])
        + nu_vis * dt / dx**2
        * (u_old[2:] - 2 * u_old[1:-1] + u_old[:-2])
    )

    # Periodic boundary condition / circle back condition
    u[0] = (
        u_old[0]
        - u_old[0] * dt / dx * (u_old[0] - u_old[-2])
        + nu_vis * dt / dx**2
        * (u_old[1] - 2 * u_old[0] + u_old[-2])
    )
    u[-1] = u[0]


# Exact solution at final time
u_analytical = np.array(
    [u_func(nt * dt, xi, nu_vis) for xi in x_grid]
)


# Plot
plt.figure(figsize=(11, 7), dpi=100)
plt.plot(x_grid, u, "o-", markersize=3, lw=1.5, label="Numerical solution")
plt.plot(x_grid, u_analytical, lw=2, label="Exact solution")

plt.xlim(0, 2 * np.pi)
plt.legend()
plt.tight_layout()

output_file = f"burgers_nx{nx}.png"
plt.savefig(output_file, dpi=150, bbox_inches="tight")
plt.show()

print(
    f"nx={nx}, nt={nt}, dx={dx:.4f}, "
    f"dt={dt:.4f}, nu={nu_vis:.4f}"
)