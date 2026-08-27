import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# Cole-Hopf analytical 

x_sym, t_sym, nu_sym = sp.symbols("x t nu")
phi_sym = (
    sp.exp(-(x_sym - 4 * t_sym) ** 2 / (4 * nu_sym * (t_sym + 1)))
    + sp.exp(-(x_sym - 4 * t_sym - 2 * sp.pi) ** 2 / (4 * nu_sym * (t_sym + 1)))
)
u_expr = -2 * nu_sym * (phi_sym.diff(x_sym) / phi_sym) + 4
u_func = sp.lambdify((t_sym, x_sym, nu_sym), u_expr, modules="numpy")


def solve_burgers(nx, nu_vis=0.5, t_final=0.5):
    dx = 2 * np.pi / (nx - 1)
    x_grid = np.linspace(0, 2 * np.pi, nx)
    u = np.array([u_func(0, x0, nu_vis) for x0 in x_grid])

    # dynamic time step 

    u_max = np.max(np.abs(u))
    dt_cfl = 0.5 * (dx / u_max)
    dt_diff = 0.4 * (dx**2 / nu_vis)

    dt = min(dt_cfl, dt_diff)
    nt = int(np.ceil(t_final / dt))
    dt = t_final / nt

    u_history = [u.copy()]

    for _ in range(nt):
        u_old = u.copy()

        # Interior update
        u[1:-1] = (
            u_old[1:-1]
            - u_old[1:-1] * dt / dx * (u_old[1:-1] - u_old[:-2])
            + nu_vis * dt / dx**2 * (u_old[2:] - 2 * u_old[1:-1] + u_old[:-2])
        )

        # Periodic boundary 
        u[0] = (
            u_old[0]
            - u_old[0] * dt / dx * (u_old[0] - u_old[-2])
            + nu_vis * dt / dx**2 * (u_old[1] - 2 * u_old[0] + u_old[-2])
        )
        u[-1] = u[0]
        u_history.append(u.copy())

    u_history = np.array(u_history)
    t_grid = np.linspace(0, t_final, nt + 1)
    l2_metric = np.sqrt(np.mean(u**2))

    return l2_metric, x_grid, t_grid, u_history, dt, nt


if __name__ == "__main__":
    nu_vis_target = 0.5
    t_end = 0.5

    # convergence study

    N = [400, 200, 100]
    results = [solve_burgers(n, nu_vis=nu_vis_target, t_final=t_end) for n in N]
    phi1, phi2, phi3 = [res[0] for res in results]

    r21 = (N[0] - 1) / (N[1] - 1)
    e21, e32 = phi2 - phi1, phi3 - phi2
    p = abs(np.log(abs(e32 / e21)) / np.log(r21))
    phi_ext = (r21**p * phi1 - phi2) / (r21**p - 1)
    gci_fine = (1.25 * abs((phi1 - phi2) / phi1)) / (r21**p - 1)

    print(f"GCI Study (nu={nu_vis_target}): p={p:.4f}, Fine GCI={gci_fine*100:.4f}%")

    # 2D comparison plot
    x1, _, u_hist1, dt1, nt1 = results[0][1:]
    u_final_num = u_hist1[-1, :]
    u_analytical = np.array([u_func(t_end, xi, nu_vis_target) for xi in x1])

    plt.figure(figsize=(11, 7), dpi=100)
    plt.plot(x1, u_final_num, "o-", markersize=3, lw=1.5, label="Numerical solution")
    plt.plot(x1, u_analytical, lw=2, label="Exact solution")
    plt.xlim(0, 2 * np.pi)
    plt.xlabel("x")
    plt.ylabel("u(x, t)")
    plt.legend()
    plt.tight_layout()

    plt.savefig(f"burgers_nx{N[0]}.png", dpi=150, bbox_inches="tight")
    plt.show()

    dx1 = 2 * np.pi / (N[0] - 1)
    print(f"nx={N[0]}, nt={nt1}, dx={dx1:.4f}, dt={dt1:.4f}, nu={nu_vis_target:.4f}")