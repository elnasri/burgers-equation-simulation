import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from burgers_simulation import solve_burgers

N1 = 400
nu_vis_target = 0.5
t_end = 0.5

_, x1, t1, u_hist1, _, _ = solve_burgers(N1, nu_vis=nu_vis_target, t_final=t_end)

fig = plt.figure(figsize=(11, 8), dpi=100)
fig.patch.set_facecolor("black")

ax = fig.add_subplot(111, projection="3d")
ax.set_facecolor("black")

X, T = np.meshgrid(x1, t1)

ax.plot_surface(
    X,
    T,
    u_hist1,
    cmap="jet",
    edgecolor="none",
    antialiased=True,
    rstride=1,
    cstride=1,
)

ax.xaxis.pane.set_facecolor((0, 0, 0, 1))
ax.yaxis.pane.set_facecolor((0, 0, 0, 1))
ax.zaxis.pane.set_facecolor((0, 0, 0, 1))

for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
    axis.pane.set_edgecolor("none")
    axis._axinfo["grid"]["color"] = (0.3, 0.3, 0.3, 0.5)

ax.set_xlabel("")
ax.set_ylabel("")
ax.set_zlabel("")

ax.xaxis.line.set_color("red")
ax.xaxis.set_tick_params(colors="red")

ax.yaxis.line.set_color("dodgerblue")
ax.yaxis.set_tick_params(colors="dodgerblue")

ax.zaxis.line.set_color("white")
ax.zaxis.set_tick_params(colors="white")

legend_elements = [
    Line2D([0], [0], color="red", lw=2, label="x axis"),
    Line2D([0], [0], color="dodgerblue", lw=2, label="t axis"),
    Line2D([0], [0], color="white", lw=2, label="u axis"),
]

ax.legend(
    handles=legend_elements,
    loc="upper left",
    facecolor="black",
    edgecolor="gray",
    fontsize=12,
    labelcolor="white",
)

ax.view_init(elev=31, azim=-60)
ax.invert_yaxis()

plt.savefig(
    f"burgers_nu{nu_vis_target}_nx{N1}.png",
    dpi=150,
    bbox_inches="tight",
    pad_inches=0.4,
    facecolor=fig.get_facecolor(),
)

plt.show()