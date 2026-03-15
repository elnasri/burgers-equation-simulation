# Plot
fig = plt.figure(figsize=(11, 8), dpi=100)
fig.patch.set_facecolor("black")

ax = fig.add_subplot(111, projection="3d")
ax.set_facecolor("black")

X, T = np.meshgrid(x_grid, t_grid)

# Plot the solution
ax.plot_surface(
    X,
    T,
    u_history,
    cmap="jet",
    edgecolor="none",
    antialiased=True,
    rstride=1,
    cstride=1
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

# simple legend for the axes
from matplotlib.lines import Line2D

legend_elements = [
    Line2D([0], [0], color="red", lw=2, label="x axis"),
    Line2D([0], [0], color="dodgerblue", lw=2, label="t axis"),
    Line2D([0], [0], color="white", lw=2, label="u axis")
]

ax.legend(
    handles=legend_elements,
    loc="upper left",
    facecolor="white",
    edgecolor="gray",
    fontsize=12
)

ax.view_init(elev=31, azim=-60)

ax.invert_yaxis()

output_file = f"burgers_nx{nx}.png"

plt.savefig(
    output_file,
    dpi=150,
    bbox_inches="tight",
    pad_inches=0.4,
    facecolor=fig.get_facecolor()
)

plt.show()