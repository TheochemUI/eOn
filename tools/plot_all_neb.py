import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splrep, splev
import glob
from cmcrameri import cm

plt.style.use("bmh")

# List all files matching the pattern 'neb_***.dat', excluding the last one
file_paths = glob.glob("neb_*.dat")
file_paths.sort()  # Sort the file paths to ensure correct order

# Set up your figure and axes
fig, ax = plt.subplots(figsize=(3.2, 2.5), dpi=200)

# Get the colormap and reverse it
cmap = cm.batlow  # sequential
# cmap = cm.roma # diverging

# Iterate over file_paths
num_files = len(file_paths)
for idx, file_path in enumerate(file_paths):
    data = np.loadtxt(file_path).T
    energy = data[2]
    rc = data[1]

    # Cubic spline interpolation
    rc_fine = np.linspace(rc.min(), rc.max(), num=100)
    spl = splrep(rc, energy, k=3)
    spl_y = splev(rc_fine, spl)

    # Color specification using colormap
    color = cmap(
        1.0 * (idx / (num_files - 1))
    )  # Normalize index to get a color from a fraction of the colormap

    # Alpha specification
    alpha = 1 if idx == 0 or idx == num_files - 1 else 0.5

    # Plotting
    if idx == num_files - 1:
        ax.plot(rc_fine, spl_y, color="red", alpha=alpha)
        ax.plot(rc, energy, linestyle="", marker="o", color="red", alpha=alpha)
    else:
        ax.plot(rc_fine, spl_y, color=color, alpha=alpha)
        ax.plot(rc, energy, linestyle="", marker="o", color=color, alpha=alpha)

# Labels and adjustments
ax.set_xlabel("Reaction Coordinate, ($\AA$)")
ax.set_ylabel("Energy, ($eV$)")
ax.minorticks_on()
ax.set_ylim(0, ax.get_ylim()[1])
ax.set_xlim(0, ax.get_xlim()[1])
ax.set_facecolor("gray")
plt.title("NEB paths over optimization steps")
plt.tight_layout(pad=0.2)
plt.subplots_adjust(left=0.20)

# Adding colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=num_files - 1))
cbar = fig.colorbar(sm, ax=ax)
cbar.set_label("File Index")

# Save and show the plot
plt.savefig("neb_energy_path_combined.pdf", transparent=True)
plt.show()
