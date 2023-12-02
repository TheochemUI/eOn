import re
import numpy as np
from pathlib import Path


def extract_interp_points(data):
    # Regular expression to find the energy type and the two-column data following it
    # https://regex101.com/r/RF6b4V/2
    # fmt: off
    pattern = (
        r".*?Iteration:.*?"
    ) + THREE_COL_NUM
    matchr = re.search(pattern, data, re.MULTILINE)
    # fmt: on

    if not matchr:
        return np.array([]), np.array([])

    # Extract and convert the data
    interpdata = matchr.group("threecolnum")
    idx_values, x_values, y_values = [], [], []
    for line in interpdata.split("\n"):
        idx, x, y = map(float, line.split())
        x_values.append(x)
        y_values.append(y)
        idx_values.append(idx)

    return np.array(x_values), np.array(y_values)

filename = "orca.interp"
interp_data = Path(filename).open().read()
print(extract_interp_points(interp_data))


# plt.style.use("bmh")
# # Set up your figure and axes
# fig, ax = plt.subplots(figsize=(3.2, 2.5), dpi=200)

# # Get the colormap and reverse it
# cmap = cm.batlow  # sequential
# # cmap = cm.roma # diverging

# # Iterate over iterations
# num_iterations = len(iterations_data)
# for idx, (iteration, data) in enumerate(iterations_data.items()):
#     data = np.array(data)
#     energy = data[:, 1]
#     distance = data[:, 0]

#     # Cubic spline interpolation
#     distance_fine = np.linspace(distance.min(), distance.max(), num=100)
#     spl = splrep(distance, energy, k=3)
#     spl_y = splev(distance_fine, spl)

#     # Color specification using colormap
#     color = cmap(1.0 * (idx / (num_iterations - 1)))

#     # Alpha specification
#     alpha = 1 if idx == 0 or idx == num_iterations - 1 else 0.5

#     # Plotting
#     if idx == num_iterations - 1:
#         ax.plot(distance_fine, spl_y, color="red", alpha=alpha)
#         ax.plot(distance, energy, linestyle="", marker="o", color="red", alpha=alpha)
#     else:
#         ax.plot(distance_fine, spl_y, color=color, alpha=alpha)
#         ax.plot(distance, energy, linestyle="", marker="o", color=color, alpha=alpha)

# # Labels and adjustments
# ax.set_xlabel("Distance (Bohr)")
# ax.set_ylabel("Energy (Eh)")
# ax.minorticks_on()
# ax.set_ylim(0, ax.get_ylim()[1])
# ax.set_xlim(0, ax.get_xlim()[1])
# ax.set_facecolor("gray")
# plt.title("NEB Paths")
# plt.tight_layout(pad=0.2)
# plt.subplots_adjust(left=0.20)

# # Adding colorbar
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=num_iterations - 1))
# cbar = fig.colorbar(sm, ax=ax)
# cbar.set_label("Iteration Index")

# # Save and show the plot
# plt.savefig("iteration_energy_path_combined.pdf", transparent=True)
# plt.show()
