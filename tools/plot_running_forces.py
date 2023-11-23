import matplotlib.pyplot as plt
import numpy as np
from cmcrameri import cm
from matplotlib.colors import Normalize

# eV/A
FORCE_NORM_CONVERGENCE = 0.025711


# Function to parse data from the file
def parse_data(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()

    data = []
    switch_iteration = None
    for line in lines:
        if "iteration" in line or "----" in line or line.strip() == "":
            continue  # Skip headers and empty lines
        if "Switched to" in line:
            switch_iteration = int(lines[lines.index(line) + 1].split()[0])
            continue
        parts = line.split()
        if len(parts) >= 5:
            try:
                iteration = int(parts[0])
                force = float(parts[2])
                maxE = float(parts[-1])
                data.append((iteration, force, maxE))
            except ValueError:
                continue  # Skip lines that don't match the expected format
    return data, switch_iteration


# Parse the data
file_path = "client_spdlog.log"  # Path to your log file
data, switch_iteration = parse_data(file_path)
is_switched = False

# Check if data is empty
if not data:
    raise ValueError("No data extracted from the file. Please check the file format.")
if switch_iteration is not None:
    is_switched = True

iterations, force_norms, max_energies = zip(*data)

# Convert to numpy arrays
iterations = np.array(iterations)
force_norms = np.array(force_norms)
max_energies = np.array(max_energies)

# Normalize max_energies for color mapping
norm = Normalize(vmin=min(max_energies), vmax=max(max_energies))
normalized_max_energies = norm(max_energies)
colors = cm.batlow(normalized_max_energies)

# Set up your figure and axes
plt.style.use("bmh")
fig, ax = plt.subplots(figsize=(10, 6), dpi=100)

# Plotting
color = cm.batlow(normalized_max_energies)
scatter = ax.scatter(iterations, force_norms, c=color, label="NEB Step")

# Annotate the first and last energies
props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)
ax.text(
    iterations[0],
    force_norms[0],
    f"{max_energies[0]:.4f}",
    color="black",
    ha="right",
    va="bottom",
)
ax.text(
    iterations[-1],
    force_norms[-1],
    f"{max_energies[-1]:.4f}",
    color=color[-1],
    ha="left",
    va="bottom",
)

# Adding a vertical line to indicate the switch iteration
if is_switched:
    ax.axvline(
        x=switch_iteration,
        color="red",
        linestyle="--",
        label="Switch to LBFGS",
        alpha=0.5,
    )

ax.axhline(
    y=FORCE_NORM_CONVERGENCE,
    color="black",
    linestyle="--",
    label="Final convergence",
    alpha=0.5,
)

# Labels and title
ax.set_xlabel("Iteration Number")
ax.set_ylabel("Perpendicular Force Norm")
# ax.set_facecolor("gray")
ax.set_title("EON NEB Calculation: Force Norm per Iteration\n Text is max energy")
ax.legend(loc="upper right")

# Add the last data point's x-value
ax.set_xticks(np.append(ax.get_xticks()[1:], [iterations.max()]))

# Adding a colorbar for maximum energy
cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cm.batlow), ax=ax)
cbar.set_label("Max Energy")
plt.tight_layout()
plt.savefig("neb_force_iteration_plot.pdf")
plt.show()
