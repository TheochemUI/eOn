import matplotlib.pyplot as plt

# Initialize variables
scg_counts = []
current_count = 0
in_section = False

# Read and process the file
with open('eon_orca_gpneb', 'r') as file:
    for line in file:
        if "optimize time:" in line:
            if in_section:
                # End of current section, save count
                scg_counts.append(current_count)
                current_count = 0
            in_section = True
        elif "SCG iteration:" in line and in_section:
            current_count += 1

# Add the last section if file doesn't end with "optimize time:"
if current_count > 0:
    scg_counts.append(current_count)

# Create histogram
plt.hist(scg_counts, bins=len(set(scg_counts)), edgecolor='black')
plt.xlabel('SCG Iteration Count')
plt.ylabel('Frequency')
plt.title('Histogram of SCG Iterations per Optimization Time Section')
plt.show()

