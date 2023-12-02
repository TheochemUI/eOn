import matplotlib.pyplot as plt
import re
from matplotlib import offsetbox

# Step 1: Read the file and extract time values
time_values = []
with open("client_spdlog.log", "r") as file:
    for line in file:
        match = re.search(r"optimize time: (\d+\.\d+)s", line)
        if match:
            time_values.append(float(match.group(1)))

# Step 2: Prepare data for plotting
cumulative_time = sum(time_values)
occurrences = list(range(1, len(time_values) + 1))

# Step 3: Plot the data and add annotations
plt.plot(occurrences, time_values, marker='o')
# Step 4: Add annotations and adjust positions
# for i, time in enumerate(time_values):
#     annotation = plt.annotate(f'{time}s',
#                               (occurrences[i], time_values[i]),
#                               xytext=(5, 5),  # Slight offset
#                               textcoords='offset points')
#     # Adjust the position of the annotation to avoid overlap
#     offsetbox.AnnotationBbox(annotation, (occurrences[i], time_values[i]))

plt.xlabel('Occurrence')
plt.ylabel('Time (s)')
plt.title(f'Optimize Time per Occurrence - Cumulative Time: {cumulative_time:.2f}s')
plt.grid(True)
plt.show()
