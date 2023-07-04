#!/usr/bin/env python

import matplotlib.pyplot as plt

lines = open("saddlesearch.dat", 'r').readlines()[1:]
all_curvatures = []
all_indices = []
translation_curvatures = []
translation_indices = []
index = 0
for line in lines:
    if "ROT" in line:
        index += 1
        all_curvatures.append(float(line.split()[5]))
        all_indices.append(index)
    else:
        translation_curvatures.append(float(line.split()[5]))
        translation_indices.append(index)
plt.plot(all_indices, all_curvatures, 'b')
plt.plot(translation_indices, translation_curvatures, 'ro')
plt.ylabel('curvature')
plt.show()
