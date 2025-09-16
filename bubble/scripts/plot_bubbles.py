#!/usr/bin/env python3
import re
import matplotlib.pyplot as plt
import sys
import numpy as np

# Define the input file
file_path = 'bubbles_and_stresses.txt'

if len(sys.argv)>1:
    file_path = sys.argv[1]

bubbles_broken = []
principal_stresses = []

bubbles_broken_re = re.compile(r"(\d+) bubbles are broken")
stress_re = re.compile(r"Value of the first principal stress: ([\d\.]+) at coordinate \(([\d\.\-\s,]+)\)")

# Read the file and extract data
with open(file_path, 'r') as file:
    for line in file:
        # Match number of bubbles broken
        bubbles_match = bubbles_broken_re.search(line)
        if bubbles_match:
            bubbles_broken.append(int(bubbles_match.group(1)))

        # Match principal stress
        stress_match = stress_re.search(line)
        if stress_match:
            principal_stresses.append(float(stress_match.group(1)))

# Plot bubbles broken as a function of principal stress
plt.figure(figsize=(8, 6))
plt.plot(principal_stresses, np.array(bubbles_broken)/max(bubbles_broken), marker='o', linestyle='-', color='b', label='Bubbles Broken')
plt.xlabel('Principal Stress', fontsize=14)
plt.ylabel('Number of Bubbles Broken', fontsize=14)

plt.grid(True)
plt.legend(fontsize=12)
plt.tight_layout()

# Show the plot
plt.savefig('plot.png')
