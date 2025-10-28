import pyvista as pv
import numpy as np
import os
from matplotlib import pyplot as plt

cur_path = os.getcwd()
file_pvd = os.path.join(cur_path, 'TestCaseBubble/TestCaseBubble.pvd')
bubble_radius = 0.400
rev_side_length = 10
bubble_center = (bubble_radius, 0., 0.)
end_point_rev = (rev_side_length/2, 0., 0.)
bubble_pressure = 1 # 1e6

data_mmm = pv.read(file_pvd)
data_simulation = data_mmm[0]  # first timestep

data_sample = data_simulation.sample_over_line(
    pointa=bubble_center, pointb=end_point_rev, progress_bar=True )

coordinate = data_sample.points[:, 0]
stress_tt = data_sample["Stress"]
stress_tt = stress_tt[:, 1]

# Analytical solution
analytical_hoop_stress = 0.5*bubble_pressure * \
    bubble_radius**3/np.power(coordinate, 3)

plt.figure()
p1 = plt.plot(coordinate, stress_tt, '-.', label="$\sigma_{YY}$, MMM")
p2 = plt.plot(coordinate, analytical_hoop_stress,
              '--', label="$\sigma_{YY}$, analytical")
plt.xlim([bubble_radius, rev_side_length/2])
plt.xlabel("Distance from the bubble center, m")
plt.ylabel("Stress, Pa")
plt.legend()
plt.grid(True)
plt.savefig('comparison.png')
plt.close()


### Comparison

abs_tol = 0.01
rel_tol = 0.1

# Calcul de l'écart mixte rel/abs

Diff = np.abs((stress_tt - analytical_hoop_stress)) - rel_tol* np.abs(analytical_hoop_stress) - abs_tol

# Définition du statut
Status = np.where( Diff <= 0, 'OK', 'MISMATCH')

# Check if the 'Status' column contains any "MISMATCH"
if (Status == 'MISMATCH').any():
    print("There is at least one MISMATCH in the data!")
    print(Diff)
    print(Status)
else:
    print("Check PASS.")
