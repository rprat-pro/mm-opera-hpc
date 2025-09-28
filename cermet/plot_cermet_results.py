import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Lecture du fichier .res (colonnes séparées par espaces) ---
data_res = np.loadtxt("cermet.res")  # par défaut espace = séparateur

x_res = data_res[:, 0]  # 1ère colonne
y_res = data_res[:, 6]  # 7ème colonne (index 6 car Python commence à 0)

# --- Lecture du fichier CSV ---
data_csv = pd.read_csv("results/cermet_cast3m.csv", sep=";")
x_csv = data_csv.iloc[:, 0]  # 1ère colonne
y_csv = data_csv.iloc[:, 6]  # 7ème colonne

# --- Tracé ---
plt.figure(figsize=(8, 6))
plt.plot(x_csv, y_csv, color="orange", label="Cast3m", linestyle="-")
plt.plot(x_res, y_res, color="blue", linestyle='-', label="MFEM/MGIS")

plt.xlabel("Time (s)")
plt.ylabel("Szz")
plt.title("Cermet test case with 5 crystals")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plot_cermet.png")
