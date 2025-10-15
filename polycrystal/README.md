# Polycrystal

## Problem definition

This test case illustrates the simulation of a Representative Volume Element (`RVE`) of a polycrystal made of uranium dioxide (`UO₂`). The objective is to study the mechanical response of the material under an uniaxial loading.

In addition to the mechanical analysis, this example demonstrates how to set up a fixed-point algorithm to handle the nonlinearities associated with crystalline plasticity at the grain scale.

- Boundary conditions: periodic boundary conditions are applied on the RVE faces. The loading is imposed in one direction, ensuring compatibility and equilibrium across periodic faces.
- Constitutive law: UO₂ crystalline plasticity law[^2].
- Finite element order: 1 (linear interpolation).
- Finite element space: H1.
- Simulation duration: 200 s.
- Number of time steps: 600.
- [Crystal] Constitutive law: UO₂ crystalline plasticity law[^2].
	- Young Modulus = 222.e9
	- Poisson ratio = 0.27
	- Shear Modulus = 54.e9 

## Mesh generation
This section explains how to generate a sample mesh with Merope.

Before running the script, make sure that the environment variable `MEROPE_DIR` is properly loaded:

Then, you can generate the mesh in two steps:


```
source ${MEROPE_DIR}/Env_Merope.sh
python3 mesh/5crystals.py # generate 5crystals.geo
gmsh -3 5crystals.geo # generate 5crystals.msh
```

You will obtain a 3D mesh (5crystals.msh) of a polycrystalline sample with 5 grains.

The geometry of the RVE is generated using the Mérope [^1] toolkit. For this example, a cermet with 5 crystals is built with a metalic interface.

Make sure to load the `MEROPE` environment before running the mesh generation script:

```
source ${MEROPE_DIR}/Env_Merope.sh
python3 mesh/5crystals.py   # generates 5crystals.geo
gmsh -3 5crystals.geo       # generates 5crystals.msh
```

### Mesh generation options

The following parameters are set in the `5crystals.py` script:

```
L = [1, 1, 1]        # Dimensions of the RVE box
nbSpheres = 5        # Number of grains (polycrystal composed of 5 crystals)
distMin = 0.4        # Minimum distance between sphere centers
randomSeed = 0       # Random seed for reproducibility
MeshOrder = 1        # Polynomial order of elements
MeshSize = 0.05      # Target mesh size
```

The resulting polycrystal is composed of 5 grains.

### Mesh Polycrystal composed of 5 crystals

![Polycrystal input mesh](doc/5crystalsGmsh.png)


## Simulation options


The main executable for this test case is uniaxial-polycrystal. Its command-line options are:

```
Usage: ./uniaxial-polycrystal --help 
```

| Option                                      | Type   | Default                      | Description                                                                                                                                   |
| ------------------------------------------- | ------ | ---------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `-h, --help`                                | —      | —                            | Print the help message and exit.                                                                                                              |
| `-m <string>, --mesh <string>`              | string | `mesh/5crystals.msh`         | Mesh file to use.                                                                                                                             |
| `-f <string>, --vect <string>`              | string | `mesh/vectors_5crystals.txt` | Vector file to use.                                                                                                                           |
| `-l <string>, --library <string>`           | string | `src/libBehaviour.so`        | Material library.                                                                                                                             |
| `-b <string>, --behaviour <string>`         | string | `Mono_UO2_Cosh_Jaco3`        | Mechanical behaviour.                                                                                                                         |
| `-o <int>, --order <int>`                   | int    | `1`                          | Finite element order (polynomial degree).                                                                                                     |
| `-r <int>, --refinement <int>`              | int    | `0`                          | Refinement level of the mesh (default = 0).                                                                                                   |
| `-v <int>, --verbosity-level <int>`         | int    | `0`                          | Verbosity level of the output.                                                                                                                |
| `-d <double>, --duration <double>`          | double | `200`                        | Duration of the simulation.                                                                                                                   |
| `-n <int>, --nstep <int>`                   | int    | `600`                        | Number of simulation steps.                                                                                                                   |
| `--linear-solver <string>`                  | string | `HyprePCG`                   | Linear solver to be used.                                                                                                                     |
| `--linear-solver-preconditioner <string>`   | string | `HypreBoomerAMG`             | Preconditioner for the linear solver. Use `none` to disable.                                                                                  |
| `--macroscopic-stress-output-file <string>` | string | `uniaxial-polycrystal.res`   | Output file containing the evolution of the diagonal components of the deformation gradient and the diagonal components of the Cauchy stress. |
| `--enable-post-processings`                 | bool   | `false`                      | Execute post-processing steps (default).                                                                                                      |
| `--disable-post-processings`                | bool   | —                            | Do not execute post-processing steps.                                                                                                         |
| `--enable-export-von-Mises-stress`          | bool   | `false`                      | Export the von Mises stress.                                                                                                                  |
| `--disable-export-von-Mises-stress`         | bool   | —                            | Do not export the von Mises stress (default).                                                                                                 |
| `--enable-export-first_eigen_stress`        | bool   | `false`                      | Export the first eigen stress.                                                                                                                |
| `--disable-export-first_eigen_stress`       | bool   | —                            | Do not export the first eigen stress (default).                                                                                               |

Note: To generate the grain orientation vectors, use the randomVectorGeneration tool provided in the distribution. This ensures a consistent and physically realistic initialization of crystallographic orientations.

![RVE of Polycrystal of UO2 with 5 crystals](doc/5crystals.png)

## Results && Post processings

You can run the simulation in parallel using MPI.

```
mpirun -n 16 ./uniaxial-polycrystal
```

## Check Results

By default, the simulation generates the file `uniaxial-polycrystal.res` when running:

```
mpirun -n 16 ./uniaxial-polycrystal
```

### Plot and Compare

To visualize and compare the results, run the following Python script:

```
python3 plot_polycrystal_results.py
```

This script generates a figure named: `plot_polycrystal.png`


![](doc/plot_polycrystal.png)

### Check the Values

To verify the simulation results, run the following Python script:

```
python3 check_polycrystal_restults.py
```

The expected output is: `Check PASS`.

An example of the detailed output:

```
      Time     MFEM/MGIS       CAST3M  RelDiff_% Status
0      1.0  6.041066e+07   63100000.0   4.451762     OK
1      2.0  7.737121e+07   79000000.0   2.105167     OK
2      3.0  8.327457e+07   84300000.0   1.231384     OK
3      4.0  8.583679e+07   86600000.0   0.889139     OK
4      5.0  8.730071e+07   87900000.0   0.686468     OK
..     ...           ...          ...        ...    ...
595  199.0  1.062465e+08  106000000.0  -0.231979     OK
596  199.0  1.062465e+08  106000000.0  -0.231979     OK
597  200.0  1.062661e+08  106000000.0  -0.250424     OK
598  200.0  1.062661e+08  106000000.0  -0.250424     OK
599  200.0  1.062661e+08  106000000.0  -0.250424     OK

[600 rows x 5 columns]
Check PASS.
```

This table shows the comparison between the simulated Cauchy stress values and the reference Cast3M results, along with the relative difference and a status check.

## References

[^1]: JOSIEN, Marc. Mérope: A microstructure generator for simulation of heterogeneous materials. Journal of Computational Science, 2024, vol. 81, p. 102359.
[^2]: PORTELETTE, Luc, AMODEO, Jonathan, MADEC, Ronan, et al. Crystal viscoplastic modeling of UO2 single crystal. Journal of Nuclear Materials, 2018, vol. 510, p. 635-643.
