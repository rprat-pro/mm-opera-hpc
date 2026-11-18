# Polycrystal

## Problem definition

This test case illustrates the simulation of a Representative Volume Element (RVE) of a polycrystal made of uranium dioxide (UO₂). The objective is to study the mechanical response of the material under a uniaxial loading (see the description of the boundary conditions below). In addition to the mechanical analysis, this example implements a fixed-point algorithm enabling the simulation of a uniaxial compression/tensile test with periodic boundary conditions.	

- **Boundary conditions:** periodic boundary conditions are applied on the RVE faces. The loading is imposed in one direction, ensuring compatibility and equilibrium across periodic faces. More precisely, the axial component Fzz of the macroscopic deformation gradient is imposed. The off-diagonal components of the macroscopic deformation gradient are set to zero. The components Fxx and Fyy are the unknowns; they are determined via the fixed-point algorithm imposing zero values for the components Sxx and Syy of the macroscopic Cauchy stress tensor. The main result used for verification is the evolution of the axial component Szz of the Cauchy stress as a function of time, which is a stress-strain curve as Fzz increases linearly with time.
- **Finite element order:** 1 (linear interpolation).
- **Finite element space:** H1.
- **Simulation duration:** 200 s.
- **Number of time steps:** 200.
- **[Crystal] Constitutive law:** The UO₂ crystal plasticity law used in the example is detailed in reference [^2]; the corresponding MFront file is `law/Mono_UO2_CosH_Jaco3.mfront`. In the case of uranium dioxide, the crystal symmetry is cubic and the corresponding orthotropic elastic properties used in the crystal plasticity law are:
	- Young's modulus = 222.e9 Pa
	- Poisson's ratio = 0.27
	- Shear modulus = 54.e9 Pa

The orthotropic basis of each grain is given as an input of the crystal plasticity law. This orthotropic basis is precomputed with the Euler angles of the grain. The fixed-point algorithm uses the homogenized elastic properties of the polycrystal in order to predict the displacement gradient needed to converge toward a uniaxial tensile test. These macroscopic properties are derived from the single-crystal elastic constants given above with the Kröner self-consistent estimate for an isotropic polycrystal (see `src/MacroscopicElasticMaterialProperties.cxx` in the repository).

## Mesh generation

This section explains how to generate a sample mesh with the Mérope toolkit [^1].

Before running the script, make sure that the environment variable
`MEROPE_DIR` is properly set.

Then, you can generate the mesh in two steps:

```
source ${MEROPE_DIR}/Env_Merope.sh
python3 mesh/5crystals.py # generate 5crystals.geo
gmsh -3 5crystals.geo # generate 5crystals.msh
```

You will obtain a 3D mesh (5crystals.msh) of a polycrystalline sample
with 5 grains.

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

### Mesh of the polycrystal composed of 5 crystals

![Polycrystal input mesh](doc/5crystalsGmsh.png)


## Simulation options


The main executable for this test case is uniaxial-polycrystal. Its
command-line options are:

```
Usage: ./uniaxial-polycrystal --help 
```

| Option                                      | Type   | Default                      | Description                                                                                                                                   |
| ------------------------------------------- | ------ | ---------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `-h, --help`                                | —      | —                            | Print the help message and exit.                                                                                                              |
| `-m <string>, --mesh <string>`              | string | `mesh/5crystals.msh`         | Mesh file to use.                                                                                                                             |
| `-mm <string>, --mesh_mode <string>`        | string | `FromScratch`                | Mesh read mode (`FromScratch`, or `Restart` for a split mesh).                                                                                |
| `-f <string>, --vect <string>`              | string | `mesh/vectors_5crystals.txt` | Vector file to use.                                                                                                                           |
| `-l <string>, --library <string>`           | string | `src/libBehaviour.so`        | Material library.                                                                                                                             |
| `-b <string>, --behaviour <string>`         | string | `Mono_UO2_Cosh_Jaco3`        | Mechanical behaviour.                                                                                                                         |
| `-o <int>, --order <int>`                   | int    | `1`                          | Finite element order (polynomial degree).                                                                                                     |
| `-r <int>, --refinement <int>`              | int    | `0`                          | Refinement level of the mesh (default = 0).                                                                                                   |
| `-v <int>, --verbosity-level <int>`         | int    | `0`                          | Verbosity level of the output.                                                                                                                |
| `-d <double>, --duration <double>`          | double | `200`                        | Duration of the simulation.                                                                                                                   |
| `-n <int>, --nstep <int>`                   | int    | `200`                        | Number of simulation steps.                                                                                                                   |
| `--linear-solver <string>`                  | string | `HyprePCG`                   | Linear solver to be used.                                                                                                                     |
| `--linear-solver-preconditioner <string>`   | string | `HypreBoomerAMG`             | Preconditioner for the linear solver. Use `none` to disable.                                                                                  |
| `--macroscopic-stress-output-file <string>` | string | `uniaxial-polycrystal.res`   | Output file containing the evolution of the diagonal components of the deformation gradient and the diagonal components of the Cauchy stress. |
| `--enable-post-processing`                  | bool   | `true`                       | Execute post-processing steps (default).                                                                                                      |
| `--disable-post-processing`                 | bool   | —                            | Do not execute post-processing steps.                                                                                                         |
| `--enable-export-von-Mises-stress`          | bool   | `false`                      | Export the von Mises stress.                                                                                                                  |
| `--disable-export-von-Mises-stress`         | bool   | —                            | Do not export the von Mises stress (default).                                                                                                 |
| `--enable-export-first_eigen_stress`        | bool   | `false`                      | Export the first principal stress.                                                                                                            |
| `--disable-export-first_eigen_stress`       | bool   | —                            | Do not export the first principal stress (default).                                                                                           |

Note: To generate the grain orientation vectors, use the randomVectorGeneration tool provided in the distribution. This ensures a consistent and physically realistic initialization of crystallographic orientations.

![RVE of Polycrystal of UO2 with 5 crystals](doc/5crystals.png)

## Results and post-processing

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

This script generates a figure named `plot_polycrystal.png`, shown below. In this figure we can observe that there is a good agreement between Cast3M and MFEM-MGIS results. Although Cast3M is our reference, the stress-strain curve shows an unexpected non-smooth evolution. These oscillations are associated with the time step discretization which was specified in the MFEM-MGIS simulation. For a larger number of time steps the Cast3M curve is smoother, like the MMM results shown below. A first analysis of this issue led us to conclude that the implicit formulation of the Cast3M code, based on a quasi-Newton algorithm using the elastic stiffness and line search type convergence acceleration, is not fully stable. The main conclusion of this work is that the implicit formulation of MFEM-MGIS (full Newton algorithm using the tangent stiffness) is very efficient, thanks to the quadratic convergence and parallelization, and provides a high-quality solution. For the verification process, the number of time steps has been significantly increased in order to avoid as much as possible the oscillations of the Cast3M solution.      

![](doc/plot_polycrystal.png)

### Check the Values

To verify the simulation results, run the following Python script:

```
python3 check_polycrystal_results.py
```

The expected output is: `Check PASS`.

An example of the detailed output:

```
      Time     MFEM/MGIS       CAST3M  RelDiff_% Status
0      1.0  6.041066e+07   63121000.0   4.486524     OK
1      2.0  7.737121e+07   78996000.0   2.099997     OK
2      3.0  8.327457e+07   84282000.0   1.209769     OK
3      4.0  8.583679e+07   86623000.0   0.915934     OK
4      5.0  8.730071e+07   87898000.0   0.684177     OK
..     ...           ...          ...        ...    ...
195  196.0  1.061861e+08  105790000.0  -0.373003     OK
196  197.0  1.062065e+08  105800000.0  -0.382700     OK
197  198.0  1.062266e+08  105800000.0  -0.401578     OK
198  199.0  1.062465e+08  105810000.0  -0.410809     OK
199  200.0  1.062661e+08  105820000.0  -0.419810     OK

[200 rows x 5 columns]
Check PASS.
```

This table shows the comparison between the simulated Cauchy stress values and the reference Cast3M results, along with the relative difference and a status check.

## References

[^1]: JOSIEN, Marc. Mérope: A microstructure generator for simulation of heterogeneous materials. Journal of Computational Science, 2024, vol. 81, p. 102359.
[^2]: PORTELETTE, Luc, AMODEO, Jonathan, MADEC, Ronan, et al. Crystal viscoplastic modeling of UO2 single crystal. Journal of Nuclear Materials, 2018, vol. 510, p. 635-643.
