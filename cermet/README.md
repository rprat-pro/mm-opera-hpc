# Cermet simulation

## Short description

This simulation consists of applying a tensile load on a cermet RVE. The cermet
is a polycrystal where each grain has a material ID (from 2 to Nmat – 1) and a
different orientation. A metallic interface is present between the different
polycrystals (material ID 1)

In addition to the mechanical analysis, this example demonstrates how to set up a fixed-point algorithm to handle the nonlinearities associated with crystalline plasticity at the grain scale.

Parameters:

- Boundary conditions: periodic boundary conditions are applied on the RVE faces. The loading is imposed in one direction, ensuring compatibility and equilibrium across periodic faces.
- [Crystal] Constitutive law: UO₂ crystalline plasticity law[^2].
	- Young Modulus = 222.e9
	- Poisson ratio = 0.27
	- Shear Modulus = 54.e9 
- [Crystal] Constitutive law: Norton.
	- Young Modulus = 276e+09
	- Poisson ratio = 0.3
	- A             = 2.5e+11
	- n1            = 4.75
	- Q             = 306.27e+03
	- D0            = 1.55e-5
	- b             = 2.5e-10
- Finite element order: 1 (linear interpolation).
- Finite element space: H1.
- Simulation duration: 200 s.
- Number of time steps: 500.
- Linear solver: HyprePCG (solver) /  (precond)

## Mesh generation

This section explains how to generate a sample mesh with `Merope`.


Before running the script, make sure that the environment variable `MEROPE_DIR` is properly loaded:

Then, you can generate the mesh in two steps:

```
source ${MEROPE_DIR}/Env_Merope.sh
python3 mesh/5grains.py # generate 5grains.geo
gmsh -3 5grains.geo # generate 5grains.msh
```

You will obtain a 3D mesh (5grains.msh) of a polycrystalline sample with 5 grains.

### Options

Mesh Generation Examples

The mesh can be customized by adjusting the input parameters in the Python script.
Below are two examples:

#### Small Example

This configuration generates a small test case with 5 grains.

```
L = [1, 1, 1]
nbSpheres = 20 
distMin = 0.3
randomSeed = 0
layer=0.02
MeshOrder = 1
MeshSize = 0.05
```

![Cermet Case](doc/cermet-5grains.png)

![Cermet Case png](doc/cermet-5grains-gmsh.png)

#### Large Example

This setup generates a realistic polycrystalline mesh with:

- 250 grains
- 12,913,361 nodes
- 86,213,779 elements


```
L = [5, 5, 5]
nbSpheres = 250
distMin = 0.1
randomSeed = 0
layer=0.04
MeshOrder = 1
MeshSize = 0.02
```

## Run your simulation

### Command-line Usage


```
Usage: ./cermet [options] ...
```

| Option                                      | Type   | Default               | Description                                                                                                                                                       |
| ------------------------------------------- | ------ | --------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-h, --help`                                | —      | —                     | Print the help message and exit.                                                                                                                                  |
| `-m <string>, --mesh <string>`              | string | `mesh/5grains.msh`    | Mesh file to use.                                                                                                                                                 |
| `-o <int>, --order <int>`                   | int    | `1`                   | Finite element order (polynomial degree).                                                                                                                         |
| `-r <int>, --refinement <int>`              | int    | `0`                   | Refinement level of the mesh (default = 1).                                                                                                                       |
| `-p <int>, --post-processing <int>`         | int    | `1`                   | Run the post-processing step.                                                                                                                                     |
| `-v <int>, --verbosity-level <int>`         | int    | `0`                   | Verbosity level of the output.                                                                                                                                    |
| `-d <double>, --duration <double>`          | double | `200`                 | Duration of the simulation (default = 5).                                                                                                                         |
| `-n <int>, --nstep <int>`                   | int    | `400`                 | Number of simulation steps (default = 40).                                                                                                                        |
| `-f <string>, --file <string>`              | string | `vectors_5grains.txt` | Vector file to use.                                                                                                                                               |
| `--macroscopic-stress-output-file <string>` | string | `cermet.res`          | Main output file containing:<br>• Evolution of the diagonal components of the deformation gradient<br>• Evolution of the diagonal components of the Cauchy stress |


### Run it

You can run the simulation in parallel using MPI.
Below are two examples:

#### Basic Test

Runs a short simulation with:

- Duration = 0.5 s
- 1 timestep
- Mesh = 5grains.msh
- Refinement level = 0


```
mpirun -n 12 ./cermet --duration 0.5 --nstep 1
```

#### Full Test

Runs a longer simulation with:

- Duration = 200 s
- 400 timesteps
- Refinement level = 1
- Custom mesh (yourmesh.msh)


```
mpirun -n 12 ./cermet --duration 200 --nstep 400 -r 1 --mesh yourmesh.msh
```

## Results
