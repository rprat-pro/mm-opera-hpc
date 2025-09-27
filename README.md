# mm-opera-hpc

This directory groups together the various test cases implemented as part of the OperaHPC project.

## Installation using Spack [recommended]

```
git clone --depth=2 --branch=v1.0.1 https://github.com/spack/spack.git
export SPACK_ROOT=$PWD/spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
```

Firstly, get the mfem-mgis spack repository.

```
git clone https://github.com/rprat-pro/spack-repo-mfem-mgis.git
spack repo add $PWD/spack-repo-mfem-mgis
```

Secondly, install mfem-mgis

```
spack install mfem-mgis@1.0.3
```

Thirdly, load mfem-mgis

```
spack load mfem-mgis
```

Create a build directory, configure the project with CMake, build it, and install.

```
git clone https://github.com/rprat-pro/mm-opera-hpc.git
cd mm-opera-hpc/
mkdir build && cd build
cmake .. -DCMAKE_PREFIX_PATH=`spack location -i tfel`/share/tfel/cmake -DCMAKE_INSTALL_PREFIX=../install
make -j 4
ctest
```

For more details on installing mfem-mgis, particularly for installing mfem-mgis without spack (cmake, not recommended) or without the internet, please visit: https://thelfer.github.io/mfem-mgis/installation_guide/installation_guide.html

## Bubble Case 

### Short Description

The default example is the rupture of a spherical, pressurized inclusion (e.g., a gas bubble) in an elastic infinite medium.

The criterium to determine the rupture or not is based on a simple geometrical assumption, i.e., if a certain distance $d_min$ is found between the position of the center of the sphere and the location of the maximum principal stress caused by the pressure exerted by the bubble on the surrounding matrix.

### Parameters

```
Usage: ./test-bubble [options] ...
Options:
   -h, --help
	Print this help message and exit.
   -m <string>, --mesh <string>, current value: mesh/single_sphere.msh
	Mesh file to use.
   -l <string>, --library <string>, current value: src/libBehaviour.so
	Material library.
   -f <string>, --bubble-file <string>, current value: mesh/single_bubble.txt
	File containing the bubbles.
   -o <int>, --order <int>, current value: 1
	Finite element order (polynomial degree).
   -r <int>, --refinement <int>, current value: 0
	refinement level of the mesh, default = 0
   -p <int>, --post-processing <int>, current value: 1
	run post processing step
   -v <int>, --verbosity-level <int>, current value: 0
	choose the verbosity level
```

### Verification against the analytical solution

The problem of a pressurized spherical inclusion in an infinite, elastic medium has a closed solution for the expressions of the stresses as a function of the distance from the sphere center. 

$$
\sigma_{\theta\theta}(r) = \dfrac{p_{in}*R_b^3}{2r^3}
$$

where $p_{in}$ is the internal pressure, $R_b$ the bubble radius, and reminding that the expression is holding for $r>R_b$.

The comparison obtained running the test-case considering the FE mesh available `bubble/mesh/single_sphere.msh` and the analytical solution is showed below, and can be generated thanks to the script available in the `verification/bubble` folder.

![Bubble Case](/img/bubble/comparison_analytical_mmm.png)

### Output With One Bubble

![Bubble Case](/img/bubble/bubble.png)

## Polycrystal Case

### Short Description

This example is under construction


### Parameters

```
Usage: ./uniaxial-polycristal [options] ...
Options:
   -h, --help
	Print this help message and exit.
   -m <string>, --mesh <string>, current value: periodic-cube.msh
	Mesh file to use.
   -f <string>, --vect <string>, current value: periodic-cube-vecteurs.txt
	Vector file to use.
   -l <string>, --library <string>, current value: src/libBehaviour.so
	Material library.
   -o <int>, --order <int>, current value: 2
	Finite element order (polynomial degree).
   -r <int>, --refinement <int>, current value: 0
	refinement level of the mesh, default = 0
   -p <int>, --post-processing <int>, current value: 1
	run post processing step
   -v <int>, --verbosity-level <int>, current value: 0
	choose the verbosity level
```

### Output Example

![Bubble Case](/img/polycrystal/polycrystal.png)
