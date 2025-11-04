# mm-opera-hpc

This directory groups together the various test cases implemented as part of the OperaHPC project.

## Installation using Spack [recommended]

```bash
git clone --depth=2 --branch=v1.0.1 https://github.com/spack/spack.git
export SPACK_ROOT=$PWD/spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
```

Firstly, get the mfem-mgis spack repository.

```bash
git clone https://github.com/rprat-pro/spack-repo-mfem-mgis.git
spack repo add $PWD/spack-repo-mfem-mgis
```

Secondly, install mfem-mgis.

```bash
spack install mfem-mgis@1.0.3
```

Thirdly, load mfem-mgis.

```bash
spack load mfem-mgis
```

Create a build directory, configure the project with CMake, build it, and install.

```bash
git clone https://github.com/rprat-pro/mm-opera-hpc.git
cd mm-opera-hpc/
mkdir build && cd build
spack load tfel
cmake .. -DCMAKE_PREFIX_PATH=`spack location -i tfel`/share/tfel/cmake -DCMAKE_INSTALL_PREFIX=../install
make -j 4
ctest
```

For more details on installing mfem-mgis, particularly for installing mfem-mgis without spack (cmake, not recommended) or without the internet, please visit: https://thelfer.github.io/mfem-mgis/installation_guide/installation_guide.html

## Bubble Case 

The applicative study case is the elastic computation of the stress field induced by an internal pressure prescribed in spherical porosities (e.g., a fission gas bubble for fuel) distributed randomly in a periodic RVE. The boundary conditions are null macroscopic displacement gradient corresponding to a macroscopic strain equal to zero. The periodic surface of the RVE is defined according a Voronoi tessellation principle in order to avoid truncated porosities and ensure a high quality of the periodic solution near the boundary. In the general case the null periodic displacement boundary condition prevents the volume expansion with a non-zero macroscopic hydrostatic stress. 

![Bubble Case](/img/bubble/bubbles.png)

More details in [bubble/README.md](./bubble/README.md)

## Polycrystal Case

This test case illustrates the simulation of a Representative Volume Element (RVE) of a polycrystal made of uranium dioxide (`UO₂`). The objective is to study the mechanical viscoplastic response of the material under a macroscopic uniaxial loading. In addition to the mechanical analysis, this example implements a fixed-point algorithm to converge toward a macroscopic uniaxial tensile state. The example given in the repository of this test case corresponds to a simplification with only 5 grains in order to enable an easy comparison with a reference solution given by the Cast3M code. 

![Polycristal Case](img/polycrystal/polycrystal.png)

More details in [polycrystal/README.md](./polycrystal/README.md)

## Cermet Case

This example illustrates the simulation of a UO₂ polycrystalline microstructure with metallic interfaces at the grain boundaries (CERMET). The objective is to study the mechanical viscoplastic response of the material, including the strain localization in the metallic interface, under a macroscopic uniaxial loading. In addition to the mechanical analysis, this example implements a fixed-point algorithm to converge toward a macroscopic uniaxial tensile state. The example given in the repository of this test case corresponds to a simplification with only 5 grains in order to enable an easy comparison with a reference solution given by the Cast3M code.

![Cermet Case](img/cermet/cermet.png)

More details in [cermet/README.md](./cermet/README.md)
