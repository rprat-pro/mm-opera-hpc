# mm-opera-hpc

This directory groups together the various test cases implemented as part of the OperaHPC project.

## Installation using Spack

```
git clone https://github.com/spack/spack.git
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
spack install mfem-mgis
```

Thirdly, load mfem-mgis

```
spack load mfem-mgis
export MFEMMGIS_DIR=`spack location -i mfem-mgis`/share/mfem-mgis/cmake/
```

Finally, build your examples:

```
cd mfem-mgis-examples
mkdir build && cd build
cmake ..
make -j 4
```


## Bubble Case 

### Short Description

TODO : The default example is the rupture of a stressed bubble.

### Parameters

```
Usage: ./test-bubble [options] ...
Options:
   -h, --help
	Print this help message and exit.
   -m <string>, --mesh <string>, current value: mesh/mesh_sphere.msh
	Mesh file to use.
   -l <string>, --library <string>, current value: src/libBehaviour.so
	Material library.
   -f <string>, --bubble-file <string>, current value: bubbles.txt
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


### Output With On Bubble

![Bubble Case](/img/bubble.png)
