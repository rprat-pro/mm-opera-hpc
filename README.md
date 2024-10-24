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


## Test Case : 1 


TODO
