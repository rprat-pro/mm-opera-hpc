# Cermet simulation

## Mesh generation

Mesh example: `cermet-opera-hpc-mini.py`


Please load `MEROPE_DIR` before.

```
source ${MEROPE_DIR}/Env_Merope.sh
python3 mesh/cermet-opera-hpc-mini.py # generate cermet-mini.geo
gmsh -3 cermet-mini.geo # generate cermet-mini.msh
```

### Options

In your python script: 

```
L = [1, 1, 1]
nbSpheres = 20 
distMin = 0.3
randomSeed = 0
layer=0.02
MeshOrder = 1
MeshSize = 0.05
```

Example to generate 250 grains, 12913361 nodes, and 86213779 elements.


```
L = [5, 5, 5]
nbSpheres = 250
distMin = 0.1
randomSeed = 0
layer=0.04
MeshOrder = 1
MeshSize = 0.02
```

### Output Example: cermet-mini.msh

![Cermet Case](doc/cermet-mini.png)

![Cermet Case png](doc/cermet-mini-png.png)

### Output Example: cermet-big.msh

![Cermet Case png](doc/cermet-big-png.png)
