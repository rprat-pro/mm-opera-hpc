# Cermet simulation

## Mesh generation

Mesh example: cermet-opera-hpc-mini.py


Please load `MEROPE_DIR` before.


```
source ${MEROPE_DIR}/Env_Merope.sh
python3 mesh/cermet-opera-hpc-mini.py # generate cermet-mini.geo
gmsh -3 cermet-mini.geo # generate cermet-mini.msh
```


### Output Example

![Cermet Case](doc/cermet-mini.png)

