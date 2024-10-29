# Bubble case - some useful tips

## Modify the geometry for the single bubble case and mesh it

The geometry for the testcase is contained in the file `.geo` stored in the `mesh` folder. One can modify it and use it as an input for `gmsh` to generate the computational mesh for the case by 

```bash
gmsh -3 single_sphere.geo
```
NB: if the bubble center, radius, or the surface label are modified, the corresponding data stored in `single_bubble.txt` must also be changed.

## Verification against the analytical solution

The script available in `verification/bubble` can be used to compare the analytical solution to the MMM one, by just placing the script in the post-processing folder and executing it 

```bash
python3 mmm_vs_analytical.py 
```