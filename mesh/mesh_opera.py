import sac_de_billes as sdb
import merope

dimensions_side = 10  # µm
L = [dimensions_side, dimensions_side, dimensions_side]

radius = 0.5  # µm
pore_fraction = 0.10  # /
distMin = radius*1.2
randomSeed = 0


theSpheres = sdb.throwSpheres_3D(sdb.TypeAlgo.RSA, sdb.NameShape.Tore, L, randomSeed, [
    [radius, pore_fraction]], [2], distMin)


tab_phase_not_to_mesh = []
for j in range(0, len(theSpheres)):
    sphere = theSpheres[j]
    sphere.phase = 2+j
    tab_phase_not_to_mesh.append(sphere.phase)

with open("bubbles.txt", 'w') as bbl:
    for sphere in theSpheres:
        bbl.write(str(sphere.phase) + " ")
        bbl.write(str(sphere.center[0]) + " " +
                  str(sphere.center[1]) + " " +
                  str(sphere.center[2]) + " ")
        bbl.write(str(sphere.radius) + "\n")

sphInc = merope.SphereInclusions_3D()
sphInc.setLength(L)
sphInc.setSpheres(theSpheres)

mi = merope.MultiInclusions_3D()
mi.setInclusions(sphInc)
mi.setMatrixPhase(1)

meshGenerator = merope.mesh.MeshGenerator()
meshGenerator.setMeshOrder(2)
meshGenerator.setMeshSize(dimensions_side*1e-2)
meshGenerator.setMultiInclusions(mi)
meshGenerator.do_not_mesh(tab_phase_not_to_mesh)
meshGenerator.set_nameOutput(["spheres.vtk"])
meshGenerator.write("mesh_spheres.geo")
