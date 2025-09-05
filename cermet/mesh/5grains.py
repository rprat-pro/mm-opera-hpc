# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 31/01/2025
#
# Copyright : see License.txt
#
# Build a .geo file to be meshed by gmsh

import sac_de_billes
import merope

L = [1, 1, 1]
nbSpheres = 5
distMin = 0.4
randomSeed = 0
layer=0.02
MeshOrder = 1
MeshSize = 0.05

theSpheres = sac_de_billes.fillMaxRSA_3D(sac_de_billes.Tore, L, nbSpheres, randomSeed, distMin)

i = 2
for sphere in theSpheres:
    sphere.phase = i
    i = i+1


sphInc = merope.LaguerreTess_3D(L,theSpheres)
mi = merope.MultiInclusions_3D()
mi.setInclusions(sphInc)
mi.addLayer(mi.getAllIdentifiers(), 1, layer)

meshGenerator = merope.mesh.MeshGenerator()
meshGenerator.setMeshOrder(MeshOrder)
meshGenerator.setMeshSize(MeshSize)
meshGenerator.setMultiInclusions(mi)
meshGenerator.set_nameOutput(["5grains.med"])
meshGenerator.setAdimMergeDistance(1e-7)
meshGenerator.write("5grains.geo", merope.mesh.MeshMethod.OpenCascade)
