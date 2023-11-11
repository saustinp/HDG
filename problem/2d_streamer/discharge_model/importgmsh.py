import gmsh
import numpy as np
gmsh.initialize()
gmsh.merge("streamer_16k-2.msh3")

p = gmsh.model.mesh.getNodes()
t = gmsh.model.mesh.getElements()

print(len(p))
# print(t)