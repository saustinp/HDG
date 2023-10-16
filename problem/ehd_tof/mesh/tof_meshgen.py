# ------------------------------------------------------------------------------
#
#  Taken from Gmsh tutorial 1 with the point coordinates changed
#
# ------------------------------------------------------------------------------

# The Python API is entirely defined in the `gmsh.py' module (which contains the
# full documentation of all the functions in the API):
import gmsh
import sys

# Before using any functions in the Python API, Gmsh must be initialized:
gmsh.initialize()

# Next we add a new model named "t1" (if gmsh.model.add() is not called a new
# unnamed model will be created on the fly, if necessary):
gmsh.model.add("t1")

lc = 7.5e-6
gmsh.model.geo.addPoint(0, 0, 0, lc, 1)

gmsh.model.geo.addPoint(0, 1e-3, 0, lc, 2)
gmsh.model.geo.addPoint(.5e-3, 1e-3, 0, lc, 3)

p4 = gmsh.model.geo.addPoint(.5e-3, 0, 0, lc)

gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(3, 2, 2)
gmsh.model.geo.addLine(3, p4, 3)
gmsh.model.geo.addLine(4, 1, p4)

gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1)

gmsh.model.geo.addPlaneSurface([1], 1)

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [1, 2, 4], 5)
gmsh.model.addPhysicalGroup(2, [1], name = "My surface")

gmsh.model.mesh.generate(2)

gmsh.write("tof_mesh20k.msh3")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
