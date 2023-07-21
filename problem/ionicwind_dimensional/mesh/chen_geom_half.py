import gmsh
import numpy as np
from gmsh import model as gm

#############################################
case_name = 'chen_geom'
#############################################

gmsh.initialize()
gm.add(case_name)

############## CREATE PLANE GEOMETRY ##############

gm.occ.importShapes('chen_geom_accurate_tip.STEP')

# Exraction of the bbox parameters
xmin, ymin, zmin, xmax, ymax, zmax = gm.occ.getBoundingBox(2, 1)
size_scale_factor = .1/(xmax-xmin)

# Scale and rotate mesh -> This is to align the imported CAD geometry with the desired CSYS and scaling for the problem.
dimtags = gm.occ.getEntities()
gmsh.model.occ.dilate(dimtags, 0,0,0, size_scale_factor, size_scale_factor, size_scale_factor)

# gmsh.model.occ.rotate(dimtags, 0,0,0, 1, 0, 0, np.pi/2)
xmin, ymin, zmin, xmax, ymax, zmax = gm.occ.getBoundingBox(2, 1)    # Need to over-write the BBox values so that we can generate an accurate length scale later

gm.occ.synchronize()
# print(gm.getEntities())
# gmsh.fltk.run()
# exit()

xlength = xmax-xmin
ylength = ymax-ymin
zlength = zmax-zmin

ref_len = max(xlength, ylength, zlength)

# Needle tip chen_geom_accurate_tip
gm.mesh.field.add("Distance", 4)
gm.mesh.field.setNumbers(4, "CurvesList", [2])
gm.mesh.field.setNumbers(4, "PointsList", [3])
gm.mesh.field.setNumber(4, "Sampling", 96)      # Tweaked this value just to hit 90k elements

# Distance field - region 1
gmsh.model.mesh.field.add("MathEval", 14)
gmsh.model.mesh.field.setString(14, "F", "0.08*F4 + 3e-6")

# Needle body
gm.mesh.field.add("Distance", 3)
gm.mesh.field.setNumbers(3, "CurvesList", [1,11,12,13])
gmsh.model.mesh.field.add("MathEval", 12)
gmsh.model.mesh.field.setString(12, "F", "5*(F3/{})^2 + .0008".format(ref_len))

gmsh.model.mesh.field.add("MathEval", 17)
gmsh.model.mesh.field.setString(17, "F", ".003")        # To allow us to visualize the distance field after the fact

# Use the minimum of all the fields as the background mesh field:
gm.mesh.field.add("Min", 20)
gm.mesh.field.setNumbers(20, "FieldsList", [12, 14])

gmsh.model.addPhysicalGroup(1, [2,3,4,5,6], name='needle')
gmsh.model.addPhysicalGroup(1, [7], name='symmetry_axis')
gmsh.model.addPhysicalGroup(1, [8], name='outflow')
gmsh.model.addPhysicalGroup(1, [9,10,11,12], name='ground_and_cylinder')
gmsh.model.addPhysicalGroup(1, [1, 13], name='farfield')
gmsh.model.addPhysicalGroup(2, [1], name='omega')

gm.mesh.field.setAsBackgroundMesh(20)
gmsh.option.setNumber("Mesh.MeshSizeMax", .003)
gm.mesh.setSmoothing(2, 1, 10)

gm.mesh.generate(2)
gmsh.write('chen_test.msh3')
gmsh.fltk.run()
gmsh.finalize()
exit()
