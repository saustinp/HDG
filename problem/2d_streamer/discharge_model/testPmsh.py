import gmsh
import numpy as np
# in this example I demonstrate how to 
# embed points in a line

gmsh.initialize()
gmsh.model.add('embedded_points_in_interval')

# I do this in 1D, it works in 2D, 3D in the same way

# some points you want in the mesh
N=10
x_internal = np.linspace(0,1,N) 

# let's make local mesh size at each point so large that the 
# final mesh only contains the line, start point, end point
# and embedded points
lc = .1

internal_tags = []
for x in x_internal:
    # adds each point to the model
    tag = gmsh.model.geo.add_point(x, 0, 0, lc)
    # collect their tags, we need them for embedding them 
    # later in the line
    internal_tags.append(tag)

# add starting and endpoint of line
p_start = gmsh.model.geo.add_point(-1, 0, 0, lc)
p_end = gmsh.model.geo.add_point(2, 0, 0, lc)

# make line from point tags p_start and p_end
line_tag = gmsh.model.geo.add_line(p_start,p_end) 
# call synchronization to make CAD kernel aware of everything
gmsh.model.geo.synchronize()

# embed all internal points
# dim = 0 since we embed points
# inDim = 1 since target is a line
#
# from gmsh.py, calling signature:
# def embed(dim, tags, inDim, inTag): 
gmsh.model.mesh.embed(0, internal_tags,1, line_tag) 

# generate 1d mesh
# since lc is so large, it will only contain 
# start point end point, line and the embedded points
gmsh.model.mesh.generate(1)

gmsh.write("embedded_points_in_interval.msh")
gmsh.fltk.run()
gmsh.finalize()