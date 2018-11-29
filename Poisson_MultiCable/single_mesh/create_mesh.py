import pygmsh
import meshio
from dolfin import timed
from IPython import embed
outer_radius = 1.2
rubber_radius = 0.255
metal_radius = 0.2
fill_marker = 15
rubber_marker = 16
metal_marker = 17
fill = 11
metaliso = 12
isofill = 13

@timed("USER_TIMING: Creating mesh")
def create_multicable(cable_pos, res=0.1):
    geo = pygmsh.built_in.geometry.Geometry()
    num_cables = int(len(cable_pos)/2)
    metal_circles = []
    rubber_circles = []
    for i in range(num_cables):
        metal_circles.append(geo.add_circle((cable_pos[2*i],cable_pos[2*i+1],0),
                                            metal_radius,lcar=res/4))
        rubber_circles.append(geo.add_circle((cable_pos[2*i],cable_pos[2*i+1],0)
                                             , rubber_radius,
                                             lcar=res/2,
                                             holes=[metal_circles[i]]))
    fill_circle = geo.add_circle((0,0,0), outer_radius, lcar=res
                                 ,holes=metal_circles+rubber_circles)

    geo.add_physical_surface([metal.plane_surface for metal in metal_circles],
                             label=metal_marker)
    geo.add_physical_surface([rubber.plane_surface for rubber in rubber_circles]
                             , label=rubber_marker)
    geo.add_physical_surface(fill_circle.plane_surface, label=fill_marker)

    metal_lines = [metal.line_loop.lines for metal in metal_circles]
    metal_flat = [val for sublist in metal_lines for val in sublist]
    geo.add_physical_line(metal_flat, label=metaliso)

    rubber_lines = [rubber.line_loop.lines for rubber in rubber_circles]
    rubber_flat = [val for sublist in rubber_lines for val in sublist]
    geo.add_physical_line(rubber_flat, label=isofill)
    geo.add_physical_line(fill_circle.line_loop.lines, label=fill)
    
    (points, cells, point_data,
     cell_data, field_data) = pygmsh.generate_mesh(geo, prune_z_0=True,
                                                   geo_filename="test.geo")
    meshio.write("multicable.xdmf",
                 meshio.Mesh(points=points, cells={"triangle": cells["triangle"]}))
    meshio.write("mf.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
    cell_data={"line": {"name_to_read":
                        cell_data["line"]["gmsh:physical"]}}))
    meshio.write("cf.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]},
        cell_data={"triangle": {"name_to_read":
                                cell_data["triangle"]["gmsh:physical"]}}))
    
