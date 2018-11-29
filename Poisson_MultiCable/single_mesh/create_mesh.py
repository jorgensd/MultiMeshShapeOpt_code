import pygmsh
import meshio
from dolfin import timed
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
def create_multicable(x0,y0, res=0.1):
    geo = pygmsh.built_in.geometry.Geometry()
    metal_circle = geo.add_circle((x0,y0,0), metal_radius, lcar=res/4)
    rubber_circle = geo.add_circle((x0,y0,0), rubber_radius, lcar=res/2,
                                   holes=[metal_circle])
    fill_circle = geo.add_circle((0,0,0), outer_radius, lcar=res
                                 ,holes=[metal_circle,rubber_circle])
    geo.add_physical_surface(metal_circle.plane_surface, label=metal_marker)
    geo.add_physical_surface(rubber_circle.plane_surface, label=rubber_marker)
    geo.add_physical_surface(fill_circle.plane_surface, label=fill_marker)
    geo.add_physical_line(metal_circle.line_loop.lines, label=metaliso)
    geo.add_physical_line(rubber_circle.line_loop.lines, label=isofill)
    geo.add_physical_line(fill_circle.line_loop.lines, label=fill)
    
    (points, cells, point_data,
     cell_data, field_data) = pygmsh.generate_mesh(geo, prune_z_0=True)
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
    
