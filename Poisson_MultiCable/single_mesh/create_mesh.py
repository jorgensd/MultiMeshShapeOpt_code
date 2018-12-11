import pygmsh
import meshio
from dolfin import timed, Timer
from IPython import embed
outer_radius = 1.2
rubber_radius = 0.255
metal_radius = 0.2
fill_marker = 115
rubber_marker = 216
metal_marker = 317
ext = 11
metaliso = 100
isofill = 200

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
        geo.add_physical_surface(metal_circles[i].plane_surface,
                                 label=int(metal_marker+i))
        geo.add_physical_surface(rubber_circles[i].plane_surface,
                                 label=int(rubber_marker+i))
        geo.add_physical_line(metal_circles[i].line_loop.lines,
                              label=int(metaliso+i))
        geo.add_physical_line(rubber_circles[i].line_loop.lines,
                              label=int(isofill+i))

    fill_circle = geo.add_circle((0,0,0), outer_radius, lcar=res
                                 ,holes=metal_circles+rubber_circles)

    geo.add_physical_surface(fill_circle.plane_surface, label=fill_marker)

    geo.add_physical_line(fill_circle.line_loop.lines, label=ext)

    with Timer("USER_TIMING: Generate mesh") as t:
        (points, cells, point_data,
         cell_data, field_data) = pygmsh.generate_mesh(geo, prune_z_0=True,
                                                       verbose=True)#,
                                                       # geo_filename="mesh.geo")

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
if __name__ == "__main__":
    create_multicable([0.1,0.1,-0.6,-0.7])
