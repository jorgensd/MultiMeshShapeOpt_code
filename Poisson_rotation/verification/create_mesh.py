from pygmsh import generate_mesh
from pygmsh.built_in.geometry import Geometry
import meshio
import os;
from IPython import embed
os.system("mkdir -p meshes")

outer_marker = 1
inner_marker = 2
L = 2.5
H = 1.75
c_x,c_y  = L/2+0.1, H/2
r_x, r_y = 0.3, 0.1

def single_mesh(res=0.025):
    geometry = Geometry()
    c = geometry.add_point((c_x,c_y,0))    

    # Elliptic obstacle
    p1 = geometry.add_point((c_x-r_x, c_y,0))
    p2 = geometry.add_point((c_x, c_y+r_y,0))
    p3 = geometry.add_point((c_x+r_x, c_y,0))
    p4 = geometry.add_point((c_x, c_y-r_y,0))
    arc_1 = geometry.add_ellipse_arc(p1, c, p2, p2)
    arc_2 = geometry.add_ellipse_arc(p2, c, p3, p3)
    arc_3 = geometry.add_ellipse_arc(p3, c, p4, p4)
    arc_4 = geometry.add_ellipse_arc(p4, c, p1, p1)
    obstacle_loop = geometry.add_line_loop([arc_1, arc_2, arc_3, arc_4])
    obstacle = geometry.add_surface(obstacle_loop)
    
    rectangle = geometry.add_rectangle(0,L,0,H,0, res, holes=[obstacle_loop])
    geometry.add_physical_surface(obstacle, label=13)
    geometry.add_physical_surface(rectangle.surface,label=12)
    geometry.add_physical_line(obstacle_loop.lines, label=inner_marker)
    geometry.add_physical_line(rectangle.line_loop.lines, label=outer_marker)
    field = geometry.add_boundary_layer(edges_list=obstacle_loop.lines,
                                        hfar=res, hwall_n=res/2, thickness=2*res)
    geometry.add_background_field([field])

    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geometry, prune_z_0=True,
                                            geo_filename="meshes/outfile.geo")
    
    meshio.write("meshes/singlemesh.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write("meshes/mf.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))
    meshio.write("meshes/cf.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]},
        cell_data={"triangle": {"name_to_read":
                            cell_data["triangle"]["gmsh:physical"]}}))

if __name__=="__main__":
    import sys
    try:
        res = float(sys.argv[1])
    except IndexError:
        res = 0.1
    single_mesh(res)
