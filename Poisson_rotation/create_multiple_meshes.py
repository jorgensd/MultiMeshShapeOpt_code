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

def background_mesh(res=0.025):
    """
    Create rectangular background mesh for the Poisson problem
    """
    geometry = Geometry()
    rectangle = geometry.add_rectangle(0,L,0,H,0, res)

    geometry.add_physical_surface(rectangle.surface,label=12)
    geometry.add_physical_line(rectangle.line_loop.lines, label=outer_marker)

    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geometry, prune_z_0=True)
    
    meshio.write("meshes/multimesh_0.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write("meshes/mf_0.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))

def front_mesh(res=0.025):
    geometry = Geometry()
    c = geometry.add_point((c_x,c_y,0))
    mesh_r = 3*res # Width of mesh
    
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

    # Surrounding mesh
    p1 = geometry.add_point((c_x-r_x-mesh_r, c_y,0))
    p2 = geometry.add_point((c_x, c_y+r_y+mesh_r,0))
    p3 = geometry.add_point((c_x+r_x+mesh_r, c_y,0))
    p4 = geometry.add_point((c_x, c_y-r_y-mesh_r,0))
    arc_5 = geometry.add_ellipse_arc(p1, c, p2, p2)
    arc_6 = geometry.add_ellipse_arc(p2, c, p3, p3)
    arc_7 = geometry.add_ellipse_arc(p3, c, p4, p4)
    arc_8 = geometry.add_ellipse_arc(p4, c, p1, p1)
    loop = geometry.add_line_loop([arc_5, arc_6, arc_7, arc_8])
    mesh = geometry.add_plane_surface(loop, holes=[obstacle_loop])
    geometry.add_physical_surface(mesh,label=12)
    geometry.add_physical_line(loop.lines, label=outer_marker)
    geometry.add_physical_line(obstacle_loop.lines, label=inner_marker)

    # Create refined mesh around geometry
    field = geometry.add_boundary_layer(edges_list=obstacle_loop.lines, hfar=res, hwall_n=res/2, thickness=2*res)
    geometry.add_background_field([field])

    # Generate mesh
    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geometry, prune_z_0=True,
                                            geo_filename="meshes/test.geo")

    # Save mesh and mesh-function to file
    meshio.write("meshes/multimesh_1.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write("meshes/mf_1.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))

        
if __name__=="__main__":
    import sys
    try:
        res = float(sys.argv[1])
    except IndexError:
        res = 0.025
    background_mesh(res)
    front_mesh(res)
