from pygmsh import generate_mesh
from pygmsh.built_in.geometry import Geometry
import meshio
import os;
from IPython import embed

outer_marker = 1
inner_marker = 2
inflow = 1
outflow = 2
walls = 3
L = 1
H = 1
c_x,c_y  = L/2, H/2
r_x = 0.126157
width_scale = 6 # Number of cells width of front  mesh w.r.t to background mesh size

def background_mesh(res=0.025):
    """
    Create rectangular background mesh for the Poisson problem
    """
    geometry = Geometry()
    square = geometry.add_rectangle(0,L,0,H,0, res)
    geometry.add_physical_surface([square.surface],label=12)
    geometry.add_physical_line([square.line_loop.lines[3]], label=inflow)
    geometry.add_physical_line([square.line_loop.lines[1]], label=outflow)
    geometry.add_physical_line([square.line_loop.lines[0],
                                square.line_loop.lines[2]], label=walls)


    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geometry, prune_z_0=True,dim=2,
                                            geo_filename="meshes/tmp.geo")
    meshio.write("multimesh_0.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write("mf_0.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))
    import os
    os.system("mv multimesh_0.* meshes/")
    os.system("mv mf_0.* meshes/")

def front_mesh_symmetric(res=0.025):
    """
    Creates a donut mesh with symmetry line through y = c_y
    """
    geometry = Geometry()
    mesh_r = width_scale*res
    c = geometry.add_point((c_x,c_y,0))
    
    # Elliptic obstacle
    p1 = geometry.add_point((c_x-r_x, c_y,0))
    p2 = geometry.add_point((c_x+r_x, c_y,0))
    # embed()
    arc_1 = geometry.add_circle_arc(p1, c, p2)
    arc_2 = geometry.add_circle_arc(p2, c, p1)

    # Surrounding mesh
    p3 = geometry.add_point((c_x-r_x-mesh_r, c_y,0))
    p4 = geometry.add_point((c_x+r_x+mesh_r, c_y,0))
    arc_5 = geometry.add_circle_arc(p3, c, p4)
    arc_6 = geometry.add_circle_arc(p4, c, p3)
    line_front = geometry.add_line(p3, p1)
    line_back = geometry.add_line(p2, p4)    
    surface_top = geometry.add_line_loop([line_front, arc_1, line_back, -arc_5])
    surface_bottom = geometry.add_line_loop([line_front, -arc_2, line_back,
                                             arc_6])

    obstacle_loop = geometry.add_line_loop([arc_1, arc_2])
    outer_loop = geometry.add_line_loop([arc_5,arc_6])
    
    top = geometry.add_plane_surface(surface_top)
    bottom = geometry.add_plane_surface(surface_bottom)
    geometry.add_physical_surface([top, bottom],label=12)
    geometry.add_physical_line(obstacle_loop.lines, label=inner_marker)
    geometry.add_physical_line(outer_loop.lines, label=outer_marker)

    # Create refined mesh around geometry
    field = geometry.add_boundary_layer(edges_list=obstacle_loop.lines,
                                        hfar=res, hwall_n=res/3,
                                        thickness=0.1*mesh_r)
    geometry.add_background_field([field])

    # Generate mesh
    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geometry, prune_z_0=True,
                                            geo_filename="meshes/test.geo")

    # Save mesh and mesh-function to file
    meshio.write("multimesh_1.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write("mf_1.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))
    import os
    os.system("mv multimesh_1.* meshes/")
    os.system("mv mf_1.* meshes/")
def front_mesh_wedge(res=0.025):
    """
    Creates a wedged mesh with symmetry line through y = c_y
    """
    geometry = Geometry()
    mesh_r = width_scale*res
    c = geometry.add_point((c_x,c_y,0))
    
    # Elliptic obstacle
    y_fac = 1
    p1 = geometry.add_point((c_x-r_x, c_y,0),lcar=res/8)
    pt = geometry.add_point((c_x,c_y+y_fac*r_x,0),lcar=res/4)
    pb = geometry.add_point((c_x, c_y-y_fac*r_x,0),lcar=res/4)
    p2 = geometry.add_point((c_x+r_x, c_y,0),lcar=res/8)
    # embed()
    arc_1 = geometry.add_bspline([p1, pt, p2])
    arc_2 = geometry.add_bspline([p2, pb, p1])

    # Surrounding mesh
    p3 = geometry.add_point((c_x-r_x-mesh_r, c_y,0),lcar=res)
    p4 = geometry.add_point((c_x+r_x+mesh_r, c_y,0),lcar=res)
    pt_ = geometry.add_point((c_x, c_y+y_fac*r_x+mesh_r,0),lcar=res)
    pb_ = geometry.add_point((c_x, c_y-y_fac*r_x-mesh_r,0),lcar=res)

    arc_5 = geometry.add_bspline([p3, pt_, p4])
    arc_6 = geometry.add_bspline([p4, pb_, p3])
    obstacle_loop = geometry.add_line_loop([arc_1, arc_2])
    outer_loop = geometry.add_line_loop([arc_5,arc_6])
    donut = geometry.add_plane_surface(obstacle_loop, holes = [outer_loop])
    geometry.add_physical_surface([donut],label=12)
    geometry.add_physical_line(obstacle_loop.lines, label=inner_marker)
    geometry.add_physical_line(outer_loop.lines, label=outer_marker)

    # Generate mesh
    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geometry, prune_z_0=True,
                                            geo_filename="meshes/test.geo")
    
    # Save mesh and mesh-function to file
    meshio.write("multimesh_1.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write("mf_1.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))
    import os
    os.system("mv multimesh_1.* meshes/")
    os.system("mv mf_1.* meshes/")

    
def front_mesh_unsym(res=0.025):
    """
    Creates unsymmetric donut mesh
    """
    geometry = Geometry()
    mesh_r = width_scale*res
    c = geometry.add_point((c_x,c_y,0))
    
    # Elliptic obstacle
    p1 = geometry.add_point((c_x-r_x, c_y,0))
    p2 = geometry.add_point((c_x+r_x, c_y,0))
    # embed()
    arc_1 = geometry.add_circle_arc(p1, c, p2)
    arc_2 = geometry.add_circle_arc(p2, c, p1)

    # Surrounding mesh
    p3 = geometry.add_point((c_x-r_x-mesh_r, c_y,0))
    p4 = geometry.add_point((c_x+r_x+mesh_r, c_y,0))
    arc_5 = geometry.add_circle_arc(p3, c, p4)
    arc_6 = geometry.add_circle_arc(p4, c, p3)
    obstacle_loop = geometry.add_line_loop([arc_1, arc_2])
    outer_loop = geometry.add_line_loop([arc_5,arc_6])
    
    donut = geometry.add_plane_surface(obstacle_loop, holes=[outer_loop])
   
    geometry.add_physical_surface([donut],label=12)
    geometry.add_physical_line(obstacle_loop.lines, label=inner_marker)
    geometry.add_physical_line(outer_loop.lines, label=outer_marker)

    # Create refined mesh around geometry
    field = geometry.add_boundary_layer(edges_list=obstacle_loop.lines,
                                        hfar=res, hwall_n=res/3,
                                        thickness=2*res)
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




def single_mesh(res=0.025):
    """ 
    Creates a single mesh containing a circular obstacle
    """

    geometry = Geometry()
    c = geometry.add_point((c_x,c_x,0))    

    # Elliptic obstacle
    p1 = geometry.add_point((c_x-r_x, c_x,0))
    p2 = geometry.add_point((c_x, c_x+r_x,0))
    p3 = geometry.add_point((c_x+r_x, c_x,0))
    p4 = geometry.add_point((c_x, c_x-r_x,0))
    arc_1 = geometry.add_ellipse_arc(p1, c, p2, p2)
    arc_2 = geometry.add_ellipse_arc(p2, c, p3, p3)
    arc_3 = geometry.add_ellipse_arc(p3, c, p4, p4)
    arc_4 = geometry.add_ellipse_arc(p4, c, p1, p1)
    obstacle_loop = geometry.add_line_loop([arc_1, arc_2, arc_3, arc_4])

    rectangle = geometry.add_rectangle(0,L,0,H,0, res, holes=[obstacle_loop])
    flow_list = [rectangle.line_loop.lines[0], rectangle.line_loop.lines[2],
                 rectangle.line_loop.lines[3]]
    wall_list = obstacle_loop.lines
    geometry.add_physical_surface(rectangle.surface,label=12)
    geometry.add_physical_line(flow_list, label=inflow)
    geometry.add_physical_line([rectangle.line_loop.lines[1]], label=outflow)
    geometry.add_physical_line(wall_list, label=walls)
    field = geometry.add_boundary_layer(edges_list=obstacle_loop.lines,
                                        hfar=res, hwall_n=res/2, thickness=2*res)
    geometry.add_background_field([field])

    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geometry, prune_z_0=True,
                                            geo_filename="meshes/singlemesh.geo")
    
    meshio.write("meshes/singlemesh.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write("meshes/mf.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))

        
if __name__=="__main__":
    import sys
    try:
        res = float(sys.argv[1])
    except IndexError:
        res = 0.01
    background_mesh(res)
    # front_mesh_symmetric(res)
    # front_mesh_unsym(res)
    front_mesh_wedge(res)
    # single_mesh(res)
