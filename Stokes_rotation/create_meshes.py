from pygmsh import generate_mesh
from pygmsh.built_in.geometry import Geometry
import meshio
import os;
from IPython import embed
os.system("mkdir -p meshes")

inlet0_marker = 1
inlet1_marker = 2
outlet_marker = 3
obstacle_marker = 4
wall_marker = 5
inlet0 = [(0,0.15,0), (0,0.25,0)] # Coordinates of inlet
inlet1 = [(0,0.73,0), (0,0.83,0)]
outlet = [(0.75,1,0), (0.95,1,0)]

L = 1
H = 1
c_x,c_y  = L/2, H/2
r_x = 0.08
r_y = 0.04
width_scale = 3 # Number of cells width of front  mesh w.r.t to background mesh size

def background_mesh(res=0.025):
    """
    Create rectangular background mesh for the Poisson problem
    """
    geo = Geometry()
    # Create points for square with two inlets and one outlet
    points = []
    points.append(geo.add_point((0,0,0), res))
    points.append(geo.add_point((1,0,0), res))
    points.append(geo.add_point((1,1,0), res))
    points.append(geo.add_point(outlet[1], res))
    points.append(geo.add_point(outlet[0], res))
    points.append(geo.add_point((0,1,0), res))
    points.append(geo.add_point(inlet1[1], res))
    points.append(geo.add_point(inlet1[0], res))
    points.append(geo.add_point(inlet0[1], res))
    points.append(geo.add_point(inlet0[0], res))
    # Connect points as lines
    lines = []
    for i in range(len(points)-1):
        lines.append(geo.add_line(points[i],points[i+1]))
    lines.append(geo.add_line(points[-1], points[0]))
    # Create geometry
    line_loop = geo.add_line_loop(lines)
    square = geo.add_plane_surface(line_loop)
    
    # Create cell and facet function
    wall_lines = lines[0:3]+lines[4:6]+[lines[7]]+[lines[9]]
    inlet0_lines = [lines[8]]
    inlet1_lines = [lines[6]]
    outlet_lines = [lines[3]]
    geo.add_physical_surface([square],label=12)
    geo.add_physical_line(inlet0_lines, label=inlet0_marker)
    geo.add_physical_line(inlet1_lines, label=inlet1_marker)
    geo.add_physical_line(outlet_lines, label=outlet_marker)
    geo.add_physical_line(wall_lines, label=wall_marker)

    # Create mesh
    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geo, prune_z_0=True)#,
    #geo_filename="meshes/tmp.geo")

    meshio.write("meshes/multimesh_0.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write("meshes/mf_0.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))

def front_mesh(res=0.025):
    geo = Geometry()
    mesh_r = width_scale*res
    c = geo.add_point((c_x,c_y,0))
    
    p1 = geo.add_point((c_x-r_x, c_y,0))
    p2 = geo.add_point((c_x, c_y+r_y,0))
    p3 = geo.add_point((c_x+r_x, c_y,0))
    p4 = geo.add_point((c_x, c_y-r_y,0))
    arc_1 = geo.add_ellipse_arc(p1, c, p2, p2)
    arc_2 = geo.add_ellipse_arc(p2, c, p3, p3)
    arc_3 = geo.add_ellipse_arc(p3, c, p4, p4)
    arc_4 = geo.add_ellipse_arc(p4, c, p1, p1)
    obstacle_loop = geo.add_line_loop([arc_1, arc_2, arc_3, arc_4])

    p1 = geo.add_point((c_x-r_x-mesh_r, c_y,0))
    p2 = geo.add_point((c_x, c_y+r_y+mesh_r,0))
    p3 = geo.add_point((c_x+r_x+mesh_r, c_y,0))
    p4 = geo.add_point((c_x, c_y-r_y-mesh_r,0))
    arc_5 = geo.add_ellipse_arc(p1, c, p2, p2)
    arc_6 = geo.add_ellipse_arc(p2, c, p3, p3)
    arc_7 = geo.add_ellipse_arc(p3, c, p4, p4)
    arc_8 = geo.add_ellipse_arc(p4, c, p1, p1)
    loop = geo.add_line_loop([arc_5, arc_6, arc_7, arc_8])
    mesh = geo.add_plane_surface(loop, holes=[obstacle_loop])

    geo.add_physical_surface(mesh,label=12)
    geo.add_physical_line(obstacle_loop.lines, label=obstacle_marker)

    # Create refined mesh around geo
    field = geo.add_boundary_layer(edges_list=obstacle_loop.lines,
                                        hfar=res, hwall_n=res/3,
                                        thickness=2*res)
    geo.add_background_field([field])

    # Generate mesh
    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geo, prune_z_0=True)
    #,
    #                                       geo_filename="meshes/test.geo")

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
        res = 0.01
    background_mesh(res)
    front_mesh(res)
