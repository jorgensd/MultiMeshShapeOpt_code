from pygmsh.built_in.geometry import Geometry
from pygmsh import generate_mesh
import meshio
outer_marker = 2
inner_marker = 1
ball_height = 1.25

import os
mesh_folder = "notebook_meshes/"
os.system("mkdir -p {0:s}".format(mesh_folder))

def single_mesh(res=0.025):
    geometry = Geometry()
    circle = geometry.add_circle((0.5,ball_height,0),0.15, lcar=res)
    square = geometry.add_rectangle(0,1,0,1.5,0, lcar=res)
    mesh = geometry.add_plane_surface(square.line_loop,
                                      holes=[circle.line_loop])
    geometry.add_physical_surface(mesh,label=12)
    geometry.add_physical_line(square.line_loop.lines, label=outer_marker)
    geometry.add_physical_line(circle.line_loop.lines, label=inner_marker)
    
    # # Generate mesh
    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geometry, prune_z_0=True,
                                            verbose=False)

    # Save mesh and mesh-function to file
    meshio.write(mesh_folder+"singlemesh.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write(mesh_folder+"mf_singlemesh.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))

def background_mesh(res=0.025):
    geometry = Geometry()
    square = geometry.add_rectangle(0,1,0,1.5,0, lcar=res)
    mesh = geometry.add_plane_surface(square.line_loop)
    geometry.add_physical_surface(mesh,label=12)
    geometry.add_physical_line(square.line_loop.lines, label=outer_marker)
    
    # # Generate mesh
    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geometry, prune_z_0=True,
                                            verbose=False)

    # Save mesh and mesh-function to file
    meshio.write(mesh_folder+"background_mesh.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write(mesh_folder+"mf_background.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))


def front_mesh(res=0.025):
    geometry = Geometry()
    circle = geometry.add_circle((0.5,ball_height,0),0.15, lcar=res)
    outer_circle = geometry.add_circle((0.5,ball_height,0),0.25, lcar=res)
    mesh = geometry.add_plane_surface(outer_circle.line_loop,
                                      holes=[circle.line_loop])
    geometry.add_physical_surface(mesh,label=12)
    geometry.add_physical_line(outer_circle.line_loop.lines, label=outer_marker)
    geometry.add_physical_line(circle.line_loop.lines, label=inner_marker)
    
    # # Generate mesh
    (points, cells, point_data,
     cell_data, field_data) = generate_mesh(geometry, prune_z_0=True,
                                            verbose=False)

    # Save mesh and mesh-function to file
    meshio.write(mesh_folder+"front_mesh.xdmf", meshio.Mesh(
        points=points, cells={"triangle": cells["triangle"]}))
    
    meshio.write(mesh_folder+"mf_front.xdmf", meshio.Mesh(
        points=points, cells={"line": cells["line"]},
        cell_data={"line": {"name_to_read":
                            cell_data["line"]["gmsh:physical"]}}))
