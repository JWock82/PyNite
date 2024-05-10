from math import pi
from PyNite.FEModel3D import FEModel3D
from PyNite.Mesh import FrustrumMesh, CylinderMesh

t = 0.25/12
E = 29000*1000*12**2
G = 11200*1000*12**2
nu = 0.3
kx_mod=1
ky_mod=1
mesh_size = 0.5
r_shell = 7.5
r_hopper = 1.5
h_shell = 20
h_hopper = 10
center = [0, 0, 0]
axis = 'Y'
start_node = 'N1'
start_element = 'Q1'

# Create a finite element model
model = FEModel3D()

# Add steel as a material
model.add_material('Steel', E, G, nu, 490)

# Add the conical hopper mesh to the model
model.add_frustrum_mesh('MSH1', mesh_size, r_shell, r_hopper, h_hopper, t, 'Steel', kx_mod, ky_mod, center, axis)

# Generate the mesh. The mesh would automatically be generated during analysis, but we want to
# manipulate it now so we'll generate it in advance.
model.Meshes['MSH1'].generate()

# Determine how many elements make up the circumference at the top of the hopper. This will
# determine the number of elements making up the circumference of the cylindrical shell.
n_hopper = model.Meshes['MSH1'].num_quads_outer

# Find an unused node and element name to start the cylinder off with
first_node = 'N' + str(len(model.Nodes.values()) + 1)
first_element = 'Q' + str(len(model.Quads.values()) + 1)

# Add the cylindrical shell mesh
model.add_cylinder_mesh('MSH2', mesh_size, r_shell, h_shell, t, 'Steel', kx_mod, ky_mod, center,
                        axis, num_elements=n_hopper, element_type='Quad')

# The two meshes have overlapping duplicate nodes that need to be merged
model.merge_duplicate_nodes()

# This next line can be used in place of the prior line to profile the merge code and find
# innefficiencies
# cProfile.run('model.merge_duplicate_nodes()', sort='cumtime')  

# Add hydrostatic pressure to each element in the model
for element in model.Quads.values():
    Yavg = (element.i_node.Y + element.j_node.Y + element.m_node.Y + element.n_node.Y)/4
    model.add_quad_surface_pressure(element.name, 62.4*(h_shell - Yavg), case='Hydrostatic')

# Add supports at the springline
for node in model.Nodes.values():
    if round(node.Y, 10) == 0:
        model.def_support(node.name, True, True, True, False, False, False)

# Add a load combination
model.add_load_combo('1.4F', {'Hydrostatic': 1.4})

# Analyze the model
# import cProfile
# cProfile.run('model.analyze()', sort='cumtime')
model.analyze(log=True)

# Render the model
from PyNite.Rendering import Renderer
rndr = Renderer(model)
rndr.annotation_size = 0.1
rndr.render_loads = False
rndr.color_map = 'SX'
rndr.combo_name = '1.4F'
rndr.show_labels = False
rndr.render_model()