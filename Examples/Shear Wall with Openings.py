# This example demonstrates how to analyze a shear wall with openings. It
# follows Section 10.5.3 of "Masonry Structures - Behavior and Design, 2nd
# Edition" by Robert G. Drysdale, Ahmad A. Hamid, and Lawrie R. Baker. The
# solution given in that text is obtained using an approximation method that
# isn't nearly as accurate as the finite element method, so some differences
# in the final results are expected.

# Import a few libraries from PyNite that we'll need
from PyNite.FEModel3D import FEModel3D
from PyNite.Mesh import RectangleMesh, RectOpening
from math import isclose

# Create a finite element model
model = FEModel3D()

# Set material properties for the wall (2 ksi masonry)
f_m = 2000        # Masonry compressive strength (psi)
E = 900*f_m/1000  # Masonry modulus of elasticity (ksi)
nu = 0.17         # Poisson's ratio for masonry
rho = 0.0000710   # Masonry density (kci)
model.add_material('Masonry', E, 0.4*E, nu, rho)

# Choose a desired mesh size. The program will try to stick to this as best as it can.
mesh_size = 6  # in

# Set the wall's dimensions
width = 26*12   # Wall overall width (in)
height = 16*12  # Wall overall height (in)
t = 8           # Masonry thickness (in)

# Generate the rectangular mesh
# The effects of cracked masonry can be modeled by adjusting the `ky_mod` factor. For this example
# uncracked masonry will be used to match the textbook problem.
model.add_rectangle_mesh('MSH1', mesh_size, width, height, t, 'Masonry', kx_mod=1, ky_mod=1,
                         origin=[0, 0, 0], plane='XY', element_type='Rect')

# Add a 4' wide x 12' tall door opening to the mesh
model.Meshes['MSH1'].add_rect_opening(name='Door 1', x_left=2*12, y_bott=0*12, width=4*12, height=12*12)

# Add a 4' wide x 4' tall window opening to the mesh
model.Meshes['MSH1'].add_rect_opening(name='Window 1', x_left=8*12, y_bott=8*12, width=4*12, height=4*12)

# Add another 4' wide x 4' tall window opening to the mesh
model.Meshes['MSH1'].add_rect_opening(name='Window 2', x_left=14*12, y_bott=8*12, width=4*12, height=4*12)

# Add another 4' wide x 12' tall door opening to the mesh
model.Meshes['MSH1'].add_rect_opening(name='Door 2', x_left=20*12, y_bott=0*12, width=4*12, height=12*12)

# Generate the mesh now that we've defined all the openings. Pynite automatically generates all
# meshes when analyzing, but generating it early with give us the chance to work with it prior to
# analysis.
model.Meshes['MSH1'].generate()

# Shear at the top of the wall
V = 100  # kip

# The shear at the top of the wall will be distributed evently to all the
# nodes at the top of the wall. Determine how many nodes are at the top of the
# wall.
n = len([node for node in model.Nodes.values() if isclose(node.Y, height)])
v = V/n

# Add supports and loads to the nodes
for node in model.Nodes.values():

    # Determine if the node is at the base of the wall
    if isclose(node.Y, 0):
        # Fix the base of the wall
        model.def_support(node.name, True, True, True, True, True, True)
    # Determine if the node is at the top of the wall
    elif isclose(node.Y, height):
        # Add out-of-plane support (provided by the diaphragm)
        model.def_support(node.name, False, False, True, False, False, False)
        # Add a seismic shear load to the top of the wall (applied by the diaphragm)
        model.add_node_load(node.name, 'FX', v, case='E')

# Add a load combination named 'Seismic' with a factor of 1.0 applied to any loads designated as
# 'E'.
model.add_load_combo('Seismic', {'E': 1.0})

# Analyze the model
model.analyze(log=True, check_statics=True)

# Render the model and plot the `Txy` shears.
# window = render_model(model, annotation_size=1, render_loads=True, deformed_shape=True,
#                       deformed_scale=200, color_map='Txy', scalar_bar=False,
#                       combo_name='Seismic', labels=False, screenshot='console')
from PyNite.Rendering import Renderer
renderer = Renderer(model)
renderer.combo_name = 'Seismic'
renderer.color_map = 'SX'
renderer.annotation_size = 1
renderer.deformed_shape = True
renderer.deformed_scale = 200
renderer.scalar_bar = True
renderer.render_model()

# Print the maximum displacement
d_max = max([node.DX['Seismic'] for node in model.Nodes.values()])
print('Max displacement: ', d_max, 'in')
print('Expected displacement from reference text: ', 7.623/E*t, 'in')
print('Wall rigidity: ', V/d_max, 'kips/in')
print('Expected rigidity from reference text: ', 1/7.623*E*t, 'kips/in')
