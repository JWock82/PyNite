# A First Course in the Finite Element Method, 4th Edition
# Daryl L. Logan
# Example 5.8
# Units for this model are kips and inches

# Import 'FEModel3D' and 'Visualization' from 'PyNite'
from PyNite import FEModel3D
from PyNite.Rendering import Renderer

# Create a new model
frame = FEModel3D()

# Define the nodes
frame.add_node('N1', 0, 0, 0)
frame.add_node('N2', -100, 0, 0)
frame.add_node('N3', 0, 0, -100)
frame.add_node('N4', 0, -100, 0)

# Define the supports
frame.def_support('N2', True, True, True, True, True, True)
frame.def_support('N3', True, True, True, True, True, True)
frame.def_support('N4', True, True, True, True, True, True)

# Create members (all members will have the same properties in this example)
J = 50
Iy = 100
Iz = 100
A = 10
frame.add_section('MySection', A, Iy, Iz, J)

# Define a material
E = 30000
G = 10000
nu = 0.3
rho = 2.836e-4
frame.add_material('Steel', E, G, nu, rho)

frame.add_member('M1', 'N2', 'N1', 'Steel', 'MySection')
frame.add_member('M2', 'N3', 'N1', 'Steel', 'MySection')
frame.add_member('M3', 'N4', 'N1', 'Steel', 'MySection')

# Add nodal loads
frame.add_node_load('N1', 'FY', -50)
frame.add_node_load('N1', 'MX', -1000)

# Analyze the model
frame.analyze(check_statics=True)

# Render the deformed shape
rndr = Renderer(frame)
rndr.annotation_size = 5
rndr.render_loads = True
rndr.deformed_shape = True
rndr.deformed_scale = 100
rndr.render_model()

# Print the node 1 displacements
print('Node 1 deformations:')
print('Calculated values: ', frame.nodes['N1'].DX, frame.nodes['N1'].DY, frame.nodes['N1'].DZ, frame.nodes['N1'].RX, frame.nodes['N1'].RY, frame.nodes['N1'].RZ)
print('Expected values: ', 7.098e-5, -0.014, -2.352e-3, -3.996e-3, 1.78e-5, -1.033e-4)