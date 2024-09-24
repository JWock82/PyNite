# A First Course in the Finite Element Method, 4th Edition
# Daryl L. Logan
# Problem 5.58
# Units for this model are kips and inches

# Import 'FEModel3D' and 'Visualization' from 'PyNite'
from PyNite import FEModel3D
from PyNite import Visualization

# Create a new model
frame = FEModel3D()

# Define the nodes
frame.add_node('N1', 0, 0, 0)
frame.add_node('N2', 10*12, 0, 0)
frame.add_node('N3', 10*12, 0, -10*12)
frame.add_node('N4', 10*12, -20*12, -10*12)

# Define the supports
frame.def_support('N1', True, True, True, True, True, True)
frame.def_support('N4', True, True, True, True, True, True)

# Create members (all members will have the same properties in this example)
J = 100
Iy = 200
Iz = 1000
A = 100
frame.add_section('MySection', A, Iy, Iz, J)

# Define a material
E = 30000
G = 10000
nu = 0.3
rho = 2.836e-4
frame.add_material('Steel', E, G, nu, rho)

frame.add_member('M12', 'N1', 'N2', 'Steel', 'MySection')
frame.add_member('M23', 'N2', 'N3', 'Steel', 'MySection')
frame.add_member('M34', 'N3', 'N4', 'Steel', 'MySection')

# Add nodal loads
frame.add_node_load('N2', 'FY', -5)
frame.add_node_load('N2', 'MX', -100*12)
frame.add_node_load('N3', 'FZ', 40)

# Analyze the frame
frame.analyze()

print('Calculated results: ', frame.nodes['N2'].DY, frame.nodes['N3'].DZ)
print('Expected results: ', -0.063, 1.825)

# Render the model for viewing
Visualization.render_model(frame, annotation_size=5, deformed_shape=True, deformed_scale=40, render_loads=True)