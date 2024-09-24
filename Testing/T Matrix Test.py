# Matrix Structural Analysis, 2nd Edition
# William McGuire, Richard H. Gallagher, Ronald D. Ziemian
# Example 5.4
# Units for this model are kilonewtons and meters

# Import 'FEModel3D' and 'Visualization' from 'PyNite'
from PyNite import FEModel3D
from PyNite import Visualization

# Create a new model
truss = FEModel3D()

# Define the nodes
truss.add_node('a', 2, 4, 8)
truss.add_node('b', 0, 0, 0)
truss.add_node('c', 8, 0, 0)
truss.add_node('d', 8, 6, 0)
truss.add_node('e', 0, 6, 0)

# Define the supports - All fixed - End releases will be used to create pins
truss.def_support('b', True, True, True, True, True, True)
truss.def_support('c', True, True, True, True, True, True)
truss.def_support('d', True, True, True, True, True, True)
truss.def_support('e', True, True, True, True, True, True)

# Define properties common to all members
E = 200000000  # kN/m^2
G = 0.4*E
truss.add_material('Steel',E, G, 0.3, 78.5)

#Define some section properties
J = 100
Iy = 100
Iz = 100

# Define members
truss.add_section('TrussSection1', 20*10**3/1000**2, Iy, Iz, J)
truss.add_member('ab', 'a', 'b', 'Steel', 'TrussSection1')

truss.add_section('TrussSection2', 30*10**3/1000**2, Iy, Iz, J)
truss.add_member('ac', 'a', 'c', 'Steel', 'TrussSection2')

truss.add_section('TrussSection3', 40*10**3/1000**2, Iy, Iz, J)
truss.add_member('ad', 'a', 'd', 'Steel', 'TrussSection3')

truss.add_member('ae', 'a', 'e', 'Steel', 'TrussSection2')

# Define member end releases
truss.def_releases('ab', False, False, False, False, True, True,
                           False, False, False, False, True, True)
truss.def_releases('ac', False, False, False, False, True, True,
                           False, False, False, False, True, True)
truss.def_releases('ad', False, False, False, False, True, True,
                           False, False, False, False, True, True)
truss.def_releases('ae', False, False, False, False, True, True,
                           False, False, False, False, True, True)

# Add nodal loads
truss.add_node_load('a', 'FY', 600)
truss.add_node_load('a', 'FX', 200)
truss.add_node_load('a', 'FZ', -800)

# Analyze the model
truss.analyze(False)

# Render the model
Visualization.RenderModel(truss, text_height=0.2, deformed_shape=True, deformed_scale=1000, render_loads=True)

# Print results
a = truss.nodes['a']
b = truss.nodes['b']
c = truss.nodes['c']
d = truss.nodes['d']
e = truss.nodes['e']

print('Node b Reactions: ', b.RxnFX, b.RxnFY, b.RxnFZ, b.RxnMX, b.RxnMY, b.RxnMZ)
print('Node c Reactions: ', c.RxnFX, c.RxnFY, c.RxnFZ, c.RxnMX, c.RxnMY, c.RxnMZ)
print('Node e Reactions: ', d.RxnFX, d.RxnFY, d.RxnFZ, d.RxnMX, d.RxnMY, d.RxnMZ)
print('Node e Reactions: ', e.RxnFX, e.RxnFY, e.RxnFZ, e.RxnMX, e.RxnMY, e.RxnMZ)
print('Node a Displacements: ', a.DX, a.DY, a.DZ)

# Results printed above matched the textbook problem results perfectly