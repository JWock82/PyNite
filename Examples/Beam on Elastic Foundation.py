# -*- coding: utf-8 -*-
"""
This example shows how to generate a beam on elastic foundation by using spring supports. All units
in this model are expressed in terms of kips (force) and inches (length).
"""

# Create a new model for the beam
from PyNite import FEModel3D
boef = FEModel3D()

# Assign the length of the beam
L = 25*12  # in

# The beam will be modeled as a single element with nodes at a 12 inch spacing. A support spring
# will be added to each node.
num_nodes = round(L/12) + 1
spacing = L/round(L/12)

# Define the spring stiffness
ks = 22.5  # k/in

# Generate the nodes
for i in range(num_nodes):

    # Add nodes spaced at 15"
    boef.add_node('N' + str(i + 1), i*spacing, 0, 0)

    # Add supports to the nodes
    if i == 0 or i == num_nodes - 1:
        # Pinned supports at beam ends
        boef.def_support('N' + str(i + 1), True, True, True, True, False, False)
    else:
        # Spring supports at all other locations
        boef.def_support_spring('N' + str(i + 1), 'DY', ks, '-')

# Define member material properties
E = 29000   # ksi
G = 11200   # ksi
boef.add_material('Steel', E, G, 0.3, 490/1000/12**3)

# Define section properties (W8x35) and add a section to the model
A = 10.3    # in^2
Iz = 127    # in^4 (strong axis)
Iy = 42.6   # in^4 (weak axis)
J = 0.769   # in^4

boef.add_section('W8x35', A, Iy, Iz, J)

# Define the member
boef.add_member('M1', 'N1', 'N' + str(num_nodes), 'Steel', 'W8x35')

# In the next few lines no load case or load combination is being specified. When this is the case,
# PyNite internally creates a default load case ('Case 1') and a default load combination
# ('Combo 1'). This can be handy for quick calculations that don't need the added complexity of
# using load combinations.

# Add a point load at midspan
if num_nodes % 2 == 0:  # Checks if there is a node at midspan or a gap between nodes
    boef.add_member_pt_load('M1', 'Fy', -40, L/2)
else:
    mid_node = 'N' + str(int(num_nodes/2))
    boef.add_node_load(mid_node, 'FY', -40)

# Analyze the model. PyNite's standard solver is most appropriate or this model since there are
# non-linear features (compression-only springs) but no large axial forces that would cause P-Delta
# effects.
boef.analyze(log=True, check_statics=True)

# Render the mdoel with the deformed shape using PyNite's buit-in renderer
from PyNite.Rendering import Renderer
renderer = Renderer(boef)
renderer.annotation_size = 1.5
renderer.deformed_shape = True
renderer.window_width = 750   # Measured in pixels
renderer.window_height = 250 
renderer.render_model()

# Find and print the largest displacement
d_min = boef.members['M1'].min_deflection('dy')
print('Minimum deflection: ', round(d_min, 4), 'in')

# Alternatively the line below could be used to get the largest nodal displacement in the model
# It will be slightly less since there is no node at the midpoint of the member.
# d_min = min([node.DY['Combo 1'] for node in boef.nodes.values()])

# Find and print the minimum moment
M_min = boef.members['M1'].min_moment('Mz')
M_max = boef.members['M1'].max_moment('Mz')
print('Minimum moment: ', round(M_min), 'k-in')
print('Maximum moment: ', round(M_max), 'k-in')