# Example of a basic 2D tension-only braced frame with gravity and lateral
# loads. Units used for the model in this example are inches and kips.

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

# Create a new finite element model
braced_frame = FEModel3D()

# Add nodes (frame is 15 ft wide x 12 ft tall)
braced_frame.add_node('N1', 0, 0, 0)
braced_frame.add_node('N2', 0, 12*12, 0)
braced_frame.add_node('N3', 15*12, 12*12, 0)
braced_frame.add_node('N4', 15*12, 0*12, 0)

# Define column properties (use W10x33 from the AISC Manual) and add section to the model:
Iy = 36.6 # in^4
Iz = 171  # in^4
J = 0.58  # in^4
A = 9.71  # in^2
braced_frame.add_section('W10x33', A, Iy, Iz, J)

# Define a material
E = 29000 # Young's modulus (ksi)
G = 11200 # Shear modulus (ksi)
nu = 0.3  # Poisson's ratio
rho = 0.490/12**2  # Density (kci)
braced_frame.add_material('Steel', E, G, nu, rho)

# Define the columns
braced_frame.add_member('Col1', 'N1', 'N2', 'Steel', 'W10x33')
braced_frame.add_member('Col2', 'N4', 'N3', 'Steel', 'W10x33')

# Define beam properties (Use W8x24) and add section to the model
Iy = 18.3 # in^4
Iz = 82.7 # in^4
J = 0.346 # in^4
A = 7.08 # in^2
braced_frame.add_section('W8x24', A, Iy, Iz, J)

# Define the beams
braced_frame.add_member('Beam', 'N2', 'N3', 'Steel', 'W8x24')
braced_frame.def_releases('Beam', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

# Define the brace properties
# We'll use a section with L/r <= 300 which is a common rule of thumb for
# tension members. We'll use L4x4x1/4. 
Iy = 3 # in^4
Iz = 3 # in^4
J = 0.0438 # in^4
A = 1.94 # in^2
braced_frame.add_section('L4x4x1/4', A, Iy, Iz, J)

# Define a brace (tension and compression - both ways)
braced_frame.add_member('Brace1', 'N1', 'N3', 'Steel', 'L4x4x1/4')

# Let's add spring supports to the base of the structure. We'll add a couple of
# extra nodes at the base of the structure that will receive the springs. The
# lengths of these springs is irrelevant, since we're defining a spring constant
# that is independent of the spring's length. Only the direction of the spring
# matters as it defines the direction of the spring's stiffness. The nodes will
# be directly below N1 and N4 in the Y-direction.
braced_frame.add_node('N1s', 0, -2*12, 0)
braced_frame.add_node('N4s', 15*12, -2*12, 0)

braced_frame.add_spring('Spring1','N1', 'N1s', 10000, tension_only=True,
                        comp_only=False)
braced_frame.add_spring('Spring2','N4', 'N4s', 10000, tension_only=False,
                        comp_only=True) # The structure would be unstable if
                                        # this was tension only

# Release the brace ends to form an axial member
braced_frame.def_releases('Brace1', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

# Springs only carry axial loads, nothing else, so we'll need to stabilize
# the column bases in the other directions. The column bases will be
# supported by the springs vertically. For the other directions (horizontally
# and about the Y-axis) we'll need to provide supports.
braced_frame.def_support('N1', support_DX=True, support_DZ=True, support_RY=True)
braced_frame.def_support('N4', support_DX=True, support_DZ=True, support_RY=True)

# Fix the nodes supporting the bottoms of the springs. Note that even though
# we're fixing these nodes, the only reactions the supports will carry will
# be in the Y-direction, due to the fact that the spring only has stiffness in
# that direction. We fix the node so that it's not free to spin or translate
# in the other directions however. If we didn't the node would be unstable and
# the model would crash. PyNite is unforgiving in this regard. Every degree of
# freedom (3 translations and 3 rotations) at every node must be stabilized so
# it's not free to move infinitely.
braced_frame.def_support('N1s', support_DX=True, support_DY=True, support_DZ=True, support_RX=True, support_RY=True, support_RZ=True)
braced_frame.def_support('N4s', support_DX=True, support_DY=True, support_DZ=True, support_RX=True, support_RY=True, support_RZ=True)

# Stabilize the frame in the global Z-direction so it doesn't tip over
# out-of-plane.
braced_frame.def_support('N2', support_DZ=True)
braced_frame.def_support('N3', support_DZ=True)

# Add self weight dead loads to the frame.
# Note that we could leave 'x1' and 'x2' undefined below and it would default
# to the full member length. Note also that the direction uses lowercase
# notations to indicate member local coordinate systems. Brace loads have been
# neglected.
braced_frame.add_member_dist_load('Beam', Direction='Fy', w1=-0.024/12,
                              w2=-0.024/12, x1=0, x2=15*12, case='D')
braced_frame.add_member_dist_load('Col1', Direction='Fx', w1=-0.033/12,
                              w2=-0.033/12, x1=0, x2=12*12, case='D')
braced_frame.add_member_dist_load('Col2', Direction='Fx', w1=-0.033/12,
                              w2=-0.033/12, x1=0, x2=12*12, case='D')

# Add nodal wind loads of 25 kips to each side of the frame. Note that the
# direction uses uppercase notation to indicate model global coordinate
# system.
braced_frame.add_node_load('N2', Direction='FX', P=25, case='W')
braced_frame.add_node_load('N3', Direction='FX', P=25, case='W')

# Create load combinations
# Note that the load combination '1.4D' has no lateral load, but does have
# gravity load. The gravity load forces the tension only spring to receive
# minor compression, which causes it to be deactivated on the first iteration.
# Once deactivated the model is unstable and an exception is thrown. This is
# normal and correct behavior. Load combination '1.4D' has been commented out,
# but you can uncomment it to see for yourself what happens.

# braced_frame.add_load_combo('1.4D', factors={'D':1.4})
braced_frame.add_load_combo('1.2D+1.0W', factors={'D':1.2, 'W':1.0})
braced_frame.add_load_combo('0.9D+1.0W', factors={'D':0.9, 'W':1.0})

# Analyze the braced frame
# P-Delta analysis could also be performed using braced_frame.analyze_PDelta().
# Generally, P-Delta analysis will have little effect on a model of a braced
# frame, as there is usually very little bending moment in the members.
braced_frame.analyze()

# Display the deformed shape of the structure magnified 50 times with the text
# height 5 model units (inches) high.
from PyNite.Rendering import Renderer
rndr = Renderer(braced_frame)
rndr.annotation_size = 5
rndr.deformed_shape = True
rndr.deformed_scale = 50
rndr.combo_name = '1.2D+1.0W'
rndr.window_width = 750
rndr.window_height = 750
rndr.render_loads = True
rndr.render_model()

# We should see upward displacement at N1 and downward displacement at N4 if
# our springs worked correctly
print('N1 displacement in Y =', braced_frame.nodes['N1'].DY['1.2D+1.0W'])
print('N4 displacement in Y =', braced_frame.nodes['N4'].DY['1.2D+1.0W'])
