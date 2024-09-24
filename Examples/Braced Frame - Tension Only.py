# Example of a basic 2D tension-only braced frame with gravity and lateral
# loads. Units used for the model in this example are inches and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

# Create a new finite element model
braced_frame = FEModel3D()

# Add nodes (frame is 15 ft wide x 12 ft tall)
braced_frame.add_node('N1', 0, 0, 0)
braced_frame.add_node('N2', 0, 12*12, 0)
braced_frame.add_node('N3', 15*12, 12*12, 0)
braced_frame.add_node('N4', 15*12, 0*12, 0)

# Define column properties (use W10x33 from the AISC Manual):
Iy = 36.6 # in^4
Iz = 171 # in^4
J = 0.58 # in^4
A = 9.71 # in^2
braced_frame.add_section('W10x33', A, Iy, Iz, J)

# Define a material
E = 29000 # ksi
G = 11400 # ksi
nu = 0.3  # Poisson's ratio
rho = 0.490/12**2  # Density (kci)
braced_frame.add_material('Steel', E, G, nu, rho)

# Define the columns
braced_frame.add_member('Col1', 'N1', 'N2', 'Steel', 'W10x33')
braced_frame.add_member('Col2', 'N4', 'N3', 'Steel', 'W10x33')

# Define beam properties (Use W8x24)
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

# Define the braces
braced_frame.add_member('Brace1', 'N1', 'N3', 'Steel', 'L4x4x1/4', tension_only=True)
braced_frame.add_member('Brace2', 'N4', 'N2', 'Steel', 'L4x4x1/4', tension_only=True)

# Release the brace ends to form an axial member
braced_frame.def_releases('Brace1', Ryi=True, Rzi=True, Ryj=True, Rzj=True)
braced_frame.def_releases('Brace2', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

# Provide pinned supports at the bases of the columns (also restrained about
# the Y-axis for stability)
braced_frame.def_support('N1', support_DX=True, support_DY=True, support_DZ=True,
                          support_RY=True)
braced_frame.def_support('N4', support_DX=True, support_DY=True, support_DZ=True,
                          support_RY=True)

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
# gravity load. The gravity load forces the braces to receive minor
# compression, which causes them to be deactivated on the first iteration.
# Once deactivated the model is unstable and an exception is thrown. This is
# normal and correct behavior. Some programs allow for slight compression in
# tension-only members to avoid this problem. Load combination '1.4D' ahs been
# commented out, but you can uncomment it to see for yourself what happens.

# braced_frame.add_load_combo('1.4D', factors={'D':1.4})
braced_frame.add_load_combo('1.2D+1.0W', factors={'D':1.2, 'W':1.0})
braced_frame.add_load_combo('0.9D+1.0W', factors={'D':0.9, 'W':1.0})

# Analyze the braced frame
# P-Delta analysis could also be performed using braced_frame.analyze_PDelta().
# Generally, P-Delta analysis will have little effect on a model of a braced
# frame, as there is usually very little bending moment in the members.
braced_frame.analyze()

# Display the deformed shape of the structure magnified 50 times with the text
# height 5 model units (inches) high
from PyNite.Visualization import Renderer
rndr = Renderer(braced_frame)
rndr.annotation_size = 5
rndr.deformed_shape = True
rndr.deformed_scale = 50
rndr.combo_name = '1.2D+1.0W'
rndr.window_width = 750
rndr.window_height = 750
rndr.render_model()

# Plot the axial load diagrams for the braces. We should see no compression on
# 'Brace2' and 64 kips on 'Brace1' if the tension-only analysis worked
# correctly.
braced_frame.members['Brace1'].plot_axial(combo_name='1.2D+1.0W')
braced_frame.members['Brace2'].plot_axial(combo_name='1.2D+1.0W')

# Report the frame reactions for the load combination '1.2D+1.0W'. We should
# see a -50 kip horizontal reaction at node 'N1', and a zero kip reaction at
# node 'N4' if the tension-only analysis worked correctly. Similarly, we
print('1.2D+1.0W: N1 reaction FY =',
      '{:.3f}'.format(braced_frame.nodes['N1'].RxnFY['1.2D+1.0W']), 'kip')
print('1.2D+1.0W: N1 reaction FX =',
      '{:.3f}'.format(braced_frame.nodes['N1'].RxnFX['1.2D+1.0W']), 'kip')
print('1.2D+1.0W: N4 reaction FY =',
      '{:.3f}'.format(braced_frame.nodes['N4'].RxnFY['1.2D+1.0W']), 'kip')
print('1.2D+1.0W: N4 reaction FX =',
      '{:.3f}'.format(braced_frame.nodes['N4'].RxnFX['1.2D+1.0W']), 'kip')