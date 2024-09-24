# Example of a basic 2D moment frame with gravity and lateral loads
# Units used for the model in this example are inches and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

# Create a new finite element model
MomentFrame = FEModel3D()

# Add nodes (frame is 15 ft wide x 12 ft tall)
MomentFrame.add_node('N1', 0, 0, 0)
MomentFrame.add_node('N2', 0, 12*12, 0)
MomentFrame.add_node('N3', 15*12, 12*12, 0)
MomentFrame.add_node('N4', 15*12, 0*12, 0)

# Define column properties (use W10x33 from the AISC Manual):
Iy = 36.6 # in^4
Iz = 171 # in^4
J = 0.58 # in^4
A = 9.71 # in^2
MomentFrame.add_section('W10x33', A, Iy, Iz, J)

# Define a material
E = 29000 # ksi
G = 11200 # ksi
nu = 0.3  # Poisson's ratio
rho = 0.490/12**3  # Density (kci)
MomentFrame.add_material('Steel', E, G, nu, rho)

# Define the columns
MomentFrame.add_member('Col1', 'N1', 'N2', 'Steel', 'W10x33')
MomentFrame.add_member('Col2', 'N4', 'N3', 'Steel', 'W10x33')

# Define beam properties (Use W8x24)
Iy = 18.3 # in^4
Iz = 82.7 # in^4
J = 0.346 # in^4
A = 7.08 # in^2
MomentFrame.add_section('W8x24', A, Iy, Iz, J)

# Define the beams
MomentFrame.add_member('Beam', 'N2', 'N3', 'Steel', 'W8x24')

# Provide fixed supports at the bases of the columns
MomentFrame.def_support('N1', support_DX=True, support_DY=True, support_DZ=True, support_RX=True, support_RY=True, support_RZ=True)
MomentFrame.def_support('N4', support_DX=True, support_DY=True, support_DZ=True, support_RX=True, support_RY=True, support_RZ=True)

# Add self weight dead loads to the frame
# Note that we could leave 'x1' and 'x2' undefined below and it would default to the full member length
# Note also that the direction uses lowercase notations to indicate member local coordinate systems
MomentFrame.add_member_dist_load('Beam', Direction='Fy', w1=-0.024/12, w2=-0.024/12, x1=0, x2=15*12, case='D')
MomentFrame.add_member_dist_load('Col1', Direction='Fx', w1=-0.033/12, w2=-0.033/12, x1=0, x2=12*12, case='D')
MomentFrame.add_member_dist_load('Col2', Direction='Fx', w1=-0.033/12, w2=-0.033/12, x1=0, x2=12*12, case='D')

# Add a nodal wind load of 10 kips at the left side of the frame
# Note that the direction uses uppercase notation to indicate model global coordinate system
MomentFrame.add_node_load('N2', Direction='FX', P=10, case='W')

# Create two load combinations
MomentFrame.add_load_combo('1.2D+1.0W', factors={'D':1.2, 'W':1.0})
MomentFrame.add_load_combo('0.9D+1.0W', factors={'D':0.9, 'W':1.0})

# Perform a P-Big Delta analysis on the frame
# Note that to capture P-little delta effects the members should idealy be broken into three segments each
MomentFrame.analyze_PDelta(log=True)

# A first-order analysis could be done with the either of the following lines instead
# The first option ignores P-Delta effects
# The second option ignores P-Delta effects and tension/compression only elements and supports
# MomentFrame.analyze()
# MomentFrame.analyze_linear(log=True)

# Display the deformed shape of the structure magnified 50 times with the text height 5 model units (inches) high
from PyNite import Visualization
Visualization.render_model(MomentFrame, annotation_size=5, deformed_shape=True, deformed_scale=50, combo_name='1.2D+1.0W')

# Plot the moment diagram for the beam
MomentFrame.members['Beam'].plot_moment('Mz', combo_name='1.2D+1.0W')

# Plot the deflection of the column
MomentFrame.members['Col1'].plot_deflection('dy', combo_name='1.2D+1.0W')

# Find the maximum shear in the first column
print('Column Shear Force:', MomentFrame.members['Col1'].max_shear('Fy', combo_name='1.2D+1.0W'), 'kip')

# Find the frame lateral drift
# Note that the deflections are stored as a dictionary in the node object by load combination, so [] are used instead of () below
print('Frame Lateral Drift:', MomentFrame.nodes['N2'].DX['1.2D+1.0W'], 'in')

# Find the maximum uplift reaction
print('Left Support Y Reaction:', MomentFrame.nodes['N1'].RxnFY['0.9D+1.0W'], 'kip')