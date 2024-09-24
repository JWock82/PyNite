# Example of a simply supported beam with a point load.
# Units used for the model in this example are inches and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D

# Import 'Visualization' for rendering the model
from PyNite import Visualization

# Create a new finite element model
simple_beam = FEModel3D()

# Add nodes (14 ft apart)
simple_beam.add_node('N1', 0, 0, 0)
simple_beam.add_node('N2', 14*12, 0, 0)

# Define a material
E = 29000       # Modulus of elasticity (ksi)
G = 11200       # Shear modulus of elasticity (ksi)
nu = 0.3        # Poisson's ratio
rho = 2.836e-4  # Density (kci)
simple_beam.add_material('Steel', E, G, nu, rho)

# Add a section with the following properties:
# Iy = 100 in^4, Iz = 150 in^4, J = 250 in^4, A = 20 in^2
simple_beam.add_section('MySection', 20, 100, 150, 250)

#Add member
simple_beam.add_member('M1', 'N1', 'N2', 'Steel', 'MySection')

# Provide simple supports
simple_beam.def_support('N1', True, True, True, True, False, False)  # Constrained for torsion at 'N1'
simple_beam.def_support('N2', True, True, True, False, False, False) # Not constrained for torsion at 'N2'

# Add a downward point load of 5 kips at the midspan of the beam
simple_beam.add_member_pt_load('M1', 'Fy', -5, 7*12, 'D') # 5 kips Dead load
simple_beam.add_member_pt_load('M1', 'Fy', -8, 7*12, 'L') # 8 kips Live load

# Add load combinations
simple_beam.add_load_combo('1.4D', {'D':1.4})
simple_beam.add_load_combo('1.2D+1.6L', {'D':1.2, 'L':1.6})

# Analyze the beam and perform a statics check
simple_beam.analyze(check_statics=True)

# Render the model
from PyNite.Visualization import Renderer
renderer = Renderer(simple_beam)
renderer.annotation_size = 10
renderer.deformed_shape = True
renderer.deformed_scale = 30
renderer.render_loads = True
renderer.combo_name = '1.2D+1.6L'
renderer.render_model()

# Print the shear, moment, and deflection diagrams
simple_beam.members['M1'].plot_shear('Fy', '1.2D+1.6L')
simple_beam.members['M1'].plot_moment('Mz', '1.2D+1.6L')
simple_beam.members['M1'].plot_deflection('dy', '1.2D+1.6L')

# Print reactions at each end of the beam
print('Left Support Reaction:', simple_beam.nodes['N1'].RxnFY['1.2D+1.6L'], 'kip')
print('Right Support Reacton:', simple_beam.nodes['N2'].RxnFY['1.2D+1.6L'], 'kip')

# Print the max/min shears and moments in the beam
print('Maximum Shear:', simple_beam.members['M1'].max_shear('Fy', '1.2D+1.6L'), 'kip')
print('Minimum Shear:', simple_beam.members['M1'].min_shear('Fy', '1.2D+1.6L'), 'kip')
print('Maximum Moment:', simple_beam.members['M1'].max_moment('Mz', '1.2D+1.6L')/12, 'kip-ft')
print('Minimum Moment:', simple_beam.members['M1'].min_moment('Mz', '1.2D+1.6L')/12, 'kip-ft')

# Print the max/min deflections in the beam
print('Maximum Deflection:', simple_beam.members['M1'].max_deflection('dy', '1.2D+1.6L'), 'in')
print('Minimum Deflection:', simple_beam.members['M1'].min_deflection('dy', '1.2D+1.6L'), 'in')

# The following lines can be uncommented to create a PDF report. Follow the instructions on the
# wiki under "Generating PDF Reports" to prevent errors. The report will be output to the PyNite
# folder unless the 'output_path' variable below is modified.

from PyNite import Reporting
Reporting.create_report(simple_beam, output_filepath='./Pynite Report.pdf', node_table=False, plate_table=False, plate_corner_forces=False, plate_center_forces=False, plate_corner_membrane=False, plate_center_membrane=False)