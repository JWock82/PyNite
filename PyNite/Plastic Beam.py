# Matrix Structural Analysis, 2nd Ed, Problem 8.6

from PyNite import FEModel3D
from PyNite.Section import SteelSection

# Create the model
plastic_beam = FEModel3D()

# Define a material
E = 29000  # ksi
G = 11200  # ksi
nu = 0.3
rho = 0.490/12**3  # kci
fy = 50  # ksi
plastic_beam.add_material('Stl_A992', E, G, nu, rho, fy)

# Define a cross-section
plastic_beam.add_steel_section('W12x65', 19.1, 20, 533, 1, 15, 96.8, 'Stl_A992')

# Add nodes
plastic_beam.add_node('N1', 0, 0, 0)
plastic_beam.add_node('N2', 8*12, 0, 0)
plastic_beam.add_node('N3', 24*12, 0, 0)

# Add supports
plastic_beam.def_support('N1', True, True, True, True, True, True)
plastic_beam.def_support('N3', False, True, True, False, False, False)

# Add a member
plastic_beam.add_member('M1', 'N1', 'N3', 'Stl_A992', 'W12x65')

# Add a load
plastic_beam.add_node_load('N3', 'FY', -0.0001, 'D')
plastic_beam.add_node_load('N2', 'FY', -0.3*325.7, 'Push')
plastic_beam.add_node_load('N3', 'FX', -1*325.7, 'Push')

# Add a load combination
plastic_beam.add_load_combo('1.4D', {'D':1.4})
plastic_beam.add_load_combo('Pushover', {'Push':0.01})

# Analyze the model
plastic_beam._not_ready_yet_analyze_pushover(log=True, check_stability=False, push_combo='Pushover', max_iter=30, sparse=True, combo_tags=None)

# Plot the moment diagram
# plastic_beam.members['M1'].plot_shear('Fy', '1.4D')
plastic_beam.members['M1'].plot_moment('Mz', '1.4D')
# plastic_beam.members['M1'].plot_deflection('dy', '1.4D')

# Render the model
# from PyNite.Visualization import Renderer
# rndr = Renderer(plastic_beam)
# rndr.combo_name = '1.4D'
# rndr.render_loads = True
# rndr.deformed_shape = True
# rndr.deformed_scale = 100
# rndr.render_model()