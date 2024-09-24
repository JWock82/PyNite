# Units used for the model in this example are inches and kips

# Import `FEModel3D` from `PyNite`
from PyNite import FEModel3D
from matplotlib import pyplot as plt
import numpy as np

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
simple_beam.add_member_pt_load('M1', 'Fy', 5, 4*12, 'D') # 5 kips Dead load
simple_beam.add_member_pt_load('M1', 'Fy', 8, 9*12, 'S') # 8 kips Live load
simple_beam.add_member_pt_load('M1', 'Fy', -8, 12*12, 'L') # 8 kips Live load

# Add load combinations
simple_beam.add_load_combo('1.4D', {'D':1.4})
simple_beam.add_load_combo('1.2D+1.6L', {'D':1.2, 'L':1.6})
simple_beam.add_load_combo('1.2D+1.6S', {'D':1.2, 'S':1.6})

# Analyze the beam and perform a statics check
simple_beam.analyze(check_statics=True)

# Visualization.render_model(simple_beam, annotation_size=10, deformed_shape=True, deformed_scale=30, render_loads=True, combo_name='1.2D+1.6L')

# Plot the shear diagram with all load cases and max/min envelope
x, M1 = simple_beam.members['M1'].moment_array("Mz", n_points=400, combo_name='1.4D')
_, M2 = simple_beam.members['M1'].moment_array("Mz", n_points=400, combo_name='1.2D+1.6L')
_, M3 = simple_beam.members['M1'].moment_array("Mz", n_points=400, combo_name='1.2D+1.6S')

max_envelope = np.maximum(np.maximum(M1, M2), M3)
min_envelope = np.minimum(np.minimum(M1, M2), M3)

plt.plot(x, np.zeros(len(x)), c="black", lw=3)
plt.plot(x, M1)
plt.plot(x, M2)
plt.plot(x, M3)
plt.plot(x, max_envelope, alpha=0.3, c="green", lw=5)
plt.plot(x, min_envelope, alpha=0.3, c="red", lw=5)