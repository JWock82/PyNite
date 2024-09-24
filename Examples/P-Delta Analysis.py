from PyNite import FEModel3D
from math import tan

# Create the cantilever model
cantilever = FEModel3D()

# Define the column and its section properties
L = 20 # ft
H = 5 # kips lateral load
P = 100 # kips axial load
I = 100/12**4 # moment of inertia (ft^4)

# Define a material
G = 11200*12**2 # shear modulus (ksf)
E = 29000*12**2 # modulus of elasticity (ksf)
nu = 0.3
rho = 0.490  # kcf
cantilever.add_material('Steel', E, G, nu, rho)

# Add nodes along the length of the column to capture P-little-delta effects
num_nodes = 6
for i in range(num_nodes):
    # Add nodes
    cantilever.add_node('N' + str(i+1), 0, i*L/(num_nodes - 1), 0)

# Add the section and member
cantilever.add_section('MySection', 10/12**2, I, I, 200/12**4)
cantilever.add_member('M1', 'N1', 'N6', 'Steel',  'MySection')

# Add a fixed support at the base of the column
cantilever.def_support('N1', True, True, True, True, True, True)

# Add a -10 kip axial load to the top of the column
cantilever.add_node_load('N6', 'FY', -P)

# Add a 5 kip lateral load to the top of the column
cantilever.add_node_load('N6', 'FX', H)

# Perform 2nd order analysis
cantilever.analyze_PDelta()

from PyNite.Visualization import Renderer
renderer = Renderer(cantilever)
renderer.annotation_size = 0.5
renderer.window_width = 750
renderer.window_height = 400
renderer.deformed_shape = True
renderer.render_model()

# The moment at the base of the column
calculated_moment = cantilever.nodes['N1'].RxnMZ['Combo 1']
calculated_moment2 = cantilever.members['M1'].plot_moment('Mz', 'Combo 1', 100)

# the deflection at the top of the column
calculated_displacement = cantilever.nodes['N6'].DX['Combo 1']*12

# Calculate the AISC benchmark problem solution:
alpha = (P*L**2/(E*I))**0.5
Mmax = H*L*(tan(alpha)/alpha)
ymax = H*L**3/(3*E*I)*(3*(tan(alpha)-alpha)/alpha**3)

# Compare the calculation results
print('Expected Moment, M_max: ', Mmax)
print('Calculated Moment, M_calc: ', calculated_moment)
print('')
print('Expected Displacement, y_max: ', ymax*12)
print('Calculated Displacement, y_calc:', calculated_displacement)