# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import enum
import unittest
from PyNite import FEModel3D
import math
import sys
from io import StringIO

class Test_AISC_Benchmark(unittest.TestCase):
    """
    AISC's benchmark problems used to determine if a second order analysis procedure is rigorous
    enough.
    """

    def setUp(self):

        # Suppress printed output temporarily
        # sys.stdout = StringIO()
        pass

    def tearDown(self):
        
        # Reset the print function to normal
        # sys.stdout = sys.__stdout__
        pass

    def test_AISC_benchmark_old(self):

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
        cantilever.add_section('BeamSection', 10/12**2, I, I, 200/12**4)

        # Add nodes along the length of the column to capture P-little-delta effects
        num_nodes = 6
        for i in range(num_nodes):
            # Add nodes
            cantilever.add_node('N' + str(i+1), 0, i*L/(num_nodes - 1), 0)
        
        # Add the member
        cantilever.add_member('M1', 'N1', 'N6', 'Steel', 'BeamSection')

        # Add a fixed support at the base of the column
        cantilever.def_support('N1', True, True, True, True, True, True)

        # Add a -10 kip axial load to the top of the column
        cantilever.add_node_load('N6', 'FY', -P)

        # Add a 5 kip lateral load to the top of the column
        cantilever.add_node_load('N6', 'FX', H)

        # Perform 2nd order analysis
        cantilever.analyze_PDelta()

        # from PyNite.Visualization import Renderer
        # renderer = Renderer(cantilever)
        # renderer.annotation_size = 0.5
        # renderer.window_width = 750
        # renderer.window_height = 750
        # renderer.deformed_shape = True
        # renderer.render_model()

        # The moment at the base of the column
        calculated_moment = cantilever.nodes['N1'].RxnMZ['Combo 1']

        # the deflection at the top of the column
        calculated_displacement = cantilever.nodes['N6'].DX['Combo 1']*12

        # Calculate the AISC benchmark problem solution:
        alpha = (P*L**2/(E*I))**0.5
        Mmax = H*L*(math.tan(alpha)/alpha)
        ymax = H*L**3/(3*E*I)*(3*(math.tan(alpha)-alpha)/alpha**3)
        
        # Compare the calculation results
        self.assertAlmostEqual(calculated_moment/Mmax, 1.0, 1)
        self.assertAlmostEqual(calculated_displacement/(ymax*12), 1.0, 1)
    
    def test_AISC_benchmark_case1(self):

        column = FEModel3D()

        column.add_node('N1', 0, 0, 0)
        column.add_node('N2', 0, 28, 0)
        
        # column.add_node('N3', 0, 28*1/4, 0)
        column.add_node('N4', 0, 28*1/2, 0)
        # column.add_node('N5', 0, 28*3/4, 0)

        column.def_support('N1', True, True, True, False, True, False)
        column.def_support('N2', True, False, True, False, False, False)

        # Define a material
        E = 29000*12**2  # ksf
        G = 11200*12**2  # ksf
        nu = 0.3
        rho = 0.490  # kcf
        column.add_material('Steel', E, G, nu, rho)

        # Define section properties
        Iy = 51.4/12**4  # ft^4
        Iz = 484/12**4   # ft^4
        A = 14.1/12**2   # ft^2
        J = 1.45/12**4   # ft^4
        column.add_section('ColumnSection', A, Iy, Iz, J)

        column.add_member('M1', 'N1', 'N2', 'Steel', 'ColumnSection')
        
        column.add_member_dist_load('M1', 'Fy', -0.200, -0.200, case='P1')
        column.add_member_dist_load('M1', 'Fy', -0.200, -0.200, case='P2')
        column.add_member_dist_load('M1', 'Fy', -0.200, -0.200, case='P3')
        column.add_member_dist_load('M1', 'Fy', -0.200, -0.200, case='P4')

        column.add_node_load('N2', 'FY', -150, 'P2')
        column.add_node_load('N2', 'FY', -300, 'P3')
        column.add_node_load('N2', 'FY', -450, 'P4')

        column.add_load_combo('Combo 1', {'P1':1})
        column.add_load_combo('Combo 2', {'P2':1})
        column.add_load_combo('Combo 3', {'P3':1})
        column.add_load_combo('Combo 4', {'P4':1})

        # from PyNite.Visualization import Renderer
        # renderer = Renderer(column)
        # renderer.annotation_size = 1
        # renderer.combo_name = 'Combo 2'
        # renderer.render_model()

        column.analyze_PDelta()

        Mmid_calculated = []
        dmid_calculated = []
        for combo in column.load_combos.values():
            Mmid_calculated.append(-column.members['M1'].moment('Mz', 14, combo.name)*12)
            dmid_calculated.append(-column.members['M1'].deflection('dy', 14, combo.name)*12)

        # Expected results per AISC
        Mmid_expected = [235, 269, 313, 375]
        dmid_expected = [0.197, 0.224, 0.261, 0.311]

        # Check that results are within 1% of expected results
        for i, val in enumerate(Mmid_calculated):
            self.assertAlmostEqual(Mmid_calculated[i]/Mmid_expected[i], 1, 2)
            self.assertAlmostEqual(dmid_calculated[i]/dmid_expected[i], 1, 2)
    
    def test_AISC_benchmark_case2(self):

        column = FEModel3D()

        column.add_node('N1', 0, 0, 0)
        column.add_node('N2', 0, 28, 0)
        
        # Add an internal node to capture P-little-delta effects
        column.add_node('N4', 0, 28*1/3, 0)
        column.add_node('N5', 0, 28*2/3, 0)

        column.def_support('N1', True, True, True, True, True, True)

        # Define a material
        E = 29000*12**2  # ksf
        G = 11200*12**2  # ksf
        nu = 0.3
        rho = 0.490  # kcf
        column.add_material('Steel', E, G, nu, rho)

        # Define section properties
        Iy = 51.4/12**4  # ft^4
        Iz = 484/12**4   # ft^4
        A = 14.1/12**2   # ft^2
        J = 1.45/12**4   # ft^4
        column.add_section('ColumnSection', A, Iy, Iz, J)

        column.add_member('M1', 'N1', 'N2', 'Steel', 'ColumnSection')
        
        column.add_node_load('N2', 'FX', 1, case='P1')
        column.add_node_load('N2', 'FX', 1, case='P2')
        column.add_node_load('N2', 'FX', 1, case='P3')
        column.add_node_load('N2', 'FX', 1, case='P4')

        column.add_node_load('N2', 'FY', -100, 'P2')
        column.add_node_load('N2', 'FY', -150, 'P3')
        column.add_node_load('N2', 'FY', -200, 'P4')

        column.add_load_combo('Combo 1', {'P1':1})
        column.add_load_combo('Combo 2', {'P2':1})
        column.add_load_combo('Combo 3', {'P3':1})
        column.add_load_combo('Combo 4', {'P4':1})

        # from PyNite.Visualization import Renderer
        # renderer = Renderer(column)
        # renderer.annotation_size = 1
        # renderer.combo_name = 'Combo 2'
        # renderer.render_model()

        column.analyze_PDelta()

        Mbase_calculated = []
        dtip_calculated = []
        for combo in column.load_combos.values():
            Mbase_calculated.append(-column.members['M1'].moment('Mz', 0, combo.name)*12)
            dtip_calculated.append(-column.members['M1'].deflection('dy', 28, combo.name)*12)

        # Expected results per AISC
        Mbase_expected = [-336, -469, -598, -848]
        dtip_expected = [0.901, 1.33, 1.75, 2.56]

        # Check that results are within 1% of expected results
        for i, val in enumerate(Mbase_calculated):
            self.assertLessEqual(abs(Mbase_calculated[i]/Mbase_expected[i] - 1), 0.03)
            self.assertLessEqual(abs(dtip_calculated[i]/dtip_expected[i] - 1), 0.05)
