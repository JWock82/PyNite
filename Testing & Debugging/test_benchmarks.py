# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
import math
import sys
from io import StringIO

class Test_AISC_Benchmark(unittest.TestCase):
    ''' AISC's benchmark problem used to determine if a second order analysis 
    procedure is rigorous enough.'''

    def setUp(self):
        """
        Sets the program to stdout.

        Args:
            self: (todo): write your description
        """
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        """
        Tear down the terminal.

        Args:
            self: (todo): write your description
        """
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_AISC_Benckmark(self):
        """
        Calculate aISCckmarkmarkmark.

        Args:
            self: (todo): write your description
        """
        # Create the cantilever model
        cantilever = FEModel3D()

        # Define the column and its properties
        L = 20 # ft
        H = 5 # kips lateral load
        P = 100 # kips axial load
        G = 11200*12**2 # shear modulus (ksf)
        E = 29000*12**2 # modulus of elasticity (ksf)
        I = 100/12**4 # moment of inertia (ft^4)

        # Break the column into several segments in order to capture P-little-delta effects
        num_segs = 5
        num_nodes = num_segs + 1
        for i in range(num_nodes):
            # Add nodes
            cantilever.AddNode(str(i+1), 0, i*L/(num_segs), 0)

        for i in range(num_segs):
            # Add members between nodes
            cantilever.AddMember(str(i+1), str(i+1), str(i+2), E, G, I, I, 200/12**4, 10/12**2)

        # Add a fixed support at the base of the column
        cantilever.DefineSupport('1', True, True, True, True, True, True)

        # Add a -10 kip axial load to the top of the column
        cantilever.AddNodeLoad(str(num_nodes), 'FY', -P)

        # Add a 5 kip lateral load to the top of the column
        cantilever.AddNodeLoad(str(num_nodes), 'FX', H)

        # Perform 2nd order analysis
        cantilever.Analyze_PDelta()

        # The moment at the base of the column
        calculated_moment = cantilever.GetNode('1').RxnMZ['Combo 1']

        # the deflection at the top of the column
        calculated_displacement = cantilever.GetNode(str(num_nodes)).DX['Combo 1']*12

        # Calculate the AISC benchmark problem solution:
        alpha = (P*L**2/(E*I))**0.5
        Mmax = H*L*(math.tan(alpha)/alpha)
        ymax = H*L**3/(3*E*I)*(3*(math.tan(alpha)-alpha)/alpha**3)
        
        # Compare the calculation results
        self.assertAlmostEqual(calculated_moment/Mmax, 1.0, 1)
        self.assertAlmostEqual(calculated_displacement/(ymax*12), 1.0, 1)
