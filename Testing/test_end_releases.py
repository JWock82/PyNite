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

class Test_End_Release(unittest.TestCase):
    ''' Test member end releases. '''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_end_release_Rz(self):
        
        # Create the model
        myModel = FEModel3D()

        # Define a material
        E = 29000  # ksi
        G = 11200  # ksi
        nu = 0.3
        rho = 0.490/12**3  # kci
        myModel.add_material('Steel', E, G, nu, rho)
        myModel.add_section('ColumnSection', 10, 100, 150, 250)

        # Add two supported nodes and one member
        myModel.add_node('N1', 0, 0, 0)
        myModel.add_node('N2', 10*12, 0, 0)
        myModel.def_support('N1', True, True, True, True, True, True)
        myModel.def_support('N2', True, True, True, True, True, True)
        myModel.add_member('M1', 'N1', 'N2', 'Steel', 'ColumnSection')

        # Release Rzi and Rzj on member M1
        myModel.def_releases('M1', False, False, False, False, False, True, \
                                    False, False, False, False, False, True)
        
        # Add a load
        myModel.add_member_dist_load('M1', 'Fy', -0.5, -0.5)
        myModel.analyze()

        # Get the resulting moments
        calculated_moment = myModel.members['M1'].min_moment('Mz')
        expected_moment = -0.5*(10*12)**2/8
        
        self.assertAlmostEqual(calculated_moment, expected_moment)
