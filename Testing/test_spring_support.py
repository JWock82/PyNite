# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2021 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
import sys
from io import StringIO

class Test_Spring_Supports(unittest.TestCase):
    ''' Tests of spring supports.'''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__
        
    def test_beam_on_elastic_foundation(self):
        '''        
        Matrix Structural Analysis, 2nd Edition
        William McGuire, Richard H. Gallagher, Ronald D. Ziemian
        Example 4.15
        Units for this model are kips and inches
        '''

        # Create a new model
        boef = FEModel3D()

        # Define nodes
        for i in range(17):

            # Add nodes spaced at 15"
            boef.add_node('N' + str(i + 1), i*15, 0, 0)

            # Add supports to the nodes
            if i == 0 or i == 16:
                boef.def_support('N' + str(i + 1), True, True, True, True, False, False)
            else:
                boef.def_support_spring('N' + str(i + 1), 'DY', 22.5, '-')

        # Define a material
        boef.add_material('Steel', 29000, 11200, 0.3, 490/1000/12**3)

        # Define section properties
        A = 10.3    # in^2
        Iz = 128.5  # in^4 (strong axis)
        Iy = 42.6   # in^4 (weak axis)
        J = 0.769   # in^4
        boef.add_section('Section', A, Iy, Iz, J)

        # Define members
        for i in range(16):

            # Add the members
            boef.add_member('M' + str(i + 1), 'N' + str(i + 1), 'N' + str(i + 2), 'Steel', 'Section')
        
        # Add a point load at midspan
        boef.add_node_load('N9', 'FY', -40)

        # Analyze the model
        boef.analyze(sparse=False)

        print(boef.members['M8'].min_moment('Mz'))
        print(boef.members['M8'].max_moment('Mz'))

        # Check that results are within 5% of the expected answer
        self.assertLess(boef.nodes['N9'].DY['Combo 1']/(-0.238) - 1, 0.05, 'Failed beam on elastic foundation test.')
        self.assertLess(-boef.members['M8'].min_moment('Mz')/547 - 1, 0.05, 'Failed beam on elastic foundation test.')