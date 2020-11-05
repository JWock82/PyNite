# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1

From:
Structural Analysis, 3rd Edition
Aslam Kassimali
Example 13.14
"""

import unittest
from PyNite import FEModel3D
import math
import sys
from io import StringIO

class Test_Support_Settlement(unittest.TestCase):
    ''' Test for support settlements. '''

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

    def test_support_settlement(self):
        """
        Test if all the available to see https : // www. python 3.

        Args:
            self: (todo): write your description
        """
        beam = FEModel3D()
        # Add nodes
        beam.AddNode('A', 0, 0, 0)
        beam.AddNode('B', 20*12, 0, 0)
        beam.AddNode('C', 40*12, 0, 0)
        beam.AddNode('D', 60*12, 0, 0)
        # Add members
        A = 20
        E = 29000
        G = 11400
        Iy = 1000
        Iz = 7800
        J = 8800
        beam.AddMember('AB', 'A', 'B', E, G, Iy, Iz, J, A)
        beam.AddMember('BC', 'B', 'C', E, G, Iy, Iz, J, A)
        beam.AddMember('CD', 'C', 'D', E, G, Iy, Iz, J, A)
        # Provide supports
        beam.DefineSupport('A', True, True, True, True, False, False)
        beam.DefineSupport('B', False, True, True, False, False, False)
        beam.DefineSupport('C', False, True, True, False, False, False)
        beam.DefineSupport('D', False, True, True, False, False, False)
        # Add a uniform load to the beam
        beam.AddMemberDistLoad('AB', 'Fy', -2/12, -2/12)
        beam.AddMemberDistLoad('BC', 'Fy', -2/12, -2/12)
        beam.AddMemberDistLoad('CD', 'Fy', -2/12, -2/12)
        # Add support settlements
        beam.AddNodeDisplacement('B', 'DY', -5/8)
        beam.AddNodeDisplacement('C', 'DY', -1.5)
        beam.AddNodeDisplacement('D', 'DY', -0.75)
        # Analyze the beam
        beam.Analyze()
        # subTest context manager prints which portion fails, if any
        correct_values = [('A', -1.098),
                          ('B',  122.373),
                          ('C', -61.451),
                          ('D',  60.176)]
        for name, value in correct_values:
            with self.subTest(node=name):
                calculated_Rxn = beam.GetNode(name).RxnFY['Combo 1']
                # Two decimal place accuracy requires +/-0.5% accuracy
                # one decimal place requires +/-5%
                self.assertAlmostEqual(calculated_Rxn/value, 1.0, 2)
