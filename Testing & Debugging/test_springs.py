# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
import sys
from io import StringIO

class Test_Spring_Elements(unittest.TestCase):
    ''' Tests of spring members.'''

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
        
    def test_spring_elements(self): 
        """
        Method to true if a list of cluster elements.

        Args:
            self: (todo): write your description
        """
        # A First Course in the Finite Element Method, 4th Edition
        # Daryl L. Logan
        # Example 2.1
        # Units for this model are pounds and inches
        system = FEModel3D()
        system.AddNode('1', 0, 0, 0)
        system.AddNode('2', 30, 0, 0)
        system.AddNode('3', 10, 0, 0)
        system.AddNode('4', 20, 0, 0)
        # Add spring members
        system.AddSpring('S1', '1', '3', 1000)
        system.AddSpring('S2', '3', '4', 2000)
        system.AddSpring('S3', '4', '2', 3000)
        # Define supports
        system.DefineSupport('1', True, True, True, True, True, True)
        system.DefineSupport('2', True, True, True, True, True, True)
        system.DefineSupport('3', False, True, True, True, True, True)
        system.DefineSupport('4', False, True, True, True, True, True)
        # Add node loads
        system.AddNodeLoad('4', 'FX', 5000)
        system.Analyze(True)
        # Check results
        # correct_values = [('3', 0.9090909090909092),
        #                   ('4', 1.3636363636363638),
        #                   ('1', -909.0909090909091),
        #                   ('2', -4090.9090909090914)]
        n3_DX = system.GetNode('3').DX['Combo 1']
        self.assertAlmostEqual(n3_DX/ 0.9090909090909092, 1.0, 2)

        n4_DX = system.GetNode('4').DX['Combo 1']
        self.assertAlmostEqual(n4_DX/1.3636363636363638, 1.0, 2)
        
        n1_rxn = system.GetNode('1').RxnFX['Combo 1']
        self.assertAlmostEqual(n1_rxn/-909.0909090909091, 1.0, 2)
        
        n2_rxn = system.GetNode('2').RxnFX['Combo 1']
        self.assertAlmostEqual(n2_rxn/-4090.9090909090914, 1.0, 2)
        