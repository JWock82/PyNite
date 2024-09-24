# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest

from numpy import True_
from PyNite import FEModel3D
import sys
from io import StringIO

class Test_Unstable(unittest.TestCase):
    ''' Tests that should raise instability errors.'''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__
        
    def test_unstable_supports(self):

        # This test checks the PyNite's ability to detect unstable support conditions
        # Units used in this test are inches, and kips
        MomentFrame = FEModel3D()

        # Add nodes (frame is 15 ft wide x 12 ft tall)
        MomentFrame.add_node("N1", 0, 0, 0)
        MomentFrame.add_node("N2", 0, 12*12, 0)
        MomentFrame.add_node("N3", 15*12, 12*12, 0)
        MomentFrame.add_node("N4", 15*12, 0*12, 0)

        # Add a material
        MomentFrame.add_material('Steel', 29000, 11400, 0.5, 490/1000/12**3)

        # Add columns with the following properties:
        # Iy = 100 in^4, Iz = 150 in^4, J = 250 in^4, A = 10 in^2
        MomentFrame.add_section('Section', 10, 100, 150, 250)
        MomentFrame.add_member("M1", "N1", "N2", 'Steel', 'Section')
        MomentFrame.add_member("M2", "N4", "N3", 'Steel', 'Section')

        # Add a beam with the following properties:
        # Iy = 100 in^4, Iz = 250 in^4, J = 250 in^4, A = 15 in^2
        MomentFrame.add_section('Section2', 15, 100, 250, 250)
        MomentFrame.add_member("M3", "N2", "N3", 'Steel', 'Section2')

        # Pin the ends of the columns
        MomentFrame.def_releases('M1', Dzi=True)
        MomentFrame.def_releases('M2', Dzi=True)

        # Provide unstable supports (unsupported in DZ)
        MomentFrame.def_support("N1", False, False, False, False, False, True)
        MomentFrame.def_support("N4", False, False, False, False, False, True)

        # Add a nodal lateral load of 50 kips at the left side of the frame
        MomentFrame.add_node_load("N2", "FX", 50)
        # Analyze the frame - we should see an error message that the structure is unstable
        with self.assertRaises(Exception):
            MomentFrame.analyze()