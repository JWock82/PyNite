# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
import sys
from io import StringIO

class Test_Unstable(unittest.TestCase):
    ''' Tests that should raise instability errors.'''

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
        
    def test_unstable_supports(self):    
        """
        Test if the test is unstable test.

        Args:
            self: (todo): write your description
        """
        # This test checks the PyNite's ability to detect unstable support conditions
        # Units used in this test are inches, and kips
        MomentFrame = FEModel3D()
        # Add nodes (frame is 15 ft wide x 12 ft tall)
        MomentFrame.AddNode("N1", 0, 0, 0)
        MomentFrame.AddNode("N2", 0, 12*12, 0)
        MomentFrame.AddNode("N3", 15*12, 12*12, 0)
        MomentFrame.AddNode("N4", 15*12, 0*12, 0)
        # Add columns with the following properties:
        # E = 29000 ksi, G = 11400 ksi, Iy = 100 in^4, Iz = 150 in^4, J = 250 in^4, A = 10 in^2
        MomentFrame.AddMember("M1", "N1", "N2", 29000, 11400, 100, 150, 250, 10)
        MomentFrame.AddMember("M2", "N4", "N3", 29000, 11400, 100, 150, 250, 10)
        # Add a beam with the following properties:
        # E = 29000 ksi, G = 11400 ksi, Iy = 100 in^4, Iz = 250 in^4, J = 250 in^4, A = 15 in^2
        MomentFrame.AddMember("M3", "N2", "N3", 29000, 11400, 100, 250, 250, 15)
        # Provide unstable supports (unsupported in DY)
        MomentFrame.DefineSupport("N1", True, False, True, True, True, True)
        # Add a nodal lateral load of 50 kips at the left side of the frame
        MomentFrame.AddNodeLoad("N2", "FX", 50)
        # Analyze the frame - we should see an error message that the structure is unstable
        with self.assertRaises(Exception):
            MomentFrame.Analyze()