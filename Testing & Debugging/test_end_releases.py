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

    def test_end_release_Rz(self):
        """
        R calculate the end of model

        Args:
            self: (todo): write your description
        """
        myModel = FEModel3D()
        # Add two supported nodes and one member
        myModel.AddNode('N1', 0, 0, 0)
        myModel.AddNode('N2', 10*12, 0, 0)
        myModel.DefineSupport('N1', True, True, True, True, True, True)
        myModel.DefineSupport('N2', True, True, True, True, True, True)
        myModel.AddMember('M1', 'N1', 'N2', 29000, 11400, 100, 150, 250, 10)
        # Release Rzi and Rzj on member M1
        myModel.DefineReleases('M1', False, False, False, False, False, True, \
                                    False, False, False, False, False, True)
        # Add a load
        myModel.AddMemberDistLoad('M1', 'Fy', -0.5, -0.5)
        myModel.Analyze()
        # Get the resulting moments
        calculated_moment = myModel.GetMember('M1').MinMoment('Mz')
        expected_moment = -0.5*(10*12)**2/8
        self.assertAlmostEqual(calculated_moment, expected_moment)
