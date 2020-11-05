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
    ''' Test torsion loads. '''

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

    def test_member_torque_load(self):
        """
        Determine the member of the member.

        Args:
            self: (todo): write your description
        """
        TorqueBeam = FEModel3D()
        # Add nodes (14 ft = 168 in apart)
        TorqueBeam.AddNode('N1', 0, 0, 0)
        TorqueBeam.AddNode('N2', 168, 0, 0)
        # Add a beam with the following properties:
        TorqueBeam.AddMember('M1', 'N1', 'N2', 29000, 11400, 100, 150, 250, 20)
        # Provide fixed supports
        TorqueBeam.DefineSupport('N1', False, True, True, True, True, True)
        TorqueBeam.DefineSupport('N2', True, True, True, True, True, True)
        # Add a point load of 5 kip-ft and 10 kip-ft at 3ft and 11 ft along the beam respectively
        TorqueBeam.AddMemberPtLoad('M1', 'Mx', 5, 3*12)
        TorqueBeam.AddMemberPtLoad('M1', 'Mx', 10, 11*12)
        TorqueBeam.Analyze(check_statics=True)
        # Support reactions (around local x axis) at each end
        left_Rxn = TorqueBeam.GetNode('N1').RxnMX['Combo 1']
        # subTest context manager prints which portion fails, if any
        with self.subTest(left_Rxn=left_Rxn):
            self.assertAlmostEqual(left_Rxn, -6.07, 2)
        
        right_Rxn = TorqueBeam.GetNode('N2').RxnMX['Combo 1']
        with self.subTest(right_Rxn=right_Rxn):
            self.assertAlmostEqual(right_Rxn, -8.93, 2)
        
        # Max/min torques on the beam
        max_torque = TorqueBeam.GetMember('M1').MaxTorsion()
        with self.subTest(max_torque=max_torque):
            self.assertAlmostEqual(max_torque, 8.93, 2)
        
        min_torque = TorqueBeam.GetMember('M1').MinTorsion()
        with self.subTest(min_torque=min_torque):
            self.assertAlmostEqual(min_torque, -6.07, 2)
