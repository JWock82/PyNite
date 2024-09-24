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
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_member_torque_load(self):
        TorqueBeam = FEModel3D()
        # Add nodes (14 ft = 168 in apart)
        TorqueBeam.add_node('N1', 0, 0, 0)
        TorqueBeam.add_node('N2', 168, 0, 0)
        # Add a material
        TorqueBeam.add_material('Steel', 29000, 11400, 0.5, 490/1000/12**3)
        #Add section
        TorqueBeam.add_section('Section', 20, 100, 150, 250)
        # Add a beam with the following properties:
        TorqueBeam.add_member('M1', 'N1', 'N2', 'Steel', 'Section')
        # Provide fixed supports
        TorqueBeam.def_support('N1', False, True, True, True, True, True)
        TorqueBeam.def_support('N2', True, True, True, True, True, True)
        # Add a point load of 5 kip-ft and 10 kip-ft at 3ft and 11 ft along the beam respectively
        TorqueBeam.add_member_pt_load('M1', 'Mx', 5, 3*12)
        TorqueBeam.add_member_pt_load('M1', 'Mx', 10, 11*12)
        TorqueBeam.analyze(check_statics=True)
        # Support reactions (around local x axis) at each end
        left_Rxn = TorqueBeam.nodes['N1'].RxnMX['Combo 1']
        # subTest context manager prints which portion fails, if any
        with self.subTest(left_Rxn=left_Rxn):
            self.assertAlmostEqual(left_Rxn, -6.07, 2)
        
        right_Rxn = TorqueBeam.nodes['N2'].RxnMX['Combo 1']
        with self.subTest(right_Rxn=right_Rxn):
            self.assertAlmostEqual(right_Rxn, -8.93, 2)
        
        # Max/min torques on the beam
        max_torque = TorqueBeam.members['M1'].max_torque()
        with self.subTest(max_torque=max_torque):
            self.assertAlmostEqual(max_torque, 8.93, 2)
        
        min_torque = TorqueBeam.members['M1'].min_torque()
        with self.subTest(min_torque=min_torque):
            self.assertAlmostEqual(min_torque, -6.07, 2)
