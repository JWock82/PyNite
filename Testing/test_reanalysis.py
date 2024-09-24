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

class Test_2D_Frame(unittest.TestCase):
    ''' Tests of analyzing 2D frames. '''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_reanalysis(self):

        beam = FEModel3D()

        # Define the nodes
        beam.add_node('N1', 0, 0, 0)
        beam.add_node('N2', 30*12, 0, 0)
        beam.add_node('N3', 60*12, 0, 0)

        # Define the supports
        beam.def_support('N1', True, True, True, True, False, False)
        beam.def_support('N2', True, True, True, False, False, False)
        beam.def_support('N3', True, True, True, False, False, False)

        # Create members (all members will have the same properties in this example)
        Iz = 24*30**2/12
        Iy = 30*24**3/12
        J = Iz + Iy
        A = 24*30
        beam.add_section('Section', A, Iy, Iz, J)

        # Define a material
        E = 57*(4500)**0.5
        G = 0.4*E
        beam.add_material('Concrete', E, G, 0.17, 0.150/12**2)
        
        beam.add_member('M1', 'N1', 'N3', 'Concrete', 'Section')

        # Add a member load
        w = -0.100*(10*12)
        beam.add_member_dist_load('M1', 'FY', w, w)

        # Analyze the model
        beam.analyze()

        d1 = beam.members['M1'].min_deflection('dy')

        # Change the moment of inertia to account for cracking
        Iz = 0.35*Iz
        beam.members['M1'].section.Iz = Iz

        beam.analyze()

        d2 = beam.members['M1'].min_deflection('dy')

        # print(d1, d2)

        self.assertNotAlmostEqual(d1/d2, 1, 3)
