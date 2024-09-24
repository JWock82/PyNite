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
    """Tests for tension/compression-only analysis"""

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_TC_members(self):

        # Create a new finite element model
        tc_model = FEModel3D()
        tc_model.add_node('N1', 0, 0, 0)
        tc_model.add_node('N2', 100, 0, 0)
        tc_model.add_node('N3', 0, 10, 0)
        tc_model.add_node('N4', 0, -10, 0)

        E = 29000 # ksi
        G = 11400 # ksi
        nu = 0.3  # Poisson's ratio
        rho = 0.490/12**2  # Density (kci)
        tc_model.add_material('Steel', E, G, nu, rho)

        Iy = 3 # in^4
        Iz = 3 # in^4
        J = 0.0438 # in^4
        A = 1.94 # in^2
        tc_model.add_section('Section', A, Iy, Iz, J)

        tc_model.add_member('both-ways', 'N1', 'N2', 'Steel', 'Section')
        tc_model.add_member('t-only top', 'N3', 'N2', 'Steel', 'Section', tension_only=True)
        tc_model.def_releases('t-only top', Ryi=True, Rzi=True, Ryj=True, Rzj=True)
        tc_model.def_releases('both-ways', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

        tc_model.def_support('N2', *[False]*2, *[True]*4)
        tc_model.def_support('N1', *[True]*6)
        tc_model.def_support('N3', *[True]*6)
        tc_model.def_support('N4', *[True]*6)
        tc_model.add_node_load('N2', 'FY', -10)

        tc_model.add_member('t-only bott', 'N4', 'N2', 'Steel', 'Section', tension_only=True)
        tc_model.def_releases('t-only bott', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

        tc_model.analyze()

        self.assertAlmostEqual(tc_model.members['t-only top'].max_axial(), -100.499, 3)
        self.assertAlmostEqual(tc_model.members['both-ways'].max_axial(), 100, 3)
        self.assertEqual(tc_model.members['t-only bott'].max_axial(), 0, 3)
        self.assertFalse(tc_model.members['t-only bott'].active['Combo 1'])
