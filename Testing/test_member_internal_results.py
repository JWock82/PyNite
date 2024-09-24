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

class TestMemberInternalResults(unittest.TestCase):
    """
    Tests for member internal results
    """

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__
    
    def test_beam_internal_forces(self):
        """
        Units for this model are kips and feet
        """

        # Define a new beam
        beam = FEModel3D()

        # Define the nodes
        beam.add_node('N1', 0, 0, 0)
        beam.add_node('N2', 10, 0, 0)

        # Define the supports
        beam.def_support('N1', True, True, True, True, False, False)
        beam.def_support('N2', True, True, True, False, False, False)

        # Define beam section proerties
        J = 400/12**4
        Iy = 200/12**4
        Iz = 200/12**4
        A = 12/12**2
        beam.add_section('Section', A, Iy, Iz, J)

        # Define a material
        E = 29000*144  # ksf
        G = 11200*144  # ksf
        nu = 0.3
        rho = 0.490  # pcf
        beam.add_material('Steel', E, G, nu, rho)

        # Create the beam
        beam.add_member('M1', 'N1', 'N2', 'Steel', 'Section')

        # Add a mid-span node
        beam.add_node('N3', 5, 0, 0)

        # Add a member distributed load along the strong axis
        beam.add_member_dist_load('M1', 'FY', -0.5, -0.5, case='D')
        beam.add_load_combo('D', {'D':1.0})

        # Add a member distributed laod along the weak axis
        beam.add_member_dist_load('M1', 'FZ', -0.5, -0.5, case='D')
        beam.add_load_combo('D', {'D':1.0})

        # Analyze the model
        beam.analyze_linear()

        # from PyNite.Visualization import Renderer
        # renderer = Renderer(beam)
        # renderer.combo_name = 'D'
        # renderer.annotation_size = 1
        # renderer.render_model()

        # Check the shear diagram
        self.assertAlmostEqual(beam.members['M1'].shear('Fy', 0, 'D'), 2.5, 2)
        self.assertAlmostEqual(beam.members['M1'].shear('Fz', 0, 'D'), 2.5, 2)
        self.assertAlmostEqual(beam.members['M1'].shear('Fy', 10, 'D'), -2.5, 2)
        self.assertAlmostEqual(beam.members['M1'].shear('Fz', 10, 'D'), -2.5, 2)
        self.assertAlmostEqual(beam.members['M1'].shear('Fy', 5, 'D'), 0, 2)
        self.assertAlmostEqual(beam.members['M1'].shear('Fz', 5, 'D'), 0, 2)

        # Check the moment diagram
        self.assertAlmostEqual(beam.members['M1'].moment('Mz', 0, 'D'), 0, 2)
        self.assertAlmostEqual(beam.members['M1'].moment('My', 0, 'D'), 0, 2)
        self.assertAlmostEqual(beam.members['M1'].moment('Mz', 5, 'D'), -6.25, 2)
        self.assertAlmostEqual(beam.members['M1'].moment('My', 5, 'D'), -6.25, 2)
        self.assertAlmostEqual(beam.members['M1'].moment('Mz', 10, 'D'), 0, 2)
        self.assertAlmostEqual(beam.members['M1'].moment('My', 10, 'D'), 0, 2)

        # Check the deflected shape
        self.assertAlmostEqual(beam.members['M1'].deflection('dy', 0, 'D')*12, 0, 5)
        self.assertAlmostEqual(beam.members['M1'].deflection('dz', 0, 'D')*12, 0, 5)
        self.assertAlmostEqual(beam.members['M1'].deflection('dz', 5, 'D')*12, 5*(-0.5)*10**4/(384*E*Iz)*12, 5)
        self.assertAlmostEqual(beam.members['M1'].deflection('dy', 5, 'D')*12, 5*(-0.5)*10**4/(384*E*Iy)*12, 5)
        self.assertAlmostEqual(beam.members['M1'].deflection('dy', 10, 'D')*12, 0, 5)
        self.assertAlmostEqual(beam.members['M1'].deflection('dz', 10, 'D')*12, 0, 5)
