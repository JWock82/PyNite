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

    def test_sloped_beam(self):
        """
        Units for this model are kips and feet
        """

        # Define a new beam
        beam = FEModel3D()

        # Define the nodes
        beam.add_node('N1', 0, 0, 0)
        beam.add_node('N2', 10, 10, 0)

        # Define the supports
        beam.def_support('N1', True, True, True, True, True, True)
        beam.def_support('N2', True, True, True, True, True, True)

        # Define a material
        beam.add_material('Steel', 29000*144, 11200*144, 0.3, 490/1000)

        # Define beam section proerties
        J = 400/12**4
        Iy = 200/12**4
        Iz = 200/12**4
        A = 12/12**2
        beam.add_section('Section', A, Iy, Iz, J)

        # Create the beam
        beam.add_member('M1', 'N1', 'N2', 'Steel', 'Section')
        
        # Hinge the ends of the beam
        beam.def_releases('M1', False, False, False, False, True, True, False, False, False, False, True, True)

        # Add a member point load
        L = (10**2 + 10**2)**0.5
        beam.add_member_pt_load('M1', 'FY', -5, L/3, 'D')
        beam.add_member_pt_load('M1', 'FY', -5, 2*L/3, 'D')
        beam.add_load_combo('1.4D', {'D':1.4})

        # Analyze the model
        beam.analyze_PDelta()

        # from PyNite.Visualization import Renderer
        # renderer = Renderer(beam)
        # renderer.combo_name = '1.4D'
        # renderer.annotation_size = 1
        # renderer.render_model()

        RY1 = beam.nodes['N1'].RxnFY['1.4D']
        RY2 = beam.nodes['N2'].RxnFY['1.4D']

        # Check the result
        self.assertAlmostEqual(RY1/(1.4*5), 1, 2)
        self.assertAlmostEqual(RY2/(1.4*5), 1, 2)
