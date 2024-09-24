# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D, Section
import sys
from io import StringIO

class TestLoads(unittest.TestCase):
    ''' Tests of member axial loads.'''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__
    
    def test_member_self_weight(self):

        cont_beam = FEModel3D()
        cont_beam.add_node('N1', 0, 0, 0)
        cont_beam.add_node('N2', 15, 0, 0)
        cont_beam.add_node('N3', 30, 0, 0)
        cont_beam.add_material('Steel', 29000*144, 11200*144, 0.3, 0.490)
        cont_beam.add_section('Section', 10/12**2, 23.3/12**4, 340/12**4, 0.569/12**4)
        cont_beam.add_member('M1', 'N1', 'N3', 'Steel', 'Section')  # W14x34
        cont_beam.def_support('N1', True, True, True, True, False, False)
        cont_beam.def_support('N2', False, True, True, False, False, False)
        cont_beam.def_support('N3', False, True, True, False, False, False)
        cont_beam.add_member_self_weight('FY', 1, 'D')
        cont_beam.add_load_combo('D', {'D':1.0})
        cont_beam.analyze()
        for member in cont_beam.members.values():
            self.assertAlmostEqual(member.DistLoads[0][1], 0.034, 4)
            self.assertAlmostEqual(member.DistLoads[0][2], 0.034, 4)

    def test_axial_distributed_load(self):

        # Units N e m
        Beam = FEModel3D()
        L = 5 # m

        # Nodes
        Beam.add_node("N1", 0, 0, 0)
        Beam.add_node("N2", L, 0, 0)

        # Define a material
        E = 2.1e11 # N/m^2
        G = 1
        nu = 0.3
        rho = 1
        Beam.add_material('Mat1', E, G, nu, rho)

        # Beams (30x50 cm)
        Iy = 0.001125 # m^4
        Iz = 0.003125 # m^4
        J = 1
        A = 0.15 # m^2
        Beam.add_section('Section', A, Iy, Iz, J)
        Beam.add_member("M1", "N1", "N2", 'Mat1', 'Section')

        # Supports
        Beam.def_support("N1", True, True, True, True, True, True)
        Beam.def_support("N2", True, True, True, False, True, True)

        # Load
        Beam.add_member_dist_load("M1", "Fx", 10, 10, 0, 5)

        # Analyze
        Beam.analyze()
        
        # Member fixed end reaction vector
        # print('M1 Displacement Vector: ', Beam.members['M1'].d())
        # print('M1 Fixed End Reaction Vector: ', Beam.members['M1'].fer())
        # Reactions
        for node_name in ('N1', 'N2'):
            with self.subTest(node=node_name):
                rxn = Beam.nodes[node_name].RxnFX['Combo 1']
                self.assertAlmostEqual(rxn/-25.0, 1.0, 2)