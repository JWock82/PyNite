# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
import sys
from io import StringIO

class Test_AxialLoads(unittest.TestCase):
    ''' Tests of member axial loads.'''

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
        
    def test_axial_distributed_load(self):    
        """
        This function to load of the axes of - axes.

        Args:
            self: (todo): write your description
        """
        # Units N e m
        Beam = FEModel3D()
        L = 5 # m
        # Nodes
        Beam.AddNode("N1", 0, 0, 0)
        Beam.AddNode("N2", L, 0, 0)
        # Beams (30x50 cm)
        E = 2.1e11 # N/m^2
        G = 1
        Iy = 0.001125 # m^4
        Iz = 0.003125 # m^4
        J = 1
        A = 0.15 # m^2
        Beam.AddMember("M1", "N1", "N2", E, G, Iy, Iz, J, A)
        # Supports
        Beam.DefineSupport("N1", True, True, True, True, True, True)
        Beam.DefineSupport("N2", True, True, True, False, True, True)
        # Load
        Beam.AddMemberDistLoad("M1", "Fx", 10, 10, 0, 5)
        # Analyze
        Beam.Analyze()
        # Member fixed end reaction vector
        # print('M1 Displacement Vector: ', Beam.GetMember('M1').d())
        # print('M1 Fixed End Reaction Vector: ', Beam.GetMember('M1').fer())
        # Reactions
        for node_name in ('N1', 'N2'):
            with self.subTest(node=node_name):
                rxn = Beam.GetNode(node_name).RxnFX['Combo 1']
                self.assertAlmostEqual(rxn/-25.0, 1.0, 2)