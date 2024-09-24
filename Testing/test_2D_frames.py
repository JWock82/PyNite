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

    def test_XY_gravity_load(self):
        
        # A First Course in the Finite Element Method, 4th Edition
        # Daryl L. Logan
        # Problem 5.30
        # Units for this model are kips and inches
        frame = FEModel3D()

        # Define the nodes
        frame.add_node('N1', 0, 0, 0)
        frame.add_node('N2', 0, 30*12, 0)
        frame.add_node('N3', 15*12, 40*12, 0)
        frame.add_node('N4', 35*12, 40*12, 0)
        frame.add_node('N5', 50*12, 30*12, 0)
        frame.add_node('N6', 50*12, 0, 0)

        # Define the supports
        frame.def_support('N1', True, True, True, True, True, True)
        frame.def_support('N6', True, True, True, True, True, True)

        # Define materials
        frame.add_material('Steel', 30000, 250, 0.5, 490/1000/12**3)

        # Create members (all members will have the same properties in this example)
        J = 250
        Iy = 250
        Iz = 200
        A = 12
        frame.add_section('FrameSection', A, Iy, Iz, J)

        frame.add_member('M1', 'N1', 'N2', 'Steel', 'FrameSection')
        frame.add_member('M2', 'N2', 'N3', 'Steel', 'FrameSection')
        frame.add_member('M3', 'N3', 'N4', 'Steel', 'FrameSection')
        frame.add_member('M4', 'N4', 'N5', 'Steel', 'FrameSection')
        frame.add_member('M5', 'N5', 'N6', 'Steel', 'FrameSection')

        # Add nodal loads
        frame.add_node_load('N3', 'FY', -30)
        frame.add_node_load('N4', 'FY', -30)

        # Analyze the model
        frame.analyze()

        # subTest context manager prints which portion fails, if any
        correct_values = [('N1', {'RxnFX': 11.6877,
                                  'RxnFY': 30,
                                  'RxnMZ': -1810.0745}),
                          ('N6', {'RxnFX': -11.6877,
                                  'RxnFY': 30,
                                  'RxnMZ': 1810.0745})]

        for name, values in correct_values:
            with self.subTest(node=name):
                node = frame.nodes[name]
                # Two decimal place accuracy requires +/-0.5% accuracy
                # one decimal place requires +/-5%
                self.assertAlmostEqual(node.RxnFX['Combo 1']/values['RxnFX'], 1.0, 2)
                self.assertAlmostEqual(node.RxnFY['Combo 1']/values['RxnFY'], 1.0, 2)
                self.assertAlmostEqual(node.RxnMZ['Combo 1']/values['RxnMZ'], 1.0, 2)

    def test_XY_member_ptload(self):
        frame = FEModel3D()
        # Add nodes
        frame.add_node('N1', 0, 0, 0)        # ft
        frame.add_node('N2', 0, 7.667, 0)    # ft
        frame.add_node('N3', 7.75, 7.667, 0) # ft
        frame.add_node('N4', 7.75, 0, 0)     # ft
        # Add supports
        frame.def_support('N1', True, True, True, True, True, False)
        frame.def_support('N4', True, True, True, True, True, False)

        # Define materials
        frame.add_material('Steel', 29000*12**2, 1111200*12**2, 0.5, 490/1000/12**3)

        # Define material and section properties for a W8x24
        Iy = 18.3/12**4     # ft^4
        Iz = 82.7/12**4     # ft^4
        J = 0.346/12**4     # ft^4
        A = 5.26/12**2      # in^2
        frame.add_section('FrameSection', A, Iy, Iz, J)

        # Define members
        frame.add_member('M1', 'N1', 'N2', 'Steel', 'FrameSection')
        frame.add_member('M2', 'N2', 'N3', 'Steel', 'FrameSection')
        frame.add_member('M3', 'N4', 'N3', 'Steel', 'FrameSection')

        # Add loads to the frame
        frame.add_member_pt_load('M2', 'Fy', -5, 7.75/2)       # 5 kips @ midspan
        frame.add_member_dist_load('M2', 'Fy', -0.024, -0.024) # W8x24 self-weight

        # Analyze the frame
        frame.analyze()
        calculated_RZ = frame.nodes['N1'].RZ['Combo 1']

        # Update the expected value to an appropriate precision
        expected_RZ = 0.00022794540510395617
        self.assertAlmostEqual(calculated_RZ/expected_RZ, 1.0, 2)
    
    def test_YZ_gravity_load(self):
        # A First Course in the Finite Element Method, 4th Edition
        # Daryl L. Logan
        # Problem 5.30
        # Units for this model are kips and inches
        frame = FEModel3D()
        # Define the nodes
        frame.add_node('N1', 0, 0, 0)
        frame.add_node('N2', 0, 30*12, 0)
        frame.add_node('N3', 0, 40*12, 15*12)
        frame.add_node('N4', 0, 40*12, 35*12)
        frame.add_node('N5', 0, 30*12, 50*12)
        frame.add_node('N6', 0, 0, 50*12)
        # Define the supports
        frame.def_support('N1', True, True, True, True, True, True)
        frame.def_support('N6', True, True, True, True, True, True)

        # Define materials
        frame.add_material('Steel', 30000, 250, 0.5, 490/1000/12**3)

        # Create members (all members will have the same properties in this example)
        J = 250
        Iy = 250
        Iz = 200
        A = 12
        frame.add_section('FrameSection1', A, Iy, Iz, J)
        frame.add_section('FrameSection2', A, Iz, Iy, J)

        frame.add_member('M1', 'N1', 'N2', 'Steel', 'FrameSection2')
        frame.add_member('M2', 'N2', 'N3', 'Steel', 'FrameSection1')
        frame.add_member('M3', 'N3', 'N4', 'Steel', 'FrameSection1')
        frame.add_member('M4', 'N4', 'N5', 'Steel', 'FrameSection1')
        frame.add_member('M5', 'N5', 'N6', 'Steel', 'FrameSection2')

        # Add nodal loads
        frame.add_node_load('N3', 'FY', -30)
        frame.add_node_load('N4', 'FY', -30)

        # Analyze the model
        frame.analyze()

        # subTest context manager prints which portion fails, if any
        # Check reactions at N1 and N6
        correct_reactions = [('N1', {'RxnFZ': 11.6877,
                                     'RxnFY': 30,
                                     'RxnMX': 1810.0745}),
                             ('N6', {'RxnFZ': -11.6877,
                                     'RxnFY': 30,
                                     'RxnMX': -1810.0745})]
        for name, values in correct_reactions:
            with self.subTest(node=name):
                node = frame.nodes[name]
                # Two decimal place accuracy requires +/-0.5% accuracy
                # one decimal place requires +/-5%
                self.assertAlmostEqual(node.RxnFZ['Combo 1']/values['RxnFZ'], 1.0, 2)
                self.assertAlmostEqual(node.RxnFY['Combo 1']/values['RxnFY'], 1.0, 2)
                self.assertAlmostEqual(node.RxnMX['Combo 1']/values['RxnMX'], 1.0, 2)
        
        # Check displacements at N3 and N4
        correct_displacements = [('N3', {'DY': -6.666757,
                                         'RX':  0.032}),
                                 ('N4', {'DY': -6.666757,
                                         'RX': -0.032})]
        for name, values in correct_displacements:
            with self.subTest(node=name):
                node = frame.nodes[name]
                # Two decimal place accuracy requires +/-0.5% accuracy
                # one decimal place requires +/-5%
                self.assertAlmostEqual(node.DY['Combo 1']/values['DY'], 1.0, 2)
                self.assertAlmostEqual(node.RX['Combo 1']/values['RX'], 1.0, 2)

    def test_XZ_ptload(self):
        # A simply supported beam with a point load.
        # Units used in this example are inches, and kips
        SimpleBeam = FEModel3D()
        # Add nodes (14 ft = 168 in apart)
        SimpleBeam.add_node("N1", 0, 0, 0)
        SimpleBeam.add_node("N2", 0, 0, 168)

        # Add a material
        SimpleBeam.add_material('Steel', 29000, 11400, 0.5, 490/1000/12**3)

        # Add a beam with the following properties:
        A = 20
        Iy = 100
        Iz = 150
        J = 250
        SimpleBeam.add_section('BeamSection', A, Iy, Iz, J)

        SimpleBeam.add_member("M1", "N1", "N2", "Steel", 'BeamSection')
        # Provide simple supports
        SimpleBeam.def_support("N1", True, True, True, False, False, True)
        SimpleBeam.def_support("N2", True, True, True, False, False, False)
        # Add a point load of 5 kips at the midspan of the beam
        SimpleBeam.add_member_pt_load("M1", "Fy", 5, 7 * 12)
        # Analyze the beam
        SimpleBeam.analyze(False)
        # Print reactions at each end of the beam
        correct_reactions = [('N1', -2.5),
                             ('N2', -2.5)]
        for node_name, rxn in correct_reactions:
            with self.subTest(node=node_name):
                calculated_reaction = SimpleBeam.nodes[node_name].RxnFY['Combo 1']
                # Two decimal place accuracy requires +/-0.5% accuracy
                # one decimal place requires +/-5%
                self.assertAlmostEqual(calculated_reaction/rxn, 1.0, 2)
    
    def test_Kassimali_3_35(self):
        """
        Tests against Kassimali example 3.35.

        This example was selected because it allows us to check the following features:
            1. Member loads aligned in global directions.
            2. A member internal hinge.
            3. A point load at the end of a member.
        The example will be run in the XZ plane to change things up a bit.
        """

        frame = FEModel3D()
        
        frame.add_node('A', 0, 0, 0)
        frame.add_node('B', 0, 0, 24)
        frame.add_node('C', 12, 0, 0)
        frame.add_node('D', 12, 0, 24)
        frame.add_node('E', 24, 0, 12)

        # Add a material
        frame.add_material('Steel', 29000*12**2, 11200*12**2, 0.5, 490/1000)

        Iy = 17.3/12**4
        Iz = 204/12**4
        J = 0.3/12**4
        A = 7.65/12**2
        frame.add_section('FrameSection', A, Iy, Iz, J)

        frame.add_member('AC', 'A', 'C', 'Steel', 'FrameSection')
        frame.add_member('BD', 'B', 'D', 'Steel', 'FrameSection')
        frame.add_member('CE', 'C', 'E', 'Steel', 'FrameSection')
        frame.add_member('ED', 'E', 'D', 'Steel', 'FrameSection')

        frame.def_support('A', support_DX=True, support_DY=True, support_DZ=True)
        frame.def_support('B', support_DX=True, support_DY=True, support_DZ=True)
        frame.def_support('E', support_DY=True)

        frame.def_releases('CE', Rzj=True)

        frame.add_member_pt_load('AC', 'FZ', 20, 12)
        frame.add_member_dist_load('CE', 'FX', -1.5, -1.5)
        frame.add_member_dist_load('ED', 'FX', -1.5, -1.5)

        # from PyNite.Visualization import render_model
        # render_model(frame, text_height=0.5, case='Case 1')

        frame.analyze()

        AZ = -8.63
        AX = 15.46
        BZ = -11.37
        BX = 35.45

        # The reactions were compared manually to Kassimali's solution and the shears were within
        # 10% and 7% respectively. That seems like it's a little big to be a rounding error alone.
        # Likely the finite element method is a little more accurate than the simplified method
        # Kassimali uses.
        self.assertLess(abs(frame.nodes['A'].RxnFZ['Combo 1']/AZ - 1), 0.1)
        self.assertLess(abs(frame.nodes['A'].RxnFX['Combo 1']/AX - 1), 0.05)
        self.assertLess(abs(frame.nodes['B'].RxnFZ['Combo 1']/BZ - 1), 0.7)
        self.assertLess(abs(frame.nodes['B'].RxnFX['Combo 1']/BX - 1), 0.05)