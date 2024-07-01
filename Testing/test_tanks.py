# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
from PyNite.Mesh import CylinderMesh
import sys
from io import StringIO

class Test_Tanks(unittest.TestCase):
    ''' Tests of analyzing plate elements. '''

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__
    
    def test_PCA_7_quad(self):
        """
        Tests against the example from Section 7 of "Circular Concrete Tanks
        without Prestressing" by PCA.
        """

        # Create a new finite element model
        tank_model = FEModel3D()

        H = 20     # Tank wall height (ft)
        D = 54     # Tank inside diameter (ft)
        R = D/2    # Tank inside radius (ft)
        t = 10/12  # Tank wall thickness (ft)

        w = 62.5   # Liquid unit weight (pcf)
        
        fc = 4000                    # Concrete compressive strength (psi)
        E = 57000*(fc)**0.5*(12**2)  # Concrete modulus of elasticity (psf)
        nu = 0.25 #0.17              # Poisson's ratio for concrete
        tank_model.add_material('Concrete', E, 0.4*E, nu, 150)

        mesh_size = 1       # Desired mesh size (ft)
        center = [0, 0, 0]  # Origin (X, Y, Z)
        axis = 'Y'          # Axis of revolution

        tank_model.add_cylinder_mesh('MSH1', mesh_size, R, H, t, 'Concrete', 1, 1, center, axis, element_type='Quad')
        tank_model.Meshes['MSH1'].generate()

        # Add hydrostatic loads to the elements
        for element in tank_model.Quads.values():

            avg_Y = (element.i_node.Y + element.j_node.Y
                   + element.m_node.Y + element.n_node.Y)/4
            
            p = (H - avg_Y)*w

            tank_model.add_quad_surface_pressure(element.name, p)

            #Use default local axes
            element.local_x_axis = None

        # Add fixed supports to the base
        for node in tank_model.Nodes.values():
            if node.Y == 0:
                tank_model.def_support(node.name, True, True, True, True, True, True)

        # Analyze the model
        tank_model.analyze()

        # Max/min moment and max hoop tension as determined by PCA.
        My_max_PCA = 14804/1.3/1.7
        My_min_PCA = -3756/1.3/1.7
        Sx_PCA = 55945/1.3/1.7

        # From Timoshenko Section 117 (p. 485)
        # The Timoshenko solution yields similar results to the PCA solution
        beta = (3*(1 - nu**2)/(R**2*t**2))**0.25  # Equation 275
        My_max_Tim = (1 - 1/(beta*H))*w*R*H*t/(12*(1 - nu**2))**0.5
        Qy_max_Tim = -(w*R*H*t)/(12*(1 - nu**2))**0.5*(2*beta - 1/H)

        My_max = max([element.moment(0, 1)[1, 0] for element in tank_model.Quads.values()])
        My_min = min([element.moment(0, 1)[1, 0] for element in tank_model.Quads.values()])
        Sx = max([element.membrane(0, 0)[0, 0] for element in tank_model.Quads.values()])*t

        # MITC4 element corner stresses are unreliable. Use the maximum
        # reaction at the base of the tank instead.
        RMy = max([node.RxnMX['Combo 1'] for node in tank_model.Nodes.values()])/mesh_size
        
        # Check that the PyNite calculated values are within 2% of expected
        # values.
        self.assertLess(abs(1 - My_max/4900), 0.02, 'Failed quad cylinder flexure test.')
        self.assertLess(abs(1 - RMy/My_max_PCA), 0.02, 'Failed quad cylinder flexure test.')
        self.assertLess(abs(1 - My_min/My_min_PCA), 0.02, 'Failed quad cylinder flexure test.')
        self.assertGreater(My_max, 0, 'Failed quad cylinder sign convention test')
        self.assertLess(abs(1 - Sx/20000), 0.02, 'Failed quad cylinder hoop tension test.')

        # Render the model
        # from PyNite.Visualization import render_model
        # render_model(tank_model, 0.25, True, 100, True, 'My', True, 'Combo 1', labels=False, screenshot=None)

    def test_PCA_7_rect(self):
        """
        Tests against the example from Section 7 of "Circular Concrete Tanks without Prestressing" by PCA.
        """

        # Create a new finite element model
        tank_model = FEModel3D()
        
        H = 20     # Tank wall height (ft)
        D = 54     # Tank inside diameter (ft)
        R = D/2    # Tank inside radius (ft)
        t = 10/12  # Tank wall thickness (ft)

        w = 62.5   # Liquid unit weight (pcf)

        fc = 4000                    # Concrete compressive strength (psi)
        E = 57000*(fc)**0.5*(12**2)  # Concrete modulus of elasticity (psf)
        nu = 0.25 #0.17                    # Poisson's ratio for concrete
        tank_model.add_material('Concrete', E, 0.4*E, nu, 150)

        mesh_size = 2       # Desired mesh size (ft)
        center = [0, 0, 0]  # Origin (X, Y, Z)
        axis = 'Y'          # Axis of revolution

        # Add a cylinder mesh to the model
        tank_model.add_cylinder_mesh('MSH1', mesh_size, R, H, t, 'Concrete', 1, 1, center, axis, element_type='Rect')
        
        # Generate the mesh prior to running so we can work with it
        tank_model.Meshes['MSH1'].generate()

        # Add hydrostatic loads to the elements
        for element in tank_model.Plates.values():

            avg_Y = (element.i_node.Y + element.j_node.Y
                   + element.m_node.Y + element.n_node.Y)/4
            
            p = (H - avg_Y)*w

            tank_model.add_plate_surface_pressure(element.name, p)
        

        # Add fixed supports to the base
        for node in tank_model.Nodes.values():
            if node.Y == 0:
                tank_model.def_support(node.name, True, True, True, True, True, True)

        # Analyze the model
        tank_model.analyze()

        # Max/min moment and max hoop tension as determined by PCA.
        My_max_PCA = 14804/1.3/1.7
        My_min_PCA = -3756/1.3/1.7
        Sx_PCA = 55945/1.3/1.7

        # From Timoshenko Section 117 (p. 485)
        # The Timoshenko solution yields similar results to the PCA solution
        beta = (3*(1 - nu**2)/(R**2*t**2))**0.25  # Equation 275
        My_max_Tim = (1 - 1/(beta*H))*w*R*H*t/(12*(1 - nu**2))**0.5
        Qy_max_Tim = -(w*R*H*t)/(12*(1 - nu**2))**0.5*(2*beta - 1/H)

        My_max = tank_model.Meshes['MSH1'].max_moment('My')
        My_min = tank_model.Meshes['MSH1'].min_moment('My')
        Sx = max([element.membrane(element.width()/2, element.height()/2)[0, 0] for element in tank_model.Plates.values()])*t

        # Check that the PyNite calculated values are within 8% of expected values. With a finer mesh the results are known to converge even closer, but 8% allows the model to run faster.
        self.assertLess(abs(1 - My_max/My_max_PCA), 0.08, 'Failed plate cylinder flexure test.')
        self.assertLess(abs(1 - My_min/My_min_PCA), 0.08, 'Failed plate cylinder flexure test.')
        self.assertGreater(My_max, 0, 'Failed plate cylinder sign convention test')
        self.assertLess(abs(1 - Sx/20000), 0.05, 'Failed plate cylinder hoop tension test.')

        # # Render the model
        # from PyNite.Visualization import render_model
        # render_model(tank_model, 0.25, True, 100, True, 'My', True, 'Combo 1', labels=False, screenshot=None)




