# References used to derive this element:
# 1. "Finite Element Procedures, 2nd Edition", Klaus-Jurgen Bathe
# 2. "A First Course in the Finite Element Method, 4th Edition", Daryl L. Logan
# 3. "Finite Element Analysis Fundamentals", Richard H. Gallagher

import numpy as np
from numpy.linalg import inv, det, norm
from math import sin, cos

class Quad3D():
    """
    An isoparametric general quadrilateral element, formulated by superimposing an isoparametric
    MITC4 bending element with an isoparametric plane stress element. Drilling stability is
    provided by adding a weak rotational spring stiffness at each node. Isotropic behavior is the
    default, but orthotropic in-plane behavior can be modeled by specifying stiffness modification
    factors for the element's local x and y axes.

    This element performs well for thick and thin plates, and for skewed plates. Element center
    stresses and corner FORCES converge rapidly; however, corner STRESSES are more representative
    of center stresses. Minor errors are introduced into the solution due to the drilling
    approximation. Orthotropic behavior is limited to acting along the plate's local axes.
    """

#%%
    def __init__(self, name, i_node, j_node, m_node, n_node, t, material_name, model, kx_mod=1.0,
                 ky_mod=1.0):

        self.name = name
        self.ID = None
        self.type = 'Quad'

        self.i_node = i_node
        self.j_node = j_node
        self.m_node = m_node
        self.n_node = n_node

        self.t = t
        self.kx_mod = kx_mod
        self.ky_mod = ky_mod

        self.pressures = []  # A list of surface pressures [pressure, case='Case 1']
    
        # Quads need a link to the model they belong to
        self.model = model

        # Get material properties for the plate from the model
        try:
            self.E = self.model.Materials[material_name].E
            self.nu = self.model.Materials[material_name].nu
        except:
            raise KeyError('Please define the material ' + str(material_name) + ' before assigning it to plates.')

#%%
    def _local_coords(self):
        '''
        Calculates or recalculates and stores the local (x, y) coordinates for each node of the
        quadrilateral.
        '''

        # Get the global coordinates for each node
        X1, Y1, Z1 = self.m_node.X, self.m_node.Y, self.m_node.Z
        X2, Y2, Z2 = self.n_node.X, self.n_node.Y, self.n_node.Z
        X3, Y3, Z3 = self.i_node.X, self.i_node.Y, self.i_node.Z
        X4, Y4, Z4 = self.j_node.X, self.j_node.Y, self.j_node.Z

        # Following Reference 1, Figure 5.26, node 3 will be used as the
        # origin of the plate's local (x, y) coordinate system. Find the
        # vector from the origin to each node.
        vector_32 = np.array([X2 - X3, Y2 - Y3, Z2 - Z3]).T
        vector_31 = np.array([X1 - X3, Y1 - Y3, Z1 - Z3]).T
        vector_34 = np.array([X4 - X3, Y4 - Y3, Z4 - Z3]).T

        # Define the plate's local x, y, and z axes
        x_axis = vector_34
        z_axis = np.cross(x_axis, vector_32)
        y_axis = np.cross(z_axis, x_axis)

        # Convert the x and y axes into unit vectors
        x_axis = x_axis/norm(x_axis)
        y_axis = y_axis/norm(y_axis)

        # Calculate the local (x, y) coordinates for each node
        self.x1 = np.dot(vector_31, x_axis)
        self.x2 = np.dot(vector_32, x_axis)
        self.x3 = 0
        self.x4 = np.dot(vector_34, x_axis)
        self.y1 = np.dot(vector_31, y_axis)
        self.y2 = np.dot(vector_32, y_axis)
        self.y3 = 0
        self.y4 = np.dot(vector_34, y_axis)

#%%
    def J(self, r, s):
        '''
        Returns the Jacobian matrix for the element
        '''
        
        # Get the local coordinates for the element
        x1, y1, x2, y2, x3, y3, x4, y4 = self.x1, self.y1, self.x2, self.y2, self.x3, self.y3, self.x4, self.y4

        # Return the Jacobian matrix
        return 1/4*np.array([[x1*(s + 1) - x2*(s + 1) + x3*(s - 1) - x4*(s - 1), y1*(s + 1) - y2*(s + 1) + y3*(s - 1) - y4*(s - 1)],
                          [x1*(r + 1) - x2*(r - 1) + x3*(r - 1) - x4*(r + 1), y1*(r + 1) - y2*(r - 1) + y3*(r - 1) - y4*(r + 1)]])

#%%
    def B_kappa(self, r, s):

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with
        # respect to r and s to derivatives with respect to x and y
        dH = np.matmul(inv(self.J(r, s)), 1/4*np.array([[1 + s, -1 - s, -1 + s,  1 - s],                 
                                                  [1 + r,  1 - r, -1 + r, -1 - r]]))
        
        # Row 1 = d(beta_x)/dx divided by the local displacement vector 'u'
        # Row 2 = d(beta_y)/dy divided by the local displacement vector 'u'
        # Row 3 = d(beta_x)/dy + d(beta_y)/dx divided by the local displacement vector 'u'
        # Note that beta_x is a function of -theta_y and beta_y is a function of +theta_x (Equations 5.99, p. 423)
        B_kappa = np.array([[0,    0,     -dH[0, 0], 0,    0,     -dH[0, 1], 0,    0,     -dH[0, 2], 0,    0,     -dH[0, 3]],
                         [0, dH[1, 0],     0,     0, dH[1, 1],     0,     0, dH[1, 2],     0,     0, dH[1, 3],     0    ],
                         [0, dH[0, 0], -dH[1, 0], 0, dH[0, 1], -dH[1, 1], 0, dH[0, 2], -dH[1, 2], 0, dH[0, 3], -dH[1, 3]]])
        
        # Below is the matrix derived from the 1984 version of the MITC4 element. It appears to be
        # the same, but with a different sign convention for the section curvatures.
        # B_kappa = np.array([[0,     0,     dH[0, 0],  0,     0,     dH[0, 1],  0,     0,     dH[0, 2],  0,     0,     dH[0, 3]],
        #                  [0,  dH[1, 0],     0,     0,  dH[1, 1],     0,     0,  dH[1, 2],     0,     0,  dH[1, 3],     0   ],
        #                  [0, -dH[0, 0], dH[1, 0],  0, -dH[0, 1], dH[1, 1],  0, -dH[0, 2], dH[1, 2],  0, -dH[0, 3], dH[1, 3]]])

        return B_kappa

#%%
    def B_gamma(self, r, s):
        '''
        Returns the [B] matrix for shear.

        This is provided for reference only and is not actually used by
        PyNite. This is the theoretical solution, but it is known to
        produce spurious shear forces. It is prone to a phenomenon called
        shear locking. Instead of this matrix, the MITC4 [B] matrix is used,
        which eliminates shear-locking and can be used for thick and thin
        plates.
        '''

        H = 1/4*np.array([(1 + r)*(1 + s), (1 - r)*(1 + s), (1 - r)*(1 - s), (1 + r)*(1 - s)])

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with respect to r and s
        # to derivatives with respect to x and y
        dH = np.matmul(inv(self.J(r, s)), 1/4*np.array([[1 + s, -1 - s, -1 + s,  1 - s],                 
                                                  [1 + r,  1 - r, -1 + r, -1 - r]]))

        # Row 1 = d(beta_x)/dx divided by the local displacement vector 'u'
        # Row 2 = d(beta_y)/dy divided by the local displacement vector 'u'
        # Row 3 = d(beta_x)/dy + d(beta_y)/dx divided by the local displacement vector 'u'
        # Note that beta_x is a function of -theta_y and beta_y is a function of +theta_x (Equations 5.99, p. 423)
        B_gamma = np.array([[dH[0, 0],   0,   H[0], dH[0, 1],   0,   H[1], dH[0, 2],   0,   H[2], dH[0, 3],   0,   H[3]],
                         [dH[1, 0], -H[0],  0,   dH[1, 1], -H[1],  0,   dH[1, 2], -H[2],  0,   dH[1, 3], -H[3],  0  ]])
        
        return B_gamma
    
    def B_gamma_MITC4(self, r, s):
        '''
        Returns the [B] matrix for shear.

        MITC stands for mixed interpolation tensoral components. MITC elements
        are used in many programs and are known to perform well for thick and
        thin plates, and for distorted plate geometries.
        '''

        # Get the local coordinates for the element
        x1, y1, x2, y2, x3, y3, x4, y4 = self.x1, self.y1, self.x2, self.y2, self.x3, self.y3, self.x4, self.y4
        x_axis = np.array([1, 0, 0]).T

        # Reference 1, Equations 5.105
        Ax = x1 - x2 - x3 + x4
        Bx = x1 - x2 + x3 - x4
        Cx = x1 + x2 - x3 - x4
        Ay = y1 - y2 - y3 + y4
        By = y1 - y2 + y3 - y4
        Cy = y1 + y2 - y3 - y4

        # Find the angles between the axes of the natural coordinate system and
        # the local x-axis.
        r_axis = np.array([(x1 + x4)/2 - (x2 + x3)/2, (y1 + y4)/2 - (y2 + y3)/2, 0]).T
        s_axis = np.array([(x1 + x2)/2 - (x3 + x4)/2, (y1 + y2)/2 - (y3 + y4)/2, 0]).T

        r_axis = r_axis/norm(r_axis)
        s_axis = s_axis/norm(s_axis)

        alpha = np.arccos(np.dot(r_axis, x_axis))
        beta = np.arccos(np.dot(s_axis, x_axis))
        # alpha = atan(Ay/Ax)
        # beta = pi/2 - atan(Cx/Cy)
        
        # Reference 1, Equations 5.103 and 5.104 (p. 426)
        det_J = det(self.J(r, s))

        gr = ((Cx + r*Bx)**2 + (Cy + r*By)**2)**0.5/(8*det_J)
        gs = ((Ax + s*Bx)**2 + (Ay + s*By)**2)**0.5/(8*det_J)

        # d      =           [    w1           theta_x1             theta_y1             w2            theta_x2              theta_y2            w3             theta_x3             theta_y3         w4             theta_x4             theta_y4      ]
        gamma_rz = gr*np.array([[(1 + s)/2, -(y1 - y2)/4*(1 + s), (x1 - x2)/4*(1 + s), -(1 + s)/2,  -(y1 - y2)/4*(1 + s), (x1 - x2)/4*(1 + s), -(1 - s)/2, -(y4 - y3)/4*(1 - s), (x4 - x3)/4*(1 - s), (1 - s)/2,  -(y4 - y3)/4*(1 - s), (x4 - x3)/4*(1 - s)]])
        gamma_sz = gs*np.array([[(1 + r)/2, -(y1 - y4)/4*(1 + r), (x1 - x4)/4*(1 + r),  (1 - r)/2,  -(y2 - y3)/4*(1 - r), (x2 - x3)/4*(1 - r), -(1 - r)/2, -(y2 - y3)/4*(1 - r), (x2 - x3)/4*(1 - r), -(1 + r)/2, -(y1 - y4)/4*(1 + r), (x1 - x4)/4*(1 + r)]])
        
        # Reference 1, Equations 5.102
        B_gamma_MITC4 = np.zeros((2, 12))
        B_gamma_MITC4[0, :] = gamma_rz*sin(beta) - gamma_sz*sin(alpha)
        B_gamma_MITC4[1, :] = -gamma_rz*cos(beta) + gamma_sz*cos(alpha)
        
        # Return the [B] matrix for shear
        return B_gamma_MITC4

#%%
    def B_m(self, r, s):

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with
        # respect to r and s to derivatives with respect to x and y
        dH = np.matmul(inv(self.J(r, s)), 1/4*np.array([[s + 1, -s - 1, s - 1, -s + 1],                 
                                                  [r + 1, -r + 1, r - 1, -r - 1]]))

        # Reference 1, Example 5.5 (page 353)
        B_m = np.array([[dH[0, 0],    0,     dH[0, 1],    0,     dH[0, 2],    0,     dH[0, 3],    0    ],
                     [   0,     dH[1, 0],    0,     dH[1, 1],    0,     dH[1, 2],    0,     dH[1, 3]],
                     [dH[1, 0], dH[0, 0], dH[1, 1], dH[0, 1], dH[1, 2], dH[0, 2], dH[1, 3], dH[0, 3]]])

        return B_m

#%%
    def Cb(self):
        '''
        Returns the stress-strain matrix for plate bending.
        '''

        # Referemce 1, Table 4.3, page 194
        nu = self.nu
        E = self.E
        h = self.t

        Cb = E*h**3/(12*(1 - nu**2))*np.array([[1,  nu,      0    ],
                                            [nu, 1,       0    ],
                                            [0,  0,  (1 - nu)/2]])
        
        return Cb

#%%
    def Cs(self):
        '''
        Returns the stress-strain matrix for shear.
        '''
        # Reference 1, Equations (5.97), page 422
        k = 5/6
        E = self.E
        h = self.t
        nu = self.nu

        Cs = E*h*k/(2*(1 + nu))*np.array([[1, 0],
                                       [0, 1]])

        return Cs

#%%
    def Cm(self):
        """
        Returns the stress-strain matrix for an isotropic or orthotropic plane stress element
        """
        
        # Apply the stiffness modification factors for each direction to obtain orthotropic
        # behavior. Stiffness modification factors of 1.0 in each direction (the default) will
        # model isotropic behavior. Orthotropic behavior is limited to the element's local
        # coordinate system.
        Ex = self.E*self.kx_mod
        Ey = self.E*self.ky_mod
        nu_xy = self.nu
        nu_yx = self.nu

        # The shear modulus will be unafected by orthotropic behavior
        # Logan, Appendix C.3, page 750
        G = self.E/(2*(1 + self.nu))

        # Gallagher, Equation 9.3, page 251
        Cm = 1/(1 - nu_xy*nu_yx)*np.array([[   Ex,    nu_yx*Ex,           0         ],
                                        [nu_xy*Ey,    Ey,              0         ],
                                        [    0,        0,     (1 - nu_xy*nu_yx)*G]])
        
        return Cm

#%%
    def k_b(self):
        '''
        Returns the local stiffness matrix for bending stresses
        '''

        Cb = self.Cb()
        Cs = self.Cs()

        # Define the gauss point for numerical integration
        gp = 1/3**0.5

        # Get the determinant of the Jacobian matrix for each gauss pointing 
        # Doing this now will save us from doing it twice below
        J1 = det(self.J(gp, gp))
        J2 = det(self.J(-gp, gp))
        J3 = det(self.J(-gp, -gp))
        J4 = det(self.J(gp, -gp))

        # Get the bending B matrices for each gauss point
        B1 = self.B_kappa(gp, gp)
        B2 = self.B_kappa(-gp, gp)
        B3 = self.B_kappa(-gp, -gp)
        B4 = self.B_kappa(gp, -gp)

        # Create the stiffness matrix with bending stiffness terms
        # See Reference 1, Equation 5.94
        k = (np.matmul(B1.T, np.matmul(Cb, B1))*J1 +
             np.matmul(B2.T, np.matmul(Cb, B2))*J2 +
             np.matmul(B3.T, np.matmul(Cb, B3))*J3 +
             np.matmul(B4.T, np.matmul(Cb, B4))*J4)

        # Get the MITC4 shear B matrices for each gauss point
        B1 = self.B_gamma_MITC4(gp, gp)
        B2 = self.B_gamma_MITC4(-gp, gp)
        B3 = self.B_gamma_MITC4(-gp, -gp)
        B4 = self.B_gamma_MITC4(gp, -gp)
        
        # Alternatively the shear B matrix below could be used. However, this matrix is prone to
        # shear locking and will overestimate the stiffness.
        # B1 = self.B_gamma(gp, gp)
        # B2 = self.B_gamma(-gp, gp)
        # B3 = self.B_gamma(-gp, -gp)
        # B4 = self.B_gamma(gp, -gp)

        # Add shear stiffness terms to the stiffness matrix
        k += (np.matmul(B1.T, np.matmul(Cs, B1))*J1 +
              np.matmul(B2.T, np.matmul(Cs, B2))*J2 +
              np.matmul(B3.T, np.matmul(Cs, B3))*J3 +
              np.matmul(B4.T, np.matmul(Cs, B4))*J4)
        
        # Following Bathe's recommendation for the drilling degree of freedom
        # from Example 4.19 in "Finite Element Procedures, 2nd Ed.", calculate
        # the drilling stiffness as 1/1000 of the smallest diagonal term in
        # the element's stiffness matrix. This is not theoretically correct,
        # but it allows the model to solve without singularities, and should
        # have a minimal effect on the final solution. Bathe recommends 1/1000
        # as a value that is weak enough but not so small that it affect the
        # results. Bathe recommends looking at all the diagonals in the
        # combined bending plus membrane stiffness matrix. Some of those terms
        # relate to translational stiffness. It seems more rational to only
        # look at the terms relating to rotational stiffness. That will be
        # PyNite's approach.
        k_rz = min(abs(k[1, 1]), abs(k[2, 2]), abs(k[4, 4]), abs(k[5, 5]),
                   abs(k[7, 7]), abs(k[8, 8]), abs(k[10, 10]), abs(k[11, 11])
                   )/1000
        
        # Initialize the expanded stiffness matrix to all zeros
        k_exp = np.zeros((24, 24))

        # Step through each term in the unexpanded stiffness matrix
        # i = Unexpanded matrix row
        for i in range(12):

            # j = Unexpanded matrix column
            for j in range(12):
                
                # Find the corresponding term in the expanded stiffness
                # matrix

                # m = Expanded matrix row
                if i in [0, 3, 6, 9]:  # indices associated with deflection in z
                    m = 2*i + 2
                if i in [1, 4, 7, 10]:  # indices associated with rotation about x
                    m = 2*i + 1
                if i in [2, 5, 8, 11]:  # indices associated with rotation about y
                    m = 2*i

                # n = Expanded matrix column
                if j in [0, 3, 6, 9]:  # indices associated with deflection in z
                    n = 2*j + 2
                if j in [1, 4, 7, 10]:  # indices associated with rotation about x
                    n = 2*j + 1
                if j in [2, 5, 8, 11]:  # indices associated with rotation about y
                    n = 2*j
                
                # Ensure the indices are integers rather than floats
                m, n = round(m), round(n)

                # Add the term from the unexpanded matrix into the expanded
                # matrix
                k_exp[m, n] = k[i, j]

        # Add the drilling degree of freedom's weak spring
        k_exp[5, 5] = k_rz
        k_exp[11, 11] = k_rz
        k_exp[17, 17] = k_rz
        k_exp[23, 23] = k_rz

        return k_exp

#%%
    def k_m(self):
        '''
        Returns the local stiffness matrix for membrane (in-plane) stresses.

        Plane stress is assumed
        '''

        t = self.t
        Cm = self.Cm()

        # Define the gauss point for numerical integration
        gp = 1/3**0.5

        # Get the membrane B matrices for each gauss point
        # Doing this now will save us from doing it twice below
        B1 = self.B_m(gp, gp)
        B2 = self.B_m(-gp, gp)
        B3 = self.B_m(-gp, -gp)
        B4 = self.B_m(gp, -gp)

        # See reference 1 at the bottom of page 353, and reference 2 page 466
        k = t*(np.matmul(B1.T, np.matmul(Cm, B1))*det(self.J(gp, gp)) +
               np.matmul(B2.T, np.matmul(Cm, B2))*det(self.J(-gp, gp)) +
               np.matmul(B3.T, np.matmul(Cm, B3))*det(self.J(-gp, -gp)) +
               np.matmul(B4.T, np.matmul(Cm, B4))*det(self.J(gp, -gp)))
        
        k_exp = np.zeros((24, 24))

        # Step through each term in the unexpanded stiffness matrix
        # i = Unexpanded matrix row
        for i in range(8):

            # j = Unexpanded matrix column
            for j in range(8):
                
                # Find the corresponding term in the expanded stiffness
                # matrix

                # m = Expanded matrix row
                if i in [0, 2, 4, 6]:  # indices associated with displacement in x
                    m = i*3
                if i in [1, 3, 5, 7]:  # indices associated with displacement in y
                    m = i*3 - 2

                # n = Expanded matrix column
                if j in [0, 2, 4, 6]:  # indices associated with displacement in x
                    n = j*3
                if j in [1, 3, 5, 7]:  # indices associated with displacement in y
                    n = j*3 - 2
                
                # Ensure the indices are integers rather than floats
                m, n = round(m), round(n)

                # Add the term from the unexpanded matrix into the expanded matrix
                k_exp[m, n] = k[i, j]
        
        return k_exp

#%%
    def k(self):
        '''
        Returns the quad element's local stiffness matrix.
        '''

        # Recalculate the local coordinate system
        self._local_coords()

        # Sum the bending and membrane stiffness matrices
        return np.add(self.k_b(), self.k_m())

#%%   
    def f(self, combo_name='Combo 1'):
        """
        Returns the quad element's local end force vector
        """
        
        # Calculate and return the plate's local end force vector
        return np.add(np.matmul(self.k(), self.d(combo_name)), self.fer(combo_name))

#%%
    def fer(self, combo_name='Combo 1'):
        '''
        Returns the quadrilateral's local fixed end reaction vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the consistent load vector for.
        '''
        
        Hw = lambda r, s : 1/4*np.array([[(1 + r)*(1 + s), 0, 0, (1 - r)*(1 + s), 0, 0, (1 - r)*(1 - s), 0, 0, (1 + r)*(1 - s), 0, 0]])

        # Initialize the fixed end reaction vector
        fer = np.zeros((12,1))

        # Get the requested load combination
        combo = self.model.LoadCombos[combo_name]

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # Initialize the element's surface pressure to zero
        p = 0

        # Loop through each load case and factor in the load combination
        for case, factor in combo.factors.items():

            # Sum the pressures
            for pressure in self.pressures:

                # Check if the current pressure corresponds to the current load case
                if pressure[1] == case:

                    # Sum the pressures
                    p -= factor*pressure[0]
        
        fer = (Hw(-gp, -gp).T*p*det(self.J(-gp, -gp))
             + Hw(-gp, gp).T*p*det(self.J(-gp, gp))
             + Hw(gp, gp).T*p*det(self.J(gp, gp))
             + Hw(gp, -gp).T*p*det(self.J(gp, -gp)))

        # Initialize the expanded vector to all zeros
        fer_exp = np.zeros((24, 1))

        # Step through each term in the unexpanded vector
        # i = Unexpanded vector row
        for i in range(12):
                
            # Find the corresponding term in the expanded vector

            # m = Expanded vector row
            if i in [0, 3, 6, 9]:   # indices associated with deflection in z
                m = 2*i + 2
            if i in [1, 4, 7, 10]:  # indices associated with rotation about x
                m = 2*i + 1
            if i in [2, 5, 8, 11]:  # indices associated with rotation about y
                m = 2*i
                
            # Ensure the index is an integer rather than a float
            m = round(m)

            # Add the term from the unexpanded vector into the expanded vector
            fer_exp[m, 0] = fer[i, 0]

        return fer_exp

#%%
    def d(self, combo_name='Combo 1'):
       """
       Returns the quad element's local displacement vector
       """

       # Calculate and return the local displacement vector
       return np.matmul(self.T(), self.D(combo_name))

#%%
    def F(self, combo_name='Combo 1'):
        """
        Returns the quad element's global force vector

        Parameters
        ----------
        combo_name : string
            The load combination to get results for.
        """
        
        # Calculate and return the global force vector
        return np.matmul(inv(self.T()), self.f(combo_name))

#%%
    def D(self, combo_name='Combo 1'):
        '''
        Returns the quad element's global displacement vector for the given
        load combination.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the displacement vector
            for (not the load combination itself).
        '''
        
        # Initialize the displacement vector
        D = np.zeros((24, 1))
        
        # Read in the global displacements from the nodes
        D[0, 0] = self.m_node.DX[combo_name]
        D[1, 0] = self.m_node.DY[combo_name]
        D[2, 0] = self.m_node.DZ[combo_name]
        D[3, 0] = self.m_node.RX[combo_name]
        D[4, 0] = self.m_node.RY[combo_name]
        D[5, 0] = self.m_node.RZ[combo_name]

        D[6, 0] = self.n_node.DX[combo_name]
        D[7, 0] = self.n_node.DY[combo_name]
        D[8, 0] = self.n_node.DZ[combo_name]
        D[9, 0] = self.n_node.RX[combo_name]
        D[10, 0] = self.n_node.RY[combo_name]
        D[11, 0] = self.n_node.RZ[combo_name]

        D[12, 0] = self.i_node.DX[combo_name]
        D[13, 0] = self.i_node.DY[combo_name]
        D[14, 0] = self.i_node.DZ[combo_name]
        D[15, 0] = self.i_node.RX[combo_name]
        D[16, 0] = self.i_node.RY[combo_name]
        D[17, 0] = self.i_node.RZ[combo_name]

        D[18, 0] = self.j_node.DX[combo_name]
        D[19, 0] = self.j_node.DY[combo_name]
        D[20, 0] = self.j_node.DZ[combo_name]
        D[21, 0] = self.j_node.RX[combo_name]
        D[22, 0] = self.j_node.RY[combo_name]
        D[23, 0] = self.j_node.RZ[combo_name]
        
        # Return the global displacement vector
        return D

#%%
    def K(self):
        '''
        Returns the quad element's global stiffness matrix
        '''

        # Get the transformation matrix
        T = self.T()

        # Calculate and return the stiffness matrix in global coordinates
        return np.matmul(np.matmul(inv(T), self.k()), T)

#%% 
    # Global fixed end reaction vector
    def FER(self, combo_name='Combo 1'):
        '''
        Returns the global fixed end reaction vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to calculate the fixed end
            reaction vector for (not the load combination itself).
        '''
        
        # Calculate and return the fixed end reaction vector
        return np.matmul(inv(self.T()), self.fer(combo_name))

#%%  
    def T(self):
        '''
        Returns the coordinate transformation matrix for the quad element.
        '''

        xi = self.i_node.X
        xj = self.j_node.X
        yi = self.i_node.Y
        yj = self.j_node.Y
        zi = self.i_node.Z
        zj = self.j_node.Z

        # Calculate the direction cosines for the local x-axis.The local x-axis will run from
        # the i-node to the j-node
        x = [xj - xi, yj - yi, zj - zi]

        # Divide the vector by its magnitude to produce a unit x-vector of
        # direction cosines
        mag = (x[0]**2 + x[1]**2 + x[2]**2)**0.5
        x = [x[0]/mag, x[1]/mag, x[2]/mag]
        
        # The local y-axis will be in the plane of the plate. Find a vector in
        # the plate's local xy plane.
        xn = self.n_node.X
        yn = self.n_node.Y
        zn = self.n_node.Z
        xy = [xn - xi, yn - yi, zn - zi]

        # Find a vector perpendicular to the plate surface to get the
        # orientation of the local z-axis.
        z = np.cross(x, xy)
        
        # Divide the z-vector by its magnitude to produce a unit z-vector of
        # direction cosines.
        mag = (z[0]**2 + z[1]**2 + z[2]**2)**0.5
        z = [z[0]/mag, z[1]/mag, z[2]/mag]

        # Calculate the local y-axis as a vector perpendicular to the local z
        # and x-axes.
        y = np.cross(z, x)
        
        # Divide the y-vector by its magnitude to produce a unit vector of
        # direction cosines.
        mag = (y[0]**2 + y[1]**2 + y[2]**2)**0.5
        y = [y[0]/mag, y[1]/mag, y[2]/mag]

        # Create the direction cosines matrix.
        dirCos = np.array([x,
                        y,
                        z])
        
        # Build the transformation matrix.
        T = np.zeros((24, 24))
        T[0:3, 0:3] = dirCos
        T[3:6, 3:6] = dirCos
        T[6:9, 6:9] = dirCos
        T[9:12, 9:12] = dirCos
        T[12:15, 12:15] = dirCos
        T[15:18, 15:18] = dirCos
        T[18:21, 18:21] = dirCos
        T[21:24, 21:24] = dirCos
        
        # Return the transformation matrix.
        return T

#%%
    def shear(self, r=0, s=0, local=True, combo_name='Combo 1'):
        '''
        Returns the interal shears at any point in the quad element.

        Internal shears are reported as a 2D array [[Qx], [Qy]] at the
        specified location in the (r, s) natural coordinate system.

        Parameters
        ----------
        r : number
            The r-coordinate. Default is 0.
        s : number
            The s-coordinate. Default is 0.
        
        Returns
        -------
        Internal shear force per unit length of the quad element.
        '''

        # Get the plate's local displacement vector
        # Slice out terms not related to plate bending
        d = self.d(combo_name)[[2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22], :]

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # Define extrapolated r and s points
        r_ex = r/gp
        s_ex = s/gp

        # Define the interpolation functions
        H = 1/4*np.array([(1 + r_ex)*(1 + s_ex), (1 - r_ex)*(1 + s_ex), (1 - r_ex)*(1 - s_ex), (1 + r_ex)*(1 - s_ex)])

        # Get the stress-strain matrix
        Cs = self.Cs()

        # Calculate the internal shears [Qx, Qy] at each gauss point
        q1 = np.matmul(Cs, np.matmul(self.B_gamma_MITC4(gp, gp), d))
        q2 = np.matmul(Cs, np.matmul(self.B_gamma_MITC4(-gp, gp), d))
        q3 = np.matmul(Cs, np.matmul(self.B_gamma_MITC4(-gp, -gp), d))
        q4 = np.matmul(Cs, np.matmul(self.B_gamma_MITC4(gp, -gp), d))

        # Extrapolate to get the value at the requested location
        Qx = H[0]*q1[0] + H[1]*q2[0] + H[2]*q3[0] + H[3]*q4[0]
        Qy = H[0]*q1[1] + H[1]*q2[1] + H[2]*q3[1] + H[3]*q4[1]

        if local:
            
            return np.array([Qx,
                             Qy])
        
        else:
            
            # Get the direction cosines for the plate's local coordinate system
            dir_cos = self.T()[:3, :3]

            return np.matmul(np.linalg.inv(dir_cos), np.array([Qx,
                                                               Qy,
                                                               [0]]))

#%%   
    def moment(self, r=0, s=0, local=True, combo_name='Combo 1'):
        '''
        Returns the interal moments at any point in the quad element.

        Internal moments are reported as a 2D array [[Mx], [My], [Mxy]] at the
        specified location in the (r, s) natural coordinate system.

        Parameters
        ----------
        r : number
            The r-coordinate. Default is 0.
        s : number
            The s-coordinate. Default is 0.
        
        Returns
        -------
        Internal moment per unit length of the quad element.
        '''

        # Get the plate's local displacement vector
        # Slice out terms not related to plate bending
        d = self.d(combo_name)[[2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22], :]

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # # Define extrapolated r and s points
        r_ex = r/gp
        s_ex = s/gp

        # Define the interpolation functions
        H = 1/4*np.array([(1 + r_ex)*(1 + s_ex), (1 - r_ex)*(1 + s_ex), (1 - r_ex)*(1 - s_ex), (1 + r_ex)*(1 - s_ex)])

        # Get the stress-strain matrix
        Cb = self.Cb()

        # Calculate the internal moments [Mx, My, Mxy] at each gauss point
        m1 = np.matmul(Cb, np.matmul(self.B_kappa(gp, gp), d))
        m2 = np.matmul(Cb, np.matmul(self.B_kappa(-gp, gp), d))
        m3 = np.matmul(Cb, np.matmul(self.B_kappa(-gp, -gp), d))
        m4 = np.matmul(Cb, np.matmul(self.B_kappa(gp, -gp), d))

        # Extrapolate to get the value at the requested location
        Mx = H[0]*m1[0] + H[1]*m2[0] + H[2]*m3[0] + H[3]*m4[0]
        My = H[0]*m1[1] + H[1]*m2[1] + H[2]*m3[1] + H[3]*m4[1]
        Mxy = H[0]*m1[2] + H[1]*m2[2] + H[2]*m3[2] + H[3]*m4[2]
        
        if local:

            return np.array([Mx,
                             My,
                             Mxy])
        
        else:
            
            # Get the direction cosines for the plate's local coordinate system
            dir_cos = self.T()[:3, :3]

            # Convert the plate flexural stresses to global coordinates
            return np.matmul(np.linalg.inv(dir_cos), np.array([Mx,
                                                               My,
                                                               [0]]))

#%%
    def membrane(self, r=0, s=0, local=True, combo_name='Combo 1'):
        
        # Get the plate's local displacement vector
        # Slice out terms not related to membrane forces
        d = self.d(combo_name)[[0, 1, 6, 7, 12, 13, 18, 19], :]

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # Define extrapolated r and s points
        r_ex = r/gp
        s_ex = s/gp

        # Define the interpolation functions
        H = 1/4*np.array([(1 + r_ex)*(1 + s_ex), (1 - r_ex)*(1 + s_ex), (1 - r_ex)*(1 - s_ex), (1 + r_ex)*(1 - s_ex)])

        # Get the stress-strain matrix
        Cm = self.Cm()
        
        # Calculate the internal stresses [Sx, Sy, Txy] at each gauss point
        s1 = np.matmul(Cm, np.matmul(self.B_m(gp, gp), d))
        s2 = np.matmul(Cm, np.matmul(self.B_m(-gp, gp), d))
        s3 = np.matmul(Cm, np.matmul(self.B_m(-gp, -gp), d))
        s4 = np.matmul(Cm, np.matmul(self.B_m(gp, -gp), d))

        # Extrapolate to get the value at the requested location
        Sx = H[0]*s1[0] + H[1]*s2[0] + H[2]*s3[0] + H[3]*s4[0]
        Sy = H[0]*s1[1] + H[1]*s2[1] + H[2]*s3[1] + H[3]*s4[1]
        Txy = H[0]*s1[2] + H[1]*s2[2] + H[2]*s3[2] + H[3]*s4[2]

        if local:

            return np.array([Sx,
                             Sy,
                             Txy])
        
        else:

            # Get the direction cosines for the plate's local coordinate system
            dir_cos = self.T()[:3, :3]

            # Convert the plate membrane stresses to global coordinates
            return np.matmul(np.linalg.inv(dir_cos), np.array([Sx,
                                                               Sy,
                                                               [0]]))

