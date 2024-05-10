import numpy as np
from numpy import zeros, array, matmul, cross, add
from numpy.linalg import inv, norm, det

#%%
class Plate3D():

    def __init__(self, name, i_node, j_node, m_node, n_node, t, material_name, model, kx_mod=1.0,
                 ky_mod=1.0):
        """
        A rectangular plate element

        Parameters
        ----------
        name : string
            A unique plate name
        i_node : Node3D
            The plate's i-node
        j_node : Node3D
            The plate's j-node
        m_node : Node3D
            The plate's m-node
        n_node : Node3D
            The plate's n-node
        t : number
            Plate thickness
        material_name : string
            The name of the plate material
        kx_mod : number
            Modification factor for stiffness in the plate's local x-direction. Default value is
            1.0, which indicates no stiffness modification (100% stiffness).
        ky_mod : number
            Modification factor for stiffness in the plate's local y-direction. Default value is
            1.0, which indicates no stiffness modification (100% stiffness).
        model : FEModel3D
            The model the plate is a part of
        """

        self.name = name
        self.ID = None
        self.type = 'Rect'

        self.i_node = i_node
        self.j_node = j_node
        self.m_node = m_node
        self.n_node = n_node

        self.t = t
        
        self.kx_mod = kx_mod
        self.ky_mod = ky_mod

        self.pressures = []  # A list of surface pressures [pressure, case='Case 1']

        # Plates need a link to the model they belong to
        self.model = model

        # Get material properties for the plate from the model
        try:
            self.E = self.model.Materials[material_name].E
            self.nu = self.model.Materials[material_name].nu
        except:
            raise KeyError('Please define the material ' + str(material_name) + ' before assigning it to plates.')
    
    def width(self):
        """
        Returns the width of the plate along its local x-axis
        """
        return self.i_node.distance(self.j_node)

    def height(self):
        """
        Returns the height of the plate along its local y-axis
        """
        return self.i_node.distance(self.n_node)
    
    def Dm(self):
        """
        Returns the plane stress constitutive matrix for orthotropic materials [Dm]
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

        # Calculate the constitutive matrix [Dm]
        Dm = 1/(1 - nu_xy*nu_yx)*array([[   Ex,    nu_yx*Ex,            0             ],
                                        [nu_xy*Ey,    Ey,               0             ],
                                        [   0,        0,     (1 - nu_xy*nu_yx)*G]])

        # Return the constitutive matrix [Dm]
        return Dm

    def Db(self):
        """
        Returns the bending constitutive matrix for orthotropic materials [Db]
        """

        # Get the plate thickness and other material parameters
        t = self.t
        nu_xy = self.nu
        nu_yx = self.nu
        Ex = self.E*self.kx_mod
        Ey = self.E*self.ky_mod
        G = self.E/(2*(1 + self.nu))

        # Calculate the constitutive matrix [Db]
        Db = t**3/(12*(1 - nu_xy*nu_yx))*array([[   Ex,    nu_yx*Ex, 0],
                                                [nu_xy*Ey,    Ey,    0],
                                                [   0,        0,     G]])

        # Return the constitutive matrix [Db]
        return Db

    def J(self, r, s):
        '''
        Returns the Jacobian matrix for the element
        '''
        
        # Get the local coordinates for the element's nodes
        x1, y1, x2, y2, x3, y3, x4, y4 = 0, 0, self.width(), 0, self.width(), self.height(), 0, self.height()

        # Return the Jacobian matrix
        return 1/4*array([[x1*(s - 1) - x2*(s - 1) + x3*(s + 1) - x4*(s + 1), y1*(s - 1) - y2*(s - 1) + y3*(s + 1) - y4*(s + 1)],
                          [x1*(r - 1) - x2*(r + 1) + x3*(r + 1) - x4*(r - 1), y1*(r - 1) - y2*(r + 1) + y3*(r + 1) - y4*(r - 1)]])
    
    def B_m(self, r, s):

        # Differentiate the interpolation functions
        # Row 1 = interpolation functions differentiated with respect to x
        # Row 2 = interpolation functions differentiated with respect to y
        # Note that the inverse of the Jacobian converts from derivatives with
        # respect to r and s to derivatives with respect to x and y
        dH = matmul(inv(self.J(r, s)), 1/4*array([[s - 1, -s + 1, s + 1, -s - 1],                 
                                                  [r - 1, -r - 1, r + 1, -r + 1]]))

        # Reference 1, Example 5.5 (page 353)
        B_m = array([[dH[0, 0],    0,     dH[0, 1],    0,     dH[0, 2],    0,     dH[0, 3],    0    ],
                     [   0,     dH[1, 0],    0,     dH[1, 1],    0,     dH[1, 2],    0,     dH[1, 3]],
                     [dH[1, 0], dH[0, 0], dH[1, 1], dH[0, 1], dH[1, 2], dH[0, 2], dH[1, 3], dH[0, 3]]])
        
        return B_m
    
    def k(self):
        """
        returns the plate's local stiffness matrix
        """
        return add(self.k_b(), self.k_m())
    
    def k_m(self):
        '''
        Returns the local stiffness matrix for membrane (in-plane) stresses.

        Plane stress is assumed
        '''

        t = self.t
        Dm = self.Dm()

        # Define the gauss point for numerical integration
        gp = 1/3**0.5

        # Get the membrane B matrices for each gauss point
        # Doing this now will save us from doing it twice below
        B1 = self.B_m(-gp, -gp)
        B2 = self.B_m(gp, -gp)
        B3 = self.B_m(gp, gp)
        B4 = self.B_m(-gp, gp)

        # See reference 1 at the bottom of page 353, and reference 2 page 466
        k = t*(matmul(B1.T, matmul(Dm, B1))*det(self.J(-gp, -gp)) +
               matmul(B2.T, matmul(Dm, B2))*det(self.J(gp, -gp)) +
               matmul(B3.T, matmul(Dm, B3))*det(self.J(gp, gp)) +
               matmul(B4.T, matmul(Dm, B4))*det(self.J(-gp, gp)))
        
        k_exp = zeros((24, 24))

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
    
    def k_b(self):
        """
        Returns the local stiffness matrix for bending
        """

        b = self.width()/2
        c = self.height()/2

        Ex = self.E
        Ey = self.E
        nu_xy = self.nu
        nu_yx = self.nu
        G = self.E/(2*(1 + self.nu))
        t = self.t

        # Stiffness matrix for plate bending. This matrix was derived using a jupyter notebook. The
        # notebook can be found in the `Derivations`` folder of this project.
        k = t**3/12*array([[(-Ex*nu_yx*b**2*c**2/4 - Ex*c**4 - Ey*nu_xy*b**2*c**2/4 - Ey*b**4 + 7*G*nu_xy*nu_yx*b**2*c**2/5 - 7*G*b**2*c**2/5)/(b**3*c**3*(nu_xy*nu_yx - 1)), (-Ex*nu_yx*c**2/2 - Ey*b**2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (Ex*c**2 + Ey*nu_xy*b**2/2 - G*nu_xy*nu_yx*b**2/5 + G*b**2/5)/(b**2*c*(nu_xy*nu_yx - 1)), (5*Ex*nu_yx*b**2*c**2 + 20*Ex*c**4 + 5*Ey*nu_xy*b**2*c**2 - 10*Ey*b**4 - 28*G*nu_xy*nu_yx*b**2*c**2 + 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (Ex*nu_yx*c**2/2 - Ey*b**2/2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (5*Ex*c**2 - G*nu_xy*nu_yx*b**2 + G*b**2)/(5*b**2*c*(nu_xy*nu_yx - 1)), (-5*Ex*nu_yx*b**2*c**2 + 10*Ex*c**4 - 5*Ey*nu_xy*b**2*c**2 + 10*Ey*b**4 + 28*G*nu_xy*nu_yx*b**2*c**2 - 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (-Ey*b**2/2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (Ex*c**2/2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), (5*Ex*nu_yx*b**2*c**2 - 10*Ex*c**4 + 5*Ey*nu_xy*b**2*c**2 + 20*Ey*b**4 - 28*G*nu_xy*nu_yx*b**2*c**2 + 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (-5*Ey*b**2 + G*nu_xy*nu_yx*c**2 - G*c**2)/(5*b*c**2*(nu_xy*nu_yx - 1)), (Ex*c**2/2 - Ey*nu_xy*b**2/2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1))],
                           [(-Ey*nu_xy*c**2/2 - Ey*b**2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), 4*(-5*Ey*b**2 + 2*G*nu_xy*nu_yx*c**2 - 2*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), Ey*nu_xy/(nu_xy*nu_yx - 1), (Ey*nu_xy*c**2/2 - Ey*b**2/2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), 2*(-5*Ey*b**2 - 4*G*nu_xy*nu_yx*c**2 + 4*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0, (Ey*b**2/2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (-5*Ey*b**2 + 2*G*nu_xy*nu_yx*c**2 - 2*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0, (5*Ey*b**2 - G*nu_xy*nu_yx*c**2 + G*c**2)/(5*b*c**2*(nu_xy*nu_yx - 1)), 2*(-5*Ey*b**2 - G*nu_xy*nu_yx*c**2 + G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0],
                           [(Ex*nu_yx*b**2/2 + Ex*c**2 - G*nu_xy*nu_yx*b**2/5 + G*b**2/5)/(b**2*c*(nu_xy*nu_yx - 1)), Ex*nu_yx/(nu_xy*nu_yx - 1), 4*(-5*Ex*c**2 + 2*G*nu_xy*nu_yx*b**2 - 2*G*b**2)/(15*b*c*(nu_xy*nu_yx - 1)), (-5*Ex*c**2 + G*nu_xy*nu_yx*b**2 - G*b**2)/(5*b**2*c*(nu_xy*nu_yx - 1)), 0, 2*(-5*Ex*c**2 - G*nu_xy*nu_yx*b**2 + G*b**2)/(15*b*c*(nu_xy*nu_yx - 1)), -(Ex*c**2/2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), 0, (-5*Ex*c**2 + 2*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1)), (-Ex*nu_yx*b**2/2 + Ex*c**2/2 + G*nu_xy*nu_yx*b**2/5 - G*b**2/5)/(b**2*c*(nu_xy*nu_yx - 1)), 0, -(10*Ex*c**2 + 8*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1))],
                           [(5*Ex*nu_yx*b**2*c**2 + 20*Ex*c**4 + 5*Ey*nu_xy*b**2*c**2 - 10*Ey*b**4 - 28*G*nu_xy*nu_yx*b**2*c**2 + 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (Ex*nu_yx*c**2/2 - Ey*b**2/2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (-5*Ex*c**2 + G*nu_xy*nu_yx*b**2 - G*b**2)/(5*b**2*c*(nu_xy*nu_yx - 1)), (-Ex*nu_yx*b**2*c**2/4 - Ex*c**4 - Ey*nu_xy*b**2*c**2/4 - Ey*b**4 + 7*G*nu_xy*nu_yx*b**2*c**2/5 - 7*G*b**2*c**2/5)/(b**3*c**3*(nu_xy*nu_yx - 1)), (-Ex*nu_yx*c**2/2 - Ey*b**2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (-Ex*c**2 - Ey*nu_xy*b**2/2 + G*nu_xy*nu_yx*b**2/5 - G*b**2/5)/(b**2*c*(nu_xy*nu_yx - 1)), (5*Ex*nu_yx*b**2*c**2 - 10*Ex*c**4 + 5*Ey*nu_xy*b**2*c**2 + 20*Ey*b**4 - 28*G*nu_xy*nu_yx*b**2*c**2 + 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (-5*Ey*b**2 + G*nu_xy*nu_yx*c**2 - G*c**2)/(5*b*c**2*(nu_xy*nu_yx - 1)), (-Ex*c**2/2 + Ey*nu_xy*b**2/2 - G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), (-5*Ex*nu_yx*b**2*c**2 + 10*Ex*c**4 - 5*Ey*nu_xy*b**2*c**2 + 10*Ey*b**4 + 28*G*nu_xy*nu_yx*b**2*c**2 - 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (-Ey*b**2/2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), -(Ex*c**2/2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1))],
                           [(Ey*nu_xy*c**2/2 - Ey*b**2/2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), 2*(-5*Ey*b**2 - 4*G*nu_xy*nu_yx*c**2 + 4*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0, (-Ey*nu_xy*c**2/2 - Ey*b**2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), 4*(-5*Ey*b**2 + 2*G*nu_xy*nu_yx*c**2 - 2*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), -Ey*nu_xy/(nu_xy*nu_yx - 1), (5*Ey*b**2 - G*nu_xy*nu_yx*c**2 + G*c**2)/(5*b*c**2*(nu_xy*nu_yx - 1)), 2*(-5*Ey*b**2 - G*nu_xy*nu_yx*c**2 + G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0, (Ey*b**2/2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (-5*Ey*b**2 + 2*G*nu_xy*nu_yx*c**2 - 2*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0],
                           [(5*Ex*c**2 - G*nu_xy*nu_yx*b**2 + G*b**2)/(5*b**2*c*(nu_xy*nu_yx - 1)), 0, 2*(-5*Ex*c**2 - G*nu_xy*nu_yx*b**2 + G*b**2)/(15*b*c*(nu_xy*nu_yx - 1)), (-Ex*nu_yx*b**2/2 - Ex*c**2 + G*nu_xy*nu_yx*b**2/5 - G*b**2/5)/(b**2*c*(nu_xy*nu_yx - 1)), -Ex*nu_yx/(nu_xy*nu_yx - 1), 4*(-5*Ex*c**2 + 2*G*nu_xy*nu_yx*b**2 - 2*G*b**2)/(15*b*c*(nu_xy*nu_yx - 1)), (Ex*nu_yx*b**2/2 - Ex*c**2/2 - G*nu_xy*nu_yx*b**2/5 + G*b**2/5)/(b**2*c*(nu_xy*nu_yx - 1)), 0, -(10*Ex*c**2 + 8*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1)), (Ex*c**2/2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), 0, (-5*Ex*c**2 + 2*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1))],
                           [(-5*Ex*nu_yx*b**2*c**2 + 10*Ex*c**4 - 5*Ey*nu_xy*b**2*c**2 + 10*Ey*b**4 + 28*G*nu_xy*nu_yx*b**2*c**2 - 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (Ey*b**2/2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), -(Ex*c**2/2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), (5*Ex*nu_yx*b**2*c**2 - 10*Ex*c**4 + 5*Ey*nu_xy*b**2*c**2 + 20*Ey*b**4 - 28*G*nu_xy*nu_yx*b**2*c**2 + 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (5*Ey*b**2 - G*nu_xy*nu_yx*c**2 + G*c**2)/(5*b*c**2*(nu_xy*nu_yx - 1)), (-5*Ex*c**2 - 25*Ey*nu_xy*b**2 + 2*b**2*(15*Ey*nu_xy - G*nu_xy*nu_yx + G))/(10*b**2*c*(nu_xy*nu_yx - 1)), (-Ex*nu_yx*b**2*c**2/4 - Ex*c**4 - Ey*nu_xy*b**2*c**2/4 - Ey*b**4 + 7*G*nu_xy*nu_yx*b**2*c**2/5 - 7*G*b**2*c**2/5)/(b**3*c**3*(nu_xy*nu_yx - 1)), (Ex*nu_yx*c**2/2 + Ey*b**2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (-Ex*c**2 - Ey*nu_xy*b**2/2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), (5*Ex*nu_yx*b**2*c**2 + 20*Ex*c**4 + 5*Ey*nu_xy*b**2*c**2 - 10*Ey*b**4 - 28*G*nu_xy*nu_yx*b**2*c**2 + 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (-Ex*nu_yx*c**2/2 + Ey*b**2/2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (-Ex*c**2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1))],
                           [(-Ey*b**2/2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (-5*Ey*b**2 + 2*G*nu_xy*nu_yx*c**2 - 2*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0, (-5*Ey*b**2 + G*nu_xy*nu_yx*c**2 - G*c**2)/(5*b*c**2*(nu_xy*nu_yx - 1)), 2*(-5*Ey*b**2 - G*nu_xy*nu_yx*c**2 + G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0, (Ey*nu_xy*c**2/2 + Ey*b**2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), 4*(-5*Ey*b**2 + 2*G*nu_xy*nu_yx*c**2 - 2*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), Ey*nu_xy/(nu_xy*nu_yx - 1), (-Ey*nu_xy*c**2/2 + Ey*b**2/2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), 2*(-5*Ey*b**2 - 4*G*nu_xy*nu_yx*c**2 + 4*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0],
                           [(Ex*c**2/2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), 0, (-5*Ex*c**2 + 2*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1)), (Ex*nu_yx*b**2/2 - Ex*c**2/2 - G*nu_xy*nu_yx*b**2/5 + G*b**2/5)/(b**2*c*(nu_xy*nu_yx - 1)), 0, -(10*Ex*c**2 + 8*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1)), (-Ex*nu_yx*b**2/2 - Ex*c**2 + G*nu_xy*nu_yx*b**2/5 - G*b**2/5)/(b**2*c*(nu_xy*nu_yx - 1)), Ex*nu_yx/(nu_xy*nu_yx - 1), 4*(-5*Ex*c**2 + 2*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1)), (Ex*c**2 - G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), 0, -(10*Ex*c**2 + 2*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1))],
                           [(5*Ex*nu_yx*b**2*c**2 - 10*Ex*c**4 + 5*Ey*nu_xy*b**2*c**2 + 20*Ey*b**4 - 28*G*nu_xy*nu_yx*b**2*c**2 + 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (5*Ey*b**2 - G*nu_xy*nu_yx*c**2 + G*c**2)/(5*b*c**2*(nu_xy*nu_yx - 1)), (5*Ex*c**2 + 25*Ey*nu_xy*b**2 - 2*b**2*(15*Ey*nu_xy - G*nu_xy*nu_yx + G))/(10*b**2*c*(nu_xy*nu_yx - 1)), (-5*Ex*nu_yx*b**2*c**2 + 10*Ex*c**4 - 5*Ey*nu_xy*b**2*c**2 + 10*Ey*b**4 + 28*G*nu_xy*nu_yx*b**2*c**2 - 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (Ey*b**2/2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (Ex*c**2/2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), (5*Ex*nu_yx*b**2*c**2 + 20*Ex*c**4 + 5*Ey*nu_xy*b**2*c**2 - 10*Ey*b**4 - 28*G*nu_xy*nu_yx*b**2*c**2 + 28*G*b**2*c**2)/(20*b**3*c**3*(nu_xy*nu_yx - 1)), (-Ex*nu_yx*c**2/2 + Ey*b**2/2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (Ex*c**2 - G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), (-Ex*nu_yx*b**2*c**2/4 - Ex*c**4 - Ey*nu_xy*b**2*c**2/4 - Ey*b**4 + 7*G*nu_xy*nu_yx*b**2*c**2/5 - 7*G*b**2*c**2/5)/(b**3*c**3*(nu_xy*nu_yx - 1)), (Ex*nu_yx*c**2/2 + Ey*b**2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (Ex*c**2 + Ey*nu_xy*b**2/2 - G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1))],
                           [(-5*Ey*b**2 + G*nu_xy*nu_yx*c**2 - G*c**2)/(5*b*c**2*(nu_xy*nu_yx - 1)), 2*(-5*Ey*b**2 - G*nu_xy*nu_yx*c**2 + G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0, (-Ey*b**2/2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), (-5*Ey*b**2 + 2*G*nu_xy*nu_yx*c**2 - 2*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0, (-Ey*nu_xy*c**2/2 + Ey*b**2/2 + G*nu_xy*nu_yx*c**2/5 - G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), 2*(-5*Ey*b**2 - 4*G*nu_xy*nu_yx*c**2 + 4*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), 0, (Ey*nu_xy*c**2/2 + Ey*b**2 - G*nu_xy*nu_yx*c**2/5 + G*c**2/5)/(b*c**2*(nu_xy*nu_yx - 1)), 4*(-5*Ey*b**2 + 2*G*nu_xy*nu_yx*c**2 - 2*G*c**2)/(15*b*c*(nu_xy*nu_yx - 1)), -Ey*nu_xy/(nu_xy*nu_yx - 1)],
                           [(-Ex*nu_yx*b**2/2 + Ex*c**2/2 + G*nu_xy*nu_yx*b**2/5 - G*b**2/5)/(b**2*c*(nu_xy*nu_yx - 1)), 0, -(10*Ex*c**2 + 8*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1)), -(Ex*c**2/2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), 0, (-5*Ex*c**2 + 2*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1)), (-Ex*c**2 + G*b**2*(nu_xy*nu_yx - 1)/5)/(b**2*c*(nu_xy*nu_yx - 1)), 0, -(10*Ex*c**2 + 2*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1)), (Ex*nu_yx*b**2/2 + Ex*c**2 - G*nu_xy*nu_yx*b**2/5 + G*b**2/5)/(b**2*c*(nu_xy*nu_yx - 1)), -Ex*nu_yx/(nu_xy*nu_yx - 1), 4*(-5*Ex*c**2 + 2*G*b**2*(nu_xy*nu_yx - 1))/(15*b*c*(nu_xy*nu_yx - 1))]])
        
        # Calculate the stiffness of a weak spring for the drilling degree of freedom (rotation
        # about local z). We'll set the weak spring to be 1000 times weaker than any of the other
        # rotational stiffnesses in the matrix.
        k_rz = min(abs(k[1, 1]), abs(k[2, 2]), abs(k[4, 4]), abs(k[5, 5]),
                   abs(k[7, 7]), abs(k[8, 8]), abs(k[10, 10]), abs(k[11, 11])
                   )/1000

        # The matrix currently only holds terms related to bending action. We need to expand it to
        # with placeholders for all the degrees of freedom so it can be directly added to the
        # membrane stiffness matrix later on.

        # Initialize the expanded stiffness matrix to all zeros
        k_exp = zeros((24, 24))

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
        
        # Return the local stiffness matrix
        return k_exp

    def f(self, combo_name='Combo 1'):
        """
        Returns the plate's local end force vector
        """
        
        # Calculate and return the plate's local end force vector
        return add(matmul(self.k(), self.d(combo_name)), self.fer(combo_name))

    def fer(self, combo_name='Combo 1'):
        """
        Returns the rectangle's local fixed end reaction vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the load vector for.
        """
        
        # Initialize the fixed end reaction vector
        fer = zeros((12, 1))

        # Get the requested load combination
        combo = self.model.LoadCombos[combo_name]

        # Initialize the element's surface pressure to zero
        p = 0
        
        # Loop through each load case and factor in the load combination 
        for case, factor in combo.factors.items():

            # Sum the pressures
            for pressure in self.pressures:

                # Check if the current pressure corresponds to the current load case
                if pressure[1] == case:

                    # Sum the pressures multiplied by their load factors
                    p += factor*pressure[0]
        
        b = self.width()/2
        c = self.height()/2
        
        fer = -4*p*c*b*array([[1/4], [c/12], [-b/12], [1/4], [c/12], [b/12], [1/4], [-c/12], [b/12], [1/4], [-c/12], [-b/12]])

        # At this point `fer` only contains terms for the degrees of freedom
        # associated with membrane action. Expand `fer` to include zero terms for
        # the degrees of freedom related to bending action. This will allow
        # the bending and membrane vectors to be summed directly
        # later on. `numpy` has an `insert` function that can be used to
        # insert rows and columns of zeros one at a time, but it is very slow
        # as it makes a temporary copy of the vector term by term each time
        # it's called. The algorithm used here accomplishes the same thing
        # much faster. Terms are copied only once.

        # Initialize the expanded vector to all zeros
        fer_exp = zeros((24, 1))

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
        
    def d(self, combo_name='Combo 1'):
       """
       Returns the plate's local displacement vector
       """

       # Calculate and return the local displacement vector
       return matmul(self.T(), self.D(combo_name))

    def F(self, combo_name='Combo 1'):
        """
        Returns the plate's global nodal force vector
        """

        # Calculate and return the global force vector
        return matmul(inv(self.T()), self.f(combo_name))

    def D(self, combo_name='Combo 1'):
        """
        Returns the plate's global displacement vector for the given load combination.
        """
        
        # Initialize the displacement vector
        D = zeros((24, 1))
        
        # Read in the global displacements from the nodes
        D.itemset((0, 0), self.i_node.DX[combo_name])
        D.itemset((1, 0), self.i_node.DY[combo_name])
        D.itemset((2, 0), self.i_node.DZ[combo_name])
        D.itemset((3, 0), self.i_node.RX[combo_name])
        D.itemset((4, 0), self.i_node.RY[combo_name])
        D.itemset((5, 0), self.i_node.RZ[combo_name])

        D.itemset((6, 0), self.j_node.DX[combo_name])
        D.itemset((7, 0), self.j_node.DY[combo_name])
        D.itemset((8, 0), self.j_node.DZ[combo_name])
        D.itemset((9, 0), self.j_node.RX[combo_name])
        D.itemset((10, 0), self.j_node.RY[combo_name])
        D.itemset((11, 0), self.j_node.RZ[combo_name])

        D.itemset((12, 0), self.m_node.DX[combo_name])
        D.itemset((13, 0), self.m_node.DY[combo_name])
        D.itemset((14, 0), self.m_node.DZ[combo_name])
        D.itemset((15, 0), self.m_node.RX[combo_name])
        D.itemset((16, 0), self.m_node.RY[combo_name])
        D.itemset((17, 0), self.m_node.RZ[combo_name])

        D.itemset((18, 0), self.n_node.DX[combo_name])
        D.itemset((19, 0), self.n_node.DY[combo_name])
        D.itemset((20, 0), self.n_node.DZ[combo_name])
        D.itemset((21, 0), self.n_node.RX[combo_name])
        D.itemset((22, 0), self.n_node.RY[combo_name])
        D.itemset((23, 0), self.n_node.RZ[combo_name])
        
        # Return the global displacement vector
        return D
 
    def T(self):
        """
        Returns the plate's transformation matrix
        """

        # Calculate the direction cosines for the local x-axis
        # The local x-axis will run from the i-node to the j-node
        xi = self.i_node.X
        xj = self.j_node.X
        yi = self.i_node.Y
        yj = self.j_node.Y
        zi = self.i_node.Z
        zj = self.j_node.Z
        x = [(xj - xi), (yj - yi), (zj - zi)]
        x = x/norm(x)
        
        # The local y-axis will be in the plane of the plate
        # Find a vector in the plate's local xy plane
        xn = self.n_node.X
        yn = self.n_node.Y
        zn = self.n_node.Z
        xy = [xn - xi, yn - yi, zn - zi]

        # Find a vector perpendicular to the plate surface to get the orientation of the local z-axis
        z = cross(x, xy)
        
        # Divide the vector by its magnitude to produce a unit z-vector of direction cosines
        z = z/norm(z)

        # Calculate the local y-axis as a vector perpendicular to the local z and x-axes
        y = cross(z, x)
        
        # Divide the z-vector by its magnitude to produce a unit vector of direction cosines
        y = y/norm(y)

        # Create the direction cosines matrix
        dirCos = array([x, y, z])
        
        # Build the transformation matrix
        transMatrix = zeros((24, 24))
        transMatrix[0:3, 0:3] = dirCos
        transMatrix[3:6, 3:6] = dirCos
        transMatrix[6:9, 6:9] = dirCos
        transMatrix[9:12, 9:12] = dirCos
        transMatrix[12:15, 12:15] = dirCos
        transMatrix[15:18, 15:18] = dirCos
        transMatrix[18:21, 18:21] = dirCos
        transMatrix[21:24, 21:24] = dirCos
        
        return transMatrix

    def K(self):
        """
        Returns the plate's global stiffness matrix
        """

        # Calculate and return the stiffness matrix in global coordinates
        return matmul(matmul(inv(self.T()), self.k()), self.T())

    def FER(self, combo_name='Combo 1'):
        """
        Returns the global fixed end reaction vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to calculate the fixed end
            reaction vector for (not the load combination itself).
        """
        
        # Calculate and return the fixed end reaction vector
        return matmul(inv(self.T()), self.fer(combo_name))

    def _C(self):
        """
        Returns the plate's displacement coefficient matrix [C]
        """

        # Find the local x and y coordinates at each node
        xi = 0
        yi = 0
        xj = self.width()
        yj = 0
        xm = xj
        ym = self.height()
        xn = 0
        yn = ym

        # Calculate the [C] coefficient matrix
        C = array([[1, xi, yi, xi**2, xi*yi, yi**2, xi**3, xi**2*yi, xi*yi**2, yi**3, xi**3*yi, xi*yi**3],
                    [0, 0, 1, 0, xi, 2*yi, 0, xi**2, 2*xi*yi, 3*yi**2, xi**3, 3*xi*yi**2],
                    [0, -1, 0, -2*xi, -yi, 0, -3*xi**2, -2*xi*yi, -yi**2, 0, -3*xi**2*yi, -yi**3],
                    
                    [1, xj, yj, xj**2, xj*yj, yj**2, xj**3, xj**2*yj, xj*yj**2, yj**3, xj**3*yj, xj*yj**3],
                    [0, 0, 1, 0, xj, 2*yj, 0, xj**2, 2*xj*yj, 3*yj**2, xj**3, 3*xj*yj**2],
                    [0, -1, 0, -2*xj, -yj, 0, -3*xj**2, -2*xj*yj, -yj**2, 0, -3*xj**2*yj, -yj**3],

                    [1, xm, ym, xm**2, xm*ym, ym**2, xm**3, xm**2*ym, xm*ym**2, ym**3, xm**3*ym, xm*ym**3],
                    [0, 0, 1, 0, xm, 2*ym, 0, xm**2, 2*xm*ym, 3*ym**2, xm**3, 3*xm*ym**2],
                    [0, -1, 0, -2*xm, -ym, 0, -3*xm**2, -2*xm*ym, -ym**2, 0, -3*xm**2*ym, -ym**3],

                    [1, xn, yn, xn**2, xn*yn, yn**2, xn**3, xn**2*yn, xn*yn**2, yn**3, xn**3*yn, xn*yn**3],
                    [0, 0, 1, 0, xn, 2*yn, 0, xn**2, 2*xn*yn, 3*yn**2, xn**3, 3*xn*yn**2],
                    [0, -1, 0, -2*xn, -yn, 0, -3*xn**2, -2*xn*yn, -yn**2, 0, -3*xn**2*yn, -yn**3]])
        
        # Return the coefficient matrix
        return C

    def _Q(self, x, y):
        """
        Calculates and returns the plate curvature coefficient matrix [Q] at a given point (x, y)
        in the plate's local system.
        """

        # Calculate the [Q] coefficient matrix
        Q =  array([[0, 0, 0, -2, 0, 0, -6*x, -2*y, 0, 0, -6*x*y, 0],
                    [0, 0, 0, 0, 0, -2, 0, 0, -2*x, -6*y, 0, -6*x*y],
                    [0, 0, 0, 0, -2, 0, 0, -4*x, -4*y, 0, -6*x**2, -6*y**2]])
 
        # Return the [Q] coefficient matrix
        return Q

    def _a(self, combo_name='Combo 1'):
        """
        Returns the vector of plate bending constants for the displacement function.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the vector for
        """

        # Get the plate's local displacement vector
        # Slice out terms not related to plate bending
        d = self.d(combo_name)[[2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22], :]

        # Return the plate bending constants
        return inv(self._C()) @ d

    def moment(self, x, y, local=True, combo_name='Combo 1'):
        """
        Returns the internal moments (Mx, My, and Mxy) at any point (x, y) in the plate's local
        coordinate system

        Parameters
        ----------
        x : number
            The x-coordinate in the plate's local coordinate system.
        y : number
            The y-coordinate in the plate's local coordinate system.
        combo_name : string
            The name of the load combination to evaluate. The default is 'Combo 1'.

        """
        
        # Calculate and return internal moments
        # A negative sign will be applied to change the sign convention to match that of
        # Pynite's quadrilateral elements.
        Mx, My, Mxy = -self.Db() @ self._Q(x, y) @ self._a(combo_name)

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
 
    def shear(self, x, y, local=True, combo_name='Combo 1'):
        """
        Returns the internal shears (Qx and Qy) at any point (x, y) in the plate's local
        coordinate system

        Parameters
        ----------
        x : number
            The x-coordinate in the plate's local coordinate system.
        y : number
            The y-coordinate in the plate's local coordinate system.
        combo_name : string
            The name of the load combination to evaluate. The default is 'Combo 1'.

        """

        # Store matrices into local variables for quicker access
        Db = self.Db()
        a = self._a(combo_name)

        # Calculate the derivatives of the plate moments needed to compute shears
        dMx_dx = (Db @ array([[0, 0, 0, 0, 0, 0, -6, 0, 0, 0, -6*y, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, -6*y],
                             [0, 0, 0, 0, 0, 0, 0, -4, 0, 0, -12*x, 0]]) @ a)[0]

        dMxy_dy = (Db @ array([[0, 0, 0, 0, 0, 0, 0, -2, 0, 0, -6*x, 0],
                              [0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 0, -6*x],
                              [0, 0, 0, 0, 0, 0, 0, 0, -4, 0, 0, -12*y]]) @ a)[2]
        
        dMy_dy = (Db @ array([[0, 0, 0, 0, 0, 0, 0, -2, 0, 0, -6*x, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 0, -6*x],
                             [0, 0, 0, 0, 0, 0, 0, 0, -4, 0, 0, -12*y]]) @ a)[1]

        dMxy_dx = (Db @ array([[0, 0, 0, 0, 0, 0, -6, 0, 0, 0, -6*y, 0],
                              [0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, -6*y],
                              [0, 0, 0, 0, 0, 0, 0, -4, 0, 0, -12*x, 0]]) @ a)[2]
        
        # Calculate internal shears
        Qx = (dMx_dx + dMxy_dy)[0]
        Qy = (dMy_dy + dMxy_dx)[0]

        # Return internal shears
        if local:

            return array([[Qx], 
                          [Qy],])
        
        else:

            # Get the direction cosines for the plate's local coordinate system
            dir_cos = self.T()[:3, :3]

            # Convert the plate shear stresses to global coordinates
            return np.matmul(np.linalg.inv(dir_cos), np.array([[Qx],
                                                               [Qy],
                                                               [0]]))

    def membrane(self, x, y, local=True, combo_name='Combo 1'):
        
        # Convert the (x, y) coordinates to (r, x) coordinates
        r = -1 + 2*x/self.width()
        s = -1 + 2*y/self.height()

        # Get the plate's local displacement vector
        # Slice out terms not related to membrane forces
        d = self.d(combo_name)[[0, 1, 6, 7, 12, 13, 18, 19], :]

        # Define the gauss point used for numerical integration
        gp = 1/3**0.5

        # Define extrapolated r and s points
        r_ex = r/gp
        s_ex = s/gp

        # Define the interpolation functions
        H = 1/4*array([(1 - r_ex)*(1 - s_ex), (1 + r_ex)*(1 - s_ex), (1 + r_ex)*(1 + s_ex), (1 - r_ex)*(1 + s_ex)])

        # Get the stress-strain matrix
        Dm = self.Dm()
        
        # Calculate the internal stresses [Sx, Sy, Txy] at each gauss point
        s1 = matmul(Dm, matmul(self.B_m(-gp, -gp), d))
        s2 = matmul(Dm, matmul(self.B_m(gp, -gp), d))
        s3 = matmul(Dm, matmul(self.B_m(gp, gp), d))
        s4 = matmul(Dm, matmul(self.B_m(-gp, gp), d))

        # Extrapolate to get the value at the requested location
        Sx = H[0]*s1[0] + H[1]*s2[0] + H[2]*s3[0] + H[3]*s4[0]
        Sy = H[0]*s1[1] + H[1]*s2[1] + H[2]*s3[1] + H[3]*s4[1]
        Txy = H[0]*s1[2] + H[1]*s2[2] + H[2]*s3[2] + H[3]*s4[2]

        if local:

            return array([Sx,
                          Sy,
                          Txy])

        else:

            # Get the direction cosines for the plate's local coordinate system
            dir_cos = self.T()[:3, :3]

            # Convert the plate membrane stresses to global coordinates
            return np.matmul(np.linalg.inv(dir_cos), np.array([Sx,
                                                               Sy,
                                                               [0]]))
