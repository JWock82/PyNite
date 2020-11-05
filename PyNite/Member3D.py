# %%
from numpy import zeros, matrix, transpose, add, subtract, matmul, insert, cross, divide
from numpy.linalg import inv
from math import isclose
from PyNite.BeamSegZ import BeamSegZ
from PyNite.BeamSegY import BeamSegY
import PyNite.FixedEndReactions
from PyNite.LoadCombo import LoadCombo

# %%
class Member3D():
    '''
    A class representing a 3D frame element in a finite element model.
    '''

    # '__plt' is used to store the 'pyplot' from matplotlib once it gets imported. Setting it to 'None' for now allows
    # us to defer importing it until it's actually needed.
    __plt = None

#%%
    def __init__(self, Name, iNode, jNode, E, G, Iy, Iz, J, A, auxNode=None, 
                 LoadCombos={'Combo 1':LoadCombo('Combo 1', 
                             factors={'Case 1':1.0})},
                 tension_only=False, comp_only=False):
        '''
        Initializes a new member.
        '''
        self.Name = Name    # A unique name for the member given by the user
        self.ID = None      # Unique index number for the member assigned by the program
        self.iNode = iNode  # The element's i-node
        self.jNode = jNode  # The element's j-node
        self.E = E  # The modulus of elasticity of the element
        self.G = G  # The shear modulus of the element
        self.Iy = Iy  # The y-axis moment of inertia
        self.Iz = Iz  # The z-axis moment of inertia
        self.J = J  # The torsional constant
        self.A = A  # The cross-sectional area
        self.auxNode = auxNode # Optional auxiliary node used to define the member's local z-axis
        self.PtLoads = []   # A list of point loads & moments applied to the element (Direction, P, x, case='Case 1') or (Direction, M, x, case='Case 1')
        self.DistLoads = [] # A list of linear distributed loads applied to the element (Direction, w1, w2, x1, x2, case='Case 1')
        self.SegmentsZ = [] # A list of mathematically continuous beam segments for z-bending
        self.SegmentsY = [] # A list of mathematically continuous beam segments for y-bending
        self.SegmentsX = [] # A list of mathematically continuous beam segments for torsion
        self.Releases = [False, False, False, False, False, False, False, False, False, False, False, False]
        self.LoadCombos = LoadCombos # The dictionary of load combinations in the model this member belongs to
        self.tension_only = tension_only # Indicates whether the member is tension-only
        self.comp_only = comp_only # Indicates whether the member is compression-only

        # Members need to track whether they are active or not for any given load combination.
        # They may become inactive for a load combination during a tension/compression-only
        # analysis. This dictionary will be used when the model is solved.
        self.active = {} # Key = load combo name, Value = True or False
        
        # The 'Member3D' object will store results for one load combination at a time. To reduce repetative calculations
        # the '__solved_combo' variable will be used to track whether the member needs to be resegmented before running
        # calculations for any given load combination.
        self.__solved_combo = None # The current solved load combination

#%%
    def L(self):
        '''
        Returns the length of the member.
        '''

        # Get the i-node and the j-node for the member
        iNode = self.iNode
        jNode = self.jNode

        # Return the distance between the two nodes
        return ((jNode.X-iNode.X)**2+(jNode.Y-iNode.Y)**2+(jNode.Z-iNode.Z)**2)**0.5

#%%
    def AuxList(self):
        '''
        Builds lists of unreleased and released degree of freedom indices for the member.

        Returns
        -------
        R1_indices : list
            A list of the indices for the unreleased DOFs
        R2_indices : list
            A list of the indices for the released DOFs
        '''

        R1_indices = []
        R2_indices = []
        for i in range(12):
            if self.Releases[i] == False:
                R1_indices.append(i)
            else:
                R2_indices.append(i)
        
        return R1_indices, R2_indices

#%%
    def k(self):
        '''
        Returns the condensed (and expanded) local stiffness matrix for the member.
        '''

        # Partition the local stiffness matrix as 4 submatrices in
        # preparation for static condensation
        k11, k12, k21, k22 = self.__Partition(self.__k_Unc())
               
        # Calculate the condensed local stiffness matrix
        k_Condensed = subtract(k11, matmul(matmul(k12, inv(k22)), k21))
        
        # Expand the condensed local stiffness matrix
        i=0
        for DOF in self.Releases:
            
            if DOF == True:
                k_Condensed = insert(k_Condensed, i, 0, axis = 0)
                k_Condensed = insert(k_Condensed, i, 0, axis = 1)
                
            i += 1

        # Return the local stiffness matrix, with end releases applied
        return k_Condensed

#%%
    def __k_Unc(self):
        '''
        Returns the uncondensed local stiffness matrix for the member.
        '''

        # Get the properties needed to form the local stiffness matrix
        E = self.E
        G = self.G
        Iy = self.Iy
        Iz = self.Iz
        J = self.J
        A = self.A
        L = self.L()
        
        # Create the uncondensed local stiffness matrix
        k = matrix([[A*E/L,  0,             0,             0,      0,            0,            -A*E/L, 0,             0,             0,      0,            0],
                    [0,      12*E*Iz/L**3,  0,             0,      0,            6*E*Iz/L**2,  0,      -12*E*Iz/L**3, 0,             0,      0,            6*E*Iz/L**2],
                    [0,      0,             12*E*Iy/L**3,  0,      -6*E*Iy/L**2, 0,            0,      0,             -12*E*Iy/L**3, 0,      -6*E*Iy/L**2, 0],
                    [0,      0,             0,             G*J/L,  0,            0,            0,      0,             0,             -G*J/L, 0,            0],
                    [0,      0,             -6*E*Iy/L**2,  0,      4*E*Iy/L,     0,            0,      0,             6*E*Iy/L**2,   0,      2*E*Iy/L,     0],
                    [0,      6*E*Iz/L**2,   0,             0,      0,            4*E*Iz/L,     0,      -6*E*Iz/L**2,  0,             0,      0,            2*E*Iz/L],
                    [-A*E/L, 0,             0,             0,      0,            0,            A*E/L,  0,             0,             0,      0,            0],
                    [0,      -12*E*Iz/L**3, 0,             0,      0,            -6*E*Iz/L**2, 0,      12*E*Iz/L**3,  0,             0,      0,            -6*E*Iz/L**2],
                    [0,      0,             -12*E*Iy/L**3, 0,      6*E*Iy/L**2,  0,            0,      0,             12*E*Iy/L**3,  0,      6*E*Iy/L**2,  0],
                    [0,      0,             0,             -G*J/L, 0,            0,            0,      0,             0,             G*J/L,  0,            0],
                    [0,      0,             -6*E*Iy/L**2,  0,      2*E*Iy/L,     0,            0,      0,             6*E*Iy/L**2,   0,      4*E*Iy/L,     0],
                    [0,      6*E*Iz/L**2,   0,             0,      0,            2*E*Iz/L,     0,      -6*E*Iz/L**2,  0,             0,      0,            4*E*Iz/L]])
        
        # Return the uncondensed local stiffness matrix
        return k

#%%
    def kg(self, P=0):
        '''
        Returns the condensed (expanded) local geometric stiffness matrix for the member.

        Parameters
        ----------
        P : number, optional
            The axial force acting on the member (compression = +, tension = -)
        '''

        # Get the properties needed to form the local geometric stiffness matrix
        Ip = self.Iy + self.Iz
        A = self.A
        L = self.L()
        
        # Create the uncondensed local geometric stiffness matrix
        kg = matrix([[0, 0,    0,     0,     0,         0,         0, 0,     0,    0,     0,         0],
                     [0, 6/5,  0,     0,     0,         L/10,      0, -6/5,  0,    0,     0,         L/10],
                     [0, 0,    6/5,   0,     -L/10,     0,         0, 0,     -6/5, 0,     -L/10,     0],
                     [0, 0,    0,     Ip/A,  0,         0,         0, 0,     0,    -Ip/A, 0,         0],
                     [0, 0,    -L/10, 0,     2*L**2/15, 0,         0, 0,     L/10, 0,     -L**2/30,  0],
                     [0, L/10, 0,     0,     0,         2*L**2/15, 0, -L/10, 0,    0,     0,         -L**2/30],
                     [0, 0,    0,     0,     0,         0,         0, 0,     0,    0,     0,         0],
                     [0, -6/5, 0,     0,     0,         -L/10,     0, 6/5,   0,    0,     0,         -L/10],
                     [0, 0,    -6/5,  0,     L/10,      0,         0, 0,     6/5,  0,     L/10,      0],
                     [0, 0,    0,     -Ip/A, 0,         0,         0, 0,     0,    Ip/A,  0,         0],
                     [0, 0,    -L/10, 0,     -L**2/30,  0,         0, 0,     L/10, 0,     2*L**2/15, 0],
                     [0, L/10, 0,     0,     0,         -L**2/30,  0, -L/10, 0,    0,     0,         2*L**2/15]])
        
        kg = kg*P/L

        # Partition the geometric stiffness matrix as 4 submatrices in
        # preparation for static condensation
        kg11, kg12, kg21, kg22 = self.__Partition(kg)
               
        # Calculate the condensed local geometric stiffness matrix
        # Note that a matrix of zeros cannot be inverted, so if P is 0 an error will occur
        if isclose(P, 0.0):
            kg_Condensed = zeros(kg11.shape)
        else:
            kg_Condensed = subtract(kg11, matmul(matmul(kg12, inv(kg22)), kg21))
        
        # Expand the condensed local geometric stiffness matrix
        i=0
        for DOF in self.Releases:
            
            if DOF == True:
                kg_Condensed = insert(kg_Condensed, i, 0, axis = 0)
                kg_Condensed = insert(kg_Condensed, i, 0, axis = 1)
                
            i += 1

        # Return the local geomtric stiffness matrix, with end releases applied
        return kg_Condensed
    
#%%
    def fer(self, combo_name='Combo 1'):
        '''
        Returns the condensed (and expanded) local fixed end reaction vector for the member for the given load combination.

        Parameters
        ----------
        combo : LoadCombo
            The load combination to construct the fixed end reaction vector for.
        '''
        
        # Get the lists of unreleased and released degree of freedom indices
        R1_indices, R2_indices = self.AuxList()

        # Partition the local stiffness matrix and local fixed end reaction vector
        k11, k12, k21, k22 = self.__Partition(self.__k_Unc())
        fer1, fer2 = self.__Partition(self.__fer_Unc(combo_name))
        
        # Calculate the condensed fixed end reaction vector
        ferCondensed = subtract(fer1, matmul(matmul(k12, inv(k22)), fer2))
        
        # Expand the condensed fixed end reaction vector
        i=0
        for DOF in self.Releases:
            
            if DOF == True:
                ferCondensed = insert(ferCondensed, i, 0, axis = 0)
                
            i += 1
        
        # Return the fixed end reaction vector        
        return ferCondensed
    
#%%
    def __fer_Unc(self, combo_name='Combo 1'):
        '''
        Returns the member's local fixed end reaction vector, ignoring the effects of end releases.
        Needed to apply the slope-deflection equation properly.
        '''
        
        # Initialize the fixed end reaction vector
        fer = zeros((12,1))

        # Get the requested load combination
        combo = self.LoadCombos[combo_name]

        # Loop through each load case and factor in the load combination
        for case, factor in combo.factors.items():
        
            # Sum the fixed end reactions for the point loads & moments
            for ptLoad in self.PtLoads:

                # Check if the current point load corresponds to the current load case
                if ptLoad[3] == case:

                    if ptLoad[0] == 'Fx':
                        fer = add(fer, PyNite.FixedEndReactions.FER_AxialPtLoad(factor*ptLoad[1], ptLoad[2], self.L()))
                    elif ptLoad[0] == 'Fy':
                        fer = add(fer, PyNite.FixedEndReactions.FER_PtLoad(factor*ptLoad[1], ptLoad[2], self.L(), 'Fy'))
                    elif ptLoad[0] == 'Fz':
                        fer = add(fer, PyNite.FixedEndReactions.FER_PtLoad(factor*ptLoad[1], ptLoad[2], self.L(), 'Fz'))
                    elif ptLoad[0] == 'Mx':
                        fer = add(fer, PyNite.FixedEndReactions.FER_Torque(factor*ptLoad[1], ptLoad[2], self.L()))
                    elif ptLoad[0] == 'My':
                        fer = add(fer, PyNite.FixedEndReactions.FER_Moment(factor*ptLoad[1], ptLoad[2], self.L(), 'My'))
                    elif ptLoad[0] == 'Mz':     
                        fer = add(fer, PyNite.FixedEndReactions.FER_Moment(factor*ptLoad[1], ptLoad[2], self.L(), 'Mz'))
                
            # Sum the fixed end reactions for the distributed loads
            for distLoad in self.DistLoads:

                if distLoad[5] == case:

                    if distLoad[0] == 'Fx':
                        fer = add(fer, PyNite.FixedEndReactions.FER_AxialLinLoad(factor*distLoad[1], factor*distLoad[2], distLoad[3], distLoad[4], self.L()))
                    else:
                        fer = add(fer, PyNite.FixedEndReactions.FER_LinLoad(factor*distLoad[1], factor*distLoad[2], distLoad[3], distLoad[4], self.L(), distLoad[0]))
        
        # Return the fixed end reaction vector, uncondensed
        return fer

#%%
    def __Partition(self, unp_matrix):
        '''
        Partitions a matrix into sub-matrices based on unreleased and released degree of freedom indices.
        '''

        # Create auxiliary lists of released/unreleased DOFs
        R1_indices, R2_indices = self.AuxList()

        # Partition the matrix by slicing
        if unp_matrix.shape[1] == 1:
            m1 = unp_matrix[R1_indices, :]
            m2 = unp_matrix[R2_indices, :]
            return m1, m2
        else:
            m11 = unp_matrix[R1_indices, :][:, R1_indices]
            m12 = unp_matrix[R1_indices, :][:, R2_indices]
            m21 = unp_matrix[R2_indices, :][:, R1_indices]
            m22 = unp_matrix[R2_indices, :][:, R2_indices]
            return  m11, m12, m21, m22

#%%   
    def f(self, combo_name='Combo 1'):
        '''
        Returns the member's local end force vector for the given load combination.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to calculate the local end force vector for (not the load combination itself).
        '''
        
        # Calculate and return the member's local end force vector
        return add(matmul(self.k(), self.d(combo_name)), self.fer(combo_name))

#%%
    def d(self, combo_name='Combo 1'):
        '''
        Returns the member's local displacement vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to construct the displacement vector for (not the load combination itself).
        '''
        
        # Calculate and return the local displacement vector
        return matmul(self.T(), self.D(combo_name))
        
#%%  
    # Transformation matrix
    def T(self):
        '''
        Returns the transformation matrix for the member.
        '''

        x1 = self.iNode.X
        x2 = self.jNode.X
        y1 = self.iNode.Y
        y2 = self.jNode.Y
        z1 = self.iNode.Z
        z2 = self.jNode.Z
        L = self.L()
        
        # Calculate the direction cosines for the local x-axis
        x = [(x2-x1)/L, (y2-y1)/L, (z2-z1)/L]
        
        # Calculate the direction cosines for the local z-axis based on the auxiliary node
        if self.auxNode != None:
            
            xa = self.auxNode.X
            ya = self.auxNode.Y
            za = self.auxNode.Z
            
            # Define a vector in the local xz plane using the auxiliary point 
            z = [xa-x1, ya-y1, za-z1]
            
            # Find the direction cosines for the local y-axis
            y = cross(z, x)
            y = divide(y, (y[0]**2 + y[1]**2 + y[2]**2)**0.5)
            
            # Ensure the z-axis is perpendicular to the x-axis.
            # If the vector from the i-node to the auxiliary node is not perpendicular to the member, this will ensure the local coordinate system is orthogonal
            z = cross(x, y)
            
            # Turn the z-vector into a unit vector of direction cosines
            z = divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)
        
        # If no auxiliary node is specified the program will determine the member's local z-axis automatically
        else:
            
            # Calculate the remaining direction cosines. The local z-axis will be kept parallel to the global XZ plane in all cases
            # Vertical members
            if isclose(x1, x2) and isclose(z1, z2):
                
                # For vertical members, keep the local y-axis in the XY plane to make 2D problems easier to solve in the XY plane
                if y2 > y1:
                    y = [-1, 0, 0]
                    z = [0, 0, 1]
                else:
                    y = [1, 0, 0]
                    z = [0, 0, 1]

            # Horizontal members
            elif isclose(y1, y2):
            
                # Find a vector in the direction of the local z-axis by taking the cross-product
                # of the local x-axis and the local y-axis. This vector will be perpendicular to
                # both the local x-axis and the local y-axis.
                y = [0, 1, 0]
                z = cross(x, y)

                # Divide the z-vector by its magnitude to produce a unit vector of direction cosines
                z = divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)

            # Members neither vertical or horizontal
            else:

                # Find the projection of x on the global XZ plane
                proj = [x2-x1, 0, z2-z1]

                # Find a vector in the direction of the local z-axis by taking the cross-product
                # of the local x-axis and its projection on a plane parallel to the XZ plane. This
                # produces a vector perpendicular to both the local x-axis and its projection. This
                # vector will always be horizontal since it's parallel to the XZ plane. The order
                # in which the vectors are 'crossed' has been selected to ensure the y-axis always
                # has an upward component (i.e. the top of the beam is always on top).
                if y2 > y1:
                    z = cross(proj, x)
                else:
                    z = cross(x, proj)

                # Divide the z-vector by its magnitude to produce a unit vector of direction cosines
                z = divide(z, (z[0]**2 + z[1]**2 + z[2]**2)**0.5)
                
                # Find the direction cosines for the local y-axis
                y = cross(z, x)
                y = divide(y, (y[0]**2 + y[1]**2 + y[2]**2)**0.5)

        # Create the direction cosines matrix
        dirCos = matrix([x, y, z])
      
        # Build the transformation matrix
        transMatrix = zeros((12, 12))
        transMatrix[0:3, 0:3] = dirCos
        transMatrix[3:6, 3:6] = dirCos
        transMatrix[6:9, 6:9] = dirCos
        transMatrix[9:12, 9:12] = dirCos
        
        return transMatrix

#%%
    # Member global stiffness matrix
    def K(self):
        """
        Compute the transpose.

        Args:
            self: (todo): write your description
        """
        
        # Calculate and return the stiffness matrix in global coordinates
        return matmul(matmul(transpose(self.T()), self.k()), self.T())

#%%
    # Member global geometric stiffness matrix
    def Kg(self, P=0):
        """
        Return the superposition matrix.

        Args:
            self: (todo): write your description
            P: (int): write your description
        """
        
        # Calculate and return the geometric stiffness matrix in global coordinates
        return matmul(matmul(transpose(self.T()), self.kg(P)), self.T())

#%%
    def F(self, combo_name='Combo 1'):
        '''
        Returns the member's global end force vector for the given load combination.
        '''
        
        # Calculate and return the global force vector
        return matmul(inv(self.T()), self.f(combo_name))
    
#%% 
    # Global fixed end reaction vector
    def FER(self, combo_name='Combo 1'):
        '''
        Returns the global fixed end reaction vector

        Parameters
        ----------
        combo_name : string
            The name of the load combination to calculate the fixed end reaction vector for (not the load combination itself).
        '''
        
        # Calculate and return the fixed end reaction vector
        return matmul(inv(self.T()), self.fer(combo_name))

#%%
    def D(self, combo_name='Combo 1'):
        '''
        Returns the member's global displacement vector.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to construct the global
            displacement vector for (not the load combination itelf).
        '''
        
        # Initialize the displacement vector
        D = zeros((12, 1))
        
        # Read in the global displacements from the nodes
        # Apply axial displacements only if the member is active
        if self.active[combo_name] == True:
            D[0, 0] = self.iNode.DX[combo_name]
            D[6, 0] = self.jNode.DX[combo_name]

        # Apply the remaining displacements
        D[1, 0] = self.iNode.DY[combo_name]
        D[2, 0] = self.iNode.DZ[combo_name]
        D[3, 0] = self.iNode.RX[combo_name]
        D[4, 0] = self.iNode.RY[combo_name]
        D[5, 0] = self.iNode.RZ[combo_name]
        D[7, 0] = self.jNode.DY[combo_name]
        D[8, 0] = self.jNode.DZ[combo_name]
        D[9, 0] = self.jNode.RX[combo_name]
        D[10, 0] = self.jNode.RY[combo_name]
        D[11, 0] = self.jNode.RZ[combo_name]      

        # Return the global displacement vector
        return D

#%%
    def Shear(self, Direction, x, combo_name='Combo 1'):
        '''
        Returns the shear at a point along the member's length.
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the shear. Must be one of the following:
                'Fy' = Shear acting on the local y-axis.
                'Fz' = Shear acting on the local z-axis.
        x : number
            The location at which to find the shear.
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]

        # Check which direction is of interest
        if Direction == 'Fy':
            
            # Check which segment 'x' falls on
            for segment in self.SegmentsZ:
                if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                    return segment.Shear(x - segment.x1)
                
            if isclose(x, self.L()):  
                lastIndex = len(self.SegmentsZ) - 1
                return self.SegmentsZ[lastIndex].Shear(x - self.SegmentsZ[lastIndex].x1)
                
        elif Direction == 'Fz':
            
            for segment in self.SegmentsY:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return segment.Shear(x - segment.x1)
                
            if isclose(x, self.L()):
                
                lastIndex = len(self.SegmentsY) - 1
                return self.SegmentsY[lastIndex].Shear(x - self.SegmentsY[lastIndex].x1)
            
#%%
    def MaxShear(self, Direction, combo_name='Combo 1'):
        '''
        Returns the maximum shear in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the maximum shear. Must be one of the following:
                'Fy' = Shear acting on the local y-axis
                'Fz' = Shear acting on the local z-axis
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        if Direction == 'Fy':
            
            Vmax = self.SegmentsZ[0].Shear(0)

            for segment in self.SegmentsZ:
                
                if segment.MaxShear() > Vmax:
                    
                    Vmax = segment.MaxShear()
                    
        if Direction == 'Fz':
            
            Vmax = self.SegmentsY[0].Shear(0)

            for segment in self.SegmentsY:
                
                if segment.MaxShear() > Vmax:
                    
                    Vmax = segment.MaxShear()
        
        return Vmax
    
#%%
    def MinShear(self, Direction, combo_name='Combo 1'):
        '''
        Returns the minimum shear in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the minimum shear. Must be one of the following:
                'Fy' = Shear acting on the local y-axis
                'Fz' = Shear acting on the local z-axis
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]   
        
        if Direction == 'Fy':
            
            Vmin = self.SegmentsZ[0].Shear(0)

            for segment in self.SegmentsZ:
                
                if segment.MinShear() < Vmin:
                    
                    Vmin = segment.MinShear()
                    
        if Direction == 'Fz':
            
            Vmin = self.SegmentsY[0].Shear(0)

            for segment in self.SegmentsY:
                
                if segment.MinShear() < Vmin:
                    
                    Vmin = segment.MinShear()
        
        return Vmin
    
#%%
    def PlotShear(self, Direction, combo_name='Combo 1'):
        '''
        Plots the shear diagram for the member
        
        Parameters
        ----------
        Direction : string
            The direction to plot the shear for.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt

        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x = []
        V = []
        
        # Calculate the shear diagram
        for i in range(21):
            x.append(self.L()/20*i)
            V.append(self.Shear(Direction, self.L()/20*i, combo_name))

        Member3D.__plt.plot(x, V)
        Member3D.__plt.ylabel('Shear')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.Name + '\n' + combo_name)
        Member3D.__plt.show()    
        
#%%
    def Moment(self, Direction, x, combo_name='Combo 1'):
        '''
        Returns the moment at a point along the member's length
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the moment. Must be one of the following:
                'My' = Moment about the local y-axis.
                'Mz' = moment about the local z-axis.
        x : number
            The location at which to find the moment.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        # Check which axis is of interest
        if Direction == 'My':
            
            # Check which segment 'x' falls on
            for segment in self.SegmentsY:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return segment.Moment(x - segment.x1)
                
            if isclose(x, self.L()):
                
                lastIndex = len(self.SegmentsY) - 1
                return self.SegmentsY[lastIndex].Moment(x - self.SegmentsY[lastIndex].x1)
                
        elif Direction == 'Mz':
            
            for segment in self.SegmentsZ:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return segment.Moment(x - segment.x1)
                
            if isclose(x, self.L()):
                
                lastIndex = len(self.SegmentsZ) - 1
                return self.SegmentsZ[lastIndex].Moment(x - self.SegmentsZ[lastIndex].x1)
            
#%%
    def MaxMoment(self, Direction, combo_name='Combo 1'):
        '''
        Returns the maximum moment in the member for the given direction.
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the maximum moment. Must be one of the following:
                'My' = Moment about the local y-axis.
                'Mz' = Moment about the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        if Direction == 'Mz':
            
            Mmax = self.SegmentsZ[0].Moment(0)

            for segment in self.SegmentsZ:
                
                if segment.MaxMoment() > Mmax:
                    
                    Mmax = segment.MaxMoment()
                    
        if Direction == 'My':
            
            Mmax = self.SegmentsY[0].Moment(0)

            for segment in self.SegmentsY:
                
                if segment.MaxMoment() > Mmax:
                    
                    Mmax = segment.MaxMoment()
        
        return Mmax

#%%
    def MinMoment(self, Direction, combo_name='Combo 1'):
        '''
        Returns the minimum moment in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the minimum moment. Must be one of the following:
                'My' = Moment about the local y-axis.
                'Mz' = Moment about the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)   
            self.__solved_combo = self.LoadCombos[combo_name]
        
        if Direction == 'Mz':
            
            Mmin = self.SegmentsZ[0].Moment(0)

            for segment in self.SegmentsZ:
                
                if segment.MinMoment() < Mmin:
                    
                    Mmin = segment.MinMoment()
                    
        if Direction == 'My':
            
            Mmin = self.SegmentsY[0].Moment(0)

            for segment in self.SegmentsY:
                
                if segment.MinMoment() < Mmin:
                    
                    Mmin = segment.MinMoment()
        
        return Mmin

#%%
    def PlotMoment(self, Direction, combo_name='Combo 1'):
        '''
        Plots the moment diagram for the member
        
        Parameters
        ----------
        Direction : string
            The direction to plot the moment for.
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
                
        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt
        
        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x = []
        M = []
        
        # Calculate the moment diagram
        for i in range(21):
            
            x.append(self.L()/20*i)
            M.append(self.Moment(Direction, self.L()/20*i, combo_name))

        Member3D.__plt.plot(x, M)
        Member3D.__plt.ylabel('Moment')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.Name + '\n' + combo_name)
        Member3D.__plt.show()
       
#%%
    def Torsion(self, x, combo_name='Combo 1'):
        '''
        Returns the torsional moment at a point along the member's length
        
        Parameters
        ----------
        x : number
            The location at which to find the torque
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
            
        # Check which segment 'x' falls on
        for segment in self.SegmentsX:
            if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                return segment.Torsion()
                
            if isclose(x, self.L()):  
                lastIndex = len(self.SegmentsX) - 1
                return self.SegmentsX[lastIndex].Torsion()

#%%
    def MaxTorsion(self, combo_name='Combo 1'):
        '''
        Returns the maximum torsional moment in the member.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]       
        
        Tmax = self.SegmentsX[0].Torsion()   
        
        for segment in self.SegmentsX:

            if segment.MaxTorsion() > Tmax:
                    
                Tmax = segment.MaxTorsion()
        
        return Tmax
    
#%%
    def MinTorsion(self, combo_name='Combo 1'):
        '''
        Returns the minimum torsional moment in the member.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        Tmin = self.SegmentsX[0].Torsion()
            
        for segment in self.SegmentsX:
                
            if segment.MinTorsion() < Tmin:
                    
                Tmin = segment.MinTorsion()
        
        return Tmin

#%%
    def PlotTorsion(self, combo_name='Combo 1'):
        '''
        Plots the axial force diagram for the member.
        
        Paramters
        ---------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt

        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x = []
        T = []
        
        # Calculate the torsional moment diagram
        for i in range(21):
            x.append(self.L()/20*i)
            T.append(self.Torsion(self.L()/20*i, combo_name))

        Member3D.__plt.plot(x, T)
        Member3D.__plt.ylabel('Torsional Moment (Warping Torsion Not Included)') # Torsion results are for pure torsion. Torsional warping has not been considered
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.Name + '\n' + combo_name)
        Member3D.__plt.show()   
        
#%%
    def Axial(self, x, combo_name='Combo 1'):
        '''
        Returns the axial force at a point along the member's length.
        
        Parameters
        ----------
        x : number
            The location at which to find the axial force.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
            
        # Check which segment 'x' falls on
        for segment in self.SegmentsZ:
            if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                return segment.Axial(x - segment.x1)
                
            if isclose(x, self.L()):  
                lastIndex = len(self.SegmentsZ) - 1
                return self.SegmentsZ[lastIndex].Axial(x - self.SegmentsZ[lastIndex].x1)

#%%
    def MaxAxial(self, combo_name='Combo 1'):
        '''
        Returns the maximum axial force in the member

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        Pmax = self.SegmentsZ[0].Axial(0)   
        
        for segment in self.SegmentsZ:

            if segment.MaxAxial() > Pmax:
                    
                Pmax = segment.MaxAxial()
        
        return Pmax
    
#%%
    def MinAxial(self, combo_name='Combo 1'):
        '''
        Returns the minimum axial force in the member.
        
        Paramters
        ---------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        Pmin = self.SegmentsZ[0].Axial(0)
            
        for segment in self.SegmentsZ:
                
            if segment.MinAxial() < Pmin:
                    
                Pmin = segment.MinAxial()
        
        return Pmin
    
#%%
    def PlotAxial(self, combo_name='Combo 1'):
        '''
        Plots the axial force diagram for the member.
        
        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt

        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x = []
        P = []
        
        # Calculate the axial force diagram
        for i in range(21):
            x.append(self.L()/20*i)
            P.append(self.Axial(self.L()/20*i, combo_name))

        Member3D.__plt.plot(x, P)
        Member3D.__plt.ylabel('Axial Force')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.Name + '\n' + combo_name)
        Member3D.__plt.show()    
                        
#%%
    def Deflection(self, Direction, x, combo_name='Combo 1'):
        '''
        Returns the deflection at a point along the member's length.
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the deflection. Must be one of the following:
                'dx' = Deflection in the local x-axis.
                'dy' = Deflection in the local y-axis.
                'dz' = Deflection in the local z-axis.
        x : number
            The location at which to find the deflection.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        # Check which axis is of interest
        if Direction == 'dx':
            
            # Check which segment 'x' falls on
            for segment in self.SegmentsZ:
                
                if round(x, 10) >= round(segment.x1, 10) and round(x, 10) < round(segment.x2, 10):
                    return segment.AxialDeflection(x - segment.x1)
                
            if isclose(x, self.L()):
                
                lastIndex = len(self.SegmentsZ) - 1
                return self.SegmentsZ[lastIndex].AxialDeflection(x - self.SegmentsZ[lastIndex].x1)

        elif Direction == 'dy':
            
            # Check which segment 'x' falls on
            for segment in self.SegmentsZ:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return segment.Deflection(x - segment.x1)
                
            if isclose(x, self.L()):
                
                lastIndex = len(self.SegmentsZ) - 1
                return self.SegmentsZ[lastIndex].Deflection(x - self.SegmentsZ[lastIndex].x1)
                
        elif Direction == 'dz':
            
            for segment in self.SegmentsY:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return segment.Deflection(x - segment.x1)
                
            if isclose(x, self.L()):
                
                lastIndex = len(self.SegmentsY) - 1
                return self.SegmentsY[lastIndex].Deflection(x - self.SegmentsY[lastIndex].x1) 

#%%
    def MaxDeflection(self, Direction, combo_name='Combo 1'):
        '''
        Returns the maximum deflection in the member.
        
        Parameters
        ----------
        Direction : {'dy', 'dz'}
            The direction in which to find the maximum deflection.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        # Initialize the maximum deflection
        dmax = self.Deflection(Direction, 0, combo_name)
        
        # Check the deflection at 100 locations along the member and find the largest value
        for i in range(100):
            d = self.Deflection(Direction, self.L()*i/99, combo_name)
            if d > dmax:
                dmax = d
        
        # Return the largest value
        return dmax
    
#%%
    def MinDeflection(self, Direction, combo_name='Combo 1'):
        '''
        Returns the minimum deflection in the member.
        
        Parameters
        ----------
        Direction : {'dy', 'dz'}
            The direction in which to find the minimum deflection.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        # Initialize the minimum deflection
        dmin = self.Deflection(Direction, 0, combo_name)
        
        # Check the deflection at 100 locations along the member and find the smallest value
        for i in range(100):
            d = self.Deflection(Direction, self.L()*i/99, combo_name)
            if d < dmin:
                dmin = d
        
        # Return the smallest value
        return dmin
              
#%%
    def PlotDeflection(self, Direction, combo_name='Combo 1'):
        '''
        Plots the deflection diagram for the member
        
        Parameters
        ----------
        Direction : {'dy', 'dz'}
            The direction in which to plot the deflection.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
                
        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt
        
        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x = []
        d = []
        
        # Calculate the deflection diagram
        for i in range(21):
            
            x.append(self.L()/20*i)
            d.append(self.Deflection(Direction, self.L()/20*i, combo_name))

        Member3D.__plt.plot(x, d)
        Member3D.__plt.ylabel('Deflection')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.Name + '\n' + combo_name)
        Member3D.__plt.show()
    
#%%   
    def RelativeDeflection(self, Direction, x, combo_name='Combo 1'):
        '''
        Returns the relative deflection at a point along the member's length
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the relative deflection. Must be one of the following:
                'dy' = Deflection in the local y-axis
                'dz' = Deflection in the local x-axis
        x : number
            The location at which to find the relative deflection
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
        
        d = self.d(self.LoadCombos[combo_name])
        dyi = d[1,0]
        dyj = d[7,0]
        dzi = d[2,0]
        dzj = d[8,0]
        L = self.L()
       
        # Check which axis is of interest
        if Direction == 'dy':
            
            # Check which segment 'x' falls on
            for segment in self.SegmentsZ:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return (segment.Deflection(x - segment.x1)) - (dyi + (dyj-dyi)/L*x)
                
            if isclose(x, self.L()):
                
                lastIndex = len(self.SegmentsZ) - 1
                return (self.SegmentsZ[lastIndex].Deflection(x - self.SegmentsZ[lastIndex].x1))-dyj
                
        elif Direction == 'dz':
            
            for segment in self.SegmentsY:
                
                if round(x,10) >= round(segment.x1,10) and round(x,10) < round(segment.x2,10):
                    
                    return (segment.Deflection(x - segment.x1)) - (dzi + (dzj-dzi)/L*x)
                
            if isclose(x, self.L()):
                
                lastIndex = len(self.SegmentsY) - 1
                return (self.SegmentsY[lastIndex].Deflection(x - self.SegmentsY[lastIndex].x1)) - dzj

#%%
    def PlotRelativeDeflection(self, Direction, combo_name='Combo 1'):
        '''
        Plots the deflection diagram for the member
        
        Parameters
        ----------
        Direction : {'dy', 'dz'}
            The direction in which to plot the deflection.
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        '''
        
        # Segment the member if necessary
        if self.__solved_combo == None or combo_name != self.__solved_combo.name:
            self.SegmentMember(combo_name)
            self.__solved_combo = self.LoadCombos[combo_name]
                
        # Import 'pyplot' if not already done
        if Member3D.__plt is None:
            from matplotlib import pyplot as plt
            Member3D.__plt = plt
        
        fig, ax = Member3D.__plt.subplots()
        ax.axhline(0, color='black', lw=1)
        ax.grid()
        
        x = []
        d_relative = []
        
        # Calculate the relative deflection diagram
        for i in range(21):
            
            x.append(self.L()/20*i)
            d_relative.append(self.RelativeDeflection(Direction, self.L()/20*i, combo_name))

        Member3D.__plt.plot(x, d_relative)
        Member3D.__plt.ylabel('Relative Deflection')
        Member3D.__plt.xlabel('Location')
        Member3D.__plt.title('Member ' + self.Name + '\n' + combo_name)
        Member3D.__plt.show()   
        
#%%    
    # Divides the element up into mathematically continuous segments along each axis
    def SegmentMember(self, combo_name='Combo 1'):
        """
        A method that creates a member of the member

        Args:
            self: (todo): write your description
            combo_name: (str): write your description
        """
        
        # Get the member's length and stiffness properties
        L = self.L()
        E = self.E
        A = self.A
        Iz = self.Iz
        Iy = self.Iy
        SegmentsZ = self.SegmentsZ
        SegmentsY = self.SegmentsY
        SegmentsX = self.SegmentsX
        
        # Get the load combination to segment the member for
        combo = self.LoadCombos[combo_name]

        # Create a list of discontinuity locations
        disconts = [0, L] # Member ends
        
        for load in self.PtLoads: 
            disconts.append(load[2]) # Point load locations
        
        for load in self.DistLoads: 
            disconts.append(load[3]) # Distributed load start locations
            disconts.append(load[4]) # Distributed load end locations
        
        # Sort the list and eliminate duplicate values
        disconts = sorted(set(disconts))
        
        # Clear out old data from any previous analyses
        SegmentsZ.clear()
        SegmentsY.clear()
        SegmentsX.clear()
        
        # Create a list of mathematically continuous segments for each direction
        for index in range(len(disconts) - 1):
            
            # z-direction segments (bending about local z-axis)
            newSeg = BeamSegZ()           # Create the new segment
            newSeg.x1 = disconts[index]   # Segment start location
            newSeg.x2 = disconts[index+1] # Segment end location
            newSeg.EI = E*Iz              # Segment flexural stiffness
            newSeg.EA = E*A               # Segment axial stiffness
            SegmentsZ.append(newSeg)      # Add the segment to the list
            
            # y-direction segments (bending about local y-axis)
            newSeg = BeamSegY()           # Create the new segment
            newSeg.x1 = disconts[index]   # Segment start location
            newSeg.x2 = disconts[index+1] # Segment end location
            newSeg.EI = E*Iy              # Segment flexural stiffness
            newSeg.EA = E*A               # Segment axial stiffness
            SegmentsY.append(newSeg)      # Add the segment to the list
            
            # x-direction segments (for torsional moment)
            newSeg = BeamSegZ()           # Create the new segment
            newSeg.x1 = disconts[index]   # Segment start location
            newSeg.x2 = disconts[index+1] # Segment end location
            newSeg.EA = E*A               # Segment axial stiffness
            SegmentsX.append(newSeg)      # Add the segment to the list
        
        # Get the member local end forces, local fixed end reactions, and local displacements
        f = self.f(combo_name)           # Member local end force vector
        fer = self.__fer_Unc(combo_name) # Member local fixed end reaction vector
        d = self.d(combo_name)           # Member local displacement vector
        
        # Get the local deflections and calculate the slope at the start of the member
        # Note 1: The slope may not be available directly from the local displacement vector if member end releases have been used,
        #         so slope-deflection has been applied to solve for it.
        # Note 2: The traditional slope-deflection equations assume a sign convention opposite of what PyNite uses for moments about
        #         the local y-axis, so a negative value has been applied to those values specifically.
        m1z = f[5, 0]       # local z-axis moment at start of member
        m2z = f[11, 0]      # local z-axis moment at end of member
        m1y = -f[4, 0]      # local y-axis moment at start of member
        m2y = -f[10, 0]     # local y-axis moment at end of member
        fem1z = fer[5, 0]   # local z-axis fixed end moment at start of member
        fem2z = fer[11, 0]  # local z-axis fixed end moment at end of member
        fem1y = -fer[4, 0]  # local y-axis fixed end moment at start of member
        fem2y = -fer[10, 0] # local y-axis fixed end moment at end of member
        delta1y = d[1, 0]   # local y displacement at start of member
        delta2y = d[7, 0]   # local y displacement at end of member
        delta1z = d[2, 0]   # local z displacement at start of member
        delta2z = d[8, 0]   # local z displacement at end of member
        SegmentsZ[0].delta1 = delta1y
        SegmentsY[0].delta1 = delta1z
        SegmentsZ[0].theta1 = 1/3*((m1z - fem1z)*L/(E*Iz) - (m2z - fem2z)*L/(2*E*Iz) + 3*(delta2y - delta1y)/L)
        SegmentsY[0].theta1 = -1/3*((m1y - fem1y)*L/(E*Iy) - (m2y - fem2y)*L/(2*E*Iy) + 3*(delta2z - delta1z)/L)

        # Add the axial deflection at the start of the member
        SegmentsZ[0].delta_x1 = d[0, 0]
        SegmentsY[0].delta_x1 = d[0, 0]
        SegmentsX[0].delta_x1 = d[0, 0]
        
        # Add loads to each segment
        for i in range(len(SegmentsZ)):
            
            # Get the starting point of the segment
            x = SegmentsZ[i].x1
            
            # Initialize the distributed loads on the segment to zero
            SegmentsZ[i].w1 = 0
            SegmentsZ[i].w2 = 0
            SegmentsZ[i].p1 = 0
            SegmentsZ[i].p2 = 0
            SegmentsY[i].w1 = 0
            SegmentsY[i].w2 = 0
            SegmentsY[i].p1 = 0
            SegmentsY[i].p2 = 0
            
            # Initialize the slope and displacement at the start of the segment
            if i > 0: # The first segment has already been initialized
                SegmentsZ[i].theta1 = SegmentsZ[i-1].Slope(SegmentsZ[i-1].Length())
                SegmentsZ[i].delta1 = SegmentsZ[i-1].Deflection(SegmentsZ[i-1].Length())
                SegmentsZ[i].delta_x1 = SegmentsZ[i-1].AxialDeflection(SegmentsZ[i-1].Length())
                SegmentsY[i].theta1 = SegmentsY[i-1].Slope(SegmentsY[i-1].Length())
                SegmentsY[i].delta1 = SegmentsY[i-1].Deflection(SegmentsY[i-1].Length())
                SegmentsY[i].delta_x1 = SegmentsY[i-1].AxialDeflection(SegmentsY[i-1].Length())
                
            # Add the effects of the beam end forces to the segment
            SegmentsZ[i].P1 = f[0, 0]
            SegmentsZ[i].V1 = f[1, 0]
            SegmentsZ[i].M1 = f[5, 0] - f[1, 0]*x
            SegmentsY[i].P1 = f[0, 0]
            SegmentsY[i].V1 = f[2, 0]
            SegmentsY[i].M1 = f[4, 0] + f[2, 0]*x
            SegmentsX[i].T1 = f[3, 0]
            
            for case, factor in combo.factors.items():
            
                # Add effects of point loads occuring prior to this segment
                for ptLoad in self.PtLoads:
                    
                    if round(ptLoad[2],10) <= round(x,10) and case == ptLoad[3]:
                    
                        if ptLoad[0] == 'Fx':
                            SegmentsZ[i].P1 += factor*ptLoad[1]
                        elif ptLoad[0] == 'Fy':
                            SegmentsZ[i].V1 += factor*ptLoad[1]
                            SegmentsZ[i].M1 -= factor*ptLoad[1]*(x - ptLoad[2])
                        elif ptLoad[0] == 'Fz':
                            SegmentsY[i].V1 += factor*ptLoad[1]
                            SegmentsY[i].M1 += factor*ptLoad[1]*(x - ptLoad[2])
                        elif ptLoad[0] == 'Mx':
                            SegmentsX[i].T1 += factor*ptLoad[1]    
                        elif ptLoad[0] == 'My':
                            SegmentsY[i].M1 += factor*ptLoad[1]
                        elif ptLoad[0] == 'Mz':
                            SegmentsZ[i].M1 += factor*ptLoad[1]
            
                # Add distributed loads to the segment
                for distLoad in self.DistLoads:
                    
                    if case == distLoad[5]:
                    
                        # Get the parameters for the distributed load
                        Direction = distLoad[0]
                        w1 = factor*distLoad[1]
                        w2 = factor*distLoad[2]
                        x1 = distLoad[3]
                        x2 = distLoad[4]
            
                        # Determine if the load affects the segment
                        if round(x1,10) <= round(x,10):
                    
                            if Direction == 'Fx':
                        
                                # Determine if the load ends after the start of the segment
                                if round(x2,10) > round(x,10):
                                                
                                    # Break up the load and place it on the segment
                                    # Note that 'w1' and 'w2' are really 'p1' and 'p2' here
                                    SegmentsZ[i].p1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                                    SegmentsZ[i].p2 += (w2 - w1)/(x2 - x1)*(SegmentsZ[i].x2 - x1) + w1
                                    SegmentsY[i].p1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                                    SegmentsY[i].p2 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x2 - x1) + w1
                            
                                    # Calculate the magnitude of the load at the start of the segment
                                    w2 = w1+(w2-w1)/(x2-x1)*(x-x1)
                                    x2 = x
                        
                                # Calculate the axial force at the start of the segment
                                SegmentsZ[i].P1 -= (w1 + w2)/2*(x2 - x1)
                                SegmentsY[i].P1 -= (w1 + w2)/2*(x2 - x1)
                    
                            elif Direction == 'Fy':
                        
                                # Determine if the load ends after the start of the segment
                                if round(x2,10) > round(x,10):
                                                
                                    # Break up the load and place it on the segment
                                    SegmentsZ[i].w1 += (w2 - w1)/(x2 - x1)*(x - x1) + w1
                                    SegmentsZ[i].w2 += (w2 - w1)/(x2 - x1)*(SegmentsZ[i].x2 - x1) + w1
                            
                                    # Calculate the magnitude of the load at the start of the segment
                                    # This will be used as the 'x2' value for the load prior to the start of the segment
                                    w2 = w1 + (w2 - w1)/(x2 - x1)*(x - x1)
                                    x2 = x
                        
                                # Calculate the shear and moment at the start of the segment due to the load
                                SegmentsZ[i].V1 += (w1 + w2)/2*(x2 - x1)
                                SegmentsZ[i].M1 -= (x1 - x2)*(2*w1*x1 - 3*w1*x + w1*x2 + w2*x1 - 3*w2*x + 2*w2*x2)/6
                    
                            elif Direction == 'Fz':
                        
                                # Determine if the load ends after the start of the segment
                                if round(x2,10) > round(x,10):
                                                
                                    # Break up the load and place it on the segment
                                    SegmentsY[i].w1 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x1 - x1) + w1
                                    SegmentsY[i].w2 += (w2 - w1)/(x2 - x1)*(SegmentsY[i].x2 - x1) + w1
                            
                                    # Calculate the magnitude of the load at the start of the segment
                                    w2 = w1 + (w2 - w1)/(x2 - x1)*(x - x1)
                                    x2 = x
                        
                                # Calculate the shear and moment at the start of the segment due to the load
                                SegmentsY[i].V1 += (w1 + w2)/2*(x2 - x1)
                                SegmentsY[i].M1 += (x1 - x2)*(2*w1*x1 - 3*w1*x + w1*x2 + w2*x1 - 3*w2*x + 2*w2*x2)/6
