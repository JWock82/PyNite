from numpy import array, dot, cross
from numpy.linalg import norm
from math import isclose, acos

from PyNite.Member3D import Member3D

class PhysMember(Member3D):
    """
    A physical member.

    Physical members can detect internal nodes and subdivide themselves into sub-members at those
    nodes.
    """

    def __init__(self, model, name, i_node, j_node, material_name, section_name, aux_node=None,
                 tension_only=False, comp_only=False):
        
        super().__init__(model, name, i_node, j_node, material_name, section_name, aux_node, tension_only, comp_only)
        self.sub_members = {}

    def descritize(self):
        """
        Subdivides the physical member into sub-members at each node along the physical member
        """

        # Clear out any old sub_members
        self.sub_members = {}

        # Start a new list of nodes along the member
        int_nodes = []

        # Create a vector from the i-node to the j-node
        Xi, Yi, Zi = self.i_node.X, self.i_node.Y, self.i_node.Z
        Xj, Yj, Zj = self.j_node.X, self.j_node.Y, self.j_node.Z
        vector_ij = array([Xj-Xi, Yj-Yi, Zj-Zi])

        # Add the i-node and j-node to the list
        int_nodes.append([self.i_node, 0])
        int_nodes.append([self.j_node, norm(vector_ij)])

        # Step through each node in the model
        for node in self.model.nodes.values():

            # Check each node in the model (except the i and j-nodes)
            if node is not self.i_node and node is not self.j_node:

                # Create a vector from the i-node to the current node
                X, Y, Z = node.X, node.Y, node.Z
                vector_in = array([X-Xi, Y-Yi, Z-Zi])

                # Calculate the angle between the two vectors
                angle = acos(round(dot(vector_in, vector_ij)/(norm(vector_in)*norm(vector_ij)), 10))

                # Determine if the node is colinear with the member
                if isclose(angle, 0):

                    # Determine if the node is on the member
                    if norm(vector_in) < norm(vector_ij):

                        # Add the node to the list of intermediate nodes
                        int_nodes.append([node, norm(vector_in)])
        
        # Create a list of sorted intermediate nodes by distance from the i-node
        int_nodes = sorted(int_nodes, key=lambda x: x[1])

        # Break up the member into sub-members at each intermediate node
        for i in range(len(int_nodes) - 1):

            # Generate the sub-member's name (physical member name + a, b, c, etc.)
            name = self.name + chr(i+97)

            # Find the i and j nodes for the sub-member, and their positions along the physical
            # member's local x-axis
            i_node = int_nodes[i][0]
            j_node = int_nodes[i+1][0]
            xi = int_nodes[i][1]
            xj = int_nodes[i+1][1]

            # Create a new sub-member
            new_sub_member = Member3D(self.model, name, i_node, j_node, self.material.name, self.section.name, self.auxNode, self.tension_only, self.comp_only)
            
            # Flag the sub-member as active
            for combo_name in self.model.load_combos.keys():
                new_sub_member.active[combo_name] = True

            # Apply end releases if applicable
            if i == 0:
                new_sub_member.Releases[0:6] = self.Releases[0:6]
            if i == len(int_nodes) - 2:
                new_sub_member.Releases[6:12] = self.Releases[6:12]

            # Add distributed to the sub-member
            for dist_load in self.DistLoads:
                
                # Find the start and end points of the distributed load in the physical member's
                # local coordinate system
                x1_load = dist_load[3]
                x2_load = dist_load[4]

                # Determine if the distributed load should be applied to this segment
                if x1_load <= xj and x2_load > xi: 
                    
                    direction = dist_load[0]
                    w1 = dist_load[1]
                    w2 = dist_load[2]
                    case = dist_load[5]

                    # Equation describing the load as a function of x
                    w = lambda x: (w2 - w1)/(x2_load - x1_load)*(x - x1_load) + w1

                    # Chop up the distributed load for the sub-member
                    if x1_load > xi:
                        x1 = x1_load - xi
                    else:
                        x1 = 0
                        w1 = w(xi)
                    
                    if x2_load < xj:
                        x2 = x2_load - xi
                    else:
                        x2 = xj - xi
                        w2 = w(xj)

                    # Add the load to the sub-member
                    new_sub_member.DistLoads.append([direction, w1, w2, x1, x2, case])

            # Add point loads to the sub-member
            for pt_load in self.PtLoads:
                
                direction = pt_load[0]
                P = pt_load[1]
                x = pt_load[2]
                case = pt_load[3]

                # Determine if the point load should be applied to this segment
                if x >= xi and x < xj or (isclose(x, xj) and isclose(xj, self.L())):

                    x = x - xi
                    
                    # Add the load to the sub-member
                    new_sub_member.PtLoads.append([direction, P, x, case])

            # Add the new sub-member to the sub-member dictionary for this physical member
            self.sub_members[name] = new_sub_member
    
    def shear(self, Direction, x, combo_name='Combo 1'):
        """
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
        """
        
        member, x_mod = self.find_member(x)
        return member.shear(Direction, x_mod, combo_name)
    
    def max_shear(self, Direction, combo_name='Combo 1'):
        """
        Returns the maximum shear in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the maximum shear. Must be one of the following:
                'Fy' = Shear acting on the local y-axis
                'Fz' = Shear acting on the local z-axis
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        """

        Vmax = None
        for member in self.sub_members.values():
            V = member.max_shear(Direction, combo_name)
            if Vmax is None or V > Vmax:
                Vmax = V
        return Vmax
    
    def min_shear(self, Direction, combo_name='Combo 1'):
        """
        Returns the minimum shear in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the minimum shear. Must be one of the following:
                'Fy' = Shear acting on the local y-axis
                'Fz' = Shear acting on the local z-axis
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """
        
        Vmin = None
        for member in self.sub_members.values():
            V = member.min_shear(Direction, combo_name)
            if Vmin is None or V < Vmin:
                Vmin = V
        return Vmin
   
    def moment(self, Direction, x, combo_name='Combo 1'):
        """
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
        """
        
        member, x_mod = self.find_member(x)
        return member.moment(Direction, x_mod, combo_name)
    
    def max_moment(self, Direction, combo_name='Combo 1'):
        """
        Returns the maximum moment in the member for the given direction.
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the maximum moment. Must be one of the following:
                'My' = Moment about the local y-axis.
                'Mz' = Moment about the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the combination itself).
        """

        Mmax = None
        for member in self.sub_members.values():
            M = member.max_moment(Direction, combo_name)
            if Mmax is None or M > Mmax:
                Mmax = M
        return Mmax
    
    def min_moment(self, Direction, combo_name='Combo 1'):
        """
        Returns the minimum moment in the member for the given direction
        
        Parameters
        ----------
        Direction : string
            The direction in which to find the minimum moment. Must be one of the following:
                'My' = Moment about the local y-axis.
                'Mz' = Moment about the local z-axis.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """
        
        Mmin = None
        for member in self.sub_members.values():
            M = member.min_moment(Direction, combo_name)
            if Mmin is None or M < Mmin:
                Mmin = M
        return Mmin

    def torque(self, x, combo_name='Combo 1'):
        """
        Returns the torsional moment at a point along the member's length
        
        Parameters
        ----------
        x : number
            The location at which to find the torque
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """
        
        member, x_mod = self.find_member(x)
        return member.torque(x_mod, combo_name)
    
    def max_torque(self, combo_name='Combo 1'):
        
        Tmax = None
        for member in self.sub_members.values():
            T = member.max_torque(combo_name)
            if Tmax is None or T > Tmax:
                Tmax = T
        return Tmax
    
    def min_torque(self, combo_name='Combo 1'):
        """
        Returns the minimum torsional moment in the member.

        Parameters
        ----------
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        Tmin = None
        for member in self.sub_members.values():
            T = member.min_torque(combo_name)
            if Tmin is None or T < Tmin:
                Tmin = T
        return Tmin

    def axial(self, x, combo_name='Combo 1'):
        """
        Returns the axial force at a point along the member's length.
        
        Parameters
        ----------
        x : number
            The location at which to find the axial force.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        member, x_mod = self.find_member(x)
        return member.axial(x_mod, combo_name)
    
    def max_axial(self, combo_name='Combo 1'):
        
        Pmax = None
        for member in self.sub_members.values():
            P = member.max_axial(combo_name)
            if Pmax is None or P > Pmax:
                Pmax = P
        return Pmax
    
    def min_axial(self, combo_name='Combo 1'):

        Pmin = None
        for member in self.sub_members.values():
            P = member.min_axial(combo_name)
            if Pmin is None or P < Pmin:
                Pmin = P
        return Pmin

    def deflection(self, Direction, x, combo_name='Combo 1'):
        """
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
        """

        member, x_mod = self.find_member(x)
        return member.deflection(Direction, x_mod, combo_name)

    def max_deflection(self, Direction, combo_name='Combo 1'):
        """
        Returns the maximum deflection in the member.
        
        Parameters
        ----------
        Direction : {'dy', 'dz'}
            The direction in which to find the maximum deflection.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        dmax = None
        for member in self.sub_members.values():
            d = member.max_deflection(Direction, combo_name)
            if dmax is None or d > dmax:
                dmax = d
        return dmax
    
    def min_deflection(self, Direction, combo_name='Combo 1'):
        """
        Returns the minimum deflection in the member.
        
        Parameters
        ----------
        Direction : {'dy', 'dz'}
            The direction in which to find the minimum deflection.
        combo_name : string
            The name of the load combination to get the results for (not the load combination itself).
        """

        dmin = None
        for member in self.sub_members.values():
            d = member.min_deflection(Direction, combo_name)
            if dmin is None or d < dmin:
                dmin = d
        return dmin

    def rel_deflection(self, Direction, x, combo_name='Combo 1'):
        """
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
        """
        
        member, x_mod = self.find_member(x)
        return member.rel_deflection(Direction, x_mod, combo_name)

    def find_member(self, x):
        """
        Returns the sub-member that the physical member's local point 'x' lies on, and 'x' modified for that sub-member's
        local coordinate system.
        """

        # Initialize a summation of sub-member lengths
        L = 0

        # Step through each sub-member (in order from start to end)
        for i, member in enumerate(self.sub_members.values()):

            # Sum the sub-member's length
            L += member.L()

            # Check if 'x' lies on this sub-member
            if x < L or (isclose(x, L) and i == len(self.sub_members.values()) - 1):

                # Return the sub-member, and a modified value for 'x' relative to the sub-member's
                # i-node
                return member, x - (L - member.L())

                # Exit the 'for' loop
                break
        else:
            raise ValueError(f"Location x={x} does not lie on this member")
            