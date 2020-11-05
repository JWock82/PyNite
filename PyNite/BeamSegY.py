from PyNite.BeamSegZ import BeamSegZ

# %%
class BeamSegY(BeamSegZ):

#%% 
    # Returns the moment at a location on the segment
    def Moment(self, x):
        '''
        Returns the moment at a location on the segment.

        Parameters
        ----------
        x : number
            Location (relative to start of segment) where moment is to be calculated
        '''

        V1 = self.V1
        M1 = self.M1
        w1 = self.w1
        w2 = self.w2
        L = self.Length()
        
        return M1 + V1*x + w1*x**2/2 + x**3*(-w1 + w2)/(6*L)

    
#%%
    def Slope(self, x):
        """
        Returns the slope at a point on the segment
        
        Parameters
        ----------
        x : number
          Location (relative to start of segment) where slope is to be calculated
        EI : number
          Flexural stiffness of the segment
        
        Returns
        -------
        Slope : number
          The slope of the segment (radians) at location "x"
        
        Notes
        -----
        Any unit system may be used as long as the units are consistent with each other
        
        """  
        
        V1 = self.V1
        M1 = self.M1
        w1 = self.w1
        w2 = self.w2
        theta1 = self.theta1
        L = self.Length()
        EI = self.EI
        
        return theta1 - (M1*x + V1*x**2/2 + w1*x**3/6 + x**4*(-w1 + w2)/(24*L))/(EI)

#%%
    # Returns the deflection at a location on the segment
    def Deflection(self, x):
        """
        Evaluate of a x and b.

        Args:
            self: (todo): write your description
            x: (todo): write your description
        """
        
        V1 = self.V1
        M1 = self.M1
        w1 = self.w1
        w2 = self.w2
        theta1 = self.theta1
        delta1 = self.delta1
        L = self.Length()
        EI = self.EI
        
        return delta1 - theta1*x + M1*x**2/(2*EI) + V1*x**3/(6*EI) + w1*x**4/(24*EI) - x**5*(w1 - w2)/(120*EI*L)

    
#%%
    # Returns the maximum moment in the segment
    def MaxMoment(self):
        """
        Calculate the maximum of the maximum number of a - axis.

        Args:
            self: (todo): write your description
        """
        
        w1 = self.w1
        w2 = self.w2
        V1 = self.V1
        L = self.Length()
        
        # Find the quadratic equation parameters
        a = (w2-w1)/(2*L)
        b = w1
        c = V1
    
        # Determine possible locations of maximum moment
        if a == 0:
            if b != 0:
                x1 = -c/b
            else:
                x1 = 0
            x2 = 0
        elif b**2-4*a*c < 0:
            x1 = 0
            x2 = 0
        else:
            x1 = (-b+(b**2-4*a*c)**0.5)/(2*a)
            x2 = (-b-(b**2-4*a*c)**0.5)/(2*a)
    
        x3 = 0
        x4 = L
    
        if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
            x1 = 0
    
        if round(x2, 10) < 0 or round(x2, 10) > round(L, 10):
            x2 = 0
    
        # Find the moment at each location of interest
        M1 = self.Moment(x1)
        M2 = self.Moment(x2)
        M3 = self.Moment(x3)
        M4 = self.Moment(x4)
    
        # Return the maximum moment
        return max(M1, M2, M3, M4)

#%%
    # Returns the minimum moment in the segment
    def MinMoment(self):
        """
        Calculate the minimum point of the line.

        Args:
            self: (todo): write your description
        """
        
        w1 = self.w1
        w2 = self.w2
        V1 = self.V1
        L = self.Length()
        
        # Find the quadratic equation parameters
        a = (w2-w1)/(2*L)
        b = w1
        c = V1
    
        # Determine possible locations of minimum moment
        if a == 0:
            if b != 0:
                x1 = -c/b
            else:
                x1 = 0
            x2 = 0
        elif b**2-4*a*c < 0:
            x1 = 0
            x2 = 0
        else:
            x1 = (-b+(b**2-4*a*c)**0.5)/(2*a)
            x2 = (-b-(b**2-4*a*c)**0.5)/(2*a)
    
        x3 = 0
        x4 = L
    
        if round(x1, 10) < 0 or round(x1, 10) > round(L, 10):
            x1 = 0
    
        if round(x2, 10) < 0 or round(x2, 10) > round(L, 10):
            x2 = 0
    
        # Find the moment at each location of interest
        M1 = self.Moment(x1)
        M2 = self.Moment(x2)
        M3 = self.Moment(x3)
        M4 = self.Moment(x4)
    
        # Return the minimum moment
        return min(M1, M2, M3, M4)