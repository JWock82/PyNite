class LoadCombo():
    '''
    A class representing a load combination
    '''

    def __init__(self, name, combo_type=None, factors={}):
        """
        Initialize a new factors object.

        Args:
            self: (todo): write your description
            name: (str): write your description
            combo_type: (todo): write your description
            factors: (float): write your description
        """
        
        self.name = name             # A unique user-defined name for the load combination
        self.combo_type = combo_type # Used to track what type of load combination (e.g. strength or serviceability)
                                     # The 'combo_type' is just a placeholder for future features for now
        self.factors = factors       # A dictionary containing each load case name and associated load factor
    
    def AddLoadCase(self, case_name, factor):
        '''
        Adds a load case with its associated load factor
        '''

        self.factors[case_name] = factor
    
    def DeleteLoadCase(self, case_name, factor):
        '''
        Deletes a load case with its associated load factor
        '''

        del self.factors[case_name]