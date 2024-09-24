from typing import Optional

class Material():
    """
    A class representing a material assigned to a Member3D, Plate or Quad in a finite element model.

    This class stores all properties related to the physical material of the element
    """
    def __init__(self, model, name:str, E:float, G:float, nu:float, rho:float, fy:Optional[float] = None):

        self.model = model
        self.name = name
        self.E = E
        self.G = G
        self.nu = nu
        self.rho = rho
        self.fy = fy