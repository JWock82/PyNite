# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 20:58:03 2017

@author: D. Craig Brinck, SE
"""
# %%
from numpy import zeros

# %%
def FER_PtLoad(P, x, L, Direction):
    """
    Returns the fixed end reaction vector for a point load
    
    Parameters:
    -----------
    P : number
        The magnitude of the point load
    x : number
        The location of the point load relative to the start of the member
    L : number
        The length of the member
    Direction : string
        The direction of the point load. Must be one of the following:
            "Fy" = Force on the member's local y-axis
            "Fz" = Force on the member's local z-axis
    """
    # Define variables
    b = L - x
    
    # Create the fixed end reaction vector
    FER = zeros((12, 1))
    
    # Populate the fixed end reaction vector
    if Direction == "Fy":
        FER.itemset((1, 0), -P*b**2*(L+2*x)/L**3)
        FER.itemset((5, 0), -P*x*b**2/L**2)
        FER.itemset((7, 0), -P*x**2*(L+2*b)/L**3)
        FER.itemset((11, 0), P*x**2*b/L**2)
    elif Direction == "Fz":
        FER.itemset((2, 0), -P*b**2*(L+2*x)/L**3)
        FER.itemset((4, 0), P*x*b**2/L**2)
        FER.itemset((8, 0), -P*x**2*(L+2*b)/L**3)
        FER.itemset((10, 0), -P*x**2*b/L**2)
        
    return FER
    
# %%
def FER_Moment(M, x, L, Direction):
    """
    Returns the fixed end reaction vector for a concentrated moment
    
    Parameters
    ----------
    M : number
        The magnitude of the moment
    x : number
        The location of the moment relative to the start of the member
    Direction : string
        The direction of the moment. Must be one of the following:
            "My" = Moment applied about the local y-axis
            "Mz" = Moment applied about the local z-axis
    """
    
    # Define variables
    b = L - x
    
    # Create the fixed end reaction vector
    FER = zeros((12, 1))
    
    # Populate the fixed end reaction vector
    if Direction == "Mz":
        FER.itemset((1, 0), 6*M*x*b/L**3)
        FER.itemset((5, 0), M*b*(2*x-b)/L**2)
        FER.itemset((7, 0), -M*x*b/L**3)
        FER.itemset((11, 0), M*x*(2*b-x)/L**2)
    elif Direction == "My":
        FER.itemset((2, 0), -6*M*x*b/L**3)
        FER.itemset((4, 0), M*b*(2*x-b)/L**2)
        FER.itemset((8, 0), M*x*b/L**3)
        FER.itemset((10, 0), M*x*(2*b-x)/L**2)
    return FER
    
# %%
# Returns the fixed end reaction vector for a linear distributed load
def FER_LinLoad(w1, w2, x1, x2, L, Direction):
    """
    R calculates the coefficient of a liquid in [ 2 ] _1.

    Args:
        w1: (str): write your description
        w2: (str): write your description
        x1: (str): write your description
        x2: (str): write your description
        L: (str): write your description
        Direction: (str): write your description
    """
        
    # Create the fixed end reaction vector
    FER = zeros((12, 1))
    
    # Populate the fixed end reaction vector
    if Direction == 'Fy':
        FER.itemset((1, 0), (x1 - x2)*(10*L**3*w1 + 10*L**3*w2 - 15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 8*w1*x1**3 + 6*w1*x1**2*x2 + 4*w1*x1*x2**2 + 2*w1*x2**3 + 2*w2*x1**3 + 4*w2*x1**2*x2 + 6*w2*x1*x2**2 + 8*w2*x2**3)/(20*L**3))
        FER.itemset((5, 0), (x1 - x2)*(20*L**2*w1*x1 + 10*L**2*w1*x2 + 10*L**2*w2*x1 + 20*L**2*w2*x2 - 30*L*w1*x1**2 - 20*L*w1*x1*x2 - 10*L*w1*x2**2 - 10*L*w2*x1**2 - 20*L*w2*x1*x2 - 30*L*w2*x2**2 + 12*w1*x1**3 + 9*w1*x1**2*x2 + 6*w1*x1*x2**2 + 3*w1*x2**3 + 3*w2*x1**3 + 6*w2*x1**2*x2 + 9*w2*x1*x2**2 + 12*w2*x2**3)/(60*L**2))
        FER.itemset((7, 0), -(x1 - x2)*(-15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 8*w1*x1**3 + 6*w1*x1**2*x2 + 4*w1*x1*x2**2 + 2*w1*x2**3 + 2*w2*x1**3 + 4*w2*x1**2*x2 + 6*w2*x1*x2**2 + 8*w2*x2**3)/(20*L**3))
        FER.itemset((11, 0), (x1 - x2)*(-15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 12*w1*x1**3 + 9*w1*x1**2*x2 + 6*w1*x1*x2**2 + 3*w1*x2**3 + 3*w2*x1**3 + 6*w2*x1**2*x2 + 9*w2*x1*x2**2 + 12*w2*x2**3)/(60*L**2))
    elif Direction == 'Fz':
        FER.itemset((2, 0), (x1 - x2)*(10*L**3*w1 + 10*L**3*w2 - 15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 8*w1*x1**3 + 6*w1*x1**2*x2 + 4*w1*x1*x2**2 + 2*w1*x2**3 + 2*w2*x1**3 + 4*w2*x1**2*x2 + 6*w2*x1*x2**2 + 8*w2*x2**3)/(20*L**3))
        FER.itemset((4, 0), -(x1 - x2)*(20*L**2*w1*x1 + 10*L**2*w1*x2 + 10*L**2*w2*x1 + 20*L**2*w2*x2 - 30*L*w1*x1**2 - 20*L*w1*x1*x2 - 10*L*w1*x2**2 - 10*L*w2*x1**2 - 20*L*w2*x1*x2 - 30*L*w2*x2**2 + 12*w1*x1**3 + 9*w1*x1**2*x2 + 6*w1*x1*x2**2 + 3*w1*x2**3 + 3*w2*x1**3 + 6*w2*x1**2*x2 + 9*w2*x1*x2**2 + 12*w2*x2**3)/(60*L**2))
        FER.itemset((8, 0), -(x1 - x2)*(-15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 8*w1*x1**3 + 6*w1*x1**2*x2 + 4*w1*x1*x2**2 + 2*w1*x2**3 + 2*w2*x1**3 + 4*w2*x1**2*x2 + 6*w2*x1*x2**2 + 8*w2*x2**3)/(20*L**3))
        FER.itemset((10, 0), -(x1 - x2)*(-15*L*w1*x1**2 - 10*L*w1*x1*x2 - 5*L*w1*x2**2 - 5*L*w2*x1**2 - 10*L*w2*x1*x2 - 15*L*w2*x2**2 + 12*w1*x1**3 + 9*w1*x1**2*x2 + 6*w1*x1*x2**2 + 3*w1*x2**3 + 3*w2*x1**3 + 6*w2*x1**2*x2 + 9*w2*x1*x2**2 + 12*w2*x2**3)/(60*L**2))

    return FER

# %%
# Returns the fixed end reaction vector for an axial point load
def FER_AxialPtLoad(P, x, L):
    """
    R calculates the value of a pipe pt

    Args:
        P: (todo): write your description
        x: (todo): write your description
        L: (todo): write your description
    """
    
    # Create the fixed end reaction vector
    FER = zeros((12, 1))
    
    # Populate the fixed end reaction vector
    FER.itemset((0, 0), -P*(L-x)/L)
    FER.itemset((6, 0), -P*x/L)
    
    return FER

# %%
# Returns the fixed end reaction vector for a distributed axial load
def FER_AxialLinLoad(p1, p2, x1, x2, L):
    """
    Calculate the adjustment of two - phase

    Args:
        p1: (str): write your description
        p2: (str): write your description
        x1: (str): write your description
        x2: (str): write your description
        L: (str): write your description
    """
    
    # Create the fixed end reaction vector
    FER = zeros((12, 1))
    
    # Populate the fixed end reaction vector
    FER.itemset((0, 0), 1/(6*L)*(x1-x2)*(3*L*p1+3*L*p2-2*p1*x1-p1*x2-p2*x1-2*p2*x2))
    FER.itemset((6, 0), 1/(6*L)*(x1-x2)*(2*p1*x1+p1*x2+p2*x1+2*p2*x2))
    
    return FER
 
def FER_Torque(T, x, L):
    """
    Returns the fixed end reaction vector for a concentrated torque

    Parameters
    ----------
    T : number
        The magnitude of the torque
    x : number
        The location of the torque relative to the start of the member
    """

    # Create the fixed end reaction vector
    FER = zeros((12, 1))
    
    # Populate the fixed end reaction vector
    FER.itemset((3, 0), -T*(L - x)/L)
    FER.itemset((9, 0), -T*x/L)

    return FER
