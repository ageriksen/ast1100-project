"""
Rewrite of first code to clarify
"""

##############################################
# imports
import numpy as np
import matplotlib.pyplot as plt
import MySolarSystem as M
import LeapFrog as lf
import scipy as scp


##############################################
# functions

def e_theta(theta, r):
    """
    transforms desired angle off of radial 
    vector to angle off of 1. axis.
    returns radians. arcatn2 input is (y, x)
    """
    theta = theta + np.arctan2(r[1], r[0])

def launchPosition(r, R, e_theta, theta):
    """
    position along surface of planet 0 to
    launch from.
    r is the position vector of the planet, 
    R is the radius of the planet, 
    e_theta is function above
    theta is angle from r in radians. 
    """
    return r + R*e_theta(theta, r)

def accelerate(r, m): 
    """
    calculates acceleration on body
    from gravitational pull from 
    body of mass m at a distance r
    """
    r_ = np.linalg.norm(r) # lenght of distance
    return -( G*m*r )/( r_**3 )

