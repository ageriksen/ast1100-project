"""
script to area for a terminal velocity of 
3 m/s at surface of planet 1
"""

def Area_needed(vt, rho0, M, R, m):
    return (2./9) * (G*M*m) / (R**2 * rho0 * vt**2)
def input():
    print 'Area_needed takes:'
    print 'vt, rho0, M, R, m'

if __name__ == '__main__':
    import numpy as np
    import MySolarSystem as MSS
    G = 6.674e-11 # N * m^2 * kg^-2 | grav const.
    M_sol = 1.9891e30 #kg in 1 solar mass
    myss = MSS.Myseed()
    M = MSS.p_mass(myss, 1)*M_sol
    print 'converted to kg', M
    R = MSS.p_radius(myss, 1)*1e3
    print 'converted to m', R
    rho0 = MSS.atmosphere(myss, 1)
    vt = 3 # m/s terminal velocity at surface
    m = 1100 # kg | mass of satellite
    area = Area_needed(vt, rho0, M, R, m)
    print 'area needed: ', area
"""
run example:
    myseed 80101 <type 'int'>                                         
    homeplanet is planet 0
    mass of planets, [solar masses]
    converted to kg 1.79733516331e+25
    radi of planets, [km]
    converted to m 9327651.84279
    Atmospheric density at surface, [kg/m^3]
    area needed:  265.825542201
"""   
