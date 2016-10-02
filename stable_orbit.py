"""
aims to calculate the stable orbit velocity of a sattelite around a planet
"""

def stable(M,  R, nr, dim='SI'):
    # takes inn mass of planet M, grav. const. G, 
    # planet radius R. 
    from numpy import sqrt
    AU = 1.496e11 # km in 1 AU
    M_sol = 1.988e30 # kg in 1 Solar mass
    if dim == 'SI':
        G = 6.67e-11 # N*m^2*kg^-2
        if M < 100:
            print 'M needs to be in kg'
            M = M*M_sol
        elif R < 1:
            print 'R needs to be in km'
            R = R*AU
        v_stab = sqrt((M*G)/R)
        print """
stable orbit velocity is: %g m/s for planet %g
              """ %(v_stab, nr)
        return v_stab
    elif dim == 'AU':
        from numpy import pi
        G = 4*pi**2 # AU^3*yr^-2M_sol^-1
        if M > 100:
            print 'M needs to be in Solar masses'
            M = M/M_sol
        elif R > 1:
            print 'R needs to be in AU'
            R = R/AU
        v_stab = sqrt((M*G)/R)
        print """
stable orbit velocity is: %g AU/yr for planet %g
              """ %(v_stab, nr)
    else: 
        print 'call needs dim to be either "SI" or "AU"'

if __name__ == '__main__':
    from numpy import pi
    from MySolarSystem import *
    syst = Myseed()
    nr = 1 # which planet is interesting
    M = p_mass(syst, nr)
    R = p_radius(syst, nr)
    stable(M, R, nr)
    stable(M,R,nr, 'AU')
    stable(M,R,nr,3)
