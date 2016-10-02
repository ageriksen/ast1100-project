"""
This module aims to generate a miniature rocket engine and simulate 1 nanosecond activity.
"""

def generate(low, high, nr, T, m_p, seed='NaN'):
    from numpy import random, sqrt
    # function aims to generate arrays of n rows of particles and 3 dimensional collumns
    # inside a box of length L on each side.
    if type(seed) != 'NaN':
        random.seed(100)
    # There is an equal chance of finding any one particle anywhere in the box, so 
    # the position vector is a uniform probability.
    r = random.uniform(low,high,(nr,3))

    # gaussian velocity. random directions means mean value mu is 0, sigma o is 
    # calculated through a boltzsmann distribution. This results in a Gaussian distribution. 
    # m_p is the mass of a particle, k is the boltzmann constant and T is the temperature in
    # Kelvin. 
    k = 1.38064852e-23 # J*K^-1 Boltzamann constant  
    mu = 0 # average velocity.
    o = sqrt((k*T)/m_p) # sigma, the standard deviation
    v = random.normal(mu,o,(nr,3))
    return r, v

def test_distribution(r,v, figname_pos=0,figname_vel=0):
    from matplotlib.pyplot import plot, show, figure, bar, title, savefig 
    from numpy import histogram
    hist, bins = histogram(r, 'auto')
    width = 0.7*(bins[1] - bins[0])
    center = (bins[:-1]+bins[1:])/2
    bar(center, hist, align='center', width=width)
    title('position')
    if figname_pos != 0:
        savefig(figname_pos)
    show()
    hist, bins = histogram(v, 'auto')
    width = 0.7*(bins[1] - bins[0])
    center = (bins[:-1]+bins[1:])/2
    bar(center, hist, align='center', width=width)
    title('velocity')
    if figname_vel != 0:
        savefig(figname_vel)
    show()


def test_values(v, N,T, m_p):
    from numpy import sum, sqrt, pi
    k = 1.38064852e-23 # J*K^-1 Boltzamann constant  
    Ek_ana = (3./2)*N*k*T # analytical kinetic energy
    Ek_num = (0.5*m_p)*sum( v**2 ) # numerical kinetic energy
    print '---------------------------------------------------------'
    print 'kinetic energy'
    print 'numerical', Ek_num
    print 'analytical', Ek_ana
    print 'relative error', abs( (Ek_num - Ek_ana)/Ek_ana )
    print '---------------------------------------------------------'

    absv_ana = 2.*N*sqrt( (2.*k*T)/(pi*m_p) ) # scaled upaverage absolute velocity in gas.
    absv_num = sum( sqrt( sum(v**2, axis=1) ) )
    print 'velocity'
    print 'analytical', absv_ana
    print 'numerical', absv_num
    print 'relative error', abs( (Ek_num - Ek_ana)/Ek_ana )
    print '---------------------------------------------------------'


def simulate(Dt, n, r, v, L, m_p, l_s, l_e):
    from numpy import where, nonzero
    # function to simulate rocket boost produced by small rocket engine. 
    # The loop checks whether or not a particle has crashed with 
    dt = Dt/n
    p_out = 0
    n_out = 0
    p_z = 0
    for timestep in range(n):
        for xyz in range(3):
            roof = where(r[:,xyz] > L, 1,0)
            floor = where(r[:,xyz] < 0, 1,0)
            crash = roof+floor
            if xyz == 2:
                # placing hole in floor of z
                n_x = where(
                (r[nonzero(floor),0]>l_s) & (r[nonzero(floor),0]<l_e), 1,0)
                n_y = where(
                (r[nonzero(floor),1]>l_s) & (r[nonzero(floor),1]<l_e), 1,0)
                n_pout = n_x*n_y
                n_pin = where(n_x*n_y > 0, 0,1)
		print 'n_pout'
		print n_pout
		print 'nonzero'
		print nonzero(n_pout)
		print 'v[:,2]'
		print v[:,2]
                p_out += 2.*m_p*sum( 
		abs( 
		v[nonzero( n_pout ), 2] 
		) 
		)
                n_out += sum(n_pout)
                r[nonzero(n_pout), xyz] += L
                r[nonzero(n_pin),xyz] = - r[nonzero(n_pin),xyz]
                r[nonzero(roof), xyz] = 2*L - r[nonzero(roof),xyz]
                v[nonzero(n_pin),xyz] = -v[nonzero(n_pin),xyz]
            else:
                v[nonzero(crash),xyz] = -v[nonzero(crash),xyz]
                r[nonzero(roof), xyz] = 2*L - r[nonzero(roof),xyz]
                r[nonzero(floor), xyz] = -r[nonzero(floor), xyz]
        p_z += 2.*m_p*sum(abs(v[nonzero(roof),2]))
        r += v*dt
    print '---------------------------------------------------------'
    print 'particles out', n_out
    print 'momentumloss', p_out
    print '---------------------------------------------------------'
    return p_out, n_out
#loop uses: n, r, L, l_s, l_e, m_p


def pressure(L, p, Dt, N, T):


    k = 1.38064852e-23 # J*K^-1 Boltzamann constant  


    print '---------------------------------------------------------'
    print 'comparing numerical pressure and analytical'
    A = L**2 # m^2
    P_num = p/(A*Dt) # Pa, pascal = Kg*m^-1*s^-2
    print 'numerical', p_num

    V = L**3
    P_ana = (N*k*T)/(V)
    print 'analytical', p_ana

    print 'relative error:', (abs(p_num - p_ana) / p_ana)
    print '---------------------------------------------------------'

def velocityChange(p, Dt, DT, m):
    dp = p/Dt
    Dp = dp*DT
    Dv = Dp/m
    print '---------------------------------------------------------'
    print 'speed gain of mass per box', Dv, 'm/s'
    print '---------------------------------------------------------'
    return Dp, Dv


def fuelConsumption(n, m_p, Dt, DT, n_box):
    dm = n*m_p/Dt
    Dm_box = dm*DT
    Dm = Dm_box*n_box
    print '---------------------------------------------------------'
    print (
    'fuel consumption of mass over time DT ' + str(DT) + 
    '\n for nr of engines ' + str(n_box) + ': \n', Dm
    )
    print '---------------------------------------------------------'




"""
version 0.1 | 180916
note: finish documentation
"""
