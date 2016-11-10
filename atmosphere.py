"""
script to simulate atmosphere density
and temperature of planet
"""
import numpy as np
import MySolarSystem as MSS
import matplotlib.pyplot as plt

np.seterr(all='raise')

G = 6.674e-11 #N*m^2*kg^-2 grav. const. 
mu = 30 # an estimate of mean molecular mass
u = 1.660539040e-27 # kg, rough mass per nucleon
mH = 1.007825*u# hydrogen mass
k = 1.38064852e-23#J*K^-1 Boltzmann's constant
gamma = 1.4 # adiabatic power
M_sol = 1.98855e30 #kg in 1 solar mass

myss = MSS.Myseed()

M = MSS.p_mass(myss, 1) * M_sol # converting to kg
R = MSS.p_radius(myss, 1)*1e3

nr = 1000000 #1000000 #1e6 number of steps
r = np.linspace(R,R+600e3, nr) #m | distance from surface to end 
#dr = 200e3/nr
P = np.zeros(nr)
T = np.zeros(nr)
rho = np.zeros(nr)
T[0] = 281 #K temperature at surface. Read from ssview. estimates lie at ~279K
P[0] = myss.rho0[1] # Atmospheric density at surface
adia = P[0]**(1-gamma)*T[0]**gamma #adiabatic constant



for i in range(nr-1):
    if i == 1:
        print i 
        print 'P[i]', type(P[i]), P[i]
        print 'T[i]', type(T[i]), T[i]
        print 'r[i]', type(r[i]), r[i]
        print 'M', type(M), M
    dr = r[i+1]-r[i]
    P[i+1] = P[i] + ((P[i] * M) / (T[i] * r[i]**2))*( -4 * np.pi * G * (mu * mH / k)**2 )*dr
    if T[i] > T[0]*0.5:
        T[i+1] = (adia/(P[i+1]**(1-gamma)))**(1/gamma)
    else: 
        T[i+1] = T[i]
    M = M + (dr**2*P[i+1]/T[i+1])*dr
    ii = i

print ii 
print 'P[ii]', type(P[ii]), P[ii]
print 'T[ii]', type(T[ii]), T[ii]
print 'r[ii]', type(r[ii]), r[ii]
print 'M', type(M), M
rho = P*mu*mH/k/T
print 'rho[ii]', type(rho[ii]), rho[ii]

r_ = r-R
r_ = r_*1e-3
plt.plot(r_, rho)
plt.title('atmospheric density')
plt.xlabel('km above surface')
plt.ylabel('N*m^-2')
plt.show()

# saving to arrays:
path = '/home/anderger/Documents/ast1100-project/atmosphere_val'
np.save(path + '/Temperature.npy', T)
np.save(path + '/Density.npy', rho)
np.save(path + '/Position.npy', r)

"""
run example:
    homeplanet is planet 0
    mass of planets, [solar masses]
    radi of planets, [km]
    1
    P[i] <type 'numpy.float64'> 1.40867059095
    T[i] <type 'numpy.float64'> 280.999607365
    r[i] <type 'numpy.float64'> 9327652.44279
    M <type 'numpy.float64'> 1.79683818763e+25
    999998
    P[ii] <type 'numpy.float64'> 0.00042087020128
    T[ii] <type 'numpy.float64'> 140.499938588
    r[ii] <type 'numpy.float64'> 9927651.24279
    M <type 'numpy.float64'> 1.79683818763e+25
    rho[ii] <type 'numpy.float64'> 1.0892921753e-08
"""
