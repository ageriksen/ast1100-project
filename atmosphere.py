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
k = 3.37e-27 # planck's constant
gamma = 1.4 # adiabatic power
M_sol = 1.98855e30 #kg in 1 solar mass

myss = MSS.Myseed()
#homeplanet is planet 0
#mass of planets, [solar masses]
#radi of planets, [km]
#P[0] 1.40867748005
#T[0] 281.0
#r[0] 9327.65184279
#M 9.03592158922e-06

M = MSS.p_mass(myss, 1) * M_sol # converting to kg
R = MSS.p_radius(myss, 1)*1e3

nr = 1000000 #1e6 number of steps
r = np.linspace(R,200e3, nr) #m | distance from surface to end 
P = np.zeros(nr)
T = np.zeros(nr)
rho = np.zeros(nr)
T[0] = 281 #K temperature at surface. Read from ssview. estimates lie at ~279K
P[0] = myss.rho0[1] # Atmospheric density at surface
adia = P[0]**(1-gamma)*T[0]**gamma #adiabatic constant


print 'P[0]', P[0]
print 'T[0]', T[0]
print 'r[0]', r[0]
print 'M', M

for i in range(nr-1):
    print i # stops at 19
    print 'P[i]', type(P[i]), P[i]
    print 'T[i]', type(T[i]), T[i]
    print 'r[i]', type(r[i]), r[i]
    print 'M', type(M), M
    P[i+1] = P[i]*M/(T[i]*r[i]**2)
    P[i+1]*(-4*np.pi*G*mu**2*mH**2/k**2)
    if T[i] > T[0]*0.5:
        T[i+1] = (adia/(P[i+1]**(1-gamma)))**(1/gamma)
    else: 
        T[i+1] = T[i]
    M += (r[i+1]-r[i])**2*(P[i+1]/T[i+1])

rho = P*mu*mH/k/T

plt.plot(r, rho)
plt.title('atmospheric density')
plt.show()
