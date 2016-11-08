"""
script to simulate atmosphere density
and temperature of planet
"""
import numpy as np
import MySolarSystem as M
import matplotlib.pyplot as plt

np.seterr(all='raise')

G = 4*np.pi #placeholder grav. const. 
mu = 30 # an estimate of mean molecular mass
mH = 1e-31 # placeholder hydrogen mass
k = 3.37e-27 # planck's constant
gamma = 1.4 # adiabatic power

myss = M.Myseed()
M = myss.mass[1]
R = myss.radius[1]


nr = 1000000 #1e6 number of steps
r = np.linspace(0,200e3, nr)
P = np.zeros(nr)
T = np.zeros(nr)
rho = np.zeros(nr)
T[0] = 281 #K
P[0] = myss.rho0[1]
adia = P[0]**(1-gamma)*T[0]**gamma #adiabatic constant


print 'P[0]', P[0]
print 'T[0]', T[0]
print 'r[0]', r[0]
print 'M', M

for i in range(nr-1):
    M += r[i]**2*(P[i]/T[i])
    P[i+1] = P[i]/(T[i]*r[i]**2)*M
    if T[i] > T[0]:
        T[i+1] = adia/P[i+1]
    else: 
        T[i+1] = T[i]

P*(-4*np.pi*G*mu**2*mH**2/k**2)
rho = P*mu*mH/k/T

plt.plot(r, rho)
plt.title('atmospheric density')
plt.show()
