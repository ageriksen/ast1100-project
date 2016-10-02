"""
Aims to produce positions of planets in Solar system w/ only star as relevant 
and in origin. 9
"""

from numpy import(
linalg, pi, zeros, linspace)
from matplotlib.pyplot import(
plot, xlabel, ylabel, title, hold, show, axis, legend)
import MySolarSystem as M
from LeapFrog import *

G = 4.*pi**2 # AU-based gravitational constant.
T = 25 #yrs of sim.
n = 20000 # timesteps per year
dt = 1./n # time per timestep in yrs
syst = M.Myseed() # retrieving seed and establishinf SolarSystem instance
m = M.starmass(syst) # mass of star
N = M.Nplanets(syst) # number of planets to calculate
r = zeros((n*T,N,2)) # positional vector for timesteps, planets in xy-plane
v = zeros((n*T,N,2)) # velocity vector

r[0] = [
[ syst.x0[i], syst.y0[i] ] for i in range(N) ]
v[0] = [
[ syst.vx0[i], syst.vy0[i] ] for i in range(N) ]

# definint acceleration function
def acceleration(r, G, m):
        from numpy import linalg
        r_ = linalg.norm(r)
        a = -(G*m*r)/(r_**3)
        return a

# advancing v a half-step.
for j in range(N): # need to step each planet
    v[0,j,:] = v0_half(
        v[0,j,:],acceleration(r[0,j,:], G, m), dt)
# integrating planet positions through LeapFrog method
for i in range(n*T-1):
    for j in range(N):
        r[i+1,j,:], v[i+1,j,:] = LeapFrog(
        r[i,j,:], v[i,j,:], acceleration, dt, var=[G, m])

# ploting position values
plot(r[:,0,0],r[:,0,1], label=('planet 0'))
plot(r[0,0,0],r[0,0,1], 'o', label=('r0 planet 0'))
hold('on')
axis('equal')
for nr in range(1,N):
    plot(r[:,nr,0],r[:,nr,1], label=('planet '+str(nr)))
    plot(r[0,nr,0],r[0,nr,1],'o', label=('r0 planet '+str(nr)))
legend()
show()

#to fit array into prescripted code, swap axes.
r = r.swapaxes(0,2)
v = v.swapaxes(0,2)
# making time-array t for functions
t = linspace(0,T, n*T)
# calling prescripted functions to check planetpositons and
# write a file
syst.checkPlanetPositions(r[:,:,::100],float(T),n/100.)
syst.orbitXml(r[:,:,::100], t[::100])
