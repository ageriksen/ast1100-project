"""
"""
####################################################################################
    #imports
from numpy import(
    linalg, pi, zeros, linspace, load, array, sum)
from matplotlib.pyplot import(
    plot, xlabel, ylabel, title, 
    hold, show, axis, legend, figure)
import MySolarSystem as M
from LeapFrog import *

####################################################################################
    #setting funcitons
def accelerate(r,  m):
    """
    Acceleration from gravitational forces on object from 
    the 3 largest planets in the solar system and the star
    """
    r_ = linalg.norm(r)
    return -(G*m*r)/(r_**3)   
    
####################################################################################
    #constants
G = 4*pi**2 # AU**3 * yr**-2 * M_s**-1 | gravitational const.
T = 25 # yrs | time span of simulation
n = 20000 # timesteps per year
dt = 1./n # time per timestep

####################################################################################
    #gathering solar system information
syst = M.Myseed() # instancing AST1100SS wit my seed. 
M_s = M.starmass(syst) # starmass in solar masses
N = M.Nplanets(syst) # number of planets in solar system
r0 = [ [syst.x0[i], syst.y0[i] ] for i in range(N) ]
v0 = [ [syst.vx0[i], syst.vy0[i] ] for i in range(N) ]

####################################################################################
    #defining arrays
r = zeros((n*T,N,2)) # initial pos. for planets in xy-plane
v = zeros((n*T,N,2)) # initial velocity. 

####################################################################################
    #filling arrays
r[0] = r0
v[0] = v0

#half-step velocity:
for j in range(N):
    v[0,j,:] = v0_half(v[0,j,:], accelerate(r[0,j,:], M_s),dt)


####################################################################################
    #integrating

for i in range(T*n-1):
    for j in range(N):
        r[i+1,j,:], v[i+1, j, :] =  LeapFrog(
            r[i,j,:], v[i,j,:], accelerate, dt, M_s)

####################################################################################
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


####################################################################################
    # checking results with AST1100SolarSystem and making movie file and position
    # array.
#to fit array into prescripted code, swap axes.
r = r.swapaxes(0,2)
v = v.swapaxes(0,2)
# making time-array t for functions
t = linspace(0,T, n*T)
# calling prescripted functions to check planetpositons and
# write a file
syst.checkPlanetPositions(r[:,:],float(T),n)
#syst.orbitXml(r[:,:,::100], t[::100])
