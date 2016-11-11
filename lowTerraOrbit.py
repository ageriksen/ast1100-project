"""
script to propell the satellite to a low orbit 
to capture images of the planet. 
initial values from AST1100SolarSystem

Initial relative position = ( 0.000899197293603 ,  0.0 ,  0 )
Initial relative velocity = ( 1.65867319879e-06 ,  0.629853245105 ,  0 )
----------------------------------------------------------
Initial relative position (SI units) = ( 134509692.386 ,  0.0 ,  0 )
Initial relative velocity (SI units) = ( 0.00786262094791 ,  2985.69804027 ,  0 )
Satellite is now receiving commands. This might take some minutes...
Automatic orientation. Time: 1.0
Pos:  [  1.34509692e+08   2.98569804e+03   0.00000000e+00]
Vel:  [ -5.84332194e-02   2.98569804e+03   0.00000000e+00]
Automatic orientation. Time: 2.0
Pos:  [  1.34509692e+08   5.97139608e+03   0.00000000e+00]
Vel:  [ -1.24729060e-01   2.98569804e+03   0.00000000e+00]
No landing detected at time =  2.0 s. Exiting.
"""
import numpy as np
import MySolarSystem as MMS
import scipy.interpolate as scp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def v_s(): #stable orbit velocity
    return np.sqrt( (M * G) / Rs )


def e_theta(): 
    # function for cartesian coordinates of polar
    theta = theta + np.arctan2()
    phi = phi + np.arctan2()

G = 6.674e-11 #N * m^2 * kg^-2 gravitational constant
AU = 149597871e3 #m in an AU
M_sol = 1.9891e30 #kg in 1 solar mass
m = 1100 #kg mass of satellite

path = '/home/anderger/Documents/ast1100-project/atmosphere_val/'
T = np.load(path+'Temperature.npy')# Temperature array as a function of distance
rho = np.load(path+'Density.npy') # temperature array as func of distance
r_dense = np.load(path+'Position.npy') # position vector for density
rho_func = scp.interpolate.interp1d(r_dense, rho)

t0 = 1 #s after init call, 
dt = 50#s timestep

myss = MMS.Myseed()
R = MMS.p_radius(myss, 1) # radius of planet
M = MMS.p_mass(myss, 1) # mass of planet
R *= 1e3 # converting R to m
M *= M_sol # converting mass to kg

# initial posittion relative to planet
r0 = np.array((1.34509692e+08, 2.98569804e+03, 0.00000000e+00))
print 'r0', type(r0), r0
# initial velocity relative to planet
v0 = np.array((-5.84332194e-02, 2.98569804e+03, 0.00000000e+00))
print 'v0', type(v0), v0

nr = 1000000#1e6
r = np.zeros((nr, 3))
v = np.zeros((nr, 3))
r[0] = r0
v[0] = v0
print 'r0', r[0]
print 'v0', v[0]

#Leapfrog integration
# v0 to vhalf:
r_ = np.linalg.norm(r) 
a = -G* ( M/r_**3 ) * r[0]
v[0] = v[0] + 0.5*a*dt

# integrater
t = 0
first = 0
second = 0
for i in range(nr-1):
    r[i+1] = r[i] + v[i]*dt
    r_ = np.linalg.norm(r[i+1])
    ag = -G*r[i+1]*M/r_**3
    if r_ < r_dense[-1]:
        if first < 1:
            print 'pretty low'
        first = 1
        rh = rho_func(r[i+1])
        v_ = np.linalg.norm(v[i])
        ad = -0.5*rh*v[i]*A*v_**3
        ratio = np.linalg.norm(ag)/np.linalg.norm(Fdm)
        if ratio < 1000:
            a = ag + ad
        else: 
            a = ag
    else: 
        if second < 1:
            print 'we use only gravity'
            print ag
        second = 1
        a = ag
    v[i+1] = v[i] + a*dt
    end = i
    t += dt

print 'time passed', t, 'seconds'
print 'at this time, position and velocity are:'
print 'r', r[end]
print 'v', v[end]

## making a 3D plot of orbital path
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(r[:,0], r[:,1], r[:,2], label='satellite')
ax.legend()
plt.show()
