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

# constants
G = 6.674e-11 #N * m^2 * kg^-2 gravitational constant
AU = 149597871e3 #m in an AU
M_sol = 1.9891e30 #kg in 1 solar mass
m = 1100 #kg mass of satellite

#from area_vterminal
A = 266 #m^2 area of parachute. 

#density and temperature
path = '/home/anderger/Documents/ast1100-project/atmosphere_val/'
T = np.load(path+'Temperature.npy')# Temperature array as a function of distance
rho = np.load(path+'Density.npy') # temperature array as func of distance
r_dense = np.load(path+'Position.npy') # position vector for density
rho_func = scp.interpolate.interp1d(r_dense, rho)

# radius, mass from AST1100solarsystem
myss = MMS.Myseed()
R = MMS.p_radius(myss, 1) # radius of planet
M = MMS.p_mass(myss, 1) # mass of planet
R *= 1e3 # converting R to m
M *= M_sol # converting mass to kg

# initial posittion relative to planet
r0 = np.array(( 134509692.386 ,  0.0 ,  0))
print 'r0', type(r0), r0
# initial velocity relative to planet
v0 = np.array((0.00786262094791 ,  2985.69804027 ,  0))
print 'v0', type(v0), v0

nr = 100000
r = np.zeros((nr, 3))
v = np.zeros((nr, 3))
r[0] = r0
v[0] = v0
print 'r0', r[0]
print 'v0', v[0]

#Leapfrog integration
t0 = 1 #s after init call, 
dt = 60#s timestep
t = 0
first = 0
second = 0
# v0 to vhalf:
r_ = np.linalg.norm(r[0]) 
a = -G* ( M/r_**3 ) * r[0]
v[0] = v[0] + 0.5*a*dt
print 'initial velocity' 
print v[0]

print 'integrating'
print '===================================='
# integrater
for i in range(nr-1):
    Dv = 0
    # checking for nan type objects
    if any(np.isnan(r[i])):
        print 'error in pos'
        print 'here is the error, at time ', t, 'index ', i
        print 'at a distance ', r[i+1]
        break
    if any(np.isnan(v[i])):
        print 'error in vel'
        print 'here is the error, at time ', t, 'index ', i
        print 'at a distance ', r[i+1]
        break

    #updating position
    r[i+1] = r[i] + v[i]*dt

    #position reliant variables
    r_ = np.linalg.norm(r[i+1])
    v_ = np.linalg.norm(v[i])
    ag = -G*r[i+1]*M/r_**3

    # checks to include drag
    if r_ < r_dense[-1]:
        ad = (-0.5*rho[i]*v[i]*A*v_**3)/m
        ratio = np.linalg.norm(ag)/np.linalg.norm(ad)
        if ratio < 1000:
            if first < 1:
                print 'too low'
                first = 1
            a = ag + ad
        else: 
            if first < 1:
                print 'still high enough'
                first = 1
            a = ag
    else: 
        if second < 1:
            print 'we use only gravity'
            print ag
            second = 1
        a = ag

    # checks to boost 
    e_v = (v[i]/v_)
    if -5e4 < r[i+1,1] < 5e4 and 0 < r[i+1,0]: 
        dv = -0.3*v_
        #th = np.arctan2(v[i,0], v[i,1])
        #Dv = np.zeros(3)
        #Dv[0] = dv*np.cos(th)
        #Dv[1] = dv*np.sin(th)
        Dv = dv*e_v
    elif -5e4 < r[i+1,1] < 5e4 and r[i+1,0] < 0: 
        v_so = np.sqrt( (M * G) / (r_) )
        dv = v_so - v_
        #th = np.arctan2(v[i,0], v[i,1])
        #Dv = np.zeros(3)
        #Dv[0] = dv*np.cos(th)
        #Dv[1] = dv*np.sin(th)
        Dv = dv*(v[i]/v_)

    # estimating new velocity
    v[i+1] = v[i] + a*dt + Dv

    # updating times and indexes to record last entry.
    end = i
    t += dt

print 'time passed', t, 'seconds'
print 'at this time, index ', end, ', position and velocity are:'
print 'r', r[end]
print 'v', v[end]

## making a 3D plot of orbital path
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")
#hold('on')

#initial position of satelite:
print 'now to plot the shit'
#ax.plot([r[end,0]], [r[end,1]], [r[end,2]], color='g', marker='o', label='end')
ax.plot([r[0,0]], [r[0,1]], [r[0,2]], color='g', marker='o', label='beginning')

##plotting a sphere for the planet.
theta, phi = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = R*np.cos(theta)*np.sin(phi)
y = R*np.sin(theta)*np.sin(phi)
z = R*np.cos(phi)
ax.plot_wireframe(x, y, z, color="r")

#plotting satellite position
ax.plot(r[:,0], r[:,1], r[:,2], label='satellite')

# plotting vectors for points of boost
#ax.quiver(v_boost)

ax.legend()
plt.show()


"""
for testing: 

    print 'at time ', t, ', index ', i, ' we need send in:'
    print 'r(',t,'): ', r[i]
    print 'v: ', v[i]
    print 'dt: ', dt
    print 'norm r', np.linalg.norm(r[i])
    print 'Dv: ', v[i]*dt


    print 'r[i+1]', r[i+1]
    print 'new norm', r_

    print 'with grav. acceleration', ag

    print 'acceleration is: ', a
"""
