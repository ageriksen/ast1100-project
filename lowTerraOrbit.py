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

def set_axes_equal(ax):
        '''
        Make axes of 3D plot have equal scale so that spheres appear as spheres,
        cubes as cubes, etc..  This is one possible solution to Matplotlib's
        ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

        Input
        ax: a matplotlib axis, e.g., as output from plt.gca().
        '''

        x_limits = ax.get_xlim3d()
        y_limits = ax.get_ylim3d()
        z_limits = ax.get_zlim3d()

        x_range = abs(x_limits[1] - x_limits[0])
        x_middle = np.mean(x_limits)
        y_range = abs(y_limits[1] - y_limits[0])
        y_middle = np.mean(y_limits)
        z_range = abs(z_limits[1] - z_limits[0])
        z_middle = np.mean(z_limits)

        # The plot bounding box is a sphere in the sense of the infinity
        # norm, hence I call half the max range the plot radius.
        plot_radius = 0.5*max([x_range, y_range, z_range])

        ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
        ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
        ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


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
#rho = rho[::-1]
#r_dense = r_dense[::-1]
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

nr = 500000
r = np.zeros((nr, 3))
v = np.zeros((nr, 3))
r[0] = r0
v[0] = v0
print 'r0', r[0]
print 'v0', v[0]

#Leapfrog integration
dt = 10#s timestep
t = 1 # seconds after "init" call
r_min = 0

# v0 to vhalf:
r_ = np.linalg.norm(r[0]) 
a = -G* ( M/r_**3 ) * r[0]
v[0] = v[0] + 0.5*a*dt

print 'integrating'
print '===================================='
# integrater
for i in range( nr - 1 ):
    r[i+1] = r[i] + v[i]*dt
    r_ = np.linalg.norm(r[i+1])
    v_ = np.linalg.norm(v[i])
    a = -G*r[i+1]*M*r_**(-3)
    v_so = np.sqrt( (M*G)/(r_) )
    theta = np.arctan2(r[i+1, 0], r[i+1, 1])
    Dv = 0
    e_v = (v[i]/v_)
    # at 1 extreme, lower vel. 
    # 1/2 round later, return to stable 
    # orbit velocity. 
    if r_min == 0:
        if abs(theta) < np.pi/20:
            Dv = -0.01*v_*e_v
        elif abs(theta - np.pi) < np.pi/20:
            dv = v_so - v_
            Dv = dv*e_v
    elif r_min != 0:
        if abs(theta) < np.pi/20:
            rxv = np.cross(r[i+1], v[i])
            e_rxv = rxv/(np.linalg.norm(rxv))
            Dv = 0.01*v_*e_rxv
    
    if r_ <= r_dense[-1]:
        ii = np.nonzero(np.in1d(r_dense ,  r_))
        if type(ii[0]) == int:
            rh = rho[ii[0]]
            ad = 0.5*rh*A*v_**2*m**(-1)
            ratio = ag / ad
    else: 
        ratio = 10000


    if ratio <= 3000:
        r_min = r_
    elif ratio < 1000:
        print 'too low again.' 
        print 'index', i
        break


    if ratio > 3000:
        if abs(theta) < np.pi/12:
            Dv = -0.01*v_*e_v
        elif abs(theta - np.pi) < np.pi/8:
            print 'hey ho'
            dv = v_so - v_
            Dv = dv*e_v
        v[i+1] = v[i] + a*dt + Dv
    elif ratio <= 3000:
        if abs(theta) < np.pi/12:
           rxv = np.cross(r[i+1], v[i])
           print 'r x v: ', rxv
           e_rxv = rxv / np.linalg.norm(rxv)
           Dv = 0.01*v_*e_v
        v[i+1] = v[i] + a*dt + Dv
        v1_ = np.linalg.norm(v[i+1])
        if v1_ != v_so:
           dv = v_so - v1_
           e_v1 = v[i+1]/v1_
           v[i+1] = v[i+1] + dv*e_v1
    end = i
    t += dt

print 'time passed', t, 'seconds'
print 'at this time, index ', end, ', position and velocity are:'
print 'r', r[end]
print 'v', v[end]

absV = np.linalg.norm(v, axis=1)
absR = np.linalg.norm(r, axis=1)
print 'absR', np.shape(absR)

#======================================================================#
## making a 3D plot of orbital path
fig = plt.figure()
ax = fig.gca(projection='3d')
#hold('on')

#initial position of satelite:
ax.plot([r[end,0]], [r[end,1]], [r[end,2]], color='g', marker='o', label='end')
ax.plot([r[0,0]], [r[0,1]], [r[0,2]], color='g', marker='o', label='beginning')

##plotting a sphere for the planet.
theta, phi = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = R*np.cos(theta)*np.sin(phi)
y = R*np.sin(theta)*np.sin(phi)
z = R*np.cos(phi)
ax.plot_wireframe(x, y, z, color="r")

#plotting satellite position
ax.plot(r[:,0], r[:,1], r[:,2], label='satellite')

ax.legend()
set_axes_equal(ax)

# plotting 2D ===================
plt.figure()
plt.plot(absR, label='absolute distance')
plt.legend()

plt.figure()
plt.plot(absV, label='absolute velocity')
plt.legend()

plt.show()
