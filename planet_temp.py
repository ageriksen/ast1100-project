"""
Short code to calculate and compare read surface
temperature of planet 1
"""
import MySolarSystem as M
import math
myss = M.Myseed()
T_s = M.startemp(myss)
R_s = M.starradius(myss)
AU = 149597871# km in 1 AU
r_px = M.x0(myss, 1)*AU
r_py = M.y0(myss, 1)*AU
r_p = math.sqrt(r_px**2 + r_py**2)
T_p = math.sqrt(R_s/(2*r_p))*T_s
T_p_read = 281 #k

print 'surface temperature of planet 1 is'
print 'calculated: ', T_p
print 'previously measured: ', T_p_read
print 'with a relative error: ', abs(T_p - T_p_read)/T_p_read

f = open('planet_temp.txt', 'r+')
f.write(str(T_p_read))
f.seek(0)
f.readline(1)
f.close()
"""
run example:

myseed 80101 <type 'int'>
homeplanet is planet 0
surface temperature of star, [K]
radius of star, [km]
initial x-position of planets, [AU]
initial y-position of planets, [AU]
surface temperature of planet 1 is
calculated:  279.869889372
previously measured:  281 
with a relative error:  0.00402174600664

"""
