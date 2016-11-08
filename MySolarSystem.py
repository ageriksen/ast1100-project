"""
this class aims to collect and distribute information of my solar system.
"""

def Myseed():
	from AST1100SolarSystem import AST1100SolarSystem
	with open('myseed.txt', 'r') as f:
		myseed = int(f.readline())
	#check to see if myseed is imported as correct object.
	print 'myseed', myseed, type(myseed)
	print 'homeplanet is planet 0'
	myStarSystem = AST1100SolarSystem(myseed)
	return myStarSystem
	
def starmass(starsyst):
	print 'mass of star, [solar masses]'
	return starsyst.starMass

def starradius(starsyst):
	print 'radius of star, [km]'
	return	starsyst.starRadius

def Nplanets(starsyst):
	N = starsyst.numberOfPlanets
	print 'number of planets', N
	return N	

def startemp(starsyst):
	print 'surface temperature of star, [K]'
	return starsyst.temperature

def semimajor(starsyst, nr='NaN'):
	print 'semi-masjor axis of planets, [AU]'
	if nr == 'NaN':
		print 'no index given'
		return starsyst.a
	else: 
		return starsyst.a[nr]

def eccentricity(starsyst, nr='NaN'):
	print 'eccentricity of planets'
	if nr == 'NaN':
		print 'no index given'
		return starsyst.e
	else:
		return starsyst.e[nr]

def init_angle(starsyst, nr='NaN'):
	if nr == 'NaN':
		print 'no index given'
		w = starsyst.omega
		print "initial angle of planets' orbits, rad:"
		print w
		return w
	else:
		w = starsyst.omega[nr]
		print "initial angle of planets' orbits, rad:"
		print w
		return w

def axisangle(starsyst, nr='NaN'):
	if nr == 'NaN':
		print 'no index given'
		psi = starsyst.psi
		print  'angle of semi-major axes for planets:'
		print psi
		return psi
	else:
		psi = starsyst.psi[nr]
		print  'angle of semi-major axes for planets:'
		print psi
		return psi


def p_mass(starsyst, nr='NaN'):
	print 'mass of planets, [solar masses]'
	if nr == 'NaN':
		print 'no index given'
		return starsyst.mass
	else:
		return starsyst.mass[nr]

def p_radius(starsyst, nr='NaN'):
	print 'radi of planets, [km]'
	if nr =='NaN':
		return starsyst.radius
	else:
		return starsyst.radius[nr]

def period(starsyst, nr='NaN'):
	print 'rotational period of planets, [earth days]'
	if nr == 'NaN':
		print 'no index given'
		return starsyst.period 
	else: 
		return starsyst.period[nr]

def x0(starsyst, nr='NaN'):
	print 'initial x-position of planets, [AU]'
	if nr == 'NaN':
		print 'no index given'
		return starsyst.x0
	else:
		return starsyst.x0[nr]

def y0(starsyst, nr='NaN'):
	print 'initial y-position of planets, [AU]'
	if nr == 'NaN':
		print 'no index given'
		return starsyst.y0
	else:
		return starsyst.y0[nr]

def vx0(starsyst, nr='NaN'):
	print 'initial x-velocityof planets, [AU/Yr]'
	if nr == 'NaN':
		print 'no index given'
		return starsyst.vx0
	else:
		return starsyst.vx0[nr]

def vy0(starsyst, nr='NaN'):
	print 'initial y-velocityof planets, [AU/Yr]'
	if nr == 'NaN':
		print 'no index given'
		return starsyst.vy0
	else:
		return starsyst.vy0[nr]

def atmosphere(starsyst, nr='NaN'):
	print 'Atmospheric density at surface, [kg/m^3]'
	if nr == 'NaN':
		print 'no index given'
		return starsyst.rho0
	else:
		return starsyst.rho0[nr]



if __name__ == '__main__':
	M = Myseed()
	print starmass(M)
	print starradius(M)
	print Nplanets(M)
	print startemp(M)
	print semimajor(M)
	print eccentricity(M)
	init_angle(M)
	axisangle(M)
	p_mass(M,0)
	p_mass(M)
	p_radius(M,0)
	p_radius(M)
	print period(M)
	print x0(M)
	print y0(M)
	print vx0(M)
	print vy0(M)
	print atmosphere(M)



"""
notes:
	should change the prints in case of nr != NaN to include whitespace
version 1.0
edit history:
added edit history, version number and notes | 140916(ddmmyy)
"""
