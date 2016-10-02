#import here... 
#use AST1100SolarSystemViewer for standard exercises 
#use AST1100SolarSystem for alternative project

system = AST1100SolarSystem(seed)
#system = AST1100SolarSystem(seed,  hasMoons=False) #if the moons are problematic 

T = #number of years 
N =  #steps per year

times = linspace(0, T, int(N*T) )
pos_computed = zeros((2, system.numberOfPlanets, int(T*N)))

#compute the positions here! 
#

#use this line to perform check of your simulations. 
#if passed, will generate a numpy array file with perfect positions. 
system.checkPlanetPositions(pos_computed, T, N)


#if the file is too large, use something like this to only output 1/100 frames
#system.orbitXml(pos_computed[:,:,::100], times[::100]) 

system.orbitXml(pos_computed, times)