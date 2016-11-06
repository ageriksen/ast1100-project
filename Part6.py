
from AST1100ShortcutSystem import AST1100SolarSystem

infile = open('myseed.txt', "rf")
seed = int(infile.readline())
myss = AST1100SolarSystem(seed) 

# in the interest of time, chi-squared of gasses is 
# more or less omitted. thus. 
mu = 30

#myss.landOnPlanet(1, 'lander.txt')
