
from AST1100ShortcutSystem import AST1100SolarSystem

infile = open('myseed.txt', "rf")
seed = int(infile.readline())
myss = AST1100SolarSystem(seed) 
myss.landOnPlanet(1, 'lander.txt')
#try: update AST1100SolarSystem.

