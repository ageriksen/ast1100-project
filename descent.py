"""
This script is to guide the sattelite down to a low-terran orbit. 
"""
from AST1100ShortcutSystem import AST1100SolarSystem

infile = open('myseed.txt', "rf")
seed = int(infile.readline())
myss = AST1100SolarSystem(seed) 

myss.landOnPlanet(1, 'lander.txt')
