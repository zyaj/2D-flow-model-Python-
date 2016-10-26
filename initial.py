# inital.py
# provide the initial condition for the 2D flow simulation
# It will set initial water depth, free surface, salinity
# and velocity field
# Yun Zhang
# Stanford University
# ver 1.0 10/10/2014

import numpy as np

# set initial water depth
def ReturnDepth(x):
	return 1

# set initial free surface
def ReturnFreeSurface(x):
	return 0

# set inital salinity field
def ReturnScalar(x):
	return 0

# set initial velocity field
def ReturnHorizontalVelocity(x):
	return 0

# set veritical velocity field
def ReturnVerticalVelocity(x):
	return 0


# set bottom roughness 
def ReturnBottomRoughness(x):
	return 0.00001

