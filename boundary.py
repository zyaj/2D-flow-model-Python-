# boundary.py
# setup h and u boundary
# and s boundary value
# Yun Zhang 
# Stanford University
# ver 1.0 10/10/2014
import numpy as np

# set up free surface boundary
def BoundaryFreeSurface(n,dt):
	amp=1
	T=3600
	t=(n+1)*dt
	return amp*np.sin(2*3.1415926/T*t)

# set up velocity boundary
def BoundaryVelocity(x):
	return 0

# set up upstream scalar boundary
def UpstreamBoundaryScalar(x=0):
	return 0

# set up downstream scalar boundary for tidal bc
def DownstreamBoundaryScalar(x=0):
	return 1
