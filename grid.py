# grid.py
# generate grid information
# which will be used in following calculation
# Yun Zhang
# Stanford University
# ver1.0 10/10/2014

import numpy as np
import initial as ic

# main grid generation function
def GridGeneration(L,nc,nk=1):
	# section 1 generate horizontal grid
	# L=river length, nc=number of cells
	# xc=cell center xe=edge center
	# d=cell center depth
	dx=float(L)/nc
	xc=np.linspace(dx/2,(L-dx/2),nc)
	xc=xc.reshape((nc,1))
	xe=np.linspace(0,L,(nc+1))
	xe=xe.reshape((nc+1,1))
        d=np.zeros(np.shape(xc))
	for i in range(nc):
		d[i]=ic.ReturnDepth(xc[i])
        # section 2 generate vertical grid
        globalzb,globalzt,Nk,dz=InitialVerticalGrid(xc,d,nk)
        
        return xc,xe,d,globalzb,globalzt,Nk,dz

# initially generate vertical grid for xc
# set the active layer for each cell
# set the bottom cell height for each cell
def InitialVerticalGrid(xc,d,nk):
	dmax=d.max()
	nc=len(xc)
	dz=dmax/nk
	# zt=cell top elevation
	zt=np.linspace(0,(dmax-dz),nk)
	zt=zt.reshape((nk,1))
	zt=zt.T
	zt=-zt
	# zb=cell bottom elevation
        zb=np.linspace(dz,dmax,nk);
	zb=zb.reshape((nk,1))
	zb=zb.T
	zb=-zb
	Z=0.5*(zb+zt)
	# the number of active layer for each cell
	Nk=np.zeros(np.shape(xc))
	for i in range(len(xc)):
		for k in range(nk):
			if zb[0][k]<=-d[i]:
				Nk[i]=k+1
				break
	return zb,zt,Nk,dz


