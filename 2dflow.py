import grid
import phys
import initial
import boundary
import numpy as np
import matplotlib.pyplot as plt
import time
import scalar as sc
import sys

if __name__ == "__main__":
	try:
		if sys.argv[1]=='default':
			L=10000;
			nc=100;
			nk=10;
			dt=1;
			theta=0.55;
			ntime=1350;
		else:
			L=sys.argv[1]
			nc=sys.argv[2]
			nk=sys.argv[3]
			dt=sys.argv[4]
			theta=sys.argv[6]
			ntime=sys.argv[5]
		xc,xe,d,zb,zt,Nk,dz=grid.GridGeneration(L,nc,nk)
		h,u,w,s,cdb=phys.Solve(ntime,dt,nc,xc,xe,d,zb,zt,nk,Nk,dz,theta,1)

	except (ValueError, IndexError):
		print"Error in command line input"
		print"Run as: python 2dflow.py <L> <nc> <nk> <dt> <theta> <ntime>"
		print"L: channel length; nc: number of horizontal cells; nk: number of vertical layers"
		print"dt: time step size; ntime: number of time steps"
		print"theta: theta value for theta method. 1 for fully implicit (recommanded value 0.55)"
		print"Also run as: python 2dflow.py default"
		print"The program will run the default problem shown in the write up"
		

