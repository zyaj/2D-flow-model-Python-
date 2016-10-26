# phys.py
# the main function to calculate velocity and
# free surface field for each time step
# also call the scalar function in scalar.py
# Yun Zhang
# Stanford University
# ver 1.0 10/10/2014

import numpy as np
import boundary as bc
import initial as ic
import time
import scalar as sc
import matplotlib.pyplot as plt

# the main function to solve every time step
def Solve(ntime,dt,nc,xc,xe,d,zb,zt,nk,Nk,dz,theta,plot):
	# initial physical variable
	[h,u,w,s,z0b]=InitialPhysicalVariables(xc,xe,d,dt,nk)
	# initial vertical grid setup
	[dzz,ctop]=UpdateDzz(xc,h,d,zb,zt,Nk,dz,nk)
        ntout=int(0.05*ntime)
	if ntout==0:
		ntout=1
	for n in range(ntime):
               
                if n%ntout==0:
			t=time.time()
		# set flux height
		[dzf,D,etop,Nke]=SetFluxHeight(dzz,ctop,Nk,u,h,nk)
		
		# set drag coefficient
		cdb=SetDragCoefficient(dzf,Nke,z0b,u)

		# calculate velocity and free surface
		[unew,hnew]=CalculateUandH(n,dt,xc,xe,dzf,D,dzz,Nk,Nke,etop,ctop,nk,u,h,cdb,theta)

		# update dzz and ctop
		[dzznew,ctopnew]=UpdateDzz(xc,hnew,d,zb,zt,Nk,dz,nk)

		# calculate vertical velocity
		wnew=CalculateVerticalVelocity(dt,xc,dzf,Nk,ctopnew,unew,nk)

		# scalar transport 
		snew=sc.ScalarTransport(dt,hnew,h,s,u,unew,w,wnew,theta,dzf,dzz,dzznew,xc,ctop,ctopnew,nk,Nk)

		# plot results
		if plot==1 and n==ntime-1:
			PlotResults(hnew,unew,wnew,snew,zb,zt,xc,xe,d,ctopnew,Nk)
		# update u w s dzz ctop
		u=unew
		w=wnew
		s=snew
		h=hnew
		dzz=dzznew
		ctop=ctopnew

		if n%ntout==0:
			if n==0:
				nn=n+1
			else:
				nn=n
			elapsed=time.time()-t
			print nn,'of',ntime,'finished,there are',elapsed*(ntime-n),'seconds left.'

	return h,u,w,s,cdb

# initialize physical variables
def InitialPhysicalVariables(xc,xe,d,dt,nk):
	nc=xc.shape[0]
	h=np.zeros((nc,1))
	h[0]=bc.BoundaryFreeSurface(0,dt)
	z0b=np.zeros((nc,1))
	s=np.zeros((nc,nk))
	u=np.zeros((nc+1,nk))
	w=np.zeros((nc,(nk+1)))
	for i in range(1,nc):
		h[i]=ic.ReturnFreeSurface(xc[i])
		z0b[i]=ic.ReturnBottomRoughness(xc[i])
		for j in range(nk):
			s[i][j]=ic.ReturnScalar(xc[i])
			w[i][j]=ic.ReturnVerticalVelocity(xc[i])
	for j in range(nk):
		s[0][j]=bc.DownstreamBoundaryScalar()
		for j in range(nk):
			if i<nc:
				u[i][j]=ic.ReturnHorizontalVelocity(xe[i])
			else:
				u[i][j]=bc.BoundaryVelocity(xe[i])
	return h,u,w,s,z0b

# calculate the index of top cell and dzz for each layer
def UpdateDzz(xc,h,d,zb,zt,Nk,dz,nk):
	ctop=np.zeros(np.shape(xc))
	dzz=np.zeros((len(xc),nk))
	for i in range(len(xc)):
		switchk=0
		for j in range(Nk[i]):
			dzz[i][j]=0
			if h[i]>=zb[0][j] and switchk==0:
				ctop[i]=j
				switchk=1
                                dzz[i][j]=-zb[0][j]+h[i];
			elif h[i]>=zb[0][j] and switchk==1:
				dzz[i][j]=dz

			if j==Nk[i]-1 and ctop[i]<(Nk[i]-1):
				dzz[i][j]=zt[0][j]+d[i]
		        elif j==Nk[i]-1 and ctop[i]==(Nk[i]-1):
				dzz[i][j]=h[i]+d[i]
			if Nk[i]==1:
				dzz[i][j]=h[i]+d[i]
	return dzz,ctop

# calculate flux height 
def SetFluxHeight(dzz,ctop,Nk,u,h,nk):
	dzf=np.zeros((len(h)+1,nk))
	D=np.zeros((len(h)+1,1));
	etop=np.zeros((len(h)+1,1))
	Nke=np.zeros((len(h)+1,1))
	# set etop
	
	for i in range(1,len(h)+1):
		nc1=i-1
		nc2=i
		if i==len(h):
			nc2=i-1
		etop[i]=max(ctop[nc1],ctop[nc2])
		Nke[i]=min(Nk[nc1],Nk[nc2])

	# set dzf
	for i in range(1,len(h)+1):
		D[i]=0
		for j in range(etop[i],Nke[i]):
			nc1=i-1
			nc2=i
			if i==len(h):
				nc2=i-1
			if u[i][j]==0:
				dzf[i][j]=min(dzz[nc1][j],dzz[nc2][j])
			if u[i][j]>0:
				dzf[i][j]=dzz[nc1][j]
			if u[i][j]<0:
				dzf[i][j]=dzz[nc2][j]
			D[i]=D[i]+dzf[i][j]
	
	return dzf,D,etop,Nke

# set drag coefficient 
def SetDragCoefficient(dzf,Nke,z0b,u):
	cdb=np.zeros((len(z0b)+1,1))

	for i in range(1,len(z0b)):
		cdb[i]=pow(np.log(0.5*dzf[i][Nke[i][0]-1]/(z0b[i-1]+z0b[i])*2)/0.41,-2);
	return cdb

# calculate velocity and free surface
def CalculateUandH(n,dt,xc,xe,dzf,D,dzz,Nk,Nke,etop,ctop,nk,u,h,cdb,theta):
	# variables
	nu=1e-4 # momentum diffusivity
	nc=len(xc)
        dx=xc[1]-xc[0] # distance in two cell center
	g=9.81 # gravity
	#theta=0.55 # theta method
	utmp=np.zeros(np.shape(u))
	# section for velocity
	for i in range(1,nc):
		# for each computational edge
	        A=np.zeros((nk,nk))
		b=np.zeros((nk,1))
		# beyond etop all zero
		if Nke[i]-etop[i]>1:
			for j in range(0,etop[i]):
				A[j][j]=1
				b[j]=0
			# below nke all zero
			for j in range(Nke[i],nk):
				A[j][j]=1
				b[j]=0
                	# inner part
			if Nke[i]-etop[i]>2:
				for j in range(etop[i]+1,Nke[i]-1):
					A[j][j]=1+nu*dt*(2/(dzf[i][j-1]+dzf[i][j])+2/(dzf[i][j]+dzf[i][j+1]))/dzf[i][j]
					A[j][j-1]=-nu*dt*2/(dzf[i][j-1]+dzf[i][j])/dzf[i][j]
					A[j][j+1]=-nu*dt*2/(dzf[i][j]+dzf[i][j+1])/dzf[i][j]
					b[j]=u[i][j]-g*dt*(h[i]-h[i-1])/dx*(1-theta)
			# etop
			j=etop[i][0]
			A[j][j]=1+nu*dt*2/(dzf[i][j]+dzf[i][j+1])/dzf[i][j]
			A[j][j+1]=-nu*dt*2/(dzf[i][j]+dzf[i][j+1])/dzf[i][j]
			b[j]=u[i][j]-g*dt*(h[i]-h[i-1])/dx*(1-theta)
			# Nke
			j=Nke[i][0]-1
			A[j][j]=1+dt*(nu*2/(dzf[i][j-1]+dzf[i][j])+cdb[i]*abs(u[i][j]))/dzf[i][j]
			A[j][j-1]=-nu*dt*2/(dzf[i][j-1]+dzf[i][j])/dzf[i][j]
			b[j]=u[i][j]-g*dt*(h[i]-h[i-1])/dx*(1-theta)
		else:
			for k in range(0,Nke[i]-1):
				A[k][k]=1.0
				b[k]=0
			for k in range(Nke[i],nk):
				A[k][k]=1.0
				b[k]=0
			j=Nke[i][0]-1
			A[j][j]=1+dt*cdb[i]*abs(u[i][j])/dzf[i][j]
			b[j]=u[i][j]-g*dt*(h[i]-h[i-1])/dx*(1-theta)
	#	print A,nk 
	#	print b
		# solve A/b to get utmp
		#print A
		#print 'nke',Nke[i][0],'etop',etop[i][0] 
	        uu=np.linalg.solve(A,b)
		utmp[i]=uu.T
	for j in range(etop[nc],Nke[nc]):
		utmp[nc][j]=bc.BoundaryVelocity(xe[nc])
	# calculate free surface
        A=np.zeros((nc,nc))
	b=np.zeros((nc,1))
	# tidal boundary condition
	A[0][0]=1
	b[0]=bc.BoundaryFreeSurface(n,dt)
	for i in range(1,nc-1):
		A[i][i]=dx+theta*theta*dt*dt*g/dx*(D[i]+D[i+1])
		A[i][i-1]=-theta*theta*dt*dt*g/dx*D[i]
		A[i][i+1]=-theta*theta*dt*dt*g/dx*D[i+1]
		b[i]=dx*h[i]
		for j in range(etop[i],Nke[i]):
			b[i]=b[i]+dt*(theta*utmp[i][j]+(1-theta)*u[i][j])*dzf[i][j]
		for j in range(etop[i+1],Nke[i+1]):
			b[i]=b[i]-dt*(theta*utmp[i+1][j]+(1-theta)*u[i+1][j])*dzf[i+1][j]

	A[nc-1][nc-1]=dx+theta*theta*dt*dt*g/dx*D[nc-1]
	A[nc-1][nc-2]=-theta*theta*dt*dt*g/dx*D[nc-1]
	b[nc-1]=dx*h[nc-1]
	for j in range(etop[nc-1],Nke[nc-1]):
		b[nc-1]=b[nc-1]+dt*((1-theta)*u[nc-1][j]+theta*utmp[nc-1][j])*dzf[nc-1][j]
	for j in range(etop[nc],Nke[nc]):
		b[nc-1]=b[nc-1]+dt*((1-theta)*u[nc][j]+theta*utmp[nc][j])*dzf[nc][j]

	hh=np.linalg.solve(A,b)
	hnew=hh
        # Set up boundaries
	#h[0]=bc.BoundaryFreeSurface(n,dt)
	# finish velocity
	for i in range(1,nc):
		for j in range(etop[i],Nke[i]):
			utmp[i][j]=utmp[i][j]-theta*g*dt/dx*(hnew[i]-hnew[i-1])

	#print '32131',hnew[1]-h[1]-(D[1]*utmp[1][0]*theta+D[1]*u[1][0]*(1-theta)-D[2]*utmp[2][0]*theta-D[2]*u[2][0]*(1-theta))

	return utmp,hnew

# Calculate Vertical Velocity
def CalculateVerticalVelocity(dt,xc,dzf,Nk,ctopnew,unew,nk):
        nc=len(xc)
	dx=xc[1]-xc[0]
	wnew=np.zeros((nc,nk+1))
	for i in range(1,nc):
		for j in range(Nk[i]-1,ctopnew[i]-1,-1):
			wnew[i][j]=wnew[i][j+1]+unew[i][j]*dzf[i][j]/dx-unew[i+1][j]*dzf[i+1][j]/dx
	# set vertical velocity for boundary cell
	for j in range(ctopnew[0],Nk[i]):
		wnew[0][j]=wnew[1][j]
	return wnew

def PlotResults(h,unew,wnew,snew,zb,zt,xc,xe,d,ctopnew,Nk):
	nk=len(zb.T)
	nc=len(xc)
	ze=np.zeros((1,nk+1))
	for i in range(nk):
		ze[0][i]=zt[0][i]

	ze[0][nk]=zb[0][nk-1]
        xbound=np.zeros((nk+1,nc+1))
	zbound=np.zeros((nc+1,nk+1))
	for i in range(nc+1):
	        #print np.shape(zbound),np.shape(ze)
		zbound[i]=ze
	for i in range(nk+1):
		#print np.shape(xbound),np.shape(xe)
		xbound[i]=xe.T
	zbound=zbound.T
        uc=np.zeros((nc,nk))
	wc=np.zeros((nk,nc))
	sc=snew
	wtmp=wnew.T
	for i in range(nc):
		for j in range(0,ctopnew[i]):
			uc[i][j]=0#np.nan
			sc[i][j]=0#np.nan	
		for j in range(ctopnew[i],Nk[i]):
			uc[i][j]=(unew[i][j]+unew[i+1][j])/2
		for j in range(Nk[i],nk):
			uc[i][j]=0#np.nan
			sc[i][j]=0#np.nan
	for j in range(0,ctopnew[0]):
		uc[0][j]=0#np.nan
	for j in range(ctopnew[0],Nk[0]):
		uc[0][j]=unew[1][j]
	for j in range(Nk[0],nk):
		uc[0][j]=0#np.nan

	sc=sc.T
	uc=uc.T
	for i in range(nk):
		for j in range(nc):
			if i<ctopnew[j]:
				wc[i][j]=0#np.nan
			elif i>=ctopnew[j] and i<Nk[j]:
				wc[i][j]=(wtmp[i][j]+wtmp[i+1][j])/2
			else:
				wc[i][j]=0#np.nan
	
	#print 'u',uc
	#print 'w',wc
	# unew 
	plt.subplot(4,1,1)
        z_min,z_max=uc.min(),uc.max()
	plt.pcolor(xbound,zbound,uc,vmin=z_min,vmax=z_max)
        plt.tick_params(\
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')
	plt.title('horizontal velocity (m/s)')
	plt.colorbar()
	#plt.show()
	 
	# wnew
	plt.subplot(4,1,2)
	z_min,z_max=wc.min(),wc.max()
	plt.pcolor(xbound,zbound,wc,vmin=z_min,vmax=z_max)
        plt.tick_params(\
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')
	plt.title('vertical velocity (m/s)')
	plt.colorbar()
	#plt.show()

 	# h
	plt.subplot(4,1,4)
	plt.plot(xc,h)
	plt.plot(xc,-d)
	plt.title('free surface (m)')
	plt.xlabel('x (m)')
	plt.colorbar()
	#plt.show()

	# s
	plt.subplot(4,1,3)
	z_min,z_max=sc.min(),sc.max()
	plt.pcolor(xbound,zbound,sc,vmin=z_min,vmax=z_max)
        plt.tick_params(\
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')
	plt.title('scalar concentration')
	#plt.axis([xbound.min(),xbound.max(),zbound.min(),zbound.max()])
	plt.colorbar() 	
	plt.show()

