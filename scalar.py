# scalar.py
# calculate scalar transport based on 
# the old and new velocity field
# Yun Zhang
# Stanford University
# ver 1.0 10/10/2014

import numpy as np
import boundary as bc

def ScalarTransport(dt,hnew,h,s,u,unew,w,wnew,theta,dzf,dzz,dzznew,xc,ctop,ctopnew,nk,Nk):
	snew=np.zeros(np.shape(s))
	nc=len(xc)
	dx=xc[1][0]-xc[0][0]
	nu_s=1e-4
	# calculate scalar field for each cell
	for i in range(1,nc):
		if ctop[i]<ctopnew[i]:
			top=ctopnew[i]
		else:
			top=ctop[i]
		A=np.zeros((nk,nk))
		b=np.zeros((nk,1))
		if Nk[i]-top>1:
			for j in range(0,top):
				b[j]=0
				A[j][j]=1.0
			for j in range(Nk[i],nk):
				b[j]=0
				A[j][j]=1.0
			# interior point
			for j in range(top,Nk[i]):
				if u[i][j]>0:
					sleft=s[i-1][j]
				else:
					sleft=s[i][j]
				if i!=nc-1:
					if u[i+1][j]>0:
						sright=s[i][j]
					else:
						sright=s[i+1][j]			
				else:
					if u[i+1][j]>0:
						sright=s[i][j]
					else:
						sright=bc.UpstreamBoundaryScalar()
	
				if j>top and j<Nk[i]-1:
					if w[i][j]>0:
						stop=s[i][j-1]
					else:
						stop=s[i][j]
					if w[i][j+1]>0:
						sbot=s[i][j+1]
					else:
						sbot=s[i][j]
				elif j==top:
					if ctop[i]<ctopnew[i]:
						if w[i][j]>0:
							stop=s[i][j-1]
						else:
							stop=s[i][j]
					else:
						stop=0

					if w[i][j+1]>0:
						sbot=s[i][j+1]
					else:
						sbot=s[i][j]
				elif j==Nk[i]-1:
					if w[i][j]>0:
						stop=s[i][j-1]
					else:
						stop=s[i][j]
					sbot=0

				# old s
				b[j]=dzz[i][j]*s[i][j]
				# horizontal advection
				b[j]=b[j]+dt*dzf[i][j]*((1-theta)*u[i][j]+theta*unew[i][j])*sleft/dx-dt*dzf[i+1][j]*((1-theta)*u[i+1][j]+theta*u[i+1][j])*sright/dx
				# vertical advection
				b[j]=b[j]-dt*(1-theta)*(w[i][j]*stop-w[i][j+1]*sbot)
				
				# vertical diffusion
				if j>top and j<Nk[i]-1:
					b[j]=b[j]+(1-theta)*nu_s*dt*(2*(s[i][j-1]-s[i][j])/(dzz[i][j]+dzz[i][j-1])-2*(s[i][j]-s[i][j+1])/(dzz[i][j]+dzz[i][j+1]))
				elif j==top:
					if ctop[i]<ctopnew[i]:
						b[j]=b[j]+(1-theta)*nu_s*dt*(2*(s[i][j-1]-s[i][j])/(dzz[i][j]+dzz[i][j-1])-2*(s[i][j]-s[i][j+1])/(dzz[i][j]+dzz[i][j+1]))
					else:
						b[j]=b[j]+(1-theta)*nu_s*dt*(-2*(s[i][j]-s[i][j+1])/(dzz[i][j]+dzz[i][j+1]))
				elif j==Nk[i]-1:
					b[j]=b[j]+(1-theta)*nu_s*dt*(2*(s[i][j-1]-s[i][j])/(dzz[i][j]+dzz[i][j-1]))
				if j>top and j<Nk[i]-1:
					A[j][j]=dzznew[i][j]+dt*nu_s*2*(1/(dzznew[i][j-1]+dzznew[i][j])+1/(dzznew[i][j+1]+dzznew[i][j]))
					A[j][j-1]=-dt*nu_s*2/(dzznew[i][j-1]+dzznew[i][j])
					A[j][j+1]=-dt*nu_s*2/(dzznew[i][j]+dzznew[i][j+1])
		   			if wnew[i][j]>0:
						A[j][j]=A[j][j]+dt*theta*wnew[i][j]
					else:
						A[j][j-1]=A[j][j-1]+dt*theta*wnew[i][j]
					if wnew[i][j+1]>0:
						A[j][j+1]=A[j][j+1]-dt*theta*wnew[i][j+1]
					else:
						A[j][j]=A[j][j]-dt*theta*wnew[i][j+1]
				elif j==top:
					A[j][j]=dzznew[i][j]+dt*nu_s*2*(1/(dzznew[i][j+1]+dzznew[i][j]))
					A[j][j+1]=-dt*nu_s*2/(dzznew[i][j]+dzznew[i][j+1])
					if wnew[i][j+1]>0:
						A[j][j+1]=A[j][j+1]-dt*theta*wnew[i][j+1]
					else:
						A[j][j]=A[j][j]-dt*theta*wnew[i][j+1]
				elif j==Nk[i]-1:
					A[j][j]=dzznew[i][j]+dt*nu_s*2*(1/(dzznew[i][j-1]+dzznew[i][j]))
					A[j][j-1]=-dt*nu_s*2/(dzznew[i][j]+dzznew[i][j-1])
					if wnew[i][j]>0:
						A[j][j]=A[j][j]+dt*theta*wnew[i][j]
					else:
						A[j][j-1]=A[j][j-1]+dt*theta*wnew[i][j]
		else:

			for j in range(0,Nk[i]-1):
				b[j]=0
				A[j][j]=1.0
			for j in range(Nk[i],nk):
				b[j]=0
				A[j][j]=1.0 
			# for only one layer
			j=Nk[i][0]-1
  			if u[i][j]>0:
				sleft=s[i-1][j]
			else:
				sleft=s[i][j]
			if i!=nc-1:
				if u[i+1][j]>0:
					sright=s[i][j]
				else:
					sright=s[i+1][j]			
			else:
				if u[i+1][j]>0:
					sright=s[i][j]
				else:
					sright=bc.UpstreamBoundaryScalar()
			stop=0
			sbot=0
			if ctop[i]<ctopnew[i]:	
				if w[i][j]>0:
					stop=s[i][j-1]
				else:
					stop=s[i][j]

			# old s
			b[j]=dzz[i][j]*s[i][j]
			# horizontal advection
			b[j]=b[j]+dt*dzf[i][j]*((1-theta)*u[i][j]+theta*unew[i][j])*sleft/dx-dt*dzf[i+1][j]*((1-theta)*u[i+1][j]+theta*unew[i+1][j])*sright/dx
			# vertical advection
			b[j]=b[j]-dt*(1-theta)*(w[i][j]*stop-w[i][j+1]*sbot)
			if ctop[i]<ctopnew[i]:
					b[j]=b[j]+(1-theta)*nu_s*dt*2*(s[i][j-1]-s[i][j])/(dzz[i][j]+dzz[i][j-1])				
			A[j][j]=dzznew[i][j]
		#if i==1:
		#	print 'dzf is ', dzznew[i]
		#	print 'A is ', A
		#	print 'b is ', b
		#	print 'top',stop,'bot',sbot,nu_s
		#	print 'u',u[1][0],u[2][0],'unew',unew[1][0],unew[2][0]
		#	print h[1]+(u[1][0]*(1-theta)+unew[1][0]*theta)*dzf[1][0]-(u[2][0]*(1-theta)+theta*unew[2][0])*dzf[2][0]-hnew[1]
		ss=np.linalg.solve(A,b)
		snew[i]=ss.T
		if ctopnew[i]<ctop[i]:
			for j in range(ctopnew[i],ctop[i]):
				snew[i][j]=snew[i][ctop[i][0]]
	# set boundary for scalar
	for j in range(ctopnew[0],Nk[0]):
		snew[0][j]=bc.DownstreamBoundaryScalar()
	return snew







