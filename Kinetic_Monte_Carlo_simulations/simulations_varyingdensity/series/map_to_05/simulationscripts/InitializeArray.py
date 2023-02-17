from __future__ import division, unicode_literals, print_function
import numpy as np

#This function initializes the positions of linkers on the channel surface, at the start of the simulation. It is called by rollingbead_f.py 

def InitializeArray(par):

	#Initialize a class with linkers and create sub-parameters:
	#Every linker is stored as this class.
	#
	#p0:		x-y coordinates on bottom surface
	#p1:		x-y-z coordinates on bead surface (for attached linkers)
	#connected: boolean array listing attached linkers
	#edge: x-y coordinates of edge of array (at the front)


	class linker:
		def __init__(self,p0,p1,connected,edge):
			self.p0=p0
			self.p1=p1
			self.connected=connected
			self.edge=edge

	if par['attachmentdirection']==0: 					#If the linkers are attached normal to bottom surface
		Rmax = np.sqrt(np.square(par['a'])-np.square(par['h']-par['Lmax']))
	else:												#If the linkers are attached normal to particle surface (this was chosen in all simulations for the paper.)
		amax = np.arccos(np.divide(par['h'],(par['a']+par['Lmax'])))
		Rmax = par['h']*np.tan(amax)					#Compute Rmax, the maximum distance linkerpositions will be placed compared to the center of the particle

	dr = 1/np.sqrt(par['rho']) 							#Compute the average distance between linkers based on the imposed linker density rho. Rmax. This is for computational efficiency, because linkers placed further than the bond contour length Lmax would be very unlikely to form.

	#Create an array of linker coordinates (origin at bead center)
	rj = np.arange(start=-Rmax-dr,stop=Rmax,step=dr)

	p=np.array([])
	edge=np.array([])

	for j in range(len(rj)):

		xj = rj
		yj = rj[j]*np.ones(len(xj))

		pj = np.transpose(np.vstack([xj,yj]))
		pj = pj[np.square(pj[:,0])+np.square(pj[:,1])<np.square(Rmax),:]	#Evaluate if the linkers to be placed are not further than Rmax from the particle center.

		if pj.any():

			if not p.any():
				p = pj[:]
			else:
				p=np.vstack([p,pj])

		if np.shape(pj)[0]>0:

			#Also compute the position of the edges, the linker positions furthest away from the particle. This edge is used throughout the simulation to determine if new linker positions have to be generated.
			if not edge.any():
				edge = pj[-1,:]
			else:
				edge = np.vstack([edge,pj[-1,:]])

	#Add extra tethers on the right (positive x direction) so that linkers can still be formed if the particle moves somewhat.
	if par['extralinkers']>0: #create extra
		nadd=par['extralinkers']
		for i in range(len(edge)): #add nadd tethers in each row

			xadd = (edge[i,0]+ dr*np.transpose(range(1,nadd+1))[np.newaxis]).T
			yadd = edge[i,1]*np.ones((nadd,1))
			newps=np.hstack([xadd,yadd])

			p = np.vstack([p,newps])	#Add the new linkers to the array p
			edge[i,0]=xadd[-1]			#Update the edge positions now that the new linkers have been added.

	#add random noise to linker positions, using the random number generator. 
	ddr = par['linkerpositionnoise']*dr
	p = p + ddr*(1-2*np.random.rand(np.shape(p)[0],np.shape(p)[1]))

	#Now that all the linker positions have been placed on the channel surface, the next step is to randomly allow some of these linkers to form, using exp(-Delta G) as a weighing factor.
	#This ensures that linkers that have to stretch further to form bonds are less likely to do so.

	#To this end, we first compute the bond length L for every potential bond that could be formed:
	if par['attachmentdirection']==0: 					#If the attachment direction is normal to the channel surface bottom surface

		L = par['h'] - np.sqrt(np.square(par['a'])-(np.square(p[:,0])+np.square(p[:,1])))

	else:												#If the attachment direction is normal to the particle surface (used in all simulations in the paper.)
		r = np.sqrt(np.square(p[:,0]) + np.square(p[:,1]))
		r[r==0] = 1e-10									#if r =0, set to a very low value to avoid instabilities due to divisions by 0
		alpha = np.arctan(r/par['h'])
		L = np.divide(r,np.sin(alpha))-par['a']

	x = np.divide(L,par['Lmax']) #Relative extension of the polymers relative to contour length
	#Randomly bind linkers with probability K/(1+K) with K=(kon/koff)^0*exp(-DeltaG_el/kT)

	if par['polymermodel']==1: #if we assume an ideal chain force-extension model
		dG = np.divide(3*np.square(L),par['Lmax']*par['lK'])
		
	if par['polymermodel']==2: #if we assume a worm-like chain force-extension model (this was chosen in all simulations in this article.)
		dG = (2/par['lK'])*(np.divide(0.25*np.square(x),1-x)+0.5*np.square(x)-(0.8/3.15)*np.power(x,3.15))


	dG[x>0.99]=np.inf #remove overstretched chains
	
	#Now that we know the extensional free energy of every potential bond, we calculate the on-rate for every bond as follows:
	K0 = par['kon0']/par['koff0']
	K = K0*np.exp(-dG)						

	#And we choose which bonds to connect!
	connected = np.random.rand(np.shape(p)[0])< np.divide(K,1+K)		#connect with the probability K/(1+K) connected, forms a list of booleans that specify which of the linkers is connected

	#next, actually connect the tethers by setting the positions p1 on the bead that they connect to!
	p1 = np.zeros((np.shape(p)[0],3))

	if par['attachmentdirection']==0:	#if attachment direction is set normal to the surface
		L=L[np.newaxis].T
		p1[connected,:] = np.hstack([p[connected,:],L[connected,:]])

	else:								#If attachment direction is set normal to the particle surface
		
		r1 = par['a']*np.sin(alpha)

		p1[connected,0] = np.divide(np.multiply(p[connected,0],r1[connected]),r[connected])
		p1[connected,1] = np.divide(np.multiply(p[connected,1],r1[connected]),r[connected])

		z1 = par['h']-par['a']*np.cos(alpha)

		p1[connected,2]=z1[connected]

	#Finally: Now that we have chosen which linkers are connected, we store all the bound and unbound linkers in the linker class, and return this to the rollingbead_f function!
	linkers = linker(p0=p,p1=p1,connected=connected,edge=edge)
	return linkers
