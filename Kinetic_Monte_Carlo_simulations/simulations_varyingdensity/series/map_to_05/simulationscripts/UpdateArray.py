from __future__ import division, unicode_literals, print_function
import numpy as np

#This function is called by rollingbead_f.py, and updates the positions of the linkers on the channel surface, in case the particle has moved during the previous step.

def UpdateArray(linkers,F,par):

	if par['attachmentdirection'] == 0:																#In case the attachment direction is normal to bottom surface
		Rmax = np.sqrt(np.square(par['a'])-np.square(par['h']-par['Lmax']))

	else:																							#In case the attachment direction is normal to particle surface. this setting was used in all simulations used in the paper.
		amax = np.arccos(np.divide(par['h'],par['a']+par['Lmax']))
		Rmax = np.multiply(par['h'],np.tan(amax))

	dr = 1/np.sqrt(par['rho'])																		#typical distance between linker positions on the channel surface. Depends on the density rho on the surface

	#First, we remove linker positions on the left of the particle in case the particle has moved forward during the previous step. To be removed, bonds have to satisfy a few conditions:
	if np.min(linkers.p0[:,0])<(-1)*(Rmax+3*dr):
		ii1 = np.where(linkers.p0[:,0]<0)[0]														#1. Find nodes with an x position smaller than 0
		ii2 = np.where((np.square(linkers.p0[:,0])+np.square(linkers.p0[:,1]))>np.square(Rmax))[0]	#2. Find nodes at least Rmax removed from the center of the bead
		ii = np.intersect1d(ii1,ii2)																#3. Find the intersect of conditions ii1 and ii2, these are nodes to remove
		ikeep = np.setdiff1d(np.arange(np.shape(linkers.p0)[0]),ii)									#Any of the remaining nodes are nodes to keep

		linkers.p0 = linkers.p0[ikeep,:]
		linkers.p1 = linkers.p1[ikeep,:]

		linkers.connected = linkers.connected[ikeep]
		F = F[ikeep]

	#Add extra linkers to the right of the particle when the particle has moved:
	if np.max(linkers.edge[:,0])<(Rmax+dr): 														#The condition to add more linkers is when the edge has moved to within Rmax+dr distance from the particle center of mass.
		ddr = par['linkerpositionnoise']*dr

		if par['extralinkers']>0:
			nadd = par['extralinkers']

			for i in range(len(linkers.edge[:,0])): 												#add nadd unbound linkers in each row, using a loop, making sure to displace each one of the randomly according to linkerpositionnoise.
				xadd = linkers.edge[i,0] + dr* np.array(range(1,nadd+1))[np.newaxis].T+ddr* (1-2*np.random.rand(nadd,1))
				yadd = linkers.edge[i,1] *np.ones((nadd,1)) +ddr* (1-2*np.random.rand(nadd,1))

				newlinkersp0s = np.hstack((xadd,yadd))												
				linkers.p0 = np.vstack((linkers.p0,newlinkersp0s))
				linkers.edge[i,0]=xadd[-1]

				newconnections = np.zeros((nadd,1),dtype=bool)
				linkers.connected = np.append(linkers.connected,newconnections)

				newp1array = np.zeros((nadd,3),dtype=float)
				linkers.p1 = np.vstack((linkers.p1,newp1array))

				newFarray = np.zeros((nadd,1),dtype=float)
				F = np.vstack((F,newFarray))
				
	#Return the updated array of linkers.
	return(linkers,F)