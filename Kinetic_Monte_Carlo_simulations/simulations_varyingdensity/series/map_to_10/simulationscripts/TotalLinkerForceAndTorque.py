from __future__ import division, unicode_literals, print_function
import numpy as np

def TotalLinkerForceAndTorque(linkers,par):

	#This function returns the force Fx=[Fx] and total force F=|F| and Torque around y, all normalized by kT
	#also returns the relative stretch of each linker x

	#calculate total force and torque on linkers

	#First, we calculate the extension L on each linker:
	L=np.zeros(np.shape(linkers.p0[:,0]))[np.newaxis].T
	r = np.zeros((np.shape(linkers.p0[linkers.connected,0])[0],3))
	r[:,0:2] = linkers.p0[linkers.connected,:] - linkers.p1[linkers.connected,0:2]
	r[:,2] = (-1)*linkers.p1[linkers.connected,2]
	L[linkers.connected,0] = np.sqrt(np.sum(np.square(r),axis=1))							#Calculate the stretch of each linker

	x = L/par['Lmax']																		#Calculate the relative stretch of each linker (compared to contour length Lmax)

	if par['polymermodel']==1:																#If the freely jointed chain force extension model is used
		F = np.divide(np.multiply(3,L), par['Lmax']*par['lK'])

	else:																					#Else, if the worm-like chain force extension model is used (This setting was used during all the simulations in the article)

		F = (2/par['lK']) * (np.divide(0.25,np.square(1-x))-0.25+x-0.8*np.power(x,2.15))
	
																							#If the bonds are stretched beyond 99% of their contour length, we set the extension force to a linear interpolation between 0.99 extension and 0.989 extension to avoid computational instabilities
																							
		F1 = (np.divide(0.25,np.square(1-0.99))-0.25+0.99-0.8*np.power(0.99,2.15))
		F2 = (np.divide(0.25,np.square(1-0.989))-0.25+0.989-0.8*np.power(0.989,2.15))

		#linear:
		a = (F1-F2)/(0.001)
		b = F1 - np.multiply(a,0.99)
		F[x>0.99] = (2/par['lK'])*(np.multiply(a,x[x>0.99])+b) 

	#Compute the x components and z components of the Forces (this is done to compute the torque.)
	Fx = np.zeros((np.shape(F)[0],1))
	Fz = np.zeros((np.shape(F)[0],1))
	Fx[linkers.connected,0] = np.divide(np.multiply(F[linkers.connected,0],r[:,0]),L[linkers.connected,0])
	Fz[linkers.connected,0] = np.divide(np.multiply(F[linkers.connected,0],r[:,2]),L[linkers.connected,0])

	#Finally, compute the torque every bond applies to the particle.
	#Torque: Ty = (z-h) * fx-x*fz
	T = np.multiply(linkers.p1[:,2][np.newaxis].T-par['h'],Fx) - np.multiply(linkers.p1[:,0][np.newaxis].T,Fz)

	#Return the array of Forces, Torques, relative extensions, and x components of the force of all the linkers..
	return(F,T,x,Fx) #Return 