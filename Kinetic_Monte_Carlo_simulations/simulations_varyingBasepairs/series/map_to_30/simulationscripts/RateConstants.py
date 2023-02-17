from __future__ import division, unicode_literals, print_function
import numpy as np
#This script is called by rollingbead_f.py, and it computes the on-rates kon for all unbound linkers and the off-rates koff for all bound linkers.
#These rate constants are then used by the KMC algorithm to decide which step to take.

def RateConstants(linkers,par,F):
	k = np.zeros((np.shape(linkers.p0)[0],1))

	#off-rates
	k[linkers.connected]=par['koff0']*np.exp(par['del']*F[linkers.connected]) 					#This is the Bell-equation that allows one to compute the off rates as a function of the extension force

	#on-rates
	if par['attachmentdirection']==0:															#In case the attachment direction is normal to the surface
		R=par['h']-np.sqrt(np.square(par['a']) - (np.square(linkers.p0[linkers.connected!=True,0]) + linkers.p0[linkers.connected!=True,1]))
	else:																						#In case the attachment direction is normal to the particle surface.
		r=np.sqrt(np.square(linkers.p0[linkers.connected!=True,0]) + np.square(linkers.p0[linkers.connected!=True,1])) #Compute R, the length a bond would have if it was drawn from the channel surface to the particle surface.
		r[r==0] = 1e-10 																		#If r is practically zero, set r to a very low value, to avoid divisions by zero
		alpha=np.arctan(r/par['h'])
		R=np.divide(r,np.sin(alpha))-par['a']

	if par['polymermodel']==1: 																	#If the freely-jointed chain behavior is used.
		dG=3*np.square(R)/(par['Lmax']*par['lK'])												#Extensional change in free energy delta G
	else: 																						#If the worm-like chain model is used.
		x=R/par['Lmax']																			#Compute the relative extension every non-formed bond would have if it was formed

		Fi = (2/par['lK'])*(np.divide(0.25,np.square(1-x))-0.25+x-0.8*np.power(x,2.15))			#Compute the extension force using the WLC force extension model
		dG = (2/par['lK'])*(np.divide(0.25*np.square(x),(1-x))+0.5*np.square(x)-(0.8/3.15)*np.power(x,3.15)) #Compute the extension energy delta G

		A=Fi*par['del']-dG
		A[A>0]=0
		A[x>1]=-np.inf

	kon = par['kon0']*np.exp(A)																	#Compute the on-rate kon for every unbound linker
	kon.shape = (np.shape(kon)[0],1)

	k[linkers.connected!=True]=kon

	#return the array with all rate constants, for both the formation and dissociation of the bonds.
	return(k)
