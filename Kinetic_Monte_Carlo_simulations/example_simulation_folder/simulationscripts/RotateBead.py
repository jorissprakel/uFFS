from __future__ import division, unicode_literals, print_function
import numpy as np

#This script is called by Equilibrate.py when the particle is rotated during the gradient descent approach.
#The script recalculates the positions r1 of every bound linker on the particle as the particle rotates.

def RotateBead(r1,par,theta):

	#Rotate the bead over an angle theta (clockwise).

	c = np.cos(theta)
	s = np.sin(theta)

	r = np.copy(r1)

	r[:,0] = r1[:,0]*c - (par['h']-r1[:,2])*s
	r[:,2] = par['h'] - (r1[:,0]*s + (par['h']-r1[:,2])*c)

	return(r)