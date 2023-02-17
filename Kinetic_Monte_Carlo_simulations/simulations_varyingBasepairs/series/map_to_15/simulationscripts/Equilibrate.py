from __future__ import division, unicode_literals, print_function
import numpy as np
import TotalLinkerForceAndTorque
import RotateBead
import copy

#This script equilibrates the energy by minimizing the forces and torques on the particle due to the shear force and torque and due to the forces and torques imposed by the linkers.
#This is done using a gradient descent algorithm, in which the bead is rotated, until the energy is minimized

def Equilibrate(linkers,par,FX,TY):

	#Equilibrate the forces and torques by demanding that sum(fx)+Fx=0 and sum(T)+Ty=0
	#or only sum(T)+Ty=0 in case there is no slip
	#with Fx and Ty the applied force and torque
	#update linkers.p1 (bead rotation) and linkers.p0 (surface movement)
	#also return movement and angle

	class linker:
		def __init__(self,p0,p1,connected,edge):
			self.p0=p0
			self.p1=p1
			self.connected=connected
			self.edge=edge

	
	#Parameters used in the gradient descent algorithm
	alphafactor = 0.8	#reduction factor of steplength
	TOL = 1e-6			#Tolerance. Stop the gradient descent if the maximum torque is lower than this value.
	SMALL = 1e-7		#Size of the small bead rotation and translation steps taken during the optimisation
	MAXIT = 10000		#maximum number of iterations. System returns an error if MAXIT is reached, but the energy is not yet minimized

	#First, check the magnitude of the force compared to the force at maximum stretch of a linker
	if par['polymermodel']==1: 		#Freely Jointed Chain
		Fmax = 3/par['lK']
	else:							#Worm-like Chain (used in the simulations in the article)
		Fmax = (2/par['lK']) * (np.divide(0.25,np.square(1-0.9999))-0.25+0.9999-0.8*np.power(0.9999,2.15))
	if FX>Fmax:
		warn='shear force exceeds maximum linker stretch force'

	#Next, we calculate the force, torque, relative extension, and x-component of the force of every bond
	F,T,x,Fx = TotalLinkerForceAndTorque.TotalLinkerForceAndTorque(linkers,par)


	#Next, we equilibrate the force and torque. In the simulations used in the paper we set 'slipconditions' to 1, in which we both equilibrate the force and torque. It is also possible to set the slipconditions to 0, in which we only equilibrate the torque.
	if par['slipconditions'] == 0:		#no slip, only equilibrate the torque
		#initial value of u = [dQ]
		u = [0]
		f = np.sum(T[linkers.connected])+TY
		it = 0 								#iteration counter
		outcon = True 						#out condition is true

		#Start iterating, until the maximum torque is lower than the tolerance
		while np.linalg.norm(f) > TOL * np.max(TY,0):
			
			it = it+1
			if it > MAXIT: #Stop and provide an error message if the maximum number of iterations has been reached
				print("Maximum iterations reached")
				outcon = False
				break

			#Get Jacobian
			linkers1 = copy.deepcopy(linkers)
			linkers1.p1[linkers.connected,:] = RotateBead.RotateBead(linkers1.p1[linkers.connected,:],par,SMALL)

			linkers1.p0[:,0] = linkers1.p0[:,0] - par['a']*SMALL
			F1,T1,x1,Fx1 = TotalLinkerForceAndTorque.TotalLinkerForceAndTorque(linkers1,par)

			dTdQ = np.divide((np.sum(T1[linkers.connected]) - np.sum(T[linkers.connected])),SMALL)
			J=np.copy(dTdQ)

			if np.abs(dTdQ)<1e-10:
				print("Singular Matrix\n")
				outcon = False
				break

			Du = -f/J

			#Update the positions of the linkers. If the bead has rotated due to the small step taken during the gradient descent, this changes the positions p1 of the linker connections on the bead, which have to be updated.

			linkersold = copy.deepcopy(linkers)
			linkers.p0[:,0] = linkers.p0[:,0] - Du*par['a']
			linkers.p1[linkers.connected,:] = RotateBead.RotateBead(linkers.p1[linkers.connected,:],par,Du)
			linkers.edge[:,0] = linkers.edge[:,0] - Du*par['a']

			#Again, calculate the total force, torque, extension and x component of the force.
			F,T,x,Fx = TotalLinkerForceAndTorque.TotalLinkerForceAndTorque(linkers,par)
			alpha = 1

			#It can happen that due to the rotation, some linkers end up being stretched beyond there contour length Lmax. To prevent this, the step size is reduced step by step with a factor alphafactor until this no longer happens in the following while loop.
			while np.max(x)>1:
				#reduce step length if linkers end up overstretched
				alpha = alpha * alphafactor
				linkers = copy.deepcopy(linkersold)
				Du = Du*alpha 						#reduce step size
				linkers.p0[:,0] = linkers.p0[:,0] - Du*par['a']
				linkers.p1[linkers.connected,:] = RotateBead.RotateBead(linkers.p1[linkers.connected,:],par,Du)
				linkers.edge[:,0] = linkers.edge[:,0] - Du*par['a']
				
				F,T,x,Fx = TotalLinkerForceAndTorque.TotalLinkerForceAndTorque(linkers,par)

				if alpha<0.01:
					print("alpha too small: alpha = "+str(alpha))

			#After making the small rotation step, update the rotation and update the total torque.
			u = u+Du
			f = np.sum(T[linkers.connected])+TY

		#When the system is minimized, update the rotation and particle position.
		dx = u * par['a']
		dQ = u

	else:  #if par['slipconditions'] == 1 (all simulations in the article) we balance both Force and Torque, by allowing displacement steps and rotation steps during the optimization.
		#set the initial value of u = [DX,dQ]

		u= np.array([0,0])[np.newaxis].T
		f= np.array([np.sum(Fx[linkers.connected])+FX,np.sum(T[linkers.connected])+TY])[np.newaxis].T

		it = 0 #Set iteration counter
		outcon = True

		while np.linalg.norm(f)>TOL*np.max(FX,0):
		
			it = it + 1
			if it>MAXIT: #if the maximum nr of iterations is reached but the system is still not optimized, the simulation crashes.
				print("Maximum iterations reached")
				outcon = False
				break

			#Get Jacobian
			linkers1 = copy.deepcopy(linkers)

			#Apply ininetesimal movement dx
			linkers1.p0[:,0] = linkers1.p0[:,0]-SMALL

			F1,T1,x1,Fx1 = TotalLinkerForceAndTorque.TotalLinkerForceAndTorque(linkers1,par)				#Calculate the total force and torque on the slightly modified linkers
			dFdx = np.divide(np.sum(Fx1[linkers.connected]) - np.sum(Fx[linkers.connected]),SMALL)			#Calculate the slope in force at this slight modification
			dTdx = np.divide(np.sum(T1[linkers.connected]) - np.sum(T[linkers.connected]),SMALL)         	#Calculate the slope in torque at this slight modification

			#return to old value																			#Return the linkers p0 values to it's original values
			linkers1 = copy.deepcopy(linkers)

			linkers1.p1[linkers.connected,:] = RotateBead.RotateBead(linkers1.p1[linkers.connected,:],par,SMALL)	#Rotate the bead slightly (changin the p1 positions of the linkers)

			F1,T1,x1,Fx1 = TotalLinkerForceAndTorque.TotalLinkerForceAndTorque(linkers1,par)				#Calculate the total force and torque after this slight rotation
			dFdQ = np.divide(np.sum(Fx1[linkers.connected]) - np.sum(Fx[linkers.connected]),SMALL)			#Calculate the slope in force upon this slight rotation
			dTdQ = np.divide(np.sum(T1[linkers.connected]) - np.sum(T[linkers.connected]),SMALL)			#Calculate the slope in torque upon this slight rotation

			J = np.array([[dFdx,dFdQ],[dTdx,dTdQ]])															#Calculate the Matrix in changes of force and torque, upon the rotation and relocalisation

			if np.abs(np.linalg.det(J))<1e-10 or np.reciprocal(np.linalg.cond(J))<1e-10:			
				print("Singular Matrix\n")
				outcon=False
				break

			Du = np.linalg.lstsq((-1)*J,f,rcond=None)[0]													#Calculate the change in position and rotation necessary to decrease the total force f

			linkersold = copy.deepcopy(linkers)																#make another copy of the old set of linkers

			linkers.p0[:,0] = linkers.p0[:,0] - Du[0]																#Change the position of the linkers according to the step change in position
			linkers.p1[linkers.connected,:] = RotateBead.RotateBead(linkers.p1[linkers.connected,:],par,Du[1])		#Change the rotation of the linkers according to the step change in angle
			
			linkers.edge[:,0] = linkers.edge[:,0]-Du[0]																#Change the position of the edge according to the step change in position
			F,T,x,Fx = TotalLinkerForceAndTorque.TotalLinkerForceAndTorque(linkers,par)								#Calculate the total force and torque in the new position
			alpha = 1

			while np.max(x) > 1: 																					# If any of the linkers linker gets stretched beyond it's contour length, reduce the step length, to prevent overstretching	
				alpha = alpha * alphafactor																			#reduce the steplength by a factor alpha factor
				linkers = copy.deepcopy(linkersold)																	#make another copy of the linkers
				Du = Du*alpha 																						#Reduce the change in position and rotation according to the alpha factor
				linkers.p0[:,0] = linkers.p0[:,0] - Du[0]
				linkers.p1[linkers.connected,:] = RotateBead.RotateBead(linkers.p1[linkers.connected,:],par,Du[1]) 	#Change the rotation of the bead 
				linkers.edge[:,0] = linkers.edge[:,0] - Du[0] 														#Change the location of the edge 
				F,T,x,Fx = TotalLinkerForceAndTorque.TotalLinkerForceAndTorque(linkers,par)

				if alpha < 0.01:																					#If a viable step has not been reached when the step size becomes very small, return an error message
					print("Error:alpha too small, no valid step size could be found: alpha = "+str(alpha))	

			u = u + Du																								#Update the position and rotation of the bead
			f = np.array([np.sum(Fx[linkers.connected])+FX, np.sum(T[linkers.connected])+TY])[np.newaxis].T 		#Update the total forces

		#After Equilibration, store change in particle location dx and change in particle rotation dQ.
		dx = u[0][0]
		dQ = u[1][0]

	return(linkers,dx,dQ,F,outcon)