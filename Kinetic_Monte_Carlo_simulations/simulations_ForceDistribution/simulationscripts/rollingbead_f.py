from __future__ import division, unicode_literals, print_function

import numpy as np
import InitializeArray					
import ShearForceAndTorque				
import Equilibrate						
import TotalLinkerForceAndTorque
import RotateBead
import UpdateArray						
import RateConstants
import pandas as pd
import pickle as pk					

def start_simulation(par,path_to_simfolder):

	MaxSteps = par['MaxSteps']
	print_steps = par['print_steps']

	pi = np.pi

	#1. Initialize the simulation
	linkers = InitializeArray.InitializeArray(par)			#Calls the initialize script to determine the initial linkers that are bound, using the parameters defined above

	xc = 0									#Sets the initial position of the center of the bead
	theta = 0								#Sets the initial rotation angle of the bead

	#2. Apply Shear force on the particle and equilibrate the position and rotation angle
	if par['stress']>0:
		alphaF,alphaT = ShearForceAndTorque.ShearForceAndTorque(par['h'],par['a'])	#Calls the ShearForceAndTorque function to calculate the Shear force and the torque applied to a particle of this height and this size
		Fx = 6*pi*par['stress']*par['a']*par['h']*alphaF				#Calculate the horizontal (forward) component of the force on the bead
		Ty = 4*pi*par['stress']*(par['a']**3)*alphaT					#Calculate the y component of the torque on the bead
		linkers,dx,dQ,F,outcon = Equilibrate.Equilibrate(linkers,par,Fx,Ty)	#Calculate the equilibrate function to find the equilibrium position and rotation angle of the bead

		xc = xc+dx														#update the position of the bead
		theta = theta+dQ												#update the rotation angle of the bead
	t=0	
	Nass = 0															#Keeps track of the number of bonds formed during the simulation
	Ndiss = 0															#Keeps track of the number of bonds broken during the simulation


	#initialize arrays for storing all the data
	steps = np.zeros((MaxSteps,1))
	N_linkers = np.zeros((MaxSteps,1))
	X_positions = np.zeros((MaxSteps,1))
	Q_angles = np.zeros((MaxSteps,1))
	T = np.zeros((MaxSteps,1))
	Nass_lst = np.zeros((MaxSteps,1))
	Ndiss_lst = np.zeros((MaxSteps,1))

	#Set all arrays to zero
	steps[:]=np.nan
	N_linkers[:]=np.nan
	X_positions[:]=np.nan
	Q_angles[:]=np.nan
	T[:]=np.nan
	Nass_lst[:]=np.nan
	Ndiss_lst[:]=np.nan

	RandNrs = np.random.rand(MaxSteps,2)								#Draw an array of random numbers
	
	#3. Start Simulation
	for step in range(MaxSteps):

		if step % 10 == 0:
			steps[int(step/10)]=step
			N_linkers[int(step/10)]=np.sum(linkers.connected)
			X_positions[int(step/10)]=xc
			Q_angles[int(step/10)]=theta
			T[int(step/10)]=t
			Nass_lst[int(step/10)]=Nass
			Ndiss_lst[int(step/10)]=Ndiss

		if np.sum(linkers.connected)<2:
			print("\nAll LINKERS BROKEN\n")
			break
		
		#I. Add or remove attachment points (in case the bead has moved)
		linkers,F = UpdateArray.UpdateArray(linkers,F,par)

		if step%print_steps == 0:
			message = 'Number of steps= '+str(step)+', time = '+str(t)+', number of attachments= '+str(np.sum(linkers.connected))+', particle position= '+str(xc)

			#Store the linker positions into a pickle!
			linkersframe=np.hstack((linkers.p0,linkers.p1,linkers.connected[np.newaxis].T))
			linkerpath = path_to_simfolder+"linkers_traj/"+str(step)+".pickle"
			with open(linkerpath, 'wb') as pickle_file:
				pk.dump(linkersframe,pickle_file)

			edgeframe = linkers.edge[np.newaxis].T
			edgepath = path_to_simfolder+"edge_traj/"+str(step)+".pickle"
			with open(edgepath, 'wb') as pickle_file:
				pk.dump(edgeframe,pickle_file)

		k = RateConstants.RateConstants(linkers,par,F)						#calculate the rate constants (on and off) for each linker
		k[k==np.inf]=1e10;													#In case the rate constant approaches infinity, set a very high value of 1e10 for computational feasibility

		#II. Run the Gillespy KMC Algorithm
		dt = np.divide(-1*np.log(RandNrs[step,0]),np.sum(k))		#Draw the waiting time until next reaction

		#Next, find out which of the reactions takes place. First, remove the reactions with zero propensity, as these may mess up the simulation, and have no chance of occuring
		valid_indx = k>0										#Find the indices referring to valid reactions
		valid_k = k[valid_indx]									#use only the valid reactions
		reactions_valid = np.where(valid_indx==True)[0]			#find the indices of the reactions which are valid

		#Next, construct the intervals
		selection_intv1 = np.cumsum(valid_k)					#Cumulative sum
		selection_intv1 = np.divide(selection_intv1,selection_intv1[-1]) #Normalize the array by dividing it by the last element
		selection_ind = np.where(selection_intv1>RandNrs[step,1])[0][0]

		j = reactions_valid[selection_ind]						#Selected tether to modify

		#Now, apply the modification!
		if linkers.connected[j]==False:							#No tether present? connect it!
			linkers.connected[j]=True
			Nass = Nass+1

			ri2 = np.square(linkers.p0[j,0])+np.square(linkers.p0[j,1])

			if ri2 == 0:
				ri2 = 1e-14
			if par['attachmentdirection']==0:					#If the set attachmentdirection is plain vertical:
				zi = par['h']-np.sqrt(np.square(par['a'])-ri2)
				linkers.p1[j,:]=np.append(linkers.p0[j,:],zi)

			else:												#If the set attachment direction is not vertical, but normal to the particle surface:
				ri1 = np.sqrt(ri2)
				alpha = np.arctan(ri1/par['h'])
				r1 = par['a']*np.sin(alpha)
				linkers.p1[j,0:2]=np.multiply(linkers.p0[j,:],(r1/ri1))
				z1=par['h']-par['a']*np.cos(alpha)
				linkers.p1[j,2]=z1

		else:													#Tether present? rupture it!
			linkers.connected[j]=False
			linkers.p1[j,:]=[0.,0.,0.]
			Ndiss=Ndiss+1

		#Equilibrate the forces again
		linkers,dx,dQ,F,outcon = Equilibrate.Equilibrate(linkers,par,Fx,Ty) #outcon?
		if not outcon:
			break

		xc=xc+dx												#Update the position of the particle
		theta = theta+dQ										#Update the angle of the particle
		t=t+dt													#Update the time
		
	return(outcon,steps,T,N_linkers,X_positions,Q_angles,Nass_lst,Ndiss_lst)

