from __future__ import division, unicode_literals, print_function
import numpy as np
import rollingbead_f
import pandas
import pickle
import sys
import os

path_to_simfolder = sys.argv[1]									#Specify the path to the location of the input file as argumnent 1
inputfile_name = sys.argv[2]									#Specify the name of the inputfile (typically input.txt) as argument 2
if not os.path.exists(path_to_simfolder+"/linkers_traj/"):
	os.makedirs(path_to_simfolder+"linkers_traj/")

if not os.path.exists(path_to_simfolder+"/edge_traj/"):
	os.makedirs(path_to_simfolder+"edge_traj/")

#Read in the parameters from the inputfile
par={}
inputfile = open(path_to_simfolder+inputfile_name,'r')
for line in inputfile:
	name,value = line.split(" = ")
	if name in ['attachmentdirection','extralinkers','slipconditions','polymermodel','MaxSteps','print_steps','seed']:
		value = int(value)
	elif name in ['trajectoryfile','parametersfile','logfile']:
		value = str(value)[:-1]
	else:
		value = float(value)
	par[name]=value

np.random.seed(par['seed'])

outcon,steps,T,N_linkers,X_positions,Q_angles,Nass_lst,Ndiss_lst=rollingbead_f.start_simulation(par,path_to_simfolder)	#Performs the KMC simulations
StorageMatrix = np.hstack((steps,T,N_linkers,X_positions,Q_angles,Nass_lst,Ndiss_lst))									#Once the simulation is finished, all the log.csv parameters (position, time etc. are stored in an array)

steps = steps[~np.isnan(steps)]
datalength = len(steps)
StorageMatrix = StorageMatrix[0:datalength,:]

#Check whether the particle dissociated from the surface during the simulation 
if (outcon or N[:-1]<3):
	dissociatedlst = np.ones((len(steps),1))
else:
	dissociatedlst = np.zeros((len(steps),1))

#Finally, the simulation results are stored as log.csv.
StorageMatrix = np.hstack((StorageMatrix,dissociatedlst))
np.savetxt(path_to_simfolder+"log.csv",StorageMatrix,delimiter=',',header ='nsteps,time,N_linkers,x_bead,theta_bead,nr_Associated,nr_Detached,detached during simulation(1=yes;0=bead remained attached')
