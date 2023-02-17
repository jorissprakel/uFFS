import numpy as np
import os
import generate_simulationslist
import generate_input
import random

########    MAIN PARAMETERS THAT CAN BE VARIED THROUGHOUT THE SIMULATION SERIES    ##################

flowrates = np.array([24,23.5,23,22.5,22,21.5,21,20.5,20,19.5,19,18.5,18,17.5,17,16.5,16,15.5,15,14.5,14,13.5,13,12.5,12,11.5,11,10.5,10])		#Applied flowrates

basepairs = np.array([15])											#nr of basepairs simulated
koff0_values = np.array([3.16e-5])									#[s-1] Bond off-rates in the absence of Force, values from Strunz et al. 1999
del_values = np.array([0.00175])									#[um] Bond activation length, values from Strunz et al. 1999
rho_values = [900]													#[linkers/um2] density of linkers on the channel surface

flowrate_to_stress = 6.67440969										#Conversion factor from flowrate Q (uL/min) to stress (kbt/um3). This factor follows from sigma = eta*6Q/h^2w (Main txt, Eq3.)
stresses = np.round(np.multiply(flowrates,flowrate_to_stress),2)

nr_simulations = 20													#Number of repeats generated for every simulation
MaxSteps = 10000000													#Number of simulation steps per simulation
print_steps = 9999999												#Number of simulation steps after which an update of the linker position and extension is stored. print_steps can be set to MaxSteps-1 for simulations that only store information on the dissociation time, but no information on the linker position.


database = {}														#Create a database to store all the variable information
database.update({"stresses":stresses})								#Update all the variable information
unitvector = np.ones(len(database["stresses"]))
nr_repeats = unitvector*nr_simulations
nr_repeats = nr_repeats.astype("int")
database.update({"MaxSteps":unitvector*MaxSteps})
database.update({"print_steps":unitvector*print_steps})
database.update({"basepairs":basepairs})
database.update({"koff0":koff0_values})
database.update({"delvalues":del_values})
database.update({"rhovalues":rho_values})

##################################### ########################
#############    OTHER SIMULATION PARAMETERS NOT VARIED THROUGHOUT THE SIMULATION SERIES    ##########
parameters_standard = {}							#Create a parameters dictionary, this dictionary should contain all parameters used during the simulation, and "None" values for the parameters specified in the database above
parameters_standard['MaxSteps'] = None				
parameters_standard['print_steps'] = None

parameters_standard['kT'] = 4.04737e-21				#Value of the constant kBT
parameters_standard['a'] = 0.5*4.34 				#[um] Particle radius 
parameters_standard['stress'] = None  				#Shear stress at the channel wall =eta*gammadot
parameters_standard['h']=parameters_standard['a']	#Height of the center of the bead. In these simulations set to a: the particles touch the channel surface exactly
parameters_standard['rho']= None					#linker density
													
parameters_standard['linkerpositionnoise'] = 0.1 	#relative fluctuations in linker positions, defined as |dr|/r with dr variation, and r average distance between linkers
parameters_standard['attachmentdirection'] = 1 		#0: initial linkers are vertical; 1: normal to bead surface
parameters_standard['extralinkers']=15				#include extra tethers on the right so new bonds can be formed if the bead moves. By standard, this is set to 15.

parameters_standard['slipconditions'] = 1 			#0: no slip, pure rolling: only balance torque, and assume Dx=theta*a
													#1: slip, roll plus slip: balance horzontal force and torque. This was set to 1 for all simulations in this paper.
#Set linker mechanical properties
parameters_standard['Lmax'] = 0.065					#[um] contour length of the DNA linker 
parameters_standard['lK'] = 0.050					#[um] Persistence length DNA linker 
parameters_standard['polymermodel'] = 2				#Sets polymer extension model 1 = ideal chain, 2 = WLC. WLC was used in all simulations in this paper

#Bond properties (for now just slip bonds)
parameters_standard['kon0']	= 1e6 					#[s-1] Base bond formation rate kon in absence of force
parameters_standard['koff0'] = None					#[s-1] koff0 in absence of force
parameters_standard['del'] = None					#[um] activation length (in Bell equation) 

parameters_standard['trajectoryfile']="trajectory.csv" #Name of trajectory file stored. Trajectory stores all the linker positions and extensions every print_steps simulatio nsteps
parameters_standard['parametersfile']="parameters.txt" #Name of parameter file stored. Stores all the parameters defined above, for later reference.
parameters_standard['logfile']="logfile.txt"		   #Name of the log file stored. This file stores all parameters such as the number of steps, time passed, number of active linkers, the x-position of the particle, the rotation angle theta of the particle, and the number of formed and broken bonds

generate_simulationslist.generate_simulationslist(basepairs=database["basepairs"],stresses=database["stresses"],nr_simulations=nr_simulations) #Generates a list of simulations to be run, for the SLURM processor to loop over.

#The for-loops indicates which parameters are looped over to form simulation inputs. In this case, these are: The nr of basepairs, the number of flowrates (=stress values), and the nr of repeats of each measurement
for bpcount,basepairs in enumerate(database["basepairs"]):
	for stressindex,stressvalue in enumerate(database['stresses']):
		for repeatnr in range(nr_repeats[stressindex]):
			os.makedirs("../simulation/BPs_"+str(basepairs)+"/stress_"+str(stressvalue)+"/rep_"+str(repeatnr)+"/") #For each simulation, a folder is created in the ../simulation directory

			#Fill in the simulation parameters from database into the standard parameters from parameters_standard to generate the final list of input parameters
			parameters_measurement = parameters_standard.copy()
			parameters_measurement.update({"seed":random.randint(0,100000)}) #Seed the random number generator
			parameters_measurement.update({"stress":stressvalue})
			parameters_measurement.update({"MaxSteps":MaxSteps})
			parameters_measurement.update({"print_steps":print_steps})
			parameters_measurement.update({"koff0":koff0_values[bpcount]})
			parameters_measurement.update({"del":del_values[bpcount]})
			parameters_measurement.update({"rho":rho_values[bpcount]})

			#Finally, write the final list of input parameters in the appropriated simulation folder
			generate_input.generate_input(parameterlib=parameters_measurement,outputfilename="../simulation/BPs_"+str(basepairs)+"/stress_"+str(stressvalue)+"/rep_"+str(repeatnr)+"/input.txt")
			