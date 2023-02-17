import numpy as np
import os
import generate_simulationslist
import generate_input
import random

########MAIN PARAMETERS OF THE SIMULATION SERIES##############

flowrates = np.array([25,24.7,24.4,24.1,23.8,23.5,23.2,22.9,22.6,22.3,22.0,21.7,21.4,21.1,20.8,20.5,20.2,20.0,19.7,19.4,19.1,18.8,18.5,18.2,17.9,17.6,17.3,17.0,16.7,16.4,16.1,15.8,15.5,15.2,14.9,14.6,14.3,14.0])

basepairs = np.array([30])
koff0_values = np.array([1e-12])
del_values = np.array([0.0028])
rho_values = [900]

flowrate_to_stress = 6.67440969						#kbt/um3 with channel dimensions of 50 um height, 1500 um wide, 2 cm long, in water at 20 C
stresses = np.round(np.multiply(flowrates,flowrate_to_stress),2)

nr_simulations = 20									#Number of repeated simulations per stress
MaxSteps = 10000000									#Number of simulation steps 
print_steps = 9999999								#Print an update on the state of the simulation every x steps

database = {}
database.update({"stresses":stresses})
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
#Other parameters not changed often:

parameters_standard = {}
parameters_standard['MaxSteps'] = None
parameters_standard['print_steps'] = None

parameters_standard['kT'] = 4.04737e-21				#kBT
parameters_standard['a'] = 0.5*4.34 				#Particle radius [um]
parameters_standard['stress'] = None  				#Shear stress at the channel wall =eta*gammadot
parameters_standard['h']=parameters_standard['a']	#Height of the center of the bead (must be larger than a)
parameters_standard['rho']= None					#linkers per um2 on the bottom surface
													#1: slip, roll plus slip: balance horzontal force and torque
parameters_standard['linkerpositionnoise'] = 0.1 	#relative fluctuations in linker positions, |dr|/r with dr variation, and r average distance between linkers
													#fluctuation is only on position on bead!
parameters_standard['attachmentdirection'] = 1 		#0: initial linkers are vertical; 1: normal to bead surface
parameters_standard['extralinkers']=15				#include extra tethers on the right so bead can roll

parameters_standard['slipconditions'] = 1 			#0: no slip, pure rolling: only balance torque, and assume Dx=theta*a

#Set linker mechanical properties
parameters_standard['Lmax'] = 0.065					#contour length of the DNA linker [um]
parameters_standard['lK'] = 0.050					#Kuhn length DNA linker [um]
parameters_standard['polymermodel'] = 2				#1 = ideal chain, 2 = WLC

#Bond properties (for now just slip bonds)
parameters_standard['kon0']	= 1e6 					#kon in absence of force
parameters_standard['koff0'] = None					#koff0 in absence of force
parameters_standard['del'] = None					#activation length (in Bell equation) [um]

parameters_standard['trajectoryfile']="trajectory.csv"
parameters_standard['parametersfile']="parameters.txt"
parameters_standard['logfile']="logfile.txt"

generate_simulationslist.generate_simulationslist(basepairs=database["basepairs"],stresses=database["stresses"],nr_simulations=nr_simulations)

for bpcount,basepairs in enumerate(database["basepairs"]):
	for stressindex,stressvalue in enumerate(database['stresses']):
		for repeatnr in range(nr_repeats[stressindex]):
			os.makedirs("../simulation/BPs_"+str(basepairs)+"/stress_"+str(stressvalue)+"/rep_"+str(repeatnr)+"/")

			parameters_measurement = parameters_standard.copy()
			parameters_measurement.update({"seed":random.randint(0,100000)})
			parameters_measurement.update({"stress":stressvalue})
			parameters_measurement.update({"MaxSteps":MaxSteps})
			parameters_measurement.update({"print_steps":print_steps})
			parameters_measurement.update({"koff0":koff0_values[bpcount]})
			parameters_measurement.update({"del":del_values[bpcount]})
			parameters_measurement.update({"rho":rho_values[bpcount]})

			generate_input.generate_input(parameterlib=parameters_measurement,outputfilename="../simulation/BPs_"+str(basepairs)+"/stress_"+str(stressvalue)+"/rep_"+str(repeatnr)+"/input.txt")
