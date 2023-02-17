import numpy as np
import os
import generate_simulationslist
import generate_input
import random

########MAIN PARAMETERS OF THE SIMULATION SERIES##############

flowrates = np.array([49.0,48.4,47.8,47.2,46.6,46.0,45.4,44.8,44.2,43.6,43.0,42.4,41.8,41.2,40.6,40.0,39.4,38.8,38.2,37.6,37.0,36.4,35.8,35.2,34.6,34.0,33.4,32.8,32.2,31.6,31.0,30.4,29.8,29.2,28.6])			#uL/min

bonds_per_um2 = np.array([2000])

flowrate_to_stress = 6.67440969						#kbt/um3 with channel dimensions of 50 um height, 1500 um wide, 2 cm long, in water at 20 C
stresses = np.round(np.multiply(flowrates,flowrate_to_stress),2)

nr_simulations = 20									#Number of repeated simulations per stress
MaxSteps = 10000000									#Number of simulation steps 
print_steps = 9999999									#Print an update on the state of the simulation every x steps

database = {}
database.update({"stresses":stresses})
unitvector = np.ones(len(database["stresses"]))

nr_repeats = unitvector*nr_simulations
nr_repeats = nr_repeats.astype("int")

database.update({"MaxSteps":unitvector*MaxSteps})
database.update({"print_steps":unitvector*print_steps})
database.update({"bonds_per_um2":bonds_per_um2})

##################################### ########################
#Other parameters not changed often:

parameters_standard = {}
parameters_standard['MaxSteps'] = None
parameters_standard['print_steps'] = None

parameters_standard['kT'] = 4.04737e-21					#kBT
parameters_standard['a'] = 0.5*4.34 				#Particle radius [um]
parameters_standard['stress'] = None  				#Shear stress at the channel wall =eta*gammadot
parameters_standard['h']=parameters_standard['a']	#Height of the center of the bead (must be larger than a)
parameters_standard['rho']= None					##linkers per um2 on the bottom surface
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
parameters_standard['koff0'] = 3.16e-10				#koff0 in absence of force
parameters_standard['del'] = 0.00245				#activation length (in Bell equation) [um]

parameters_standard['trajectoryfile']="trajectory.csv"
parameters_standard['parametersfile']="parameters.txt"
parameters_standard['logfile']="logfile.txt"

generate_simulationslist.generate_simulationslist(bonddensity=database["bonds_per_um2"],stresses=database["stresses"],nr_simulations=nr_simulations)

for bonddensitycount,bonddensity in enumerate(database['bonds_per_um2']):
	for stressindex,stressvalue in enumerate(database['stresses']):
		for repeatnr in range(nr_repeats[stressindex]):
			os.makedirs("../simulation/density_"+str(bonddensity)+"/stress_"+str(stressvalue)+"/rep_"+str(repeatnr)+"/")

			parameters_measurement = parameters_standard.copy()
			parameters_measurement.update({"seed":random.randint(0,100000)})
			parameters_measurement.update({"stress":stressvalue})
			parameters_measurement.update({"MaxSteps":MaxSteps})
			parameters_measurement.update({"print_steps":print_steps})
			parameters_measurement.update({"rho":bonddensity})

			generate_input.generate_input(parameterlib=parameters_measurement,outputfilename="../simulation/density_"+str(bonddensity)+"/stress_"+str(stressvalue)+"/rep_"+str(repeatnr)+"/input.txt")
			