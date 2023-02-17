import numpy as np
import os
import generate_simulationslist
import generate_input
import random

########MAIN PARAMETERS OF THE SIMULATION SERIES##############


flowrates = np.array([19.0])
flowrate_to_stress = 6.67440969
stresses = np.round(np.multiply(flowrates,flowrate_to_stress),2)

nr_simulations = 20									#Number of repeated simulations per stress
MaxSteps = 100000									#Number of simulation steps 
print_steps = 10									#Print an update on the linker force and positions every print_steps
dilutionfactors = [1]								#Multiplied by the parameter "rho" to determine the number of linkers per um2 as dilutionfactor*rho

database = {}
database.update({"stresses":stresses})
unitvector = np.ones(len(database["stresses"]))

nr_repeats = unitvector*nr_simulations
nr_repeats = nr_repeats.astype("int")

database.update({"MaxSteps":unitvector*MaxSteps})
database.update({"print_steps":unitvector*print_steps})
database.update({"dilutionfactors":dilutionfactors})

##################################### ########################
#Other parameters not changed often:

parameters_standard = {}
parameters_standard['MaxSteps'] = None
parameters_standard['print_steps'] = None

parameters_standard['kT'] = 4.1e-21					#kBT
parameters_standard['a'] = 0.5*4.34 				#Particle radius [um]
parameters_standard['stress'] = None  				#Shear stress at the channel wall =eta*gammadot
parameters_standard['h']=parameters_standard['a'] 	#Height of the center of the bead (must be larger than a)
parameters_standard['rho']= 900 					#linkers per um2 on the bottom surface
													#1: slip, roll plus slip: balance horzontal force and torque
parameters_standard['dilutionfactor'] = None													
parameters_standard['linkerpositionnoise'] = 0.1 	#relative fluctuations in linker positions, |dr|/r with dr variation, and r average distance between linkers
													#fluctuation is only on position on bead!
parameters_standard['attachmentdirection'] = 1 		#0: initial linkers are vertical; 1: normal to bead surface
parameters_standard['extralinkers']=15				#include extra tethers on the right so bead can roll

parameters_standard['slipconditions'] = 1 			#0: no slip, pure rolling: only balance torque, and assume Dx=theta*a

#Set linker mechanical properties
parameters_standard['Lmax'] = 0.065					#contour length of the DNA linker [um]
parameters_standard['lK'] = 0.050					#Kuhn length DNA linker [um]
parameters_standard['polymermodel'] = 2				#1 = ideal chain, 2 = WLC

#Bond properties
parameters_standard['kon0']	= 1e6 					#kon in absence of force
parameters_standard['koff0'] = 1e-12				#parameters for the 30 BP bond. [s-1] This is koff0 in absence of force
parameters_standard['del'] = 0.0028					#parameters for the 30 BP bond. [um] activation length

parameters_standard['trajectoryfile']="trajectory.csv"
parameters_standard['parametersfile']="parameters.txt"
parameters_standard['logfile']="logfile.txt"

generate_simulationslist.generate_simulationslist(dilutionfactors=database['dilutionfactors'],stresses=database["stresses"],nr_simulations=nr_simulations)

for dilutioncount,dilutionfactor in enumerate(database['dilutionfactors']):
	for stressindex,stressvalue in enumerate(database['stresses']):
		for repeatnr in range(nr_repeats[stressindex]):
			os.makedirs("../simulation/dil_"+str(dilutionfactor)+"/stress_"+str(stressvalue)+"/rep_"+str(repeatnr)+"/")

			parameters_measurement = parameters_standard.copy()
			parameters_measurement.update({"seed":random.randint(0,100000)})
			parameters_measurement.update({"stress":stressvalue})
			parameters_measurement.update({"MaxSteps":MaxSteps})
			parameters_measurement.update({"print_steps":print_steps})
			parameters_measurement.update({"dilutionfactor":dilutionfactor})

			generate_input.generate_input(parameterlib=parameters_measurement,outputfilename="../simulation/dil_"+str(dilutionfactor)+"/stress_"+str(stressvalue)+"/rep_"+str(repeatnr)+"/input.txt")
