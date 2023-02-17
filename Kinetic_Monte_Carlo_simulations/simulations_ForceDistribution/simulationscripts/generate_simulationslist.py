#This file takes a number of looping parameters, in this case basepairs, stresses applied, and nr of repeats per simultion.
#Next, a simulationslist.txt file is written, Which can be read by the SLURM task manager to start all the simulations.

def generate_simulationslist(basepairs,stresses,nr_simulations):

	print(nr_simulations)
	simulationslst = range(0,nr_simulations)

	lstfile = open('simulationslist.txt','w')
	for BP in basepairs:
		for elem in stresses:
			for nr in simulationslst:
				lstfile.write(str(BP)+' '+str(elem)+' '+str(nr)+'\n')
	
	lstfile.close()
