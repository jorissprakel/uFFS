# Microfluidic Force Spectroscopy (uFFS) Data Analysis and Simulation Scripts.

This depository contains the scripts used in the paper "Rapid Molecular Mechanotyping using Microfluidic Force Spectroscopy", 
by Martijn van Galen, Annemarie Bok, Taieesa Peshkovsky, Jasper van der Gucht, Bauke Albada and Joris Sprakel.

Please note that the scripts are annotated with comments, to instruct the user.
Both scripts used for the uFFS experimental data analysis and scripts used to perform and analyse the Kinetic Monte Carlo Simulations are provided:

-------------------------- uFFS experimental data analysis scripts --------------------

The uFFS_experimental_data_analysis/ folder contains a sample dataset of one measurement, and the data analysis scripts used in the article.
To analyse a uFFS dataset, the following analysis scripts should be run in this order:

1.preprocess_and_track: 	This script locates the dataset, carries out the pre-processing as described in SI section S5, tracks the particles in the dataset and stores the trajectories.

2.overlay_traces:		This script makes plots in which the particle trajectories are overlayed on the dataset. This allows the user to check if the particle tracking went well.

3.compute_dissociationtimes:	This script Computes the dissociation time of each particle from the trajectories

4.Plot_fit_Ptcurves:		This script plots the P-t curves using the dissociation times and fits them with a stretched exponential function to compute the lifetimes. It also computes the applied shear force from the flowrate data.

5.Plot_Force_lifetimecurves:	This script plots the force-lifetime curves, using the shearforces and lifetimes determined by 4.Plot_fit_Ptcurves.py.

Please note that the provided sample dataset is not a full measurement, so scripts 4 and 5 do not work on this dataset.


----------------------------- Kinetic Monte Carlo Simulations -------------------------

Contents:

-example_simulation_folder: test example of a typical simulation series

-simulations_varyingdensity: Contains all scripts required to perform the simulation series with varying bond density (Fig. 5B)

-simulations_varyingBasepairs: Contains all scripts required to perform the simulation series with varying numbers of basepairs (Fig. 5C)

-simulations_ForceDistribution: Contains all scripts required to perform the simulation series for the force distribution heatmaps (Fig. 5D,E)

Here follows a description on how to perform the simulations, along with an overview of all the scripts in the folders mentioned above.
Each simulation series contains the following subfolders:


Simulation_series_name/
	simulationscripts/
	simulation/
	analysis/

All scripts required to run the simulations are stored in simulationscripts/
During the simulations, raw simulation results will be stored in the simulation/ folder,
The analysisscripts/ folder contains all the scripts required to analyse the raw data.
and processed results after running the analysisscripts will be stored in the analysis/ folder.

Here, we will go over all the simulation scripts in the simulationscripts/ folder required to 1. generate the input parameters for a simulation series and 2. run the simulations. 

--- 1. Generate input parameters for the simulations ---

	-generate_measurearchitecture.py: This file allows the user to specify all the simulation parameters, and generates a number of input files that can be read by the run_simulation.py script. In the top of the file, a list of variable parameters is specified: The single-bond mechanochemical parameters (dx,koff,0), the applied flowrates (=shear stress), and the number of repeats. Below that, the remaining constant parameters is specified, such as the size of the particle, kT, etc. When run using "python generate_measurearchitecture.py", this script creates the ../simulation/ folder and generates an input.txt file containing the parameters for each individual simulation. Furthermore, a simulationslist.txt file is generated, which is a list of all the simulations that need to be performed.

	-generate_input.py: This script is called by generate_measurementarchitecture.py

	-generate_simulationslist.py: This script is called by generate_measurearchitecture.py

--- 2. Running the simulations ---
The main simulation script is run_simulation.py
This script is used to start a simulation. An individual simulation can be run with this script, using "python ../simulation/path_to_simulationfolder/ input.txt"

-rollingbead_f.py:
This is the main KMC simulation function, which is called by run_simulation.py. An explanation of all the steps is given in the script. rollingbead_f.py calls several other scripts that perform part of the simulation:

	-InitializeArray.py: 		This function creates initial positions for the linkers on the channel surface. These positions are chosen such that the bond density rho is met. Next, the function randomly allows some linkers to form.
	
	-ShearForceAndTorque.py: 	This function computes the Force and Torque on the particle caused by the fluid shear.
	
	-Equilibrate.py:		This function equilibrates the free energy of the system by rotating and translating the bead using a grandient descent approach
	
	-RotateBead.py: 	This function is called by Equilibrate.py. It is used to recalculate the positions p1 of the linkers on the particle as the particle rotates.
	
	-TotalLinkerForceAndTorque.py:	This function is called by both Equilibrate.py and rollingbead_f.py. It computes the total extensional force F and torque T all bound linkers apply to the particle, the relative extension x, and finaly the x component of the force F. 
	
	-UpdateArray.py:		This function is called by rollingbead_f.py at the start of each simulation loop. It creates and removes positions of linkers on the channel surface in case the bead has moved during the previous step. Parts of UpdateArray.py are similar to InitializeArray.py
	
	-RateConstants.py:		This function is called by rollingbead_f.py before each kinetic monte carlo step, to compute the rate constants for formation and dissociation of every linkers.

The batch of simulations generated by generate_measurearchitecture.py can be run using the SLURM workload manager "https://slurm.schedmd.com/documentation.html", by executing the initiate_jobarrayparallel.sh file.
This is done using the command line: "sbatch initiate_jobarrayparallel.sh simulationslist.txt"

The initiate_jobarrayparallel.sh script reads in the list of simulations specified in simulationslist.txt, and then runs all simulations. Please make sure a sufficient number of tasks (--array command) is specified at the top of the initiate_jobarrayparallel.sh script.

Alternatively, individual simulations can be run separately by running python from the command line, with the command: "python run_simulation.py ../simulation/dil_X/stressX/repX/ input.txt".

Finally, we go over all the scripts used to analyse the simulated data. These scripts are found in the analysisscripts/ folder of each simulation series.

	-plot_store_force_lifetimecurve.py:	This function reads simulation results (.log files) and plots the force-lifetime curve.
	
	-overlay_force_lifetimecurve.py:	This function overlays the force-lifetime over the experimental uFFS data, to generate the figures shown in fig. 5B and 5C.
	-plot_heatmap_forcedistribution_averaged.py:	This script is only present in the simulations_ForceDistribution/ folder. It computes the heatmaps of the force distributions on the linkers, as shown in figures 5D and 5E of the uFFS paper.
	
