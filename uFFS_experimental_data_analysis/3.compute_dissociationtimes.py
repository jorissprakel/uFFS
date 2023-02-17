from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import scipy.stats
import scipy.optimize as fit
import os
import pims
import trackpy as tp
from scipy.ndimage import gaussian_filter
from matplotlib.ticker import AutoMinorLocator,Locator

#This script takes the particle trajectories, and computes the dissociation time of every particle.

datafolder = 'sampledata_notfulldataset/dataset/' #Specify the location of the data again.
framespersecond = 40 							#Frames per second used in the datacollection (used to calculate velocities)
channelwidth_um = 895							#width of the channel in micrometers. This is used to convert distances in pixels to distances in um.

compute_MSD_before = True 						#Calculate the mean square displacement in the time prior to turning on the flowrate, using the lag time below.
MSD_lagtime = 1 								#seconds, the lag time used in the MSD calculations
time_before_flow = 5 							#seconds, marks the end of the time window until the flow rate was turned on. 
dissociationthreshold_MSD = 20 					#um^2 threshold MSD value used to determine whether a particle is detached or not

#List all the measurements in the datafolder and loop over them.
datapoints = os.listdir(datafolder)
for dataindex,datapoint in enumerate(datapoints):

	print("--- Analysing measurement: "+datafolder+datapoint+" ---")

	trackingparameters=np.loadtxt(datafolder+datapoint+"/analysis/trackingparameters.csv",skiprows=1,delimiter=',')
	framestep = int(trackingparameters[0])
	time_per_frame = framestep/framespersecond #time difference between two frames during this analysis

	channelleft = int(trackingparameters[1])
	channelright = int(trackingparameters[2])
	channelwidth_px = channelright-channelleft
	um_per_px = channelwidth_um/channelwidth_px #Calculate the width of one pixel in um, for converting distances in pixels to distances in um.

	particleradius = float(trackingparameters[4])
	
	#Read the particle trajectories.
	t = pd.read_pickle(datafolder+datapoint+"/analysis/trajectories")

	#Analyse only those particles that were also found in frame 0 of the movie (The particles that were stuck to the channel at the beginning, before the flowrate is turned on.)
	t_particles = t.loc[t['frame']==0] 								
	particlelist = t_particles.get('particle').to_numpy()
	t_analyse = t.loc[t['particle'].isin(particlelist)]


	particlesfound = int(t_analyse.max()['particle'])
	analysedcounter =0
	percentage_old = None

	xpositions = []
	ypositions = []
	averagevelocities = []
	dissociationtimes = []
	MSDvalues = []

	#loop over all particles found in this movie
	for pindex in range(particlesfound):
		percentage = str(round(analysedcounter/particlesfound*100))		
		percentage_old = percentage

		ptraj=t_analyse.loc[t_analyse['particle']==pindex]
		xvals = ptraj.get('x').to_numpy()
		yvals = ptraj.get('y').to_numpy()
		frames = ptraj.get('frame').to_numpy()

		analysedcounter+=1

		#Use numpy arrays to calculate all the instantaneous velocities for this particle in all frames in one go.
		xvals_offset = xvals[1:]
		yvals_offset = yvals[1:]
		frames_offset = frames[1:]

		xdisplacement = xvals_offset-xvals[:-1]
		ydisplacement = yvals_offset-yvals[:-1]
		framedifference = frames_offset-frames[:-1]

		xsquared = np.power(xdisplacement,2)
		ysquared = np.power(ydisplacement,2)

		absdisplacement = np.sqrt(xsquared+ysquared)

		absdisplacement_um = np.multiply(absdisplacement,um_per_px)

		framedifference_s = np.multiply(framedifference,time_per_frame)
		times_s = np.multiply(frames,time_per_frame)
		velocity_um_per_s = np.divide(absdisplacement_um,framedifference_s)
		
		#Store the x and y position to keep track of where the particle was located
		xpos_um = np.multiply(xvals[0],um_per_px)
		ypos_um = np.multiply(yvals[0],um_per_px)
		
		#Determine the particles displacement versus initial point:
		xdispl_vinitial = xvals[:]-xvals[0]
		ydispl_vinitial = yvals[:]-yvals[0]
		xdisplvi_sq = np.power(xdispl_vinitial,2)
		ydisplvi_sq = np.power(ydispl_vinitial,2)
		absdisplacement_vi = np.sqrt(xdisplvi_sq+ydisplvi_sq)
		absdisplacement_vi_um = np.multiply(absdisplacement_vi,um_per_px)

		#The next part calculates the MSD values, if requested.
		if compute_MSD_before == True:

			xvals_um = np.multiply(xvals,um_per_px)
			yvals_um = np.multiply(yvals,um_per_px)

			xvals_um_shifted = xvals_um[int(MSD_lagtime/time_per_frame):] #Shift the positions by the set lag time*frames per second.
			yvals_um_shifted = yvals_um[int(MSD_lagtime/time_per_frame):] #Shift the positions by the set lag time*frames per second.

			length_shifted = len(xvals_um_shifted)

			xvals_um_forMSD = xvals_um[:length_shifted]
			yvals_um_forMSD = yvals_um[:length_shifted]

			displacements_x = np.absolute(xvals_um_shifted-xvals_um_forMSD)
			displacements_y = np.absolute(yvals_um_shifted-yvals_um_forMSD)

			totaldisplacement = np.sqrt(np.power(displacements_x,2)+np.power(displacements_y,2))
			squaredisplacement = np.power(totaldisplacement,2)
			
			times_MSDwindow = times_s[:length_shifted]
			
			#Store the MSD value before the flow is turned on, by averaging the squared displacements before flow.
			MSD_before_flowrate = np.average(squaredisplacement[:int(time_before_flow/time_per_frame-MSD_lagtime/time_per_frame)])
			MSDvalues.append(MSD_before_flowrate)

		times_s_corrected = times_s[:]
		absdisplacement_um_corrected=absdisplacement_um[:]
		velocity_um_per_s_corrected=velocity_um_per_s[:]

		detachthreshold = dissociationthreshold_MSD
		squaredisplacement_flipped = np.flip(squaredisplacement)
		inversedetachmentpoint = np.where(squaredisplacement_flipped<detachthreshold)[0]#The first point where the MSD is below the threshold
		
		if len(inversedetachmentpoint)==0:
			detachmentpoint_corr=0

		else:
			inversedetachmentpoint=inversedetachmentpoint[0]
			detachmentpoint = len(absdisplacement_um)-inversedetachmentpoint #the inverse of this threshold point
			detachmentpoint_corr = int(detachmentpoint+(MSD_lagtime/time_per_frame)) #Correct for the lag time to find the actual detachment point.

		if detachmentpoint_corr>len(absdisplacement_um)-1:
			detachmentpoint_corr=len(absdisplacement_um)-1#if the detachmentpoint projected in this way is greater than the measured time, set it to the maximum measured time: i.e. this indicates that the particle did not dissociate over the course of the measurement.

		dissociationtime=times_s_corrected[int(detachmentpoint_corr)]

		#Store the information about this particle in an array.
		dissociationtimes.append(times_s_corrected[detachmentpoint_corr])
		averagevelocities.append(np.average(velocity_um_per_s_corrected[detachmentpoint_corr:]))
		xpositions.append(xpos_um)
		ypositions.append(ypos_um)

	#Store the calculated dissociationtimes and MSD values in a csv file.
	resultmatrix=np.transpose(np.vstack((xpositions,ypositions,dissociationtimes,averagevelocities,MSDvalues)))
	np.savetxt(datafolder+datapoint+"/analysis/traceanalysisresults.csv",resultmatrix,delimiter=",",header ='xposition_particle (um),yposition_particle (um),dissociationtime (s),averagevelocity_afterdissociation(um/s),MSDvalues(um2)')
