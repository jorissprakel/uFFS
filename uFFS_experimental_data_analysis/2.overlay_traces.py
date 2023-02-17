from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import scipy.optimize as fit
import os
import time
import pims
import trackpy as tp
from scipy.ndimage import gaussian_filter
from matplotlib.ticker import AutoMinorLocator,Locator

#######		This script overlays the tracked particles located by 1.preprocess_and_track on top of the raw image files. This way, the user can check if the particle tracking worked properly.		#######
######		The overlayed traces are stored in datafolder/datapoint/analysis/traces/

datafolder = 'sampledata_notfulldataset/dataset/'			#Location of the datafolder
#resultfolder = 'sampledata_notfulldataset/results/'			#Location of the result folder where the overlayed traces should end up.

nr_of_frames = 1634											#Number of frames in the video.
framespersecond = 40 										#Frames per second used in the datacollection (used to calculate velocities)
frames_trace=(0,1634)										#Select which frames where analysed (by default all of them)
channelwidth_um = 895 										#width of the channel in micrometers. This is used to convert distances in pixels to distances in um.

datapoints = os.listdir(datafolder)
for dataindex,datapoint in enumerate(datapoints):

	#Load in the frames
	try:
		frames = pims.open(datafolder+datapoint+'/rawdata/*.bmp')
		print("\n---- Dataset loaded succesfully, Continuing with analysis ----")
	except:
		print("\n----  Dataset not found, exiting... ---- ")
		quit()

	#Load in the tracking parameters
	trackingparameters=np.loadtxt(datafolder+datapoint+"/analysis/trackingparameters.csv",skiprows=1,delimiter=',')
	framestep = int(trackingparameters[0])
	channelleft = int(trackingparameters[1])
	channelright = int(trackingparameters[2])

	channelwidth_px = channelright-channelleft
	selectionleft = channelleft
	selectionright = channelright

	um_per_pixel = channelwidth_um/(channelright-channelleft)
	frames = pims.process.crop(frames,((0,25),(selectionleft,len(frames[0][0,:])-selectionright)))

	tracepath = datafolder+datapoint+"/analysis/traces/"
	if not os.path.exists(tracepath):
		os.makedirs(tracepath)

	splitdatapoints = datapoint.split("-")
	flowrate = splitdatapoints[0][:-5] #uL/min
	
	#Load in the trajectories
	t = pd.read_pickle(datafolder+datapoint+"/analysis/trajectories")
	t_particles_frame0 = t.loc[t['frame']==frames_trace[0]] #Find the particles present at frame 0
	particlelist_frame0 = t_particles_frame0.get('particle').to_numpy()

	#Plot a figure to store the overlayed image.
	figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (12,10))
	plt.subplots_adjust(left= 0.05,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.07, right = 0.93)
	ax.set_aspect(1)

	print(framestep)
	print(nr_of_frames)
	print(tracepath)

	#Loop over the frames, overlay the trajectories on the images, and save the images.
	tracestepper=0
	for framenr in range(frames_trace[0],frames_trace[1],framestep):
		framestepper = framenr

		ax.imshow(frames[framenr],cmap="Greys_r")

		if not framenr == 0:
			tracestepper = round(framenr/framestep)

		#Find all the particles in this frame, that also existed in frame 0. These are the particles we are interested in in the uFFS method.
		t_particles_infr0 = t.loc[t['particle'].isin(particlelist_frame0)]
		t_particles_infr_infr0 = t_particles_infr0.loc[t_particles_infr0['frame']==tracestepper]
		t_particles_framelst_infr0 = t_particles_infr0.loc[t_particles_infr0['frame']<=tracestepper]

		xpositions = t_particles_infr_infr0.get('x').to_numpy()
		ypositions = t_particles_infr_infr0.get('y').to_numpy()

		#overlay the tracked particles using a scatter
		ax.scatter(xpositions,ypositions,s=80,facecolors='none',edgecolors='firebrick')
		
		#Save the figure
		plt.savefig(tracepath+str(framenr)+".png")

		ax.cla()

	plt.close()
