from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import scipy
import time
import pims
import trackpy as tp
from scipy.ndimage import gaussian_filter

######		This script loads in the raw uFFS data, preprocesses it before particle tracking, and performs the particle tracking 		########

#######		Fill in the data settings		########
datafolder = 'sampledata_notfulldataset/dataset/'	#Location of the dataset. The sample dataset included here contains only 100 frames, the datasets used in the experiment contained 1634
framespersecond = 40 						#frames per second used while collecting the data
nr_of_frames = 1634							#Total number of frames in the experiment
filename = "/rawdata/"						#fill in the start of the name of every frame here, with a * behind it.
filetype = "*.bmp"							#extension of the image file (.bmp,.png, etc.)
analyze_frames=(0,nr_of_frames)				#Tell which frames to track, all of them by default.
channelwidth_microns = 895					#um, width of the field of view

###		Particle tracking parameters for Trackpy (see http://soft-matter.github.io/trackpy/v0.5.0/ for more information)	###
particleradius_um = 2.17 					#um, radius of the particles used
masksize = 9 								#odd integer number. Mask size in pixels set by the localisation algorithm
minmass = 1 								#Minimal intensity for accepting localised spots
maxsize = 6 								#Maximal size

###     Preset parameters for the trajectory linking:
memory = 9 									#Number of frames a particle may be unseen before it stops being tracked
distance = 3 								#Distance in pixels a particle is allowed to move between two frames
framestep = 1 								#Analyze every x frames of the movie 1 by default

###		Preprocessing Parameters 			########
thresholdintensity = 35						#Background threshold used in the pre-processing. This value has to be larger than the most intense picture in the background of the image.
Gaussianblurwidth = 1						#Width of the Gaussian filter used in the pre-processing.

datasets = os.listdir(datafolder)

borders_left = []
borders_right = []
for dataindex,dataset in enumerate(datasets):

	loadpath=datafolder+dataset+filename
	storepath=datafolder+dataset+'/analysis/'

	if not os.path.exists(storepath):
		os.makedirs(storepath)

	print("\n---- Loading designated dataset ----")
	#Try to load in the images:
	try:
		frames = pims.open(loadpath+filetype)

		print("\n---- Dataset loaded succesfully, Continuing with analysis ----")
	except:
		print("\n----  Dataset not found, exiting... ---- ")
		quit()

	frames=frames[analyze_frames[0]:analyze_frames[1]:framestep]   # 

	#Set the size of the window to analyse, using channelleft and channelright. This is the entire width of the image by default (and in all the datasets used in the paper), but can be reduced by the user if so desired.
	channelleft = 0
	channelright = len(frames[0][0,:])-1

	channelwidth_px = channelright-channelleft
	micron_per_pixel = np.divide(channelwidth_microns,channelwidth_px)

	selectionleft = channelleft
	selectionright = channelright
	borders_left.append(selectionleft)
	borders_right.append(len(frames[0][0,:])-selectionright)

#Now that the edges of the field of view have been specified, load in the data
for dataindex,dataset in enumerate(datasets):

	loadpath=datafolder+dataset+filename
	storepath=datafolder+dataset+'/analysis/'
	
	if not os.path.exists(storepath):
		os.makedirs(storepath)

	print("\n---- Loading designated dataset ----")

	try:
		frames = pims.open(loadpath+filetype)
		print("\n---- Dataset loaded succesfully, Continuing with analysis ----")
	except:
		print("\n----  Dataset not found, exiting... ---- ")
		quit()

	frames=frames[analyze_frames[0]:analyze_frames[1]:framestep]
	frames = pims.process.crop(frames,((0,25),(borders_left[dataindex],borders_right[dataindex])))

	#change to float32 datatype:
	frames = np.array(frames)
	frames = frames.astype('float32')
	print(frames)

	#########		First, perform the processing steps required before particle tracking can be done. 			##########

	#Subtract the selected boundaries, in case a smaller window was chosen using channelleft and channelright above (this was not done in the data-analysis used in the paper):
	print("\n---- 1. Cropping the selected boundaries from the images ----")
	time.sleep(1)

	figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (16,10))
	plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.19, right = 0.93)
	midframe = round(len(frames)/2)

	#Save an example of a raw cropped image.
	ax.imshow(frames[0],cmap="Greys_r")
	plt.savefig(storepath+"1.channelscropped.png",dpi=300)
	plt.close()

	#Perform the background threshold processing step
	print("\n ---- 2. Performing background thresholding at threshold intensity "+str(thresholdintensity))

	shape = (len(frames),len(frames[0][:,0]),len(frames[0][0,:]))
	frames_thresholded = np.zeros(shape)
	frames_thresholded=frames_thresholded.astype('float32')
	for frindex,frame in enumerate(frames):
		above_threshold = np.where(frame>thresholdintensity)

		frames_thresholded[frindex]=np.where(frame>thresholdintensity,frame,thresholdintensity)

	#Store a picture of a thresholded image.
	figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (16,10))
	plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.19, right = 0.93)
	ax.imshow(frames_thresholded[0],cmap="Greys_r")
	plt.savefig(storepath+"2.frames_thresholded.png",dpi=300)
	plt.close()

	#Perform the Gaussian blur step.
	print("\n---- 3. Performing a Gaussian blur of width "+str(Gaussianblurwidth)+" px ----")#Make a plot of the loaded raw dataset, showing the first frame, middle frame, and last frame of the series. Plot is stored in the storepath.

	frames_gaussianblurred = np.zeros(shape)
	frames_gaussianblurred = frames_gaussianblurred.astype('float32')
	
	#For each frame in the movie, perform the gaussian blur.
	for fridx,frame in enumerate(frames_thresholded):
		frames_gaussianblurred[fridx]=scipy.ndimage.filters.gaussian_filter(frames_thresholded[fridx],Gaussianblurwidth)
	
	#Store a picture of a gaussian blurred image
	figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (16,10))
	plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.19, right = 0.93)
	midframe = round(len(frames_gaussianblurred)/2)
	ax.imshow(frames_gaussianblurred[0],cmap="Greys_r")
	plt.savefig(storepath+"3.gaussian_blurred.png",dpi=300)
	plt.close()

	#Rescale the blurred image, so that the intensity goes between 0 and 255.
	print("---- 4. Rescaling the blurred images ----")
	arrayminimum = np.amin(frames_gaussianblurred)
	arraymaximum = np.amax(frames_gaussianblurred)
	arrayrange = arraymaximum - arrayminimum

	frames_rescaled = np.zeros(shape)
	frames_rescaled = frames.astype('float32')

	#Perform the rescaling
	for fridx,frame in enumerate(frames_gaussianblurred):
		frame = np.subtract(frame,arrayminimum)
		frame = np.divide(frame,arrayrange)
		frame = np.multiply(frame,255)
		frames_rescaled[fridx]=frame

	#Store a picture of the rescaled image.
	figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (16,10))
	plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.19, right = 0.93)
	ax.imshow(frames_rescaled[0],cmap="Greys_r")
	plt.savefig(storepath+"4.frames_rescaled.png",dpi=300)
	plt.close()

	#########		Now that all the processing has been done, start with the particle tracking. 			##########

	print("---- 5. Localizing features ----")
	print("\nParameters used by the localisation algorithm:\n-mask diameter: "+str(masksize)+" pixels\n-minimal mass: "+str(minmass)+"\n-maxsize: "+str(maxsize)+"\n-radius: "+str(particleradius_um)+" um")

	time.sleep(2)
	#The tp.batch function runs the tracking, using the parameters specified at the top of the script.
	f = tp.batch(frames_rescaled[:], processes=1,diameter=masksize, invert=False,noise_size=0,minmass=minmass,maxsize=maxsize) 

	#Plot a picture of a frame with the localized features.
	print("---- 6. Plotting located features ----")
	figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (16,10))
	plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.19, right = 0.93)
	ax.imshow(frames[0],cmap="Greys_r")
	positions = f.loc[f['frame']==0]
	xpositions = positions.get('x').to_numpy()
	ypositions = positions.get('y').to_numpy()
	ax.scatter(xpositions,ypositions, c="r",s = 0.1)
	plt.savefig(storepath+"6.features_localized.png",dpi=300)
	plt.close()

	#########		After tracking the particles, link the tracked particles into trajectories, using the tp.link function. 		#########
	t= tp.link(f,distance,memory=memory)

	#Store the linked trajectories.
	t.to_pickle(path=storepath+"trajectories")

	#Finally, store the tracking parameters used in the particle tracking algorithm for future reference.
	trackingparameters = np.array([framestep,channelleft,channelright,particleradius_um,masksize,minmass,maxsize,memory,distance,Gaussianblurwidth])
	np.savetxt(storepath+"trackingparameters.csv",trackingparameters,delimiter='',header = 'framestep (analysed every x frames),channelleft(px),channelright(px),particleradius (um),masksize (px), minmass (I), maxsize(px), memory (frames),distance(px),Gaussianblur width (px)')
