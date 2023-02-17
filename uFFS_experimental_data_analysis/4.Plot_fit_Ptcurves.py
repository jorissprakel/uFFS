from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize as fit
import os
from matplotlib.ticker import AutoMinorLocator,Locator
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

#This script takes the dissociation times computed by 3.compute_dissociationtimes.py, and constructs P_t curves.
#Note, this script does not work on the test dataset submitted here, because this is not a full dataset and only contains the first 5 frames.
#Next, the P_t curves are fit using a stretched exponential curve to obtain force-lifetime information
#The force-lifetime information is then stored separately.

plasma = cm.get_cmap('plasma', 255)

datafolder = 'sampledata_notfulldataset/dataset/'
resultfolder = 'sampledata_notfulldataset/results/'

framespersecond = 40 				#frames per second used while collecting the data
nr_of_frames = 1634					#Total number of frames in the measurement
fitting = False 						#Fit the Pt curves? yes or no.

flowrateaveragetimewindow = 10 		#s, length of the time window used to determine the average flow rate from the flowrate data.

#Set the lowest and highest flowrate of the measurement series. This is only used to set the limits of the colour map in the Pt plot
minflowrate = 9
maxflowrate = 15

nr_timepoints = 800 				#Number of timebins for obtaining the survivalprobability
scaling = "log"

force_per_flowrate=4.018261 		#pN/uLs-1 Conversion factor from flowrate Q(uL/min) to shear force Fshear(pN). Calculated via equations 3 and 4 in the main text of the uFFS article.
flowratestarttime = 6 #s, rough time estimate before the flowrate was turned on

#Initial guesses for the fitting parameters (tau,alpha,b), and one set of guesses per measurement
guesses=[[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
[28,3,0.7],
]

def stretchedexponential_withoffset(t,tau,alpha,b):
	#Note: This is a stretched exponential fit that aims to separate the time periods before and after the flowrate is turned on. 
	#We let the fit algorithm determine the flowrate turn on time b.
	#Then, we fit the equation as: P(t) = 1 before t=b and P(t) = exp((-t/tau)^alpha) after t=b, as follows:
	f = np.ones(np.shape(t)[0])
	f = np.where(t>b,np.exp(-1*((t-b)/tau)**alpha),f)
	return(f)

#Initiate a matplotlib figure.
figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (8,4.2))
plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.19, right = 0.73)

ax.tick_params(axis='y',**{'which':'major','width':2.2,'length':6,'direction':'in','left':True,'right':True})
ax.tick_params(axis='y',**{'which':'minor','width':1.5,'length':3,'direction':'in','left':True,'right':True})	
ax.tick_params(axis='x',**{'which':'major','width':1.9,'length':6.5,'direction':'in','bottom':True,'top':True})	
ax.tick_params(axis='x',**{'which':'minor','width':1.3,'length':3.5,'direction':'in','bottom':True,'top':True})
[line.set_linewidth(1.6) for line in ax.spines.values()]

datapoints = os.listdir(datafolder)


targetflowrates = []
targetforces = []
forces = []
forces_stdevs = []

lifetimes = []
alphas = []
bvalues = []

fiterrs_tau = [] #fit errors
fiterrs_alpha = []
fiterrs_bvalue = []
for dataindex,datapoint in enumerate(datapoints):

	#Read out the target flowrate from the datapoint name. The datapoint should be named according to the target flow rate, for example: 12.75ulmin, 18.00ulmin, etc.
	splitdatapoints = datapoint.split("-")
	flowrate = float(splitdatapoints[0][:-5]) #uL/min
	
	#Load the analysed data produced by 3.compute_dissociationtimes.py
	analysisresults = np.loadtxt(datafolder+datapoint+"/analysis/traceanalysisresults.csv",delimiter=",",skiprows=1)
	xpositions=analysisresults[:,0]
	ypositions=analysisresults[:,1]
	dissociationtimes = analysisresults[:,2]
	v_average_afterdissociation = analysisresults[:,3]
	MSD_values = analysisresults[:,4]		
	particles = np.ones(np.shape(MSD_values)[0],dtype=bool)

	xpositions=xpositions[particles]
	ypositions=ypositions[particles]
	v_average_afterdissociation = v_average_afterdissociation[particles]
	dissociationtimes = dissociationtimes[particles]

	#Obtain and fit the survival probability of the particles.
	#For every timebin, calculate how many particles remained attach until at least this time. This forms the P_t curve.
	timeseries= np.linspace(start=0,stop=nr_of_frames/framespersecond,num=nr_timepoints)
	nr_particles = len(xpositions)
	P_t = np.array(None)
	for tindex,time in enumerate(timeseries):
		P_t=np.append(P_t,np.sum(dissociationtimes>time))
	P_t = P_t[1:]

	#Determine the average flowrate and the error in the flowrate from the flowrate data stored during the experiment:
	flowrateinfo = np.loadtxt(datafolder+datapoint+"/flowrate_data.txt",delimiter ='\t',usecols=(0,10),skiprows=1)
	times = flowrateinfo[:,0]
	Qvalues = flowrateinfo[:,1] #flowrate values

	#Determine the window for computing the average flowrate. (beyond t=6 s in this case)
	Qaverage_minx_idx = int(min(np.where(times>flowratestarttime)[0]))
	if max(times)>flowrateaveragetimewindow:
		Qaverage_maxx_idx = int(min(np.where(times>flowratestarttime+flowrateaveragetimewindow)[0]))
	else:
		Qaverage_maxx_idx=len(times)-1
		
	#Compute the average flowrate and flowrate standard deviation in this window.
	Qvalues_average=np.average(Qvalues[Qaverage_minx_idx:Qaverage_maxx_idx])
	Qvalues_stdev = np.std(Qvalues[Qaverage_minx_idx:Qaverage_maxx_idx])

	#Compute the force in these windows.
	force = Qvalues_average*force_per_flowrate
	force_stdev=Qvalues_stdev*force_per_flowrate
	
	#Plot the flowrate profile.				
	figure3,ax3=plt.subplots(nrows=1,ncols=1,figsize = (12,3.8),linewidth=2)
	plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.19, right = 0.93)
	ax3.tick_params(axis='y',**{'which':'major','width':2.2,'length':6,'direction':'in','left':True,'right':True})
	ax3.tick_params(axis='y',**{'which':'minor','width':1.5,'length':3,'direction':'in','left':True,'right':True})	
	ax3.tick_params(axis='x',**{'which':'major','width':1.9,'length':6.5,'direction':'in','bottom':True,'top':True})	
	ax3.tick_params(axis='x',**{'which':'minor','width':1.3,'length':3.5,'direction':'in','bottom':True,'top':True})
	ax3.plot(times,Qvalues)
	ax3.vlines(0,ymin=min(Qvalues),ymax=max(Qvalues),colors='k',linewidth=1.6)
	ax3.hlines(y=Qvalues_average,xmin=times[Qaverage_minx_idx],xmax=times[Qaverage_maxx_idx],colors='firebrick')
	ax3.set_xlabel(r'$t\ \mathregular{[s]}$',fontsize=22)
	ax3.set_ylabel(r'$Q\ \mathregular{[uL/min]}$',fontsize = 22,rotation="vertical")
	figure3.savefig(datafolder+datapoint+"/analysis/flowrateresults.png",dpi=600)
	plt.close(figure3)

	#Normalize the P_T curves so that they start at 1
	P_t = np.divide(P_t,P_t[0])

	#Compute the color of the current P_t curve and plot the p_t curve. Curves are colorcoded according to the flowrate.
	colorvalue = plasma(((flowrate-minflowrate)/(maxflowrate-minflowrate)))
	ax.scatter(timeseries,P_t,s=5,color=colorvalue,alpha=1,edgecolors=None,label=str(round(force))+" pN - N = "+str(nr_particles))

	#Next, fit the P_t curves with the stretched exponential function.
	if fitting:
		fitparams,cov = fit.curve_fit(stretchedexponential_withoffset,timeseries,P_t,p0=guesses[0])		
		ax.plot(timeseries,stretchedexponential_withoffset(timeseries,fitparams[0],fitparams[1],fitparams[2]),color=colorvalue,alpha=1)
		lifetime = fitparams[0]
		alpha = fitparams[1]
		bvalue = fitparams[2]

		#Store the obtained fit parameters in arrays.
		targetflowrates.append(flowrate)
		targetforces.append(flowrate*force_per_flowrate)#Target force in pN
		forces.append(force)
		lifetimes.append(lifetime)
		alphas.append(alpha)
		bvalues.append(bvalue)

		forces_stdevs.append(force_stdev)
		fiterrs_tau.append(np.sqrt(cov[0,0]))
		fiterrs_alpha.append(np.sqrt(cov[1,1]))
		fiterrs_bvalue.append(np.sqrt(cov[2,2]))

#Change some more axis parameters.
plt.yscale(scaling)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.set_xlabel(r'$t\ \mathregular{[s]}$',fontsize=22)
ax.set_ylabel(r'$P(t)$',fontsize = 22,rotation="vertical")
ax.set_xlim(left=0,right = max(timeseries))

if scaling == 'linear':
	ax.set_ylim(bottom=0,top=1.5)

elif scaling == 'log':
	ax.set_ylim(bottom=0.05,top=2)
	
for tick in ax.xaxis.get_major_ticks():
 	tick.label.set_fontsize(17)
for tick in ax.yaxis.get_major_ticks():
 	tick.label.set_fontsize(17)

legend = ax.legend(frameon=False,labelspacing=0.3,loc="upper right",bbox_to_anchor=(1.5,1),fontsize = 14,handlelength=0,markerscale=3)

#Save the P-t curve figure in the resultfolder.
filename = resultfolder+"particlesurvival_oneplot_d1_"
filename += scaling
filename += ".svg"
plt.savefig(filename,dpi=600,Transparent=True)
plt.close()

#Save the force-lifetime data
resultmatrix = np.transpose(np.vstack((forces,forces_stdevs,lifetimes,alphas,bvalues,fiterrs_tau,fiterrs_alpha,fiterrs_bvalue,targetforces,targetflowrates)))
np.savetxt(resultfolder+"force-lifetimecurve.csv",resultmatrix,delimiter=",",header="F(pN),F_stdev(pN),tau(s),alpha(-),bvalue(flowrateoffset-s),covariance_tau,covariance_alpha,covariance-bvalue,target force (pN),target flowrate (uL/min),fitfunction_used=stoppedflow_plot_fit_Pt_MSD_simplified")
