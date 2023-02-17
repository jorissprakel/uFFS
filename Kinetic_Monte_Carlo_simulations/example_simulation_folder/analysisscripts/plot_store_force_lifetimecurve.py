from __future__ import division, unicode_literals, print_function

import numpy as np
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator,Locator
from matplotlib import cm
from matplotlib.colors import rgb2hex
import pickle as pk
import os
import scipy.optimize as fit

#This script reads the log.csv files in all the simulation data in the simulation/ folder, calculates the average dissociation times per shear stress, and plots a force-lifetime curve.
#This curve can be fitted with an exponential decay curve if desired.


#exponential fit function.
def exponential_logfit_dx(F,dx,tau0):
	F = np.exp(F)
	F = np.multiply(F,10**(-12))
	taus=tau0*np.exp((-1)*F*dx/(kbT))
	taus = np.log(taus)
	return taus

inferno = cm.get_cmap('inferno', 255)

datafolder = "../simulation"

#Here, fill in the same parameters specified in generate_measurearchitecture.py
flowrates = np.array([24,23.5,23,22.5,22,21.5,21,20.5,20,19.5,19,18.5,18,17.5,17,16.5,16,15.5,15,14.5,14,13.5,13,12.5,12,11.5,11,10.5,10])
basepairs = np.array([15,25,30])
koff0_values = np.array([3.16e-5,3.16e-10,1e-12])
del_values = np.array([0.00175,0.00245,0.0028])
flowrate_to_stress = 6.67440969					
stressvalues = np.round(np.multiply(flowrates,flowrate_to_stress),2)
nr_repeats = 10

fitting = True					#Set to True if you want to fit the simulated force-lifetime curve with an experimental fit.
guesses=[1.7e-9,1e10]			#Initial guesses for the exponential fit ([dx, tau0])
global kbT

plotcolors = ['steelblue','firebrick','forestgreen'] #List with colours used in the plots.

#generate and format a matplotlib figure.
figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (6.4,5.4))
plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.20, right = 0.95)

ax.tick_params(axis='y',**{'which':'major','width':2.2,'length':6,'direction':'in','left':True,'right':True})
ax.tick_params(axis='y',**{'which':'minor','width':1.5,'length':3,'direction':'in','left':True,'right':True})	
ax.tick_params(axis='x',**{'which':'major','width':1.9,'length':6.5,'direction':'in','bottom':True,'top':True})	
ax.tick_params(axis='x',**{'which':'minor','width':1.3,'length':3.5,'direction':'in','bottom':True,'top':True})
[x.set_linewidth(1.6) for x in ax.spines.values()]

tau0values_determined = []
dxvalues_determined = []
basepairvalues_imposed = []

for bpcount, basepairs in enumerate(basepairs):				#Loop over all the DNA constructs simulated
	basepairname = "BPs_"+str(basepairs)

	if not os.path.exists("../analysis/"+basepairname):
		os.makedirs("../analysis/"+basepairname)

	stresses = []
	shearforces_pN = []
	lifetime_averages = []
	lifetime_stdevs = []
	
	for stressindex,stress in enumerate(stressvalues):		#Loop over all the stresses simulated
		stressvaluename = "stress_"+str(stress)
		lifetimes_atstress=[]
		detachmentsfound = 0
		for repeat in range(nr_repeats):					#Loop over all repeats for each simulation.
			
			logfilename=datafolder+"/"+basepairname+"/"+stressvaluename+"/rep_"+str(repeat)+"/log.csv"		
			#Read in the simulation parameters from the input file used for the simulations
			par={}
			inputfile = open(datafolder+"/"+basepairname+"/"+stressvaluename+"/rep_"+str(repeat)+"/input.txt",'r')	#Load in the simulation parameters from the input file.
			for line in inputfile:
				name,value = line.split(" = ")
				if name in ['attachmentdirection','extralinkers','slipconditions','polymermodel','MaxSteps','print_steps','seed']:
					value = int(value)
				elif name in ['trajectoryfile','parametersfile','logfile']:
					value = str(value)[:-1]
				else:
					value = float(value)
				par[name]=value

			shearforcevalue = 32*par['a']**2*stress #in kBT/um 			
			shearforcevalue = shearforcevalue*10**6 #in kBT/m
			kbT = par['kT']
			shearforcevalue = shearforcevalue*par['kT'] #in J/m = N
			shearforcevalue_pN = shearforcevalue*10**(12)					#Convert the shear stress into shear forces [pN]

			if os.path.exists(logfilename):

				logfile = np.loadtxt(logfilename,delimiter=",",skiprows=1)	#Load in the log file
					
				if logfile.ndim > 1:
					time = logfile[:,1]
					detached = logfile[0,7]

					if detached:	#If the current simulation detached during the course of the simulation, append the detachment time
						lifetimes_atstress.append(time[-1])
						detachmentsfound+=1
		
		if detachmentsfound>0:
			stresses.append(stress)
			shearforces_pN.append(shearforcevalue_pN)
			lifetime_averages.append(np.average(lifetimes_atstress))		#Calculate the average detachment time per stress
			lifetime_stdevs.append(np.std(lifetimes_atstress))				#Calculate the standard deviation of the detachment time.

	colorvalue = plotcolors[bpcount]

	ax.semilogy(shearforces_pN,lifetime_averages,'ks-',color=colorvalue,linewidth=0,markersize=8,markeredgewidth=2,markerfacecolor='white',label = str(basepairs)+" BPs")

	curve_for_storage=np.vstack((shearforces_pN,lifetime_averages,lifetime_stdevs)).T
	np.savetxt("../analysis/"+basepairname+"/force_lifetimecurve.csv",curve_for_storage,delimiter = ",", header = "F_shear[pN],tau_average[s],tau_stdevs[s]") #Save the force-lifetime curve.

	#Next, we fit the force-lifetime curve
	if fitting:
		if len(lifetime_averages)>3:	#In case enough force-lifetime datapoints have been found to carry out a proper fit
			
			logforces = np.log(shearforces_pN)
			loglifetimes = np.log(lifetime_averages)

			#Fit using scipy.optimize
			fitparams,covariance = fit.curve_fit(exponential_logfit_dx,logforces,loglifetimes,p0=guesses,bounds=((0,10**3),(10*10**(-9),10**17)))

			#obtain the optimized dx and tau0 values from the fit.
			dx = fitparams[0]
			tau0=fitparams[1]

			dxvalues_determined.append(dx)
			tau0values_determined.append(tau0)
			basepairvalues_imposed.append(basepairs)

			fitforces = np.linspace(0,max(shearforces_pN)*1.2,200)
			logfitforces = np.log(fitforces)

			answer = exponential_logfit_dx(logfitforces,dx,tau0)
			answer = np.exp(answer)

			ax.plot(fitforces,answer,c=colorvalue) 				#Plot the fit

for tick in ax.xaxis.get_major_ticks():							#Adjust tick parameters
 	tick.label.set_fontsize(17) 

for tick in ax.yaxis.get_major_ticks():							#Adjust tick parameters
 	tick.label.set_fontsize(17) 

ax.xaxis.set_minor_locator(AutoMinorLocator(2))					#Set the minor tick interval
ax.yaxis.set_minor_locator(AutoMinorLocator(2))					#Set the minor tick interval

ax.set_yscale('log')											#Logarithmic y axis
ax.set_ylim(bottom=1e-6, top = 1e2)								#Axis limits
ax.set_xlim(left = -0.1, right = 100)							#Axis limits
	
ax.set_xlabel(r'$F_{shear}\ [\mathregular{pN}]$',fontsize=22)	#Axis labels
ax.set_ylabel(r'$\tau \ [\mathregular{s}]$',fontsize = 22,rotation="vertical") #Axis labels

legend = ax.legend(frameon=False,labelspacing=0.3,loc="upper right",fontsize = 14,handlelength=0,markerscale=1)	#legend
plt.savefig("force_lifetimecurve_log.png",dpi=300)				#Store the plot
plt.show()
plt.close()

basepairvalues_imposed = np.array(basepairvalues_imposed)		
dxvalues_determined = np.array(dxvalues_determined)*1e9	#convert the activation length to nm
tau0values_determined = np.array(tau0values_determined)

storagematrix = np.hstack((basepairvalues_imposed[:,np.newaxis],dxvalues_determined[:,np.newaxis],tau0values_determined[:,np.newaxis]))
storagefile = "parameters_information.csv"						#Store the optimized fit parameters in storagefile
np.savetxt(storagefile,storagematrix,delimiter=',',header="basepairs(imposed - # BPs),dx(determined - nm),tau0(determined - s)")