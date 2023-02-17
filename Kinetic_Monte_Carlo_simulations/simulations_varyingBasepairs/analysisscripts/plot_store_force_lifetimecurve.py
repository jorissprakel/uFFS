from __future__ import division, unicode_literals, print_function

import numpy as np
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator,Locator
from matplotlib import cm
from matplotlib.colors import rgb2hex
import pickle as pk
import os
import scipy.optimize as fit

#fit function
def exponential_logfit_dx(F,dx,tau0):
	F = np.exp(F)
	F = np.multiply(F,10**(-12))
	taus=tau0*np.exp((-1)*F*dx/(kbT))
	taus = np.log(taus)
	return taus

inferno = cm.get_cmap('inferno', 255)

datafolder = "../series/"

colors=[
'steelblue',
'firebrick',
'forestgreen'
]

legends=[
'sim. 15BP',
'sim. 25BP',
'sim. 30BP'
]

#minimumtime = 1e-5 #The lowest time taken into account for the fits

nr_repeats = 20
fitting = True
guesses=[1.7e-9,1e10]
global kbT

figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (6.4,5.4))
plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.20, right = 0.95)

ax.tick_params(axis='y',**{'which':'major','width':2.2,'length':6,'direction':'in','left':True,'right':True})
ax.tick_params(axis='y',**{'which':'minor','width':1.5,'length':3,'direction':'in','left':True,'right':True})	
ax.tick_params(axis='x',**{'which':'major','width':1.9,'length':6.5,'direction':'in','bottom':True,'top':True})	
ax.tick_params(axis='x',**{'which':'minor','width':1.3,'length':3.5,'direction':'in','bottom':True,'top':True})
[x.set_linewidth(1.6) for x in ax.spines.values()]

tau0values_determined = []
dxvalues_determined = []
nr_basepairs_imposed = []

seriesnames=os.listdir(datafolder)

for seriesindex,seriesname in enumerate(seriesnames):
	bpvalues = os.listdir(datafolder+seriesname+"/simulation/")

	for bpindex,bpname in enumerate(bpvalues):
		print(bpname)
		bps = float(bpname.split("_")[1])

		if not os.path.exists(datafolder+seriesname+"/analysis/"+bpname):
			os.makedirs(datafolder+seriesname+"/analysis/"+bpname)

		stresses = []
		shearforces_pN = []

		lifetime_averages = []
		lifetime_stdevs = []

		stressvaluenames = os.listdir(datafolder+seriesname+"/simulation/"+bpname+"/")
		stressvalues = []
		for i in stressvaluenames:
			stressvalues.append(float(i.split("_")[1]))

		for stressindex,stress in enumerate(stressvalues):
			stressvaluename = stressvaluenames[stressindex]
		
			lifetimes_atstress=[]
			detachmentsfound = 0
			for repeat in range(nr_repeats):


				logfilename=datafolder+seriesname+"/simulation/"+bpname+"/"+stressvaluename+"/rep_"+str(repeat)+"/log.csv"

				#Read in the simulation parameters from the input file used for the simulations
				par={}
				inputfile = open(datafolder+seriesname+"/simulation/"+bpname+"/"+stressvaluename+"/rep_"+str(repeat)+"/input.txt",'r')
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
				shearforcevalue_pN = shearforcevalue*10**(12)

				if os.path.exists(logfilename):

					logfile = np.loadtxt(logfilename,delimiter=",",skiprows=1)

					if logfile.ndim > 1:
						time = logfile[:,1]
						detached = logfile[0,7]

						if detached:	#If the current simulation detached during the course of the simulation, append the detachment time
							lifetimes_atstress.append(time[-1])
							detachmentsfound+=1
			
			if detachmentsfound>0:
				stresses.append(stress)
				shearforces_pN.append(shearforcevalue_pN)
				lifetime_averages.append(np.average(lifetimes_atstress))
				lifetime_stdevs.append(np.std(lifetimes_atstress))

	ax.scatter(shearforces_pN,lifetime_averages,color=colors[seriesindex],s=20,label = legends[seriesindex])
	
	curve_for_storage=np.vstack((shearforces_pN,lifetime_averages,lifetime_stdevs)).T
	np.savetxt(datafolder+seriesname+"/analysis/"+bpname+"/force_lifetimecurve.csv",curve_for_storage,delimiter = ",", header = "F_shear[pN],tau_average[s],tau_stdevs[s]")

	if fitting:
		if len(lifetime_averages)>3:
			
			logforces = np.log(shearforces_pN)
			loglifetimes = np.log(lifetime_averages)
			fitparams,covariance = fit.curve_fit(exponential_logfit_dx,logforces,loglifetimes,p0=guesses,bounds=((0,10**3),(10*10**(-9),10**17)))

			dx = fitparams[0]
			tau0=fitparams[1]

			dxvalues_determined.append(dx)
			tau0values_determined.append(tau0)
			nr_basepairs_imposed.append(bps)

			fitforces = np.linspace(0,max(shearforces_pN)*1.2,200)
			logfitforces = np.log(fitforces)

			answer = exponential_logfit_dx(logfitforces,dx,tau0)
			answer = np.exp(answer)

			ax.plot(fitforces,answer,c=colors[seriesindex])

for tick in ax.xaxis.get_major_ticks():
 	tick.label.set_fontsize(17) 
	
for tick in ax.yaxis.get_major_ticks():
 	tick.label.set_fontsize(17) 
	
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))

ax.set_yscale('log')
ax.set_ylim(top = 1e3, bottom=1e-3)
ax.set_xlim(left = -0.1, right = 130)

ax.set_xlabel(r'$F_{shear}\ [\mathregular{pN}]$',fontsize=22)
ax.set_ylabel(r'$\tau \ [\mathregular{s}]$',fontsize = 22,rotation="vertical")

legend = ax.legend(frameon=False,labelspacing=0.3,loc="upper right",fontsize = 14,handlelength=0,markerscale=1)
plt.savefig("force_lifetimecurve_log.png",dpi=300)

plt.close()

nr_basepairs_imposed = np.array(nr_basepairs_imposed)		
dxvalues_determined = np.array(dxvalues_determined)*1e9		#convert the activation length from m to nm
tau0values_determined = np.array(tau0values_determined)

storagematrix = np.hstack((nr_basepairs_imposed[:,np.newaxis],dxvalues_determined[:,np.newaxis],tau0values_determined[:,np.newaxis]))
storagefile = "parameters_information.csv"

np.savetxt(storagefile,storagematrix,delimiter=',',header="nr bps(imposed - #),dx(determined - nm),tau0(determined - s)")
