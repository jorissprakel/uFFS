from __future__ import division, unicode_literals, print_function

import numpy as np
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator,Locator
from matplotlib import cm
from matplotlib.colors import rgb2hex
import os


#This script overlays the simulated force-lifetime plots produced by plot_store_force_lifetimecurve.py on top of the experimental uFFS experimental results, to generate plots such as in Figure 5B and 5C.


############### Parameters ####################

resultfolders = [ #Here, put in the result folders of the experiments (after analysis with the uFFS data analysis scripts)
#Please make sure these strings properly link the user to the result folder of the experimental data, as below:
'../../../../uffsexperiments/Data/22-11-18-uFFSexp_RAT15_5pc_highresistance/results/',
'../../../../uffsexperiments/Data/22-11-11-uFFSexp_RAT25_5pc_highresistance/results/',
'../../../../uffsexperiments/Data/22-11-04-uFFSexp_RAT30_5pc_highresistance/results/'
]

colors=[ #colors to give the experimental curves
'steelblue',
'firebrick',
'forestgreen'
]

legends_exp=[#Legends to give to the experimental curves
'exp. 15bp',
'exp. 25bp',
'exp. 30bp'
]

legends_sim=[#Legends to give to the simulated curves.
'sim. 15bp',
'sim. 25bp',
'sim. 30bp'
]

############### Exponential curve ##############

def exponential_logfit_dx(F,dx,tau0):
	F = np.exp(F)
	F = np.multiply(F,10**(-12))
	taus=tau0*np.exp((-1)*F*dx/(kbT))
	taus = np.log(taus)
	return taus

global kbT
kbT = 4.1e-21 

############### - First, import and plot the EXPERIMENTAL force-lifetime data - ###############

#Generate an empty matplotlib figure and format it.
inferno = cm.get_cmap('inferno', 255)
figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (6.4,5.4))
plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.20, right = 0.95)

ax.tick_params(axis='y',**{'which':'major','width':2.2,'length':6,'direction':'in','left':True,'right':True})
ax.tick_params(axis='y',**{'which':'minor','width':1.5,'length':3,'direction':'in','left':True,'right':True})	
ax.tick_params(axis='x',**{'which':'major','width':1.9,'length':6.5,'direction':'in','bottom':True,'top':True})	
ax.tick_params(axis='x',**{'which':'minor','width':1.3,'length':3.5,'direction':'in','bottom':True,'top':True})
[x.set_linewidth(1.6) for x in ax.spines.values()]

for dataindex,resultfolder in enumerate(resultfolders):
	if dataindex in range(len(resultfolders)):

		#For each experimental dataset, load in the force-lifetime data:
		forcelifetimedata = np.genfromtxt(resultfolder+"force-lifetimecurve.csv",skip_header=1,delimiter=",")
		forces = forcelifetimedata[:,0]
		forces_stdev = forcelifetimedata[:,1]
		lifetimes = forcelifetimedata[:,2]

		#Plot the force-lifetime data.
		ax.scatter(forces,lifetimes,color=colors[dataindex],s=120,edgecolors='k',linewidths=1.75,label=legends_exp[dataindex],alpha=0.8)
		ax.errorbar(forces,lifetimes,xerr=forces_stdev,fmt='None',ecolor=colors[dataindex],elinewidth=1.75,size=10,alpha=1)


############### - Next, import and plot the SIMULATED force-lifetime data - ###############

analysisfolder = "../analysis/"

iterators = os.listdir(analysisfolder)
iteratorvalues = []
for i,it in enumerate(iterators):
	iteratorvalues.append(int(it.split("_")[1]))
iterator_start = it.split("_")[0]
iteratorvalues = np.sort(iteratorvalues)

#import the fit results obtained from plot_store_force_lifetimecruve.py:
fitparameters = np.loadtxt("parameters_information.csv",skiprows=1,delimiter=",")
iterator_parameters = fitparameters[:,0]
dxvalues_determined = fitparameters[:,1]
tau0values_determined = fitparameters[:,2]

#Loop over the iterators that have been found
for iterindex, iterator in enumerate(iteratorvalues):
	iteratorname = iterator_start+"_"+str(iterator)

	colorvalue = colors[iterindex]

	#Read in the force-lifetime curves stored by plot_store_force_lifetimecurve.py
	force_lifetimedata = np.loadtxt(analysisfolder+iteratorname+"/force_lifetimecurve.csv",skiprows=1,delimiter=",")
	shearforces_pN = force_lifetimedata[:,0]
	lifetime_averages = force_lifetimedata[:,1]
	lifetime_stdevs = force_lifetimedata[:,2]

	#Plot these force-lifetime curves in the same figure as the experimental data.
	ax.semilogy(shearforces_pN,lifetime_averages,'ks-',color=colors[iterindex],linewidth=0,markersize=8,markeredgewidth=2,markerfacecolor='white',label = legends_sim[iterindex])

	#Reproduce the exponential fit produced by plot_store_force_lifetimecurve.py	
	dxvalue = dxvalues_determined[iterindex]
	tau0value = tau0values_determined[iterindex]

	fitforces = np.linspace(0,max(shearforces_pN)*1.2,200)
	logfitforces = np.log(fitforces)

	answer = exponential_logfit_dx(logfitforces,dxvalue/1e9,tau0value)
	answer = np.exp(answer)

	#Plot the exponential fit too.
	ax.plot(fitforces,answer,c=colors[iterindex])

for tick in ax.xaxis.get_major_ticks():								#Tick parameters
 	tick.label.set_fontsize(19) 

for tick in ax.yaxis.get_major_ticks():								#Tick parameters
 	tick.label.set_fontsize(19) 

ax.xaxis.set_minor_locator(AutoMinorLocator(5))						#Minor tick locations
ax.yaxis.set_minor_locator(AutoMinorLocator(5))						#Minor tick locations

ax.set_yscale('log')												#Logarithmic y axis
ax.set_ylim(top = 1e3,bottom = 1e-3)								#Axis limits
ax.set_xlim(left = -0.1, right = 130)								#Axis limits

ax.set_xlabel(r'$F_{shear}\ [\mathregular{pN}]$',fontsize=25)		#Axis labels
ax.set_ylabel(r'$\tau \ [\mathregular{s}]$',fontsize = 25,rotation="vertical") #Axis labels

legend = ax.legend(frameon=False,labelspacing=0.3,loc="upper right",fontsize = 18,handlelength=0,markerscale=1)	#Plot the lgend
plt.savefig("force_lifetimecurve_log.png",dpi=300)					#Save the figure
plt.show()															#Show the figure.
plt.close()

