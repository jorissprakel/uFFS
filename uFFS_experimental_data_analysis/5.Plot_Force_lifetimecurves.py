from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import scipy.optimize as fit
from matplotlib.ticker import AutoMinorLocator,Locator
import os

#This script takes the force-lifetime curves computed by 4.Plot_fit_Ptcurves.py of multiple analysed datasets
#The force-lifetime curves are then plotted in one figure to plot figures (3A and 3B) from the uFFS paper.

resultfolders = [#Insert the series of dataset names here.
'../Data/22-11-10-uFFSexp_RAT15_5pc_highresistance/results/',
'../Data/22-11-18-uFFSexp_RAT15_5pc_highresistance/results/',
'../Data/22-11-11-uFFSexp_RAT25_5pc_highresistance/results/',
'../Data/22-11-04-uFFSexp_RAT30_5pc_highresistance/results/',
]

colors=[#Plot colors for every dataset
'steelblue',
'navy',
'firebrick',
'forestgreen'
]

legends=[#Legends
'15BP-1',
'15BP-2',
'25BP',
'30BP',
]

#Build a matplotlib figure
figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (5.4,4.8))
plt.subplots_adjust(left= 0.18,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.20, right = 0.95)

ax.tick_params(axis='y',**{'which':'major','width':2.2,'length':6,'direction':'in','left':True,'right':True})
ax.tick_params(axis='y',**{'which':'minor','width':1.5,'length':3,'direction':'in','left':True,'right':True})	
ax.tick_params(axis='x',**{'which':'major','width':1.9,'length':6.5,'direction':'in','bottom':True,'top':True})	
ax.tick_params(axis='x',**{'which':'minor','width':1.3,'length':3.5,'direction':'in','bottom':True,'top':True})
[x.set_linewidth(1.6) for x in ax.spines.values()]

#Loop over the datasets
for dataindex,resultfolder in enumerate(resultfolders):

	if dataindex in range(len(resultfolders)):

		#Load the force-lifetime curve
		forcelifetimedata = np.genfromtxt(resultfolder+"force-lifetimecurve.csv",skip_header=1,delimiter=",")
		
		print(resultfolder+"force-lifetimecurve.csv")
		forces = forcelifetimedata[:,0]
		forces_stdev = forcelifetimedata[:,1]
		lifetimes = forcelifetimedata[:,2]

		#Plot the force-lifetime data
		ax.scatter(forces,lifetimes,color=colors[dataindex],s=120,edgecolors='k',linewidths=1.75,label=legends[dataindex],alpha=0.8)
		ax.errorbar(forces,lifetimes,xerr=forces_stdev,fmt='None',ecolor=colors[dataindex],elinewidth=2,size=10,alpha=1)

#Change some further plot parameters.
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_yscale("log")

for tick in ax.xaxis.get_major_ticks():
 	tick.label.set_fontsize(15) 
for tick in ax.yaxis.get_major_ticks():
 	tick.label.set_fontsize(15) 

plt.legend(frameon=False,labelspacing=0.3,handlelength = 0,loc="upper right",fontsize = 14)
ax.set_xlabel(r'$F_{shear}\ \mathregular{[pN]}$',fontsize=20)
ax.set_ylabel(r'$\tau\ \mathregular{[s]}$',fontsize = 20,rotation="vertical")
ax.set_xlim(left=0,right=120)
ax.set_ylim(bottom=1e-2,top=1e3)

#Save the force-lifetime curves
plt.savefig("force_lifetimecurve.png",dpi=600)
plt.close()