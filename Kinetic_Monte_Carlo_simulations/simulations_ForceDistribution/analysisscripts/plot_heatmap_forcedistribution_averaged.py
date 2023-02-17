from __future__ import division, unicode_literals, print_function

import numpy as np
import TotalLinkerForceAndTorque
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator,Locator
import matplotlib.ticker
from matplotlib import cm
from matplotlib.colors import rgb2hex
import matplotlib.colors as colors
import pickle as pk
import os

#This script plots a heatmap of the force distribution, as found in figures 5D and 5E of the uFFS paper. 
#For every simulation snapshot saved during the simulations (/linkers_traj/), the position p0 of the linkers on the channel surface is determined and the extension force F is calculated.
#Using this information, a heatmap is then generated. These heatmaps are averaged across all repeats of the same simulation and normalized.
#

colormap = cm.get_cmap('plasma_r', 255)

heatmap_dimensions = (150,23) 				#nr of datapoints in the xy axes for the heat
heatmap_size = ((-0.7,0.42),(-0.56,0.56))		#boundaries of the heatmap in um (xmin,xmax),(ymin,ymax) 

cutoff_count = 1							#Only count those bins that have at least cutoff_count registered forces. 1 by default (any bins are counted)
a =	0.5*4.34								#um, Set particle radius
T = 293.15 									#K, Temperature
Force_kTperum_to_pN = 4.0454*10**(-3)		

#update heatmap_size according to particle radius, to find the size in um
heatmap_size = np.multiply(heatmap_size,a)		#boundaries of the heatmap in um (xmin,xmax),(ymin,ymax) 

#Redefine the class linker to be able to load in the linker information from the trajectory files.
class linker:
	def __init__(self,p0,p1,connected,edge):
		self.p0=p0
		self.p1=p1
		self.connected=connected
		self.edge=edge

datafolder = "../simulation"
dilutions = os.listdir(datafolder)

counter = 0

#Loop over all the dilution factors.
for dindex, dilution in enumerate(dilutions):
	stresses = os.listdir(datafolder+"/"+dilution+"/")

	#Loop over all the applied stresses.
	for stressindex,stress in enumerate(stresses):
		repeats = os.listdir(datafolder+"/"+dilution+"/"+stress)

		#Initialize an array with linker counts and an array to store the force data of the linkers.
		counts = np.zeros((heatmap_dimensions[0],heatmap_dimensions[1]))
		forcedata = np.zeros((heatmap_dimensions[0],heatmap_dimensions[1]))

		#Define the edges of the bins in um, to create the grid to store the heatmap in.
		distancebinedges_x = np.linspace(start=heatmap_size[0][0],stop=heatmap_size[0][1],num=heatmap_dimensions[0]+1)
		distancebinedges_y = np.linspace(start=heatmap_size[1][0],stop=heatmap_size[1][1],num=heatmap_dimensions[1]+1)

		#Loop over all the repeats of the same simulation.
		for repindex,repeat in enumerate(repeats):
			datapath=datafolder+"/"+dilution+"/"+stress+"/"+repeat+"/"
			analysisfolder = "../analysis/"+dilution+"/"+stress+"/heatmap_forcedistribution/"

			if not os.path.exists(analysisfolder):
				os.makedirs(analysisfolder)

			#Read in the simulation parameters from the input file used for the simulations
			par={}
			inputfile = open(datapath+"input.txt",'r')
			for line in inputfile:
				name,value = line.split(" = ")
				if name in ['attachmentdirection','extralinkers','slipconditions','polymermodel','MaxSteps','print_steps','seed']:
					value = int(value)
				elif name in ['trajectoryfile','parametersfile','logfile']:
					value = str(value)[:-1]
				else:
					value = float(value)
				par[name]=value

			#Find the pickle files containing the trajectories
			picklefiles = os.listdir(datapath+"/linkers_traj/")
			logfile = np.loadtxt(datapath+"/log.csv",delimiter=',',skiprows=1)
			steps=logfile[:,0]
			beadpositions = logfile[:,3]

			#Loop over all the trajectory files
			for pkindex,pickle in enumerate(picklefiles):

				print(pkindex)
				totalnr = len(picklefiles)*len(repeats)*len(stresses)*len(dilutions)	
				percentage = (counter/totalnr)*100

				counter+=1
				linkerpath = picklefiles[pkindex]
				stepvalue = int(linkerpath[:-7])

				#Load the information on all the linkers in this trajectory file.
				with open(datapath+"/linkers_traj/"+linkerpath, 'rb') as pickle_file:
					linkerframe=pk.load(pickle_file)
				with open(datapath+"/edge_traj/"+linkerpath, 'rb') as pickle_file:
					edgeframe=pk.load(pickle_file)

				linkers = linker(p0=linkerframe[:,0:2],p1=linkerframe[:,2:5],connected=linkerframe[:,5].astype('bool'),edge=edgeframe[:])

				#Compute the force and location of every linker
				p0s_connected = linkers.p0[linkers.connected==True]
				p0s_connected = np.multiply(p0s_connected,par['a'])
				linkers_connected =linkerframe[linkerframe[:,5]==1]
				F,T,x,Fx=TotalLinkerForceAndTorque.TotalLinkerForceAndTorque(linkers,par)

				#convert force to pN
				F = np.multiply(F,Force_kTperum_to_pN)
				Forces_connected = F[linkers.connected==True]

				#Now, for every linker, figure out in which heatmap bin it fits, and add its force to this bin:
				for bondindex, p0 in enumerate(p0s_connected):
					if p0s_connected[bondindex,0]>heatmap_size[0][0] and p0s_connected[bondindex,0]<heatmap_size[0][1] and p0s_connected[bondindex,1]>heatmap_size[1][0] and p0s_connected[bondindex,1]<heatmap_size[1][1]:
						xbin = np.where(distancebinedges_x>p0s_connected[bondindex,0])[0][0]-1
						ybin = np.where(distancebinedges_y>p0s_connected[bondindex,1])[0][0]-1

						counts[xbin,ybin]+=1													#counts = nr of bonds within this heatmap bin
						forcedata[xbin,ybin]=forcedata[xbin,ybin]+Forces_connected[bondindex]	#force data of this heatmap bin.

		#After the loop, normalize the heatmap bins by dividing by the number of counts per bin
		forcedata[counts<cutoff_count]=0
		forcedata_normalized=np.divide(forcedata,counts).T
		forcedata_normalized_plot = np.log(forcedata_normalized)

		#Plot the heatmap using matplotlib.
		figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (5.6,3.7))
		plt.subplots_adjust(left= 0.15,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.17, right = 0.91)

		ax.tick_params(axis='y',**{'which':'major','width':1.8,'length':4,'direction':'in','left':True,'right':True})
		ax.tick_params(axis='y',**{'which':'minor','width':1.1,'length':1.5,'direction':'in','left':True,'right':True})	
		ax.tick_params(axis='x',**{'which':'major','width':1.3,'length':4,'direction':'in','bottom':True,'top':True})	
		ax.tick_params(axis='x',**{'which':'minor','width':0.9,'length':1.75,'direction':'in','bottom':True,'top':True})
		[axisline.set_linewidth(1.2) for axisline in ax.spines.values()]

		plot=ax.imshow(forcedata_normalized,norm=colors.LogNorm(vmin=1e-4,vmax=forcedata_normalized[np.isfinite(forcedata_normalized)].max()),cmap=colormap,extent=(heatmap_size[0][0],heatmap_size[0][1],heatmap_size[1][0],heatmap_size[1][1]))
		
		cbar=figure.colorbar(plot,ax=ax)
		
		cbar.set_label(r'$F\ \mathregular{[pN]}$',fontsize=22)
		cbar.ax.tick_params(axis='y',**{'which':'major','width':1.6,'length':6,'direction':'out','left':False,'right':True})
		cbar.ax.tick_params(axis='y',**{'which':'minor','width':0.8,'length':3,'direction':'out','left':False,'right':True})
		[axisline.set_linewidth(1.6) for axisline in cbar.ax.spines.values()]

		cbar.ax.tick_params(labelsize=17)
		
		for tick in ax.xaxis.get_major_ticks():
		 	tick.label.set_fontsize(13) 

		for tick in ax.yaxis.get_major_ticks():
		 	tick.label.set_fontsize(13) 

		ax.set_xlim(left = -0.64*par['a'],right=0.4*par['a'])
		ax.set_ylim(bottom = -0.6*par['a'], top = 0.6*par['a'])

		ax.set_xlabel(r'$x\ \mathregular{[\mu m]}$',fontsize=18)
		ax.set_ylabel(r'$y\ \mathregular{[\mu m]}$',fontsize = 18,rotation="vertical")

		#Store the heatmap (Fig. 5D in the paper.)
		plt.savefig(analysisfolder+"force_heatmap.png",dpi=300)
		plt.close()

		#Next, we plot the cross-section through the heatmap (Fig. 5E in the paper.)
		figure,ax=plt.subplots(nrows=1,ncols=1,figsize = (5.6,3.8))
		plt.subplots_adjust(left= 0.23,wspace=0.38, hspace=0.38,top=0.92,bottom = 0.19, right = 0.93)

		ax.tick_params(axis='y',**{'which':'major','width':2.2,'length':6,'direction':'in','left':True,'right':True})
		ax.tick_params(axis='y',**{'which':'minor','width':1.5,'length':3,'direction':'in','left':True,'right':True})	
		ax.tick_params(axis='x',**{'which':'major','width':1.9,'length':6.5,'direction':'in','bottom':True,'top':True})	
		ax.tick_params(axis='x',**{'which':'minor','width':1.3,'length':3.5,'direction':'in','bottom':True,'top':True})

		lines = [12] #Choose which vertical line to take. In our case, line 12 is the center cross-section.
		distancebinmiddles_x = distancebinedges_x[:-1]+0.5*(distancebinedges_x[1]-distancebinedges_x[0])

		#Calculate a predicted force-distance function based on the assumption that every linker is connected to the nearest position on the surface (See Section S1 in the SI.)
		xvalues_predicted = np.linspace(start = distancebinedges_x[0], stop = distancebinedges_x[-1], num = 200)
		Lvalues_predicted = np.sqrt(par['h']**2+np.square(xvalues_predicted))-par['a']
		extense_predicted = Lvalues_predicted/par['Lmax']
		Forces_predicted = (2/par['lK']) * (np.divide(0.25,np.square(1-extense_predicted))-0.25+extense_predicted-0.8*np.power(extense_predicted,2.15))

		Forces_predicted = np.multiply(Forces_predicted,Force_kTperum_to_pN)
		xvalues_predicted = np.multiply(xvalues_predicted,par['a'])

		#Plot the cross-section, and the force-distance function through the cross-section.
		for line in lines:
			ax.scatter(distancebinmiddles_x,forcedata_normalized[line,:],c=forcedata_normalized_plot[line,:],cmap=colormap)
			ax.plot(xvalues_predicted[xvalues_predicted>-0.54*par['a']],Forces_predicted[xvalues_predicted>-0.54*par['a']],'k',linewidth = 1.5,label="Model")

		loc = matplotlib.ticker.MultipleLocator(base=0.5) # this locator puts ticks at regular intervals
		ax.xaxis.set_major_locator(loc)

		for tick in ax.xaxis.get_major_ticks():
		 	tick.label.set_fontsize(20) 
			
		for tick in ax.yaxis.get_major_ticks():
		 	tick.label.set_fontsize(20) 
			
		ax.xaxis.set_minor_locator(AutoMinorLocator(2))

		ax.set_yscale('log')
		ax.set_ylim(bottom=1e-4, top=1e2)
		ax.set_xlim(left = -0.64*par['a'],right=0.4*par['a'])
		ax.set_xlabel(r'$x\ \mathregular{[\mu m]}$',fontsize=23)
		ax.set_ylabel(r'$F\ \mathregular{[pN]}$',fontsize = 23,rotation="vertical")

		#Store the cross-section
		plt.savefig(analysisfolder+"force_lines_through_heatmap.png",dpi=300)
		plt.close()
