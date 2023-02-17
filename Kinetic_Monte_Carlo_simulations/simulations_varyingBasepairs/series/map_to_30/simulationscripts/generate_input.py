
import random

#This file takes the library of parameters provided by generate_measurearchitecture.py and writes an input file for the simulation.
def generate_input(parameterlib,outputfilename):

	print(outputfilename)
	outputfile = open(outputfilename,'w')
	for elem in parameterlib:
		outputfile.write(elem+" = "+str(parameterlib[elem])+"\n")

	outputfile.close()
