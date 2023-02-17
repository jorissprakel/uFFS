
import random

def generate_input(parameterlib,outputfilename):

	print(outputfilename)
	outputfile = open(outputfilename,'w')
	for elem in parameterlib:
		outputfile.write(elem+" = "+str(parameterlib[elem])+"\n")

	outputfile.close()
	#print("--- oxDNA input file written to: "+outputfilename+" ---")
