from __future__ import division, unicode_literals, print_function
import numpy as np

#This script computes the shear force and torque on the particle due to the shear flow.
#This is done via interpolation from Goldman et al. (1967) Slow Viscous Motion of a Sphere Parallel to a Plane Wall .2. Couette Flow

def ShearForceAndTorque(h,a):
	#Find alpha_F=Fx/6pi*eta*a*h*shearrate
	#and alpha_T=T/4pi*eta*a^3*shear rate
	#From Goldman et al. by interpolation
	h_a=np.flip([10.0677,3.7622,2.3524,1.5431,1.1276,1.0453,1.005004,1.003202,1]);
	aF=np.flip([1.0587,1.1671,1.278,1.4391,1.616,1.6682,1.6969,1.6982,1.7005]);
	aT=np.flip([.99981,.99711,.99010,.97419,.95374,.94769,.94442,.94427,.94399]);

	if np.divide(h,a)>10.0677:
		alphaF=1
		alphaT=1

	else:
		alphaF=np.interp(h/a,h_a,aF)
		alphaT=np.interp(h/a,h_a,aT)

	return(alphaF,alphaT)