sys.path.append('../..')
from problems.functionals import *
	
"""
theta: [0] gain [1] read noise [2] dark current
"""
def fp_hlva(theta, x):
	#define interest params:
	gain = theta[0]
	rn = theta[1]
	dc = theta[2]
	#define parameters:
	tau = x["tau"] #seconds, exposure time

	QoI = math.sqrt(rn**2 + (tau*dc)**2)
	return QoI
	
"""
theta: [0] gain [1] read noise [2] dark current
"""
def fp_qe_hlva(theta, x):
	#define interest params:
	gain = theta[0]
	rn = theta[1]
	dc = theta[2]
	qe = theta[3]
	#define parameters:
	tau = x["tau"] #seconds, exposure time
	
	#we need to turn qe into a single parameter
	#qe is a Functional
	list_wavelength_qe = qe.to_array(10000) #something suitably big
	qes = [wave_qe[1] for wave_qe in list_wavelength_qe]
	avg_qe = np.mean(qes)

	QoI = math.sqrt(rn**2 + (tau*dc)**2) / avg_qe
	return QoI
	
"""
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def fp_cost(d, x):
	#define design vars:

	

	return 0
	