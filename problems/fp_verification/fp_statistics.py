
	
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
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def fp_cost(d, x):
	#define design vars:

	

	return 0
	