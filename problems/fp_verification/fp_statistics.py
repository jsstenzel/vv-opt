
	
"""
theta: [0] gain [1] read noise [2] dark current
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def fp_hlva(theta):
	#define interest params:
	gain = theta[0]
	rn = theta[1]
	dc = theta[2]
    #define parameters:
    tau = 1800 #seconds, exposure time
	
    QoI = math.sqrt(rn**2 + (tau*dc)**2)
	return QoI
	
	
"""
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def fp_cost(d):
	#define design vars:

	

	return y
	