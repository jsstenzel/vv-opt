
"""
Full matrix experiment model
"""
def fp_likelihood_fn(theta, d):
	#define interest params:
	#gain = theta[0]
	#rn = theta[1]
	#dc = theta[2]
    #define design vars:
    
    y1 = read_noise_exp([theta[0],theta[1]], d)
    y2 = dark_current_exp([theta[0],theta[1],theta[2]], d)
    y3 = gain_exp([theta[0],theta[1],theta[2]], d)
    
    return [y1, y2, y3]


"""
theta: [0] gain [1] read noise [2] dark current
d: [0] number exposures
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def read_noise_exp(_theta, _d):
	#define interest params:
	gain = _theta[0]
	rn = _theta[1]
	#define parameters:
	nx = 
	ny = 
	sigma_stray = 
	sigma_dc = 
	t = .1 #100 ms exposure, best for both noise and 
	
	sigma_si = gain * math.sqrt(rn**2 + (sigma_dc*t)**2 + sigma_stray**2)
	
	n = nx * ny * _d[0]
	err = (2/n) * sigma_si**4
    random = 
	
	y = sigma_si + random
	return y
	
	
"""
theta: [0] gain [1] read noise [2] dark current
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def dark_current_exp(_theta, _d):
	#define interest params:
	gain = _theta[0]
	rn = _theta[1]
	dc = _theta[2]
	#define parameters:

	

	return y
	
	
"""
theta: [0] gain [1] read noise [2] dark current
d: 
x parameters: sigma_stray, sigma_dc, nx, ny
"""
def gain_exp(_theta, _d):
	#define interest params:
	gain = _theta[0]
	rn = _theta[1]
	dc = _theta[2]
	#define parameters:

	

	return y
	