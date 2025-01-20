import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv

sys.path.append('../..')
from problems.problem_definition import *
from approx.gaussian_process import *

#_wave_min = 350
_wave_min = 690
_wave_max = 975
_bandpass = _wave_max - _wave_min

wasach_llamas2200_red = [[690.000000,0.422154],[692.850000,0.437670],[695.700000,0.453143],[698.550000,0.468548],[701.400000,0.483856],[704.250000,0.499042],[707.100000,0.514082],[709.950000,0.528950],[712.800000,0.543623],[715.650000,0.558077],[718.500000,0.572292],[721.350000,0.586246],[724.200000,0.599919],[727.050000,0.613292],[729.900000,0.626346],[732.750000,0.639064],[735.600000,0.651429],[738.450000,0.663427],[741.300000,0.675042],[744.150000,0.686262],[747.000000,0.697074],[749.850000,0.707466],[752.700000,0.717427],[755.550000,0.726949],[758.400000,0.736022],[761.250000,0.744639],[764.100000,0.752793],[766.950000,0.760478],[769.800000,0.767689],[772.650000,0.774422],[775.500000,0.780673],[778.350000,0.786441],[781.200000,0.791724],[784.050000,0.796521],[786.900000,0.800831],[789.750000,0.804656],[792.600000,0.807997],[795.450000,0.810857],[798.300000,0.813237],[801.150000,0.815142],[804.000000,0.816576],[806.850000,0.817543],[809.700000,0.818048],[812.550000,0.818098],[815.400000,0.817699],[818.250000,0.816858],[821.100000,0.815581],[823.950000,0.813877],[826.800000,0.811754],[829.650000,0.809220],[832.500000,0.806284],[835.350000,0.802956],[838.200000,0.799245],[841.050000,0.795161],[843.900000,0.790713],[846.750000,0.785913],[849.600000,0.780770],[852.450000,0.775296],[855.300000,0.769501],[858.150000,0.763397],[861.000000,0.756993],[863.850000,0.750303],[866.700000,0.743337],[869.550000,0.736106],[872.400000,0.728622],[875.250000,0.720896],[878.100000,0.712940],[880.950000,0.704765],[883.800000,0.696382],[886.650000,0.687804],[889.500000,0.679041],[892.350000,0.670104],[895.200000,0.661005],[898.050000,0.651755],[900.900000,0.642364],[903.750000,0.632844],[906.600000,0.623205],[909.450000,0.613457],[912.300000,0.603612],[915.150000,0.593678],[918.000000,0.583667],[920.850000,0.573589],[923.700000,0.563451],[926.550000,0.553266],[929.400000,0.543040],[932.250000,0.532785],[935.100000,0.522507],[937.950000,0.512217],[940.800000,0.501922],[943.650000,0.491631],[946.500000,0.481351],[949.350000,0.471090],[952.200000,0.460857],[955.050000,0.450657],[957.900000,0.440498],[960.750000,0.430387],[963.600000,0.420330],[966.450000,0.410334],[969.300000,0.400404],[972.150000,0.390547],[975.000000,0.380767]]

def prior_mean_vph_red(t):
	f = scipy.interpolate.interp1d(np.array(wasach_llamas2200_red).T[0], np.array(wasach_llamas2200_red).T[1], kind='quadratic')
	return f(t)

vph_red_var = .05
vph_red_ls = (0.05*_bandpass)**2
vph_red_pts = [pt[0] for pt in wasach_llamas2200_red]
prior_gp_vph_red = ["gp_expquad", [vph_red_var, vph_red_ls, vph_red_pts, prior_mean_vph_red]]  
#variance, ls, prior_pts, mean_fn


theta_defs = [
				["vph_thru_red", prior_gp_vph_red, "continuous"],
			]
#need to update with range somehow? These can't be negative

y_defs = 	[	
				"y_vph_red_pts", #expands?
			]

d_defs = 	[
				["d_vph_n_pts", ['uniform', [0,_bandpass*2]], "discrete"],
			]
	

x_defs = 	[
				["wave_min", ["nonrandom", [_wave_min]], "continuous", _wave_min],
				["wave_max", ["nonrandom", [_wave_max]], "continuous", _wave_max],
				["measurement_stddev", ["nonrandom", [.01]], "continuous", .01]
			]


def eta(theta, d, x, err=True):
	vph_thru_red = theta["vph_thru_red"]
	d_vph_n_pts = d["d_vph_n_pts"]
	wave_min = x["wave_min"]
	wave_max = x["wave_max"]
	if err:
		measurement_stddev = x["measurement_stddev"]
	else:
		measurement_stddev = 0
	
	#choose the measurement points
	measurement_pts = np.linspace(wave_min, wave_max, num=d_vph_n_pts)
	
	#make the measurements, assuming that there is one y_i for each measurement point ki
	y_thru = vph_thru_red.eval_gp_cond(measurement_pts, measurement_stddev)
	
	#apply the 0..1 boundaries
	for i,yi in enumerate(y_thru):
		if yi < 0:
			y_thru[i] = 0
		if yi > 1:
			y_thru[i] = 1
		
	return y_thru

def H(theta, x):
	vph_thru_red = theta["vph_thru_red"]
	
	#Take the average of theta over k:
	throughput = vph_thru_red.evaluate()
	avg = np.mean(throughput)
	
	return avg

def Gamma(d, x):
	d_vph_n_pts = theta["d_vph_n_pts"]
	return 100 + d_vph_n_pts #the more measurements, the more expensive, linearly


#_dim_d, _dim_theta, _dim_y, _dim_x, _eta, _H, _G, _x_default, _priors)
gp_test = ProblemDefinition(eta, H, Gamma, theta_defs, y_defs, d_defs, x_defs)


def update_gp_problem(problem, d):
	d_masked = [(math.floor(dd) if problem.d_masks[i]=='discrete' else dd) for i,dd in enumerate(d)]
	d_dict = dict(zip(problem.d_names, d_masked))
	num_vph = d_dict["d_vph_n_pts"]
	
	y_vph = ["y_vph_red_pts"]
	
	mod_y_defs = []
	for yname in problem.y_names:
		if yname in y_vph:
			new_y = [yname+"_"+str(i) for i in range(num_vph)]
			mod_y_defs += new_y
		else:
			mod_y_defs.append(yname)
	
	#problem_d = ProblemDefinition(
	#	problem.eta, 
	#	problem.H, 
	#	problem.G, 
	#	problem.theta_defs, 
	#	mod_y_defs, 
	#	problem.d_defs, 
	#	problem.x_defs
	#)
	
	problem_d = problem
	problem_d.dim_y = len(mod_y_defs)
	problem_d.y_names = mod_y_defs
		
	return problem_d
	
	
if __name__ == "__main__":
	#configure the problem for d
	#print(gp_test)
	d = [20]
	print('\n',d,'\n')
	gp_test_d = update_gp_problem(gp_test, d)
	print(gp_test_d)
	
	#check the prior?
	theta = gp_test_d.prior_rvs(1)
	gp = theta[0].evaluate()
	plt.plot(np.array(wasach_llamas2200_red).T[0], np.array(wasach_llamas2200_red).T[1], 'b')
	plt.plot(gp_test_d.priors[0][1][2], gp, 'g')
	plt.show()
	
	#take a measurement
	y = gp_test_d.eta(theta, d)
	print(y)
	plt.plot(np.array(wasach_llamas2200_red).T[0], np.array(wasach_llamas2200_red).T[1], 'b') #annoying way to plot the mean prior
	plt.plot(gp_test_d.priors[0][1][2], gp, 'g') #annoying way to plot the smapled prior
	plt.plot(np.linspace(_wave_min, _wave_max, num=d[0]), y, 'r') #annoying way to plot the measurement
	plt.show()
	
	#check QoI
	q = gp_test_d.H(theta)
	print(q)
