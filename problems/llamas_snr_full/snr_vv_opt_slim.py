import sys
import os
import scipy.stats
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv
import dill

sys.path.append('../..')
#llamas
from problems.problem_definition import *
from problems.llamas_snr_full.snr_cost_model import *
#analysis
from obed.obed_gbi import *
#from obed.pdf_estimation import *
from inference.bn_modeling import *
from inference.bn_evaluation import *
from uq.uncertainty_propagation import *
from opt.nsga import *

#This version of the optimization code does not include or pass around the problem object
#this is intended to cut down on pickling time & memory incurred by multiprocessing

def construct_slim_llamas_snr_problem():
	empty_prior = ["uniform", [0,1]]
	theta_defs = [                             #mean, variance
						["gain_red", empty_prior, "continuous"],
						["gain_gre", empty_prior, "continuous"],
						["gain_blu", empty_prior, "continuous"],
						["rn_red", empty_prior, "continuous"],
						["rn_gre", empty_prior, "continuous"],
						["rn_blu", empty_prior, "continuous"],
						["dc_red", empty_prior, "continuous"],
						["dc_gre", empty_prior, "continuous"],
						["dc_blu", empty_prior, "continuous"],
						["qe_red_t", empty_prior, "continuous"],
						["qe_gre_t", empty_prior, "continuous"],
						["qe_blu_t", empty_prior, "continuous"],
						["vph_red_t", empty_prior, "continuous"],
						["vph_gre_t", empty_prior, "continuous"],
						["vph_blu_t", empty_prior, "continuous"],
						["sl_t", empty_prior, "continuous"],
						["bg_t", empty_prior, "continuous"],
						["coll_t", empty_prior, "continuous"],
						["red_l1_t", empty_prior, "continuous"],
						["red_l2_t", empty_prior, "continuous"],
						["red_l3_t", empty_prior, "continuous"],
						["red_l4_t", empty_prior, "continuous"],
						["red_l5_t", empty_prior, "continuous"],
						["red_l6_t", empty_prior, "continuous"],
						["red_l7_t", empty_prior, "continuous"],
						["gre_l1_t", empty_prior, "continuous"],
						["gre_l2_t", empty_prior, "continuous"],
						["gre_l3_t", empty_prior, "continuous"],
						["gre_l4_t", empty_prior, "continuous"],
						["gre_l5_t", empty_prior, "continuous"],
						["gre_l6_t", empty_prior, "continuous"],
						["gre_l7_t", empty_prior, "continuous"],
						["blu_l1_t", empty_prior, "continuous"],
						["blu_l2_t", empty_prior, "continuous"],
						["blu_l3_t", empty_prior, "continuous"],
						["blu_l4_t", empty_prior, "continuous"],
						["blu_l5_t", empty_prior, "continuous"],
						["blu_l6_t", empty_prior, "continuous"],
						["blu_l7_t", empty_prior, "continuous"],
						["blu_l8_t", empty_prior, "continuous"],
						["fiber_frd", empty_prior, "continuous"]
					]
					
	y_defs = [	
				"y_gain_red", 
				"y_gain_gre", 
				"y_gain_blu", 
				"y_rn_red", 
				"y_rn_gre", 
				"y_rn_blu", 
				"y_dc_red", 
				"y_dc_gre", 
				"y_dc_blu", 
				"y_qe_red_t0", 
				"y_qe_red_t1", 
				"y_qe_red_t2", 
				"y_qe_red_t3", 
				"y_qe_red_t4", 
				"y_qe_gre_t0", 
				"y_qe_gre_t1", 
				"y_qe_gre_t2", 
				"y_qe_gre_t3", 
				"y_qe_gre_t4", 
				"y_qe_blu_t0", 
				"y_qe_blu_t1", 
				"y_qe_blu_t2", 
				"y_qe_blu_t3", 
				"y_qe_blu_t4", 
				"y_vph_red_p0", 
				"y_vph_red_p1", 
				"y_vph_red_p2", 
				"y_vph_gre_p0", 
				"y_vph_gre_p1", 
				"y_vph_gre_p2",
				"y_vph_blu_p0", 
				"y_vph_blu_p1", 
				"y_vph_blu_p2", 
				"y_sl_t0", 
				"y_sl_t1", 
				"y_sl_t2", 
				"y_sl_t3", 
				"y_bg_t0", 
				"y_bg_t1", 
				"y_bg_t2", 
				"y_bg_t3", 
				"y_coll_t",
				"y_red_cam_t",
				"y_gre_cam_t",
				"y_blu_cam_t",
				"y_frd"
			]
	
	d_defs_discrete = [
					["t_gain", ['uniform', [0, 600]], "discrete"], #gain
					["I_gain", ['uniform', [0, 100]], "discrete"],    #gain
					["n_meas_rn", ['uniform', [0, 50]], "discrete"],  #rn
					["d_num", ['uniform', [2, 25]], "discrete"],      #dc
					["d_max", ['uniform', [.1, 12000]], "discretized100"], #dc
					["d_pow", ['uniform', [0.01,3]], "discretized100"],      #dc
					
					["n_qe", ['uniform', [0, 100]], "discrete"],   #qe
					["t_qe", ['uniform', [0, 600]], "discrete"],#qe
					
					["d_vph_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
					["d_dichroic_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
					["d_coll_n_pts", ['uniform', [0,_bandpass*10]], "discrete"],
					["d_redcam_n_pts", ['uniform', [0,(_wave_max-_wave_redgreen)*10]], "discrete"], #i guess??
					["d_greencam_n_pts", ['uniform', [0,(_wave_redgreen-_wave_greenblue)*10]], "discrete"], #i guess??
					["d_bluecam_n_pts", ['uniform', [0,(_wave_greenblue-_wave_min)*10]], "discrete"], #i guess??
					["d_frd_n_meas", ['uniform', [0,2400]], "discrete"],
				]
	
	x_defs = [
				#rn exp
				["t_rn", [], "continuous", .1], #100 ms exposure
				#dc exp
				["t_0", [], "continuous", 0.1], #100ms baseline exposure assumed
				#cost
				["t_gain_setup", [], "continuous", 1200], #rough estimate based on experience
				["t_gain_buffer", [], "continuous", 5], #rough estimate based on experience
				["t_rn_buffer", [], "continuous", 5], #rough estimate based on experience
				["t_dc_buffer", [], "continuous", 5], #rough estimate based on experience
				["fp_testbed_setup", [], "continuous", 1800], #rough estimate based on experience
				["t_qe_setup", [], "continuous", 1200], #WAG just assuming all setups are 1hr
				["t_qe_buffer", [], "continuous", 60], #WAG - need to adjust photodiode or sensor or something, so say about a minute
				["t_vph_setup", [], "continuous", 1200], #WAG just assuming all setups are 1hr
				["t_vph_per_pt", [], "continuous", 5*60], #WAG - individual measurements -5 min each with buffer
				["t_dichroic_setup", [], "continuous", 1200], #WAG just assuming all setups are 1hr
				["t_dichroic_per_pt", [], "continuous", 0.01], #WAG - this captures tradeoff between spectral resolution and integration time
				["t_coll_setup", [], "continuous", 1200], #WAG just assuming all setups are 1hr
				["t_coll_per_pt", [], "continuous", 0.05], #WAG - this captures tradeoff between spectral resolution and integration time
				["t_camera_test_setup", [], "continuous", 1800], #WAG this setup is probably more complicated
				["t_camera_per_pt", [], "continuous", 5*60], #WAG - individual measurements - 5 min each with buffer
				["t_frd_setup", [], "continuous", 1200], #WAG just assuming all setups are 1hr
				["t_frd_test", [], "continuous", 600], #WAG i think this test was probably pretty fiddly, since we only tested 10
				#["C_engineer", [], "continuous", 0.00694444444] #WAG $/s, from $25/hr
				["C_engineer", [], "continuous", 1] #just count time
	]
	
	eta = None
	H = None
	Gamma = snr_cost
	llamas_snr = ProblemDefinition(eta, H, Gamma, theta_defs, y_defs, d_defs_discrete, x_defs)
	print("llamas_snr_full SLIM problem constructed.",flush=True)
	return llamas_snr
	
	
def vv_OPT(problem, gmm_file, ysamples_file, design_pts, epsilon, util_err, do_hrs, do_min, threads, popSize, nMC, displayFreq=10):
	#Load the GMM and presampled y from file
	print("Loading GMM and presamples...",flush=True)
	gmm20 = bn_load_gmm(gmm_file)
	presampled_ylist = bn_load_y(problem, ysamples_file, doPrint=False, doDiagnostic=False)

	nmcp4_1results = [
		[10,0,0,7,240.098,0.8173,6,73,0,5579,1451,0, 0,17,3],
		[10,0,0,7,360.097,2.6412,6,37,0,5110,1529,1, 0,19,3],
		[10,0,0,7,360.097,2.6412,6,65,0,5110,1529,1, 0,19,3],
		[11,0,1,6,1080.09,2.7608,6,14,4,5105,1499,0,13,15,7],
		[10,39,0,7,360.097,2.6412,6,65,2,5139,1529,1, 0,19,11],
		[10,2,0,3,240.098,0.6379,6,60,0,5015,1446,0, 0,127,3],
		[10,1,0,2,360.097,1.9236,6,35,11,4941,1906,0, 0,141,3],
		[10,0,0,7,240.098,1.9236,6,66,12,5325,1440,0, 0,124,16],
		[10,2,0,3,480.096,0.6379,6,69,1,5015,1936,2, 0,194,3],
		[13,38,0,12,840.093,2.8206,11,86,8,5041,1459,1,34,140,3],
		[10,41,0,12,1680.09,2.8206,12,86,8,4971,1459,16,60,140,0],
		[10,0,2,12,240.098,2.9402,14,90,12,5325,1440,2,117,124,16],
		[11,39,3,12,2040.08,2.7309,14,189,74,5042,2956,109,135,258,5],
		[46,33,2,12,1560.09,1.9236,14,104,225,4986,2586,34,91,228,4],
		[46,39,2,12,1920.08,2.8206,14,190,225,5039,2417,34,133,253,5],
		[25,39,2,12,1080.09,2.9701,14,170,224,4965,1741,4,643,228,11],
		[25,39,0,12,1080.09,2.7309,14,155,350,4965,1396,4,643,228,11],
		[13,2,2,12,1680.09,2.8505,14,190,550,4940,2408,77,747,251,7],
		[16,35,2,12,1080.09,2.9402,14,170,469,4934,1743,691,870,270,3],
		[26,28,2,14,1560.09,2.8804,14,192,178,4967,2599,1653,805,321,72],
		[27,34,0,12,1560.09,2.7608,14,188,558,4994,2417,1726,753,253,6],
		[26,28,2,14,1560.09,2.8804,14,192,500,4967,2599,1792,805,321,6],
		[13,28,0,12,1560.09,2.9103,14,190,543,4964,2599,1792,753,322,6],
		[10,40,0,12,1920.08,2.0432,14,190,543,5023,2337,1690,843,257,2201],
	]

	nmcp4_2results = [
		[10,0,0,7,240.098,2.4917,5,64,0,5361,1279,0,0,17,3],
		[10,0,0,7,240.098,0.8173,5,73,0,5361,1451,0,0,17,3],
		[10,0,0,7,360.097,2.6412,6,35,0,5110,1529,0,0,19,3],
		[10,36,0,7,360.097,2.6412,6,72,0,5368,1528,0,0,19,10],
		[10,0,0,7,360.097,2.6412,5,42,5,5110,1529,1,17,19,3],
		[10,0,0,2,240.098,2.6412,6,57,0,5072,1529,0,0,121,3],
		[13,28,0,12,1080.09,2.8206,11,86,12,4965,1405,1,60,140,2],
		[13,38,0,12,840.093,2.8206,13,155,14,5107,1459,6,60,140,3],
		[10,0,0,2,240.098,2.6412,5,58,3,5110,1530,0,1,326,0],
		[13,4,0,13,240.098,1.9535,13,63,178,4961,1396,0,0,335,18],
		[9,4,8,13,1080.09,2.6711,11,85,172,5042,1525,0,34,318,1],
		[46,33,0,12,1560.09,2.0432,6,104,225,5105,1479,33,14,248,3],
		[45,32,1,12,1560.09,2.6412,14,104,225,4986,2586,34,91,228,4],
		[25,39,0,12,1080.09,2.7309,14,155,221,4965,1396,4,643,228,11],
		[27,39,2,12,1680.09,2.9701,14,170,224,4965,2418,4,643,252,11],
		[26,34,0,11,1560.09,2.7608,14,188,169,4996,2417,1726,753,253,6],
		[12,34,0,12,1560.09,2.7608,11,184,558,4994,2403,1726,753,253,6],
		[25,32,2,6,1560.09,2.8505,14,193,558,4967,2528,1889,803,320,1],
		[26,40,0,14,1920.08,2.0432,14,190,543,5023,2337,1796,843,257,2201],
		[0.1, #t_gain                            0-length exposure for gain exp
                                                0, #I_gain                              no measurements for gain exp
                                                0, #n_meas_rn                   no measurements for rn exp
                                                0, #d_num                               no measurements for dc exp
                                                1, #d_max                               dc exp filer
                                                0.01, #d_pow                            dc exp filler
                                                0, #n_qe                                no measurements for qe exp
                                                0.1, #t_qe                              0-length exposure for qe exp
                                                0, #d_vph_n_pts                 no measurements for vph exp
                                                0, #d_dichroic_n_pts    no measurements for dichroic exp
                                                0, #d_coll_n_pts                no measurements for collimator exp
                                                0, #d_redcam_n_pts              no measurements for camera exp
                                                0, #d_greencam_n_pts            no measurements for camera exp
                                                0, #d_bluecam_n_pts             no measurements for camera exp
                                                0],  #d_frd_n_meas                no measurements for FRD exp
	[600, #t_gain
                                                100, #I_gain
                                                50, #n_meas_rn
                                                25, #d_num
                                                12000, #d_max
                                                3, #d_pow
                                                100, #n_qe
                                                600, #t_qe
                                                (975-350)*10, #d_vph_n_pts
                                                (975-350)*10, #d_dichroic_n_pts
                                                (975-350)*10, #d_coll_n_pts
                                                (975-690)*10, #d_redcam_n_pts
                                                (690-480)*10, #d_greencam_n_pts
                                                (480-350)*10, #d_bluecam_n_pts
                                                2400]  #d_frd_n_meas
	]
	#nmcp4_2results_aug = nmcp4_2results + d0_aug + d0_aug + dmax_aug + dmax_aug  

	nmcp4_3results = [
		[0,0,0,0,0.1,0.01,0,0,0,0,0,0,0,0,0],
		[0,0,0,2,0.1,0.01,0,0,0,0,0,0,0,0,0],
		[0,0,2,2,0.1,0.0698,0,0,0,0,47,0,0,0,0],
		[0,0,0,2,0.1,0.0698,0,0,0,1,33,0,0,0,0],
		[22,0,0,2,0.1,0.1595,0,0,0,4982,4,0,0,0,0],
		[22,0,0,2,0.1,0.0698,0,0,0,4982,126,0,0,0,0],
		[0,0,0,2,0.1,0.1595,0,0,0,4942,2374,0,0,0,0],
		[0,0,0,2,120.099,0.0698,0,0,0,4945,2451,0,0,0,0],
		[0,0,0,2,240.098,0.0698,0,56,0,4945,2497,0,0,0,0],
		[24,0,0,2,240.098,0.0997,0,1,0,5075,1274,0,0,2,0],
		[0,6,0,2,840.093,0.0997,0,1,0,5168,1683,0,0,0,0],
		[4,6,0,2,840.093,0.0997,0,1,0,5168,1683,0,0,0,0],
		[0,30,0,11,240.098,0.0698,0,62,0,5385,2497,0,0,0,0],
		[10,0,0,7,360.097,2.6412,6,35,0,5110,1529,0,0,19,3],
		[37,35,0,6,360.097,2.7907,6,35,0,5177,1595,0,0,19,3],
		[22,0,0,6,120.099,1.9236,5,60,0,5188,1268,0,0,120,1],
		[0,0,3,6,120.099,1.9236,5,60,0,5189,1677,0,0,120,1],
		[22,0,0,6,240.098,2.4917,5,64,0,4943,1267,0,0,120,3],
		[45,33,0,12,240.098,2.0432,6,57,0,4965,1479,0,0,218,2],
		[13,38,0,12,1080.09,2.8206,13,155,5,4964,1459,5,32,141,2],
		[22,42,0,11,1200.09,2.7309,14,155,5,4961,1280,0,425,140,11],
		[26,32,0,11,1080.09,1.9834,11,105,159,4952,2421,34,88,227,3],
		[9,39,0,11,1080.09,2.7309,6,162,9,5088,1523,4,636,243,10],
		[45,32,1,11,1560.09,2.6412,12,105,225,4957,2586,34,81,227,4],
		[25,39,0,12,1080.09,2.7309,14,155,221,4965,1396,4,643,228,11],
		[26,31,17,10,1200.09,1.9834,12,188,400,4942,2417,1726,753,253,75],
		[26,29,16,11,1200.09,2.7309,13,188,400,4942,2417,1726,753,253,77],
		[12,34,0,12,1560.09,2.7608,11,184,558,4994,2403,1726,753,253,6]	
	]		

	nmcp4_final = [
		[2,0,0,2,0.1,1.9834,0,0,0,0,0,0,0,0,0],
		[0,0,0,2,0.1,1.8339,0,4,0,4294,0,0,0,0,0],
		[61,0,0,2,0.1,1.804,0,0,0,4297,0,0,0,0,0],
		[64,0,1,2,0.1,1.804,0,5,0,4313,0,0,0,0,0],
		[62,0,0,2,0.1,1.9535,0,0,0,4319,0,0,0,0,0],
		[0,0,0,2,0.1,0.01,0,5,0,4353,0,0,0,0,0],
		[61,0,0,2,0.1,1.9535,0,0,0,4358,0,0,0,0,0],
		[101,0,0,2,0.1,1.9535,0,0,0,4360,0,0,0,0,0],
		[9,0,0,2,0.1,1.9535,0,0,0,4487,0,0,0,0,0],
		[67,0,0,2,0.1,1.804,0,0,0,4519,0,0,0,0,0],
		[1,0,1,2,0.1,1.9834,0,0,0,4520,0,0,0,0,0],
		[1,0,0,2,0.1,2.0133,0,0,0,4529,0,0,0,0,0],
		[1,0,0,2,0.1,2.0133,0,0,0,4564,0,0,0,0,0],
		[64,0,0,2,0.1,1.9834,0,5,0,4611,0,0,0,0,0],
		[63,0,1,2,0.1,1.9834,0,5,0,4624,0,0,0,0,0],
		[0,0,0,2,0.1,0.0698,0,5,0,4645,0,0,0,0,0],
		[101,0,0,2,0.1,2.2824,0,0,0,4731,0,0,0,0,0],
		[0,0,0,2,0.1,0.1296,0,5,0,4905,0,0,0,0,0],
		[20,0,0,2,0.1,0.01,0,0,0,4914,0,0,0,0,0],
		[64,0,0,2,0.1,1.9834,0,0,0,4931,0,0,0,0,0],
		[1,0,0,2,0.1,0.0698,0,0,0,4952,0,0,0,0,0],
		[62,0,0,2,0.1,2.0133,0,19,0,4958,0,0,0,0,0],
		[62,0,0,2,0.1,2.0133,0,5,0,4964,0,0,0,0,0],
		[64,0,0,2,0.1,1.8339,0,5,0,5070,0,0,0,0,0],
		[21,0,0,2,0.1,0.1595,0,0,0,5112,0,0,0,0,0],
		[4,0,1,10,240.098,1.9834,0,55,0,5151,1677,0,0,0,0],
		[0,0,1,10,240.098,1.9834,0,2,0,5261,1832,0,0,0,0],
		[6,0,1,10,240.098,0.4585,0,57,0,5151,1677,0,0,0,0],
		[21,0,0,6,0.1,1.9236,0,59,0,5121,707,0,0,107,0],
		[21,0,0,6,0.1,1.9236,0,59,0,5121,1345,0,0,107,0],
		[21,0,0,2,480.096,0.1595,0,60,0,4912,0,0,0,118,1],
		[22,0,0,2,0.1,1.9236,5,60,0,5071,1217,0,0,113,0],
		[22,0,0,6,0.1,1.9236,5,59,0,4887,1177,0,0,114,0],
		[22,0,0,6,0.1,1.9236,5,59,0,4947,1267,0,0,115,0],
		[22,1,0,6,1080.09,1.9236,5,60,0,4950,1268,0,0,120,1],
		[26,32,0,9,1080.09,1.9834,11,170,159,4892,2285,6,88,227,3],
		[26,32,0,11,1080.09,1.9834,11,105,159,4952,2421,34,88,227,3],
		[24,0,0,10,720.094,1.8339,0,178,10,4906,1333,3,642,228,11],
		[25,0,0,10,720.094,0.01,0,178,10,4906,1333,4,643,228,11],
		[10,31,0,2,120.099,1.804,0,181,179,4955,3225,3,614,219,2],
		[25,0,0,14,720.094,0.01,0,178,171,4906,3318,4,643,228,11],
		[10,0,0,14,0.1,2.2525,0,146,253,4915,1765,11,641,213,3],
		[10,0,0,14,0.1,2.2525,12,184,253,4913,1554,11,733,249,3],
		[25,30,0,13,360.097,1.9834,11,202,1,4963,1791,1553,630,227,2],
		[11,30,16,13,1080.09,1.9834,11,181,175,4963,3318,1631,644,248,2],
		[26,31,15,10,1080.09,1.804,12,187,280,4954,3106,1541,651,253,2]
	]

	costs, utilities, designs = nsga2_obed_bn(	
					n_threads=threads, 
					prob=problem, 
					hours=do_hrs, 
					minutes=do_min, 
					popSize=popSize, 
					nSkip=2, 
					tolDelta=epsilon,#utility_conf95/1.96, 
					nPeriod=5, 
					nMonteCarlo=nMC, 
					GMM=gmm20, 
					Ylist=presampled_ylist,
					displayFreq=displayFreq,
					initial_pop=nmcp4_final
				)
	plot_nsga2(costs, utilities, design_pts, util_err=util_err, showPlot=True, savePlot=False, logPlotXY=[False,False])

if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--run', metavar='string', required=True, help='Function to run for this vvopt analysis')
	parser.add_argument('-n', type=int, default=0, help='Number of iterations to give to the function')
	parser.add_argument('--filename', metavar='string', default="SA_QoI", help='Base name to five to SA_QoI_sample')
	args = parser.parse_args()
	
	###Problem Definition
	_wave_min = 350.0
	_wave_max = 975.0
	_wave_redgreen = 690.0
	_wave_greenblue = 480.0
	_bandpass = _wave_max - _wave_min

	d_historical = [
						20,   #t_gain
						30,   #I_gain
						1,	#n_meas_rn
						8,	#d_num
						9600, #d_max
						2,	 #d_pow   #approx, 
						12,  #n_qe #questionable; kind of wasn't measured
						2,  #t_qe
						3, #d_vph_n_pts
						801, #d_dichroic_n_pts
						1501, #d_coll_n_pts #MIT Run 9-2730		
						13, #d_redcam_n_pts #protoLLAMAS camera testing
						0, #d_greencam_n_pts #protoLLAMAS camera testing
						13, #d_bluecam_n_pts #protoLLAMAS camera testing
						10 #d_frd_n_meas #from evaluating_cleaving_through_bigger.xlsx
					]
	d_0 = [
						0.1, #t_gain 				0-length exposure for gain exp
						0, #I_gain 				no measurements for gain exp
						0, #n_meas_rn 			no measurements for rn exp
						0, #d_num 				no measurements for dc exp
						1, #d_max 				dc exp filer
						0.01, #d_pow 				dc exp filler
						0, #n_qe 				no measurements for qe exp
						0.1, #t_qe 				0-length exposure for qe exp
						0, #d_vph_n_pts 		no measurements for vph exp
						0, #d_dichroic_n_pts 	no measurements for dichroic exp
						0, #d_coll_n_pts		no measurements for collimator exp
						0, #d_redcam_n_pts		no measurements for camera exp
						0, #d_greencam_n_pts		no measurements for camera exp
						0, #d_bluecam_n_pts		no measurements for camera exp
						0  #d_frd_n_meas 		no measurements for FRD exp
					]
	d_min = [
						0.1, #t_gain 				0-length exposure for gain exp
						1, #I_gain 				no measurements for gain exp
						1, #n_meas_rn 			no measurements for rn exp
						1, #d_num 				no measurements for dc exp
						1, #d_max 				dc exp filer
						0.01, #d_pow 				more short experiments
						1, #n_qe 				no measurements for qe exp
						0.1, #t_qe 				0-length exposure for qe exp
						1, #d_vph_n_pts 		no measurements for vph exp
						1, #d_dichroic_n_pts 	no measurements for dichroic exp
						1, #d_coll_n_pts		no measurements for collimator exp
						1, #d_redcam_n_pts		no measurements for camera exp
						1, #d_greencam_n_pts		no measurements for camera exp
						1, #d_bluecam_n_pts		no measurements for camera exp
						1  #d_frd_n_meas 		no measurements for FRD exp
					]

	d_max = [
						600, #t_gain
						100, #I_gain
						50, #n_meas_rn
						25, #d_num
						12000, #d_max
						3, #d_pow
						100, #n_qe
						600, #t_qe
						_bandpass*10, #d_vph_n_pts
						_bandpass*10, #d_dichroic_n_pts
						_bandpass*10, #d_coll_n_pts
						(_wave_max-_wave_redgreen)*10, #d_redcam_n_pts
						(_wave_redgreen-_wave_greenblue)*10, #d_greencam_n_pts
						(_wave_greenblue-_wave_min)*10, #d_bluecam_n_pts
						2400  #d_frd_n_meas
					]
	d_med = [dd/2 for dd in d_max]
	d_med[5] = 1 #d_pow
	
	problem = construct_slim_llamas_snr_problem()
	#instead of constructing the problem, im just constructing what nsga needs up top
	
	design_pts = [
		[problem.G(d_historical), 0.004240541527302059, "d_hist", 5.384503718405341e-05],
		[problem.G(d_0), 0.0047856764435419575, "d_0", 5.384503718405341e-05],
		[problem.G(d_max), 0.0022957744137691916, "d_max", 5.384503718405341e-05],
		[problem.G(d_med), 0.004294863943612242, "d_med", 5.384503718405341e-05],
		#[problem.G(d_min), 0.004772754162991483, "d_min", 5.384503718405341e-05]
	]
	
	if args.run == "OPT_test":
		vv_OPT(
			problem,
			gmm_file="ncomp_testing/BN_model_1639027_ncomp20.pkl", 
			ysamples_file="BN_samples_1639027.csv", 
			design_pts=design_pts,
			epsilon=0.01, #probably the ncomp20 model is better than this
			util_err=0.001,
			do_hrs = 0,
			do_min = 1,
			threads = 1 if args.n==0 else args.n,
			popSize=10,# if args.n==0 else args.n,
			nMC=10,
			displayFreq=10
		)

	elif args.run == "OPT_nmc_p4":
		conf95 = 0.0006962967583258008
		std_frac = conf95 / (1.96*(0.0047856764435419575 - 0.0022957744137691916))
		vv_OPT(
			problem,
			gmm_file="ncomp_testing/BN_model_1639027_ncomp200.pkl", 
			ysamples_file="BN_samples_1639027.csv", 
			design_pts=design_pts,
			epsilon=0.001,
			util_err=conf95,
			do_hrs = 11 if args.filename == "timed" else 0,
			do_min = 30 if args.filename == "timed" else 0,
			threads = 25,#1 if args.n == 0 else args.n,
			popSize=50,#30 if args.n==0 else args.n,
			nMC=10**4,
			displayFreq=5
		)

	elif args.run == "OPT_nmc_p5":
		conf95 = 0.00019389166166701476
		std_frac = conf95 / (1.96*(0.0047856764435419575 - 0.0022957744137691916))
		vv_OPT(
			problem,
			gmm_file="ncomp_testing/BN_model_1639027_ncomp200.pkl", 
			ysamples_file="BN_samples_1639027.csv", 
			design_pts=design_pts,
			epsilon=0.01,#std_frac, #rough analogue for nMC=10^5
			util_err=conf95,
			do_hrs = 11 if args.filename == "timed" else 0,
			do_min = 30 if args.filename == "timed" else 0,
			threads = 25,# if args.n == 0 else args.n,
			popSize=50,
			nMC=10**5,
			displayFreq=1
		)
		
	elif args.run == "OPT_nmc_p6":
		#prepare initializaiton samples based on previous run + prepared designs
	    
		conf95 = 5.384503718405341e-05
		std_frac = conf95 / (1.96*(0.0047856764435419575 - 0.0022957744137691916))
		vv_OPT(
			problem,
			gmm_file="ncomp_testing/BN_model_1639027_ncomp200.pkl", 
			ysamples_file="BN_samples_1639027.csv", 
			design_pts=design_pts,
			epsilon=std_frac, #the result for half the 95% confidence interval width for nMC=10^6
			util_err=conf95,
			do_hrs = 11 if args.filename == "timed" else 0,
			do_min = 30 if args.filename == "timed" else 0,
			threads = 8 if args.n == 0 else args.n,
			popSize=40,
			nMC=10**6,
			displayFreq=5
		)
	
	else:
		print("I dont recognize the command",args.run)
