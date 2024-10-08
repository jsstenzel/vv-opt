#This details an obed script that can be called regardless of the problem
#takes problem as input

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial

from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
#from pymoo.visualization.scatter import Scatter
#from pymoo.operators.crossover.pntx import TwoPointCrossover
#from pymoo.operators.mutation.bitflip import BitflipMutation
#from pymoo.operators.sampling.rnd import BinaryRandomSampling

#parallel:
from pymoo.core.problem import ElementwiseProblem
from multiprocessing.pool import ThreadPool
from pymoo.core.problem import *
from pymoo.factory import get_termination
from pymoo.termination.default import DefaultMultiObjectiveTermination
from pymoo.termination.robust import RobustTermination
from pymoo.termination.ftol import MultiObjectiveSpaceTermination

sys.path.append('..')
from problems.problem_definition import *

#remove this later, not appropriate, this should be outside and passed in from a wrapper fn
from obed.obed_gbi import *


def ngsa2_problem(prob, nGenerations=100, popSize=100, nMonteCarlo=10**5, nGMM=10**5):
	###1. Define the utility fn as an NGSA-II problem
	###Define the pymoo Problem class
	class VerificationProblem(Problem):
		param = {}
	
		def __init__(self, aux_matrix=None):
			d_lowers = [d_dist[1][0] for d_dist in prob.d_dists]
			d_uppers = [d_dist[1][1] for d_dist in prob.d_dists]
			param_types = [float if (dmask=="continuous") else int for dmask in prob.d_masks]

			#translate the key problem parameters here
			#2 objectives, utility and cost
			#For now, we won't consider that cost or utility or anything else are constrained
			super().__init__(n_var=prob.dim_d, n_obj=2, n_constr=0, xl=d_lowers, xu=d_uppers, vtype=param_types)
			
		def _evaluate(self, X, out, *args, **kwargs):
			print(flush=True, end="")
			utility_list = []
			cost_list = []

			for num, design in enumerate(X):
				#print(design, flush=True)
				#Construct dictionaries
				
				U, _ = U_varH_gbi(design, prob, n_mc=nMonteCarlo, n_gmm=nGMM, ncomp=5, doPrint=False)
				cost = np.log(prob.G(design))
				
				utility_list.append(U)
				cost_list.append(cost)
		
			out["F"] = np.column_stack([cost_list, utility_list])
			#out["G"] = np.column_stack([...,...,...])
			
	ngsa_problem = VerificationProblem()
	
	samples = [d for d in prob.sample_d(nGenerations-2)]
	samples.append([d_dist[1][0] for d_dist in prob.d_dists])
	samples.append([d_dist[1][1] for d_dist in prob.d_dists])	
	samples = np.array(samples)

	###2. Run NGSA-II
	algorithm = NSGA2(pop_size=popSize,
					  sampling=samples,
					  #crossover=TwoPointCrossover(), #what
					  #mutation=BitflipMutation(), #what
					  eliminate_duplicates=True)

	res = minimize(ngsa_problem, #minimize the Hvar, minimize the cost
				   algorithm,
				   #termination,
				   ('n_gen', nGenerations),
				   #seed=1,
				   verbose=True)
	
	#flip cost back
	pareto_costs = [f[0] for f in res.F]
	pareto_utilities = [f[1] for f in res.F]
	pareto_designs = res.X
	pareto_costs, pareto_utilities, pareto_designs = zip(*sorted(zip(pareto_costs, pareto_utilities, pareto_designs)))
	
	print(pareto_costs)
	print(pareto_utilities)
	print(res.X) #all of the pareto front points
	#print(res.history) #idk, returns None
		
	return pareto_costs, pareto_utilities, res.X
	
def plot_ngsa2(pareto_costs, pareto_utilities, design_pts, showPlot=False, savePlot=False, logPlotXY=[False,False]):
	costs, utilities = zip(*sorted(zip(pareto_costs, pareto_utilities)))
	
	plt.scatter(costs, utilities, s=80, facecolors='none', edgecolors='r')
	plt.plot(costs, utilities, c='k', alpha=0.7)
	plt.xlabel("log cost")
	plt.ylabel(r"Var $[p(Q|\mathbf{y})]$")
	
	if logPlotXY[0]:
		plt.xscale('log')
	if logPlotXY[1]:
		plt.yscale('log')
		
	for pt in design_pts:
		plt.scatter(pt[0],pt[1], s=80, facecolors='b', edgecolors='b')
		plt.annotate("  "+str(pt[2]), (pt[0], pt[1]))
		
	if savePlot:
		#name = "ngsa_" + str(nGenerations) + '_' + str(popSize) + '_' + str(nMonteCarlo) + '_' + str(nGMM) +'.png'
		plt.savefig("ngsa.png", format='png', bbox_inches='tight')
			
	if showPlot:
		plt.show()
		
	#plt.clf()
	#plt.close()	
	

def ngsa2_problem_parallel(n_threads, prob, hours, minutes, popSize, nMonteCarlo, nGMM):
	###1. Define the utility fn as an NGSA-II problem
	###Define the pymoo Problem class
	class VerificationProblemSingle(ElementwiseProblem):
		param = {}
	
		def __init__(self, **kwargs):
			d_lowers = [d_dist[1][0] for d_dist in prob.d_dists]
			d_uppers = [d_dist[1][1] for d_dist in prob.d_dists]
			param_types = [float if (dmask=="continuous") else int for dmask in prob.d_masks]

			#translate the key problem parameters here
			#2 objectives, utility and cost
			#For now, we won't consider that cost or utility or anything else are constrained
			super().__init__(n_var=prob.dim_d, n_obj=2, n_constr=0, xl=d_lowers, xu=d_uppers, vtype=param_types, **kwargs)
			
		def _evaluate(self, design, out, *args, **kwargs):
			print(flush=True, end="")
			utility_list = []
			cost_list = []

			#print(design, flush=True)
			#Construct dictionaries
			
			U, _ = U_varH_gbi(design, prob, n_mc=nMonteCarlo, n_gmm=nGMM, ncomp=5, doPrint=False)
			cost = np.log(prob.G(design))
			
			out["F"] = np.column_stack([cost, U])
			#out["G"] = np.column_stack([...,...,...])
			
	# initialize the thread pool and create the runner
	pool = ThreadPool(n_threads)
	runner = StarmapParallelization(pool.starmap)

	# define the problem by passing the starmap interface of the thread pool
	ngsa_problem = VerificationProblemSingle(elementwise_runner=runner)

	###2. Run NGSA-II
	algorithm = NSGA2(pop_size=popSize,
					  #sampling=BinaryRandomSampling(), #what
					  #crossover=TwoPointCrossover(), #what
					  #mutation=BitflipMutation(), #what
					  eliminate_duplicates=True)

	if hours==0 and minutes==0:
		#termination = DefaultMultiObjectiveTermination(
		#	f_tol=0.0025,
		#	period=30,
		#	n_max_gen=500
		#)	
		termination = RobustTermination(MultiObjectiveSpaceTermination(tol=0.05, n_skip=5), period=10)
	else:
		time_string = f"{hours:02}"+":"+f"{minutes:02}"+":00"
		print(time_string)
		termination=get_termination("time", time_string)
		
	res = minimize(ngsa_problem,
				   algorithm,
				   termination,
				   #seed=1,
				   verbose=True)
				   
	pool.close()
	
	pareto_costs = [f[0] for f in res.F]
	pareto_utilities = [f[1] for f in res.F]
	pareto_designs = res.X
	pareto_costs, pareto_utilities, pareto_designs = zip(*sorted(zip(pareto_costs, pareto_utilities, pareto_designs)))
	
	print(pareto_costs)
	print(pareto_utilities)
	print(res.X) #all of the pareto front points
	#print(res.history) #idk, returns None
		
	return pareto_costs, pareto_utilities, res.X
