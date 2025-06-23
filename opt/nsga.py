#This details an obed script that can be called regardless of the problem
#TODO someday rework and generalize these methods to be independentof OBED approach, and even independent of the thing to be optimized

import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from functools import partial
from tabulate import tabulate

from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.core.population import Population
#from pymoo.visualization.scatter import Scatter
#from pymoo.operators.crossover.pntx import TwoPointCrossover
from pymoo.operators.mutation.pm import PM
from pymoo.operators.mutation.rm import ChoiceRandomMutation
#from pymoo.operators.sampling.rnd import BinaryRandomSampling
from pymoo.util.display.display import Display
from pymoo.util.display.multi import MultiObjectiveOutput
from pymoo.operators.repair.rounding import RoundingRepair

#parallel:
from pymoo.core.problem import ElementwiseProblem
from multiprocessing import Pool
from pymoo.core.problem import *
from pymoo.factory import get_termination
from pymoo.termination.default import DefaultMultiObjectiveTermination
from pymoo.termination.robust import RobustTermination
from pymoo.termination.ftol import MultiObjectiveSpaceTermination

sys.path.append('..')
from problems.problem_definition import *

#remove this later, not appropriate, this should be outside and passed in from a wrapper fn
from obed.obed_gbi import *

""" TODO fold this in as an option for nsga2_problem_parallel
def nsga2_problem(prob, nGenerations=100, popSize=100, nMonteCarlo=10**5, nGMM=10**5):
	###1. Define the utility fn as an nsga-II problem
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
			
	nsga_problem = VerificationProblem()
	
	samples = [d for d in prob.sample_d(nGenerations-2)]
	samples.append([d_dist[1][0] for d_dist in prob.d_dists])
	samples.append([d_dist[1][1] for d_dist in prob.d_dists])	
	samples = np.array(samples)

	###2. Run nsga-II
	algorithm = NSGA2(pop_size=popSize,
					  sampling=samples,
					  #crossover=TwoPointCrossover(), #what
					  #mutation=BitflipMutation(), #what
					  eliminate_duplicates=True)

	res = minimize(nsga_problem, #minimize the Hvar, minimize the cost
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
"""
	
def plot_nsga2(pareto_costs, pareto_utilities, design_pts, util_err=None, showPlot=False, savePlot=False, logPlotXY=[False,False]):
	costs, utilities = zip(*sorted(zip(pareto_costs, pareto_utilities)))
	
	if not util_err:
		plt.scatter(costs, utilities, s=80, facecolors='none', edgecolors='r')
	else:
		plt.errorbar(costs, utilities, yerr=util_err, fmt='.', c='r', capsize=5)
	plt.plot(costs, utilities, c='k', alpha=0.7)
	plt.xlabel("total testbed time [s]")
	plt.ylabel(r"Var $[p(Q|\mathbf{y})]$")
	
	if logPlotXY[0]:
		plt.xscale('log')
	if logPlotXY[1]:
		plt.yscale('log')
		
	for pt in design_pts:
		plt.annotate("  "+str(pt[2]), (pt[0], pt[1]))		
		if not len(pt) == 4:
			plt.scatter(pt[0],pt[1], s=80, facecolors='b', edgecolors='b')
		else:
			plt.errorbar(pt[0], pt[1], yerr=pt[3], fmt='.', c='b', capsize=5)
		
	if savePlot:
		#name = "nsga_" + str(nGenerations) + '_' + str(popSize) + '_' + str(nMonteCarlo) + '_' + str(nGMM) +'.png'
		plt.savefig("nsga.png", format='png', bbox_inches='tight')
			
	if showPlot:
		plt.show()
		
	#plt.clf()
	#plt.close()	
	
def design_print(costs, utils, designs, prob):
	headers = ["Cost", "Utility"]+list(prob.d_names)
	#print(res.history) #idk, returns None
	table = []
	for cost, util, design in zip(costs, utils, designs):
		design_exact = prob._mask(design, prob.d_masks, prob.d_dists)
		row = [cost] + [util] + design_exact
		table.append(row)
	print(tabulate(table,headers),flush=True)

#This function works for LLAMAS problem
"""
Optimization parameters:
n_threads, 
hours, 
minutes, 

Genetic Algorithm parameters:
popSize - number of running samples in the GA
nSkip - number of generations that are skipped between evals of the optimization stopping criterion
tolDelta -  The tolerance of the optimization stopping criterion, representing the maximum threshold difference of solution sets across generations
nPeriod - The number of generations that are considered for each evaluation of the optimization stopping criterion

OBED parameters:
prob - the vv-opt problem
nMonteCarlo - number of Monte Carlo iterations per U(d) evaluation 
GMM - the Bayesian network model which is conditioned on y and d for Q
Ylist - a list of presampled y data to pull from in Monte Carlo simulation
"""
# Define a custom display class to print the population at each generation
class PeriodicPopulationDisplay(Display):
	def __init__(self, displayFreq, prob, output=MultiObjectiveOutput(), progress=False, verbose=True):
		self.displayFreq = displayFreq
		self.prob = prob
		super().__init__(output=output, progress=progress, verbose=verbose)

	def update(self, algorithm, **kwargs):
		super().update(algorithm, **kwargs)
		if self.displayFreq != 0:
			if algorithm.n_gen % self.displayFreq == 0 or algorithm.n_gen == 1:  # Print every generation (change to desired frequency)
				F = algorithm.pop.get('F')
				costs = [f[0] for f in F]
				utilities = [f[1] for f in F]
				print(f"Generation {algorithm.n_gen}:")
				design_print(costs, utilities, algorithm.pop.get('X'), self.prob)
				print("-" * 20)
					
class VerificationProblemSingle(ElementwiseProblem):
	param = {}

	def __init__(self, prob, GMM, Ylist, nMonteCarlo, **kwargs):
		self.prob = prob
		self.GMM = GMM
		self.Ylist = Ylist
		self.nMonteCarlo = nMonteCarlo
	
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
		
		U, _ = U_varH_gbi_joint_presampled(design, None, gmm_qyd=self.GMM, presampled_ylist=self.Ylist, n_mc=self.nMonteCarlo, doPrint=False)
		cost = self.prob.G(design)
		
		out["F"] = np.column_stack([cost, U])
		#out["G"] = np.column_stack([...,...,...])
			
def nsga2_obed_bn(n_threads, prob, hours, minutes, popSize, nSkip, tolDelta, nPeriod, nMonteCarlo, GMM, Ylist, displayFreq=10, initial_pop=[], doMutation="default"):	
	# initialize the thread pool and create the runner
	pool = Pool(n_threads)
	runner = StarmapParallelization(pool.starmap)

	# define the problem by passing the starmap interface of the thread pool
	nsga_problem = VerificationProblemSingle(prob=prob, GMM=GMM, Ylist=Ylist, nMonteCarlo=nMonteCarlo, elementwise_runner=runner)

	#set parameters
	if doMutation=="RM":
		mutation = PM(prob=1.0, eta=1.0, vtype=int, repair=RoundingRepair()) #ChoiceRandomMutation()	
	else:
		mutation = PM(eta=20)

	###2. Run nsga-II
	#Allow initialization with a pre-sampled population
	#https://pymoo.org/customization/initialization.html
	#it can be a different size from nmc, but not smaller if we're eliminating duplicates
	if len(initial_pop) == 0:
		algorithm = NSGA2(pop_size=popSize,
						  #sampling=BinaryRandomSampling(), #what
						  #crossover=TwoPointCrossover(), #what
						  mutation=mutation,
						  eliminate_duplicates=True)
	elif len(initial_pop) < popSize:
		#draw more random samples so that we're up to pop = nmc
		diff = popSize - len(initial_pop)
		random_sample = list(prob.sample_d(diff)) if diff>1 else [prob.sample_d(diff)]
		greater_pop = initial_pop + random_sample
		pop_obj = Population.new("X", greater_pop)
		algorithm = NSGA2(pop_size=popSize,
			mutation=mutation,
			sampling=pop_obj,
			eliminate_duplicates=True)
	elif len(initial_pop) == popSize:
		pop_obj = Population.new("X", initial_pop)
		algorithm = NSGA2(pop_size=popSize,
			mutation=mutation,			
			sampling=pop_obj,
			eliminate_duplicates=True)
	else: #len(initial_pop) > popSize
		pop_obj = Population.new("X", initial_pop)
		algorithm = NSGA2(pop_size=popSize,
			mutation=mutation,
			sampling=pop_obj,
			eliminate_duplicates=False)

	if hours==0 and minutes==0:
		#termination = DefaultMultiObjectiveTermination(
		#	f_tol=0.0025,
		#	period=30,
		#	n_max_gen=500
		#)	
		termination = RobustTermination(MultiObjectiveSpaceTermination(tol=tolDelta, n_skip=nSkip), period=nPeriod)
	else:
		time_string = f"{hours:02}"+":"+f"{minutes:02}"+":00"
		print(time_string)
		termination=get_termination("time", time_string)
		
	res = minimize(nsga_problem,
				   algorithm,
				   termination,
				   #seed=1,
				   display=PeriodicPopulationDisplay(displayFreq, prob),
				   verbose=True)
				   
	pool.close()
	
	pareto_costs = [f[0] for f in res.F]
	pareto_utilities = [f[1] for f in res.F]
	pareto_designs = res.X
	pareto_costs, pareto_utilities, pareto_designs = zip(*sorted(zip(pareto_costs, pareto_utilities, pareto_designs)))
	
	#do a final print of the optimal designs
	design_print(pareto_costs, pareto_utilities, pareto_designs, prob)
		
	return pareto_costs, pareto_utilities, pareto_designs
