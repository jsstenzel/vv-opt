import sys
import os
import scipy.stats
import scipy.interpolate
import numpy as np
import matplotlib.pyplot as plt

import pymc

sys.path.append('../..')
#from obed.obed import *

#I need a nice simple interface to deal with these general messy things
#leave uncertainty handling to the prior fns!
class Functional:
	def __init__(self, *args, ):
		if len(args) == 1:
			xs = [point[0] for point in args[0]]
			ys = [point[1] for point in args[0]]
			self.size = len(args[0])
		elif len(args) == 2:
			xs = args[0]
			ys = args[1]
			if len(xs) != len(ys):
				raise ValueError("List lengths in Functional definition don't match!", len(xs), "and". len(ys))
			self.size = len(args[0])
		else:
			raise ValueError("Wrong number of Functional arguments!")

		self.xpts = xs
		self.ypts = ys
		
		self.xmin = xs[0]
		self.xmax = xs[-1]
		self.ymin = min(ys)
		self.ymax = max(ys)
		
		self.bspline = None
		
	def set_xlim(self,xmin,xmax):
		#risky
		self.xmin = xmin
		self.xmax = xmax
	
	def set_ylim(self,ymin,ymax):
		self.ymin = ymin
		self.ymax = ymax
	
	def spline_interp(self, order):
		self.bspline = scipy.interpolate.make_interp_spline(self.xpts, self.ypts, k=order, bc_type=None, axis=0, check_finite=True)
		
	def idx(self, i):
		return [self.xpts[i], self.ypts[i]]
		
	def f(self, x):
		#enforce x boundary condition strongly
		if x < self.xmin or x > self.xmax:
			raise ValueError("Functional evaluated outside domain!")
		if self.bspline == None:
			self.spline_interp(order=1)
			
		y = self.bspline(x)#extrapolate=None)
		#enforce y boundary condition softly
		if y > self.ymax:
			y = self.ymax
		elif y < self.ymin:
			y = self.ymin
		return float(y)
		
	def to_array(self, steps):
		#convert to array
		list_funct = []
		for _x in np.linspace(self.xmin, self.xmax, steps):
			_y = self.f(_x)
			#enforce boundary condition
			if _y > self.ymax:
				_y = self.ymax
			elif _y < self.ymin:
				_y = self.ymin
			list_funct.append([_x, _y])
		return list_funct
		
	def plot(self, linecolor='k', pointcolor='k', numpoints=10000, show=True):
		plt.xlim(self.xmin, self.xmax)
		plt.ylim(self.ymin, self.ymax)
		
		#scatter the data
		plt.scatter(self.xpts, self.ypts, c=pointcolor)
		
		#plot the bspline if there is one
		if self.bspline != None:
			_x = np.linspace(self.xmin, self.xmax, numpoints)
			_y = [self.f(point) for point in _x]
			plt.plot(_x, _y, c=linecolor)
		
		if show:
			plt.show()


def noise_to_functional(funct, noise):
	xs = funct.xpts
	ys = [y + scipy.stats.norm.rvs(size=1, loc=0, scale=noise)[0] for y in funct.ypts]
	
	sample = Functional(xs, ys)
	sample.set_xlim(funct.xmin, funct.xmax)
	sample.set_ylim(funct.ymin, funct.ymax)
	sample.spline_interp(funct.bspline.k)
	
	return sample
	
	
	
########## Now I need to do Gaussian process
#https://www.pymc.io/projects/docs/en/v3/Gaussian_Processes.html

if __name__ == "__main__":
	with pymc.Model() as model:
		cov_func = pymc.gp.cov.ExpQuad(input_dim=1, ls=0.1)
		theta_gp = pymc.gp.Latent(cov_func=cov_func)
		
		X = np.linspace(0, 1, 10)[:, None]
		f = theta_gp.prior("f", X)
		print(f)