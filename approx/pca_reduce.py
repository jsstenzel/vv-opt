import sys
import math
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA

sys.path.append('..')
from inference.bn_modeling import *

def standardize(data):
	print("Standardizing...",flush=True)
	data_mean = np.mean(data, axis=0)
	data_std = np.std(data, axis=0)
	data_standard = np.array([list((dd - data_mean)/data_std) for i,dd in enumerate(data)])
	return data_standard

def pca_analysis(problem, bn_samples, do_subset=0, doPlot=True):
	q, yd = bn_load_samples(problem, savefile=bn_samples, do_subset=do_subset, doPrint=True, doDiagnostic=True)
	data = [[qi]+ydi for qi,ydi in zip(q,yd)]
	ndim = len(data[0])
	print("Dimension of the ydq space:",ndim,flush=True)
	
	#standardize data
	data_standard = standardize(data)
	
	pca = PCA(n_components=ndim) #do PCA with full dimension, to decide how to chop things down
	pca.fit(data_standard)
	print(pca.explained_variance_)
	print(pca.explained_variance_ratio_)
	print(pca.singular_values_)
	
	if doPlot:
		# Plot explained variance ratio
		cumulative_variance_ratio = np.cumsum(pca.explained_variance_ratio_)
		plt.plot(cumulative_variance_ratio, marker='o', color='darkgreen')
		plt.xlabel('Number of Principal Components')
		plt.ylabel('Cumulative Explained Variance Ratio')
		plt.title('Cumulative Explained Variance Ratio by Principal Components')
		plt.show()
		
def pca_train(problem, bn_samples, ncomp, savefile, do_subset=0):
	q, yd = bn_load_samples(problem, savefile=bn_samples, do_subset=do_subset, doPrint=True, doDiagnostic=True)
	data = [[qi]+ydi for qi,ydi in zip(q,yd)]
	ndim = len(data[0])
	if ncomp >= ndim:
		print("Dimension of the ydq space:",ndim,"but ncomp =",ncomp,flush=True)
		sys.exit()
	
	#standardize data
	data_standard = standardize(data)
	
	pca = PCA(n_components=ncomp)
	pca.fit(data_standard)
	print(pca.explained_variance_)
	print(pca.explained_variance_ratio_)
	print(pca.singular_values_)
	
	with open(savefile, 'wb') as file:
		pickle.dump(pca, file)
	
	return pca
	
def pca_load(pca_file):
	filename = pca_file if pca_file.endswith('.pkl') else pca_file+'.pkl'
	
	with open(filename, 'rb') as file:
		pca = pickle.load(file)

	return pca